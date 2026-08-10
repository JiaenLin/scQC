"""The task graph: which adapter runs when, and which gate judges the result.

This is the seam between the two halves of scQC. The adapters produce numbers and know nothing
about policy; the modules judge numbers and never touch data. Every join between them is made
here and nowhere else, so there is exactly one place to look when asking "what actually decides
whether this run continues".

THE ORDER IS NOT A PREFERENCE

Doublet scoring precedes quality filtering because scDblFinder's documentation requires it. The
light floor precedes doublet scoring because a detector given near-empty droplets models noise.
The cell-call gate precedes everything downstream because a population lost there cannot be
recovered by any later step. Each of those is a constraint, not a taste, and reordering them
silently changes what the deliverable means.

WHAT A STEP MAY NOT DO

A step may not decide policy. If a threshold is needed, it comes from `decisions` or from a
module's `derive`; a step that picks its own number has moved an adjudicated parameter into code
where nobody will look for it.
"""

from __future__ import annotations

from pathlib import Path

from .pipeline import step_module
from .task import Refusal, Task, TaskFailure


# --------------------------------------------------------------------------------------------
# helpers


def _tables(p) -> Path:
    return p.results / "tables"


def _objects(p) -> Path:
    return p.results / "objects"


def _tool(p, name: str, default: str) -> str:
    """Resolve a tool path from the project config, falling back to the bare name on PATH."""
    return str(p.decisions.get("tools", {}).get(name, default))


def _design(samples: list[dict], max_levels: int = 6) -> dict:
    """Design factors discovered from the samplesheet, never declared.

    Three exclusions, and the third is the one that is easy to miss:

      * a column that is CONSTANT separates nothing;
      * a column with more than `max_levels` levels is not a factor this pipeline can test;
      * a column with ONE SAMPLE PER LEVEL is an identifier, not a factor. A replicate id or a
        library barcode in a four-sample cohort has four levels, passes a bare `<= 6` test, and
        then every differential check computes a ratio between single libraries - arithmetic
        with no evidence in it, reported in the same words as a real design differential.

    So a factor must leave at least one level holding more than one sample.
    """
    skip = {"sample", "platform", "species", "reference", "assay",
            "fastq_r1", "fastq_r2", "matrix"}
    out: dict = {}
    if not samples:
        return out
    n = len(samples)
    ceiling = min(max_levels, max(n - 1, 1))
    for col in samples[0]:
        if col in skip:
            continue
        vals = {r.get(col) for r in samples if str(r.get(col) or "").strip()}
        if 2 <= len(vals) <= ceiling:
            out[col] = {r["sample"]: r[col] for r in samples if r.get(col)}
    return out


# --------------------------------------------------------------------------------------------
# step 0 - ingest


def _ingest(task, pipeline, log):
    from adapters import matrix as mx

    row = task.params["row"]
    sample = row["sample"]
    ing = step_module("ingest")
    registry = ing.read_registry(pipeline.project.parent / "references" / "_registry"
                                 / "registry.tsv")
    if not registry:
        registry = ing.read_registry(Path(__file__).resolve().parents[1]
                                     / "references" / "_registry" / "registry.tsv")

    errs = ing.validate_row(row, registry)
    if errs:
        raise Refusal(f"samplesheet row for {sample} is not usable:\n"
                      + "\n".join(f"    - {e}" for e in errs))


    # The adapter supplies this callable precisely so the calling convention lives in one place:
    # plan_one() calls stats_fn(path) and then verify(name=path, **result), so the result must be
    # verify's keyword arguments WITHOUT `name` and without the extra keys summary_stats returns.
    # Assembling it here instead duplicated that convention and got it wrong - the returned
    # payload was discarded and the JSON re-read under a different shape, so step 0 could not run.
    #
    # MEASURED IN THE ANALYSIS INTERPRETER, NOT THIS ONE.
    #
    # `ingest_stats_fn` measures in-process, which requires pandas and anndata in whatever
    # interpreter is running the orchestrator. That is the one thing this pipeline can never
    # assume: the aligner, the denoiser and the analysis stack have incompatible pins and live in
    # separate environments, so `scqc` itself routinely runs under a bare system python. Step 0
    # died on `No module named 'pandas'` for exactly that reason - the adapter had the
    # out-of-process route all along (build_stats_argv / run_summary_stats) and the engine
    # called the in-process one.
    #
    # `--python` names the interpreter that has the analysis stack; it is required here rather
    # than defaulted, because falling back to this process would make step 0 succeed on a
    # developer laptop and fail on every cluster.
    python_exe = task.params.get("python_exe")
    if not python_exe:
        raise Refusal(
            "00_ingest: no analysis interpreter was given. Measuring a matrix needs pandas and "
            "anndata, which the process running scqc is not required to have. Pass --python.")

    def stats_fn(matrix_path):
        res = mx.run_summary_stats(
            matrix=matrix_path,
            out_json=pipeline.work / f"{sample}_ingest_stats.json",
            log=pipeline.work / f"{sample}_ingest_stats.log",
            python_exe=python_exe,
            expected_genes=task.params.get("expected_genes"),
            tmp_dir=pipeline.work / f"{sample}_extract",
            executor=pipeline.executor)
        # run_summary_stats returns {'outputs','metrics','versions'}; plan_one wants verify()'s
        # keyword arguments. verify_kwargs is the same selector the in-process path uses, so the
        # two routes cannot drift into accepting different keys.
        return mx.verify_kwargs(res["metrics"], include_name=False)

    plan = ing.plan_one(row, registry, stats_fn)
    print(f"    {plan}")
    if plan.mode == "blocked":
        raise Refusal(f"{sample} cannot be ingested: {plan.reason}\n"
                      f"    A blocked sample is not skipped - a cohort missing a library is a "
                      f"different cohort, and continuing would not say so.")
    # Step 0 decides; it writes nothing of its own, so it promises no outputs.
    outs: list = []
    return {"outputs": outs,
            "metrics": {"mode": plan.mode, "processor": plan.processor,
                        "reason": plan.reason},
            "versions": {}}


# --------------------------------------------------------------------------------------------
# step 1 - ambient, and its cohort audit


def _ambient(task, pipeline, log):
    from adapters import cellbender as cbd

    sample = task.sample
    raw = task.params["raw"]
    out_h5 = _objects(pipeline) / f"{sample}_cellbender.h5"
    res = cbd.run_remove_background(
        sample=sample, input_path=raw, output_h5=out_h5,
        exe=task.params["exe"], env_bin=task.params.get("env_bin"),
        fpr=task.params.get("fpr", 0.0),
        learning_rate=task.params.get("learning_rate"),
        device=task.params.get("device", "cuda"),
        log=log, executor=pipeline.executor)
    return res


def _ambient_supplied(task, pipeline, log):
    """Validate and place matrices that were denoised elsewhere. Corrects nothing.

    The provenance goes through `plan_ambient(..., supplied=...)`, the same function that decides
    whether CellBender runs, so an unattributed object cannot enter the pipeline through a side
    door that the ordinary path would have refused.

    The object is COPIED to where step 1's output would have been rather than referenced in
    place. A symlink into someone else's results directory makes this run's deliverable depend on
    a file this run does not own, and the provenance record then describes a matrix that can
    change underneath it.
    """
    import shutil

    am = step_module("ambient")
    supplied = task.params["supplied"]
    assay = task.params["assay"]
    outs, records = [], {}
    for s, rec in sorted(supplied.items()):
        src = Path(rec["path"])
        if not src.exists():
            raise Refusal(
                f"01_ambient ({s}): the supplied ambient matrix {src} does not exist. A missing "
                f"input is not an uncorrected sample - it is a run that cannot start.")
        try:
            plan = am.plan_ambient(s, assay.get(s), supplied={
                k: rec[k] for k in ("tool", "version", "params", "produced_by")})
        except Exception as e:                                            # noqa: BLE001
            raise Refusal(f"01_ambient ({s}): {e}") from None
        print(f"    {plan}")
        dst = _objects(pipeline) / f"{s}_ambient.h5"
        if dst.resolve() != src.resolve():
            shutil.copy2(src, dst)
        outs.append(str(dst))
        records[s] = {**plan.supplied, "source": str(src), "state": plan.state}

    import json
    p = _tables(pipeline) / "ambient_supplied.json"
    p.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"outputs": outs + [str(p)],
            "metrics": {"supplied": len(records), "corrected_here": 0},
            "versions": {}}


def _ambient_audit(task, pipeline, log):
    from adapters import cellbender as cbd

    aa = step_module("audit_ambient")
    supplied = task.params.get("supplied") or {}
    rows, per_gene = [], []
    for s, r in task.params["per_sample"].items():
        m = cbd.parse_metrics(r["h5"], r["raw"])
        rows.append({"sample": s,
                     "fraction_removed_overall": m["fraction_removed_overall"],
                     "genes_fully_removed": m["genes_fully_removed"]})
        for g in m.get("per_gene", []):
            per_gene.append({"sample": s, **g})

    # No pandas. This runs in the ORCHESTRATOR's interpreter, which on a cluster is a bare
    # python - the analysis stack lives in --python and the denoiser in its own environment.
    # `import pandas` here killed step 1 with ModuleNotFoundError on a cohort where the audit
    # had nothing to compute anyway. modules/01_ambient/audit_ambient.py takes plain rows.
    import csv
    out = _tables(pipeline) / "ambient_summary.csv"
    cols = ["sample", "fraction_removed_overall", "genes_fully_removed"]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})

    if supplied:
        # NOT audited, and said so rather than passed over. Measuring what ambient correction
        # removed needs BOTH the raw and the denoised counts; a supplied matrix arrives without
        # its raw, so the fraction removed is not merely unknown to this run - it is
        # unmeasurable by it. An audit that silently covers 0 libraries and reports "ok" is the
        # worst available outcome, because it reads exactly like an audit that passed.
        names = ", ".join(sorted(supplied))
        print(f"    NOT AUDITED: {len(supplied)} library(ies) arrived already corrected "
              f"({names}).")
        print("      The fraction removed cannot be measured without the raw counts. Their "
              "provenance is in tables/ambient_supplied.json.")
        if not rows:
            pipeline.gate("01_ambient", [], "NOT RUN")
            return {"outputs": [str(out)],
                    "metrics": {"libraries": 0, "supplied_not_audited": len(supplied)},
                    "versions": {}}

    findings = aa.audit(rows, per_gene, task.params["design"])
    pipeline.gate("01_ambient", findings, aa.verdict(findings))
    return {"outputs": [str(out)],
            "metrics": {"libraries": len(rows), "supplied_not_audited": len(supplied)},
            "versions": {}}


# --------------------------------------------------------------------------------------------
# step 2 - the cell-call gate


def _cellcall(task, pipeline, log):
    """Compare the aligner's cell call with the denoiser's, over the barcodes themselves.

    `calls` was read straight from task.params and nothing ever put it there, so this step had
    never run. It is built here from each caller's barcode list, because `lost` - aligner cells
    the denoiser did not call - is a SET DIFFERENCE. Two callers can agree on a total and
    disagree about which cells; the gate turns on the difference, not the totals.
    """
    from adapters import matrix as mx

    cg = step_module("cellcall_gate")
    paths = task.params.get("call_paths") or {}
    calls, missing = {}, []
    for s in task.params["samples"]:
        p = paths.get(s) or {}
        a_path, c_path = p.get("aligner"), p.get("cellbender")
        if not a_path or not c_path:
            missing.append(
                f"{s}: " + ", ".join(
                    x for x in (("aligner_cells" if not a_path else None),
                                ("cellbender_barcodes" if not c_path else None)) if x))
            continue
        a = set(mx.called_barcodes(a_path))
        c = set(mx.called_barcodes(c_path))
        calls[s] = {"aligner": len(a), "cellbender": len(c), "lost": len(a - c)}
        print(f"    {s:<14} aligner {len(a):>7,}   denoiser {len(c):>7,}   "
              f"lost {len(a - c):>6,}  ({100 * len(a - c) / max(len(a), 1):.2f}%)")

    if missing:
        # REFUSED, not skipped. This step exists to catch a population lost between two callers,
        # and a population lost here cannot be recovered by anything downstream. Running the rest
        # of the pipeline with the check silently absent is the outcome it was written to prevent.
        raise Refusal(
            "02_cells cannot compare cell calls - the samplesheet does not say where they are:\n"
            + "\n".join(f"    - {m}" for m in missing)
            + "\n    `aligner_cells` is the aligner's filtered matrix directory (CeleScope "
              "outs/filtered, CellRanger filtered_feature_bc_matrix).\n"
              "    `cellbender_barcodes` is CellBender's <stem>_cell_barcodes.csv; it is required "
              "when the denoised object was SUPPLIED, because this run never produced one.")

    findings = cg.gate(calls, task.params["design"])
    pipeline.gate("02_cells", findings, cg.verdict(findings))
    out = _tables(pipeline) / "cell_calls.csv"
    import csv
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "aligner", "ambient", "lost"])
        for s, c in sorted(calls.items()):
            w.writerow([s, c["aligner"], c["cellbender"], c["lost"]])
    return {"outputs": [str(out)], "metrics": {"libraries": len(calls)}, "versions": {}}


# --------------------------------------------------------------------------------------------
# step 5 - thresholds derived from measured valleys


def _quality(task, pipeline, log):
    q = step_module("quality")
    valleys = task.params["valleys"]        # [{"sample":..,"value":..,"bimodal":..}, ...]
    metric = task.params["metric"]
    objs = [q.Valley(v["sample"], metric, float(v["value"]), bool(v["bimodal"]))
            for v in valleys]
    try:
        prop = q.derive(objs, metric, light_floor=task.params.get("light_floor"))
    except Exception as e:                                            # noqa: BLE001
        # A refusal here is a verdict about the valleys, and it stops the run - but it is
        # reported as a refusal, not as a crash, because the two mean different things to a
        # reader and only one of them is a bug.
        raise Refusal(f"05_quality ({metric}): {e}") from None
    print(f"    {prop}")
    out = _tables(pipeline) / f"valleys_{metric}.csv"
    import csv
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "metric", "valley", "bimodal"])
        for v in valleys:
            w.writerow([v["sample"], metric, v["value"], v["bimodal"]])
    return {"outputs": [str(out)],
            "metrics": {"proposed": prop.constant, "spread": prop.spread},
            "versions": {}}


# --------------------------------------------------------------------------------------------
# the report - always last, and always produced


def _report(task, pipeline, log):
    from report.build import build_report

    payload = pipeline.report_payload(stopped=None)
    payload.update(task.params.get("extra", {}))
    out_html = pipeline.results / "reports" / "qc_report.html"
    out_json = pipeline.results / "reports" / "report.json"
    build_report(payload, out_html, out_json)
    return {"outputs": [str(out_html), str(out_json)],
            "metrics": {"findings": len(pipeline.findings)}, "versions": {}}


# --------------------------------------------------------------------------------------------
# graph construction


def build_tasks(pipeline, python_exe: str, tools: dict) -> list[Task]:
    """Assemble the graph for this project and mode.

    In evidence mode the apply task is NOT PLACED IN THE GRAPH. Not disabled, not guarded by a
    flag - absent, so there is no code path from `--mode evidence` to a deletion.
    """
    tasks: list[Task] = []
    samples = pipeline.samples
    design = _design(samples)
    if not design:
        print("  NOTE: no design factor found in the samplesheet. Every differential check "
              "will report NOT CHECKED,\n        which is its own outcome and does not read "
              "as a pass.")

    for row in samples:
        s = row["sample"]
        tasks.append(Task(
            key=f"00_ingest/{s}", step="00_ingest", sample=s, fn=_ingest,
            inputs=tuple(x for x in (row.get("matrix"), row.get("fastq_r1")) if x),
            params={"row": row, "python_exe": python_exe,
                    "expected_genes": row.get("expected_genes")},
        ))

    # The remaining per-sample steps are added by the caller once step 0 has decided whether a
    # matrix is accepted or must be rebuilt: a graph that assumes the answer would either skip a
    # needed alignment or run one that was not needed.
    tasks.append(Task(
        key="report", step="report", fn=_report,
        needs=tuple(t.key for t in tasks),
        params={"extra": {}},
    ))
    return tasks


# --------------------------------------------------------------------------------------------
# step 0b - alignment, when step 0 decided a matrix must be rebuilt


def _align(task, pipeline, log):
    row = task.params["row"]
    proc = (task.params.get("processor") or "").lower()
    tools = task.params["tools"]
    work = pipeline.work / f"{task.sample}_align"
    if proc == "celescope":
        from adapters import celescope as cs
        # `threads`, not `thread`. The two adapters were called with keywords neither accepts,
        # so this branch raised TypeError before any of its guards could run - the alignment
        # path was dead and every check inside it unreachable.
        return cs.run_celescope(
            sample=task.sample, fastq_dir=row["fastq_r1"], genome_dir=row["reference"],
            chemistry=row.get("chemistry"), work_dir=work, log=log,
            threads=task.cpus, env_bin=tools.get("celescope_bin"),
            exe=tools.get("celescope", "multi_rna"), executor=pipeline.executor)
    if proc == "cellranger":
        from adapters import cellranger as cr
        return cr.run_cellranger(
            sample=task.sample, fastq_dirs=[row["fastq_r1"]], transcriptome=row["reference"],
            work_dir=work, log=log, localcores=task.cpus, localmem_gb=task.memory_gb,
            exe=tools.get("cellranger", "cellranger"), executor=pipeline.executor)
    raise TaskFailure(
        f"{task.sample}: no processor for platform {row.get('platform')!r}. Step 0 decided this "
        f"sample must be rebuilt from FASTQ but named no tool to do it with, which is a defect "
        f"in the plan rather than an option missing here.")


# --------------------------------------------------------------------------------------------
# step 4 - doublet scoring, and its cohort health check


def _doublets(task, pipeline, log):
    from adapters import doublets as db

    p = task.params
    mtx = pipeline.work / f"{p['sample']}_dbl_mtx"
    out_csv = _tables(pipeline) / f"{p['sample']}_doublets.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    export = db.export_matrix(p["h5"], mtx, min_umi=int(p["light_floor"]))
    # The never-scored population, named from what ExportedMatrix actually carries.
    #
    # This was `getattr(export, "unscored", None)`, and ExportedMatrix has no `.unscored` - it
    # has `below_floor` and `not_selected`, two different reasons a barcode was never handed to
    # the detector, deliberately kept apart. So the expression evaluated to None on every run and
    # the never-scored population was never recorded at all: DoubletCalls.unscored stayed None,
    # which means "nobody said what was left out" and supports strictly fewer claims than the
    # truth did.
    #
    # Both reasons are unioned here because both end up UNKNOWN, which is what the doublet rate's
    # denominator needs. They stay separable on `export` for anyone asking which threshold did it.
    unscored = tuple(export.below_floor) + tuple(export.not_selected)
    return db.run_scdblfinder(
        rscript=p["rscript"], mtx_dir=mtx, out_csv=out_csv,
        dbr=p.get("dbr"), dbr_sd=p.get("dbr_sd"), seed=int(p.get("seed", 0)),
        log=log, executor=pipeline.executor, sample=p["sample"],
        unscored=unscored)


def _doublet_health(task, pipeline, log):
    from adapters import doublets as db

    dh = step_module("doublet_health")
    rates, unscored_total = {}, 0
    for s in task.params["samples"]:
        calls = db.read_calls(_tables(pipeline) / f"{s}_doublets.csv", sample=s)
        # TWO BUGS LIVED IN THE TWO LINES THIS REPLACES, and each hid the other.
        #
        # `getattr(calls, "scored", None) or {}` - DoubletCalls has no `.scored`; it IS the
        # mapping. The default turned a wrong attribute name into an empty result, so every
        # library reported no rate and the step refused with "no library produced a doublet
        # rate". The refusal was honest and the cause was invisible: a missing attribute must
        # raise, not resolve to nothing.
        #
        # `sum(1 for v in scored.values() if v)` - each value is a (score, is_doublet) TUPLE, and
        # a non-empty tuple is always truthy. Had the first bug not emptied the mapping, this
        # would have counted every scored barcode as a doublet and reported a 100% rate. The
        # health module would then have refused for "rate variance at zero", which is a real
        # check firing on a fabricated number - the worst of the available outcomes.
        n_scored = len(calls)
        n_pos = sum(1 for _score, is_doublet in calls.values() if is_doublet)
        # A barcode below the light floor was never scored. It is UNKNOWN, counted as unknown
        # rather than folded into the denominator as a negative.
        rates[s] = (n_pos / n_scored) if n_scored else None
        unscored_total += len(calls.unscored or ()) if calls.unscored is not None else 0
        print(f"    {s:<14} scored {n_scored:>7,}   doublets {n_pos:>6,}   "
              f"{100 * n_pos / max(n_scored, 1):.2f}%")

    known = {s: r for s, r in rates.items() if r is not None}
    if not known:
        raise Refusal("04_doublets: no library produced a doublet rate, so the health check was "
                      "not run. NOT CHECKED is its own outcome and does not read as a pass.")
    findings = dh.health(known, task.params["design"], n_kept_unscored=unscored_total,
                         detector_name="scDblFinder", reproducible=True)
    pipeline.gate("04_doublets", findings, dh.verdict(findings))
    return {"outputs": [], "metrics": {"libraries": len(known),
                                       "never_scored": unscored_total}, "versions": {}}


# --------------------------------------------------------------------------------------------
# step 5 - measure a valley per library, propose one cohort constant


def _scanpy(pipeline, op, h5, prefix, params, log, python_exe):
    from adapters import scanpy_ops as so
    return so.run_scanpy_op(op=op, h5ad_in=h5, out_prefix=prefix, params=params,
                            log=log, executor=pipeline.executor, python_exe=python_exe)


def _quality_stage(task, pipeline, log):
    valleys, mito_stats = [], {}
    for s in task.params["samples"]:
        res = _scanpy(pipeline, "valley",
                      pipeline.results / "objects" / f"{s}_ambient.h5",
                      pipeline.work / f"{s}_qc", {"metric": "umi"},
                      log, task.params["python_exe"])
        m = res.get("metrics", {})
        valleys.append({"sample": s, "value": m.get("valley"), "bimodal": m.get("bimodal")})
        q = m.get("mito_quartiles")
        if q:
            mito_stats[s] = q

    unknown = [v["sample"] for v in valleys if v["value"] is None or v["bimodal"] is None]
    if unknown:
        raise Refusal(
            f"05_quality: no valley was established for {', '.join(unknown)}. A library with no "
            f"measured valley cannot contribute to a cohort constant, and treating its absence "
            f"as agreement lets the other libraries decide on its behalf.")
    sub = Task(key=task.key, step=task.step, fn=_quality,
               params={"valleys": valleys, "metric": "umi",
                       "light_floor": task.params.get("light_floor")})
    out = _quality(sub, pipeline, log)
    out = _mito_ceiling_stage(task, pipeline, mito_stats, out)
    return out


def _mito_ceiling_stage(task, pipeline, mito_stats, out):
    """Derive the per-library mitochondrial ceiling alongside the count floors.

    Runs in the SAME step as the floors and off the SAME per-library pass, because the one thing
    that must be true of all three thresholds is that they describe the same population. Kept a
    separate function only so its refusals name the ceiling rather than the floors.
    """
    q = step_module("quality")
    samples = list(task.params["samples"])
    missing = [s for s in samples if s not in mito_stats]
    if missing:
        # Not defaulted and not skipped. A library with no mitochondrial summary is not a library
        # with no mitochondrial contamination, and filtering it on the cohort's other ceilings is
        # the borrowed-threshold failure this whole module exists to prevent.
        raise Refusal(
            f"05_quality (mitochondrial ceiling): no quartile summary for "
            f"{', '.join(missing)}. Either obs carries no pct_counts_mt - check the mt_prefix "
            f"matched anything - or the library has too few cells to place a quartile. Refusing "
            f"rather than deriving a ceiling for it from the other libraries.")

    assay = task.params.get("assay", "snrna")
    bounds = task.params.get("mito_bounds")
    declared_by = task.params.get("mito_bound_declared_by")
    try:
        d = q.derive_mito_ceiling_from_quartiles(
            mito_stats, assay=assay,
            bounds=tuple(bounds) if bounds else None, declared_by=declared_by)
    except Exception as e:                                                # noqa: BLE001
        raise Refusal(f"05_quality (mitochondrial ceiling): {e}") from None

    lo, hi = d["bounds"]
    ceil = d["ceilings"]
    print(f"    mitochondrial ceiling: per library "
          f"{min(m.ceiling for m in ceil.values()):.2f}-{max(m.ceiling for m in ceil.values()):.2f}% "
          f"(bound {lo}-{hi}%, {assay})")
    for n in d["notes"]:
        print(f"      {n}")

    import csv
    p = _tables(pipeline) / "mito_ceiling_per_sample.csv"
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "n", "median", "q1", "q3", "iqr", "derived", "ceiling", "clamped",
                    "iqr_mult", "assay", "bound_lo", "bound_hi"])
        for s in samples:
            m = ceil[s]
            w.writerow([s, m.n, f"{m.median:.6f}", f"{m.q1:.6f}", f"{m.q3:.6f}",
                        f"{m.iqr:.6f}", f"{m.derived:.6f}", f"{m.ceiling:.6f}", m.clamped,
                        d["mult"], assay, lo, hi])

    out = dict(out)
    out["outputs"] = list(out.get("outputs", [])) + [str(p)]
    out["metrics"] = dict(out.get("metrics", {}))
    out["metrics"].update({
        "mito_ceiling_lo": min(m.ceiling for m in ceil.values()),
        "mito_ceiling_hi": max(m.ceiling for m in ceil.values()),
        "mito_bound_binds": sum(1 for m in ceil.values() if m.clamped)})
    return out


# --------------------------------------------------------------------------------------------
# step 6 - cluster, profile, flag


def _cluster(task, pipeline, log):
    p = task.params
    return _scanpy(pipeline, "cluster",
                   pipeline.results / "objects" / f"{p['sample']}_ambient.h5",
                   pipeline.work / f"{p['sample']}_clusters",
                   {"resolution": p["resolution"], "seed": p["seed"]},
                   log, p["python_exe"])


def _cluster_flags(task, pipeline, log):
    import csv as _csv

    cf = step_module("cluster_flags")
    d = (task.params.get("decisions") or {}).get("cluster_check") or {}
    thr = cf.Thresholds(
        a_umi_frac=float(d.get("a_umi_fraction", 0.5)),
        b_pct_mt=float(d.get("b_mito_pct", 15.0)),
        c_uninformative=float(d.get("c_uninformative_pct", 50.0)),
        d_doublet=float(d.get("d_doublet_pct", 70.0)),
        source=("decisions.yml" if d else "PROPOSED from the cohort - not approved"))

    rows: list = []
    for f in sorted(pipeline.work.glob("*profile*.csv")):
        with open(f, encoding="utf-8", newline="") as fh:
            rows.extend(list(_csv.DictReader(fh)))
    if not rows:
        raise Refusal("06_cluster_check: no cluster profile was produced, so no cluster was "
                      "examined. That is not the same as no cluster being flagged.")

    flagged = cf.apply_flags(rows, thr)
    out = _tables(pipeline) / "cluster_profile.csv"
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=sorted({k for r in flagged for k in r}))
        w.writeheader()
        w.writerows(flagged)
    return {"outputs": [str(out)], "metrics": {"clusters": len(flagged)}, "versions": {}}


# --------------------------------------------------------------------------------------------
# step 7 - the only step that removes anything


def _apply(task, pipeline, log):
    import csv as _csv

    from .decisions import action_string, validate

    ap = step_module("apply")
    required = {
        "quality.umi_floor": "the UMI floor",
        "quality.gene_floor": "the gene floor",
        "quality.mito_ceiling_pct": "the mitochondrial ceiling",
    }
    resolved = validate(task.params["decisions"], required)
    apply_block = (task.params["decisions"].get("apply") or {})
    action = apply_block.get("action") or action_string(resolved)

    rows: list = []
    prof = _tables(pipeline) / "cluster_profile.csv"
    if prof.exists():
        with open(prof, encoding="utf-8", newline="") as fh:
            rows = list(_csv.DictReader(fh))
    for f in ap.preflight(rows, kept_total=0):
        pipeline.findings.append({"step": "07_apply", "check": f.check, "severity": f.severity,
                                  "message": f.message, "detail": []})

    raise Refusal(
        "07_apply is not wired to a written object.\n"
        f"    The decisions file is complete and authorises: {action}\n"
        "    The pre-flight ran and its findings are in the report. What does NOT exist yet is\n"
        "    the combined object this step would filter, so nothing is written.\n"
        "    Refusing rather than producing a deliverable that would not carry its own removal\n"
        "    ledger - an unrecoverable removal is the one outcome this pipeline exists to make\n"
        "    impossible.")


# --------------------------------------------------------------------------------------------
# What each step does, and what it cannot establish.
#
# The second half is required by docs/REPORT_DESIGN.md, and the report marks its absence as a
# defect rather than omitting the block - an omitted limit reads as no limit. Both halves live
# here, beside the code that produces the numbers, so a step and its caveat cannot drift apart.

STEP_TEXT = {
    "00_ingest": (
        "Validates the samplesheet and decides, per library, whether the supplied matrix is "
        "accepted or must be rebuilt from FASTQ.",
        "Whether the matrix is CORRECT. It establishes only that the values are raw counts and "
        "that empty droplets are still present - a matrix can be raw in both senses and still "
        "have been produced against the wrong reference."),
    "00_align": (
        "Rebuilds an unfiltered count matrix from FASTQ with the declared processor.",
        "Whether the chemistry declared for a library is the chemistry it was made with. A "
        "mismatch produces a near-empty matrix rather than an error, so a plausible cell count "
        "is not evidence the declaration was right."),
    "01_ambient": (
        "Denoises each library and audits the removal for degeneracy and for evenness across "
        "the design.",
        "Whether the denoising is CORRECT. Every check is cohort-relative or a differential; "
        "none can tell a well-removed ambient transcript from a badly-removed real one, because "
        "that requires knowing which cells should have expressed it."),
    "02_cells": (
        "Compares the aligner's cell call with the denoiser's and gates the loss.",
        "Whether the cells LOST were real. It measures how many and how unevenly, never what "
        "they were; only annotation could say that."),
    "04_doublets": (
        "Scores doublets per library above the light floor and checks the calls for health.",
        "Whether a called doublet IS one. It checks that the rate is a measurement rather than "
        "the prior, and that it does not fall unevenly across the design. Barcodes below the "
        "light floor were never scored and are reported as unknown, not as singlets."),
    "05_quality": (
        "Derives the count floors and the mitochondrial ceiling, in one pass over the same "
        "population. The floors are the density valley measured per library and proposed as ONE "
        "cohort constant; the ceiling is each library's OWN upper Tukey fence, Q3 + 1.5*IQR, "
        "bounded by a declared statement about what a nucleus can be.",
        "Whether either threshold is RIGHT. For the floors it establishes that the distribution "
        "has two modes and that the cut sits between them and inside plausible bounds - a tight "
        "real population can be refused by the dispersion test and a large enough artifact can "
        "pass it. For the ceiling it establishes only that each library was cut at its own "
        "outlier fence and that the result is even across the design; whether a "
        "mitochondria-high POPULATION is damage or a mitochondria-rich cell type needs an "
        "identity this pipeline does not establish."),
    "06_cluster_check": (
        "Clusters each library and flags clusters by depth, mitochondrial content, marker "
        "informativeness and doublet fraction.",
        "Whether a flagged cluster is TECHNICAL. A cluster flagged for mitochondrial content "
        "cannot be told from a mitochondria-rich cell type without an identity, which this "
        "pipeline does not establish."),
    "07_apply": (
        "Applies the recorded decisions and writes the deliverable with its removal ledger.",
        "Whether the removal was JUSTIFIED. It establishes that every removal was authorised in "
        "the operator's own words and that each removed observation remains recoverable with "
        "the criterion that removed it."),
    "report": (
        "Assembles this document from what the run recorded.",
        "Anything the run did not measure. A section reading NOT STATED is a gap in the "
        "pipeline, not a finding about the data."),
}


def step_text(step: str) -> tuple:
    """(what it does, what it cannot establish) for a step key, or two empty strings."""
    return STEP_TEXT.get(step, ("", ""))


__all__ = [
    "STEP_TEXT", "step_text", "_design", "_ingest", "_align", "_ambient", "_ambient_supplied",
    "_ambient_audit", "_cellcall",
    "_doublets", "_doublet_health", "_quality", "_quality_stage", "_cluster",
    "_cluster_flags", "_apply", "_report",
]
