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

import json
from pathlib import Path

from .pipeline import step_module
from .task import Refusal, Task, TaskFailure


# --------------------------------------------------------------------------------------------
# helpers


#: How many (rank, count) pairs step 0 keeps per library for figure F1. FIXED, not a threshold:
#: the pairs are downsampled log-uniformly and the curve is read on log-log axes, so this changes
#: how heavy the file is and nothing about where the knee falls. No decision reads it.
RANK_POINTS = 2000


def _tables(p) -> Path:
    return p.results / "tables"


def _objects(p) -> Path:
    return p.results / "objects"


def _promote(pipeline, result: dict, suffix: str, dest_name: str, *, what: str):
    """Copy one output an op wrote under `scratch` into `results/tables/`, and return where.

    WHY A COPY AND NOT A DIFFERENT out_prefix. The scanpy ops write four or five files at one
    prefix and only some of them are evidence. `<sample>.valleys.json` and
    `<sample>.called_barcodes.csv` are intermediates this pipeline reads and nothing else opens;
    the density curve and the embedding are the data two figures are drawn FROM, and rule 1 of
    `report/collect.py` is that a figure comes from a named file a reader can open. Pointing the
    whole prefix at `tables/` would publish the intermediates too and invite them to be quoted.

    WHY IT EXISTS AT ALL. The density curve was already being written, correctly, on every run -
    and then left in a scratch directory the report never looks in. The report duly printed "step
    5 fits a KDE, takes the minimum, records the valley position, and discards the curve" about a
    curve sitting on disk. A figure's data being PRODUCED is not the same as its being KEPT, and
    the gap between the two is invisible from either end.

    Returns None when the op declared no such output - which is what a null declaration produces
    and is not an error here. The caller reports the absence.
    """
    import shutil

    src = next((Path(o) for o in (result.get("outputs") or []) if str(o).endswith(suffix)), None)
    if src is None:
        return None
    if not src.exists():
        raise TaskFailure(
            f"{what}: the op declared {src} and it is not on disk, so it cannot be published to "
            f"tables/. A declared output that is absent is a broken record of the run, not a "
            f"missing figure.")
    dst = _tables(pipeline) / dest_name
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst


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

    measured: dict = {}

    def stats_fn(matrix_path):
        res = mx.run_summary_stats(
            matrix=matrix_path,
            out_json=pipeline.work / f"{sample}_ingest_stats.json",
            log=pipeline.work / f"{sample}_ingest_stats.log",
            python_exe=python_exe,
            expected_genes=task.params.get("expected_genes"),
            tmp_dir=pipeline.scratch / f"{sample}_extract",
            # THE CURVE BEHIND THE VERDICT, measured in the pass that is already reading this
            # matrix. Step 0 opens every raw matrix to establish that the empty droplets are
            # still there; F1 is the picture of that same fact, and the only figure in the report
            # a reader can check the claim "this input is raw" against. Asking for it here costs
            # one more pass over row totals in a subprocess that has the matrix open.
            #
            # RANK_POINTS IS A FIXED PROCEDURE PARAMETER, not a threshold. The curve is read on
            # log-log axes; the pairs are downsampled log-uniformly, so 2,000 of them draw the
            # knee at the same place 200,000 would and no decision anywhere depends on the number.
            rank_points=RANK_POINTS,
            executor=pipeline.executor)
        measured["barcode_rank"] = res.get("barcode_rank")
        measured["n_barcodes"] = res["metrics"].get("n_barcodes")
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

    # Step 0 decides, and now also KEEPS what it measured while deciding. A sample rebuilt from
    # FASTQ has no supplied matrix and no curve, which is why this is conditional rather than
    # promised: the file is absent for that library and F1 draws the libraries it has, with the
    # n stated on the figure.
    # THE VERDICT ITSELF, in the shape the report's provenance block reads.
    #
    # Step 0's whole job is to establish P1 (raw VALUES) and P2 (raw DROPLETS), it did establish
    # them for every library, and the verdict went to the console and nowhere else - so the report
    # said "the raw-input verification was not recorded. A pre-filtered matrix cannot be
    # un-filtered and nothing downstream detects it." That was true of the RECORD and false of the
    # run, which is the worst pairing: the check happened and the document could not say so.
    #
    # Recorded even when it FAILED, and even when there was no matrix to check - `plan.verdict` is
    # None for a library rebuilt from FASTQ, and that is written as a name with no determination
    # rather than omitted. The report is three-valued here: NOT CHECKED is not a PASS.
    v = plan.verdict
    check = {"name": sample,
             "p1_raw_values": None if v is None else bool(v.p1_raw_values),
             "p2_raw_droplets": None if v is None else bool(v.p2_raw_droplets),
             "reasons": list(getattr(v, "reasons", []) or [])}

    outs: list = []
    curve = measured.get("barcode_rank")
    if curve:
        import csv as _csv

        out = _tables(pipeline) / f"{sample}.barcode_rank.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as fh:
            w = _csv.writer(fh)
            w.writerow(["sample", "rank", "total_counts", "n_barcodes"])
            for r, t in curve:
                w.writerow([sample, r, t, measured.get("n_barcodes")])
        outs.append(str(out))
    return {"outputs": outs,
            "metrics": {"mode": plan.mode, "processor": plan.processor,
                        "reason": plan.reason,
                        "n_rank_points": len(curve) if curve else None,
                        "input_check": check},
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
        m = cbd.parse_metrics(r["h5"], r["raw"],
                              cell_barcodes_csv=(r.get("cell_barcodes") or None))
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

    # UNMEASURABLE AND UNMEASURED MUST NOT READ THE SAME. A library with no raw counts genuinely
    # cannot have its removal fraction measured and is named here; one that was merely SUPPLIED is
    # audited like any other, because the raw counts the samplesheet names are what the fraction
    # is measured against and it makes no difference who ran the denoiser.
    no_raw = list(task.params.get("no_raw") or [])
    if no_raw:
        print(f"    NOT AUDITED: {len(no_raw)} library(ies) have no raw matrix "
              f"({', '.join(sorted(no_raw))}). The fraction removed cannot be measured without "
              f"the counts the denoiser started from.")
    if supplied:
        print(f"    {len(supplied)} library(ies) arrived already corrected and ARE audited "
              f"against their raw counts; provenance in tables/ambient_supplied.json.")
    if not rows:
        pipeline.gate("01_ambient", [], "NOT RUN")
        return {"outputs": [str(out)],
                "metrics": {"libraries": 0, "no_raw": len(no_raw),
                            "supplied": len(supplied)},
                "versions": {}}

    findings = aa.audit(rows, per_gene, task.params["design"])

    # --- IS ANY FIT UNLIKE ITS SIBLINGS? modules/01_ambient/lr_policy.py has answered this since
    # it was written and nothing had ever called it - neither it nor the two adapter entry points
    # that exist to feed it. A degenerate CellBender fit has no per-sample symptom: on its own it
    # produces a matrix, a report and a plausible cell count, and is only legible beside the other
    # libraries. An unwired cohort check is therefore not a partial safeguard, it is none.
    lr_findings, lr_outliers = _ambient_lr_findings(task, pipeline, rows, aa)
    findings += lr_findings

    pipeline.gate("01_ambient", findings, aa.verdict(findings))
    return {"outputs": [str(out)],
            "metrics": {"libraries": len(rows), "no_raw": len(no_raw),
                        "supplied": len(supplied),
                        # Carried on the result so step 2 can read it from the manifest and judge
                        # an outlying fit together with an unusual cell-call loss.
                        "lr_outliers": lr_outliers},
            "versions": {}}


def _ambient_lr_findings(task, pipeline, rows, aa) -> list:
    """Cohort-relative learning-rate assessment, as findings. Re-runs nothing.

    Detection and reporting only: where a library is unlike its siblings this says so and names
    it, and a human decides whether to re-run the denoiser at half the rate. `lr_policy` describes
    what that re-run must then prove - halve, RE-MEASURE, adopt only if the diagnostic resolves,
    and change the rate for every library rather than one - and none of that is done here.

    A library whose learning curve was not declared is reported as NOT MEASURED and named. It is
    not dropped from the cohort quietly: `assess_cohort` is a comparison between siblings, so a
    missing sibling changes what the remaining ones are being compared against.
    """
    from adapters import cellbender as cbd

    lp = step_module("lr_policy")
    per = task.params["per_sample"]
    diag, missing = {}, []
    for r in rows:
        s = r["sample"]
        spec = per.get(s) or {}
        source = spec.get("log") or spec.get("metrics_csv") or spec.get("h5")
        curve = {}
        if source:
            try:
                curve = cbd.parse_learning_curve(
                    source, metrics_csv=spec.get("metrics_csv") or None) or {}
            except Exception as e:                                        # noqa: BLE001
                curve = {"error": f"{type(e).__name__}: {e}"}
        ci = curve.get("convergence_indicator")
        if ci is None:
            missing.append(s)
            continue
        diag[s] = {"fraction_removed": r.get("fraction_removed_overall"),
                   "convergence_indicator": ci}

    import csv as _csv
    p = _tables(pipeline) / "ambient_lr_diagnostics.csv"
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh)
        w.writerow(["sample", "fraction_removed", "convergence_indicator", "measured"])
        for r in rows:
            s = r["sample"]
            d = diag.get(s) or {}
            w.writerow([s, d.get("fraction_removed", ""), d.get("convergence_indicator", ""),
                        s in diag])

    out = []
    if missing:
        out.append(aa.Finding(
            "learning rate: diagnostics not measured", "REVIEW",
            f"{len(missing)} of {len(rows)} library(ies) carry no learning curve "
            f"({', '.join(sorted(missing))}), so they took no part in the comparison. Declare "
            f"`ambient_log` or `ambient_metrics` for them. A cohort check is a comparison between "
            f"siblings, and a missing sibling changes what the rest were compared against.",
            ["looked for a declared log, then a declared metrics file, then the object itself"]))
    if len(diag) < 4:
        out.append(aa.Finding(
            "learning rate: not assessed", "REVIEW",
            f"only {len(diag)} library(ies) had comparable diagnostics, and this check does not "
            f"answer below about four. It is cohort-relative by construction: one library's "
            f"numbers can say that the optimiser complained, never that the fit is degenerate.",
            []))
        return out, []

    v = lp.assess_cohort(diag, label="as delivered")
    print(f"    learning rate: {v.action}" + (f" - {v.note}" if v.note else ""))
    for s in sorted(v.outliers or {}):
        for why in (v.outliers[s] or []):
            print(f"      OUTLIER {s}: {why}")

    if v.action == "keep_default":
        out.append(aa.Finding(
            "learning rate: is any fit unlike its siblings", "ok",
            f"no library is an outlier on either diagnostic across {len(diag)} libraries. "
            + (v.note or ""), []))
        return out, []

    # REFUSED, not noted. Continuing would deliver a filtered object built on a fit this pipeline
    # has just judged degenerate, and a degenerate fit is invisible in its own run - which is the
    # whole reason the check is cohort-relative. Nothing is re-run automatically: what the re-run
    # would have to prove is stated, and a human decides.
    out.append(aa.Finding(
        "learning rate: a fit is unlike its siblings", "REFUSE",
        f"{len(v.outliers)} of {len(diag)} libraries sit outside the cohort "
        f"({', '.join(sorted(v.outliers))}); the policy's action is {v.action!r}. "
        + (v.note + " " if v.note else "")
        + "Nothing was re-run - this pipeline detects and reports, and the decision is yours. A "
          "half-rate re-run has to HALVE, RE-MEASURE and resolve the diagnostic to count: moving "
          "a complaint from one epoch to another is not a resolution, and the rate changes for "
          "every library rather than the flagged one, or the denoising becomes a technical "
          "property varying across the design.",
        [f"{s}: {'; '.join(why)}" for s, why in sorted((v.outliers or {}).items())]))
    return out, sorted(v.outliers or {})


# --------------------------------------------------------------------------------------------
# step 2 - the cell-call gate


def _cellcall(task, pipeline, log):
    """Compare the aligner's cell call with the denoiser's, over the barcodes themselves.

    THE DENOISER'S CALL COMES FROM THE OBJECT EVERY OTHER STEP READS.

    It used to come from a `_cell_barcodes.csv` named in the samplesheet - a second artefact,
    from a run nothing verified was this one. Two denoising runs of the same library produce
    different calls in identically-shaped files, so the comparison could be against the wrong
    denoising and nothing downstream could tell.

    The failure is not hypothetical. On the calibration cohort the object's barcodes carried a
    `<sample>_` prefix and the CSV's did not: the same run, described two ways, intersecting in
    ZERO barcodes. A check that had simply intersected them would have reported every cell lost.

    So the denoiser's call is read from `<sample>_ambient.h5` - written by step 1, read by steps
    4, 5 and 6 - and the aligner's call remains an external file, because the aligner's output
    genuinely is one. A supplied `cellbender_barcodes` is no longer the source; it is CHECKED
    against the object, after normalising the sample prefix, and a disagreement stops the run.
    """
    from adapters import matrix as mx

    cg = step_module("cellcall_gate")
    paths = task.params.get("call_paths") or {}
    calls, missing, notes = {}, [], []
    for s in task.params["samples"]:
        pth = paths.get(s) or {}
        a_path = pth.get("aligner")
        if not a_path:
            missing.append(f"{s}: aligner_cells")
            continue

        # Read from the measurement task's own result, so step 2 and step 5 are looking at one
        # pass over one object rather than two reads that could disagree.
        r = pipeline.results_by_key.get(f"05_quality/{s}")
        called_csv = (getattr(r, "metrics", None) or {}).get("called_barcodes") if r else None
        if not called_csv or not Path(called_csv).exists():
            missing.append(f"{s}: no denoiser call was measured from {s}_ambient.h5")
            continue
        called = mx.called_barcodes(called_csv)

        a = {_bare(b, s) for b in mx.called_barcodes(a_path)}
        c = {_bare(b, s) for b in called}

        # A supplied CSV is a cross-check now, never the source.
        csv_path = pth.get("cellbender")
        if csv_path and Path(csv_path).exists():
            declared = {_bare(b, s) for b in mx.called_barcodes(csv_path)}
            if declared != c:
                only_csv, only_obj = declared - c, c - declared
                raise Refusal(
                    f"02_cells ({s}): the declared cellbender_barcodes disagrees with the object "
                    f"step 1 produced. {len(only_csv):,} barcode(s) are in the CSV and not the "
                    f"object; {len(only_obj):,} are in the object and not the CSV. They are not "
                    f"the same denoising run, and comparing the aligner against the wrong one "
                    f"produces a loss figure that describes nothing.\n"
                    f"    object: {len(c):,} called   CSV: {len(declared):,} called\n"
                    f"    Remove cellbender_barcodes, or point ambient_h5 at the run it came from.")
            notes.append(s)

        calls[s] = {"aligner": len(a), "cellbender": len(c), "lost": len(a - c)}
        print(f"    {s:<14} aligner {len(a):>7,}   denoiser {len(c):>7,}   "
              f"lost {len(a - c):>6,}  ({100 * len(a - c) / max(len(a), 1):.2f}%)")

    if missing:
        raise Refusal(
            "02_cells cannot compare cell calls:\n"
            + "\n".join(f"    - {m}" for m in missing)
            + "\n    `aligner_cells` is the aligner's filtered matrix directory (CeleScope "
              "outs/filtered, CellRanger filtered_feature_bc_matrix). The denoiser's call is "
              "read from <sample>_ambient.h5 and needs no declaration.")
    if notes:
        print(f"    cross-checked against a declared cellbender_barcodes for "
              f"{len(notes)} library(ies); all agree with the object")

    findings = cg.gate(calls, task.params["design"])

    # THE TWO OBSERVATIONS ABOUT ONE DENOISING, READ TOGETHER.
    #
    # Step 1 can say that a library's fit is unlike its siblings; step 2 can say that a library
    # loses an unusual share of the aligner's cells. Separately each is weak - an outlying
    # diagnostic may be harmless, and a loss may be ordinary biology - and nothing had ever put
    # them side by side, so a library that was BOTH read as two unremarkable entries in different
    # sections of the report. A degenerate fit that also discards cells is a different claim from
    # either half.
    lr = (getattr(pipeline.results_by_key.get("01_ambient_audit"), "metrics", None) or {})
    flagged = set(lr.get("lr_outliers") or [])
    if flagged and calls:
        worst = {s: 100 * c["lost"] / max(c["aligner"], 1) for s, c in calls.items()}
        both = sorted(s for s in flagged if worst.get(s, 0) > 0)
        findings.append(cg.GateFinding(
            "denoising: flagged fit and cell-call loss in the same library",
            "REFUSE" if both else "REVIEW",
            (f"{', '.join(both)} were flagged by the learning-rate check AND lose cells the "
             f"aligner called ("
             + "; ".join(f"{s} {worst[s]:.2f}%" for s in both)
             + "). Each on its own is weak evidence; together they describe a fit that is both "
               "unlike its siblings and discarding data."
             if both else
             f"{len(flagged)} library(ies) were flagged by the learning-rate check "
             f"({', '.join(sorted(flagged))}) and none of them loses any cell the aligner "
             f"called. The flag stands; it is not corroborated here."),
            sorted(f"{s}: lost {calls[s]['lost']:,} of {calls[s]['aligner']:,}"
                   for s in flagged if s in calls)))

    pipeline.gate("02_cells", findings, cg.verdict(findings))
    out = _tables(pipeline) / "cell_calls.csv"
    import csv
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "aligner", "denoiser", "lost"])
        for s, c in sorted(calls.items()):
            w.writerow([s, c["aligner"], c["cellbender"], c["lost"]])
    return {"outputs": [str(out)], "metrics": {"libraries": len(calls)}, "versions": {}}


def _bare(barcode, sample) -> str:
    """A barcode without the `<sample>_` prefix some conversions add.

    The same cell is written `AACGT...` by one tool and `lib3_AACGT...` by another, and the two
    sets then intersect in nothing while describing identical data. Normalising here means the
    comparison is between cells rather than between naming conventions - and stripping only an
    exact `<sample>_` prefix cannot merge two distinct barcodes, which a general strip could.
    """
    b = str(barcode)
    pre = f"{sample}_"
    return b[len(pre):] if b.startswith(pre) else b


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


#: Every threshold and every count this pipeline derives, one column each, with the SCOPE that
#: says who it was derived for. A cohort constant repeated down the column is not clutter: a
#: reader comparing two libraries has to be able to see, without going anywhere else, which
#: numbers differ because the libraries differ and which are the same by construction.
#:
#: `key` is looked up in the row dict built below; `scope` is one of:
#:    per library      measured or derived for that library alone
#:    cohort constant  one value proposed for every library
PER_SAMPLE_COLUMNS = (
    ("cells_aligner", "aligner cells", "per library", "02_cells"),
    ("cells_denoiser", "denoiser cells", "per library", "05_quality"),
    ("cells_lost", "cells lost at the call", "per library", "02_cells"),
    ("light_floor_umi", "light floor (UMI)", "cohort constant", "04_doublets"),
    ("doublets_scored", "scored for doublets", "per library", "04_doublets"),
    ("doublets_called", "called doublet", "per library", "04_doublets"),
    ("doublet_rate_pct", "doublet rate %", "per library", "04_doublets"),
    ("umi_valley", "UMI valley", "per library", "05_quality"),
    ("umi_bimodal", "UMI bimodal", "per library", "05_quality"),
    ("umi_floor_proposed", "UMI floor proposed", "cohort constant", "05_quality"),
    ("gene_valley", "gene valley", "per library", "05_quality"),
    ("gene_bimodal", "gene bimodal", "per library", "05_quality"),
    ("gene_floor_proposed", "gene floor proposed", "cohort constant", "05_quality"),
    ("mito_pop_floor", "mito derived above (UMI)", "cohort constant", "05_quality"),
    ("mito_pop_n", "mito derived over (n)", "per library", "05_quality"),
    ("mito_q1", "mito Q1 %", "per library", "05_quality"),
    ("mito_q3", "mito Q3 %", "per library", "05_quality"),
    ("mito_ceiling_pct", "MITO CEILING %", "per library", "05_quality"),
    ("mito_clamped", "ceiling clamped", "per library", "05_quality"),
    ("clusters", "clusters", "per library", "06_cluster_check"),
    ("clusters_flagged", "clusters flagged", "per library", "06_cluster_check"),
)


def _per_sample_thresholds(pipeline, samples: list) -> tuple:
    """One row per library, one column per threshold, written to a file and returned.

    Assembled from the RUN MANIFEST and from the per-sample tables the steps wrote, never
    recomputed here: a report that derives its own numbers can disagree with the run it describes
    and nothing would say which was right.

    A value that no step recorded stays None and is rendered as such. The point of a per-library
    table is to show which thresholds vary by library and which do not, and a blank that reads as
    a zero destroys exactly that.
    """
    import csv as _csv

    def metrics(key):
        r = pipeline.results_by_key.get(key)
        return (getattr(r, "metrics", None) or {}) if r is not None else {}

    def table(path, key="sample"):
        p = Path(path)
        if not p.exists():
            return {}
        with p.open(encoding="utf-8", newline="") as fh:
            return {r[key]: r for r in _csv.DictReader(fh) if r.get(key)}

    cohort = metrics("05_quality")
    calls = table(_tables(pipeline) / "cell_calls.csv")
    ceilings = table(_tables(pipeline) / "mito_ceiling_per_sample.csv")

    # A library that WAS profiled and had nothing flagged is a zero, not a blank. Counting only
    # the hits left the two libraries with no flagged cluster indistinguishable from the ones
    # step 6 never reached - the unknown-is-not-a-zero rule, running in the direction people
    # forget: a measured zero must not read as unknown either.
    flagged: dict = {}
    prof = _tables(pipeline) / "cluster_profile.csv"
    if prof.exists():
        with prof.open(encoding="utf-8", newline="") as fh:
            for r in _csv.DictReader(fh):
                s = r.get("sample", "")
                flagged.setdefault(s, 0)
                if str(r.get("FLAG", "")).strip().lower() == "true":
                    flagged[s] += 1

    rows = []
    for s in samples:
        q5, d4 = metrics(f"05_quality/{s}"), metrics(f"04_doublets/{s}")
        val, bim = q5.get("valleys") or {}, q5.get("bimodal") or {}
        pop, cc, ce = (q5.get("mito_population") or {}), calls.get(s, {}), ceilings.get(s, {})
        rate = d4.get("rate_over_scored")
        rows.append({
            "sample": s,
            "cells_aligner": _int(cc.get("aligner")),
            "cells_denoiser": q5.get("n_called_by_denoiser"),
            "cells_lost": _int(cc.get("lost")),
            "light_floor_umi": pop.get("floor_umi"),
            "doublets_scored": d4.get("n_scored"),
            "doublets_called": d4.get("n_called"),
            "doublet_rate_pct": round(100 * rate, 2) if rate is not None else None,
            "umi_valley": _round(val.get("umi"), 1),
            "umi_bimodal": bim.get("umi"),
            "umi_floor_proposed": cohort.get("umi_proposed"),
            "gene_valley": _round(val.get("genes"), 1),
            "gene_bimodal": bim.get("genes"),
            "gene_floor_proposed": cohort.get("genes_proposed"),
            "mito_pop_floor": pop.get("floor_umi"),
            "mito_pop_n": pop.get("n_at_or_above"),
            "mito_q1": _round(ce.get("q1"), 3),
            "mito_q3": _round(ce.get("q3"), 3),
            "mito_ceiling_pct": _round(ce.get("ceiling"), 3),
            "mito_clamped": ce.get("clamped") or None,
            "clusters": metrics(f"06_cluster/{s}").get("n_clusters"),
            "clusters_flagged": flagged.get(s),
        })

    out = _tables(pipeline) / "thresholds_per_sample.csv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = _csv.writer(fh)
        w.writerow(["sample"] + [k for k, _l, _s, _st in PER_SAMPLE_COLUMNS])
        w.writerow(["scope"] + [sc for _k, _l, sc, _st in PER_SAMPLE_COLUMNS])
        for r in rows:
            w.writerow([r["sample"]] + ["" if r.get(k) is None else r.get(k)
                                        for k, _l, _s, _st in PER_SAMPLE_COLUMNS])
    block = {"source": str(out),
             "columns": [{"key": k, "label": lab, "scope": sc, "step": st}
                         for k, lab, sc, st in PER_SAMPLE_COLUMNS],
             "rows": rows}
    return block, out


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _round(v, n):
    try:
        return round(float(v), n)
    except (TypeError, ValueError):
        return None


def _report(task, pipeline, log):
    from report.build import build_report

    # report_payload() assembles the per-library table and writes its CSV, because scqc_cli
    # rebuilds this payload after the graph finishes and would otherwise write over anything
    # added here. This step reports the file; it does not produce it.
    payload = pipeline.report_payload(stopped=None)

    # The figures come from report_payload(), like every other section - NOT assembled here.
    # Assembling them here worked and was still wrong: finish() rebuilds the payload after the
    # graph and wrote a figure-less document over the top of this one.
    figures = payload.get("figures") or {}
    figure_notes = payload.get("figure_notes") or {}
    print(f"    figures: {len(figures)} assembled ({', '.join(sorted(figures)) or 'none'})")
    for fid in sorted(figure_notes):
        print(f"      {fid} not produced: {figure_notes[fid]}")

    per_sample = payload.get("per_sample") or {}
    out_html = pipeline.results / "reports" / "qc_report.html"
    out_json = pipeline.results / "reports" / "report.json"

    # The PAYLOAD, beside the document built from it. report.json is the rendered document and
    # cannot be fed back in; without the payload, changing a caption or a figure meant running
    # the whole pipeline again to redraw it. With it, `scqc report <results>` rebuilds from the
    # run's own files in seconds and cannot reach the matrices, so no number can change.
    out_payload = pipeline.results / "reports" / "payload.json"
    out_payload.parent.mkdir(parents=True, exist_ok=True)
    out_payload.write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")

    build_report(payload, out_html, out_json)
    outputs = [str(out_html), str(out_json), str(out_payload)]
    if per_sample.get("source"):
        outputs.append(str(per_sample["source"]))
    return {"outputs": outputs,
            "metrics": {"findings": len(pipeline.findings),
                        "per_sample_rows": len(per_sample.get("rows") or []),
                        "per_sample_columns": len(per_sample.get("columns") or [])},
            "versions": {}}


# --------------------------------------------------------------------------------------------
# graph construction


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
    import csv as _csv

    from adapters import doublets as db

    p = task.params
    # Local scratch: this triple is written by Python and read straight back by R, and nothing
    # downstream opens it. On NFS it was the whole cost of step 4.
    mtx = pipeline.scratch / f"{p['sample']}_dbl_mtx"
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

    # WHAT THE LIGHT FLOOR DID, written down where it was applied.
    #
    # The floor is applied HERE, in the export, and its consequence lived only in this local
    # `export` object: `<sample>_doublets.csv` holds the SCORED barcodes and nothing else, so the
    # never-examined population existed nowhere on disk and step 3 had no record in the report at
    # all. The two reasons stay apart, as `ExportedMatrix` keeps them - below the floor is the
    # threshold doing its documented job; not selected is the caller handing in a smaller
    # population - because only the first is explained by a threshold.
    floor_csv = _tables(pipeline) / f"{p['sample']}.light_floor.csv"
    floor_csv.parent.mkdir(parents=True, exist_ok=True)
    with floor_csv.open("w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh)
        w.writerow(["sample", "light_floor_umi", "n_exported", "n_below_floor",
                    "n_not_selected", "n_never_examined"])
        w.writerow([p["sample"], int(p["light_floor"]), len(export.exported),
                    len(export.below_floor), len(export.not_selected), len(unscored)])

    res = db.run_scdblfinder(
        rscript=p["rscript"], mtx_dir=mtx, out_csv=out_csv,
        dbr=p.get("dbr"), dbr_sd=p.get("dbr_sd"), seed=int(p.get("seed", 0)),
        log=log, executor=pipeline.executor, sample=p["sample"],
        unscored=unscored)
    res = {**res, "outputs": list(res.get("outputs") or []) + [str(floor_csv)]}
    res["metrics"] = {**(res.get("metrics") or {}),
                      "light_floor_umi": int(p["light_floor"]),
                      "n_exported": len(export.exported),
                      "n_below_floor": len(export.below_floor),
                      "n_not_selected": len(export.not_selected)}
    return res


def _light_floor(task, pipeline, log):
    """Step 3's record: what the light floor left out, per library. Removes nothing.

    WHY THIS IS A TASK AND NOT A PARAMETER. `03_light_floor` is one of the report's eight steps
    and had no task at all, so every run said "no record of this step at all. Whether it ran is
    unknown, and unknown is not the same as did-not-run". That was accurate: the floor is applied
    inside step 4's export and nothing reported on it. A step in the spine with nothing filling it
    is the report's own defect, and adding STEP_TEXT does not fix it - `step_text()` is only
    consulted for steps that have tasks.

    IT REPORTS, IT DOES NOT DECIDE. The floor is DECLARED (`--light-floor`), applied in step 4,
    and this reads what that application recorded. It is deliberately not a place where the floor
    could be changed: a threshold that can be set in two places is a threshold nobody can trace.
    """
    import csv as _csv

    lf = step_module("light_floor")
    samples = list(task.params["samples"])
    rows, absent = [], []
    for s in samples:
        r = pipeline.results_by_key.get(f"04_doublets/{s}")
        m = (getattr(r, "metrics", None) or {}) if r is not None else {}
        if m.get("n_exported") is None:
            absent.append(s)
            continue
        rows.append({"sample": s, "light_floor_umi": m.get("light_floor_umi"),
                     "n_exported": m.get("n_exported"),
                     "n_below_floor": m.get("n_below_floor"),
                     "n_not_selected": m.get("n_not_selected")})
    if absent:
        raise Refusal(
            f"03_light_floor: step 4 recorded no export for {', '.join(absent)}, so what the "
            f"floor left out in that library is unknown. Reporting the cohort from the rest "
            f"would state a coverage this run did not measure.")

    out = _tables(pipeline) / "light_floor.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=["sample", "light_floor_umi", "n_exported",
                                            "n_below_floor", "n_not_selected"])
        w.writeheader()
        w.writerows(rows)

    floors = set()
    for r in rows:
        floors.add(r["light_floor_umi"])
    total_below = sum(int(r["n_below_floor"] or 0) for r in rows)
    total_other = sum(int(r["n_not_selected"] or 0) for r in rows)
    for r in rows:
        print(f"    {r['sample']:<14} exported {int(r['n_exported'] or 0):>7,}   "
              f"below floor {int(r['n_below_floor'] or 0):>7,}   "
              f"not selected {int(r['n_not_selected'] or 0):>6,}")
    print(f"    floor {sorted(floors)} - DECLARED; {total_below:,} barcode(s) below it were "
          f"never examined for doublets, which is UNKNOWN and not a singlet")

    # The module's own statement of what the floor is for, so the number and its meaning cannot
    # drift apart in the report.
    note = getattr(lf, "PURPOSE", None)
    return {"outputs": [str(out)],
            "metrics": {"light_floor_umi": sorted(floors)[0] if len(floors) == 1 else None,
                        "libraries": len(rows),
                        "n_below_floor": total_below,
                        "n_not_selected": total_other,
                        "class": "DECLARED - --light-floor; applied in step 4's export",
                        **({"purpose": note} if note else {})},
            "versions": {}}


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

    # THE RATES THIS GATE JUDGED, written down. They were returned as metrics and the task
    # declared no output, so the report said every one of them "states a value with no source
    # file" - correctly: the numbers the health check refuses or passes on were nowhere a reader
    # could open them. A gate whose evidence is not on disk cannot be re-checked.
    import csv as _csv

    out = _tables(pipeline) / "doublet_health.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh)
        w.writerow(["sample", "n_scored", "n_called_doublet", "rate_over_scored"])
        for s in sorted(rates):
            calls = db.read_calls(_tables(pipeline) / f"{s}_doublets.csv", sample=s)
            n_pos = sum(1 for _score, is_doublet in calls.values() if is_doublet)
            w.writerow([s, len(calls), n_pos, rates[s]])
    return {"outputs": [str(out)],
            "metrics": {"libraries": len(known), "never_scored": unscored_total},
            "versions": {}}


def _doublet_sweep(task, pipeline, log):
    """Re-score ONE library at every swept dbr.sd. Writes no calls anyone downstream reads.

    THIS DOES NOT CHANGE THE DELIVERABLE. The applied calls are the ones `_doublets` wrote at the
    declared setting; these go to scratch and exist so figure F5 can show whether the rate the
    pipeline applied moved with the data or with the prior. A flat rate across libraries that
    differ 2.5-fold in size is the prior's flatness, and one point cannot reveal it.

    PER LIBRARY, NOT PER COHORT, AND THE ADAPTER IS STILL CALLED ONCE PER SAMPLE. `db.sweep()`
    loops over samples internally; handing it all ten would run thirty scDblFinder invocations in
    a single task, serially, on a machine with a scheduler outside it. One sample at a time keeps
    every cross-setting check the adapter makes - same seed, one version of the tool, no label
    reused - and lets the ten run at once.
    """
    from adapters import doublets as db

    p = task.params
    s = p["sample"]
    # Its OWN export, not the one `_doublets` left in scratch. Depending on another task's
    # intermediate makes this step's correctness depend on when scratch is cleaned, and the
    # export is cheap beside thirty seconds of xgboost.
    mtx = pipeline.scratch / f"{s}_dbl_sweep_mtx"
    export = db.export_matrix(p["h5"], mtx, min_umi=int(p["light_floor"]))
    res = db.sweep(list(p["settings"]), p["rscript"], {s: mtx},
                   pipeline.scratch / f"{s}_dbl_sweep", p["dbr"], int(p.get("seed", 0)),
                   pipeline.work, pipeline.executor,
                   unscored={s: export.unscored},
                   # The deepest-decile alarm needs the per-barcode depth, and the export has just
                   # measured it. Computed per library here; whether the COHORT figure can be
                   # assembled from ten of them is settled at the barrier, which refuses to.
                   umi_by_barcode={s: export.umi_by_barcode})
    per_setting = (res.get("metrics") or {}).get("per_setting") or {}

    # THE SWEPT RATES AS A FILE, not as a nested blob in the manifest.
    #
    # This returned `per_setting` - a mapping of mappings - as a metric, and declared no output.
    # The report turned each metric into a finding whose source is the task's first output, so
    # every library produced two findings reading "states a value with no source file", and the
    # "value" was a dict nobody could read anyway. One table per library, exactly the pairing
    # `<sample>.percell.csv` has with `removal_ledger.csv`.
    import csv as _csv

    out = _tables(pipeline) / f"{s}.doublet_sweep.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    scalars: dict = {}
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh)
        w.writerow(["sample", "setting", "dbr_sd_value", "n_scored", "n_called",
                    "rate_over_scored"])
        for label in sorted(per_setting):
            e = per_setting[label]
            rate = (e.get("per_sample_rate_over_scored") or {}).get(s)
            n_scored = (e.get("per_sample_n_scored") or {}).get(s)
            n_called = (e.get("per_sample_n_called") or {}).get(s)
            w.writerow([s, label, e.get("dbr_sd_value"), n_scored, n_called, rate])
            print(f"    {s:<14} dbr.sd={label:<8} "
                  + (f"{100 * rate:.2f}% of scored" if rate is not None else "NO RATE"))
            # Scalars only. A metric is rendered as one line of a findings table, and a nested
            # mapping there is a value no reader reads and no reviewer can check.
            if rate is not None:
                scalars[f"rate_at_dbr_sd_{label}"] = rate
            if n_scored is not None:
                scalars["n_scored"] = n_scored
    return {"outputs": [str(out)],
            "metrics": {"settings": sorted(per_setting), **scalars},
            "versions": res.get("versions", {})}


def _doublet_sweep_stage(task, pipeline, log):
    """The barrier: one row per (library, setting) for F5, and what the sweep recommends.

    THE DEEP-DECILE ARM IS NOT EVALUATED HERE, AND THAT IS RECORDED RATHER THAN PAPERED OVER.

    `doublet.recommend()` rejects the fully-free setting when it calls half of the deepest UMI
    decile - a cohort quantity, over every library's barcodes pooled. Ten per-library deciles are
    not that number and no combination of them is: the worst of them is a different statistic, and
    the mean of them is a different statistic again. Both would print. So `deep_decile_rate` is
    left None, which `recommend()` is documented to read as NO EVIDENCE rather than as
    reassurance, and the table says which arm was not evaluated.
    """
    import csv as _csv

    d = step_module("doublet")
    samples = list(task.params["samples"])
    merged: dict = {}
    absent = []
    # FROM THE FILES THE TASKS DECLARED, taken out of the manifest rather than found by a glob -
    # the rule `_cluster_flags` states and for the same two reasons: a glob cannot tell this run's
    # output from an earlier one's, and it silently accepts nine files where ten were swept.
    for s in samples:
        r = pipeline.results_by_key.get(f"04_doublet_sweep/{s}")
        paths = [o for o in (getattr(r, "outputs", None) or [])
                 if str(o).endswith(".doublet_sweep.csv")]
        if not paths:
            absent.append(s)
            continue
        with open(paths[0], encoding="utf-8", newline="") as fh:
            rows = list(_csv.DictReader(fh))
        if not rows:
            absent.append(s)
            continue
        for row in rows:
            label = str(row["setting"])
            slot = merged.setdefault(label, {"dbr_sd_value": row.get("dbr_sd_value") or None,
                                             "rate": {}, "n_scored": {}, "n_called": {}})
            rate = row.get("rate_over_scored")
            slot["rate"][s] = float(rate) if str(rate).strip() else None
            slot["n_scored"][s] = row.get("n_scored")
            slot["n_called"][s] = row.get("n_called")
    if absent:
        raise Refusal(
            f"04_doublet_sweep: no swept rate was recorded for {', '.join(absent)}. A sweep "
            f"missing a library is not a sweep of this cohort, and the spread across the "
            f"remainder would still print as one.")
    if not merged:
        raise Refusal("04_doublet_sweep: no setting produced a rate, so nothing was swept.")

    applied = task.params.get("dbr_sd_applied")
    out = _tables(pipeline) / "doublet_sweep.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh)
        w.writerow(["sample", "setting", "dbr_sd_value", "n_scored", "n_called",
                    "rate_over_scored", "dbr_sd_applied"])
        for label in sorted(merged):
            slot = merged[label]
            for s in samples:
                w.writerow([s, label, slot["dbr_sd_value"], slot["n_scored"].get(s),
                            slot["n_called"].get(s), slot["rate"].get(s), applied])

    results = [d.SweepResult(setting=label,
                             per_sample_rate={s: v for s, v in slot["rate"].items()
                                              if v is not None},
                             deep_decile_rate=None)
               for label, slot in sorted(merged.items())]
    for r in results:
        print(f"    {r}")
    rec = d.recommend(results)
    print(f"    {rec}")
    print("    the deepest-decile arm was NOT evaluated - it is a cohort quantity and this "
          "sweep measured it per library; see _doublet_sweep_stage")
    # Every metric a scalar or a plain string: the report renders each as one row of a findings
    # table, and `rec.rejected` is a mapping whose values are paragraphs. Flattened rather than
    # dropped, because WHY a setting was rejected is the useful half of a recommendation.
    rejected_text = "; ".join(f"{k}: {v}" for k, v in sorted(rec.rejected.items())) or "none"
    return {"outputs": [str(out)],
            "metrics": {"settings": ", ".join(sorted(merged)), "samples": len(samples),
                        "recommended": rec.setting, "reason": rec.reason,
                        "rejected": rejected_text,
                        "deep_decile_evaluated": False,
                        "applied": applied},
            "versions": {}}


# --------------------------------------------------------------------------------------------
# step 5 - measure a valley per library, propose one cohort constant


def _scanpy(pipeline, op, h5, prefix, params, log, python_exe):
    from adapters import scanpy_ops as so
    return so.run_scanpy_op(op=op, h5ad_in=h5, out_prefix=prefix, params=params,
                            log=log, executor=pipeline.executor, python_exe=python_exe)


#: The valley op measures BOTH count axes in one pass over each object. Step 5 derives a floor per
#: metric, and a metric silently not measured yields a cohort whose floor for it was never
#: proposed - which reads identically to one that was proposed and accepted.
VALLEY_METRICS = ("umi", "genes")


def _require_gene_patterns(step: str, p: dict) -> None:
    """Refuse a step that needs the gene-class patterns and was not given them.

    Shared by step 5 and step 6 because both measure the mitochondrial percentage, and a step that
    got its own copy of this check is a step that will one day not have it.
    """
    for k in ("mt_prefix", "ribo_pattern"):
        if not p.get(k):
            raise Refusal(
                f"{step} ({p['sample']}): `{k}` is not declared. It is species-specific and "
                f"nothing here will guess it: mouse and human differ in case alone for the "
                f"mitochondrial prefix, and a wrong one gives every cell pct_counts_mt == 0, "
                f"which is indistinguishable from a clean library and passes every mitochondrial "
                f"gate.\n"
                f"    mouse: mt_prefix=mt-  ribo_pattern=^Rp[sl]\n"
                f"    human: mt_prefix=MT-  ribo_pattern=^RP[SL]\n"
                f"    Add them as samplesheet columns. Matching is case-sensitive by design.")


def _quality_measure(task, pipeline, log):
    """Measure ONE library: both count valleys and the mitochondrial quartiles, in one pass.

    Split out of _quality_stage, which looped over every library inside a single task. Concurrency
    releases TASKS, not iterations, so ten libraries were measured one after another whatever
    --jobs said - and the measurement is the expensive half of step 5: reading the object and
    running calculate_qc_metrics over a 39k x 34k sparse matrix. The KDE that follows is cheap.
    """
    p = task.params
    _require_gene_patterns("05_quality", p)
    s = p["sample"]
    res = _scanpy(pipeline, "valley",
                  pipeline.results / "objects" / f"{s}_ambient.h5",
                  pipeline.scratch / f"{s}_qc",
                  {"sample": s, "metrics": list(VALLEY_METRICS),
                   "mt_prefix": p["mt_prefix"], "ribo_pattern": p["ribo_pattern"],
                   # The mitochondrial quartiles are taken above this floor, below the
                   # derivation ceiling, and the VALLEYS take neither. One pass, two populations,
                   # deliberately - see _op_valley for the measurement behind both restrictions.
                   # `mito_derivation_max` is passed only when this run declares one; absent, the
                   # adapter's own default applies and records itself in the ceiling table, so
                   # there is no route by which the applied value goes unrecorded.
                   "mito_floor_umi": p["light_floor"],
                   # Passed only when the samplesheet names one. Absent, no nuclear fraction is
                   # read, none is written, and the run is what it was before this existed.
                   **({"cellreads_stats": p["cellreads_stats"],
                       "nf_antisense": p.get("nf_antisense", False)}
                      if p.get("cellreads_stats") else {}),
                   **({"mito_derivation_max": p["mito_derivation_max"]}
                      if p.get("mito_derivation_max") is not None else {})},
                  log, p["python_exe"])
    m = res.get("metrics", {}) or {}
    got, bim = m.get("valleys") or {}, m.get("bimodal") or {}
    for metric in VALLEY_METRICS:
        print(f"    {s:<14}{metric:<7}valley {str(got.get(metric)):>8}   "
              f"bimodal {bim.get(metric)}")
    # Carried on the task RESULT so the barrier reads it from the manifest rather than from a
    # shared mutable the workers would race on.
    # The path to the denoiser's call travels with the result, so step 2 consumes this one pass
    # rather than opening the object again.
    called_csv = next((o for o in res.get("outputs", []) if str(o).endswith(".called_barcodes.csv")),
                      None)
    # THE CURVE THE VALLEY WAS FOUND ON, published rather than left in scratch. This is F6's
    # only possible input: `valleys_umi.csv` carries the position of the minimum and nothing
    # about the shape it sits in, and the whole question F6 answers - does the cut fall where the
    # data separates - is about the shape.
    density = _promote(pipeline, res, ".valley_density.csv", f"{s}.valley_density.csv",
                       what=f"05_quality ({s})")
    # PROMOTED to tables/, because step 7 reads it and a criterion whose input lives only in
    # scratch is one nobody can re-check the filter against.
    nf_csv = _promote(pipeline, res, ".nuclear_fraction.csv", f"{s}.nuclear_fraction.csv",
                      what=f"05_quality ({s})") if m.get("nf_median") is not None else None
    outputs = list(res.get("outputs", [])) + ([density] if density else []) \
        + ([nf_csv] if nf_csv else [])
    return {"outputs": outputs,
            "metrics": {"valleys": got, "bimodal": bim,
                        "mito_quartiles": m.get("mito_quartiles"),
                        # The ceiling's population travels with the ceiling. A threshold whose
                        # derivation population is not recorded cannot be compared with anyone
                        # else's, which is the entire disagreement this run went and settled.
                        "mito_population": m.get("mito_population"),
                        "called_barcodes": str(called_csv) if called_csv else None,
                        "n_called_by_denoiser": m.get("n_called_by_denoiser"),
                        # The nuclear fraction travels as METRICS, like everything else the
                        # barrier reads, so no worker writes to shared state.
                        "nf_median": m.get("nf_median"),
                        "nf_n_in_median": m.get("nf_n_in_median"),
                        "nf_n_joined": m.get("nf_n_joined"),
                        "nf_n_defined": m.get("nf_n_defined"),
                        "nf_join_pct": m.get("nf_join_pct"),
                        "nf_source": m.get("nf_source"),
                        "nf_columns": m.get("nf_columns"),
                        "nf_csv": str(nf_csv) if nf_csv else None},
            "versions": res.get("versions", {})}


def _quality_stage(task, pipeline, log):
    """The barrier: one cohort constant per count axis, plus the per-library mito ceilings.

    Reads what the ten `05_quality/<sample>` tasks measured. Nothing is measured here, so this is
    cheap and runs once - which is the point of splitting it out.
    """
    valleys = {m: [] for m in VALLEY_METRICS}
    mito_stats = {}
    mito_pop = {}
    nf_per = {}
    missing = []
    for s in task.params["samples"]:
        r = pipeline.results_by_key.get(f"05_quality/{s}")
        met = (getattr(r, "metrics", None) or {}) if r is not None else {}
        got, bim = met.get("valleys") or {}, met.get("bimodal") or {}
        if not got:
            missing.append(s)
            continue
        for metric in VALLEY_METRICS:
            valleys[metric].append({"sample": s, "value": got.get(metric),
                                    "bimodal": bim.get(metric)})
        q = met.get("mito_quartiles")
        if q:
            mito_stats[s] = q
        mito_pop[s] = met.get("mito_population") or {}
        nf_per[s] = {k: met.get(k) for k in
                     ("nf_median", "nf_n_in_median", "nf_median_floor_umi", "nf_n_joined",
                      "nf_n_defined", "nf_join_pct", "nf_columns", "nf_source", "nf_csv")}
    if missing:
        raise Refusal(
            f"05_quality: no measurement was recorded for {', '.join(missing)}. A library whose "
            f"valley was never measured cannot contribute to a cohort constant, and treating its "
            f"absence as agreement lets the other libraries decide on its behalf.")

    out = {"outputs": [], "metrics": {}, "versions": {}}
    for metric in VALLEY_METRICS:
        unknown = [v["sample"] for v in valleys[metric]
                   if v["value"] is None or v["bimodal"] is None]
        if unknown:
            raise Refusal(
                f"05_quality ({metric}): no valley was established for {', '.join(unknown)}. A "
                f"library with no measured valley cannot contribute to a cohort constant, and "
                f"treating its absence as agreement lets the other libraries decide on its "
                f"behalf.")
        sub = Task(key=task.key, step=task.step, fn=_quality,
                   params={"valleys": valleys[metric], "metric": metric,
                           "light_floor": task.params.get("light_floor")})
        r = _quality(sub, pipeline, log)
        out["outputs"] += list(r.get("outputs", []))
        for k, v in (r.get("metrics") or {}).items():
            out["metrics"][f"{metric}_{k}"] = v

    out = _mito_ceiling_stage(task, pipeline, mito_stats, out, mito_pop)

    # ----------------------------------------------- the cohort nuclear-fraction floor
    #
    # AFTER the ceiling, and deliberately: the trigger is the LOWER MITOCHONDRIAL BOUND that stage
    # just applied, read from the metric it published rather than re-derived here. The two cannot
    # drift apart, and this criterion adds no threshold anyone has to justify separately.
    #
    # ONE NUMBER FOR THE COHORT, the median of the per-library medians - equal weight per LIBRARY.
    # Pooling every cell instead lets a 16,000-cell library outvote a 7,500-cell one, and the floor
    # becomes a statement about how many nuclei each animal yielded rather than about nuclear
    # fractions. A PER-LIBRARY floor would be worse: it ADAPTS, so the library whose fractions run
    # low gets the loosest floor, keeps what a stricter library loses, and the criterion applies a
    # different standard in each arm of the design - differential by construction.
    samples = list(task.params["samples"])
    nf_medians = {s: (nf_per.get(s) or {}).get("nf_median") for s in samples}
    nf_have = {s: v for s, v in nf_medians.items() if v is not None}
    nf_floor = nf_trigger = None
    if nf_have:
        if len(nf_have) != len(samples):
            absent = sorted(set(samples) - set(nf_have))
            raise Refusal(
                f"05_quality (nuclear fraction): {len(nf_have)} of {len(samples)} libraries have "
                f"one and {', '.join(absent)} do not. A floor pooled over some libraries and "
                f"applied to all of them filters the rest on a boundary measured entirely in "
                f"other animals, and which libraries had a source declared is not a random subset "
                f"of a design. Declare a source for every library, or for none.")
        from adapters import nuclear_fraction as _nfmod
        nf_floor = _nfmod.median(list(nf_have.values()))
        nf_trigger = out["metrics"].get("mito_bound_lo")
        if nf_trigger is None:
            raise Refusal(
                "05_quality (nuclear fraction): the mitochondrial bound was not published, so "
                "there is no trigger. The trigger is that bound by design; inventing one here "
                "would be a second threshold that can drift away from the ceiling's own.")
        nf_trigger = float(nf_trigger)
        print(f"    nuclear fraction: per-library medians "
              f"{min(nf_have.values()):.4f}-{max(nf_have.values()):.4f}; cohort floor "
              f"{nf_floor:.4f} (median of {len(nf_have)} library medians, equal weight)")
        print(f"      joint criterion armed: mt > {nf_trigger:g}% AND nf < {nf_floor:.4f} "
              f"- ADDITIVE to the ceiling; the trigger is the declared lower bound")
        import csv as _csvnf
        pnf = _tables(pipeline) / "nuclear_fraction.csv"
        with open(pnf, "w", newline="", encoding="utf-8") as fh:
            w = _csvnf.writer(fh)
            w.writerow(["sample", "nf_median", "n_in_median", "median_floor_umi", "n_joined",
                        "n_defined", "join_pct", "columns", "source",
                        "cohort_floor", "trigger_pct"])
            for s in samples:
                r = nf_per.get(s) or {}
                w.writerow([s,
                            f"{r.get('nf_median'):.6f}" if r.get("nf_median") is not None else "",
                            r.get("nf_n_in_median"), r.get("nf_median_floor_umi"),
                            r.get("nf_n_joined"), r.get("nf_n_defined"),
                            f"{r.get('nf_join_pct'):.2f}" if r.get("nf_join_pct") is not None
                            else "", r.get("nf_columns"), r.get("nf_source"),
                            f"{nf_floor:.6f}", f"{nf_trigger:g}"])
        print(f"      per-library medians and the cohort floor: {pnf}")
        out["outputs"] = list(out.get("outputs", [])) + [str(pnf)]

    # PUBLISHED on the barrier's metrics - the channel `_apply_thresholds` already reads. Null
    # where no library declared a source: step 7 then sees three nulls and records the criterion
    # as NOT EVALUATED, which is a different fact from it firing on nothing.
    out["metrics"]["nf_floor"] = nf_floor
    out["metrics"]["nf_trigger_pct"] = nf_trigger
    out["metrics"]["nf_csv_by_sample"] = ({s: (nf_per.get(s) or {}).get("nf_csv")
                                           for s in samples} if nf_floor is not None else {})
    return out


def _mito_ceiling_stage(task, pipeline, mito_stats, out, mito_pop=None):
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
        # `pop_floor_umi` and `pop_n_all` are not decoration. A ceiling means nothing without the
        # population it was taken over: two implementations of the same published rule, applied to
        # the same libraries, disagreed on 7 of 10 ceilings for no other reason than this column
        # being absent on both sides, and neither could tell which one was different.
        # `derived` is the APPLIED fence (MAD route); `tukey` is the independent cross-check and is
        # never applied. Both are written, with the skew ratio that explains any gap between them,
        # so a reader can see whether the two routes agreed instead of being told that they did.
        # `pop_derivation_max_pct` and the two counts beside it complete the population record the
        # previous three columns began. The ceiling is derived over cells above `pop_floor_umi`
        # AND below `pop_derivation_max_pct`; `pop_max_above_floor` is the true observed maximum
        # taken BEFORE that cut, so the one number that shows a library is full of ambient is not
        # the number the derivation population erased.
        w.writerow(["sample", "n", "median", "q1", "q3", "iqr", "mad", "smad", "skew_iqr_over_smad",
                    "derived", "tukey_crosscheck", "tukey_over_derived", "ceiling", "clamped",
                    "mad_k", "iqr_mult", "assay", "bound_lo", "bound_hi",
                    "pop_floor_umi", "pop_n_all_called", "pop_derivation_max_pct",
                    "pop_n_above_floor", "pop_n_excluded_high_mito", "pop_max_above_floor"])
        for s in samples:
            m = ceil[s]
            pp = (mito_pop or {}).get(s) or {}
            mx = pp.get("max_above_floor")
            w.writerow([s, m.n, f"{m.median:.6f}", f"{m.q1:.6f}", f"{m.q3:.6f}",
                        f"{m.iqr:.6f}", f"{m.mad:.6f}", f"{m.smad:.6f}",
                        "" if m.skew_ratio is None else f"{m.skew_ratio:.6f}",
                        f"{m.derived:.6f}", f"{m.tukey:.6f}",
                        "" if m.cross_check is None else f"{m.cross_check:.6f}",
                        f"{m.ceiling:.6f}", m.clamped,
                        d.get("k", ""), d["mult"], assay, lo, hi,
                        pp.get("floor_umi", ""), pp.get("n_all_with_a_value", ""),
                        pp.get("derivation_max_pct", ""), pp.get("n_above_floor", ""),
                        pp.get("n_excluded_high_mito", ""),
                        "" if mx is None else f"{float(mx):.6f}"])

    out = dict(out)
    out["outputs"] = list(out.get("outputs", [])) + [str(p)]
    out["metrics"] = dict(out.get("metrics", {}))
    # k, WHERE IT CAME FROM, AND WHETHER THE BOUND DECIDED IT - carried into the metrics because
    # the report's parameter table reads them. It used to state "DERIVED, k derived from each
    # library's Tukey fence" as fixed text, which stopped being true for snRNA on 2026-08-13 and
    # would have had the report describing a filter the run did not apply.
    out["metrics"].update({
        "mito_ceiling_lo": min(m.ceiling for m in ceil.values()),
        "mito_ceiling_hi": max(m.ceiling for m in ceil.values()),
        "mito_bound_binds": sum(1 for m in ceil.values() if m.clamped),
        "mito_k": d.get("k"),
        "mito_k_source": d.get("k_source"),
        "mito_provenance": d.get("provenance"),
        "mito_assay": assay,
        "mito_bound_lo": lo,
        "mito_bound_hi": hi,
        "mito_derivation_max_pct": next(
            (pp.get("derivation_max_pct") for pp in (mito_pop or {}).values()
             if pp and pp.get("derivation_max_pct") is not None), None)})
    return out


# --------------------------------------------------------------------------------------------
# step 6 - cluster, profile, flag


def _population_for_cluster(pipeline, sample):
    """The cells step 6 is specified to cluster: quality-filtered, doublets NOT applied.

    # rule-one: no-removal - this reads a table and a task's metrics and returns four numbers.

    `cluster_flags.py` requires "the step-5 object: quality-filtered". The object step 6 opens is
    the denoised FULL DROPLET matrix, so without this the flags describe empty droplets.

    WHERE EACH NUMBER COMES FROM, and why it is not one place:

      * the per-library mitochondrial ceiling from `mito_ceiling_per_sample.csv`, which the
        05_quality task writes before step 6 runs;
      * the two cohort count floors from that task's METRICS, not from
        `thresholds_per_sample.csv`. That table is written by a LATER task - measured on
        2026-08-11 at 17:29 against step 6 at 17:19 - so reading it here found nothing. The
        floors are already in `results_by_key["05_quality"].metrics` as `umi_proposed` and
        `genes_proposed`; the same lookup is used for the decisions layer.

    The doublet criterion is deliberately absent: criterion D is only computable before a
    removal, since afterwards every cluster is 0% doublet by construction.

    REFUSES rather than returning None. A silent fallback here clustered the full droplet matrix
    on 2026-08-11 while the run reported success, which is the failure this whole lookup exists
    to prevent.
    """
    import csv as _csv

    tdir = _tables(pipeline)
    ceiling_path = tdir / "mito_ceiling_per_sample.csv"
    ceiling = None
    try:
        with open(ceiling_path, newline="", encoding="utf-8") as fh:
            for row in _csv.DictReader(fh):
                if row.get("sample") == sample:
                    ceiling = float(row["ceiling"])
                    break
    except (OSError, KeyError, ValueError) as e:
        raise Refusal(
            f"06_cluster_check ({sample}): could not read the mitochondrial ceiling from "
            f"{ceiling_path} ({type(e).__name__}: {e}). Step 6 must cluster the cells that reach "
            f"the deliverable; clustering the droplet matrix instead is the defect this lookup "
            f"exists to fix, so it stops rather than falling back to it.") from None
    if ceiling is None:
        raise Refusal(
            f"06_cluster_check ({sample}): {ceiling_path} carries no row for this library, so "
            f"the population it should be clustered on is unknown.")

    cohort = (getattr(pipeline.results_by_key.get("05_quality"), "metrics", None) or {})
    umi, gene = cohort.get("umi_proposed"), cohort.get("genes_proposed")
    if umi is None or gene is None:
        raise Refusal(
            f"06_cluster_check ({sample}): the 05_quality task recorded no "
            f"{'umi_proposed' if umi is None else 'genes_proposed'} in its metrics, so the count "
            f"floors the deliverable will use are not knowable here. Do NOT read them from "
            f"thresholds_per_sample.csv - that table is written after this step runs.")

    return {"cell_call_key": "cellbender_cell", "umi_floor": float(umi),
            "gene_floor": float(gene), "mito_ceiling": float(ceiling)}


def _cluster(task, pipeline, log):
    p = task.params
    _require_gene_patterns("06_cluster_check", p)
    res = _scanpy(pipeline, "cluster",
                  pipeline.results / "objects" / f"{p['sample']}_ambient.h5",
                  pipeline.scratch / f"{p['sample']}_clusters",
                  {"sample": p["sample"],
                   "resolution": p["resolution"], "seed": p["seed"],
                   # Profiled beside the applied resolution, into sibling files. Forwarded
                   # explicitly rather than by passing `p` through: the op's parameters are a
                   # declared contract, and a task parameter that never reaches it is a feature
                   # that silently does nothing while every test of the wiring still passes.
                   "extra_resolutions": p.get("extra_resolutions"),
                   # Needed to measure this object at all: it is the denoised one, so nothing
                   # has computed total_counts or pct_counts_mt on it, and cluster() refuses to
                   # normalise over an object whose depth was never measured.
                   "mt_prefix": p["mt_prefix"], "ribo_pattern": p["ribo_pattern"],
                   # The doublet calls are ATTACHED here and applied nowhere. Step 6's criterion
                   # D is only computable before a removal - afterwards every cluster is 0%
                   # doublet by construction - and this is the last point at which that holds.
                   "doublet_csv": p["doublet_csv"],
                   "doublet_key": "doublet_class",
                   "doublet_positive": "doublet",
                   # The population cluster_flags.py specifies. Step 5 has already derived every
                   # floor and this library's ceiling, so they are declared here rather than
                   # re-derived: step 6 must cluster the cells that reach the deliverable, not
                   # the droplet matrix they were selected from. The doublet criterion is
                   # deliberately absent - see the op.
                   "population": _population_for_cluster(pipeline, p["sample"]),
                   # THE COORDINATES F10 AND F11 ARE DRAWN ON, over a WIDER population than the
                   # clustering. `cell_called` is not a looser version of the clustering
                   # population, it is the only one that can answer F11: a projection built
                   # after the floors and the ceiling have been applied contains no nucleus they
                   # removed, so "did the removed nuclei leave as a population?" could only ever
                   # come back "there were none to look at". The clustering itself is unchanged.
                   "embedding": {"population": "cell_called"}},
                  log, p["python_exe"])
    # Published for the report, per `_promote`. A run whose embedding was declared null writes no
    # such file and the report says the coordinates were never asked for.
    emb = _promote(pipeline, res, ".embedding.csv", f"{p['sample']}.embedding.csv",
                   what=f"06_cluster ({p['sample']})")
    if emb:
        res = {**res, "outputs": list(res.get("outputs", [])) + [emb]}
    return res


def _cluster_flags(task, pipeline, log):
    import csv as _csv

    cf = step_module("cluster_flags")

    # THE PROFILES COME FROM THE MANIFEST, NOT FROM A GLOB.
    #
    # This globbed `pipeline.work` for `*profile*.csv` while the profiles are written under
    # `pipeline.scratch`, so it matched nothing and the step refused with "no cluster profile was
    # produced" - an honest refusal about something that had in fact been produced. A glob also
    # cannot tell this run's output from a previous one's, and it silently accepts nine files
    # where ten libraries were profiled. Each task reported the file it wrote and the orchestrator
    # checked that file exists, so ask it.
    profiles, absent = [], []
    for s in task.params["samples"]:
        r = pipeline.results_by_key.get(f"06_cluster/{s}")
        got = [o for o in (getattr(r, "outputs", None) or [])
               if str(o).endswith(".cluster_profile.csv")]
        profiles.extend(got)
        if not got:
            absent.append(s)
    if absent:
        raise Refusal(
            f"06_cluster_check: no cluster profile was recorded for {', '.join(absent)}. A "
            f"library that was never profiled has no flagged clusters and no unflagged ones, and "
            f"folding it in as neither lets the other libraries decide on its behalf.")

    rows: list = []
    for f in profiles:
        try:
            rows.extend(cf.read_profile_csv(f))
        except cf.ClusterRefusal as e:
            raise Refusal(f"06_cluster_check: {e}") from None
    if not rows:
        raise Refusal("06_cluster_check: no cluster profile was produced, so no cluster was "
                      "examined. That is not the same as no cluster being flagged.")

    # THRESHOLDS ARE PROPOSED FROM THIS COHORT, and anything declared overrides them by name.
    #
    # They used to be `d.get("a_umi_fraction", 0.5)` and friends - the calibration cohort's own
    # numbers, hardcoded here, under a `source` string that read "PROPOSED from the cohort". The
    # module says in its own header that these do not transfer: B is the p95 of one cohort's
    # cluster mito and C works only where that distribution was bimodal. A number carried from
    # another dataset while labelled as measured from this one is the worst of both.
    try:
        proposed = cf.propose(rows)
    except cf.ClusterRefusal as e:
        raise Refusal(f"06_cluster_check: {e}") from None
    d = (task.params.get("decisions") or {}).get("cluster_check") or {}
    declared = {k: float(d[y]) for k, y in (("a_umi_frac", "a_umi_fraction"),
                                            ("b_pct_mt", "b_mito_pct"),
                                            ("c_uninformative", "c_uninformative_pct"))
                if d.get(y) is not None}
    if declared:
        base = {"a_umi_frac": proposed.a_umi_frac, "b_pct_mt": proposed.b_pct_mt,
                "c_uninformative": proposed.c_uninformative}
        thr = cf.Thresholds(**{**base, **declared},
                            source=(f"decisions.yml declares {', '.join(sorted(declared))}; "
                                    f"the rest is {proposed.source}"))
    else:
        thr = proposed
    print(f"    thresholds: {thr}")

    flagged = cf.apply_flags(rows, thr)
    out = _tables(pipeline) / "cluster_profile.csv"
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=sorted({k for r in flagged.rows for k in r}))
        w.writeheader()
        w.writerows(flagged.rows)

    # WHAT POPULATION WAS CLUSTERED, recorded because it is not the one this step's module
    # specifies. cluster_flags.py opens by saying the check runs on the quality-filtered object
    # with the doublet flags attached; no such object exists here, and in evidence mode none can,
    # because nothing is removed. So the clustering is over the denoised object as delivered -
    # empty droplets included - and a table that does not say so reads exactly like one measured
    # on nuclei. Measured, not judged: a cluster whose median barcode carries no counts is a
    # statement about the input, and no threshold of ours decided it.
    empty = [r for r in rows if r.get("median_umi") in (None, 0, 0.0)]
    cells = sum(int(r.get("n") or 0) for r in rows)
    if empty:
        pipeline.record_findings("06_cluster_check", [{
            "step": "06_cluster_check", "check": "population clustered", "severity": "REVIEW",
            "message": (
                f"{len(empty)} of {len(rows)} clusters have a median UMI of zero, holding "
                f"{sum(int(r.get('n') or 0) for r in empty):,} of {cells:,} barcodes. This step "
                f"clustered the denoised object as delivered: no count floor and no cell call "
                f"have been applied to it, because evidence mode removes nothing and this "
                f"pipeline builds no quality-filtered object for it to read. The flags below "
                f"therefore describe a mixture of nuclei and empty droplets rather than the "
                f"population modules/06_cluster_check/cluster_flags.py specifies, and the "
                f"thresholds proposed from it inherit that."),
            "detail": ["clusters per library: " + ", ".join(
                f"{s}={sum(1 for r in rows if r.get('sample') == s)}"
                for s in task.params["samples"])]}])

    fired, unknown = flagged.counts(), flagged.unknown_counts()
    for k in ("A", "B", "C", "FLAG", "WATCH"):
        print(f"    {k:<6}fired {fired[k]:>4}   not evaluated {unknown[k]:>4}   "
              f"of {len(flagged.rows)}")
    # Both numbers, always. "How many fired" alone cannot distinguish a criterion that cleared
    # every cluster from one that was never computed on any of them.
    metrics = {"clusters": len(flagged.rows), "libraries": len(task.params["samples"]),
               "thresholds": str(thr)}
    metrics.update({f"{k}_fired": v for k, v in fired.items()})
    metrics.update({f"{k}_not_evaluated": v for k, v in unknown.items()})

    # THE EXTRA RESOLUTIONS, FLAGGED THE SAME WAY, INTO SIBLING TABLES.
    #
    # `cluster_profile.csv` above is untouched and remains the only table step 7 and every
    # downstream stage read. These carry the SAME columns and the SAME rule, and differ from it
    # in the resolution and in nothing else - which is the whole point of having them: a flag
    # that appears at one resolution and not another is a statement about the clustering, and one
    # that appears at all of them is a statement about the cells.
    #
    # Thresholds are proposed WITHIN each resolution, exactly as the default's are proposed
    # within its own rows. Carrying the default's thresholds across would compare a resolution
    # against a cut-point derived from a different number of clusters.
    extras_written, extras_failed = [], {}
    by_res: dict = {}
    for s in task.params["samples"]:
        r = pipeline.results_by_key.get(f"06_cluster/{s}")
        for o in (getattr(r, "outputs", None) or []):
            name = str(o)
            if ".cluster_profile.res" in name and name.endswith(".csv"):
                by_res.setdefault(name.rsplit(".cluster_profile.res", 1)[1][:-4], []).append(o)
    for res_txt, files in sorted(by_res.items()):
        try:
            erows: list = []
            for f in files:
                erows.extend(cf.read_profile_csv(f))
            ethr = cf.propose(erows)
            eflag = cf.apply_flags(erows, ethr)
            eout = _tables(pipeline) / f"cluster_profile.res{res_txt}.csv"
            with open(eout, "w", encoding="utf-8", newline="") as fh:
                w = _csv.DictWriter(fh, fieldnames=sorted({k for r in eflag.rows for k in r}))
                w.writeheader()
                w.writerows(eflag.rows)
            efired = eflag.counts()
            print(f"    resolution {res_txt}: {len(eflag.rows)} clusters, "
                  f"FLAG {efired['FLAG']}, WATCH {efired['WATCH']}  -> {eout.name}")
            extras_written.append(str(eout))
            metrics[f"res{res_txt}_clusters"] = len(eflag.rows)
            metrics.update({f"res{res_txt}_{k}_fired": v for k, v in efired.items()})
        except Exception as e:                              # noqa: BLE001 - named, not hidden
            # An extra resolution is a diagnostic. It must not be able to cost the run its
            # default flags, and it must not fail silently either.
            extras_failed[res_txt] = f"{type(e).__name__}: {e}"
            print(f"    resolution {res_txt}: FAILED - {type(e).__name__}: {e}")
    if extras_failed:
        metrics["extra_resolutions_failed"] = extras_failed
    return {"outputs": [str(out)] + extras_written, "metrics": metrics, "versions": {}}


# --------------------------------------------------------------------------------------------
# step 7 - the only step that removes anything


#: `mito_ceiling_pct: per_library` in decisions.yml means "apply what step 5 derived for each
#: library", rather than one number for all ten. It is a word and not a number because that is
#: what is being approved: a single cohort ceiling and a per-library fence are different
#: decisions, and on the calibration cohort the difference is a removal that falls 2.35x more
#: heavily on one arm of the design than the other against 1.23x - the rule-one Q3 quantity.
PER_LIBRARY = "per_library"


def _ceilings_for(pipeline, samples, declared):
    """(ceiling per sample, how it was decided). Refuses rather than defaulting.

    A single declared number is applied to every library. `per_library` reads step 5's own table,
    which is the only place the derived fences exist.
    """
    if str(declared).strip().lower() != PER_LIBRARY:
        try:
            v = float(declared)
        except (TypeError, ValueError):
            raise Refusal(
                f"07_apply: quality.mito_ceiling_pct is {declared!r}, which is neither a number "
                f"nor {PER_LIBRARY!r}. There is no default: one cohort ceiling and a per-library "
                f"fence remove different populations.") from None
        return {s: v for s in samples}, f"one declared ceiling of {v:g}% for every library"

    import csv as _csv
    p = _tables(pipeline) / "mito_ceiling_per_sample.csv"
    if not p.exists():
        raise Refusal(
            f"07_apply: {PER_LIBRARY} was approved but {p} does not exist, so there is nothing "
            f"to apply. Run step 5 in evidence mode first - the ceilings are derived there.")
    with p.open(encoding="utf-8", newline="") as fh:
        got = {r["sample"]: float(r["ceiling"]) for r in _csv.DictReader(fh) if r.get("ceiling")}
    missing = [s for s in samples if s not in got]
    if missing:
        raise Refusal(
            f"07_apply: {p} has no ceiling for {', '.join(missing)}. A library filtered on "
            f"another library's ceiling is filtered on a threshold nobody derived for it.")
    lo, hi = min(got[s] for s in samples), max(got[s] for s in samples)
    return {s: got[s] for s in samples}, (f"the per-library fence step 5 derived, "
                                          f"{lo:.2f}-{hi:.2f}% across {len(samples)} libraries")


#: What step 7 cannot proceed without, and what each one is for. Read by the measure tasks and by
#: the apply task, so an incomplete decisions file is refused identically wherever it is noticed.
APPLY_REQUIRES = {
    "quality.umi_floor": "the UMI floor",
    "quality.gene_floor": "the gene floor",
    "quality.mito_ceiling_pct": "the mitochondrial ceiling",
}


def _attested(block):
    """The value of a decisions entry, if it is properly attested. None otherwise.

    A number with no `approved_by` and no `verbatim` is not a decision, it is a number somebody
    typed - so it does not count as ADJUDICATED and the derived value is used instead. Half an
    approval must not read as a whole one.
    """
    if not isinstance(block, dict):
        return None
    v = block.get("value")
    if v is None or str(v).strip() == "":
        return None
    if not str(block.get("approved_by") or "").strip():
        return None
    if not str(block.get("verbatim") or "").strip():
        return None
    return v


def _cohort_metric(pipeline, key):
    """A value the 05_quality barrier published, or None. Never a default.

    Step 7 reads what step 5 derived through the same channel `_apply_thresholds` uses, so no
    threshold can reach the filter by a route the rest of the run did not travel.
    """
    return (getattr(pipeline.results_by_key.get("05_quality"), "metrics", None) or {}).get(key)


def _nf_csv_for(pipeline, sample):
    """This library's per-barcode nuclear-fraction table, or None where none was written."""
    return (_cohort_metric(pipeline, "nf_csv_by_sample") or {}).get(sample)


def _apply_thresholds(task, pipeline, samples):
    """(values, who decided each, ceiling per sample, how the ceilings were decided).

    A THRESHOLD IS EITHER ADJUDICATED OR DERIVED, AND THE RESULT SAYS WHICH. A decisions file is
    optional: without one the pipeline applies what it derived in evidence mode, which is the
    honest default because those values are a function of the data rather than of anybody's
    preference. What must never happen is a derived value being reported as a chosen one, so the
    class travels with the value into the ledger, the object and the report.

    Resolved when the task RUNS, not when the graph is built: with `per_library` the ceilings come
    from step 5's table, and on a clean run that file does not exist yet while the graph is being
    assembled. Both step-7 tasks come through here, so they cannot resolve differently.
    """
    q = (task.params.get("decisions") or {}).get("quality") or {}
    cohort = (getattr(pipeline.results_by_key.get("05_quality"), "metrics", None) or {})
    derived = {
        "quality.umi_floor": cohort.get("umi_proposed"),
        "quality.gene_floor": cohort.get("genes_proposed"),
        "quality.mito_ceiling_pct": PER_LIBRARY,
    }
    values, classes = {}, {}
    for dotted in APPLY_REQUIRES:
        leaf = dotted.split(".", 1)[1]
        declared = _attested(q.get(leaf))
        if declared is not None:
            values[dotted], classes[dotted] = declared, "ADJUDICATED"
            continue
        if derived[dotted] is None:
            raise Refusal(
                f"07_apply: {dotted} is neither declared nor derived. Step 5 proposes it in "
                f"evidence mode and this run has no record of that, so there is nothing to "
                f"apply - and a threshold nobody chose and nothing measured is not a default, it "
                f"is a guess.")
        values[dotted], classes[dotted] = derived[dotted], "DERIVED"
    ceilings, basis = _ceilings_for(pipeline, samples, values["quality.mito_ceiling_pct"])
    return values, classes, ceilings, basis


def _apply_measure(task, pipeline, log):
    """Measure every applied criterion for ONE library. Removes nothing, writes no object."""
    p = task.params
    _require_gene_patterns("07_apply", p)
    s = p["sample"]
    resolved, _classes, ceilings, _basis = _apply_thresholds(task, pipeline, [s])
    # RESULTS, NOT SCRATCH. This table is the audit trail of the removal: every barcode the
    # library held, the four measured values, and which of the five criteria fired on each. The
    # ledger names only what LEFT, so without this there is no released record of what stayed or
    # of how close a retained nucleus came to a threshold. It is the file a reader re-checks the
    # filter with, and an intermediate is not something anyone is invited to open.
    return _scanpy(pipeline, "apply_measure",
                   pipeline.results / "objects" / f"{s}_ambient.h5",
                   _tables(pipeline) / f"{s}",
                   {"sample": s,
                    "mt_prefix": p["mt_prefix"], "ribo_pattern": p["ribo_pattern"],
                    "umi_floor": resolved["quality.umi_floor"],
                    "gene_floor": resolved["quality.gene_floor"],
                    "mito_ceiling_pct": ceilings[s],
                    "light_floor": p["light_floor"],
                    # ALWAYS PRESENT, null where the criterion is not armed. The op refuses a
                    # MISSING key precisely so that "not evaluated" and "wired wrong" cannot
                    # produce the same per-cell table.
                    "nf_csv": _nf_csv_for(pipeline, s),
                    "nf_floor": _cohort_metric(pipeline, "nf_floor"),
                    "nf_trigger_pct": _cohort_metric(pipeline, "nf_trigger_pct"),
                    "doublet_csv": p["doublet_csv"],
                    "doublet_key": "doublet_class", "doublet_positive": "doublet"},
                   log, p["python_exe"])


#: Columns of `removal_by_criterion.csv`, in order. Long format - one row per library per
#: criterion, plus one `ALL` row per criterion for the cohort - because a wide table gains a
#: column every time a criterion is added and every reader of it has to be changed.
_BREAKDOWN_COLUMNS = ("sample", "criterion", "n_in", "n_fired", "n_sole", "n_removed_any",
                      "pct_of_library", "pct_sole_of_library")

#: The literal a cohort row carries in its `sample` column. A row whose sample is blank reads as
#: a library whose name was lost, which is a different claim from a total over all of them.
BREAKDOWN_ALL = "ALL"


def _removal_breakdown(percell, criteria, masks, removed_mask) -> list:
    """Per library and per criterion: how many it fired on, and how many it removed ALONE.

    # rule-one: no-removal - this counts an already-decided removal and returns rows.

    `n_fired` counts every observation the criterion fired on, so the criteria overlap and their
    counts sum to more than `n_removed_any`. `n_sole` counts the ones NO other criterion would
    have removed, and that is the number that says whether a threshold did any work: a criterion
    with a large total and a sole count of zero removed nothing the others were not already
    removing, and changing it would change nothing at all.

    Both are needed and neither substitutes for the other. Reporting only the total makes every
    criterion look load-bearing; reporting only the sole count makes a criterion that agrees with
    its neighbours look inert when it may be the reason they agree.

    Built from `masks` rather than by re-reading the ledger. The ledger is written from these same
    booleans, so re-deriving the counts from the file would be a second route to one number with
    no cross-check attached - it can only agree, or disagree with nothing to say which is right.
    """
    named = {str(r.get("sample") or "") for r in percell}
    if BREAKDOWN_ALL in named:
        # A library actually called ALL would have its own row silently merged into the cohort
        # total, and the total would then be double-counted into itself. Refused rather than
        # renamed: renaming a library in one table and not the others is worse than stopping.
        raise Refusal(
            f"07_apply: a library is named {BREAKDOWN_ALL!r}, which is the label the per-criterion "
            f"breakdown uses for the cohort total. Its row and the total cannot be told apart. "
            f"Rename the library in the samplesheet.")

    per: dict = {}
    for i, row in enumerate(percell):
        s = row.get("sample") or BREAKDOWN_ALL
        d = per.setdefault(s, {"n_in": 0, "n_removed": 0,
                               "fired": {c: 0 for c in criteria},
                               "sole": {c: 0 for c in criteria}})
        d["n_in"] += 1
        if removed_mask[i]:
            d["n_removed"] += 1
        hit = [c for c in criteria if masks[c][i]]
        for c in hit:
            d["fired"][c] += 1
        if len(hit) == 1:
            d["sole"][hit[0]] += 1

    total = {"n_in": sum(d["n_in"] for d in per.values()),
             "n_removed": sum(d["n_removed"] for d in per.values()),
             "fired": {c: sum(d["fired"][c] for d in per.values()) for c in criteria},
             "sole": {c: sum(d["sole"][c] for d in per.values()) for c in criteria}}

    def pct(n, d):
        return round(100.0 * n / d, 4) if d else ""

    rows = []
    for s in list(per) + [BREAKDOWN_ALL]:
        d = total if s == BREAKDOWN_ALL else per[s]
        for c in criteria:
            rows.append({"sample": s, "criterion": c, "n_in": d["n_in"],
                         "n_fired": d["fired"][c], "n_sole": d["sole"][c],
                         "n_removed_any": d["n_removed"],
                         "pct_of_library": pct(d["fired"][c], d["n_in"]),
                         "pct_sole_of_library": pct(d["sole"][c], d["n_in"])})
    return rows


def _apply(task, pipeline, log):
    import csv as _csv

    from .decisions import action_string

    ap = step_module("apply")
    samples = list(task.params["samples"])
    resolved, classes, _ceilings, ceiling_basis = _apply_thresholds(task, pipeline, samples)
    apply_block = (task.params["decisions"].get("apply") or {})
    # WHO AUTHORISED THIS, RECORDED RATHER THAN DEMANDED.
    #
    # The pipeline no longer refuses without an operator's approval. It does not need to: nothing
    # here overwrites anything, the inputs are untouched and every run writes under its own
    # content-addressed directory, so a removal is recoverable by construction and the refusal was
    # guarding a hazard the layout has removed.
    #
    # What the approval ALSO did was attribute the thresholds, and that is not replaced by a
    # directory name - so it is recorded instead. Every threshold carries ADJUDICATED or DERIVED
    # into the ledger, the object and the report. A run with no decisions file applies what the
    # data proposed and says so; it never reports a proposal as a decision.
    #
    # When a decisions file IS present and complete, the gate still runs, and it still refuses on
    # words that do not match the action they were given for.
    verbatim = str(apply_block.get("verbatim") or "")
    declared_action = str(apply_block.get("action") or "").strip()
    action = action_string(resolved)
    adjudicated = [k for k, v in classes.items() if v == "ADJUDICATED"]
    authorised_by = (f"ADJUDICATED: {', '.join(sorted(adjudicated))}" if adjudicated
                     else "DERIVED - no decisions file; the thresholds are what this pipeline "
                          "measured, not what anyone chose")

    # Read through step 6's own reader, not csv.DictReader. `preflight` compares pct_doublet
    # against a float and tests FLAG with `is True`; from a raw DictReader every cell is a string,
    # so the first raised TypeError and the second silently found no flagged cluster in a table
    # of 126 of them - it reported that none of the deliverable sits in one.
    cf = step_module("cluster_flags")
    rows: list = []
    prof = _tables(pipeline) / "cluster_profile.csv"
    if prof.exists():
        try:
            rows = cf.read_profile_csv(prof)
        except cf.ClusterRefusal as e:
            raise Refusal(f"07_apply: step 6's profile cannot be read: {e}") from None
    # --- what the per-library measurements found. One table per library, covering every
    # barcode it saw - kept and removed alike - because a record of only the survivors cannot
    # be audited against anything.
    percell: list = []
    absent = []
    for s in samples:
        r = pipeline.results_by_key.get(f"07_measure/{s}")
        got = [o for o in (getattr(r, "outputs", None) or []) if str(o).endswith(".percell.csv")]
        if not got:
            absent.append(s)
            continue
        with open(got[0], encoding="utf-8", newline="") as fh:
            percell.extend(dict(row, _src=str(got[0])) for row in _csv.DictReader(fh))
    if absent:
        raise Refusal(
            f"07_apply: no per-cell criteria table was recorded for {', '.join(absent)}. A "
            f"library whose criteria were never measured cannot be filtered, and delivering the "
            f"other nine as though it had been is a cohort missing a library that says so "
            f"nowhere.")

    def flag(row, col):
        return str(row.get(col, "")).strip().lower() == "true"

    criteria = list(so_apply_criteria())
    n_in = len(percell)
    ids = [f"{r['sample']}|{r['barcode']}" for r in percell]
    masks = {c: [flag(r, c) for r in percell] for c in criteria}
    removed_mask = [flag(r, "removed") for r in percell]
    n_removed = sum(1 for x in removed_mask if x)

    # The pre-flight is told the real retained total, not a placeholder: its whole job is to say
    # what share of the DELIVERABLE sits in a cluster step 6 flagged.
    pipeline.record_findings("07_apply", [
        {"step": "07_apply", "check": f.check, "severity": f.severity,
         "message": f.message, "detail": list(f.detail or [])}
        for f in ap.preflight(rows, kept_total=n_in - n_removed)])

    # --- the ledger, BEFORE the gate. One row per removed barcode, listing every criterion that
    # fired on it, so what left can be named afterwards and re-read from the input.
    record = ap.build_removal_record(ids, masks)
    ledger = _tables(pipeline) / "removal_ledger.csv"
    ap.write_removal_record(record, ledger)

    # --- and the same removal counted, PER LIBRARY AND PER CRITERION. F9 draws this for the
    # cohort; the cohort is not where the question lives. A criterion that removes 2% of one
    # library and 42% of another has put a technical gradient exactly where the biology is
    # measured, and a single bar cannot show it. Built from `masks` - the same booleans the
    # ledger and the gate were built from - so the table and the removal cannot disagree.
    breakdown = _removal_breakdown(percell, criteria, masks, removed_mask)
    bpath = _tables(pipeline) / "removal_by_criterion.csv"
    with open(bpath, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=list(_BREAKDOWN_COLUMNS))
        w.writeheader()
        for r in breakdown:
            w.writerow(r)

    # --- the removal is checked and recorded. The arithmetic and the ledger are verified the
    # same way either route; only the approval is conditional.
    try:
        if adjudicated and verbatim.strip():
            if not declared_action:
                raise Refusal(
                    "07_apply: a decisions file declares thresholds and supplies approval words, "
                    "but apply.action is empty. The action is what the approval is matched "
                    "against; with it blank the gate would compare a value with itself.\n"
                    f"    This run's thresholds authorise: {action}")
            kept = ap.apply_removal(
                n_in=n_in, removed_mask_sum=n_removed, action=action, user_verbatim=verbatim,
                approvals={declared_action: verbatim},
                record=record, record_path=str(ledger))
            authorised_by = f"APPROVED by {apply_block.get('approved_by') or 'unnamed'}"
        else:
            kept = ap.record_removal(
                n_in=n_in, removed_mask_sum=n_removed, action=action,
                record=record, record_path=str(ledger), authorised_by=authorised_by)
    except ap.ApplyRefusal as e:
        raise Refusal(f"07_apply: {e}") from None

    # --- the audit, on the completed decision and still before anything is written.
    au = step_module("audit_removal")
    table = {c: [flag(r, c) for r in percell] for c in criteria}
    table.update({
        "cellbender_cell": [flag(r, "cellbender_cell") for r in percell],
        "keep": [flag(r, "keep") for r in percell],
        "removed": removed_mask,
        "total_counts": [float(r["total_counts"]) if r.get("total_counts") else 0.0
                         for r in percell],
        "doublet_scored": [flag(r, "doublet_scored") for r in percell],
    })
    findings = au.audit(table, criteria=criteria, scored_col="doublet_scored",
                        light_floor=task.params["light_floor"],
                        quality_floor=float(resolved["quality.umi_floor"]),
                        doublet_criterion="fail_doublet")
    pipeline.gate("07_apply", findings, au.verdict(findings))

    # --- and only now is a matrix touched.
    # `kept` is the APPROVED COUNT from the gate above and is checked against the object that
    # gets written, so it must survive this loop. Naming the per-library barcode list `kept`
    # too shadowed it, and the check below - the one that proves the delivered object holds the
    # population the gate approved - died on int(<list>) instead of comparing anything. It had
    # never run: nothing had previously reached a successful write for it to run on.
    keep_lists = []
    for s in samples:
        kept_barcodes = [r["barcode"] for r in percell if r["sample"] == s and flag(r, "keep")]
        kp = pipeline.scratch / f"{s}_apply.keep.txt"
        kp.write_text("\n".join(kept_barcodes) + "\n", encoding="utf-8")
        keep_lists.append({
            "sample": s, "keep_csv": str(kp),
            "h5ad": str(pipeline.results / "objects" / f"{s}_ambient.h5"),
            "annotations_csv": _deliverable_annotations(
                pipeline, ap, s, kept_barcodes, rows,
                [r for r in percell if r["sample"] == s])})
    # The run's own identity, carried onto the delivered objects so a reader holding one can say
    # which run produced it without finding this directory again. Resolved here because the
    # adapter runs out of process and has no pipeline to ask. Anything unresolvable is passed
    # EMPTY rather than guessed: a declaration naming a run key the object did not come from is
    # worse than one naming none.
    res = _scanpy(pipeline, "apply_write", pipeline.results / "objects" / f"{samples[0]}_ambient.h5",
                  pipeline.results / "objects" / "cohort",
                  {"libraries": keep_lists, "provenance": _run_identity(pipeline)},
                  log, task.params["python_exe"])

    delivered = (res.get("metrics") or {}).get("n_delivered")
    if delivered != int(kept):
        raise Refusal(
            f"07_apply: the gate approved {int(kept):,} observations and the object written "
            f"holds {delivered:,}. The deliverable and the ledger describe different "
            f"populations, and the one on disk is not the one that was approved.")
    print(f"    delivered {delivered:,} of {n_in:,}; {n_removed:,} removed, ledger at {ledger}")
    print(f"    ceilings: {ceiling_basis}")
    written = res.get("metrics") or {}
    print(f"    per-library objects: {len(written.get('per_sample_objects') or [])}"
          f"   obs carries: {', '.join(written.get('obs_columns') or [])}")
    return {"outputs": ([str(p) for p in res.get("outputs", [])]
                        + [str(ledger), str(bpath)]),
            "metrics": {"n_in": n_in, "n_delivered": delivered, "n_removed": n_removed,
                        "action": action, "ceilings": ceiling_basis,
                        "authorised_by": authorised_by,
                        "per_sample_objects": written.get("per_sample_objects") or [],
                        "obs_columns": written.get("obs_columns") or [],
                        "removal_breakdown": str(bpath),
                        **{f"n_{c}": sum(1 for x in masks[c] if x) for c in criteria}},
            "versions": res.get("versions", {})}


def _run_identity(pipeline) -> dict:
    """Run key, commit and version, for stamping onto the delivered objects.

    Every field degrades to an empty string rather than to a plausible value. A compute node
    often has no git, and `git_provenance` already returns a marker instead of a hash there;
    passing that marker on as though it were one would put a fiction in the object's own record.
    """
    from . import provenance as _prov

    repo = Path(__file__).resolve().parents[1]
    try:
        git = _prov.git_provenance(repo)
    except Exception:                                                     # noqa: BLE001
        git = {}
    commit = str(git.get("commit") or "")
    if commit and not all(c in "0123456789abcdef" for c in commit.lower()):
        commit = ""                       # "not a git checkout" is not a commit
    try:
        version = (repo / "VERSION").read_text(encoding="utf-8").strip()
    except OSError:
        version = ""
    return {"run_key": str(getattr(pipeline, "run_key", "") or ""),
            "commit": commit, "version": version}


def _deliverable_annotations(pipeline, ap, sample, kept, profile, percell) -> str:
    """Everything the delivered object should carry per barcode, for one library.

    The values the criteria were evaluated on, and what step 6 found about the cluster the
    nucleus sits in.

    Step 6 clusters each library and flags each cluster, and until this existed none of it reached
    the deliverable. The per-cell assignment was computed and discarded, so anything downstream
    had to re-cluster to ask a question step 6 had already answered - and re-clustering a
    filtered object does not even give the same answer, because criterion D is a tautology once
    the doublets are gone.

    The join is barcode -> cluster, from step 6's own labels, then cluster -> flags through
    `apply.annotate_kept`, which is the module's own function for exactly this and had never been
    called by anything.

    A barcode step 6 did not label, or a cluster absent from the profile, yields empty cells. An
    empty cell is not `False`: "this cluster was not flagged" and "this barcode was never
    examined" are different facts, and a deliverable that conflates them says a cluster check
    happened where none did.
    """
    import csv as _csv

    r = pipeline.results_by_key.get(f"06_cluster/{sample}")
    labels_path = next((o for o in (getattr(r, "outputs", None) or [])
                        if str(o).endswith(".cluster_labels.csv")), None)
    if labels_path is None:
        raise Refusal(
            f"07_apply ({sample}): step 6 recorded no per-barcode cluster labels, so its result "
            f"cannot be carried onto the nuclei it describes. The deliverable would be written "
            f"having never been cluster-checked, and nothing in it would say so.")

    with open(labels_path, encoding="utf-8", newline="") as fh:
        of_barcode = {str(row["barcode"]): row.get("cluster") or None
                      for row in _csv.DictReader(fh)}

    obs_rows = [{"barcode": b, "sample": sample, "cluster": of_barcode.get(b)} for b in kept]
    ap.annotate_kept(obs_rows, profile, cluster_key="cluster", sample_key="sample")

    # THE VALUES THE FILTER READ, CARRIED ONTO THE OBJECT IT PRODUCED.
    #
    # The delivered matrix used to arrive with no total_counts, no gene count, no mitochondrial
    # percentage and no doublet call - `qc_metrics()` is computed inside each op, on a freshly
    # loaded object, and was never written back - so a reader could not see why any nucleus had
    # survived, and had to recompute all four to find out. Every one of them was already in hand.
    #
    # Taken from the per-cell table rather than recomputed here, so the object provably carries
    # the numbers the criteria were evaluated on. Recomputing would give the same answers and
    # would not be able to prove it.
    measured = {str(r["barcode"]): r for r in (percell or [])}
    for row in obs_rows:
        m = measured.get(str(row["barcode"])) or {}
        for col in ("total_counts", "n_genes", "pct_counts_mt", "doublet_score",
                    "doublet_class"):
            v = m.get(col)
            row[col] = None if v is None or str(v).strip() == "" else v

    out = pipeline.scratch / f"{sample}_apply.annotations.csv"
    fields = ["barcode", "total_counts", "n_genes", "pct_counts_mt",
              "doublet_score", "doublet_class",
              "cluster", "cluster_FLAG", "cluster_WATCH",
              "cluster_pct_doublet", "cluster_median_pct_mt"]
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in obs_rows:
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fields})
    return str(out)


def so_apply_criteria():
    """The applied criteria, named by the adapter that writes them rather than repeated here."""
    from adapters.scanpy_ops import APPLY_CRITERIA
    return APPLY_CRITERIA


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
    "03_light_floor": (
        "Sets the depth below which a barcode is not handed to the doublet detector, because "
        "scDblFinder documents a 200-UMI floor to avoid erroring on near-empty droplets.",
        "Anything about QUALITY. It is a technical precondition of one tool and not a filter: "
        "nothing is removed here, and a barcode below the floor is reported as never examined "
        "for doublets - which is UNKNOWN, not a singlet. Reading it as a quality threshold is "
        "the mistake the separate name exists to prevent."),
    "04_doublets": (
        "Scores doublets per library above the light floor and checks the calls for health.",
        "Whether a called doublet IS one. It checks that the rate is a measurement rather than "
        "the prior, and that it does not fall unevenly across the design. Barcodes below the "
        "light floor were never scored and are reported as unknown, not as singlets."),
    "05_quality": (
        "Derives the count floors and the mitochondrial ceiling, in one pass over the same "
        "object and deliberately not the same population. The floors are the density valley "
        "measured per library and proposed as ONE cohort constant, over every barcode; the "
        "ceiling is each library's own median + k*1.4826*MAD taken over the barcodes above the "
        "light floor, with k DERIVED as the cohort median of the multiple at which each "
        "library's Tukey fence sits, and the whole bounded by a declared statement about what a "
        "nucleus can be. Tukey is carried as an independent second derivation and is never "
        "applied.",
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
    "_doublets", "_doublet_health", "_doublet_sweep", "_doublet_sweep_stage", "_light_floor",
    "_quality", "_quality_measure", "_quality_stage", "_cluster",
    "_cluster_flags", "_apply", "_apply_measure", "_ceilings_for", "_report",
]
