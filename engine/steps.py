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


def _tables(p) -> Path:
    return p.results / "tables"


def _objects(p) -> Path:
    return p.results / "objects"


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
            tmp_dir=pipeline.scratch / f"{sample}_extract",
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
                   # The mitochondrial quartiles are taken above this floor and the VALLEYS are
                   # not. One pass, two populations, deliberately - see _op_valley for the
                   # measurement behind it.
                   "mito_floor_umi": p["light_floor"]},
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
    return {"outputs": list(res.get("outputs", [])),
            "metrics": {"valleys": got, "bimodal": bim,
                        "mito_quartiles": m.get("mito_quartiles"),
                        # The ceiling's population travels with the ceiling. A threshold whose
                        # derivation population is not recorded cannot be compared with anyone
                        # else's, which is the entire disagreement this run went and settled.
                        "mito_population": m.get("mito_population"),
                        "called_barcodes": str(called_csv) if called_csv else None,
                        "n_called_by_denoiser": m.get("n_called_by_denoiser")},
            "versions": res.get("versions", {})}


def _quality_stage(task, pipeline, log):
    """The barrier: one cohort constant per count axis, plus the per-library mito ceilings.

    Reads what the ten `05_quality/<sample>` tasks measured. Nothing is measured here, so this is
    cheap and runs once - which is the point of splitting it out.
    """
    valleys = {m: [] for m in VALLEY_METRICS}
    mito_stats = {}
    mito_pop = {}
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

    return _mito_ceiling_stage(task, pipeline, mito_stats, out, mito_pop)


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
        w.writerow(["sample", "n", "median", "q1", "q3", "iqr", "derived", "ceiling", "clamped",
                    "iqr_mult", "assay", "bound_lo", "bound_hi",
                    "pop_floor_umi", "pop_n_all_called"])
        for s in samples:
            m = ceil[s]
            pp = (mito_pop or {}).get(s) or {}
            w.writerow([s, m.n, f"{m.median:.6f}", f"{m.q1:.6f}", f"{m.q3:.6f}",
                        f"{m.iqr:.6f}", f"{m.derived:.6f}", f"{m.ceiling:.6f}", m.clamped,
                        d["mult"], assay, lo, hi,
                        pp.get("floor_umi", ""), pp.get("n_all_with_a_value", "")])

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
    _require_gene_patterns("06_cluster_check", p)
    return _scanpy(pipeline, "cluster",
                   pipeline.results / "objects" / f"{p['sample']}_ambient.h5",
                   pipeline.scratch / f"{p['sample']}_clusters",
                   {"sample": p["sample"],
                    "resolution": p["resolution"], "seed": p["seed"],
                    # Needed to measure this object at all: it is the denoised one, so nothing
                    # has computed total_counts or pct_counts_mt on it, and cluster() refuses to
                    # normalise over an object whose depth was never measured.
                    "mt_prefix": p["mt_prefix"], "ribo_pattern": p["ribo_pattern"],
                    # The doublet calls are ATTACHED here and applied nowhere. Step 6's criterion
                    # D is only computable before a removal - afterwards every cluster is 0%
                    # doublet by construction - and this is the last point at which that holds.
                    "doublet_csv": p["doublet_csv"],
                    "doublet_key": "doublet_class",
                    "doublet_positive": "doublet"},
                   log, p["python_exe"])


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
                                            ("c_uninformative", "c_uninformative_pct"),
                                            ("d_doublet", "d_doublet_pct"))
                if d.get(y) is not None}
    if declared:
        base = {"a_umi_frac": proposed.a_umi_frac, "b_pct_mt": proposed.b_pct_mt,
                "c_uninformative": proposed.c_uninformative, "d_doublet": proposed.d_doublet}
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
        pipeline.findings.append({
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
                for s in task.params["samples"])]})

    fired, unknown = flagged.counts(), flagged.unknown_counts()
    for k in ("A", "B", "C", "D", "FLAG", "WATCH"):
        print(f"    {k:<6}fired {fired[k]:>4}   not evaluated {unknown[k]:>4}   "
              f"of {len(flagged.rows)}")
    # Both numbers, always. "How many fired" alone cannot distinguish a criterion that cleared
    # every cluster from one that was never computed on any of them.
    metrics = {"clusters": len(flagged.rows), "libraries": len(task.params["samples"]),
               "thresholds": str(thr)}
    metrics.update({f"{k}_fired": v for k, v in fired.items()})
    metrics.update({f"{k}_not_evaluated": v for k, v in unknown.items()})
    return {"outputs": [str(out)], "metrics": metrics, "versions": {}}


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
                    "doublet_csv": p["doublet_csv"],
                    "doublet_key": "doublet_class", "doublet_positive": "doublet"},
                   log, p["python_exe"])


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
    for f in ap.preflight(rows, kept_total=n_in - n_removed):
        pipeline.findings.append({"step": "07_apply", "check": f.check, "severity": f.severity,
                                  "message": f.message, "detail": list(f.detail or [])})

    # --- the ledger, BEFORE the gate. One row per removed barcode, listing every criterion that
    # fired on it, so what left can be named afterwards and re-read from the input.
    record = ap.build_removal_record(ids, masks)
    ledger = _tables(pipeline) / "removal_ledger.csv"
    ap.write_removal_record(record, ledger)

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
    res = _scanpy(pipeline, "apply_write", pipeline.results / "objects" / f"{samples[0]}_ambient.h5",
                  pipeline.results / "objects" / "cohort",
                  {"libraries": keep_lists}, log, task.params["python_exe"])

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
    return {"outputs": [str(p) for p in res.get("outputs", [])] + [str(ledger)],
            "metrics": {"n_in": n_in, "n_delivered": delivered, "n_removed": n_removed,
                        "action": action, "ceilings": ceiling_basis,
                        "authorised_by": authorised_by,
                        "per_sample_objects": written.get("per_sample_objects") or [],
                        "obs_columns": written.get("obs_columns") or [],
                        **{f"n_{c}": sum(1 for x in masks[c] if x) for c in criteria}},
            "versions": res.get("versions", {})}


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
    "_doublets", "_doublet_health", "_quality", "_quality_measure", "_quality_stage", "_cluster",
    "_cluster_flags", "_apply", "_apply_measure", "_ceilings_for", "_report",
]
