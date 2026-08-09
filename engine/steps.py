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
    stats_fn = mx.ingest_stats_fn(
        expected_genes=task.params.get("expected_genes"),
        tmp_dir=pipeline.work / f"{sample}_extract")

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


def _ambient_audit(task, pipeline, log):
    from adapters import cellbender as cbd

    aa = step_module("audit_ambient")
    rows, per_gene = [], []
    for s, r in task.params["per_sample"].items():
        m = cbd.parse_metrics(r["h5"], r["raw"])
        rows.append({"sample": s,
                     "fraction_removed_overall": m["fraction_removed_overall"],
                     "genes_fully_removed": m["genes_fully_removed"]})
        for g in m.get("per_gene", []):
            per_gene.append({"sample": s, **g})

    import pandas as pd
    summ = pd.DataFrame(rows)
    genes = pd.DataFrame(per_gene) if per_gene else None
    out = _tables(pipeline) / "ambient_summary.csv"
    summ.to_csv(out, index=False)

    findings = aa.audit(summ, genes, task.params["design"])
    pipeline.gate("01_ambient", findings, aa.verdict(findings))
    return {"outputs": [str(out)], "metrics": {"libraries": len(rows)}, "versions": {}}


# --------------------------------------------------------------------------------------------
# step 2 - the cell-call gate


def _cellcall(task, pipeline, log):
    cg = step_module("cellcall_gate")
    calls = task.params["calls"]
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

    payload = pipeline.payload(stopped=None)
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


__all__ = ["build_tasks", "_design", "_ingest", "_ambient", "_ambient_audit",
           "_cellcall", "_quality", "_report"]
