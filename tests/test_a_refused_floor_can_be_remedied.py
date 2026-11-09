# Calls step 5's barrier and the threshold resolver of steps 6 and 7 on a Pipeline in a temporary
# project, with the per-library measurements put in place by hand; reads no matrix.
"""A floor whose derivation refuses can still be decided, by a person (harness autodev W020, B1-042).

Step 5 refused - no library bimodal, a valley outside the bounds - BEFORE the decisions file was
read, so the remedy its own message named ("choose the cut explicitly and record it") could not be
applied: the project's rule "never block the remedy" (docs/PRINCIPLES.md section 3), broken by its
own step. Checked here, as the register's CHECK states it:

  1  every library's UMI valley a shoulder, and quality.umi_floor = 400 attested in decisions.yml:
     step 5 completes, publishes 400 as ADJUDICATED with a REVIEW finding, and steps 6 and 7 apply
     400 as ADJUDICATED
  2  the same without it: step 5 completes with a REVIEW finding and no UMI floor; step 6's resolver
     goes on without one (evidence completes); step 7's refuses, naming quality.umi_floor
  3  a derivation that succeeds is untouched by the decisions file
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from engine import steps  # noqa: E402
from engine.executor import make_executor  # noqa: E402
from engine.pipeline import Pipeline  # noqa: E402
from engine.task import Refusal, Status, Task, TaskResult  # noqa: E402

fails: list[str] = []
SAMPLES = ["N1", "N2", "N3", "N4"]
ATTESTED = {"quality": {"umi_floor": {"value": 400, "approved_by": "the analyst",
                                      "verbatim": "the shoulders sit near 400, so the cut goes there"}}}


def barrier(td: Path, decisions, umi_bimodal: bool):
    pipe = Pipeline(project=td / "proj", mode="apply", executor=make_executor("local"),
                    samples=[{"sample": s, "assay": "snrna"} for s in SAMPLES],
                    decisions=decisions, jobs=1)
    for i, s in enumerate(SAMPLES):
        pipe.results_by_key[f"05_quality/{s}"] = TaskResult(
            f"05_quality/{s}", Status.DONE, step="05_quality", sample=s,
            metrics={"valleys": {"umi": 300.0 + 20 * i, "genes": 220.0 + 10 * i},
                     "bimodal": {"umi": umi_bimodal, "genes": True},
                     "mito_quartiles": None, "mito_population": {}})
    # the ceilings are not this test's subject: step 5's ceiling stage is left out, and the table
    # steps 6 and 7 read the ceilings from is written by hand
    real = steps._mito_ceiling_stage
    steps._mito_ceiling_stage = lambda task, pipeline, mito_stats, out, mito_pop=None: out
    try:
        t = Task(key="05_quality", step="05_quality", fn=steps._quality_stage,
                 params={"samples": SAMPLES, "light_floor": 150, "decisions": decisions,
                         "assay": {s: "snrna" for s in SAMPLES},
                         "call_by": {s: "denoiser" for s in SAMPLES}})
        out = steps._quality_stage(t, pipe, None)
    finally:
        steps._mito_ceiling_stage = real
    pipe.results_by_key["05_quality"] = TaskResult("05_quality", Status.DONE, step="05_quality",
                                                   metrics=out["metrics"])
    (pipe.results / "tables" / "mito_ceiling_per_sample.csv").write_text(
        "sample,ceiling\n" + "".join(f"{s},8.0\n" for s in SAMPLES))
    return pipe, out


# ---- 1. attested
with tempfile.TemporaryDirectory() as td:
    try:
        pipe, out = barrier(Path(td), ATTESTED, umi_bimodal=False)
        m = out["metrics"]
        if (m.get("umi_proposed"), m.get("umi_provenance")) != (400, "adjudicated"):
            fails.append(f"1: step 5 published umi {m.get('umi_proposed')!r} as "
                         f"{m.get('umi_provenance')!r}, not 400 adjudicated")
        if not any(f.get("check") == "umi floor: derivation refused"
                   and f.get("severity") == "REVIEW" for f in pipe.findings):
            fails.append("1: the refused derivation was not recorded as a REVIEW finding")
        for where in ("06_cluster_check (N1)", "07_apply"):
            v, c, _ce, _b = steps._resolve_thresholds(pipe, ATTESTED, SAMPLES, step=where,
                                                      allow_refused=where.startswith("06"))
            if (v["quality.umi_floor"], c["quality.umi_floor"]) != (400, "ADJUDICATED"):
                fails.append(f"1: {where} resolves the UMI floor as {v['quality.umi_floor']!r} "
                             f"{c['quality.umi_floor']!r}")
    except Refusal as e:
        fails.append(f"1: step 5 refused although a floor was attested: {str(e)[:200]}")
print("1. all shoulders + an attested 400: ADJUDICATED 400 through steps 5, 6 and 7")

# ---- 2. not attested
with tempfile.TemporaryDirectory() as td:
    try:
        pipe, out = barrier(Path(td), {}, umi_bimodal=False)
        m = out["metrics"]
        if (m.get("umi_proposed"), m.get("umi_provenance")) != (None, "refused"):
            fails.append(f"2: step 5 published umi {m.get('umi_proposed')!r} "
                         f"{m.get('umi_provenance')!r}, not None refused")
        pop = steps._population_for_cluster(pipe, "N1", "cellbender_cell", {})
        if pop.get("umi_floor") is not None or pop.get("gene_floor") is None:
            fails.append(f"2: step 6's population is {pop}")
        try:
            steps._resolve_thresholds(pipe, {}, SAMPLES, step="07_apply")
            fails.append("2: step 7 applied a floor that was neither declared nor derived")
        except Refusal as e:
            if "quality.umi_floor" not in str(e) or "decisions.yml" not in str(e):
                fails.append(f"2: step 7's refusal does not name the remedy: {str(e)[:200]}")
    except Refusal as e:
        fails.append(f"2: step 5 still refuses before its remedy can be read: {str(e)[:200]}")
print("2. all shoulders, nothing attested: steps 5 and 6 go on, step 7 refuses naming the remedy")

# ---- 3. a derivation that succeeds is untouched
with tempfile.TemporaryDirectory() as td:
    pipe, out = barrier(Path(td), ATTESTED, umi_bimodal=True)
    m = out["metrics"]
    if m.get("umi_provenance") == "adjudicated" or m.get("umi_proposed") == 400:
        fails.append(f"3: a successful derivation was replaced: {m.get('umi_proposed')!r}")
print("3. a derivation that succeeds is what step 5 publishes")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("a refused floor no longer blocks its remedy")
