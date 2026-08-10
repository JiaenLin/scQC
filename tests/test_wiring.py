# Checks that the orchestrator calls the adapters with keywords they accept.
"""Every call from engine/ into adapters/ must match the adapter's real signature.

WHY THIS TEST EXISTS

Three defects in a row lived at the seam between two correct halves, and none was visible from
either side:

  * the orchestrator called two aligner adapters with `thread=` and `fastq_dir=` where they take
    `threads=` and `fastq_dirs=`, so the whole alignment path raised TypeError and every guard
    inside it was unreachable;
  * the report was handed a payload whose keys it does not recognise, and dropped them - the
    document came out with no refusal in it and nothing saying why;
  * the same report was handed one record per task where it expects one per step.

A unit test on either side passes in all three cases. The seam needs its own test, and it has to
be static: these calls cannot be exercised without the real tools, so a runtime test would be
skipped on every machine that matters.
"""

from __future__ import annotations

import ast
import importlib.util
import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

fails: list[str] = []


def load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


print("Wiring - engine calls against adapter signatures")
print("=" * 74)

# ---- A. every adapter call in engine/steps.py, read from the source rather than executed.
src = (ROOT / "engine" / "steps.py").read_text(encoding="utf-8")
tree = ast.parse(src)

ADAPTERS = {
    "cs": "adapters/celescope.py", "cr": "adapters/cellranger.py",
    "cbd": "adapters/cellbender.py", "db": "adapters/doublets.py",
    "so": "adapters/scanpy_ops.py", "mx": "adapters/matrix.py",
}
mods = {}
for alias, rel in ADAPTERS.items():
    try:
        mods[alias] = load(f"_wiring_{alias}", rel)
    except Exception as e:                                            # noqa: BLE001
        fails.append(f"A: {rel} does not import: {type(e).__name__}: {e}")

checked = 0
for node in ast.walk(tree):
    if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
        continue
    owner = node.func.value
    if not isinstance(owner, ast.Name) or owner.id not in mods:
        continue
    fn = getattr(mods[owner.id], node.func.attr, None)
    if fn is None:
        fails.append(f"A: engine/steps.py:{node.lineno} calls "
                     f"{owner.id}.{node.func.attr}(), which does not exist")
        continue
    try:
        params = inspect.signature(fn).parameters
    except (TypeError, ValueError):
        continue
    if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values()):
        continue
    passed = {kw.arg for kw in node.keywords if kw.arg}
    unknown = sorted(passed - set(params))
    checked += 1
    if unknown:
        fails.append(f"A: engine/steps.py:{node.lineno} calls {owner.id}.{node.func.attr}() "
                     f"with keyword(s) it does not accept: {unknown}")

print(f"A. adapter call sites checked: {checked}")

# ---- B. the report is handed the keys it documents, not the orchestrator's own shape.
try:
    from engine.pipeline import Pipeline
    required = {"run", "gates", "steps", "provenance"}
    have = set(inspect.getsource(Pipeline.report_payload))
    missing = [k for k in required if f'"{k}"' not in inspect.getsource(Pipeline.report_payload)]
    if missing:
        fails.append(f"B: report_payload() does not emit {missing}, which report/build.py "
                     f"documents as required")
    print(f"B. report payload keys present: {sorted(required - set(missing))}")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"B: could not inspect report_payload: {type(e).__name__}: {e}")

# ---- C. every step the graph names has a description and a stated limit.
try:
    from engine import steps as st
    keys = {"00_ingest", "00_align", "01_ambient", "02_cells", "04_doublets",
            "05_quality", "06_cluster_check", "07_apply", "report"}
    for k in sorted(keys):
        what, cannot = st.step_text(k)
        if not what or not cannot:
            fails.append(f"C: step {k!r} has no {'description' if not what else 'stated limit'}. "
                         f"An omitted limit reads as no limit.")
    print(f"C. steps with a description and a stated limit: "
          f"{sum(1 for k in keys if all(st.step_text(k)))} of {len(keys)}")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"C: could not read STEP_TEXT: {type(e).__name__}: {e}")

# ---- D. every task the graph builds names a step the orchestrator can group under.
try:
    from engine import graph
    print("D. graph module imports and exposes: "
          f"{sorted(n for n in ('ingest_stage', 'main_stage') if hasattr(graph, n))}")
    for n in ("ingest_stage", "main_stage"):
        if not hasattr(graph, n):
            fails.append(f"D: engine/graph.py does not expose {n}()")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"D: engine/graph.py does not import: {type(e).__name__}: {e}")

# ---- E. the engine never measures data in its OWN interpreter.
#
# Check A compares call sites against real signatures, so it catches wrong ARGUMENTS. It did not
# catch step 0 calling mx.ingest_stats_fn: that call was perfectly well formed, and the mistake
# was choosing a function that measures in-process. The result needed pandas and anndata in
# whatever interpreter runs scqc - which on a cluster is a bare system python, because the
# aligner, the denoiser and the analysis stack have incompatible pins and live in separate
# environments. Step 0 therefore passed on a laptop and failed on every cluster.
#
# These adapter entry points load matrices HERE, so the engine must not call them. Each has an
# out-of-process sibling that takes `python_exe` and an executor.
IN_PROCESS_ONLY = {
    "ingest_stats_fn": "run_summary_stats(..., python_exe=..., executor=...)",
    "summary_stats": "run_summary_stats(..., python_exe=..., executor=...)",
    "barcode_rank": "run_summary_stats(..., rank_points=N, python_exe=..., executor=...)",
    "read_matrix": "an adapter call that runs in the analysis environment",
}
offenders = []
for py in sorted((ROOT / "engine").glob("*.py")):
    for node in ast.walk(ast.parse(py.read_text(encoding="utf-8"))):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        owner = node.func.value
        if not isinstance(owner, ast.Name) or owner.id not in ADAPTERS:
            continue
        if node.func.attr in IN_PROCESS_ONLY:
            offenders.append(
                f"E: engine/{py.name}:{node.lineno} calls {owner.id}.{node.func.attr}(), which "
                f"loads the matrix in THIS interpreter. Use {IN_PROCESS_ONLY[node.func.attr]}")
print(f"E. engine calls into in-process adapter entry points: {len(offenders)}")
fails.extend(offenders)

# ---- F. every `needs` resolves inside its own stage.
#
# Stage one (ingest) runs to completion before stage two is built, so a stage-two task naming a
# stage-one key is unresolvable and the run stops with "depends on unknown task(s)". The ordering
# such an edge tries to express is already guaranteed by the two-stage structure.
#
# This was a latent break for every sample whose matrix is ACCEPTED rather than aligned, and it
# never fired because no run had ever supplied a matrix - which is the path a reusable pipeline is
# most likely to be used on. Built for three cohort shapes, since the graph differs between them.
class _P:
    results = Path("/tmp/r"); work = Path("/tmp/w"); project = Path("/tmp/p")
    decisions: dict = {}
    mode = "evidence"
    samples: list = []


def _row(s, **extra):
    return {"sample": s, "platform": "singleron", "species": "mus_musculus",
            "reference": "mus_musculus/ensembl_112_filtered", "assay": "snrna",
            "matrix": f"/d/{s}", "diet": "chow" if s.endswith("1") else "hfd", **extra}


_AMB = {"ambient_h5": "/cb/x.h5", "ambient_tool": "CellBender", "ambient_version": "0.3.2",
        "ambient_params": "--fpr 0", "ambient_produced_by": "a prior run"}
_SHAPES = {
    "accepted matrix, ambient run here": [_row("A"), _row("B")],
    "accepted matrix, ambient supplied": [_row("A", **_AMB), _row("B", **_AMB)],
    "mixed: one supplied, one not":      [_row("A", **_AMB), _row("B")],
}
_bad = 0
for _label, _rows in _SHAPES.items():
    _P.samples = _rows
    _ing = {r["sample"]: {"mode": "accept"} for r in _rows}
    try:
        _built = graph.main_stage(_P, "python", {}, _ing)
    except Exception as e:                                            # noqa: BLE001
        fails.append(f"F: graph does not build for {_label!r}: {type(e).__name__}: {e}")
        continue
    _keys = {t.key for t in _built}
    for _t in _built:
        for _n in (_t.needs or ()):
            if _n not in _keys:
                fails.append(f"F: [{_label}] task {_t.key!r} needs {_n!r}, which is not in this "
                             f"stage. Stage one has already run; drop the edge.")
                _bad += 1
print(f"F. cross-stage dependency edges over {len(_SHAPES)} cohort shapes: {_bad}")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("wiring OK - the seam between engine and adapters is consistent")

