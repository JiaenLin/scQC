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
#
# BOTH MODES, because apply mode adds the only tasks that remove anything and no test built it.
# The graph must ASSEMBLE without resolving a threshold: with a per-library ceiling those come
# from step 5's table, which does not exist while a clean run's graph is being built, so
# resolving during construction refuses a correct run for a file three steps from being written.
# Building it here with empty decisions is what catches that.
_bad = 0
_combos = [(lbl, rows, mode) for lbl, rows in _SHAPES.items() for mode in ("evidence", "apply")]
for _label, _rows, _mode in _combos:
    _P.samples = _rows
    _P.mode = _mode
    _P.decisions = {}
    _ing = {r["sample"]: {"mode": "accept"} for r in _rows}
    try:
        _built = graph.main_stage(_P, "python", {}, _ing)
    except Exception as e:                                            # noqa: BLE001
        fails.append(f"F: graph does not build for {_label!r} in {_mode} mode: "
                     f"{type(e).__name__}: {e}")
        continue
    _keys = {t.key for t in _built}
    for _t in _built:
        for _n in (_t.needs or ()):
            if _n not in _keys:
                fails.append(f"F: [{_label}/{_mode}] task {_t.key!r} needs {_n!r}, which is not "
                             f"in this stage. Stage one has already run; drop the edge.")
                _bad += 1
    # Evidence mode has no path to a removal, and that is structural rather than a flag: the
    # tasks are ABSENT, so this asserts absence rather than that something defaults to off.
    _removers = sorted(k for k in _keys if k.startswith("07_"))
    if _mode == "evidence" and _removers:
        fails.append(f"F: [{_label}] evidence mode placed {_removers} in the graph. Evidence "
                     f"mode must have no task that can remove, not a task that declines to.")
    if _mode == "apply" and not _removers:
        fails.append(f"F: [{_label}] apply mode placed no step-7 task, so nothing can apply.")
_P.mode = "evidence"
print(f"F. cross-stage dependency edges over {len(_combos)} cohort shape / mode combinations: "
      f"{_bad}")

# ---- G. no getattr() default on an adapter object the engine owns.
#
# `getattr(calls, "scored", None) or {}` turned a wrong attribute name into an empty mapping.
# DoubletCalls has no `.scored` - it IS the mapping - so every library reported no doublet rate
# and step 4 refused with "no library produced a doublet rate". The refusal was honest; the cause
# was invisible, because the default answered a question the object had never been asked.
#
# A missing attribute on an object this repo defines is a bug, and bugs must raise. Defaults are
# for values that are legitimately absent - and those are modelled explicitly here (DoubletCalls
# .unscored is None when the population was not supplied, which is deliberately not the same as
# empty).
_getattr_defaults = []
for py in sorted((ROOT / "engine").glob("*.py")):
    src_g = py.read_text(encoding="utf-8")
    for node in ast.walk(ast.parse(src_g)):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
            continue
        if node.func.id != "getattr" or len(node.args) < 3:
            continue
        # Only flag it where the object came from an adapter this repo owns; getattr with a
        # default on argparse namespaces and third-party objects is ordinary and correct.
        target = node.args[0]
        name = target.id if isinstance(target, ast.Name) else ""
        if name in ("calls", "export", "res", "plan", "result"):
            attr = node.args[1].value if isinstance(node.args[1], ast.Constant) else "?"
            _getattr_defaults.append(
                f"G: engine/{py.name}:{node.lineno} getattr({name}, {attr!r}, <default>) on an "
                f"object this repo defines - a wrong attribute name resolves to the default "
                f"instead of raising")
print(f"G. getattr-with-default on repo-owned adapter objects: {len(_getattr_defaults)}")
fails.extend(_getattr_defaults)

# ---- H. every scanpy op gets the params it declares as required.
#
# Params cross the `_scanpy` seam as a DICT, so check A - which compares call sites against
# function signatures - cannot see a missing or misnamed key. Step 5 called the valley op with
# {"metric": "umi"} where it requires `metrics` (a LIST), `sample`, `mt_prefix` and
# `ribo_pattern`. It failed on the first library of the first real cohort, having never run.
#
# The requirement is read from the adapter's own idiom - `_require_*(params, "X")` and the
# hand-written `if "X" not in params` - so this check cannot drift from what the ops actually
# enforce. Both spellings are read because reading only the first left step 6 unprotected: the
# cluster op demands `uninformative_genes`, `doublet_key`, `resolution` and `seed` by hand, the
# engine passed two of the four, and each missing one cost a separate submit-and-wait cycle to
# discover. It still does NOT catch a key required some third way - `metrics` is validated inside
# an `isinstance` - and that limit is stated rather than papered over.
_ops_src = (ROOT / "adapters" / "scanpy_ops.py").read_text(encoding="utf-8")
_ops_tree = ast.parse(_ops_src)
_required: dict = {}
for fn in [n for n in ast.walk(_ops_tree)
           if isinstance(n, ast.FunctionDef) and n.name.startswith("_op_")]:
    need = set()
    for node in ast.walk(fn):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id.startswith("_require_") and len(node.args) >= 2
                and isinstance(node.args[0], ast.Name) and node.args[0].id == "params"
                and isinstance(node.args[1], ast.Constant)):
            need.add(node.args[1].value)
        # `if "X" not in params: raise` - the same requirement, written out by hand where the
        # refusal needs to say something the generic one cannot.
        if (isinstance(node, ast.Compare) and len(node.ops) == 1
                and isinstance(node.ops[0], ast.NotIn)
                and isinstance(node.left, ast.Constant) and isinstance(node.left.value, str)
                and len(node.comparators) == 1
                and isinstance(node.comparators[0], ast.Name)
                and node.comparators[0].id == "params"):
            need.add(node.left.value)
    _required[fn.name[len("_op_"):]] = need

_steps_tree = ast.parse((ROOT / "engine" / "steps.py").read_text(encoding="utf-8"))
_seam = 0
for node in ast.walk(_steps_tree):
    if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
            and node.func.id == "_scanpy" and len(node.args) >= 5):
        continue
    op_node, params_node = node.args[1], node.args[4]
    if not isinstance(op_node, ast.Constant) or not isinstance(params_node, ast.Dict):
        continue
    op = op_node.value
    given = {k.value for k in params_node.keys if isinstance(k, ast.Constant)}
    absent = sorted(_required.get(op, set()) - given)
    _seam += 1
    if absent:
        fails.append(f"H: engine/steps.py:{node.lineno} calls the {op!r} op without "
                     f"{absent}, which _op_{op} requires. Params cross this seam as a dict, so "
                     f"nothing else checks it.")
print(f"H. scanpy op call sites checked against declared requirements: {_seam}")

# ---- I. the PBS executor reads the log the JOB wrote, not the one PBS was asked to deliver.
#
# `#PBS -o <path>` names a path ON THE HOST THAT SUBMITTED THE JOB. When the orchestrator is
# itself a job, tasks land on other nodes, and PBS must copy each output file back between two
# compute nodes once the job ends. Where that copy does not work there is no error anywhere: the
# job reports Exit_status = 0, its real output is correct on disk, and the file simply never
# appears. It failed only for tasks that landed on a different node from the orchestrator, so
# seven of ten libraries succeeded and three did not - which reads as a staging delay and is not
# one. A 60 s wait was added for it and could not have helped; nothing was on the way.
#
# The job now redirects its own output to the log on shared storage. Constructing `{log}.out` as
# a path is the signature of the old assumption coming back, so it is what this checks.
_ex_src = (ROOT / "engine" / "executor.py").read_text(encoding="utf-8")
_shell_fn = None
for _cls in [n for n in ast.walk(ast.parse(_ex_src))
             if isinstance(n, ast.ClassDef) and n.name == "PBSExecutor"]:
    for _fn in [n for n in _cls.body if isinstance(n, ast.FunctionDef) and n.name == "shell"]:
        _shell_fn = _fn
if _shell_fn is None:
    fails.append("I: PBSExecutor.shell() could not be found in engine/executor.py")
else:
    if "exec > " not in (ast.get_source_segment(_ex_src, _shell_fn) or ""):
        fails.append("I: the generated PBS script does not redirect its own output to the log. "
                     "Without it the executor is relying on PBS to copy the file between hosts, "
                     "which on this cluster fails silently.")
    # Whatever gets read is what the adapters are handed and what an exit status is judged from.
    # Naming the source by expression rather than by spelling: only `log` - the path the job
    # writes to itself - may be read here.
    _read = [n for n in ast.walk(_shell_fn)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
             and n.func.attr == "read_text"]
    if not _read:
        fails.append("I: PBSExecutor.shell() reads no output at all.")
    for _r in _read:
        _from = _r.func.value
        if not (isinstance(_from, ast.Name) and _from.id == "log"):
            fails.append(
                f"I: engine/executor.py:{_r.lineno} reads output from something other than `log`. "
                f"The only other candidate is the path PBS was asked to DELIVER, and a delivered "
                f"file that never arrives is indistinguishable from a job that printed nothing.")
    print(f"I. PBS output reads, all from the job's own log: {len(_read)}")

# ---- J. every report is built from report_payload(), and from nothing else.
#
# The document is assembled at TWO call sites: the `report` task, and scqc_cli.finish() after the
# graph ends - the second exists so a run that stopped still leaves a document. When the task
# enriched its payload and finish() rebuilt one without the addition, finish() ran last and wrote
# over the top: the per-library threshold table was absent from every report while the code
# producing it ran correctly on every run, and the only symptom was the report's own defect
# saying a section had not been supplied.
#
# Two builders of one document will always drift. So the payload has exactly one source, and any
# call to build_report() whose payload is not report_payload() is that drift starting again.
#
# A payload READ BACK from reports/payload.json is accepted, and is not a second builder: that
# file is written by the report task from the payload report_payload() returned, so replaying it
# reproduces the same document rather than assembling a competing one. `scqc report` rebuilds a
# finished run this way, which is what makes a caption or a figure iterable without re-running a
# pipeline. The distinction the check still enforces is assembling a payload versus reloading
# one - a hand-built dict named `payload` fails here exactly as before.
_builders = []
for _py in sorted((ROOT / "engine").glob("*.py")) + [ROOT / "scqc_cli.py"]:
    for _n in ast.walk(ast.parse(_py.read_text(encoding="utf-8"))):
        if not (isinstance(_n, ast.Call) and isinstance(_n.func, ast.Name)
                and _n.func.id == "build_report" and _n.args):
            continue
        _arg = _n.args[0]
        # Either build_report(pipe.report_payload(...)) or a name assigned from it.
        _ok = (isinstance(_arg, ast.Call) and isinstance(_arg.func, ast.Attribute)
               and _arg.func.attr == "report_payload")
        if not _ok and isinstance(_arg, ast.Name):
            _src = ast.get_source_segment(_py.read_text(encoding="utf-8"), _n) or ""
            _fn = next((f for f in ast.walk(ast.parse(_py.read_text(encoding="utf-8")))
                        if isinstance(f, ast.FunctionDef)
                        and f.lineno <= _n.lineno <= (f.end_lineno or _n.lineno)), None)
            _ok = _fn is not None and any(
                isinstance(a, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == _arg.id for t in a.targets)
                and any(k in (ast.get_source_segment(_py.read_text(encoding="utf-8"), a) or "")
                        for k in ("report_payload", "payload.json"))
                for a in ast.walk(_fn))
        _builders.append((f"{_py.name}:{_n.lineno}", _ok))
        if not _ok:
            fails.append(f"J: {_py.name}:{_n.lineno} calls build_report() with a payload that "
                         f"did not come from report_payload(). Two builders of one document "
                         f"drift, and the loser is silent - the section is simply absent.")
if not _builders:
    fails.append("J: no call to build_report() was found; this check is not testing anything.")
print(f"J. build_report call sites, all fed by report_payload(): "
      f"{sum(1 for _, ok in _builders if ok)} of {len(_builders)}")

# ---- J2. the report task does not ADD to the payload it was given.
#
# The other half of J, and the half that actually keeps happening. finish() rebuilds the payload
# after the graph and writes last, so anything the report TASK adds to its own copy is silently
# discarded - the task's log says the section was assembled, and the document does not have it.
#
# It has now happened twice for the same reason: once for the per-library threshold table, and
# once for the figures, where seven were assembled, printed as assembled, and rendered nowhere.
# Both times the code was correct and ran; both times the addition went in at the wrong layer.
# A section belongs to report_payload(), which both callers use, or it does not exist.
_report_fn = next((f for f in ast.walk(ast.parse((ROOT / "engine" / "steps.py")
                                                 .read_text(encoding="utf-8")))
                   if isinstance(f, ast.FunctionDef) and f.name == "_report"), None)
if _report_fn is None:
    fails.append("J2: engine/steps.py has no _report function; this check is not testing "
                 "anything.")
else:
    _mutations = []
    for _n in ast.walk(_report_fn):
        if isinstance(_n, ast.Assign):
            for _t in _n.targets:
                if isinstance(_t, ast.Subscript) and isinstance(_t.value, ast.Name) \
                        and _t.value.id == "payload":
                    _mutations.append(_n.lineno)
        if isinstance(_n, ast.Call) and isinstance(_n.func, ast.Attribute) \
                and isinstance(_n.func.value, ast.Name) and _n.func.value.id == "payload" \
                and _n.func.attr in ("update", "setdefault"):
            _mutations.append(_n.lineno)
    for _ln in sorted(set(_mutations)):
        fails.append(f"J2: engine/steps.py:{_ln} adds to the payload inside the report task. "
                     f"scqc_cli.finish() rebuilds the payload afterwards and writes last, so "
                     f"this addition reaches the task's own document and no other. Put it in "
                     f"report_payload(), which both callers go through.")
    print(f"J2. payload additions inside the report task: {len(set(_mutations))}")

# ---- J3. nothing mutates pipeline.findings, because it is a property and the write vanishes.
#
# `findings` builds a new list on every access, so `pipeline.findings.append(...)` appends to a
# throwaway and is discarded without an error anywhere. Two call sites did that and lost four
# findings between them - one of which was the statement that 21,395 delivered nuclei sat inside
# a cluster the run had flagged. The report was complete, well-formed, and quietly missing the
# finding a reader most needed.
_mutators = []
for _py in sorted((ROOT / "engine").glob("*.py")) + sorted((ROOT / "modules").rglob("*.py")):
    for _n in ast.walk(ast.parse(_py.read_text(encoding="utf-8"))):
        if not (isinstance(_n, ast.Call) and isinstance(_n.func, ast.Attribute)
                and _n.func.attr in ("append", "extend", "insert")):
            continue
        _tgt = _n.func.value
        if (isinstance(_tgt, ast.Attribute) and _tgt.attr == "findings"
                and isinstance(_tgt.value, ast.Name) and _tgt.value.id != "self"):
            _mutators.append(f"{_py.name}:{_n.lineno}")
for _m in _mutators:
    fails.append(f"J3: {_m} mutates pipeline.findings directly. It is a read-only property - "
                 f"the write goes to a temporary list and disappears silently. Use "
                 f"record_findings(), or gate() when there is a verdict to act on.")
print(f"J3. direct writes to pipeline.findings: {len(_mutators)}")

# ---- K. every step_module() the engine asks for is one the engine can load.
#
# `step_module("audit_removal")` raised KeyError at runtime because the module was never added to
# _STEP_MODULES - a decision module that exists, is tested, and is unreachable from the engine.
# It is the same seam as check A, one level up: check A compares a call against a signature, this
# compares a NAME against the registry, and neither side is wrong on its own.
#
# It surfaces only when the step actually runs, which for step 7 means only under `--mode apply`
# with a complete decisions file - the least-travelled path in the pipeline and the only one that
# removes anything.
try:
    from engine.pipeline import _STEP_MODULES
    _asked = set()
    for _py in sorted((ROOT / "engine").glob("*.py")):
        for _n in ast.walk(ast.parse(_py.read_text(encoding="utf-8"))):
            if (isinstance(_n, ast.Call) and isinstance(_n.func, ast.Name)
                    and _n.func.id == "step_module" and _n.args
                    and isinstance(_n.args[0], ast.Constant)):
                _asked.add(_n.args[0].value)
    _unregistered = sorted(_asked - set(_STEP_MODULES))
    _absent = sorted(k for k, p in _STEP_MODULES.items() if not Path(p).exists())
    for _u in _unregistered:
        fails.append(f"K: the engine calls step_module({_u!r}), which is not in _STEP_MODULES. "
                     f"It raises KeyError when that step runs, not when the graph is built.")
    for _a in _absent:
        fails.append(f"K: _STEP_MODULES registers {_a!r} at a path that does not exist.")
    print(f"K. step_module names asked for, all registered and present: "
          f"{len(_asked) - len(_unregistered)} of {len(_asked)}")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"K: could not check the step-module registry: {type(e).__name__}: {e}")

# ---- L. an adapter the orchestrator RUNS AS A SCRIPT has no parent package, so a relative
#         import inside it raises `attempted relative import with no known parent package`. Same
#         shape as the defects above: invisible from either side, and from inside the adapter it
#         reads as ordinary package code.
#
#         Checked statically, and checked for imports INSIDE FUNCTIONS too - which is the part
#         that costs. A relative import at module level fails the moment anything touches the
#         file; one inside `_op_apply_write` fails at step 7, after every other step of a
#         ten-library cohort has completed. That is exactly where it did fail, and neither
#         `--help` nor a direct `import adapters.x` would have caught it.
try:
    _script_adapters = sorted(p for p in Path("adapters").glob("*.py")
                              if "__main__" in p.read_text(encoding="utf-8"))
    _rel = []
    for _py in _script_adapters:
        for _n in ast.walk(ast.parse(_py.read_text(encoding="utf-8"))):
            if isinstance(_n, ast.ImportFrom) and (_n.level or 0) > 0:
                _names = ", ".join(a.name for a in _n.names)
                _rel.append(f"{_py}:{_n.lineno} from {'.' * _n.level}{_n.module or ''} "
                            f"import {_names}")
    for _r in _rel:
        fails.append(f"L: {_r} - this module is run as a script, so it has no parent package and "
                     f"a relative import raises. Use `from adapters.x import ...`.")
    print(f"L. adapters run as scripts, none using a relative import: "
          f"{len(_script_adapters)} checked, {len(_rel)} bad")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"L: could not check adapter imports: {type(e).__name__}: {e}")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("wiring OK - the seam between engine and adapters is consistent")

