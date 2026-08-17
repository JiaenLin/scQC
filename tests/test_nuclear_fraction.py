"""WHY THIS SUITE EXISTS

Six defects, and four of them produce a filter that runs to completion and removes the wrong cells.

ONE — AN UNDEFINED FRACTION READ AS ZERO. A barcode with no assignable reads has no measurement.
Coerced to 0.0 it becomes the *most* nuclear-poor droplet in the library and is removed by any
floor at all, on evidence nobody collected. Section A.

TWO — THE JOIN SILENTLY MATCHING NOTHING. The object spells barcodes `<sample>_ACGT…`; the
aligner's file spells them `ACGT…`. Joined directly the intersection is exactly ZERO while both
files describe the same library — which looks like a library with no signal rather than like a
mismatch. Section B.

THREE — AGGREGATE ROWS SUMMED AS CELLS. STARsolo writes rows such as `CBnotInPasslist` for reads it
could not assign. Treated as barcodes they contribute one enormous "cell" to any pooled statistic.
Section B.

FOUR — THE BIGGEST LIBRARY DECIDING THE COHORT FLOOR. Pooling every cell lets a 16,000-cell library
outvote a 7,500-cell one, and the floor becomes a statement about how many nuclei each animal
yielded. The floor is the median of the per-library medians, equal weight per LIBRARY. Section C.

FIVE — THE CRITERION REPLACING THE FENCE INSTEAD OF ADDING TO IT. It must be a strict superset of
the run without it: `fail_mito_ceiling` keeps removing exactly what it removed. Section D.

SIX — THE FEATURE ARMING ITSELF. Absent a declared source nothing is measured, nothing fires, and
the run is what it was before this existed. Section E.

WHAT THIS SUITE DOES NOT DO: it chooses no floor, asserts no threshold for any cohort, and reads no
real data. Every number below is arithmetic built to make one branch distinguishable from another.
"""

from __future__ import annotations

import pathlib
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name if cond else f"{name}  [{detail}]")
    print(f"  {'ok    ' if cond else 'FAILED'} {name}" + (f"   {detail}" if not cond else ""))


from adapters import nuclear_fraction as nf  # noqa: E402

# ------------------------------------------------------- A. undefined is a THIRD state, not zero

check("the ratio is intronic / (intronic + exonic)", nf.nuclear_fraction(30, 70) == 0.3)
check("all intronic is 1.0", nf.nuclear_fraction(5, 0) == 1.0)
check("all exonic is 0.0", nf.nuclear_fraction(0, 5) == 0.0)
check("ZERO AND ZERO IS None, not 0.0", nf.nuclear_fraction(0, 0) is None,
      "a barcode with no assignable reads has no measurement; 0.0 would invent one")
check("a non-numeric value is None rather than an exception",
      nf.nuclear_fraction("x", 1) is None)
check("median ignores None instead of raising", nf.median([0.1, None, 0.3]) == 0.2)
check("a median of nothing is None", nf.median([None, None]) is None)

# ------------------------------------------------------------- B. the join, and what it must skip

HEADER = "CB\tcbMatch\texonic\tintronic\texonicAS\tintronicAS\tmito\n"
ROWS = [
    "AAA\t9\t70\t30\t1\t99\t2",       # sense 30/100 = 0.30; with AS 129/200 = 0.645
    "BBB\t9\t20\t80\t0\t0\t1",        # nf = 0.80
    "CCC\t9\t0\t0\t0\t0\t0",          # UNDEFINED
    "CBnotInPasslist\t9\t999999\t1\t0\t0\t0",   # aggregate row - never a cell
    "DDD\t9\t50\t50\t0\t0\t0",        # present in the file, NOT wanted by the caller
]


def _stats_file(tmp):
    p = pathlib.Path(tmp) / "CellReads.stats"
    p.write_text(HEADER + "\n".join(ROWS) + "\n", encoding="utf-8")
    return p


with tempfile.TemporaryDirectory() as tmp:
    src = _stats_file(tmp)
    targets = ["S1_AAA", "S1_BBB", "S1_CCC"]          # object spelling, with the sample prefix
    vals, st = nf.read_cellreads(src, targets, sample="S1")

    check("the `<sample>_` prefix is stripped for the join", set(vals) == set(targets),
          f"joined {sorted(vals)} - a direct join would intersect in ZERO barcodes")
    check("values are keyed as the OBJECT spells them", vals["S1_AAA"] == 0.3)
    check("...for every wanted barcode", vals["S1_BBB"] == 0.8)
    check("an undefined barcode is present with value None", "S1_CCC" in vals
          and vals["S1_CCC"] is None,
          "absent and undefined are different facts and must not collapse")
    check("aggregate rows are never returned as barcodes",
          not any(k.endswith("CBnotInPasslist") for k in vals))
    check("rows the caller did not ask for are not returned", "S1_DDD" not in vals
          and "DDD" not in vals,
          "the file holds one row per DROPLET; building all of them costs ~45x the memory")
    check("sense-only by default", st["antisense"] is False)
    check("the columns summed are recorded beside the values", "intronic" in st["columns"])
    check("the join rate is reported", st["n_target"] == 3 and st["n_joined"] == 3)
    check("...and how many were DEFINED, separately", st["n_defined"] == 2)

    va, sa = nf.read_cellreads(src, ["S1_AAA"], sample="S1", antisense=True)
    check("antisense changes the answer when asked for",
          abs(va["S1_AAA"] - 0.645) < 1e-9, f"got {va['S1_AAA']}, expected 129/200")
    check("...and says so, so the number is never silently a different quantity",
          sa["antisense"] is True and "intronicAS" in sa["columns"])

    # Column ORDER is not assumed - the header is authoritative.
    p2 = pathlib.Path(tmp) / "reordered.stats"
    p2.write_text("mito\tintronic\tCB\texonic\n5\t30\tAAA\t70\n", encoding="utf-8")
    v2, _ = nf.read_cellreads(p2, ["AAA"])
    check("columns are found by NAME, not by position", v2["AAA"] == 0.3)

    p3 = pathlib.Path(tmp) / "wrong.stats"
    p3.write_text("CB\tsomething\nAAA\t1\n", encoding="utf-8")
    try:
        nf.read_cellreads(p3, ["AAA"])
        check("a file without the needed columns is refused", False, "it did not refuse")
    except nf.NuclearFractionError as e:
        check("a file without the needed columns is refused, by name",
              "intronic" in str(e) or "exonic" in str(e))

# ------------------------------------------------- C. the cohort floor weights LIBRARIES equally

small = [0.10, 0.12, 0.14]                    # median 0.12
big = [0.80] * 10_000                         # median 0.80, ten thousand cells
per_library = [nf.median(small), nf.median(big), nf.median([0.20, 0.22, 0.24])]
floor_equal = nf.median(per_library)
pooled = nf.median(small + big + [0.20, 0.22, 0.24])
check("the cohort floor is the median of the per-library medians", floor_equal == 0.22,
      f"got {floor_equal}")
check("a ten-thousand-cell library cannot outvote a three-cell one", floor_equal != pooled,
      "pooling every cell makes the floor a statement about library size")
grown = nf.median([nf.median(small), nf.median([0.80] * 500_000),
                   nf.median([0.20, 0.22, 0.24])])
check("...and growing that library changes the equal-weight floor NOT AT ALL",
      grown == floor_equal)

# ---------------------------------------------------------- D. the criterion is ADDITIVE, and (E)

_ops = (ROOT / "adapters" / "scanpy_ops.py").read_text(encoding="utf-8")
check("the criterion is registered so the ledger and percell table carry it",
      '"fail_mito_nf"' in _ops.split("APPLY_CRITERIA = ", 1)[1].split(")", 1)[0])
_body = _ops.split("the joint mitochondrial x nuclear-fraction criterion", 1)[1].split(
    "removed = np.zeros", 1)[0]
check("fail_mito_ceiling is not reassigned by it - the removal stays a strict SUPERSET",
      "fail_mito =" not in _body,
      "an additive criterion that edits the fence can remove FEWER cells than the fence alone")
check("the trigger is a comparison, never a removal on its own",
      "m > nf_trigger" in _body and "and (v is not None" in _body)
check("an unknown fraction fails nothing", "v is not None and float(v) < nf_floor" in _body,
      "None must not compare as below the floor")
check("a floor outside (0,1) is refused as a percentage mistake",
      "must lie strictly in (0, 1)" in _body)

_steps = (ROOT / "engine" / "steps.py").read_text(encoding="utf-8")
check("E: the trigger is the DECLARED lower mitochondrial bound, not a new parameter",
      'out["metrics"].get("mito_bound_lo")' in _steps,
      "a second threshold would need its own justification and could drift from the ceiling's")

# --------------------------------------------------- F. the names are bound where they are USED
#
# THIS EXISTS BECAUSE A SOURCE GREP CANNOT CATCH IT. The floor derivation was first written into
# `_mito_ceiling_stage`, which holds the bound - while the per-library medians it reads are
# collected in `_quality_stage`, its caller. Every grep in section E passed, `py_compile` passed,
# and the run died with `NameError: name 'nf_per' is not defined` eighteen minutes in, after ten
# libraries had been measured and the doublet sweep had run.
#
# A whole-module linter would catch this, and there is none in this environment. So the check is
# scoped and dependency-free: every `nf_*` name LOADED inside a function must be STORED in that
# same function.
import ast  # noqa: E402

_tree = ast.parse(_steps)
for _fn in [n for n in ast.walk(_tree) if isinstance(n, ast.FunctionDef)]:
    loaded = {n.id for n in ast.walk(_fn)
              if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)
              and n.id.startswith("nf_")}
    stored = {n.id for n in ast.walk(_fn)
              if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store)}
    stored |= {a.arg for a in _fn.args.args}
    unbound = sorted(loaded - stored)
    if loaded:
        check(f"F: every nf_* name {_fn.name}() uses is bound in {_fn.name}()",
              not unbound, f"unbound: {unbound} - the block is in the wrong function")
check("E: a partial cohort is refused rather than filtered on other animals' boundary",
      "applied to all of them filters the rest on a boundary measured entirely in" in _steps)
check("E: the floor reaches step 7 through the barrier's metrics, like every other threshold",
      '_cohort_metric(pipeline, "nf_floor")' in _steps)
check("E: absent a declared source nothing is passed and nothing arms",
      'if p.get("cellreads_stats") else {}' in _steps)

_graph = (ROOT / "engine" / "graph.py").read_text(encoding="utf-8")
check("E: the source is per-library, from the samplesheet", '"cellreads_stats"' in _graph)

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
for f in FAIL:
    print(f"  FAILED: {f}")
sys.exit(1 if FAIL else 0)
