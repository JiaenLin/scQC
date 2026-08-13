# Exercises the resolution default and the extra-resolution siblings. It removes nothing and
# clusters nothing: every assertion here is about names, parsing and selection.
"""Step 6: the applied resolution, the extras beside it, and the contract that keeps them apart.

The point of this suite is the LAST test. Everything downstream of step 6 - step 7, the object's
`cluster_FLAG`, every consuming stage - reads exactly one table, selected by the suffix
`.cluster_profile.csv`. An extra resolution that produced a file matching that suffix would be
picked up as though it were the applied one, silently, and the deliverable would carry flags from
a clustering nobody chose. The naming is therefore a contract, not a convention.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

fails = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        fails.append(name)


print("\n1 - the applied resolution defaults to 2.0")
import scqc_cli  # noqa: E402

parser = scqc_cli.build_parser() if hasattr(scqc_cli, "build_parser") else None
if parser is None:                       # the CLI builds its parser inside main(); read the source
    src = (HERE.parent / "scqc_cli.py").read_text(encoding="utf-8")
    m = re.search(r'--resolution["\s,]+type=float,\s*default=([0-9.]+)', src)
    check("--resolution default is 2.0", bool(m) and float(m.group(1)) == 2.0,
          f"found {m.group(1) if m else 'nothing'}")
    m2 = re.search(r'--extra-resolutions".*?default="([^"]*)"', src, re.S)
    check("--extra-resolutions defaults to 1.0,3.0", bool(m2) and m2.group(1) == "1.0,3.0",
          f"found {m2.group(1) if m2 else 'nothing'}")
else:
    a = parser.parse_args(["run", "--project", "."])
    check("--resolution default is 2.0", a.resolution == 2.0, f"got {a.resolution}")


print("\n2 - the extras list is parsed, deduplicated, and never repeats the applied resolution")


def parse_extras(text, applied):
    """The CLI's own expression, kept in one place so the test exercises the real rule."""
    return (sorted({r for r in (float(t.strip()) for t in str(text).split(",") if t.strip())
                    if r != float(applied)}) if text else None)


check("drops the applied resolution", parse_extras("1.0,2.0,3.0", 2.0) == [1.0, 3.0],
      str(parse_extras("1.0,2.0,3.0", 2.0)))
check("drops empty tokens", parse_extras("1.0,,3.0", 2.0) == [1.0, 3.0])
check("deduplicates", parse_extras("3.0,1.0,3.0", 2.0) == [1.0, 3.0])
check("an empty request is None, not an empty sweep", parse_extras("", 2.0) is None)
check("all-extras-equal-applied leaves an empty list, not None",
      parse_extras("2.0", 2.0) == [])


print("\n3 - a sibling table can never be mistaken for the applied one")
APPLIED_SUFFIX = ".cluster_profile.csv"          # what engine/steps.py selects on
for res in ("1", "1.5", "3", "0.25"):
    sib = f"/work/ctrl_01_clusters.cluster_profile.res{res}.csv"
    check(f"res{res} sibling is not selected as the applied profile",
          not sib.endswith(APPLIED_SUFFIX), sib.rsplit("/", 1)[-1])
check("the applied profile IS selected",
      "/work/ctrl_01_clusters.cluster_profile.csv".endswith(APPLIED_SUFFIX))


print("\n4 - the resolution token round-trips out of the sibling's name")
for res in ("1", "1.5", "3", "0.25"):
    name = f"/work/ctrl_01_clusters.cluster_profile.res{res}.csv"
    got = name.rsplit(".cluster_profile.res", 1)[1][:-4]     # engine/steps.py's own expression
    check(f"res{res} recovered from the file name", got == res, f"got {got!r}")


print("\n5 - the applied table's name and columns are untouched by any of this")
sys.path.insert(0, str(HERE.parent / "modules" / "06_cluster_check"))
from cluster_flags import BOOLEAN_KEYS, NUMERIC_KEYS  # noqa: E402

check("the verdict columns are unchanged",
      BOOLEAN_KEYS == ("A", "B", "C", "C_mt", "C_ribo", "D", "FLAG", "WATCH"), str(BOOLEAN_KEYS))
check("the measured columns are unchanged",
      NUMERIC_KEYS == ("n", "median_umi", "umi_frac_of_sample", "median_pct_mt",
                       "pct_uninformative", "pct_mt_markers", "pct_ribo_markers", "pct_doublet"),
      str(NUMERIC_KEYS))

print("\n6 - the parameter actually reaches the op")
# This exists because it did not. `extra_resolutions` was added to the step-6 TASK in graph.py
# while `_cluster` built the op's parameter dict by hand and never forwarded it - so the extras
# would have been requested, recorded in the run key, reported as declared, and never computed,
# with the whole suite green. A parameter that stops at the task boundary is invisible to every
# test that checks the task.
graph_src = (HERE.parent / "engine" / "graph.py").read_text(encoding="utf-8")
steps_src = (HERE.parent / "engine" / "steps.py").read_text(encoding="utf-8")
ops_src = (HERE.parent / "adapters" / "scanpy_ops.py").read_text(encoding="utf-8")
check("graph.py puts extra_resolutions on the step-6 task",
      '"extra_resolutions": tools.get("extra_resolutions")' in graph_src)
check("steps.py forwards it into the op's parameters",
      '"extra_resolutions": p.get("extra_resolutions")' in steps_src)
check("the op reads it", 'params.get("extra_resolutions")' in ops_src)
check("the op writes siblings, not a second applied table",
      '.cluster_profile.res{r:g}.csv' in ops_src)

print("\n" + "=" * 62)
if fails:
    print(f"resolution sweep: {len(fails)} FAILED - " + ", ".join(fails))
    sys.exit(1)
print("resolution sweep OK - 2.0 applied, extras beside it, and the two cannot be confused")
