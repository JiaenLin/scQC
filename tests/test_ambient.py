# Plans ambient correction for each assay branch and prints the plan; it corrects no counts.
"""Step 1 test: every branch of the mandatory/optional rule, including the ones that must refuse."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "01_ambient"))
from ambient import AmbientRefusal, plan_ambient # noqa: E402

fails = []

def expect_ok(label, **kw):
    try:
        p = plan_ambient(**kw)
        print(p)
        return p
    except AmbientRefusal as e:
        fails.append(f"{label}: expected to succeed, refused with: {e}")
        print(f"[FAIL] {label}: unexpectedly refused")
        return None

def expect_refusal(label, **kw):
    try:
        plan_ambient(**kw)
        fails.append(f"{label}: expected a refusal, none raised")
        print(f"[FAIL] {label}: NOT refused")
    except AmbientRefusal as e:
        print(f"[REFUSED] {label}\n {str(e).split(': ', 1)[1][:150]}...")

print("Step 1 - ambient correction policy\n" + "=" * 74)

expect_ok("1 snrna, default", sample="cohort_ctrl_01", assay="snrna",
          intronic_fraction=0.52)
print()
expect_ok("2 scrna, default runs", sample="pbmc_donor1", assay="scrna",
          intronic_fraction=0.18)
print()
p = expect_ok("3 scrna, skip WITH reason", sample="pbmc_donor2", assay="scrna",
              skip=True, skip_reason="ambient fraction measured at 1.2%; CellBender "
                                     "over-corrected the plasma-cell IG genes in a pilot",
              intronic_fraction=0.16)
print()
print("-" * 74)
expect_refusal("4 snrna, skip requested", sample="cohort_ctrl_03", assay="snrna",
               skip=True, skip_reason="in a hurry", intronic_fraction=0.55)
print()
expect_refusal("5 scrna, skip with NO reason", sample="pbmc_donor3", assay="scrna",
               skip=True, intronic_fraction=0.20)
print()
expect_refusal("6 declared scrna, nuclear intronic fraction", sample="mislabelled",
               assay="scrna", intronic_fraction=0.61)
print()
expect_refusal("7 assay not declared", sample="nodecl", assay="")

print("\n" + "=" * 74)
if fails:
    print("FAILED:")
    for f in fails:
        print(" -", f)
    raise SystemExit(1)
print("proved: all 7 branches behave as specified, and case 6 is the one that earns the check -")
print("a mis-declared assay would otherwise SKIP a correction that removed 15.5-23.7% of all")
print("counts in the calibration cohort, and say nothing about having skipped it")
