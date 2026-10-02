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
# Single cells are NOT denoised here - the PI's ruling of 2026-10-02, and DENOISE in the module.
# Case 2 read "scrna, default runs" until then, and case 5 refused a single-cell skip given no
# reason; neither branch exists now, because the skip is the pipeline's rule rather than a request.
p2 = expect_ok("2 scrna, default: NOT run", sample="pbmc_donor1", assay="scrna",
               intronic_fraction=0.18)
if p2 is not None and (p2.run or p2.state != "SKIP" or "NOT RUN" not in p2.reason):
    fails.append(f"2: a single-cell library was planned {p2.state} (run={p2.run}); the policy "
                 f"is that this pipeline never denoises one")
print()
p = expect_ok("3 scrna, a reason given anyway is kept", sample="pbmc_donor2", assay="scrna",
              skip=True, skip_reason="ambient fraction measured at 1.2% in a pilot",
              intronic_fraction=0.16)
if p is not None and (p.run or "1.2%" not in p.reason):
    fails.append("3: the reason given beside the policy was dropped")
print()
p3 = expect_ok("3b scrna, SUPPLIED is still accepted", sample="pbmc_donor4", assay="scrna",
               supplied={"tool": "CellBender", "version": "0.3.2", "params": "fpr=0.0",
                         "produced_by": "a collaborator"})
if p3 is not None and p3.state != "SUPPLIED":
    fails.append(f"3b: a supplied single-cell correction was planned {p3.state}")
print()
print("-" * 74)
expect_refusal("4 snrna, skip requested", sample="cohort_ctrl_03", assay="snrna",
               skip=True, skip_reason="in a hurry", intronic_fraction=0.55)
print()
expect_ok("5 scrna, skip with no reason: the policy is the reason", sample="pbmc_donor3",
          assay="scrna", skip=True, intronic_fraction=0.20)
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
print("proved: every branch behaves as specified, and case 6 is the one that earns the check -")
print("a mis-declared assay would otherwise SKIP a correction that removed 15.5-23.7% of all")
print("counts in the calibration cohort, and say nothing about having skipped it")
