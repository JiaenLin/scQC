# Assesses recorded ambient-correction diagnostics against the learning-rate policy and prints.
"""Step 1b test: two learning-rate experiments, and a policy that has to get both right.

The pair is the point. Run on the same tool and the same data, the two experiments reached
OPPOSITE conclusions, and the only difference between them was how many samples were looked at.
A policy that reproduces one of them and not the other has learned the wrong lesson.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "01_ambient"))
from lr_policy import assess_cohort, compare_after_halving # noqa: E402

# The cohort experiment: all ten libraries at the package default 1e-4. treat_02 is the
# degenerate fit - a clearly smaller share of counts removed than elsewhere, a convergence
# indicator several times the others'. The rest are spread across ordinary ranges (synthetic).
COHORT_DEFAULT = {
    "ctrl_01": {"fraction_removed": 0.18, "convergence_indicator": 0.5},
    "ctrl_02": {"fraction_removed": 0.16, "convergence_indicator": 0.45},
    "treat_01": {"fraction_removed": 0.17, "convergence_indicator": 0.6},
    "treat_02": {"fraction_removed": 0.08, "convergence_indicator": 1.8},  # degenerate
    "ctrl_03": {"fraction_removed": 0.165, "convergence_indicator": 0.4},
    "ctrl_04": {"fraction_removed": 0.155, "convergence_indicator": 0.5},
    "ctrl_05": {"fraction_removed": 0.13, "convergence_indicator": 0.55},
    "treat_03": {"fraction_removed": 0.14, "convergence_indicator": 0.7},
    "treat_04": {"fraction_removed": 0.125, "convergence_indicator": 0.65},
    "treat_05": {"fraction_removed": 0.115, "convergence_indicator": 0.55},
}
# At half the rate, 5e-5, the outlier becomes ordinary.
COHORT_HALF = {k: dict(v) for k, v in COHORT_DEFAULT.items()}
COHORT_HALF["treat_02"] = {"fraction_removed": 0.15, "convergence_indicator": 0.60}
COHORT_HALF["treat_03"] = {"fraction_removed": 0.135, "convergence_indicator": 0.60}

fails = []
print("Step 1b - learning-rate policy\n" + "=" * 74)

print("A. the one-sample experiment: treat_03 and nothing else")
one = {"treat_03": COHORT_DEFAULT["treat_03"]}
v = assess_cohort(one)
print(v)
if v.action != "keep_default":
    fails.append("one-sample: must not act - a degenerate fit is only visible against siblings")
print(" -> a single sample cannot see the problem: there is nothing to be an outlier against\n")

print("B. the cohort experiment: all ten at the package default")
v = assess_cohort(COHORT_DEFAULT, label="1e-4")
print(v)
if v.action != "rerun_half" or "treat_02" not in v.outliers:
    fails.append("cohort: must flag treat_02 and trigger a half-rate re-run")
print()

print("C. after halving - did it resolve?")
v2 = compare_after_halving(COHORT_DEFAULT, COHORT_HALF, flagged=list(v.outliers))
print(v2)
if v2.action != "adopt_half":
    fails.append("post-halving: treat_02 resolves, so half must be adopted cohort-wide")
print()

print("D. the case where halving fixes NOTHING - the diagnostic simply moves")
unchanged = {k: dict(x) for k, x in COHORT_DEFAULT.items()}
v3 = compare_after_halving(COHORT_DEFAULT, unchanged, flagged=["treat_02"])
print(v3)
if v3.action != "keep_default":
    fails.append("unresolved: must retain the default rather than adopt a fix that fixes nothing")

print("\n" + "=" * 74)
if fails:
    print("FAILED:")
    for f in fails:
        print(" -", f)
    raise SystemExit(1)
print("proved: all four branches behave as specified, including the two that reach opposite")
print("conclusions - and A and B differ ONLY in how many samples were examined")
