# Assesses ambient-correction diagnostics against the learning-rate policy and prints.
"""Step 1b test: two learning-rate experiments, and a policy that has to get both right.

The pair is the point. Run on the same tool and the same data, the two experiments reach
OPPOSITE conclusions, and the only difference between them is how many samples are looked at.
A policy that reproduces one of them and not the other has learned the wrong lesson.

Every value here is CONSTRUCTED - ordinary libraries one step apart around a centre, one
degenerate library far from it - and every expected outlier is COMPUTED below from the rule the
policy documents, so no value stands for any cohort's.
"""
import sys
from pathlib import Path
from statistics import NormalDist, median

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "01_ambient"))
from lr_policy import assess_cohort, compare_after_halving # noqa: E402

# THE RULE THIS TEST HOLDS THE POLICY TO, stated here rather than read from the module: a library
# is unlike its siblings when its robust z - the distance from the median in units of the MAD
# scaled to a normal SD - is above MAD_CUT, and nothing is judged below MIN_LIBRARIES.
MAD_CUT, MIN_LIBRARIES = 3.5, 4
SCALE = 1 / NormalDist().inv_cdf(3 / 4)

LIBS = ["ctrl_01", "ctrl_02", "ctrl_03", "ctrl_04", "treat_01", "treat_02", "treat_03", "treat_04"]
DEGENERATE = "treat_02"
FRACTION, FRACTION_STEP = 1 / 5, 1 / 80          # share of counts removed: centre and step
CONVERGENCE, CONVERGENCE_STEP = 1 / 2, 1 / 16    # convergence indicator: centre and step
OFFSETS = [-3, -2, -1, 0, 1, 2, 3]               # the ordinary libraries, one step apart


def ordinary(fraction_offset, convergence_offset):
    return {"fraction_removed": FRACTION + FRACTION_STEP * fraction_offset,
            "convergence_indicator": CONVERGENCE + CONVERGENCE_STEP * convergence_offset}


# At the package default 1e-4: the degenerate fit removes a quarter of the centre's share and its
# convergence indicator is four times the centre; the others are ordinary, their two offsets
# paired in opposite orders so neither diagnostic mirrors the other.
others = [s for s in LIBS if s != DEGENERATE]
COHORT_DEFAULT = {s: ordinary(a, b) for s, a, b in zip(others, OFFSETS, reversed(OFFSETS))}
COHORT_DEFAULT[DEGENERATE] = {"fraction_removed": FRACTION / 4,
                              "convergence_indicator": CONVERGENCE * 4}
# At half the rate, 5e-5, the degenerate library becomes ordinary.
COHORT_HALF = {k: dict(v) for k, v in COHORT_DEFAULT.items()}
COHORT_HALF[DEGENERATE] = ordinary(0, 0)


def expected(diag):
    """{(library, diagnostic)} the documented rule flags."""
    if len(diag) < MIN_LIBRARIES:
        return set()
    out = set()
    for d in ("fraction_removed", "convergence_indicator"):
        vals = {s: m[d] for s, m in diag.items()}
        med = median(vals.values())
        mad = median(abs(x - med) for x in vals.values())
        out |= {(s, d) for s, x in vals.items() if mad > 0 and abs(x - med) / (SCALE * mad) > MAD_CUT}
    return out


def flagged(verdict):
    """{(library, diagnostic)} the policy reports."""
    return {(s, why.split(" = ")[0]) for s, reasons in verdict.outliers.items() for why in reasons}


fails = []
print("Step 1b - learning-rate policy\n" + "=" * 74)

print("A. the one-sample experiment: one ordinary library and nothing else")
one = {"treat_03": COHORT_DEFAULT["treat_03"]}
v = assess_cohort(one)
print(v)
if v.action != "keep_default" or flagged(v):
    fails.append("one-sample: must not act - a degenerate fit is only visible against siblings")
print(" -> a single sample cannot see the problem: there is nothing to be an outlier against\n")

print(f"B. the cohort experiment: all {len(LIBS)} at the package default")
v = assess_cohort(COHORT_DEFAULT, label="1e-4")
print(v)
want = expected(COHORT_DEFAULT)
if want != {(DEGENERATE, "fraction_removed"), (DEGENERATE, "convergence_indicator")}:
    fails.append(f"cohort: the construction is wrong - the rule flags {sorted(want)}")
if flagged(v) != want or v.action != "rerun_half":
    fails.append(f"cohort: the rule flags {sorted(want)} and a half-rate re-run; the policy flagged "
                 f"{sorted(flagged(v))} and chose {v.action!r}")
print()

print("C. after halving - did it resolve?")
if expected(COHORT_HALF):
    fails.append(f"post-halving: the construction is wrong - the rule flags {sorted(expected(COHORT_HALF))}")
v2 = compare_after_halving(COHORT_DEFAULT, COHORT_HALF, flagged=list(v.outliers))
print(v2)
if v2.action != "adopt_half":
    fails.append("post-halving: the degenerate library resolves, so half must be adopted cohort-wide")
print()

print("D. the case where halving fixes NOTHING - the diagnostic simply moves")
unchanged = {k: dict(x) for k, x in COHORT_DEFAULT.items()}
v3 = compare_after_halving(COHORT_DEFAULT, unchanged, flagged=[DEGENERATE])
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
