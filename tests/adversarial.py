# Calls the pipeline's functions with hostile inputs and prints what they do; removes nothing.
"""Adversarial review: specific bug classes, tested rather than read for.

Re-reading code finds little. These are hypotheses about failure MODES, each with a concrete
input that would expose it.

Each probe must be able to FAIL. A fixture that cannot reach the failing comparison reports
"tolerated" whether or not the bug is present, which is worse than not probing at all.
"""
import sys
from pathlib import Path

# The repository root, resolved from this file. A hardcoded install path makes the tests
# runnable only on the machine they were written on.
P = Path(__file__).resolve().parents[1]
for m in ("01_ambient", "02_cells", "04_doublets", "05_quality", "06_cluster_check", "07_apply"):
    sys.path.insert(0, str(P / "modules" / m))

from cluster_flags import propose # noqa: E402
from apply import ApplyRefusal, apply_removal, preflight # noqa: E402
from cluster_flags import sweep_summary # noqa: E402
from lr_policy import assess_cohort # noqa: E402

found = []
print("Adversarial review\n" + "=" * 74)

# H1 - cluster_flags.propose(): the "bimodal midpoint" arithmetic
print("H1 propose(): does the bimodal branch compute a midpoint?")
prof = [dict(umi_frac_of_sample=1.0, median_pct_mt=2.0, pct_uninformative=0.0)] * 8
prof += [dict(umi_frac_of_sample=1.0, median_pct_mt=20.0, pct_uninformative=44.0),
         dict(umi_frac_of_sample=1.0, median_pct_mt=25.0, pct_uninformative=55.0)]
t = propose(prof)
print(f" bulk = 0, tail starts at 44 -> C proposed at {t.c_uninformative}")
if t.c_uninformative == 44.0:
    found.append("H1 CONFIRMED: propose() returns the START of the tail, not a midpoint. "
                 "`(0 + m)/2 + m/2` simplifies to `m`. The comment says midpoint; the code "
                 "returns min(nonzero). A cluster at exactly the tail minimum is flagged.")
    print(" ^ BUG: that is min(nonzero), not a midpoint between 0 and 44 (which is 22)")

# H2 - preflight(): does a MISSING doublet fraction read as zero?
print("\nH2 preflight(): does pct_doublet=None count as 0% doublet?")
rows = [dict(sample="X", cluster=1, n=1000, pct_doublet=None, median_pct_mt=5.0, FLAG=True)]
out = preflight(rows, kept_total=1000)
msg = [x for x in out if "flagged clusters" in x.check][0].message
print(f" {msg[:100]}")
if "1,000 of 1,000" in msg:
    found.append("H2 CONFIRMED: `(c.get('pct_doublet') or 0)` turns an UNKNOWN doublet "
                 "fraction into 0, so every nucleus is counted as surviving. This is the "
                 "same class as the C/FLAG bug - a blank reading as a value.")
    print(" ^ BUG: unknown treated as 0% doublet")

# H3 - apply_removal(): does whitespace in a legitimate approval cause a false refusal?
print("\nH3 apply_removal(): trailing whitespace in a recorded approval")
ACT = "some action"
try:
    apply_removal(100, 10, ACT, "CONFIRM", {ACT: "CONFIRM\n"})
    print(" tolerated")
except ApplyRefusal:
    found.append("H3 CONFIRMED: an approval stored with a trailing newline refuses a correct "
                 "CONFIRM. A gate that fires on correct behaviour gets switched off, which "
                 "costs more than the whitespace it caught.")
    print(" ^ BUG: refused a legitimate approval over a newline")

# H4 - sweep_summary(): does a None profile value crash the max() across clusters?
# TWO clusters, deliberately. max() over a single element returns that element without ever
# comparing it to anything, so a one-row fixture never reaches the comparison that fails and
# reports "tolerated" against buggy and fixed code alike. The failure needs a None AND a number.
print("\nH4 sweep_summary(): a cluster with no mitochondrial value beside one that has it")
try:
    sweep_summary({1.0: [
        dict(A=False, B=False, D=False, median_pct_mt=None, pct_doublet=1.0),
        dict(A=False, B=False, D=False, median_pct_mt=4.5, pct_doublet=2.0),
    ]})
    print(" tolerated")
except TypeError as e:
    found.append(f"H4 CONFIRMED: sweep_summary() raises on a None value - {e}")
    print(f" ^ BUG: {type(e).__name__}: {e}")

# H5 - _mad_outliers(): does a zero MAD manufacture outliers?
print("\nH5 assess_cohort(): >50% identical values give MAD = 0")
d = {f"s{i}": {"fraction_removed": 0.15, "convergence_indicator": 0.5} for i in range(9)}
d["s9"] = {"fraction_removed": 0.1501, "convergence_indicator": 0.5}
v = assess_cohort(d)
print(f" 9 identical + 1 differing by 0.0001 -> action {v.action}, "
      f"{len(v.outliers)} outlier(s)")
if v.outliers:
    found.append("H5 CONFIRMED: with MAD = 0 the guard `or 1e-12` makes any deviation "
                 "~1e12 robust SD, so a 0.07% difference is flagged as a degenerate fit.")
    print(" ^ BUG: a trivial difference becomes an extreme outlier")

print("\n" + "=" * 74)
if not found:
    print("no hypothesis confirmed")
else:
    print(f"{len(found)} CONFIRMED:")
    for f in found:
        print(" * " + f)
