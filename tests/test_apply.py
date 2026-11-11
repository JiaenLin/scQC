# Exercises the step-7 pre-flight and approval gate and prints the result; removes nothing.
"""Step 7 test: the step-6 contradiction surfaces, and the approval gate cannot be talked past."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "07_apply"))
from apply import (ApplyRefusal, annotate_kept, apply_removal, # noqa: E402
                   preflight, propose_cluster_removal)

# A SYNTHETIC profile of nine flagged clusters at one clustering resolution. The values are
# invented, not measured: what matters to this suite is the SHAPE of the profile, and inventing
# it keeps the fixture from standing in for a real cohort's per-cluster table. The shape that is
# load-bearing, and must survive any future edit to these numbers:
#
#   - exactly one cluster above the 70% cluster-doublet line, so the residual has something to
#     find and propose_cluster_removal() has exactly one victim;
#   - that cluster on ONE design level only, so the one-sided-removal warning fires;
#   - the rest below the line, so they are counted as flagged-but-retained rather than removed.
PROF = [
    dict(sample="ctrl_01", cluster=4, n=1800, median_umi=1200, median_pct_mt=19.0,
         pct_doublet=7.0, FLAG=True, WATCH=False),
    dict(sample="ctrl_01", cluster=7, n=400, median_umi=1100, median_pct_mt=22.0,
         pct_doublet=15.0, FLAG=True, WATCH=False),
    dict(sample="treat_01", cluster=18, n=190, median_umi=1000, median_pct_mt=18.5,
         pct_doublet=7.5, FLAG=True, WATCH=False),
    # the one above the line, and the only one
    dict(sample="treat_02", cluster=11, n=200, median_umi=2980, median_pct_mt=2.90,
         pct_doublet=72.50, FLAG=True, WATCH=False),
    dict(sample="ctrl_03", cluster=6, n=2000, median_umi=950, median_pct_mt=21.5,
         pct_doublet=5.0, FLAG=True, WATCH=False),
    dict(sample="ctrl_04", cluster=9, n=600, median_umi=850, median_pct_mt=20.5,
         pct_doublet=8.0, FLAG=True, WATCH=False),
    dict(sample="ctrl_05", cluster=5, n=250, median_umi=900, median_pct_mt=23.0,
         pct_doublet=10.0, FLAG=True, WATCH=False),
    dict(sample="ctrl_05", cluster=8, n=1000, median_umi=800, median_pct_mt=19.0,
         pct_doublet=4.5, FLAG=True, WATCH=False),
    dict(sample="treat_04", cluster=12, n=900, median_umi=750, median_pct_mt=12.5,
         pct_doublet=2.0, FLAG=True, WATCH=False),
]
N_ABOVE_LINE = sum(c["n"] for c in PROF if c["pct_doublet"] > 70.0)
CONDITION = {"ctrl_01": "ctrl", "ctrl_02": "ctrl", "treat_01": "treat", "treat_02": "treat",
        "ctrl_03": "ctrl", "ctrl_04": "ctrl", "ctrl_05": "ctrl", "treat_03": "treat",
        "treat_04": "treat", "treat_05": "treat"}
ACTION = "a run's filter: UMI>=400, genes>=250, mito<=20%, scDblFinder dbr.sd=0.04"  # synthetic

fails = []
print("Step 7 - apply\n" + "=" * 74)
print("A. pre-flight against a deliverable of 110,000 kept nuclei")
pf = preflight(PROF, kept_total=110000)
for x in pf:
    print(x)
resid = [x for x in pf if "residual" in x.check]
if not resid or resid[0].severity != "REVIEW":
    fails.append("A: survivors inside a >70% doublet cluster must surface as REVIEW")

print("\n" + "-" * 74)
print("B. the cluster removal the cited method would make - prepared, not taken")
p = propose_cluster_removal(PROF, design=CONDITION)
for k in ("rule", "n_removed", "status"):
    print(f" {k}: {p[k]}")
print(f" clusters: {p['clusters']}")
print(f" by design level: {p.get('by_design_level')}")
if "warning" in p:
    print(f" WARNING: {p['warning']}")
if p["n_removed"] != N_ABOVE_LINE:
    fails.append(f"B: expected {N_ABOVE_LINE} nuclei proposed, got {p['n_removed']}")
if "warning" not in p:
    fails.append("B: a removal falling on one design level must carry the one-sided warning")

print("\n" + "-" * 74)
print("C. flags carried into the deliverable")
obs = [dict(sample="treat_02", cluster=11, barcode="bc1"),
       dict(sample="ctrl_01", cluster=99, barcode="bc2")]
annotate_kept(obs, PROF)
for r in obs:
    print(f" {r['barcode']}: FLAG={r['cluster_FLAG']} "
          f"pct_doublet={r['cluster_pct_doublet']}")
if obs[0]["cluster_FLAG"] is not True or obs[1]["cluster_FLAG"] is not None:
    fails.append("C: flags must attach where the cluster is known and be None where not")

print("\n" + "-" * 74)
print("D. the approval gate")
approvals = {ACTION: "CONFIRM"}
n = apply_removal(200000, 90000, ACTION, "CONFIRM", approvals)
print(f" approved: {n:,} kept")
if n != 110000:
    fails.append("D: 200,000 - 90,000 must be 110,000")

for label, kw in (
        ("no verbatim words", dict(action=ACTION, user_verbatim="", approvals=approvals)),
        ("action text changed", dict(action=ACTION + " (revised)", user_verbatim="CONFIRM",
                                     approvals=approvals)),
        ("words from another decision", dict(action=ACTION, user_verbatim="approve this filter.",
                                             approvals=approvals))):
    try:
        apply_removal(200000, 90000, **kw)
        fails.append(f"D: {label} must refuse")
    except ApplyRefusal as e:
        print(f" [REFUSED] {label}: {str(e)[:110]}...")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: nuclei surviving inside a mostly-doublet cluster are surfaced rather than")
print("counted as clean, cluster flags reach the deliverable and stay None where the cluster")
print("is unknown, a removal falling on one design level is warned about because no ratio can")
print("check it, and the approval gate refuses every approval it was not actually given")
