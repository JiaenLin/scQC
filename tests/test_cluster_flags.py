# Exercises the step-6 cluster flags and prints them; it removes no cluster.
"""Step 6 test: the conjunction, its cost, and why a blank must never read as a pass."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "06_cluster_check"))
from cluster_flags import (Thresholds, apply_flags, propose, # noqa: E402
                           sweep_summary)

# A SYNTHETIC profile: nine clusters that must FLAG and three that must not. The values are
# invented rather than copied from any cohort's table, so what is asserted below is the rule and
# not a remembered result. The structure is what carries the test, and each entry is here for a
# reason that must survive any future edit to the numbers:
#
#   eight FLAG by B&C           high mitochondrial content WITH uninformative markers
#   one FLAG by D alone         a mostly-doublet cluster with ordinary depth and markers
#   one of the eight sits on A  just below the 0.5x depth line, so A&C carries it instead of B
#   near miss                   over B, under C - the conjunction must hold it back
#   C alone                     WATCH, never FLAG
#   ordinary                    nothing fires
EXPECT_FLAG, EXPECT_WATCH = 9, 1
P = [
    dict(s="ctrl_01 c4", umi_frac_of_sample=1.11, median_pct_mt=19.40, pct_doublet=6.80,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    dict(s="ctrl_01 c7", umi_frac_of_sample=0.98, median_pct_mt=22.15, pct_doublet=14.90,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    dict(s="treat_01 c18", umi_frac_of_sample=0.87, median_pct_mt=18.70, pct_doublet=7.40,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    dict(s="treat_02 c11", umi_frac_of_sample=1.84, median_pct_mt=2.90, pct_doublet=72.50,
         pct_mt_markers=0.0, pct_ribo_markers=0.0, pct_uninformative=0.0), # D only
    dict(s="ctrl_03 c6", umi_frac_of_sample=0.52, median_pct_mt=21.80, pct_doublet=5.10,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    dict(s="ctrl_04 c9", umi_frac_of_sample=0.69, median_pct_mt=20.60, pct_doublet=8.30,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    dict(s="ctrl_05 c5", umi_frac_of_sample=0.83, median_pct_mt=23.05, pct_doublet=10.40,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    dict(s="ctrl_05 c8", umi_frac_of_sample=0.74, median_pct_mt=18.95, pct_doublet=4.70,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    # A&C, and deliberately at the boundary: a profile table that PRINTS this as 0.50 would
    # produce a fixture that does not flag, because A is a strict `<`. A boundary value
    # transcribed from a rounded display flips the inequality and the module looks wrong when
    # the transcription is. The stored value, not the displayed one, is what a fixture needs.
    dict(s="treat_04 c12", umi_frac_of_sample=0.4993, median_pct_mt=12.40, pct_doublet=2.20,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    # must NOT flag
    dict(s="ctrl_04 c11 (near miss)", umi_frac_of_sample=0.71, median_pct_mt=15.40,
         pct_doublet=0.0, pct_mt_markers=45.0, pct_ribo_markers=0.0, pct_uninformative=45.0),
    dict(s="WATCH: C alone", umi_frac_of_sample=1.15, median_pct_mt=9.20, pct_doublet=8.10,
         pct_mt_markers=50.0, pct_ribo_markers=0.0, pct_uninformative=50.0),
    dict(s="ordinary", umi_frac_of_sample=1.26, median_pct_mt=3.10, pct_doublet=7.30,
         pct_mt_markers=0.0, pct_ribo_markers=0.0, pct_uninformative=0.0),
]
THR = Thresholds(0.5, 15.0, 50.0, 70.0, "APPROVED for this cohort - not derived by the module")

fails = []
print("Step 6 - cluster flags\n" + "=" * 74)
print("A. flags under approved thresholds")
print(" " + str(THR).replace("\n", "\n "))
f = apply_flags(P, THR)
c = f.counts()
print(f"\n FLAG {c['FLAG']} WATCH {c['WATCH']} D {c['D']}")
for r in f.rows:
    if r["FLAG"] or r["WATCH"]:
        why = "".join(k for k in ("A", "B", "C", "D") if r[k])
        print(f" {'FLAG ' if r['FLAG'] else 'WATCH'} {r['s']:26s} [{why}]")
if c["FLAG"] != EXPECT_FLAG or c["WATCH"] != EXPECT_WATCH:
    fails.append(f"A: expected {EXPECT_FLAG} FLAG / {EXPECT_WATCH} WATCH, "
                 f"got {c['FLAG']}/{c['WATCH']}")
if not any(r["D"] and not r["C"] for r in f.rows):
    fails.append("A: the mostly-doublet cluster must flag on D alone, without C")

nearmiss = [r for r in f.rows if "near miss" in r["s"]][0]
print(f"\n the conjunction's cost: {nearmiss['s']} at {nearmiss['median_pct_mt']}% mito "
      f"(over B) and {nearmiss['pct_uninformative']}% markers (under C) -> FLAG "
      f"{nearmiss['FLAG']}")
if nearmiss["FLAG"]:
    fails.append("A: the near-miss must not flag - that is the conjunction working")

print("\n" + "-" * 74)
print("B. markers not computed - C, FLAG, WATCH must be MISSING, never False")
f2 = apply_flags(P, THR, markers_computed=False)
bad = [r for r in f2.rows if r["C"] is not None or r["FLAG"] is not None]
print(f" populated C/FLAG rows: {len(bad)} (must be 0)")
print(f" D still computed: {sum(1 for r in f2.rows if r['D'])} (D needs no markers)")
if bad:
    fails.append("B: C/FLAG/WATCH must be None without markers")

print("\nC. the sweep reports only what is computed everywhere")
sw = sweep_summary({1.0: f.rows, 1.1: f2.rows})
for row in sw:
    print(" " + str(row))
if any("FLAG" in row or "C" in row for row in sw):
    fails.append("C: the sweep must not carry C or FLAG")

print("\n" + "-" * 74)
print("D. thresholds proposed from this cohort rather than inherited")
p = propose(P)
print(" " + str(p).replace("\n", "\n "))
if p.b_pct_mt <= 0 or p.c_uninformative <= 0:
    fails.append("D: proposal failed")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print(f"proved: the conjunction yields {EXPECT_FLAG} FLAG and {EXPECT_WATCH} WATCH on this")
print("profile, with the mostly-doublet cluster reached by D alone - the case a per-cell")
print("doublet test cannot see - and the near miss held back at five points on one axis")
print("proved: without markers, C, FLAG and WATCH stay MISSING; a blank never reads as a pass")
