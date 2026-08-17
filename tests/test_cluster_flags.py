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
#   one mostly-doublet cluster  ordinary depth and markers: it must NOT flag, because
#                               criterion D was removed on 2026-08-15 (see the end of
#                               this file). Kept in the fixture precisely so that a
#                               reinstated D would light it up and fail the run.
#   one of the eight sits on A  just below the 0.5x depth line, so A&C carries it instead of B
#   near miss                   over B, under C - the conjunction must hold it back
#   C alone                     WATCH, never FLAG
#   ordinary                    nothing fires
# 8, not 9: the fixture's mostly-doublet cluster no longer flags. Criterion D was removed on
# 2026-08-15 and that cluster fires nothing else - which is exactly why D went.
EXPECT_FLAG, EXPECT_WATCH = 8, 1
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
THR = Thresholds(0.5, 15.0, 50.0, "APPROVED for this cohort - not derived by the module")

fails = []
print("Step 6 - cluster flags\n" + "=" * 74)
print("A. flags under approved thresholds")
print(" " + str(THR).replace("\n", "\n "))
f = apply_flags(P, THR)
c = f.counts()
print(f"\n FLAG {c['FLAG']} WATCH {c['WATCH']}")
for r in f.rows:
    if r["FLAG"] or r["WATCH"]:
        why = "".join(k for k in ("A", "B", "C") if r[k])
        print(f" {'FLAG ' if r['FLAG'] else 'WATCH'} {r['s']:26s} [{why}]")
if c["FLAG"] != EXPECT_FLAG or c["WATCH"] != EXPECT_WATCH:
    fails.append(f"A: expected {EXPECT_FLAG} FLAG / {EXPECT_WATCH} WATCH, "
                 f"got {c['FLAG']}/{c['WATCH']}")

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

print("\n" + "-" * 74)
print("E. the profile survives being written and read back")
#
# The flags are computed in the ORCHESTRATOR, from a CSV the analysis environment wrote, because
# the two do not share an interpreter. A CSV carries no types: every cell comes back a string,
# `"0.71" < 0.5` raises TypeError, and where a comparison happens to be string-to-string it
# succeeds and sorts lexicographically. The blank cell the writer uses for "not computed" comes
# back as "", which `_unknown` does not recognise. So the round trip is part of the contract and
# is tested as one: the same profile, through a file, must produce the same verdicts.
import csv as _csv                                                        # noqa: E402
import tempfile                                                           # noqa: E402

sys.path.insert(0, str(HERE.parent))
from adapters.scanpy_ops import write_profile_csv                          # noqa: E402
from cluster_flags import ClusterRefusal, read_profile_csv                 # noqa: E402

with tempfile.TemporaryDirectory() as _d:
    _p = write_profile_csv(P, Path(_d) / "profile.csv")
    back = read_profile_csv(_p)
    strings = sorted({k for r in back for k, v in r.items()
                      if k in ("umi_frac_of_sample", "median_pct_mt", "pct_doublet",
                               "pct_uninformative", "n", "median_umi")
                      and isinstance(v, str)})
    print(f" numeric columns still carrying strings: {strings or 'none'}")
    if strings:
        fails.append(f"E: {strings} came back as text; every threshold on them would raise or "
                     f"compare lexicographically")
    blanks = sum(1 for r in back if r["median_umi"] is None)
    print(f" blank cells read as unknown (not 0.0, not ''): {blanks} of {len(back)}")
    if blanks != len(back):
        fails.append("E: a blank cell must read back as None")

    c_back = apply_flags(back, THR).counts()
    print(f" verdicts through the file: FLAG {c_back['FLAG']} WATCH {c_back['WATCH']}"
          f"   (in memory: {c['FLAG']}/{c['WATCH']})")
    if (c_back["FLAG"], c_back["WATCH"]) != (c["FLAG"], c["WATCH"]):
        fails.append("E: the round trip changed the verdicts")

    # A FLAGGED profile, written and read back, must still answer `is True`. Step 7 selects with
    # `[c for c in profile if c.get("FLAG") is True]`, and `"True" is True` is False - so a
    # verdict left as text does not raise, it reports a table of 126 flagged clusters as having
    # none, and the deliverable is described as sitting in no flagged cluster at all.
    _fp = write_profile_csv(f.rows, Path(_d) / "flagged.csv")
    _fb = read_profile_csv(_fp)
    kinds = sorted({type(r[k]).__name__ for r in _fb for k in ("A", "B", "C", "FLAG")})
    n_is_true = sum(1 for r in _fb if r.get("FLAG") is True)
    print(f" verdict types after the round trip: {kinds}")
    print(f" clusters answering `FLAG is True`: {n_is_true} (in memory: {c['FLAG']})")
    if n_is_true != c["FLAG"]:
        fails.append(f"E: {n_is_true} of {c['FLAG']} flagged clusters survive `is True` through "
                     f"the file; step 7 selects on exactly that")

    # A column that did not survive is a defect, and must not be absorbed as "unknown".
    _bad = Path(_d) / "corrupt.csv"
    with open(_p, encoding="utf-8", newline="") as _fh:
        _rows = list(_csv.DictReader(_fh))
    _rows[0]["median_pct_mt"] = "19.40%"
    with open(_bad, "w", encoding="utf-8", newline="") as _fh:
        _w = _csv.DictWriter(_fh, fieldnames=list(_rows[0]))
        _w.writeheader()
        _w.writerows(_rows)
    try:
        read_profile_csv(_bad)
        fails.append("E: a non-numeric, non-blank cell must be refused, not read as unknown")
        print(" unparseable cell: ACCEPTED - wrong")
    except ClusterRefusal as _e:
        print(f" unparseable cell refused: {str(_e).splitlines()[0][:96]}")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print(f"proved: the conjunction yields {EXPECT_FLAG} FLAG and {EXPECT_WATCH} WATCH on this")
print("profile, with the mostly-doublet cluster reached by D alone - the case a per-cell")
print("doublet test cannot see - and the near miss held back at five points on one axis")
print("proved: without markers, C, FLAG and WATCH stay MISSING; a blank never reads as a pass")


# ---------------------------------------------------------------- D is gone, and stays gone
#
# Criterion D flagged a whole CLUSTER when more than 70% of its cells were called doublets.
# Retired 2026-08-15 on PI instruction: *"just let the scDblfinder alone handle the doublet,
# remove D is good as no evidence prove they are also doublet"*.
#
# The measurement behind it: on the calibration cohort D fired on exactly one cluster in nearly two
# hundred, and that cluster fired nothing else - not low UMI, not high mitochondrial, not marker
# share. Its entire marginal contribution was a few dozen nuclei that the per-cell caller had
# DECLINED to call doublets. So
# the criterion's only effect was to remove cells by association with their neighbours, against
# per-cell evidence that they were not doublets.
#
# It is REMOVED, not defaulted off, for the reason scAnno's FLAG_SHARE was removed: a disabled
# threshold is one edit away from returning, and the next person to see the parameter will assume
# it was meant to be used. Per-cell `doublet_class` and `doublet_score` survive untouched, so a
# consumer that wants to act on doublets still can - on evidence about the cell in front of it.
def _chk(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not cond: fails.append(name)

print("\nF. criterion D is removed, and cannot come back")
_root = HERE.parent
_src = (_root / "modules/06_cluster_check/cluster_flags.py").read_text(encoding="utf-8")
_eng = (_root / "engine/steps.py").read_text(encoding="utf-8")
import re as _re
# WORD BOUNDARIES. A substring test reports `n_called_doublet` as a surviving `d_doublet`, which
# is a false alarm about a real rule - and a gate that cries wolf is a gate that gets deleted.
for _name, _text in (("cluster_flags.py", _src), ("steps.py", _eng)):
    for _tok in ("d_doublet", "D_DEFAULT"):
        _chk(f"{_name} has no {_tok}", not _re.search(rf"\b{_tok}\b", _text))
import cluster_flags as _cf
_chk("BOOLEAN_KEYS has no D", "D" not in _cf.BOOLEAN_KEYS, str(_cf.BOOLEAN_KEYS))
_chk("Thresholds has no d_doublet field",
      not hasattr(Thresholds(0.5, 1.0, 1.0), "d_doublet"))
_chk("apply_flags writes no D column",
      all("D" not in r for r in apply_flags(P, Thresholds(0.5, 1.0, 1.0)).rows))
_chk("the threshold description does not mention D",
      "D >" not in str(Thresholds(0.5, 1.0, 1.0)), str(Thresholds(0.5, 1.0, 1.0)))
