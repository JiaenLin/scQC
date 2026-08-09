# Assesses doublet-scoring coverage against the light floor and prints it; it filters nothing.
"""Step 3 test: coverage bookkeeping, the rates that need their denominator, and the refusals."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "03_light_floor"))
from light_floor import FloorRefusal, assess # noqa: E402

fails = []
print("Step 3 - light floor and coverage bookkeeping\n" + "=" * 74)

print("A. the calibration cohort's coverage, as run")
c = assess(n_total=244968, n_scored=185646, floor=200, quality_floor=350,
           max_unscored_umi=199, min_scored_umi=200, n_kept_unscored=0)
print(c)
if c.n_unscored != 59322 or round(c.pct_scored, 2) != 75.78:
    fails.append("A: the calibration cohort's coverage does not reproduce")

print("\nB. the denominator must travel with the rate (22,656 called)")
print(" " + c.rate(22656, "scored"))
print(" " + c.rate(22656, "all"))
try:
    c.rate(22656, denominator="")
    fails.append("B: a rate with no denominator must be refused")
except ValueError as e:
    print(f" [REFUSED] bare rate: {e}")

print("\n" + "-" * 74)
print("C. the light floor creeping up to the quality floor")
try:
    assess(n_total=244968, n_scored=150000, floor=350, quality_floor=350)
    fails.append("C: floor == quality floor must refuse")
except FloorRefusal as e:
    print(f"[REFUSED] {str(e)[:190]}...")

print("\nD. a nucleus scored below the floor (the scored set is not what the floor defines)")
try:
    assess(n_total=1000, n_scored=900, floor=200, quality_floor=350, min_scored_umi=87)
    fails.append("D: a sub-floor scored nucleus must refuse")
except FloorRefusal as e:
    print(f"[REFUSED] {str(e)[:170]}...")

print("\nE. retained nuclei that were never examined")
c2 = assess(n_total=244968, n_scored=185646, floor=200, quality_floor=250,
            n_kept_unscored=4127)
print(c2)
if "UNKNOWN, not negative" not in str(c2):
    fails.append("E: retained-unscored nuclei must be marked UNKNOWN, not negative")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: coverage reproduces at 185,646 of 244,968 scored, the floor boundary is clean")
print("and no retained nucleus went unscored - and a rate refuses to be quoted without its")
print("denominator, because 5.70% and 7.52% are BOTH true of the same count and comparing")
print("either against a published band means nothing until you know which one it is")
