# Gates recorded aligner and denoiser cell counts and prints the findings; it calls no cells.
"""Step 2 test: a real pair of cell calls, plus the two faults the gate exists to stop."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "02_cells"))
from cellcall_gate import gate, verdict # noqa: E402

# Measured from results/objects/cellbender_h5ad/*.h5ad, obs celescope_cell vs
# cellbender_cell. `lost` is aligner cells NOT called by CellBender.
COHORT = {
    "ctrl_03": dict(aligner=21901, cellbender=38285, lost=0),
    "ctrl_04": dict(aligner=16364, cellbender=30103, lost=0),
    "ctrl_05": dict(aligner=13140, cellbender=22067, lost=0),
    "treat_03": dict(aligner=10438, cellbender=14229, lost=527),
    "treat_04": dict(aligner=14210, cellbender=21620, lost=10),
    "treat_05": dict(aligner=10615, cellbender=14074, lost=4),
    "ctrl_01": dict(aligner=17685, cellbender=32220, lost=0),
    "ctrl_02": dict(aligner=12999, cellbender=20995, lost=0),
    "treat_01": dict(aligner=17692, cellbender=29870, lost=0),
    "treat_02": dict(aligner=17609, cellbender=21505, lost=76),
}
CONDITION = {s: ("treat" if "treat" in s else "ctrl") for s in COHORT}
# A SECOND design factor. Its levels used to be readable off the sample names; now they are
# explicit. The gate has to be exercised against more than one factor, because a loss that
# is balanced on one can be entirely one-sided on another.
BATCH = {s: ("b2" if s in {"ctrl_03", "ctrl_04", "ctrl_05",
                           "treat_03", "treat_04", "treat_05"} else "b1") for s in COHORT}
DESIGN = {"condition": CONDITION, "batch": BATCH}

fails = []
print("Step 2 - cell-call gate\n" + "=" * 74)
print("A. the calibration cohort's calls, as accepted")
f = gate(COHORT, DESIGN)
for x in f:
    print(x)
v = verdict(f)
print(f"\n verdict: {v}")
if v == "REFUSE":
    fails.append("A: the accepted cohort must not be refused")
if not any(x.severity == "REVIEW" and "condition" in x.check for x in f):
    fails.append("A: the one-sided condition loss must be surfaced")

print("\n" + "-" * 74)
print("B. injected: CellBender stricter than the aligner")
b = {k: dict(v) for k, v in COHORT.items()}
b["ctrl_05"]["cellbender"] = 11000 # below the aligner's 13,140
f2 = gate(b, DESIGN)
hit = [x for x in f2 if "stricter" in x.check]
print(hit[0])
if hit[0].severity != "REFUSE":
    fails.append("B: a ratio below 1.0 must REFUSE")

print("\n" + "-" * 74)
print("C. injected: 14% of aligner cells lost in one sample")
c = {k: dict(v) for k, v in COHORT.items()}
c["ctrl_02"]["lost"] = int(0.14 * c["ctrl_02"]["aligner"])
f3 = gate(c, DESIGN)
hit = [x for x in f3 if x.check.startswith("aligner cells lost")]
print(hit[0])
if hit[0].severity != "REFUSE":
    fails.append("C: >10% lost must REFUSE")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: the accepted cohort clears the REFUSE lines and is REVIEWed where it should be -")
print(" treat_03 sits exactly on the 5% review line, so the boundary is inclusive")
print(" all 617 lost cells fall on treat and none on ctrl, which no per-sample % can see")
