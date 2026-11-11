# Gates recorded aligner and denoiser cell counts and prints the findings; it calls no cells.
"""Step 2 test: a synthetic pair of cell calls, plus the two faults the gate exists to stop."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "02_cells"))
from cellcall_gate import gate, verdict # noqa: E402

# Synthetic calls in the shape of an accepted nuclei cohort's: CellBender calls more cells than the
# aligner everywhere, one library loses exactly 5% of the aligner's cells (ON the REVIEW line), and
# every loss falls on one arm. The losses are invented, not any cohort's. `lost` is aligner cells NOT called by CellBender.
COHORT = {
    "ctrl_03": dict(aligner=20000, cellbender=36000, lost=0),
    "ctrl_04": dict(aligner=16000, cellbender=29000, lost=0),
    "ctrl_05": dict(aligner=13000, cellbender=22000, lost=0),
    "treat_03": dict(aligner=10000, cellbender=14000, lost=500),
    "treat_04": dict(aligner=14000, cellbender=21000, lost=12),
    "treat_05": dict(aligner=10500, cellbender=14000, lost=3),
    "ctrl_01": dict(aligner=17500, cellbender=32000, lost=0),
    "ctrl_02": dict(aligner=13000, cellbender=21000, lost=0),
    "treat_01": dict(aligner=17500, cellbender=29500, lost=0),
    "treat_02": dict(aligner=17500, cellbender=21500, lost=60),
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
print("A. a nuclei cohort's calls (synthetic), as accepted")
f = gate(COHORT, DESIGN, assay="snrna")
for x in f:
    print(x)
v = verdict(f)
print(f"\n verdict: {v}")
if v == "REFUSE":
    fails.append("A: the accepted cohort must not be refused")
if not any(x.severity == "REVIEW" and "condition" in x.check for x in f):
    fails.append("A: the one-sided condition loss must be surfaced")
if not any(x.severity == "REVIEW" and x.check.startswith("aligner cells lost")
           and any(d.startswith("treat_03:") for d in x.detail) for x in f):
    fails.append("A: a loss exactly on the 5% line must be REVIEWed (the line is inclusive)")

print("\n" + "-" * 74)
print("B. injected: CellBender stricter than the aligner")
# STRICTER IS A FILTER WHEN IT COSTS CELLS (harness autodev W020, B1-103): a ratio below 1.0 with
# under REVIEW_LOST of the aligner's cells lost is REVIEW; with REVIEW_LOST or more it is REFUSED.
b = {k: dict(v) for k, v in COHORT.items()}
b["ctrl_05"]["cellbender"] = 11000 # below the aligner's 13,000
f2 = gate(b, DESIGN, assay="snrna")
hit = [x for x in f2 if "stricter" in x.check]
print(hit[0])
if hit[0].severity != "REVIEW":
    fails.append("B: a ratio below 1.0 that loses under the REVIEW line must be REVIEW")
b["ctrl_05"]["lost"] = int(0.06 * b["ctrl_05"]["aligner"])
hit = [x for x in gate(b, DESIGN, assay="snrna") if "stricter" in x.check]
print(hit[0])
if hit[0].severity != "REFUSE":
    fails.append("B: a ratio below 1.0 that loses 6% of the aligner's cells must REFUSE")

print("\n" + "-" * 74)
print("C. injected: 14% of aligner cells lost in one sample")
c = {k: dict(v) for k, v in COHORT.items()}
c["ctrl_02"]["lost"] = int(0.14 * c["ctrl_02"]["aligner"])
f3 = gate(c, DESIGN, assay="snrna")
hit = [x for x in f3 if x.check.startswith("aligner cells lost")]
print(hit[0])
if hit[0].severity != "REFUSE":
    fails.append("C: >10% lost must REFUSE")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: the accepted cohort clears the REFUSE lines and is REVIEWed where it should be -")
print(" treat_03 sits exactly on the 5% review line, so the boundary is inclusive")
print(" every lost cell falls on treat and none on ctrl, which no per-sample % can see")
