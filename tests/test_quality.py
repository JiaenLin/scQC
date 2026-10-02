# Derives count-floor proposals from recorded valleys and prints them; it applies no threshold.
"""Step 5 test: real density valleys, the bounds around them, and the cases that must refuse."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "05_quality"))
from quality import (Valley, ThresholdRefusal, derive, # noqa: E402
                     mito_ceiling_note, count_bounds, GENE_BOUNDS, UMI_BOUNDS)

S = ["ctrl_01", "ctrl_02", "treat_01", "treat_02", "ctrl_03",
     "ctrl_04", "ctrl_05", "treat_03", "treat_04", "treat_05"]
# the calibration cohort measured: UMI median 348 (274-473), genes median 264 (184-352)
UMI = [473, 300, 361, 348, 274, 336, 352, 289, 383, 344]
GENES = [352, 214, 271, 264, 184, 249, 268, 208, 291, 258]

fails = []
print("Step 5 - quality thresholds\n" + "=" * 74)
print(f"bounds: UMI {UMI_BOUNDS}, genes {GENE_BOUNDS}\n")
# A-H are the calibration cohort's NUCLEI, bounded as nuclei. The bounds are per assay since
# 2026-10-02; section I below is what a single-cell cohort gets, and why it is not this.
SN = "snrna"

print("A. the calibration cohort's real valleys")
pu = derive([Valley(s, "umi", v, True) for s, v in zip(S, UMI)], "umi", light_floor=200, assay=SN)
print(pu)
pg = derive([Valley(s, "genes", v, True) for s, v in zip(S, GENES)], "genes", assay=SN)
print(pg)
if pu.constant != 350 or pg.constant != 260:
    fails.append(f"A: expected 350/260, got {pu.constant}/{pg.constant}")
print(f" -> the calibration cohort applied 350/250; derived {pu.constant}/{pg.constant} from the same valleys")

print("\nB. the gene lower bound must NOT be the UMI one")
print(f" smallest measured gene valley is {min(GENES)}, below the UMI lower bound "
      f"{UMI_BOUNDS[SN][0]}")
if min(GENES) >= UMI_BOUNDS[SN][0]:
    fails.append("B: premise wrong")
print(" applying UMI bounds to genes would refuse a real library for being correct")

print("\n" + "-" * 74)
print("C. ONE unimodal library - classified, not refused")
# This asserted a refusal until 2026-08-10. A shallow minimum says the number is a judgement
# rather than a measurement; it does not say the number is unusable, and refusing discarded the
# information the depth test had just produced. On a real ten-library cohort six were shoulders,
# all six landed inside the bounds, and the constant they produced was the one that cohort had
# applied. What the depth test decides now is the PROVENANCE label.
try:
    p_mixed = derive([Valley(s, "umi", v, s != "ctrl_03") for s, v in zip(S, UMI)], "umi", assay=SN)
    print(f"[ACCEPTED] constant {p_mixed.constant}, provenance {p_mixed.provenance!r}, "
          f"shoulders {p_mixed.shoulders}")
    if p_mixed.provenance != "declared_informed":
        fails.append(f"C: one shoulder must class as declared_informed, got {p_mixed.provenance}")
    if p_mixed.shoulders != ("ctrl_03",):
        fails.append(f"C: the shoulder library must be NAMED, got {p_mixed.shoulders}")
    if not any("NOT a pure measurement" in n for n in p_mixed.notes):
        fails.append("C: the proposal must say it is not a pure measurement")
except ThresholdRefusal as e:
    fails.append(f"C: one shoulder must NOT refuse - {str(e)[:100]}")

print("\nC2. EVERY library unimodal - nothing to take a median OF")
try:
    derive([Valley(s, "umi", v, False) for s, v in zip(S, UMI)], "umi", assay=SN)
    fails.append("C2: an all-unimodal cohort must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:180]}...")

print("\nD. a valley above the UMI upper bound")
try:
    bad = UMI.copy(); bad[0] = 1240
    derive([Valley(s, "umi", v, True) for s, v in zip(S, bad)], "umi", assay=SN)
    fails.append("D: >1000 must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:175]}...")

print("\nE. a valley below the UMI lower bound")
try:
    bad = UMI.copy(); bad[4] = 150
    derive([Valley(s, "umi", v, True) for s, v in zip(S, bad)], "umi", assay=SN)
    fails.append("E: <200 must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:150]}...")

print("\nF. a gene valley above 600")
try:
    bad = GENES.copy(); bad[2] = 640
    derive([Valley(s, "genes", v, True) for s, v in zip(S, bad)], "genes", assay=SN)
    fails.append("F: genes >600 must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:140]}...")

print("\nG. the constant must clear the light floor")
try:
    low = [210] * 10
    derive([Valley(s, "umi", v, True) for s, v in zip(S, low)], "umi", light_floor=250, assay=SN)
    fails.append("G: constant at/below the light floor must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:165]}...")

print("\n" + "-" * 74)
print("I. single cells are bounded as single cells, never as nuclei")
# THE DEFECT THIS SECTION PINS: one pair of bounds, read off nuclear valleys, applied to every
# assay. A whole cell carries about ten times a nucleus's RNA; its valley has no reason to sit
# where a nucleus's does, and before 2026-10-02 nothing let it sit anywhere else.
if count_bounds("umi", "scrna") == count_bounds("umi", "snrna"):
    fails.append("I: single cells and nuclei share the UMI bounds")
if count_bounds("genes", "scrna") == count_bounds("genes", "snrna"):
    fails.append("I: single cells and nuclei share the gene bounds")
if count_bounds("umi", "scrna") != (500, 2000) or count_bounds("genes", "scrna") != (250, 1500):
    fails.append(f"I: the PI's declared single-cell bounds are not what the module holds: "
                 f"{count_bounds('umi', 'scrna')} / {count_bounds('genes', 'scrna')}")
# A valley of 1,200 UMI is past the nuclear ceiling and inside the whole-cell range: the same
# number must refuse on nuclei and be proposed on cells.
cells = [1180, 1240, 1100, 1310, 1150, 1205, 1190, 1260, 1220, 1170]
try:
    derive([Valley(s, "umi", v, True) for s, v in zip(S, cells)], "umi", light_floor=200, assay=SN)
    fails.append("I: a 1,200-UMI valley was accepted as a NUCLEAR floor")
except ThresholdRefusal:
    pass
try:
    pc = derive([Valley(s, "umi", v, True) for s, v in zip(S, cells)], "umi", light_floor=200,
                assay="scrna")
    print(f" whole cells, valleys 1100-1310: proposed {pc.constant} inside {pc.bounds}")
    if pc.bounds != (500, 2000):
        fails.append(f"I: a single-cell proposal carried bounds {pc.bounds}")
except ThresholdRefusal as e:
    fails.append(f"I: a 1,200-UMI whole-cell valley was refused: {str(e)[:120]}")
# ...and the nuclear calibration cohort's own valleys sit BELOW the whole-cell floor, so lending
# either assay's bounds to the other refuses a correct library.
try:
    derive([Valley(s, "umi", v, True) for s, v in zip(S, UMI)], "umi", assay="scrna")
    fails.append("I: nuclear valleys (274-473) were accepted under single-cell bounds")
except ThresholdRefusal as e:
    print(f" [REFUSED, correctly] nuclear valleys under whole-cell bounds: {str(e)[:90]}...")
# No default assay: the old single pair WAS the nuclear pair, and a default is how it returns.
try:
    derive([Valley(s, "umi", v, True) for s, v in zip(S, UMI)], "umi", assay="")
    fails.append("I: an undeclared assay was bounded anyway")
except ThresholdRefusal:
    pass
try:
    derive([Valley(s, "umi", v, True) for s, v in zip(S, UMI)], "umi")
    fails.append("I: derive() accepted a call with no assay at all")
except TypeError:
    pass

print("\n" + "-" * 74)
print("H. the mitochondrial ceiling")
print(" " + mito_ceiling_note("snrna")[:300] + "...")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: A derives 350 from the same ten valleys the cohort measured - the procedure")
print("ships and the number does not, because those valleys span 274-473 within ONE cohort")
