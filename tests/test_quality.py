# Derives count-floor proposals from recorded valleys and prints them; it applies no threshold.
"""Step 5 test: real density valleys, the bounds around them, and the cases that must refuse."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "05_quality"))
from quality import (Valley, ThresholdRefusal, derive, # noqa: E402
                     mito_ceiling_note, GENE_BOUNDS, UMI_BOUNDS)

S = ["ctrl_01", "ctrl_02", "treat_01", "treat_02", "ctrl_03",
     "ctrl_04", "ctrl_05", "treat_03", "treat_04", "treat_05"]
# the calibration cohort measured: UMI median 348 (274-473), genes median 264 (184-352)
UMI = [473, 300, 361, 348, 274, 336, 352, 289, 383, 344]
GENES = [352, 214, 271, 264, 184, 249, 268, 208, 291, 258]

fails = []
print("Step 5 - quality thresholds\n" + "=" * 74)
print(f"bounds: UMI {UMI_BOUNDS}, genes {GENE_BOUNDS}\n")

print("A. the calibration cohort's real valleys")
pu = derive([Valley(s, "umi", v, True) for s, v in zip(S, UMI)], "umi", light_floor=200)
print(pu)
pg = derive([Valley(s, "genes", v, True) for s, v in zip(S, GENES)], "genes")
print(pg)
if pu.constant != 350 or pg.constant != 260:
    fails.append(f"A: expected 350/260, got {pu.constant}/{pg.constant}")
print(f" -> the calibration cohort applied 350/250; derived {pu.constant}/{pg.constant} from the same valleys")

print("\nB. the gene lower bound must NOT be the UMI one")
print(f" smallest measured gene valley is {min(GENES)}, below the UMI lower bound "
      f"{UMI_BOUNDS[0]}")
if min(GENES) >= UMI_BOUNDS[0]:
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
    p_mixed = derive([Valley(s, "umi", v, s != "ctrl_03") for s, v in zip(S, UMI)], "umi")
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
    derive([Valley(s, "umi", v, False) for s, v in zip(S, UMI)], "umi")
    fails.append("C2: an all-unimodal cohort must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:180]}...")

print("\nD. a valley above the UMI upper bound")
try:
    bad = UMI.copy(); bad[0] = 1240
    derive([Valley(s, "umi", v, True) for s, v in zip(S, bad)], "umi")
    fails.append("D: >1000 must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:175]}...")

print("\nE. a valley below the UMI lower bound")
try:
    bad = UMI.copy(); bad[4] = 150
    derive([Valley(s, "umi", v, True) for s, v in zip(S, bad)], "umi")
    fails.append("E: <200 must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:150]}...")

print("\nF. a gene valley above 600")
try:
    bad = GENES.copy(); bad[2] = 640
    derive([Valley(s, "genes", v, True) for s, v in zip(S, bad)], "genes")
    fails.append("F: genes >600 must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:140]}...")

print("\nG. the constant must clear the light floor")
try:
    low = [210] * 10
    derive([Valley(s, "umi", v, True) for s, v in zip(S, low)], "umi", light_floor=250)
    fails.append("G: constant at/below the light floor must refuse")
except ThresholdRefusal as e:
    print(f"[REFUSED] {str(e)[:165]}...")

print("\n" + "-" * 74)
print("H. the mitochondrial ceiling")
print(" " + mito_ceiling_note()[:300] + "...")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: A derives 350 from the same ten valleys the cohort measured - the procedure")
print("ships and the number does not, because those valleys span 274-473 within ONE cohort")
