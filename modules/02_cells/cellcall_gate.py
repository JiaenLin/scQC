# This module gates the cell call: it compares the aligner's and the denoiser's calls and returns
# findings. It removes nothing and subsets nothing; it exists to refuse a cell-call boundary,
# never to apply one.
"""Step 2 - the cell-call gate: CellBender must not be stricter than the aligner.

WHY THE GATE POINTS THIS WAY

Step 2 gates the first boundary at which observations could be discarded. The gate itself removes
nothing - step 7 is the only code path that drops an observation - but the boundary it inspects is
inherited from the tools, and a mistake in it is not symmetric. A nucleus wrongly called empty
never reaches the steps that would have judged it. A droplet wrongly called a cell still has to
pass the count floors, the mitochondrial ceiling, the doublet call and the cluster check. So
permissive is the safe error, and the gate enforces exactly that asymmetry.

It also enforces a boundary set deliberately: CellBender is a DENOISER, not the cell caller. Its
cell/empty split exists here for one reason - scDblFinder requires input with "empty drops having
already been removed". The moment CellBender starts calling FEWER cells than the aligner it has
stopped being a boundary and become a filter, and a cell-selection decision has migrated into a
tool chosen for denoising.

THRESHOLDS, CALIBRATED ON A REAL COHORT RATHER THAN CHOSEN

Measured across the calibration cohort's ten libraries (CALIBRATION.md), CeleScope EmptyDrops_CR
against CellBender:

    cells called, CellBender / aligner    1.22 - 1.84    cohort 1.60
    aligner cells LOST by CellBender      0.0% - 5.0%    cohort 0.4%, 617 of 152,653

Two things follow, and they are why the gate has two tiers rather than one:

  * A 10% REFUSE line leaves the whole accepted cohort clear by a factor of two. Nothing in that
    cohort would ever have tested it, and a gate that has never fired on the only data available
    is a gate nobody knows works.
  * A 5% REVIEW line puts exactly one library on it: treat_03, at 5.0% - the same library whose
    two callers disagree over 6,183 droplets and which therefore needed a written decision
    anyway. Flagging the sample that already needed one is the gate behaving correctly, not a
    false positive.

So: REVIEW at 5%, REFUSE at 10%. The first is where a human looks; the second is where the run
stops. Neither number is round by accident - 5% is where that cohort's worst library sits, and
10% is twice it. Both are DECLARED rather than derived, so they carry to another dataset as a
starting point and not as a measurement.

THE THIRD CHECK, WHICH THE THRESHOLDS DO NOT COVER

A loss can be small and still be structured. In the calibration cohort all 617 lost cells fall on
one level of the condition factor and none on the other - 0.4% overall, and entirely on one arm of
the design. No per-sample percentage catches that; it is only visible when the losses are grouped
by the design. The removal checklist, question 3 (docs/PRINCIPLES.md), applies to the cell call
exactly as it applies to a filter.

REPORTING AND ENFORCING ARE SEPARATE CALLS

`gate()` returns findings and `verdict()` summarises them as a string; neither stops anything.
`enforce()` raises GateRefusal on a REFUSE verdict, and exists because a returned string is
ignorable by accident where an exception is not. Use it when a refusal must halt the caller;
the bundled CLI instead prints the findings and exits 2, which is the same decision made
visible. What is not safe is to read the verdict and carry on regardless.
"""

from __future__ import annotations

from dataclasses import dataclass, field

REVIEW_LOST = 0.05 # aligner cells lost by CellBender
REFUSE_LOST = 0.10
REFUSE_RATIO = 1.00 # CellBender calling fewer cells than the aligner
DESIGN_REFUSE = 3.0 # the design-differential check
# Binary floating point puts an exact 3x differential on the wrong side of the line: 0.30/0.10
# evaluates to 2.9999999999999996, so an exactly-3x loss - the one case the threshold names -
# would not fire. The comparison therefore carries a relative tolerance.
RATIO_TOL = 1e-9
# A RATIO TEST IS UNRELIABLE AT BOTH EXTREMES. Above about 33% loss in the worst arm the check
# cannot fail at all, because three times the lowest arm already exceeds 100%. Near zero it fails
# the other way: a ratio between two near-zero rates is dominated by a single library, and
# refusing a run over half a percentage point of compositional distortion is how a gate gets
# switched off. So the differential refuses only when the loss is also MATERIAL. Below this floor
# a >=3x differential is still reported - it is informative - but as REVIEW, and it does not stop
# the run. That is deliberate: a real >=3x differential passes as REVIEW whenever the worst arm
# loses under 1% of its cells.
MATERIAL_LOSS = 0.01 # worst arm must lose at least this fraction for the ratio to bind

class GateRefusal(RuntimeError):
    """Raised when the cell-call gate refuses the run."""

@dataclass
class GateFinding:
    check: str
    severity: str
    message: str
    detail: list = field(default_factory=list)

    def __str__(self) -> str:
        s = f"[{self.severity:6s}] {self.check}\n {self.message}"
        for d in self.detail:
            s += f"\n - {d}"
        return s

def _ratio_binds(hi: float, lo: float, threshold: float = DESIGN_REFUSE) -> bool:
    """True when `hi` is at least `threshold` times `lo`, the exact boundary included.

    Written as a cross-multiplication with a relative tolerance rather than `hi / lo >=
    threshold`: the division form misses an exactly-threshold differential, which is the single
    case the threshold is named for. `hi == 0` means no loss anywhere and never binds.
    """
    return hi > 0 and hi + abs(hi) * RATIO_TOL >= threshold * lo

def gate(calls, design=None) -> list:
    """calls: {sample: {"aligner": n, "cellbender": n, "lost": n}} - `lost` is aligner cells
    NOT called by CellBender. design: {factor: {sample: level}}."""
    out = []

    strict = {s: v["cellbender"] / max(v["aligner"], 1) for s, v in calls.items()}
    bad = {s: r for s, r in strict.items() if r < REFUSE_RATIO}
    out.append(GateFinding(
        "CellBender must not be stricter than the aligner",
        "REFUSE" if bad else "ok",
        (f"{len(bad)} sample(s) call FEWER cells than the aligner - the boundary has become a "
         f"filter, and a cell-selection decision has migrated into a denoiser"
         if bad else
         f"all {len(strict)} call more (ratio {min(strict.values()):.2f}-"
         f"{max(strict.values()):.2f})"),
        [f"{s}: ratio {r:.2f}" for s, r in sorted(bad.items())]))

    lost = {s: v["lost"] / max(v["aligner"], 1) for s, v in calls.items()}
    refuse = {s: f for s, f in lost.items() if f > REFUSE_LOST}
    review = {s: f for s, f in lost.items() if REVIEW_LOST <= f <= REFUSE_LOST}
    sev = "REFUSE" if refuse else ("REVIEW" if review else "ok")
    out.append(GateFinding(
        "aligner cells lost by CellBender",
        sev,
        (f"worst {100*max(lost.values()):.1f}%; REVIEW at {100*REVIEW_LOST:.0f}%, "
         f"REFUSE at {100*REFUSE_LOST:.0f}%"),
        [f"{s}: {100*f:.1f}% ({calls[s]['lost']:,} of {calls[s]['aligner']:,}) "
         f"{'REFUSE' if f > REFUSE_LOST else 'review'}"
         for s, f in sorted({**refuse, **review}.items(), key=lambda x: -x[1])]))

    # Small losses can still be one-sided. A percentage per sample cannot see this.
    if design:
        for factor, mapping in design.items():
            by = {}
            for s, v in calls.items():
                lvl = mapping.get(s)
                if lvl is None:
                    continue
                a, l = by.setdefault(lvl, [0, 0])
                by[lvl] = [a + v["aligner"], l + v["lost"]]
            if len(by) < 2:
                continue
            rate = {k: l / max(a, 1) for k, (a, l) in by.items()}
            hi, lo = max(rate.values()), min(rate.values())
            txt = " · ".join(f"{k} {100*v:.2f}%" for k, v in sorted(rate.items()))
            if lo == 0 and hi > 0:
                out.append(GateFinding(
                    f"loss differential: {factor}", "REVIEW",
                    f"all loss falls on one level and none on another ({txt}). The ratio is "
                    f"undefined rather than large, which is why a threshold on it would not "
                    f"have fired - but a one-sided loss on a design factor is exactly what "
                    f"the design-differential check exists to catch, however small"))
            else:
                r = hi / max(lo, 1e-12)
                binds = _ratio_binds(hi, lo)
                material = hi >= MATERIAL_LOSS
                if binds and material:
                    sev, note = "REFUSE", (
                        f" - REFUSED: at or above the {DESIGN_REFUSE:.0f}x design-differential "
                        f"line, and the loss is material - the worst arm loses "
                        f"{100*hi:.2f}% of its cells")
                elif binds:
                    sev, note = "REVIEW", (
                        f" ratio reaches {DESIGN_REFUSE:.0f}x but the loss is NOT material - "
                        f"the worst arm loses {100*hi:.2f}%, under the {100*MATERIAL_LOSS:.0f}% "
                        f"floor. A ratio between two near-zero rates is dominated by single "
                        f"libraries; reported, not refused")
                else:
                    sev, note = "ok", ""
                out.append(GateFinding(
                    f"loss differential: {factor}", sev,
                    f"max/min = {r:.2f}x ({txt}){note}"))
    else:
        out.append(GateFinding("loss differential", "REVIEW",
                               "no design given - a one-sided loss cannot be detected"))
    return out

def verdict(findings) -> str:
    """REPORT the findings as one severity - "REFUSE", "REVIEW" or "PASS" - and raise nothing.

    This is the reporting form, for printing and for tables. It is not the enforcing form: a
    caller that acts on the cell call uses enforce(), because a returned "REFUSE" can be ignored
    by accident.
    """
    if any(f.severity == "REFUSE" for f in findings):
        return "REFUSE"
    return "REVIEW" if any(f.severity == "REVIEW" for f in findings) else "PASS"

def enforce(findings) -> str:
    """ENFORCE the verdict: raise GateRefusal on REFUSE, otherwise return "REVIEW" or "PASS".

    Every gate here refuses rather than warns, and this is where step 2 does it. REVIEW is
    returned rather than raised on purpose - it means a human must look, not that the run is
    wrong - so a caller that wants to stop on REVIEW as well checks the returned value.
    """
    v = verdict(findings)
    if v == "REFUSE":
        raise GateRefusal(
            "the cell-call gate REFUSES this run:\n"
            + "\n".join(str(f) for f in findings if f.severity == "REFUSE"))
    return v
