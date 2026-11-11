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

They were read off a calibration cohort's libraries (CALIBRATION.md says how; the measured values
are that cohort's and live in its own records), comparing the aligner's cell call (CeleScope
EmptyDrops_CR) against CellBender's: the ratio of cells called, and the share of the aligner's
cells that CellBender LOST. Two things follow, and they are why the gate has two tiers rather
than one:

  * A 10% REFUSE line left the whole accepted cohort clear. Nothing in that cohort would ever have
    tested it, and a gate that has never fired on the only data available is a gate nobody knows
    works.
  * A 5% REVIEW line put on it the one library whose two callers disagreed over many droplets
    and which therefore needed a written decision anyway. Flagging the sample that already needed one is the gate behaving correctly, not a
    false positive.

So: REVIEW at 5%, REFUSE at 10%. The first is where a human looks; the second is where the run
stops. Neither number is round by accident - 5% is where that cohort's worst library sits, and
10% is twice it. Both are DECLARED rather than derived, so they carry to another dataset as a
starting point and not as a measurement.

THE THIRD CHECK, WHICH THE THRESHOLDS DO NOT COVER

A loss can be small and still be structured. In the calibration cohort every lost cell fell on
one level of the condition factor and none on the other - well under 1% overall, and entirely on
one arm of the design. No per-sample percentage catches that; it is only visible when the losses are grouped
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

# PER ASSAY (harness autodev W020, B1-074): the lines below were measured on nuclei (the header
# above), and a whole cell is not one. The same values on both assays, with LINES_SOURCE saying so.
REVIEW_LOST = {"snrna": 0.05, "scrna": 0.05}   # aligner cells lost by CellBender
REFUSE_LOST = {"snrna": 0.10, "scrna": 0.10}
REFUSE_RATIO = {"snrna": 1.00, "scrna": 1.00}  # CellBender calling fewer cells than the aligner
DESIGN_REFUSE = {"snrna": 3.0, "scrna": 3.0}   # the design-differential check
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
MATERIAL_LOSS = {"snrna": 0.01, "scrna": 0.01}  # worst arm must lose at least this for the ratio to bind


#: Where each per-assay line came from. The values are the same on both assays today: every one
#: was set or measured on the calibration cohort's NUCLEI, and lending it to cells is done here
#: knowingly, in a table that says so, never by a default (harness autodev W020, B1-074; the PI,
#: 2026-10-02: single cells and single nuclei are fundamentally different).
LINES_SOURCE = {"snrna": "set or measured on the calibration cohort's nuclei",
                "scrna": "calibrated on nuclei, unmeasured on cells - the same value, lent knowingly "
                         "until a single-cell cohort measures its own"}


def per_assay(table: dict, assay, what: str):
    """The assay's entry of a per-assay table. No default: an unknown assay is refused by name."""
    a = str(assay or "").strip().lower()
    if a not in table:
        raise ValueError(f"{what}: no value for assay {assay!r} (known: {sorted(table)}). The "
                         f"lines are declared per assay, and one assay's is not lent to another.")
    return table[a]

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

def _ratio_binds(hi: float, lo: float, threshold: float) -> bool:
    """True when `hi` is at least `threshold` times `lo`, the exact boundary included.

    Written as a cross-multiplication with a relative tolerance rather than `hi / lo >=
    threshold`: the division form misses an exactly-threshold differential, which is the single
    case the threshold is named for. `hi == 0` means no loss anywhere and never binds.
    """
    return hi > 0 and hi + abs(hi) * RATIO_TOL >= threshold * lo

def gate(calls, design=None, *, assay) -> list:
    """calls: {sample: {"aligner": n, "cellbender": n, "lost": n}} - `lost` is aligner cells
    NOT called by CellBender. design: {factor: {sample: level}}. `assay` selects the lines
    (per assay above); it is required and has no default."""
    out = []
    review_lost = per_assay(REVIEW_LOST, assay, "REVIEW_LOST")
    refuse_lost = per_assay(REFUSE_LOST, assay, "REFUSE_LOST")
    refuse_ratio = per_assay(REFUSE_RATIO, assay, "REFUSE_RATIO")
    design_refuse = per_assay(DESIGN_REFUSE, assay, "DESIGN_REFUSE")
    material_loss = per_assay(MATERIAL_LOSS, assay, "MATERIAL_LOSS")

    strict = {s: v["cellbender"] / max(v["aligner"], 1) for s, v in calls.items()}
    bad = {s: r for s, r in strict.items() if r < refuse_ratio}
    lost = {s: v["lost"] / max(v["aligner"], 1) for s, v in calls.items()}
    # STRICTER IS A FILTER ONLY WHEN IT COSTS CELLS (harness autodev W020, B1-103). A denoiser that
    # calls fewer cells than the aligner while losing under review_lost of the aligner's own cells
    # has dropped droplets the aligner called and kept nearly every cell - the loss check below
    # judges exactly that, and refusing on the ratio alone stopped runs the loss check passed.
    # Such a library is REVIEW; one that is stricter AND loses at least review_lost is REFUSED.
    material = {s for s in bad if lost[s] >= review_lost}
    out.append(GateFinding(
        "CellBender must not be stricter than the aligner",
        "REFUSE" if material else ("REVIEW" if bad else "ok"),
        (f"{len(material)} sample(s) call FEWER cells than the aligner and lose at least "
         f"{100*review_lost:.0f}% of the aligner's cells - the boundary has become a filter, and "
         f"a cell-selection decision has migrated into a denoiser"
         if material else
         f"{len(bad)} sample(s) call fewer cells than the aligner, but each loses under "
         f"{100*review_lost:.0f}% of the aligner's cells - reported, not refused"
         if bad else
         f"all {len(strict)} call more (ratio {min(strict.values()):.2f}-"
         f"{max(strict.values()):.2f})"),
        [f"{s}: ratio {r:.2f}, lost {100*lost[s]:.1f}%" for s, r in sorted(bad.items())]))

    refuse = {s: f for s, f in lost.items() if f > refuse_lost}
    review = {s: f for s, f in lost.items() if review_lost <= f <= refuse_lost}
    sev = "REFUSE" if refuse else ("REVIEW" if review else "ok")
    out.append(GateFinding(
        "aligner cells lost by CellBender",
        sev,
        (f"worst {100*max(lost.values()):.1f}%; REVIEW at {100*review_lost:.0f}%, "
         f"REFUSE at {100*refuse_lost:.0f}%"),
        [f"{s}: {100*f:.1f}% ({calls[s]['lost']:,} of {calls[s]['aligner']:,}) "
         f"{'REFUSE' if f > refuse_lost else 'review'}"
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
                binds = _ratio_binds(hi, lo, design_refuse)
                material = hi >= material_loss
                if binds and material:
                    sev, note = "REFUSE", (
                        f" - REFUSED: at or above the {design_refuse:.0f}x design-differential "
                        f"line, and the loss is material - the worst arm loses "
                        f"{100*hi:.2f}% of its cells")
                elif binds:
                    sev, note = "REVIEW", (
                        f" ratio reaches {design_refuse:.0f}x but the loss is NOT material - "
                        f"the worst arm loses {100*hi:.2f}%, under the {100*material_loss:.0f}% "
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
