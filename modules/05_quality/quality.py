# This module derives quality thresholds and refuses implausible ones.
# It applies nothing; the removal happens at step 7 under a recorded approval.
"""Step 5 - quality thresholds: the count floors are MEASURED, the mitochondrial ceiling is not.

THE CLEAREST CASE IN THE PIPELINE OF A PROCEDURE TRANSFERRING WHERE A NUMBER CANNOT

Across ten libraries of ONE cohort - one sample type, one platform, one operator - the measured
density valley between the debris mode and the nucleus mode ranged:

    UMI              274 - 473   (median 346)
    detected genes   184 - 352   (median 261)

A 1.7-fold spread inside a single experiment. A fixed floor of 500 UMI and 300 genes - a common
published choice, and the obvious one to reach for - sits ABOVE the measured valley in all ten of
those libraries on UMI and in nine of ten on genes, so it cuts into the nucleus mode rather than
the debris below it. No fixed floor is right for most libraries even within one study, which is
why what ships is the KDE-valley PROCEDURE and not the pair of numbers it produced here.

BIMODALITY IS TESTED, NOT ASSUMED

A density minimum exists in any smooth curve. If the distribution is unimodal the minimum is not
a valley, it is the flank of the only mode, and returning it as a threshold dresses an arbitrary
cut as a measurement. So bimodality is checked first and the module REFUSES to return a floor
when it fails - which is also why the mitochondrial ceiling cannot be derived this way at all:
its distribution is unimodal, there is no valley, and every cut is a judgement.

THE VALLEY IS MEASURED PER LIBRARY AND APPLIED AS A COHORT CONSTANT

Ten valleys were measured on the calibration cohort and ONE pair of numbers applied, 350/250.
That is deliberate. A per-library threshold makes the filter a technical property that varies
across the design - the same objection that forces one learning rate for all ten libraries at
step 1, and the same shape the design-differential check exists to catch. The per-library
values decide WHAT the constant should be; they are not themselves applied.

THE BOUNDS, AND WHY THEY ARE NOT THE THRESHOLD

`UMI_BOUNDS = (200, 1000)`. These do not choose a floor - the data does. They catch a floor that
cannot be right:

  below 200    the quality floor would sit at or under the doublet-scoring light floor,
               inverting the pipeline's ordering. Step 3 already refuses a light floor that is
               not strictly below the quality floor; this is the same constraint seen from the
               other side.
  above 1000   more than twice the largest valley measured on the calibration cohort (473), and
               twice the 500 already shown to sit above the valley in every one of its
               libraries. A KDE returning 1,200 has found something other than the
               debris/nucleus boundary.

A value inside the bounds is not thereby correct; it is merely not obviously wrong.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import median

UMI_BOUNDS = (200, 1000)
# Genes get their own bounds, 100-600, and NOT the UMI ones. The measured gene valleys ran
# 184-352, so the smallest of them is ALREADY BELOW 200 - applying the UMI lower bound to genes
# would have refused a real library for being correct. The upper bound of 600 is ~1.7x the
# largest measured valley, the same headroom the UMI bound leaves.
GENE_BOUNDS = (100, 600)
SPREAD_REVIEW = 2.0 # per-library valleys differing more than this question a constant

class ThresholdRefusal(RuntimeError):
    """Raised when a derived threshold cannot be trusted or applied."""

@dataclass
class Valley:
    sample: str
    metric: str
    value: float
    bimodal: bool
    note: str = ""

@dataclass
class Proposal:
    metric: str
    per_library: dict
    constant: int
    bounds: tuple
    notes: list = field(default_factory=list)

    @property
    def spread(self) -> float:
        """max/min across libraries, or None when a library has no positive valley.

        A missing or zero valley makes the ratio UNDEFINED. Guarding the denominator with an
        epsilon would report ~1e12 instead, which prints like a measurement and would be read
        as one; undefined has to be reported as undefined.
        """
        v = list(self.per_library.values())
        if not v or any(x is None or x <= 0 for x in v):
            return None
        return max(v) / min(v)

    def __str__(self) -> str:
        v = list(self.per_library.values())
        sp = self.spread
        sp_txt = f"{sp:.2f}x" if sp is not None else "UNDEFINED"
        s = [f"{self.metric}: per-library valley {min(v):.0f}-{max(v):.0f} "
             f"(median {median(v):.0f}, spread {sp_txt})",
             f" PROPOSED cohort constant: {self.constant} "
             f"bounds {self.bounds[0]}-{self.bounds[1]}"]
        s += [f" {n}" for n in self.notes]
        return "\n".join(s)

def derive(valleys, metric, light_floor=None) -> Proposal:
    """Turn per-library valleys into a proposed cohort constant, refusing what cannot be right."""
    bounds = UMI_BOUNDS if metric == "umi" else GENE_BOUNDS

    bad = [v for v in valleys if not v.bimodal]
    if bad:
        raise ThresholdRefusal(
            f"{metric}: {len(bad)} library(ies) are NOT bimodal "
            f"({', '.join(v.sample for v in bad)}). A density minimum exists in any smooth "
            f"curve; without two modes it is the flank of the only one, and returning it as a "
            f"threshold dresses an arbitrary cut as a measurement. Choose the cut explicitly and "
            f"record it as a judgement, as the mitochondrial ceiling is.")

    per = {v.sample: v.value for v in valleys}
    out = [v for s, v in per.items() if not (bounds[0] <= v <= bounds[1])]
    if out:
        who = ", ".join(f"{s} {v:.0f}" for s, v in per.items()
                        if not (bounds[0] <= v <= bounds[1]))
        raise ThresholdRefusal(
            f"{metric}: valley outside {bounds[0]}-{bounds[1]} in {len(out)} library(ies) "
            f"({who}). Below the lower bound the quality floor would sit at or under the "
            f"doublet-scoring light floor and invert the ordering; above the upper bound the KDE "
            f"has found something other than the debris/nucleus boundary. Neither is a threshold "
            f"to apply - look at the density first.")

    constant = int(round(median(list(per.values())) / 10.0) * 10)
    notes = []
    if not (bounds[0] <= constant <= bounds[1]):
        raise ThresholdRefusal(f"{metric}: proposed constant {constant} is outside "
                               f"{bounds[0]}-{bounds[1]}")
    if light_floor is not None and constant <= light_floor:
        raise ThresholdRefusal(
            f"{metric}: proposed constant {constant} is not above the {light_floor}-UMI light "
            f"floor. The quality filter would then apply at or before doublet scoring, which is "
            f"the ordering scDblFinder's documentation forbids.")
    p = Proposal(metric, per, constant, bounds, notes)
    # `p.spread` cannot be None here: the bounds check above has already refused any valley that
    # is not a positive number inside the bounds, so every value is > 0 by this line. The
    # property still handles None because a Proposal can be built directly by a caller that did
    # not come through derive(), and __str__ prints UNDEFINED rather than a fabricated ratio.
    # A branch for it HERE would be unreachable, and unreachable code reads as a handled case.
    if p.spread > SPREAD_REVIEW:
        notes.append(
            f"per-library valleys differ {p.spread:.2f}x - a single constant fits some libraries "
            f"better than others. It is still preferred to a per-library threshold, which would "
            f"make the filter a technical property varying across the design, but the poorest-fit "
            f"library should be looked at")
    notes.append("PROPOSED, not applied. The constant is a cohort decision and needs approval; "
                 "the per-library valleys decide what it should be, they are not applied "
                 "themselves")
    return p

def mito_ceiling_note() -> str:
    """The mitochondrial ceiling cannot be derived by this route, and says so."""
    return (
        "mitochondrial ceiling: NOT DERIVABLE. The distribution is unimodal - there is no valley, "
        "so every cut is a judgement rather than a measurement, and it is also specific to what "
        "was sampled: where the majority population is the large, RNA-rich, mitochondria-dense "
        "one, a ceiling borrowed from a different sample type removes that population and "
        "reports the removal as quality control. It is an "
        "ADJUDICATED parameter: the pipeline emits the per-library distribution, what each "
        "candidate cut removes, the design differential at each, and cluster-level medians so a "
        "reader can see whether a cut removes a POPULATION or trims a tail - then stops.")
