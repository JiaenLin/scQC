# This module removes nothing: it selects which nuclei are handed to the doublet detectors and
# records the rest as UNSCORED, not as negative. Every nucleus below the floor continues to
# step 5 and is judged there. The masks here choose an input set, not an output.
"""Step 3 - the light floor: which nuclei get scored for doublets, and the bookkeeping that must
travel with that choice.

WHAT IT IS AND IS NOT

A minimum UMI applied ONLY to select the doublet detectors' input. It is not a quality filter.
Nothing is removed from the analysis by it.

Doublet detectors build their null by SUMMING PAIRS OF OBSERVED TRANSCRIPTOMES. If the pool they
sample from contains near-empty droplets, the artificial "doublets" are debris+cell or
debris+debris, which are not doublets, and the null is calibrated on the wrong thing. Only one of
the tools puts a number on it:

    scDblFinder: "it might be necessary to remove cells with a very low coverage
                  (e.g. <200 reads) to avoid errors."

The same vignette records the hard failure behind it - a 'Size factors should be positive' error
"if you have some cells that have zero reads (or a very low read count, leading to zero after
feature selection)". So the floor is partly numerical stability, not statistical taste.

DoubletDetection and scds document no such requirement. For them the floor is a reasonable
assumption and is recorded as one.

THE ORDERING IS THE POINT

The light floor must sit BELOW the applied count floor. Filtering hard before doublet detection
is what scDblFinder's FAQ warns against - "further quality filtering should be performed
downstream of doublet detection" - and a UMI CEILING is what doi:10.1101/140848 used AS its
doublet filter. Aggressive pre-filtering does not merely lose cells; it deletes the population
the next step exists to find. So a light floor that creeps up to or past the quality floor has
silently become a quality filter, and this module refuses it.

TWO FAILURES THIS STEP CAUSES IF NOBODY WATCHES IT

  1  "NOT CALLED" IS NOT "NOT EXAMINED". Colour the never-scored nuclei the same grey as nuclei
     the tools looked at and cleared, and absence of a call reads as evidence of absence. In the
     calibration cohort that is 59,322 nuclei, 24% of the set. The unscored set is therefore
     carried as its own state and is never folded into the negatives.

  2  A RATE WITHOUT ITS DENOMINATOR IS NOT COMPARABLE. In the calibration cohort the doublet
     rate is 5.70% against all denoised cells and 7.52% against the nuclei actually scored,
     because 24% of the first denominator had never been examined. Both numbers are true;
     quoting one against a published band without saying which it is, is not. So a rate
     computed downstream carries its denominator with it and cannot be printed without one.
"""

from __future__ import annotations

from dataclasses import dataclass, field

DEFAULT_FLOOR = 200 # scDblFinder's documented value; the only one of the tools to give one
DOCUMENTED_BY = ("scDblFinder",)
ASSUMED_FOR = ("DoubletDetection", "scds")

class FloorRefusal(RuntimeError):
    """Raised when the light floor has become a quality filter."""

@dataclass
class Coverage:
    """Who was scored, who was not, and what that permits saying."""
    n_total: int
    n_scored: int
    floor: int
    max_unscored_umi: float = None
    min_scored_umi: float = None
    n_kept_unscored: int = None # retained nuclei that were never examined
    notes: list = field(default_factory=list)

    @property
    def n_unscored(self) -> int:
        return self.n_total - self.n_scored

    @property
    def pct_scored(self) -> float:
        return 100 * self.n_scored / max(self.n_total, 1)

    def rate(self, n_called: int, denominator: str = "scored") -> str:
        """A doublet rate ALWAYS carries its denominator. There is no bare-number form."""
        if denominator == "scored":
            return (f"{100 * n_called / max(self.n_scored, 1):.2f}% of nuclei SCORED "
                    f"({n_called:,} of {self.n_scored:,})")
        if denominator == "all":
            return (f"{100 * n_called / max(self.n_total, 1):.2f}% of ALL cells "
                    f"({n_called:,} of {self.n_total:,}), including "
                    f"{self.n_unscored:,} never examined")
        raise ValueError("denominator must be 'scored' or 'all' - a rate without one is not "
                         "comparable to anything, including a published figure")

    def __str__(self) -> str:
        s = [f"scored {self.n_scored:,} of {self.n_total:,} ({self.pct_scored:.2f}%) "
             f"at a {self.floor}-UMI floor",
             f" unscored: {self.n_unscored:,} - NOT negative, never examined"]
        if self.max_unscored_umi is not None and self.min_scored_umi is not None:
            clean = self.min_scored_umi >= self.floor > self.max_unscored_umi
            s.append(f" boundary: unscored max {self.max_unscored_umi:,.0f}, scored min "
                     f"{self.min_scored_umi:,.0f} - {'clean' if clean else 'NOT CLEAN'}")
        if self.n_kept_unscored is not None:
            s.append(f" retained nuclei never scored: {self.n_kept_unscored:,}"
                     + (" (the coverage gap does not reach the deliverable)"
                        if self.n_kept_unscored == 0 else
                        " <- these carry a doublet status of UNKNOWN, not negative"))
        s += [f" {n}" for n in self.notes]
        return "\n".join(s)

def check_floor(floor: int, quality_floor: int) -> list:
    """The light floor must stay below the quality floor, or it has become one."""
    notes = []
    if floor >= quality_floor:
        raise FloorRefusal(
            f"light floor {floor} is not below the quality floor {quality_floor}. It has stopped "
            f"selecting a scoring set and become a quality filter applied BEFORE doublet "
            f"detection - which is what scDblFinder's FAQ warns against, and a UMI cut is what "
            f"doi:10.1101/140848 used AS a doublet filter. Lower it, or move the quality floor.")
    if floor != DEFAULT_FLOOR:
        notes.append(f"floor {floor} departs from the documented {DEFAULT_FLOOR} "
                     f"({', '.join(DOCUMENTED_BY)}); state why")
    notes.append(f"documented by {', '.join(DOCUMENTED_BY)}; an assumption for "
                 f"{', '.join(ASSUMED_FOR)}")
    return notes

def assess(n_total, n_scored, floor, quality_floor, max_unscored_umi=None,
           min_scored_umi=None, n_kept_unscored=None) -> Coverage:
    notes = check_floor(floor, quality_floor)
    c = Coverage(n_total, n_scored, floor, max_unscored_umi, min_scored_umi,
                 n_kept_unscored, notes)
    if min_scored_umi is not None and min_scored_umi < floor:
        raise FloorRefusal(
            f"a nucleus with {min_scored_umi:,.0f} UMI was scored but the floor is {floor} - "
            f"the scored set is not the set the floor defines, so its coverage is unknown")
    return c
