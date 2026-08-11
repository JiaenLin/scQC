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
    """Turn per-library valleys into a proposed cohort constant, refusing what cannot be right.

    A SHALLOW VALLEY CHANGES WHAT THE FLOOR IS CALLED, NOT WHETHER YOU GET ONE.

    This used to REFUSE the moment any library failed the bimodality test, and that was wrong in
    a way worth naming: it discarded the very information it had just computed. The depth test
    answers "is this minimum a measurement or a judgement?", which is a question about the
    PROVENANCE of the number, not about whether the number is usable.

    On a real ten-library cohort, four had a deep valley stable across kernel bandwidths and six
    had a shallow saddle whose minimum still landed at 268-378 UMI - squarely inside the bounds,
    and agreeing with the four to within 100 UMI. Refusing all ten because six were shoulders
    would have stopped a pipeline that had just measured the right answer.

    It is also not what the field does. Published cardiac and lung atlases choose a bounded floor
    and treat the UMI curve as a CONSIDERATION - "excluded barcodes with less than 200 detected
    genes" (doi:10.1101/2020.01.21.20018358); "CellRanger count metrics ... in addition to the
    UMI curves" (doi:10.1016/j.celrep.2024.115091). None derives the floor automatically, and the
    serious per-cell QC happens after clustering, within cell types.

    THE BOUNDS ARE THE GUARD. A valley outside them is still refused, and that is what catches
    the real failure: on the same cohort, two libraries' minima wandered to ~1,040 UMI under a
    narrower kernel, which the upper bound rejects on sight. Depth is not needed to catch those.

    So `bimodal` now classifies rather than gates. Nothing is loosened - `min_valley_depth` is
    unchanged and still splits that cohort at its natural break (0.090 | 0.133) - and the
    libraries that failed it are NAMED in the proposal, so a reader is told the constant rests on
    a shoulder rather than a dip in those.
    """
    bounds = UMI_BOUNDS if metric == "umi" else GENE_BOUNDS

    shoulders = [v.sample for v in valleys if not v.bimodal]
    if shoulders and len(shoulders) == len(valleys):
        # Every library a shoulder is a different claim from some. There is then no library whose
        # minimum was ever shown to separate two modes, so the cohort constant would be a median
        # of numbers none of which was a measurement, presented as if derived.
        raise ThresholdRefusal(
            f"{metric}: NO library is bimodal ({', '.join(shoulders)}). A density minimum exists "
            f"in any smooth curve; without two modes anywhere in the cohort there is nothing to "
            f"take a median OF, and the result would be a judgement wearing a derivation's "
            f"clothes. Choose the cut explicitly and record it as one, as the mitochondrial "
            f"ceiling's bound is.")

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
    # The provenance of the number, stated where a reader will meet it. A constant resting partly
    # on shoulders is usable and is not a measurement, and the difference has to survive into the
    # report - it is the difference between "the data put the floor here" and "we put the floor
    # here, and the data did not object".
    if shoulders:
        notes.append(
            f"DECLARED, informed by the curve - NOT a pure measurement. {len(shoulders)} of "
            f"{len(valleys)} libraries have a shoulder rather than a dip at their minimum "
            f"({', '.join(sorted(shoulders))}), so for those the position is a judgement the "
            f"density did not contradict. The remaining {len(valleys) - len(shoulders)} carry a "
            f"valley deep enough to be called measured, and the constant sits inside "
            f"{bounds[0]}-{bounds[1]} either way. Record it as a judgement, not as a derivation")
    else:
        notes.append(
            f"DERIVED - every one of {len(valleys)} libraries carries a valley deep enough to "
            f"separate two modes")
    notes.append("PROPOSED, not applied. The constant is a cohort decision and needs approval; "
                 "the per-library valleys decide what it should be, they are not applied "
                 "themselves")
    p.shoulders = tuple(sorted(shoulders))
    p.provenance = "declared_informed" if shoulders else "derived"
    return p

def mito_ceiling_note() -> str:
    """What the mitochondrial ceiling is, and what remains a judgement after it is derived."""
    return (
        "mitochondrial ceiling: DERIVED PER LIBRARY, bounded by a DECLARED statement. The "
        "distribution is unimodal, so there is no valley and the count-floor route does not "
        "apply - but 'no valley' does not mean 'no derivation'. The applied fence is "
        "median + k*1.4826*MAD of each library's own distribution, where k is itself DERIVED: it "
        "is the cohort median of the multiple at which Tukey's Q3 + 1.5*IQR sits, rounded to an "
        "integer. Tukey calibrates how far out the fence belongs and is then retained per library "
        "as an INDEPENDENT CROSS-CHECK; it is never applied, "
        "because two robust routes to one number either agree - which is evidence - or disagree, "
        "which is a finding. The BOUND on the fence is declared by the analyst, because it is a "
        "statement about what a nucleus can be and not a property of this cohort. What stays "
        "ADJUDICATED is whether a mitochondria-high POPULATION is damage or a mitochondria-rich "
        "cell type - that needs an identity, so the pipeline emits cluster-level medians and "
        "stops.")


# --- the mitochondrial ceiling ------------------------------------------------------------------
#
# WHY THIS IS PER LIBRARY WHERE THE COUNT FLOORS ARE A COHORT CONSTANT
#
# The docstring above argues against per-library thresholds: they make the filter a technical
# property that varies across the design. That objection is real and it is why the count floors
# ship as one constant. It does not carry to the mitochondrial ceiling, for two reasons, and both
# are checked at runtime rather than asserted here.
#
#   1. The objection is TESTABLE and it is tested. `assess_differential` below measures what
#      fraction each design arm loses and refuses when they diverge. On the calibration cohort
#      the per-library fence gave 1.05-1.16x across every design factor, against a 3x refusal
#      line - i.e. the feared technical-property-varying-with-design did not occur, and if it
#      does occur on another cohort the derivation refuses instead of shipping it.
#   2. Unlike the valley, there is nothing to hold constant. The valley is the SAME physical
#      boundary in every library - debris against nuclei - so one number is a defensible summary
#      of ten estimates of one quantity. Mitochondrial content is not one boundary measured ten
#      times; libraries genuinely differ in it (4x in Q3 on the calibration cohort). Averaging
#      that produces a number describing no library.
#
# WHY A TUKEY FENCE AND NOT A MAD
#
# Both were computed on the calibration cohort and they agree closely. The fence is anchored at
# Q3 and has NO FREE PARAMETER; MAD needs a k, and k from 3 to 12 is a continuum with no
# principled stopping point - the choice of k silently becomes the threshold. Preferring the rule
# with nothing to tune is not a claim that it performs better. It did not, measurably.
#
# WHY THE BOUND IS DECLARED AND NOT DERIVED
#
# A fence fitted to a library whose distribution is pathological will follow it anywhere. The
# bound is the analyst's statement of what a nucleus can be - "below X the filter risks cutting
# real signal; above Y it is not a viable nucleus" - and it must be recorded in the analyst's own
# words, because it is the one part of this that the data cannot supply. It is a GUARD RAIL: if
# it binds in most libraries it has stopped being a policy and become the threshold, and the
# ceiling is then REPORTED AS BOUND-DOMINATED rather than derived, so that a declared number
# never masquerades as a derived one.

# Assay defaults. NOT thresholds - the outer limits within which a derived fence is allowed to
# land. Whole cells retain cytoplasm and legitimately carry more mitochondrial signal than nuclei.
#
# The snRNA lower bound was 5.0 until 2026-08-11. It was raised to 10.0 on a biological argument
# rather than a statistical one: on the calibration cohort, 5.0 let a library's ceiling fall to
# 6.12%, and the nuclei that cut removed were indistinguishable from the ones it kept - median
# depth 0.88-0.96x that of retained nuclei in three of the four affected libraries, i.e. ordinary
# cells, not debris. In heart especially, cardiomyocytes are the most mitochondria-rich cell type
# in the body, so a 6% ceiling preferentially removes the cell type the study is about.
#
# The floor is also what a bound is FOR. Those libraries had low fences because their preps were
# tight, not because their nuclei were biologically cleaner: across that cohort the fence varied
# MORE between two mice of the same group (3.87x) than between the design groups (2.59x), so the
# spread is a technical property. A declared floor is the right instrument against that - "however
# tight this prep looks, an 8% cardiac nucleus is not an outlier" - and it is exactly the statement
# about what a nucleus can be that the bound exists to carry.
MITO_BOUNDS = {"snrna": (10.0, 25.0), "scrna": (10.0, 30.0)}

# THE APPLIED FENCE IS MAD-BASED; TUKEY IS RETAINED AS AN INDEPENDENT SECOND DERIVATION.
#
#     applied     median + MAD_K * 1.4826 * MAD
#     cross-check Q3 + IQR_MULT * IQR
#
# 1.4826 makes the scaled MAD equal sigma for normal data, so MAD_K reads as a z-score.
#
# Why two, and why this one applied. The two estimators do NOT measure the same thing on a skewed
# distribution: Tukey is anchored at Q3 and grows with the middle-50% width, MAD is anchored at the
# median and grows with a symmetric spread. Mitochondrial percentage is bounded at 0 and strongly
# right-skewed, so IQR/(1.4826*MAD) - exactly 1.349 for a normal - ran 1.59 to 2.82 across the
# calibration cohort. That is why no single MAD_K reproduces Tukey: the k that would match it
# tracks each library's skew, from 3.44 to 6.56.
#
# k IS DERIVED, NOT DECLARED. Tukey is the calibration instrument; MAD is the applied estimator.
#
# For each library, solve for the k that would put the MAD fence exactly on the Tukey fence:
#
#     k_i = (Q3_i + mult*IQR_i - median_i) / (1.4826 * MAD_i)
#
# and take the cohort median, rounded to an integer. This is the same shape as the count floors:
# measure per library, propose ONE cohort constant. It makes k a property of the data rather than
# a number someone chose - on the calibration cohort it lands on 4 (per-library 3.44 to 6.56,
# median 4.26), and on a cohort with different skew it will land elsewhere, which is the point.
#
# Why not simply apply Tukey, if Tukey selects k. Because the two do different jobs. Tukey adapts
# without limit to a wide tail - to 25.88% on one calibration library - and adapting that far is a
# statement about that prep, not about what a nucleus can be. The MAD fence at the calibrated k
# tracks the same cohort-level scale while compressing the extremes, so the applied ceiling stays
# inside a biologically defensible range. Tukey decides HOW FAR out the fence should sit; MAD
# decides how much any single library is allowed to differ in getting there.
#
# k is rounded to an integer deliberately. The per-library k values span 1.9x, so the cohort
# summary is not precise to a decimal, and a k of 4.26 would imply a resolution the spread does
# not support.
#
# WHAT IT COSTS, measured, and recorded because it is not visible from the number: threshold
# evenness and removal evenness trade off monotonically, across every rule tested.
#
#     rule                       ceiling spread   removal differential (interaction arm)
#     cohort constant 12.65%          1.00x                 2.77x
#     MAD k=4, bound 10-25            2.05x                 1.75x   <- applied
#     MAD k=5, bound 10-25            2.46x                 1.60x
#     Tukey,   bound 10-25            2.50x                 1.37x
#     Tukey,   unbounded              4.23x                 1.09x
#
# Every step toward a more uniform THRESHOLD costs a less uniform EFFECT. The reason is that the
# quantity being fenced genuinely varies: Q3 spreads 3.91x and IQR 4.75x across those libraries,
# so an estimator whose output spreads less than that is not being more consistent, it is
# under-adapting. `ceiling_spread` reports the first property, the design differential the second,
# and a reader needs both.
MAD_SCALE = 1.4826
IQR_MULT = 1.5
# Sanity limits on the DERIVED k. Not a tuning knob: a cohort whose Tukey-implied k falls outside
# these has a tail shape this calibration does not describe, and the run says so rather than
# quietly using the edge value. k below 2 fences inside the bulk of the distribution; k above 10
# is so permissive that the bound, not the derivation, is doing all the work.
MAD_K_BOUNDS = (2, 10)
# Per-library implied k spreading more than this says no single k fits the cohort - the libraries
# differ in SHAPE, not just in scale, and a cohort constant is then a compromise rather than a
# measurement. Reported, because the calibration cohort itself spans 1.90x.
K_SPREAD_REVIEW = 2.0
# Tukey and the MAD fence differing by more than this on a library the bound did NOT clamp is
# reported. Only unclamped libraries count: two clamped values agree because the bound made them
# agree, and reading that as corroboration is confidence the data did not supply.
CROSS_CHECK_REVIEW = 1.5
# If the bound binds in more than this fraction of libraries it IS the threshold, not a rail - the
# ceiling is reported as bound-dominated rather than derived. Lower- and upper-binding are counted
# SEPARATELY because they are opposite events: an upper clamp removes MORE than the library's own
# fence asks, a lower clamp retains more. Only one of those can delete signal.
BOUND_BINDS_REVIEW = 0.5
# Applied ceilings differing by more than this across libraries are worth a look, by the same
# reasoning as SPREAD_REVIEW for valleys: one filter is then treating samples of one experiment
# very differently, and that difference is technical unless something says otherwise.
CEILING_SPREAD_REVIEW = 3.0


@dataclass
class MitoCeiling:
    """One library's derived ceiling and the quantities it came from.

    `derived` is the APPLIED fence (MAD-based). `tukey` is the independent second derivation,
    carried so a reader can see whether the two routes agreed rather than being told they did.
    """
    sample: str
    n: int
    median: float
    q1: float
    q3: float
    mad: float
    derived: float          # the raw APPLIED fence (MAD route), before the bound
    tukey: float            # the cross-check fence (Q3 + 1.5*IQR), never applied
    ceiling: float          # what would be applied
    clamped: str            # "", "lower" or "upper"

    @property
    def iqr(self) -> float:
        return self.q3 - self.q1

    @property
    def smad(self) -> float:
        """MAD on the same scale as sigma, which is what MAD_K multiplies."""
        return MAD_SCALE * self.mad

    @property
    def skew_ratio(self):
        """IQR / scaled MAD. Exactly 1.349 for a normal distribution; higher means right-skewed.

        This is why the two fences diverge, so it is reported beside them rather than left for a
        reader to reconstruct from the quartiles.
        """
        return (self.iqr / self.smad) if self.smad > 0 else None

    @property
    def cross_check(self):
        """Tukey / applied fence. 1.0 means the two independent routes landed together."""
        return (self.tukey / self.derived) if self.derived > 0 else None


def fence(values, mult: float = IQR_MULT) -> tuple:
    """Q1, Q3 and the raw upper Tukey fence of `values`. No bound applied."""
    v = sorted(float(x) for x in values if x == x and x not in (float("inf"), float("-inf")))
    if len(v) < 4:
        raise ThresholdRefusal(
            f"a Tukey fence needs at least 4 finite values; got {len(v)}. A library with too few "
            f"nuclei to place a quartile is not a library with no mitochondrial contamination.")

    def q(p):
        # Linear interpolation between order statistics - numpy's default, restated so this
        # module has no hard numpy dependency and so the definition is visible rather than
        # inherited. Quartile conventions differ by up to a whole observation on small n, and a
        # threshold that changes with the library that computed it is not reproducible.
        h = (len(v) - 1) * p
        lo = int(h)
        hi = min(lo + 1, len(v) - 1)
        return v[lo] + (h - lo) * (v[hi] - v[lo])

    q1, q3 = q(0.25), q(0.75)
    return q1, q3, q3 + mult * (q3 - q1)


def select_mad_k(stats, mult=IQR_MULT, k_bounds=MAD_K_BOUNDS) -> dict:
    """Derive the MAD multiplier from the Tukey fence. Tukey calibrates; MAD applies.

    # rule-one: no-removal - this reads a per-library summary and returns a number.

    For each library the k that would place the MAD fence exactly on the Tukey fence is
    `(tukey - median) / (1.4826 * MAD)`. The cohort value is the median of those, rounded to an
    integer, because the per-library values do not agree closely enough to justify a decimal.

    Libraries with `MAD == 0` are EXCLUDED from the calibration rather than skipped silently: the
    implied k is undefined there (division by zero), and a library that cannot contribute to the
    calibration must not be able to change it by being counted as some default.
    """
    implied, undefined = {}, []
    for s, st in stats.items():
        mad = float(st.get("mad", 0.0) or 0.0)
        med, q1, q3 = float(st["median"]), float(st["q1"]), float(st["q3"])
        if mad <= 0:
            undefined.append(s)
            continue
        implied[s] = ((q3 + mult * (q3 - q1)) - med) / (MAD_SCALE * mad)
    if not implied:
        raise ThresholdRefusal(
            "no library has a positive MAD, so the Tukey-implied k is undefined everywhere and "
            "there is nothing to calibrate against. A cohort in which more than half of every "
            "library shares one mitochondrial value is not one this fence describes.")

    vals = sorted(implied.values())
    n = len(vals)
    med_k = vals[n // 2] if n % 2 else 0.5 * (vals[n // 2 - 1] + vals[n // 2])
    k = int(round(med_k))
    lo_k, hi_k = k_bounds
    clamped = ""
    if k < lo_k:
        k, clamped = lo_k, "lower"
    elif k > hi_k:
        k, clamped = hi_k, "upper"

    spread = max(vals) / min(vals) if min(vals) > 0 else None
    notes = [f"MAD k DERIVED from the Tukey fence: per-library {min(vals):.2f}-{max(vals):.2f}, "
             f"median {med_k:.2f}, applied k = {k}"]
    if undefined:
        notes.append(f"excluded from the calibration (MAD is zero): {', '.join(sorted(undefined))}")
    if clamped:
        notes.append(
            f"REVIEW - the derived k was clamped at the {clamped} sanity limit {k_bounds}. A "
            f"cohort whose Tukey-implied k falls outside those limits has a tail shape this "
            f"calibration does not describe; the applied number is the limit, not a derivation")
    if spread is not None and spread > K_SPREAD_REVIEW:
        notes.append(
            f"REVIEW - the implied k spans {spread:.2f}x across libraries, so no single k fits "
            f"them all. The libraries differ in the SHAPE of their tail, not only its scale, and "
            f"the cohort k is a compromise. The libraries at the extremes are the ones to look at")
    return {"k": k, "median_k": med_k, "per_library": implied, "spread": spread,
            "undefined": tuple(sorted(undefined)), "clamped": clamped, "notes": notes}


def derive_mito_ceiling(per_library, assay="snrna", bounds=None, mult=IQR_MULT,
                        declared_by=None, k=None) -> dict:
    """Derive one mitochondrial ceiling per library, bounded by a DECLARED statement.

    `per_library` maps sample -> an iterable of that library's per-nucleus mitochondrial
    percentages. It must be measured on the SAME population the doublet caller sees - cells at or
    above the light floor, BEFORE any quality threshold. Deriving a fence on an already-filtered
    set measures the previous fence, not the data: on the calibration cohort that mistake moved
    the ceilings by up to 39% and capped the observed maximum at 47.8% against a true 99.5%.

    `declared_by` is the analyst's own words for why the bound is what it is. It is REQUIRED when
    the bound is not an assay default, because a bound with no recorded reasoning is indis-
    tinguishable later from a number someone tuned until the result looked right.

    In the pipeline the per-nucleus values live in a worker process and are never shipped to the
    engine; `from_quartiles` below is the entry point used there. Both routes share this
    function's bound and refusal logic so the two cannot diverge.
    """
    stats = {}
    for s, vals in per_library.items():
        v = sorted(float(x) for x in vals if x == x)
        q1, q3, _ = fence(v, mult)
        # The median by linear interpolation, matching `fence()` - not v[len//2], which is a
        # different statistic on even n and would make the applied fence disagree with the
        # cross-check for a reason that has nothing to do with the data.
        h = (len(v) - 1) * 0.5
        lo_i = int(h)
        hi_i = min(lo_i + 1, len(v) - 1)
        med = v[lo_i] + (h - lo_i) * (v[hi_i] - v[lo_i])
        dev = sorted(abs(x - med) for x in v)
        mad = dev[lo_i] + (h - lo_i) * (dev[hi_i] - dev[lo_i])
        stats[s] = {"n": len(v), "median": med, "q1": q1, "q3": q3, "mad": mad}
    return derive_mito_ceiling_from_quartiles(
        stats, assay=assay, bounds=bounds, mult=mult, declared_by=declared_by, k=k)


def derive_mito_ceiling_from_quartiles(stats, assay="snrna", bounds=None, mult=IQR_MULT,
                                       declared_by=None, k=None) -> dict:
    """As `derive_mito_ceiling`, from a precomputed per-library summary.

    `stats` maps sample -> {"n", "median", "q1", "q3", "mad"}. This is what the pipeline uses: the
    worker that already has the matrix open computes five numbers per library instead of
    returning a per-nucleus array, so nothing scales with cell count across the task boundary.

    The APPLIED fence is `median + k * 1.4826 * MAD`. Tukey's `Q3 + mult * IQR` is computed from
    the same summary and carried as an independent cross-check; it is never applied.
    """
    if assay not in MITO_BOUNDS and bounds is None:
        raise ThresholdRefusal(
            f"unknown assay {assay!r} and no explicit bounds. Known: {sorted(MITO_BOUNDS)}. "
            f"Guessing a bound for an unknown assay is guessing what a nucleus can be.")
    lo, hi = bounds if bounds is not None else MITO_BOUNDS[assay]
    if bounds is not None and not declared_by:
        raise ThresholdRefusal(
            "a non-default mitochondrial bound requires `declared_by` - the analyst's own words "
            "for why it is what it is. The bound is the one part of this the data cannot supply.")
    if not (0 <= lo < hi <= 100):
        raise ThresholdRefusal(f"mitochondrial bounds must satisfy 0 <= lo < hi <= 100; got "
                               f"{lo}-{hi}.")

    # k is DERIVED from the Tukey fence unless the caller states one. A caller-supplied k is
    # allowed - a re-run reproducing an earlier cohort needs it - but it is recorded as declared,
    # because a number that came from an argument did not come from the data.
    if k is None:
        sel = select_mad_k({s: st for s, st in stats.items()
                            if all(key in st for key in ("median", "q1", "q3", "mad"))},
                           mult=mult)
        k, k_notes, k_source = sel["k"], sel["notes"], "derived"
    else:
        sel = None
        k_notes = [f"MAD k = {k} DECLARED by the caller, not derived from this cohort's Tukey "
                   f"fences. Reproducing an earlier run is the reason to do this; choosing a "
                   f"number because the result looked better is not."]
        k_source = "declared"
    if not (k > 0):
        raise ThresholdRefusal(f"MAD multiplier k must be positive; got {k}. A fence at or below "
                               f"the median is not an upper fence.")

    out = {}
    for s, st in stats.items():
        # "mad" is required, not defaulted. A summary produced before this module applied the MAD
        # fence has no `mad` key, and silently falling back to Tukey would apply a DIFFERENT
        # threshold under the same name - the exact substitution this pipeline exists to prevent.
        missing = [key for key in ("n", "median", "q1", "q3", "mad") if key not in st]
        if missing:
            raise ThresholdRefusal(
                f"{s}: per-library summary is missing {missing}. A ceiling cannot be derived from "
                f"a partial summary, and defaulting the missing part would invent the threshold. "
                f"If 'mad' is the missing key, the summary predates the MAD fence and the run must "
                f"be recomputed rather than reinterpreted.")
        q1, q3 = float(st["q1"]), float(st["q3"])
        if not (q1 <= q3):
            raise ThresholdRefusal(f"{s}: q1 ({q1}) exceeds q3 ({q3}) - the summary is not a "
                                   f"quartile summary of anything.")
        mad = float(st["mad"])
        if mad < 0:
            raise ThresholdRefusal(f"{s}: MAD is {mad}, which is not a deviation.")
        med = float(st["median"])
        raw = med + k * MAD_SCALE * mad          # APPLIED
        tuk = q3 + mult * (q3 - q1)              # cross-check, never applied
        c = min(max(raw, lo), hi)
        out[s] = MitoCeiling(
            sample=s, n=int(st["n"]), median=med, q1=q1, q3=q3, mad=mad,
            derived=raw, tukey=tuk, ceiling=c,
            clamped="upper" if c < raw else ("lower" if c > raw else ""))

    lower = sorted(s for s, m in out.items() if m.clamped == "lower")
    upper = sorted(s for s, m in out.items() if m.clamped == "upper")
    n_clamped = len(lower) + len(upper)
    notes = [f"bound {lo}-{hi}% ({'assay default: ' + assay if bounds is None else 'DECLARED'})"]
    notes.extend(k_notes)
    if declared_by:
        notes.append(f"declared by the analyst: {declared_by!r}")

    # The two directions are reported apart because they are opposite events. An UPPER clamp means
    # the library's own fence asked for a ceiling higher than the analyst allows, so the bound
    # REMOVES nuclei the derivation would have kept - that is the direction that can delete signal.
    # A LOWER clamp means the fence was stricter than the analyst allows, so the bound RETAINS
    # nuclei the derivation would have cut. Counting them together gives one number that cannot
    # answer the only question worth asking of it: which way did the declaration push?
    notes.append(
        f"bound binds in {n_clamped} of {len(out)} libraries "
        f"(lower {len(lower)}: {', '.join(lower) if lower else 'none'}; "
        f"upper {len(upper)}: {', '.join(upper) if upper else 'none'})")

    # Reclassification, not refusal. Until 2026-08-11 this raised, on the reasoning that a bound
    # binding in most libraries IS the threshold and reporting it as derived would be false. The
    # reasoning is right and the remedy was wrong: refusing does not stop the number being wrongly
    # classified, it stops the run - and it fires hardest exactly when an analyst has deliberately
    # narrowed the bound, which is a legitimate thing to do. So the class changes instead, and it
    # travels with the result. A parameter whose class is honest can be argued with; a run that
    # will not start cannot.
    bound_dominated = bool(out) and n_clamped / len(out) > BOUND_BINDS_REVIEW
    provenance = "bound_dominated" if bound_dominated else "derived"
    if bound_dominated:
        notes.append(
            f"REVIEW - BOUND-DOMINATED, not derived. The bound {lo}-{hi}% decides the ceiling in "
            f"{n_clamped} of {len(out)} libraries, so for most of this cohort the applied number "
            f"is the DECLARATION and not that library's own fence. It is reported as such rather "
            f"than described as derived. This is usable and it is not a measurement: to make it a "
            f"derivation again, widen the bound; to keep it, say in the record that the ceiling "
            f"here is a declared constant with a per-library exception where the fence is tighter")
    else:
        notes.append("DERIVED per library, but the BOUND is declared - and what a mitochondria-high "
                     "population IS remains adjudicated; see mito_ceiling_note()")

    # Spread of the APPLIED ceiling. This is not the same property as the design differential and
    # is not implied by it: the differential asks whether removal falls evenly across the arms of
    # the design, this asks whether one filter is treating samples of one experiment alike at all.
    # A cohort can pass the differential while its ceiling varies fourfold, which is how a 6.12%
    # library and a 25.00% library sat in one deliverable without either check objecting.
    ceilings = [m.ceiling for m in out.values()]
    spread = (max(ceilings) / min(ceilings)) if ceilings and min(ceilings) > 0 else None
    if spread is not None:
        notes.append(f"applied ceiling spans {min(ceilings):.2f}-{max(ceilings):.2f}% "
                     f"({spread:.2f}x across {len(ceilings)} libraries)")
        if spread > CEILING_SPREAD_REVIEW:
            notes.append(
                f"REVIEW - the applied ceiling differs {spread:.2f}x between the least and most "
                f"permissive library. One filter is treating samples of one experiment very "
                f"differently; that is technical unless the libraries themselves differ that "
                f"much. Check whether the widest and narrowest sit in different arms of the "
                f"design before accepting it")

    # A MAD of zero means more than half this library's nuclei share one mitochondrial value, so
    # the fence collapses onto the median and every nucleus above it is called an outlier. The
    # bound catches the case where the median is below `lo`; above it, nothing else would.
    flat = sorted(s for s, m in out.items() if m.mad <= 0)
    if flat:
        notes.append(
            f"REVIEW - MAD is zero in {', '.join(flat)}. More than half of those libraries' nuclei "
            f"carry one mitochondrial value, so the fence sits ON the median and roughly half the "
            f"library is above it. That is a property of the measurement, not an outlier "
            f"population; do not apply the ceiling there without looking at the distribution")

    # The independent second derivation. Only UNCLAMPED libraries are compared: two clamped values
    # agree because the bound made them agree, and counting that as corroboration reports
    # confidence the data did not supply.
    checked = {s: m for s, m in out.items() if not m.clamped and m.cross_check is not None}
    diverged = sorted(s for s, m in checked.items()
                      if m.cross_check > CROSS_CHECK_REVIEW or m.cross_check < 1 / CROSS_CHECK_REVIEW)
    if checked:
        ratios = [m.cross_check for m in checked.values()]
        notes.append(
            f"cross-check against Tukey (Q3 + {mult}*IQR), on the {len(checked)} of {len(out)} "
            f"libraries the bound did not clamp: ratio {min(ratios):.2f}-{max(ratios):.2f}")
    else:
        notes.append(
            f"cross-check against Tukey NOT EVALUATED - the bound clamped every library, so the "
            f"two routes cannot be compared on any of them. This is not agreement")
    if diverged:
        notes.append(
            f"REVIEW - the two independent fences disagree by more than {CROSS_CHECK_REVIEW}x in "
            f"{', '.join(diverged)}. They are anchored differently (median vs Q3) and respond "
            f"differently to skew, so a large divergence says that library's tail has a shape "
            f"unlike its cohort's, not that one estimator is wrong. Look at it before accepting "
            f"the ceiling")

    return {"ceilings": out, "bounds": (lo, hi), "mult": mult, "assay": assay, "k": k,
            "k_source": k_source, "k_selection": sel,
            "declared_by": declared_by, "notes": notes, "provenance": provenance,
            "clamped_lower": tuple(lower), "clamped_upper": tuple(upper),
            "ceiling_spread": spread, "cross_check_diverged": tuple(diverged),
            "cross_check_evaluated": tuple(sorted(checked)), "mad_zero": tuple(flat)}


# The design-differential refusal line and the materiality floor below which a ratio is not a
# meaningful statistic. Same values and same reasoning as modules/04_doublets/doublet_health.py
# and modules/02_cells/cellcall_gate.py - restated rather than imported so each module's refusal
# is legible where it fires.
DESIGN_REFUSE = 3.0
MATERIAL = 0.01     # worst arm must lose at least this fraction before the ratio binds
# A threshold written in decimal is not the threshold the machine tests. This is not theoretical
# here: on the calibration cohort a candidate fixed ceiling produced a diet ratio of EXACTLY
# 3.00x, which is precisely the value a bare `>= 3.0` can miss by one unit in the last place.
ON_THE_LINE = 1e-9


def _at_least(value: float, threshold: float) -> bool:
    """`value >= threshold`, counting a value that is exactly on the line as on it."""
    return value >= threshold - abs(threshold) * ON_THE_LINE


def assess_mito_removal(per_library, ceilings, design=None) -> dict:
    """What the derived ceilings WOULD remove, and whether that is even across the design.

    `per_library` maps sample -> per-nucleus mitochondrial percentages, the same input
    `derive_mito_ceiling` saw. `ceilings` is its `{"ceilings": ...}` output or the bare mapping.
    `design` maps factor -> {sample: level}; without it the check that matters most cannot run,
    and that is reported rather than passed over.

    Refuses when one design arm loses at least DESIGN_REFUSE times the fraction another does.
    That is rule one Q3: a filter removing 53% of one arm and 6% of another has converted a
    technical property into an apparent biological difference, and nothing downstream can undo it.
    """
    cs = ceilings.get("ceilings", ceilings) if isinstance(ceilings, dict) else ceilings
    per_sample_rate, removed, total = {}, 0, 0
    for s, vals in per_library.items():
        if s not in cs:
            raise ThresholdRefusal(
                f"no derived ceiling for {s!r}. Refusing rather than filtering one library on "
                f"another library's number.")
        c = cs[s].ceiling if hasattr(cs[s], "ceiling") else float(cs[s])
        v = [float(x) for x in vals if x == x]
        n = sum(1 for x in v if x > c)
        per_sample_rate[s] = (n / len(v)) if v else 0.0
        removed += n
        total += len(v)

    out = {"per_sample_rate": per_sample_rate, "n_removed": removed, "n_total": total,
           "overall_rate": (removed / total) if total else 0.0, "arms": {}, "notes": []}

    if not design:
        out["notes"].append(
            "NO DESIGN SUPPLIED - the differential check did not run. This is the check that "
            "catches a filter converting a technical property into an apparent biological "
            "difference; its absence is not a pass.")
        return out

    worst = None
    for factor, mapping in design.items():
        levels = {}
        for s, r in per_sample_rate.items():
            lv = mapping.get(s)
            if lv is None:
                continue
            n_s = len([x for x in per_library[s] if x == x])
            levels.setdefault(str(lv), [0, 0])
            levels[str(lv)][0] += r * n_s
            levels[str(lv)][1] += n_s
        rates = {k: (a / b if b else 0.0) for k, (a, b) in levels.items()}
        if len(rates) < 2:
            continue
        lo_r, hi_r = min(rates.values()), max(rates.values())
        ratio = (hi_r / lo_r) if lo_r > 0 else None
        out["arms"][factor] = {"levels": rates, "ratio": ratio}
        if ratio is not None and (worst is None or ratio > worst[1]):
            worst = (factor, ratio, hi_r)

    if worst is not None:
        factor, ratio, hi_r = worst
        out["worst_factor"], out["worst_ratio"] = factor, ratio
        # The materiality floor exists because a ratio computed on tiny losses is arithmetic, not
        # evidence: 0.02% against 0.005% is 4x and means nothing. It cuts the other way too - at a
        # very LARGE removal the 3x test cannot fail, because three times the smallest arm exceeds
        # 100%. Both are stated so a passing check is not read as a strong one.
        if hi_r < MATERIAL:
            out["notes"].append(
                f"worst arm ({factor}) loses {100*hi_r:.3f}%, below the {100*MATERIAL:.0f}% "
                f"materiality floor - the {ratio:.2f}x ratio is arithmetic on small numbers and "
                f"does not bind.")
        elif _at_least(ratio, DESIGN_REFUSE):
            raise ThresholdRefusal(
                f"the mitochondrial ceiling removes unevenly across {factor}: "
                + ", ".join(f"{k} {100*v:.2f}%" for k, v in out['arms'][factor]['levels'].items())
                + f" - a {ratio:.2f}x differential against a {DESIGN_REFUSE}x refusal line. A "
                f"filter this uneven across a design factor has converted a technical property "
                f"into an apparent biological difference, and no downstream analysis can undo "
                f"that. Do not apply it.")
        else:
            out["notes"].append(
                f"worst design-arm ratio {ratio:.2f}x on {factor}, under the "
                f"{DESIGN_REFUSE}x refusal line.")
        if 3 * min(min(a["levels"].values()) for a in out["arms"].values()) > 1.0:
            out["notes"].append(
                "NOTE: at this removal rate the 3x test CANNOT fail - three times the lowest arm "
                "exceeds 100%. The observed ratios are informative; the test passing is not.")
    return out
