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
    """What the mitochondrial ceiling is, and what remains a judgement after it is derived."""
    return (
        "mitochondrial ceiling: DERIVED PER LIBRARY, bounded by a DECLARED statement. The "
        "distribution is unimodal, so there is no valley and the count-floor route does not "
        "apply - but 'no valley' does not mean 'no derivation'. Each library's own upper Tukey "
        "fence (Q3 + 1.5*IQR of its own distribution) is derived by mito_ceiling(); the BOUND on "
        "that fence is declared by the analyst, because it is a statement about what a nucleus "
        "can be and not a property of this cohort. What stays ADJUDICATED is whether a "
        "mitochondria-high POPULATION is damage or a mitochondria-rich cell type - that needs an "
        "identity, so the pipeline emits cluster-level medians and stops.")


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
# it binds in most libraries it has stopped being a policy and become the threshold, and
# `derive_mito_ceiling` refuses in that case rather than letting a declared number masquerade as
# a derived one.

# Assay defaults. NOT thresholds - the outer limits within which a derived fence is allowed to
# land. Whole cells retain cytoplasm and legitimately carry more mitochondrial signal than nuclei.
MITO_BOUNDS = {"snrna": (5.0, 25.0), "scrna": (10.0, 30.0)}
IQR_MULT = 1.5
# If the bound binds in more than this fraction of libraries it IS the threshold, not a rail.
BOUND_BINDS_REVIEW = 0.5


@dataclass
class MitoCeiling:
    """One library's derived ceiling and the quantities it came from."""
    sample: str
    n: int
    median: float
    q1: float
    q3: float
    derived: float          # the raw fence, before the bound
    ceiling: float          # what would be applied
    clamped: str            # "", "lower" or "upper"

    @property
    def iqr(self) -> float:
        return self.q3 - self.q1


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


def derive_mito_ceiling(per_library, assay="snrna", bounds=None, mult=IQR_MULT,
                        declared_by=None) -> dict:
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
        stats[s] = {"n": len(v), "median": v[len(v) // 2], "q1": q1, "q3": q3}
    return derive_mito_ceiling_from_quartiles(
        stats, assay=assay, bounds=bounds, mult=mult, declared_by=declared_by)


def derive_mito_ceiling_from_quartiles(stats, assay="snrna", bounds=None, mult=IQR_MULT,
                                       declared_by=None) -> dict:
    """As `derive_mito_ceiling`, from precomputed per-library quartiles.

    `stats` maps sample -> {"n", "median", "q1", "q3"}. This is what the pipeline uses: the
    worker that already has the matrix open computes four numbers per library instead of
    returning a per-nucleus array, so nothing scales with cell count across the task boundary.
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

    out = {}
    for s, st in stats.items():
        missing = [k for k in ("n", "median", "q1", "q3") if k not in st]
        if missing:
            raise ThresholdRefusal(
                f"{s}: quartile summary is missing {missing}. A ceiling cannot be derived from a "
                f"partial summary, and defaulting the missing part would invent the threshold.")
        q1, q3 = float(st["q1"]), float(st["q3"])
        if not (q1 <= q3):
            raise ThresholdRefusal(f"{s}: q1 ({q1}) exceeds q3 ({q3}) - the summary is not a "
                                   f"quartile summary of anything.")
        raw = q3 + mult * (q3 - q1)
        c = min(max(raw, lo), hi)
        out[s] = MitoCeiling(
            sample=s, n=int(st["n"]), median=float(st["median"]), q1=q1, q3=q3,
            derived=raw, ceiling=c,
            clamped="upper" if c < raw else ("lower" if c > raw else ""))

    n_clamped = sum(1 for m in out.values() if m.clamped)
    notes = [f"bound {lo}-{hi}% ({'assay default: ' + assay if bounds is None else 'DECLARED'})"]
    if declared_by:
        notes.append(f"declared by the analyst: {declared_by!r}")
    if out and n_clamped / len(out) > BOUND_BINDS_REVIEW:
        raise ThresholdRefusal(
            f"the declared bound {lo}-{hi}% binds in {n_clamped} of {len(out)} libraries. A bound "
            f"is a guard rail on a derived fence; when it binds in most libraries it IS the "
            f"threshold, and reporting it as derived would be false. Either widen the bound, or "
            f"apply it deliberately as a declared constant and say so.")
    notes.append(f"bound binds in {n_clamped} of {len(out)} libraries")
    notes.append("DERIVED per library, but the BOUND is declared - and what a mitochondria-high "
                 "population IS remains adjudicated; see mito_ceiling_note()")
    return {"ceilings": out, "bounds": (lo, hi), "mult": mult, "assay": assay,
            "declared_by": declared_by, "notes": notes}


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
