# Audits a CellBender run for abnormal removal and returns findings.
# It removes nothing, writes no object and subsets no data; every mask selects rows to count.
# Its purpose is to refuse a denoising run, never to alter one.
"""Step 1c - audit a CellBender run for abnormal removal.

WHAT AN ABNORMAL REMOVAL LOOKS LIKE, AND WHY THE TOTAL WILL NOT SHOW IT

CellBender reports a fraction removed. That number is almost useless on its own: the calibration
cohort's ten libraries span 15.5-23.7% and all ten are fine, while the one library that WAS
degenerate sat at 8.3% - inside a plausible-looking range, low rather than high. A denoiser that
has failed does not announce it; it returns a matrix, a report and a believable cell count.

So the auditor asks five questions the total cannot answer, in descending order of how much
damage a missed answer does.

  1 IS THE REMOVAL DIFFERENTIAL ACROSS THE DESIGN?
     The one that matters most and the only one that can silently become a result. If ambient
     removal takes 30% from one arm and 10% from another, a technical property has been
     converted into an apparent biological difference and nothing downstream can undo it. An
     even run measures close to 1x on every factor; the calibration cohort gives condition
     1.08x, stratum 1.07x, chemistry 1.07x and batch 1.04x. The refusal line is the
     design-differential threshold, bounded twice so that it cannot fire on a clean run:

       ratio >= threshold, worst arm material     REFUSE
       ratio >= threshold, worst arm immaterial   REVIEW - two near-zero rates, reported
       one arm at zero, another material          REFUSE - undefined is the MOST one-sided
       one arm at zero, another immaterial        REVIEW - undefined, but nothing was lost
       a sample carries no level for the factor   REVIEW - it is named, never dropped in silence
       fewer than two levels present              REVIEW - NOT CHECKED is its own outcome

  2 IS ANY SAMPLE UNLIKE ITS SIBLINGS?
     A degenerate fit is detectable only against a cohort. Cohort-relative, MAD-based.

  3 IS ANY GENE BEING GUTTED?
     Ambient removal should shave a gene, not delete it. In the calibration cohort the worst
     median per-gene removal is 34.6% and no gene exceeds 90%. A gene losing nearly all its
     counts is either genuinely ambient or a marker being destroyed, and the two are
     indistinguishable from the fraction alone - so the auditor reports the LIST, never a
     count. The removal checklist, question 1 (docs/PRINCIPLES.md): print the symbols and read
     them.

  4 IS DETECTION COLLAPSING?
     A gene detected in 40% of droplets before and 2% after has not been shaved, it has been
     removed from the analysis in all but name. Counts and detection can move independently.

  5 HOW MANY GENES VANISH ENTIRELY?
     In the calibration cohort, 12-240 per library. An order of magnitude more is a different
     kind of run.

WHAT THIS AUDITOR DOES NOT DO. It does not say the denoising is CORRECT. Every check is either
cohort-relative or a differential; none can tell a well-removed ambient transcript from a
badly-removed real one, because that requires knowing which cells should have expressed it -
which is annotation, and is not available here. It refuses runs that are unlike their siblings
or uneven across the design. That is a weaker claim than correct and is stated as one.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import median

MAD_K = 3.5
DESIGN_REFUSE = 3.0 # the design-differential check: worst arm this many times the best
# A ratio between two negligible rates is arithmetic, not evidence: arms removing 0.01% and
# 0.05% differ by 5x and neither has meaningfully touched the data. So the ratio binds only
# where the worst arm is also material; below the floor it is REPORTED and not refused. A gate
# that fires on correct behaviour gets switched off.
MATERIAL_REMOVAL = 0.01 # worst arm must remove at least this fraction before the ratio binds


# 0.30 / 0.10 is 2.9999999999999996 in binary floating point, so a bare `>= 3.0` misses the
# exact case the threshold was written for. Compare with a relative tolerance instead.
RATIO_TOL = 1e-9


def _material(rate: float) -> bool:
    """Is the worst arm removing enough for a ratio between arms to mean anything?

    Compared with the same relative tolerance the ratio uses. A bare `>=` on a float
    at an exact threshold is the defect fixed one branch above: a mean of two arms at
    0.26% and 1.74% is 0.009999999999999998, nominally 1.00%, and read as immaterial
    while the same cohort written as 1.00%/1.00% is read as material.
    """
    return rate >= MATERIAL_REMOVAL * (1 - RATIO_TOL)

GENE_GUT = 0.90 # a gene losing >=90% of its counts
DETECTION_COLLAPSE = 0.80 # detected fraction falling by >=80% of its original
TOP_LIST = 25 # how many symbols to print; the list is the finding, not the count

@dataclass
class Finding:
    check: str
    severity: str # REFUSE | REVIEW | ok
    message: str
    detail: list = field(default_factory=list)

    def __str__(self) -> str:
        tag = {"REFUSE": "REFUSE", "REVIEW": "REVIEW", "ok": "ok "}[self.severity]
        s = f"[{tag}] {self.check}\n {self.message}"
        for d in self.detail[:TOP_LIST]:
            s += f"\n - {d}"
        if len(self.detail) > TOP_LIST:
            s += f"\n ... and {len(self.detail) - TOP_LIST} more"
        return s

def _mad_out(values: dict, k: float = MAD_K) -> dict:
    vals = [v for v in values.values() if v is not None]
    if len(vals) < 4:
        return {}
    med = median(vals)
    mad = median([abs(v - med) for v in vals])
    if mad <= 0:
        # More than half the values are identical, so the robust scale is zero and `or 1e-12`
        # turned ANY deviation into ~1e12 robust SD - a 0.07% difference read as a degenerate
        # fit. With no spread to measure against there is no outlier to report.
        return {}
    return {s: (v, med, abs(v - med) / (1.4826 * mad))
            for s, v in values.items()
            if v is not None and abs(v - med) / (1.4826 * mad) > k}

def _rows(table):
    """A list of dict-like rows from either a list of dicts or a DataFrame.

    This module has no pandas dependency and must not acquire one: the process that calls it is
    the orchestrator, which on a cluster is a bare interpreter, because the aligner, the denoiser
    and the analysis stack have incompatible pins and live elsewhere. Requiring a DataFrame here
    forced `import pandas` into the engine and step 1's audit died on ModuleNotFoundError.

    Both shapes are accepted so existing callers keep working.
    """
    if table is None:
        return []
    if hasattr(table, "iterrows"):
        return [r for _, r in table.iterrows()]
    return list(table)


def _median_by(rows, key, value):
    """{key: median(value)} - the one pandas groupby this module used, done in the standard lib."""
    buckets: dict = {}
    for r in rows:
        v = r.get(value)
        if v is None:
            continue
        try:
            buckets.setdefault(r[key], []).append(float(v))
        except (TypeError, ValueError):
            continue
    return {k: median(v) for k, v in buckets.items() if v}


def audit(summary, per_gene=None, design=None) -> list:
    """summary: rows (list of dicts, or a DataFrame) with sample, fraction_removed_overall,
                genes_fully_removed.
    per_gene: optional rows with sample, symbol, fraction_removed, raw_detection_frac,
              denoised_detection_frac.
    design: optional {factor: {sample: level}} for the differential check.
    """
    out = []
    summary = _rows(summary)
    per_gene = _rows(per_gene)

    # ---- 1. differential across the design. First, because it is the one that becomes a result.
    if design:
        for factor, mapping in design.items():
            by, unmapped = {}, []
            for r in summary:
                lvl = mapping.get(r["sample"])
                if lvl is None:
                    unmapped.append(str(r["sample"]))
                    continue
                by.setdefault(lvl, []).append(r["fraction_removed_overall"])
            # A sample with no level is not evidence of anything - but dropping it in silence
            # makes a differential computed over part of the cohort read as one computed over
            # all of it. Name them.
            if unmapped:
                out.append(Finding(
                    f"design coverage: {factor} (unmapped samples)", "REVIEW",
                    f"{len(unmapped)} of {len(summary)} sample(s) carry no level for "
                    f"'{factor}', so the differential below speaks for part of the cohort "
                    f"only - READ THE LIST",
                    sorted(unmapped)))
            if len(by) < 2:
                out.append(Finding(
                    f"design differential: {factor}", "REVIEW",
                    f"NOT CHECKED - {len(by)} level(s) present among the mapped samples and a "
                    f"differential needs two. Never-checked is its own outcome and must not "
                    f"read as a pass"))
                continue
            means = {k: sum(v) / len(v) for k, v in by.items()}
            lvls = " · ".join(f"{k} {100*v:.2f}%" for k, v in sorted(means.items()))
            hi, lo = max(means.values()), min(means.values())
            if hi <= 0:
                out.append(Finding(
                    f"design differential: {factor}", "ok",
                    f"no ambient removed at any level ({lvls}) - nothing to be differential"))
            elif lo <= 0:
                # Where the denominator is zero the ratio is UNDEFINED, not infinite. Dividing
                # by a floor of 1e-12 manufactures a 1e8x differential out of arms at 0.00%
                # and 0.01% and refuses a run that removed almost nothing.
                #
                # Materiality still decides the severity, and it must: an undefined ratio is the
                # MOST one-sided a removal can be, not an exempt case. Routing every
                # zero-denominator case to REVIEW made a wholly one-sided material removal
                # milder than a merely large one - 50% against 0.0000001% refused, and the same
                # 50% against exactly 0% did not.
                if _material(hi):
                    out.append(Finding(
                        f"design differential: {factor}", "REFUSE",
                        f"one level removes nothing while another removes {100*hi:.2f}% of its "
                        f"counts ({lvls}). REFUSED: the design-differential check - the ratio is "
                        f"UNDEFINED rather than large, which is the most one-sided a removal can "
                        f"be, and it is material. A technical property has become an apparent "
                        f"biological difference"))
                else:
                    out.append(Finding(
                        f"design differential: {factor}", "REVIEW",
                        f"one level removes nothing while another does ({lvls}). The ratio is "
                        f"UNDEFINED rather than large, and the worst arm is below the "
                        f"{100*MATERIAL_REMOVAL:.0f}% materiality floor - reported, not refused. "
                        f"A one-sided removal across a design factor is still exactly what this "
                        f"check exists to catch, however small"))
            else:
                ratio = hi / lo
                exceeds = ratio >= DESIGN_REFUSE * (1 - RATIO_TOL)
                material = _material(hi)
                if exceeds and material:
                    sev, note = "REFUSE", (
                        f" REFUSED: the design-differential check - and the removal is "
                        f"material, the worst arm losing {100*hi:.2f}% of its counts. A "
                        f"technical property has become an apparent biological difference")
                elif exceeds:
                    sev, note = "REVIEW", (
                        f" ratio exceeds {DESIGN_REFUSE:.0f}x but the removal is NOT material "
                        f"- the worst arm removes {100*hi:.2f}%, under the "
                        f"{100*MATERIAL_REMOVAL:.0f}% floor. A ratio between two near-zero "
                        f"rates is dominated by single libraries; reported, not refused")
                else:
                    sev, note = "ok", ""
                out.append(Finding(
                    f"design differential: {factor}", sev,
                    f"max/min = {ratio:.3f}x ({lvls}){note}"))
    else:
        out.append(Finding("design differential", "REVIEW",
                           "no design given - the check that matters most was not run"))

    # ---- 2. cohort outliers
    fr = {r["sample"]: r["fraction_removed_overall"] for r in summary}
    o = _mad_out(fr)
    out.append(Finding(
        "cohort outlier: fraction removed",
        "REVIEW" if o else "ok",
        (f"{len(o)} sample(s) unlike the cohort" if o else
         f"all {len(fr)} within {MAD_K} robust SD "
         f"(range {100*min(fr.values()):.1f}-{100*max(fr.values()):.1f}%)"),
        [f"{s}: {100*v:.1f}% against a cohort median of {100*m:.1f}% ({z:.1f} robust SD)"
         for s, (v, m, z) in o.items()]))

    # ---- 5. genes vanishing entirely
    if "genes_fully_removed" in summary:
        gf = {r["sample"]: float(r["genes_fully_removed"]) for r in summary}
        o = _mad_out(gf)
        out.append(Finding(
            "genes removed entirely", "REVIEW" if o else "ok",
            (f"{len(o)} sample(s) unlike the cohort" if o else
             f"{int(min(gf.values()))}-{int(max(gf.values()))} per library, no outlier"),
            [f"{s}: {int(v):,} genes against a cohort median of {int(m):,} ({z:.1f} robust SD)"
             for s, (v, m, z) in o.items()]))

    if per_gene is None or not len(per_gene):
        out.append(Finding("per-gene checks", "REVIEW",
                           "no per-gene table given - gutting and detection collapse not checked"))
        return out

    # ---- 3. genes being gutted. The LIST is the finding.
    g = _median_by(per_gene, "symbol", "fraction_removed")
    gutted = sorted(k for k, v in g.items() if v >= GENE_GUT)
    out.append(Finding(
        "genes gutted", "REVIEW" if gutted else "ok",
        (f"{len(gutted)} gene(s) lose >={100*GENE_GUT:.0f}% of their counts - READ THE LIST. "
         f"Ambient and a destroyed marker are indistinguishable from the fraction alone"
         if gutted else
         (f"no gene loses >={100*GENE_GUT:.0f}%; worst median removal is "
          f"{100*max(g.values()):.1f}%" if g else "no per-gene rows to check")),
        gutted))

    # ---- 4. detection collapse
    cols = set(per_gene[0]) if per_gene else set()
    if {"raw_detection_frac", "denoised_detection_frac"} <= cols:
        raw = _median_by(per_gene, "symbol", "raw_detection_frac")
        den = _median_by(per_gene, "symbol", "denoised_detection_frac")
        # Genes barely detected to begin with are ignored: a collapse from 0.4% to 0.1% is noise
        # wearing the shape of a finding.
        drop = {k: 1 - (den[k] / raw[k]) for k in raw
                if raw[k] > 0.01 and k in den}
        coll = dict(sorted(((k, v) for k, v in drop.items() if v >= DETECTION_COLLAPSE),
                           key=lambda kv: -kv[1]))
        out.append(Finding(
            "detection collapse", "REVIEW" if len(coll) else "ok",
            (f"{len(coll)} gene(s) lose >={100*DETECTION_COLLAPSE:.0f}% of the droplets they "
             f"were detected in - shaved in counts is not the same as removed from the analysis"
             if len(coll) else
             f"no gene loses >={100*DETECTION_COLLAPSE:.0f}% of its detection; worst is "
             f"{100*max(drop.values()):.1f}%" if drop else "no gene passed the detection floor"),
            [f"{s}: detected {100*raw[s]:.1f}% -> {100*den[s]:.1f}%" for s in coll]))
    return out

def verdict(findings) -> str:
    if any(f.severity == "REFUSE" for f in findings):
        return "REFUSE"
    return "REVIEW" if any(f.severity == "REVIEW" for f in findings) else "PASS"
