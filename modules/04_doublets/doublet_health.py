# This module health-checks doublet calls: it reads the calls and returns findings.
# It removes nothing and alters no call.
"""Step 4b - the failure signatures, applied to whatever detector ran.

WHY THIS EXISTS INSTEAD OF A SECOND DETECTOR

Scrublet failed on this cohort. The failure was not subtle: its automatic threshold returned ~0%
in seven of ten libraries, and the silence was CONDITION-ASSOCIATED - a usable rate in 3 of the 5
ctrl libraries and in 0 of the 5 treated ones, which lands a technical collapse squarely on one
arm of the study's primary readout.

Nothing about a vote was required to see that. A look at the per-sample rates was enough. So the
protection against a detector failing is a check on its output, not a second detector whose own
failures are then averaged in - and this check applies to a supplied detector exactly as it
applies to the shipped one.

WHAT EACH CHECK CATCHES, AND WHAT IT CANNOT

  zero_rate      a library where the detector called (almost) nothing. Absence of a call is
                 not evidence of absence, and a threshold that collapses returns silence.
  design         the removal rate per level of every design factor, which is the check that can
                 silently become a result. Bound by materiality: a ratio between two near-zero
                 rates is dominated by single libraries, and 541 cells against 76 is arithmetic
                 rather than evidence.
  spread         per-sample variation, read in BOTH directions - too little means the rate is
                 imposed rather than measured, too much means the threshold is unstable.
  imposed        rate variance across libraries at or near zero. scds returned exactly 6.00%
                 ten times over because its call was a fixed quantile; a detector that returns
                 the same fraction of every library is not measuring any of them.
  deep_decile    the share of the deepest UMI decile called. COMPUTED, NOT INTERPRETED: where
                 the deepest nuclei are a large, RNA-rich majority population a high value is
                 alarming, and where no such population dominates it means something else. The
                 pipeline reports the number and refuses to decide which case it is in.
  coverage       nuclei that reached the deliverable without ever being scored. Their doublet
                 status is UNKNOWN, not negative.

REPORTING A REFUSAL IS NOT ENFORCING ONE

`verdict()` returns the findings as one severity, for printing and for tables. `enforce()` raises
`HealthRefusal` on a REFUSE verdict; use it where a refusal must halt the caller, because a
returned "REFUSE" is a string and a string can be ignored by accident. The bundled CLI instead
prints the findings and exits 2 - the same decision, made visible. REVIEW is returned rather than
raised on purpose: it means a human must look, not that the run is wrong.

None of these says the calls are CORRECT. There is no ground truth in this experiment - no cell
hashing, no SNP demultiplexing - and homotypic doublets are invisible to every method at every
setting, so the rate is a floor. These checks catch a detector that has failed, not one that is
wrong in a way its output does not show.
"""

from __future__ import annotations

from dataclasses import dataclass, field

ZERO_RATE = 0.005 # a library calling <0.5% has effectively returned silence
SPREAD_IMPOSED = 1.05 # rate variance this low means the fraction is being set, not measured
SPREAD_UNSTABLE = 5.0
DESIGN_REFUSE = 3.0
MATERIAL = 0.01 # worst arm must call at least this before the ratio binds
# A threshold written in decimal is not the threshold the machine tests. 0.30 / 0.10 is
# 2.9999999999999996 in binary floating point, so a differential that is EXACTLY 3x slips under a
# bare `>= 3.0` by one unit in the last place - and exactly-on-the-line is the case a reader
# checks by hand, so it is the one that must not silently pass. The ratio and spread comparisons
# go through the two helpers below, which count a value within this relative distance of a
# threshold as being on it.
ON_THE_LINE = 1e-9

def _at_least(value: float, threshold: float) -> bool:
    """`value >= threshold`, counting a value that is exactly on the line as on it."""
    return value >= threshold - abs(threshold) * ON_THE_LINE

def _at_most(value: float, threshold: float) -> bool:
    """`value <= threshold`, counting a value that is exactly on the line as on it."""
    return value <= threshold + abs(threshold) * ON_THE_LINE

class HealthRefusal(RuntimeError):
    """Raised when the calls carry a failure signature that must stop the run."""

@dataclass
class HealthFinding:
    check: str
    severity: str
    message: str
    detail: list = field(default_factory=list)

    def __str__(self) -> str:
        s = f"[{self.severity:6s}] {self.check}\n {self.message}"
        for d in self.detail:
            s += f"\n - {d}"
        return s

def health(per_sample_rate, design=None, deep_decile_rate=None, n_kept_unscored=None,
           detector_name="detector", reproducible=True) -> list:
    """per_sample_rate: {sample: fraction of that library called a doublet}.
    design: {factor: {sample: level}} - without it the check that matters most cannot run.
    deep_decile_rate: fraction of the deepest UMI decile called, if it was computed.
    n_kept_unscored: nuclei in the deliverable that were never scored.
    Returns findings. Use `verdict()` to summarise them and `enforce()` to act on them.
    """
    out = []
    rates = per_sample_rate
    v = list(rates.values())

    zero = {s: r for s, r in rates.items() if r < ZERO_RATE}
    out.append(HealthFinding(
        "silent libraries", "REFUSE" if zero else "ok",
        (f"{len(zero)} library(ies) called under {100*ZERO_RATE:.1f}% - the threshold has "
         f"collapsed and returned silence, which is not the same as finding no doublets. This "
         f"is Scrublet's failure signature: 0% in seven of ten"
         if zero else f"every library called above {100*ZERO_RATE:.1f}%"),
        [f"{s}: {100*r:.2f}%" for s, r in sorted(zero.items())]))

    if not v:
        sev, msg = "REVIEW", ("no per-sample rates were given, so the spread could not be "
                              "computed - not measured is not the same as unremarkable")
    elif min(v) <= 0:
        # A zero denominator makes max/min UNDEFINED. Guarding it with an epsilon would print
        # ~1e12, which reads as an extreme measurement rather than as a missing one - and the
        # library that produced the zero is already refused by the silent-library check above.
        sev, msg = "REVIEW", (
            "spread is UNDEFINED - at least one library called nothing, so max/min has no "
            "value. Reported as undefined rather than as a very large number: a zero "
            "denominator is a missing measurement, not an extreme one")
    else:
        spread = max(v) / min(v)
        if _at_most(spread, SPREAD_IMPOSED):
            sev, msg = "REVIEW", (
                f"spread {spread:.2f}x - the rate barely varies across libraries. A detector "
                f"returning the same fraction of every library is imposing it, not measuring "
                f"them. scds returned exactly 6.00% ten times because its call was a fixed "
                f"quantile")
        elif _at_least(spread, SPREAD_UNSTABLE):
            sev, msg = "REVIEW", (
                f"spread {spread:.2f}x - the rate varies more than {SPREAD_UNSTABLE:.0f}-fold "
                f"across libraries. Either they genuinely differ that much or the threshold is "
                f"unstable, and the output alone cannot separate those")
        else:
            sev, msg = "ok", f"spread {spread:.2f}x across libraries"
    out.append(HealthFinding("rate spread", sev, msg))

    if design:
        for factor, mapping in design.items():
            by = {}
            for s, r in rates.items():
                lvl = mapping.get(s)
                if lvl is not None:
                    by.setdefault(lvl, []).append(r)
            if len(by) < 2:
                continue
            m = {k: sum(x) / len(x) for k, x in by.items()}
            hi, lo = max(m.values()), min(m.values())
            txt = " · ".join(f"{k} {100*x:.2f}%" for k, x in sorted(m.items()))
            if lo <= 0:
                # One whole arm called nothing. The ratio is undefined rather than large, so
                # no threshold on it can fire - but every library in that arm is under the
                # zero-rate line, and the silent-library check above has already refused.
                out.append(HealthFinding(
                    f"design differential: {factor}", "REVIEW",
                    f"one level of this factor called nothing at all ({txt}), so max/min is "
                    f"UNDEFINED rather than large and no threshold on it can fire. A detector "
                    f"that goes silent on one arm of the design is the failure this check "
                    f"exists for; the refusal for it comes from the silent-library check"))
                continue
            ratio = hi / lo
            material = _at_least(hi, MATERIAL)
            if _at_least(ratio, DESIGN_REFUSE) and material:
                sev, note = "REFUSE", (
                    f" REFUSED: the differential is at or above {DESIGN_REFUSE:.0f}x and the "
                    f"rate is material - the worst arm calls {100*hi:.2f}%. A detection rate "
                    f"that tracks an arm of the design turns a technical property into an "
                    f"apparent biological difference, and nothing downstream can undo it")
            elif _at_least(ratio, DESIGN_REFUSE):
                sev, note = "REVIEW", (f" ratio is at or above {DESIGN_REFUSE:.0f}x but the "
                                       f"worst arm calls only {100*hi:.2f}%, under the "
                                       f"{100*MATERIAL:.0f}% floor - a ratio between two "
                                       f"near-zero rates is dominated by single libraries, so "
                                       f"this is reported, not refused")
            else:
                sev, note = "ok", ""
            out.append(HealthFinding(f"design differential: {factor}", sev,
                                     f"max/min = {ratio:.2f}x ({txt}){note}"))
    else:
        out.append(HealthFinding(
            "design differential", "REVIEW",
            "no design given - the check that matters most was not run. A per-sample rate "
            "cannot show that the calls track one arm of the design, and that is the failure "
            "that becomes a result"))

    if deep_decile_rate is not None:
        out.append(HealthFinding(
            "deepest UMI decile", "REVIEW",
            f"{100*deep_decile_rate:.1f}% of the deepest decile called. REPORTED, NOT "
            f"INTERPRETED: where the deepest nuclei are a large, RNA-rich majority population "
            f"this is as consistent with eating that population as with finding doublets, and "
            f"depth cannot tell them apart - within matched depth bins, called doublets and "
            f"singlets have indistinguishable genes per 1,000 UMI. What it means is the "
            f"reader's to decide"))

    if n_kept_unscored is not None:
        out.append(HealthFinding(
            "coverage of the deliverable", "REVIEW" if n_kept_unscored else "ok",
            (f"{n_kept_unscored:,} retained nuclei were never scored - their doublet status is "
             f"UNKNOWN, not negative"
             if n_kept_unscored else
             "every retained nucleus was scored; the call covers the whole deliverable")))

    if not reproducible:
        out.append(HealthFinding(
            "reproducibility", "REVIEW",
            f"{detector_name} declares itself not reproducible - a re-run will move these calls "
            f"and every downstream count inherits that"))
    return out

def verdict(findings) -> str:
    """REPORT the findings as one severity - "REFUSE", "REVIEW" or "PASS" - and raise nothing.

    This is the reporting form, for printing and for tables. It is not the enforcing form: a
    caller that acts on these calls uses enforce(), because a returned "REFUSE" can be ignored
    by accident.
    """
    if any(f.severity == "REFUSE" for f in findings):
        return "REFUSE"
    return "REVIEW" if any(f.severity == "REVIEW" for f in findings) else "PASS"

def enforce(findings, detector_name: str = "detector") -> str:
    """ENFORCE the verdict: raise HealthRefusal on REFUSE, otherwise return "REVIEW" or "PASS".

    A failed detector is not repairable downstream - its rate decides what is scored as a
    doublet at all - so this is where the run stops. REVIEW is returned rather than raised on
    purpose: it means a human must look, not that the run is wrong, so a caller that wants to
    stop on REVIEW as well checks the returned value.
    """
    v = verdict(findings)
    if v == "REFUSE":
        raise HealthRefusal(
            f"the doublet health check REFUSES the calls from {detector_name}:\n"
            + "\n".join(str(f) for f in findings if f.severity == "REFUSE"))
    return v
