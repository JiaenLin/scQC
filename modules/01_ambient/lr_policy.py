# Decides whether to re-run CellBender at a halved learning rate.
# It removes no observation and discards no output; both rates are kept.
"""Step 1b - is the default learning rate acceptable, and what to do when it is not.

THE DETECTION IS COHORT-RELATIVE. ONE SAMPLE CANNOT TELL YOU.

The same diagnostics answer different questions depending on how many libraries are in hand:

    scope       what it can establish                     what it cannot
    ----------  ----------------------------------------  ------------------------------
    one sample  that the optimiser complained, and where  whether the fit is degenerate
    the cohort  that one fit is unlike its siblings       whether any fit is correct

Testing the learning rate on a single library is not merely weak, it misleads. CellBender's own
report can advise halving because of a period of motion in the wrong direction near some epoch;
the half-rate run can repeat the same complaint a few epochs earlier; the two denoised matrices
can agree at r = 0.9999; and the half rate can call thousands FEWER cells. Every one of those
observations is equally compatible with the default being right and with it being wrong.

Beside its siblings the same evidence is legible. One library removing 8.3% of its counts where
the other nine remove 11-18%, with a convergence indicator of 1.71 against 0.35-0.73, is a
degenerate fit; at half the rate it becomes ordinary. On its own that library produced a matrix,
a report and a plausible cell count - it looked exactly like a result.

Two lessons the pipeline has to encode, and they pull in opposite directions:

  1. THE TOOL'S OWN SUGGESTION IS NOT EVIDENCE. A non-default adopted because a tool suggested
     it is an unjustified choice wearing a recommendation. So: halve, RE-MEASURE, and adopt only
     if the diagnostic actually resolves. Moving a complaint from one epoch to another is not a
     resolution.

  2. A DEGENERATE FIT IS INVISIBLE IN ITS OWN RUN. There is no per-sample number that announces
     it, which is why every check in this module is cohort-relative and why the module refuses
     to answer at all below about four libraries.

WHY THE FIX IS APPLIED COHORT-WIDE AND NOT ONLY TO THE FAILING SAMPLE

Running most libraries at one learning rate and one at another makes the denoising a technical
property that varies across the design - exactly what the design-differential check exists to
catch. If the outlier's arm is also the arm that differs biologically, no downstream analysis
can separate them. So a rate that has to change changes for everyone, and the superseded outputs
are kept.

ON convergence_indicator. This module does not assert a direction of "better" from it, because
its definition could not be located in the CellBender source. It is used only as an OUTLIER
statistic: you do not need to know what a number means to notice that one sample's is three
times every other's.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import median

# Diagnostics that are cohort-comparable. Each maps to "is this sample unlike its siblings",
# never to "is this value good", which is a claim none of them can support alone.
DIAGNOSTICS = ("fraction_removed", "convergence_indicator")
MAD_K = 3.5 # conventional robust-outlier cut; with n=10 it flags the plainly odd
HALVE = 0.5

@dataclass
class LRVerdict:
    outliers: dict = field(default_factory=dict) # sample -> [reasons]
    action: str = "keep_default" # keep_default | rerun_half | escalate
    note: str = ""

    def __str__(self) -> str:
        s = [f"action: {self.action}"]
        if self.note:
            s.append(f" {self.note}")
        for smp, why in self.outliers.items():
            s.append(f" OUTLIER {smp}")
            s += [f" - {w}" for w in why]
        return "\n".join(s)

def _mad_outliers(values: dict) -> dict:
    """Samples whose value is more than MAD_K robust deviations from the cohort median."""
    vals = [v for v in values.values() if v is not None]
    if len(vals) < 4:
        return {}
    med = median(vals)
    mad = median([abs(v - med) for v in vals])
    if mad <= 0:
        # More than half the values are identical, so the robust scale is zero. Guarding that
        # with `or 1e-12` would turn ANY deviation into ~1e12 robust SD - a 0.07% difference
        # read as a degenerate fit. With no spread to measure against there is no outlier.
        return {}
    out = {}
    for s, v in values.items():
        if v is None:
            continue
        z = abs(v - med) / (1.4826 * mad)
        if z > MAD_K:
            out[s] = (v, med, z)
    return out

def assess_cohort(diag: dict, label: str = "default") -> LRVerdict:
    """`diag` is {sample: {diagnostic: value}} for one learning rate, all samples."""
    v = LRVerdict()
    if len(diag) < 4:
        v.note = (f"only {len(diag)} samples - a degenerate fit is detected against its "
                  f"siblings, so this check is not meaningful below about four")
        return v
    for d in DIAGNOSTICS:
        vals = {s: m.get(d) for s, m in diag.items()}
        for s, (val, med, z) in _mad_outliers(vals).items():
            v.outliers.setdefault(s, []).append(
                f"{d} = {val:.4g} against a cohort median of {med:.4g} ({z:.1f} robust SD)")
    if v.outliers:
        v.action = "rerun_half"
        v.note = (f"{len(v.outliers)} sample(s) unlike the cohort at learning rate '{label}'. "
                  f"Re-running ALL samples at half - a rate that changes changes for everyone, "
                  f"or the denoising becomes a technical difference across the design.")
    else:
        v.note = f"no sample is an outlier at learning rate '{label}'; the default stands"
    return v

def compare_after_halving(before: dict, after: dict, flagged: list) -> LRVerdict:
    """Did halving actually resolve the outliers? Adopt only if it did."""
    v = LRVerdict()
    still = assess_cohort(after, label="half").outliers
    unresolved = [s for s in flagged if s in still]
    resolved = [s for s in flagged if s not in still]

    if unresolved and not resolved:
        v.action = "keep_default"
        v.outliers = {s: still[s] for s in unresolved}
        v.note = ("halving resolved NOTHING - the same samples remain outliers. The behaviour "
                  "is a property of the data, not of the learning rate, and adopting a "
                  "non-default that fixes nothing trades a defensible default for an "
                  "unjustified choice. A complaint that merely moves from one epoch to another "
                  "is the case this branch exists for.")
    elif resolved and not unresolved:
        v.action = "adopt_half"
        v.note = (f"halving resolved every flagged sample ({', '.join(resolved)}). Adopt it "
                  f"COHORT-WIDE and retain the default-rate outputs; the change is a stated "
                  f"departure from the package default and must be recorded as one.")
    elif resolved and unresolved:
        v.action = "escalate"
        v.outliers = {s: still[s] for s in unresolved}
        v.note = (f"halving resolved {len(resolved)} sample(s) and not {len(unresolved)}. That "
                  f"is not a hyperparameter answer - a rate that fixes some libraries and not "
                  f"others is describing something about those libraries. Do not halve again "
                  f"until the unresolved ones are understood.")
    else:
        v.action = "keep_default"
        v.note = "nothing was flagged; there was nothing for halving to resolve"
    return v
