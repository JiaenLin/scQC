# Doublet scoring: this module scores nuclei and writes flags, and removes nothing.
# Whether a flagged nucleus is removed is a separate, adjudicated decision taken at step 7.
"""Step 4 - doublet detection: one detector, tuned, with a contract for supplying another.

ONE TOOL, NOT A VOTE - AND THE REASON IS MEASURED

A consensus - two of three detectors must agree before a nucleus is called - is the obvious
design. Four detectors were measured on the calibration cohort's ten libraries, and the numbers
rule it out:

    scds               6.00% in ALL TEN libraries, spread 1.00x - its call is
                       `top_frac(hybrid, 0.06)`, a quantile fixed in advance
    scDblFinder        10.04-10.79%, spread 1.07x
    DoubletDetection   1.82-20.90%, spread 11.48x - "PhenoGraph does not support
                       random seeds", by its own documentation
    pairwise Jaccard   0.21-0.28

Two of the three voters carried no rate information, and the consensus inherited its 3.21x
per-sample spread from the one member that cannot be reproduced. A vote does not protect against
a bad detector; it averages one in. So this module runs scDblFinder alone.

THE PROTECTION IS A DIAGNOSTIC, NOT A SECOND OPINION

Scrublet failed on this cohort and the failure was plainly visible - 0% in seven of ten
libraries, and the silence fell on one arm of the design. No vote was needed to see that, only a
look at the per-sample rates. So one detector runs and `doublet_health.py` checks its calls for
exactly those signatures. A supplied detector gets the same check with no exemption.

WHAT IS DECLARED, AND WHY THERE IS NO DEFAULT FOR IT

`dbr`, the expected doublet rate. scDblFinder's default is the 10x Chromium loading formula,
~(cells/1000)*0.01, and its own vignette says "different protocols may create considerably more
doublets, and that this should be updated accordingly". Applied to a Singleron run it tracked
library size at r = 0.872 and, because the ctrl libraries were larger, imported a condition
differential of 19.65% vs 14.98% onto half the study's primary readout. So `dbr` comes from the
platform and the module refuses rather than falling back to a formula belonging to another one.

TUNING IS `dbr.sd`, AND A FLAT RATE IS SUSPECT RATHER THAN REASSURING

At the package default, scDblFinder returned 10.04-10.79% across libraries differing 2.5-fold in
size - a 0.75 percentage-point spread. That flatness was the prior's, not the data's: at
`dbr.sd = 1`, which the documentation describes as disabling any expectation of the number of
doublets, the same tool spans 12.81-31.69%.

This is the argument against a per-sample MAD5 mitochondrial rule - a method that equalises the
fraction removed cannot reveal a library that genuinely differs - applied to a doublet caller.
But fully freeing the prior is not the answer either: at `dbr.sd = 1` the tool calls 63.7% of the
deepest UMI decile, and where the deepest nuclei are a large, RNA-rich majority population, that
is as consistent with eating that population as with finding doublets.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

SWEEP = ("default", "dbr", "1") # dbr.sd settings; "dbr" means dbr.sd = dbr
DEEP_DECILE_ALARM = 0.50 # share of the deepest UMI decile called

class DoubletRefusal(RuntimeError):
    """Raised when the detector cannot be run safely as configured."""

@runtime_checkable
class Detector(Protocol):
    """What a detector must DECLARE before it may be used.

    Every field here can only be learned from the tool's own documentation, and one of them -
    reproducibility - is stated there plainly and read by almost nobody.
    """
    name: str
    reproducible: bool # scDblFinder yes; DoubletDetection no, by its own docs
    needs_empty_drops_removed: bool # scDblFinder: yes
    min_umi_floor: int | None # scDblFinder: 200, "to avoid errors"
    imports_rate_prior: bool # if True, dbr must be DECLARED

    def score(self, matrix, sample: str, seed: int, **kw):
        """Return (scores, calls) for the nuclei given. Must not subset the input."""

@dataclass
class SweepResult:
    setting: str
    per_sample_rate: dict # sample -> fraction called
    worst_arm_ratio: float = None
    worst_arm_factor: str = ""
    deep_decile_rate: float = None
    in_published_band: int = None
    n_samples_in_band: int = None

    @property
    def spread(self) -> float:
        """max/min across libraries, or None when the lowest library called nothing.

        A rate of zero makes the ratio UNDEFINED, not enormous. Guarding the denominator with
        an epsilon returns ~1e12, which prints like a measurement and sorts like one; an
        undefined ratio has to be reported as undefined, not manufactured.
        """
        v = list(self.per_sample_rate.values())
        if not v or min(v) <= 0:
            return None
        return max(v) / min(v)

    def __str__(self) -> str:
        v = list(self.per_sample_rate.values())
        sp = self.spread
        sp_txt = f"{sp:5.2f}x" if sp is not None else "UNDEFINED (a library called nothing)"
        s = (f"dbr.sd={self.setting:8s} cohort {100*sum(v)/len(v):6.2f}% "
             f"per-sample {100*min(v):5.2f}-{100*max(v):5.2f} spread {sp_txt}")
        if self.worst_arm_ratio is not None:
            s += f" worst arm {self.worst_arm_ratio:.2f}x ({self.worst_arm_factor})"
        if self.deep_decile_rate is not None:
            s += f" deepest decile {100*self.deep_decile_rate:.1f}%"
        if self.in_published_band is not None:
            s += f" in band {self.in_published_band}/{self.n_samples_in_band}"
        return s

@dataclass
class Recommendation:
    setting: str
    reason: str
    rejected: dict = field(default_factory=dict)

    def __str__(self) -> str:
        s = [f"recommended dbr.sd = {self.setting}", f" {self.reason}"]
        s += [f" rejected {k}: {v}" for k, v in self.rejected.items()]
        return "\n".join(s)

def check_detector(det, dbr=None, light_floor=None) -> list:
    """Refuse a configuration the detector's own declarations say is unsafe."""
    notes = []
    for f in ("name", "reproducible", "needs_empty_drops_removed", "min_umi_floor",
              "imports_rate_prior"):
        if not hasattr(det, f):
            raise DoubletRefusal(
                f"detector does not declare '{f}'. Every field in the contract can only be "
                f"learned from the tool's own documentation; a detector that will not state "
                f"them cannot be checked.")
    if det.imports_rate_prior and dbr is None:
        raise DoubletRefusal(
            f"{det.name} imports a rate prior and no dbr was DECLARED. There is no default: "
            f"scDblFinder's is the 10x loading formula, which on a non-10x cohort tracked library "
            f"size at r = 0.872 and imported a 19.65% vs 14.98% condition differential.")
    if det.min_umi_floor is not None:
        if light_floor is None:
            raise DoubletRefusal(f"{det.name} requires a floor of {det.min_umi_floor} UMI and "
                                 f"none was applied")
        if light_floor < det.min_umi_floor:
            notes.append(f"light floor {light_floor} is below the {det.min_umi_floor} "
                         f"{det.name} documents; it may error on very low-coverage nuclei")
    if not det.reproducible:
        notes.append(f"{det.name} declares itself NOT reproducible - a re-run will move the "
                     f"calls, and any downstream count inherits that")
    return notes

def recommend(results, published_band=None) -> Recommendation:
    """Pick a dbr.sd setting from the sweep. Flat is suspect; fully free may be implausible."""
    by = {r.setting: r for r in results}
    rejected = {}

    # An undefined spread is not a small one. A setting at which some library called nothing
    # is the signature `doublet_health` refuses on, so it cannot be recommended - and it must
    # be excluded before any comparison, because there is no number to compare.
    for name, r in by.items():
        if r.spread is None:
            rejected[name] = (
                "at least one library called nothing at this setting, so the per-sample spread "
                "is UNDEFINED rather than large. A silent library is a collapsed threshold, not "
                "a measurement of that library")

    base = by.get("default")
    if base and base.spread is not None and base.spread < 1.15:
        rejected["default"] = (
            f"rate is flat across libraries ({base.spread:.2f}x) despite them differing in size "
            f"- the prior is setting the threshold, not the data. This is the property that "
            f"rules out a per-sample MAD5 mitochondrial rule: equalising the fraction removed "
            f"cannot reveal a library that genuinely differs")

    free = by.get("1")
    if free and free.deep_decile_rate and free.deep_decile_rate >= DEEP_DECILE_ALARM:
        rejected["1"] = (
            f"calls {100*free.deep_decile_rate:.1f}% of the deepest UMI decile. Where the "
            f"deepest nuclei are a large, RNA-rich majority population, that is as consistent "
            f"with eating that population as with finding doublets, and depth cannot "
            f"distinguish the two")

    for cand in ("dbr", "1", "default"):
        r = by.get(cand)
        if r is None or cand in rejected:
            continue
        why = [f"rate moves with the data (spread {r.spread:.2f}x)"]
        if r.worst_arm_ratio is not None:
            why.append(f"worst design arm {r.worst_arm_ratio:.2f}x")
        if r.in_published_band is not None:
            why.append(f"{r.in_published_band}/{r.n_samples_in_band} libraries inside the "
                       f"published band")
        return Recommendation(cand, "; ".join(why), rejected)

    return Recommendation("default",
                          "every swept setting was rejected - this needs a human, not a fallback",
                          rejected)
