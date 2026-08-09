# Step 0 ingest: validates a samplesheet and decides how each sample will be ingested.
# It removes no observation; it REFUSES inputs rather than altering them.
"""Step 0 - ingest: decide, per sample, whether a supplied matrix can be used or must be rebuilt.

THE MODE IS DECIDED BY THE DATA, NOT BY A FLAG.

A supplied matrix that looks like the pipeline's input frequently is not one. A delivered count
matrix is often the aligner's `outs/filtered`, already through cell calling, and nothing about
the file says so - not its name, not its shape, not a header field. Discovering it late costs a
full reprocess from FASTQ, and before it is discovered it supports a confident and wrong
conclusion: that ambient-RNA correction is not possible with the files in hand.

So a supplied matrix is not accepted because it exists. It is accepted because it passes
`lib/verify_raw.py`, which checks two INDEPENDENT properties that both hide under the word "raw":

    P1  raw VALUES     unnormalised integers, no ceiling, no gene subsetting
    P2  raw DROPLETS   every barcode, including the empties

A matrix can pass P1 and fail P2 - that is the usual failure - and P2 is the one that matters
next: an ambient model learns the background profile FROM the empty droplets.

WHAT IS DECLARED AND WHAT IS NOT

`platform`, `species` and `reference` are DECLARED with no defaults, because the pipeline cannot
infer them and a wrong guess is silent. The doublet-rate model is platform-specific: applying one
platform's loading formula to another platform's run imported a 19.65% against 14.98% condition
differential into the calibration cohort's primary readout. Gene-class patterns are equally
specific - the regex that selects mitochondrial or ribosomal genes depends on the naming
convention of the annotation, and one written for a different reference silently matches nothing
rather than failing.

`chemistry` is NOT declared. It is detected per sample by the processor and RECORDED, because in
the calibration cohort it differed WITHIN a single cohort - one chemistry version for half the
libraries and an earlier one for the rest - and a run-level declaration would have been wrong for
half the samples without ever saying so. It is also confounded with one of that cohort's design
factors, which is what makes that factor's coefficient uninterpretable. A pipeline that hides
chemistry hides the limitation along with it.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from verify_raw import verify # noqa: E402

REQUIRED = ("sample", "platform", "species", "reference")
PROCESSOR = {"10x": "cellranger", "singleron": "celescope"}
FUTURE = {"bgi": "not implemented - declare it and the run will refuse rather than guess"}

@dataclass
class Plan:
    sample: str
    platform: str
    processor: str
    mode: str # "accept" | "run" | "blocked"
    reason: str = ""
    verdict: object = None

    def __str__(self) -> str:
        head = {"accept": "ACCEPT", "run": "RUN ", "blocked": "BLOCKED"}[self.mode]
        s = f"[{head}] {self.sample:12s} {self.platform:10s} -> {self.processor}"
        if self.reason:
            s += f"\n {self.reason}"
        if self.verdict is not None and not self.verdict.usable:
            for r in self.verdict.reasons:
                s += f"\n - {r}"
        return s

def validate_row(row: dict, registry: dict) -> list:
    """Every DECLARED field checked. Missing means fail, never a default."""
    errs = []
    for c in REQUIRED:
        if not str(row.get(c, "")).strip():
            errs.append(f"'{c}' is required and has no default - the pipeline cannot infer it")
    p = str(row.get("platform", "")).strip().lower()
    if p in FUTURE:
        errs.append(f"platform '{p}': {FUTURE[p]}")
    elif p and p not in PROCESSOR:
        errs.append(f"platform '{p}' unknown; supported: {', '.join(sorted(PROCESSOR))}")
    ref = str(row.get("reference", "")).strip()
    if ref and ref not in registry:
        errs.append(f"reference '{ref}' is not in references/_registry/registry.tsv - "
                    f"known: {', '.join(sorted(registry)) or '(registry empty)'}")
    if not str(row.get("matrix", "")).strip() and not str(row.get("fastq_r1", "")).strip():
        errs.append("neither 'matrix' nor 'fastq_r1' given - nothing to ingest")
    return errs

def plan_one(row: dict, registry: dict, stats_fn) -> Plan:
    """Decide how this sample is ingested. `stats_fn(path)` returns matrix summary stats."""
    sample = str(row.get("sample", "?")).strip()
    platform = str(row.get("platform", "")).strip().lower()
    proc = PROCESSOR.get(platform, "?")

    errs = validate_row(row, registry)
    if errs:
        return Plan(sample, platform, proc, "blocked", "; ".join(errs))

    matrix = str(row.get("matrix", "")).strip()
    has_fastq = bool(str(row.get("fastq_r1", "")).strip())

    if not matrix:
        return Plan(sample, platform, proc, "run",
                    "no matrix supplied - running the processor from FASTQ")

    if not Path(matrix).exists():
        if has_fastq:
            return Plan(sample, platform, proc, "run",
                        f"matrix '{matrix}' does not exist - falling back to FASTQ")
        return Plan(sample, platform, proc, "blocked",
                    f"matrix '{matrix}' does not exist and no FASTQ was given")

    st = stats_fn(matrix)
    v = verify(name=matrix, **st)
    if v.usable:
        return Plan(sample, platform, proc, "accept",
                    "supplied matrix verified: raw values AND raw droplets", v)
    if has_fastq:
        # NOT a warning. A matrix that fails is not usable, and continuing with it is how a
        # cell-called matrix reaches an ambient model that needs empty droplets.
        return Plan(sample, platform, proc, "run",
                    "supplied matrix REJECTED - rebuilding from FASTQ", v)
    return Plan(sample, platform, proc, "blocked",
                "supplied matrix REJECTED and no FASTQ to rebuild from", v)

def read_registry(path: Path) -> dict:
    if not path.exists():
        return {}
    out = {}
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines()):
        if not line.strip() or i == 0:
            continue
        f = line.split("\t")
        if len(f) >= 3:
            out[f"{f[0]}/{f[1]}"] = f[2]
    return out
