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

ON SINGLE CELLS A CELL-CALLED MATRIX IS ACCEPTED, AND SAID TO BE ONE (harness autodev W020, B1-046).
No ambient model runs on cells and nothing in their route reads an empty droplet, so P2 is required
on `snrna` only (`verify_raw.P2_REQUIRED`). A single-cell matrix that fails P2 alone - the authors'
filtered cells, which is how most public single-cell data is published - is ACCEPTED as a
cell-called input: the plan says so, the run records it, and nothing downstream treats its barcodes
as droplets. A FASTQ given beside it is not used: rebuilding the droplets would buy empties that
nothing on this route reads.

WHAT IS DECLARED AND WHAT IS NOT

`platform`, `species` and `reference` are DECLARED with no defaults, because the pipeline cannot
infer them and a wrong guess is silent. The platform names the aligner that would rebuild a matrix
from FASTQ, so it must be one this pipeline can align only on that route; a verified supplied matrix
from any declared platform is accepted, and the record says no processor was involved. The doublet-rate model is platform-specific: applying one
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
from verify_raw import CALLED_INPUT_NOTE, verify # noqa: E402

REQUIRED = ("sample", "platform", "species", "reference", "assay")
#: The assays this pipeline has a policy for. `modules/01_ambient` (ASSAYS) and `modules/05_quality`
#: (MITO_BOUNDS) key on the same words, and tests/test_ingest.py holds the three equal.
#:
#: REQUIRED SINCE 2026-10-02. It was optional, and the mitochondrial step read it as
#: `params.get("assay", "snrna")` from a task never handed it - so every cohort was bounded as
#: single-nucleus whatever its samplesheet said, and no test could tell, because every cohort this
#: pipeline had run WAS single-nucleus. A cell or a nucleus is a fact about the experiment; it has
#: no default.
ASSAYS = ("snrna", "scrna")
#: The aligner that REBUILDS a matrix from FASTQ, per platform. Read on the FASTQ route only
#: (harness autodev W020, B1-048): a verified supplied matrix needs no aligner, so any declared
#: platform is accepted with one and recorded `processor: none (supplied matrix)` when it has none.
PROCESSOR = {"10x": "cellranger", "singleron": "celescope"}
FUTURE = {"bgi": "not implemented - declare it and the run will refuse rather than guess"}
#: What the record says when a supplied matrix was accepted on a platform this pipeline cannot
#: align: nothing was, or would have been, rebuilt.
NO_PROCESSOR = "none (supplied matrix)"


def _processor_refusal(p: str) -> str:
    """Why the FASTQ route cannot run for platform `p`, or '' when it can."""
    if p in FUTURE:
        return f"platform '{p}': {FUTURE[p]}"
    if p not in PROCESSOR:
        return (f"platform '{p}' has no aligner here ({', '.join(sorted(PROCESSOR))} do), so "
                f"nothing can be rebuilt from FASTQ - supply a verified matrix instead")
    return ""

@dataclass
class Plan:
    sample: str
    platform: str
    processor: str
    mode: str # "accept" | "run" | "blocked"
    reason: str = ""
    verdict: object = None

    @property
    def input_kind(self):
        """'raw_droplets', 'called_cells', or None where no supplied matrix was verified."""
        return None if self.verdict is None else self.verdict.input_kind

    def __str__(self) -> str:
        head = {"accept": "ACCEPT", "run": "RUN ", "blocked": "BLOCKED"}[self.mode]
        s = f"[{head}] {self.sample:12s} {self.platform:10s} -> {self.processor}"
        if self.reason:
            s += f"\n {self.reason}"
        if self.verdict is not None and (not self.verdict.usable or self.verdict.called_input):
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
    # THE PLATFORM SELECTS AN ALIGNER, AND ONLY THE FASTQ ROUTE RUNS ONE (B1-048). It is checked
    # here only where the route is already certain - no matrix, FASTQ given; a supplied matrix's
    # route is the verdict's, and plan_one checks the platform if that verdict sends it to FASTQ.
    if p and not str(row.get("matrix", "")).strip() and str(row.get("fastq_r1", "")).strip():
        why = _processor_refusal(p)
        if why:
            errs.append(why)
    a = str(row.get("assay", "")).strip().lower()
    if a and a not in ASSAYS:
        errs.append(f"assay '{a}' unknown; supported: {', '.join(ASSAYS)} - a cell and a nucleus "
                    f"get different mitochondrial bounds and a different ambient policy, and an "
                    f"assay outside these has neither")
    ref = str(row.get("reference", "")).strip()
    if ref and ref not in registry:
        errs.append(f"reference '{ref}' is not in the registry (the project's `--registry`) - "
                    f"known: {', '.join(sorted(registry)) or '(registry empty)'}")
    if not str(row.get("matrix", "")).strip() and not str(row.get("fastq_r1", "")).strip():
        errs.append("neither 'matrix' nor 'fastq_r1' given - nothing to ingest")
    return errs

def plan_one(row: dict, registry: dict, stats_fn) -> Plan:
    """Decide how this sample is ingested. `stats_fn(path)` returns matrix summary stats."""
    sample = str(row.get("sample", "?")).strip()
    platform = str(row.get("platform", "")).strip().lower()
    proc = PROCESSOR.get(platform, NO_PROCESSOR)

    errs = validate_row(row, registry)
    if errs:
        return Plan(sample, platform, proc, "blocked", "; ".join(errs))

    matrix = str(row.get("matrix", "")).strip()
    has_fastq = bool(str(row.get("fastq_r1", "")).strip())

    def fastq_route(why: str, v=None) -> Plan:
        # Every route to the aligner comes through here, so the platform is checked on each.
        refusal = _processor_refusal(platform)
        if refusal:
            return Plan(sample, platform, proc, "blocked", f"{why}; {refusal}", v)
        return Plan(sample, platform, proc, "run", why, v)

    if not matrix:
        return fastq_route("no matrix supplied - running the processor from FASTQ")

    if not Path(matrix).exists():
        if has_fastq:
            return fastq_route(f"matrix '{matrix}' does not exist - falling back to FASTQ")
        return Plan(sample, platform, proc, "blocked",
                    f"matrix '{matrix}' does not exist and no FASTQ was given")

    st = stats_fn(matrix)
    assay = str(row.get("assay", "")).strip().lower()
    v = verify(name=matrix, assay=assay, **st)
    if v.usable and v.called_input:
        # ACCEPTED AND NAMED. Not a pass of P2 - the verdict keeps p2_raw_droplets False - and not
        # a silent one: the plan, the run's input record and STATUS.json all say cell-called.
        return Plan(sample, platform, proc, "accept",
                    f"supplied matrix verified: raw values; CELL-CALLED, no empty droplets - "
                    f"accepted on {assay}, where nothing reads an empty droplet ({CALLED_INPUT_NOTE})"
                    + (". The FASTQ given was not used: rebuilding the droplets buys empties nothing "
                       "on this route reads" if has_fastq else ""), v)
    if v.usable:
        return Plan(sample, platform, proc, "accept",
                    "supplied matrix verified: raw values AND raw droplets", v)
    if has_fastq:
        # NOT a warning. A matrix that fails is not usable, and continuing with it is how a
        # cell-called matrix reaches an ambient model that needs empty droplets.
        return fastq_route("supplied matrix REJECTED - rebuilding from FASTQ", v)
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
