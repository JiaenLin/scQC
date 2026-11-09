"""What one observation is called in the report, by the cohort's declared assay.

Every caption, axis label and step text said "nucleus" - written when the only cohort was nuclei -
so a single-cell run's report described an object the run never had. Single cells are not single
nuclei (single-cell-harness RULES, rule 2), and the report is where a reader meets the difference
first. The words come from the samplesheet's declared `assay`, read once, and travel with the
figure data and the payload.

An assay that is missing, mixed or unknown is named neutrally - "barcode" - never as either. A run
can only reach a report with no single assay when it was refused at ingest or at the graph, and
`scqc report` rebuilds without a samplesheet; neither may be described as nuclei by default.
"""
from __future__ import annotations

#: (one, many) per declared assay.
UNIT = {"snrna": ("nucleus", "nuclei"), "scrna": ("cell", "cells")}
#: For a cohort with no single declared assay. Every observation is a barcode, whatever it holds.
NEUTRAL = ("barcode", "barcodes")


def cohort_assay(rows) -> str | None:
    """The one assay the rows declare, normalised as ingest normalises it, or None."""
    kinds = {str((r or {}).get("assay") or "").strip().lower() for r in (rows or [])}
    kinds.discard("")
    return kinds.pop() if len(kinds) == 1 else None


def unit_words(assay) -> tuple:
    """(one, many) for this assay; the neutral pair when it is not one this report knows."""
    return UNIT.get(str(assay or "").strip().lower(), NEUTRAL)


def fill(text: str, assay) -> str:
    """`text` with `{one}` and `{many}` replaced by this assay's words. Nothing else is touched,
    so a caption that carries a literal brace is not broken by being formatted."""
    one, many = unit_words(assay)
    return str(text).replace("{one}", one).replace("{many}", many)
