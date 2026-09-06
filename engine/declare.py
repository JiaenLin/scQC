"""What this tool declares about itself, in one place, for a reader or a host that will not read
the code: what it needs, what it provides, what it was shown, what its numbers cannot establish,
and which version of the NUMBERS this is.

`scqc describe --json` prints it. The fields follow the plugin contract of the harness this tool
is orchestrated by (needs / provides / sees / cannot_show / state_version / gates / escapes), so
an adapter is a field mapping rather than a reading of the source.

STATE_VERSION versions the numbers, not the code. It goes up whenever the same samplesheet and
the same parameters would produce a different deliverable — a changed derivation, a changed
default, a changed criterion. The commit already keys the run directory (engine/runkey.py); this
is the coarser, human-legible statement that the numbers moved, and a cached materialisation
anywhere downstream is keyed on it.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

STATE_VERSION = 1

# What the pipeline was SHOWN beyond the matrices: the design, through the samplesheet's
# factor columns, which the differential checks are computed over. Never the labels — this
# tool runs before annotation and cannot see them.
SEES = ["design"]

NEEDS = ["matrix/{counts}", "column/{sample}"]

# One mask per criterion, never one opaque mask: every removed observation is paired with the
# criteria that fired (`adapters/apply.build_removal_record`), so a later split is a move and
# not a rewrite.
# ONE NAME PER REASON A BARCODE MAY BE SET ASIDE. Written as a literal tuple rather than inlined
# into the comprehension below so that it can be READ - by a person, and by `sch dev check`,
# which parses this file to answer "is the new criterion registered?". A registry built inside a
# comprehension is invisible to anything but the interpreter, and a check that cannot see it
# reports "no criteria" when it means "I cannot tell", which is a different and worse answer.
CRITERIA = ("fail_not_cellbender_cell", "fail_umi_floor", "fail_gene_floor",
            "fail_mito_ceiling", "fail_doublet", "fail_mito_nf")

PROVIDES = [f"mask/{c}" for c in CRITERIA] + [
    "column/total_counts", "column/n_genes", "column/pct_counts_mt", "column/pct_counts_ribo",
    "column/nuclear_fraction", "column/doublet_score", "column/doublet_class",
    "column/cluster_FLAG", "table/removal_ledger", "table/removal_by_criterion",
    "table/thresholds_per_sample", "table/*.percell", "object/*.filtered"]

# The gates, each a module that returns PASS / REVIEW / REFUSE from findings it computed with the
# same functions the CSV subcommands expose — one instrument, two consumers.
GATES = {
    # `verdicts` is what each gate can return, in the closed vocabulary; `as_written` is the word
    # the module itself uses, so a reader can grep for it. Two gates refuse by raising rather than
    # by returning a word, and one writes FAIL where the vocabulary says REFUSE — declared, not
    # papered over: the test asserts the words are in the modules.
    "00_ingest": {"module": "modules/00_ingest/ingest.py", "probe": "scqc verify",
                  "verdicts": ["PASS", "REFUSE"], "as_written": ["REFUSES"]},
    "01_ambient": {"module": "modules/01_ambient/audit_ambient.py", "probe": None,
                   "verdicts": ["PASS", "REVIEW", "REFUSE"], "as_written": ["REFUSE", "REVIEW", "ok"]},
    "02_cells": {"module": "modules/02_cells/cellcall_gate.py", "probe": "scqc gate-cells",
                 "verdicts": ["PASS", "REVIEW", "REFUSE"], "as_written": ["REFUSE", "REVIEW"]},
    "04_doublets": {"module": "modules/04_doublets/doublet_health.py", "probe": "scqc doublet-health",
                    "verdicts": ["PASS", "REVIEW", "REFUSE"], "as_written": ["REFUSE", "REVIEW"]},
    "05_quality": {"module": "modules/05_quality/quality.py", "probe": "scqc quality",
                   "verdicts": ["PASS", "REFUSE"], "as_written": ["ThresholdRefusal"]},
    "06_cluster_check": {"module": "modules/06_cluster_check/cluster_flags.py", "probe": "scqc cluster-preflight",
                         "verdicts": ["PASS", "REFUSE"], "as_written": ["ClusterRefusal"]},
    "07_apply": {"module": "modules/07_apply/audit_removal.py", "probe": None,
                 "verdicts": ["PASS", "REVIEW", "REFUSE"], "as_written": ["FAIL", "REVIEW", "ok"],
                 "note": "the module writes FAIL for REFUSE"},
}
VERDICTS = ["PASS", "REVIEW", "REFUSE"]

# How a refusal is lifted, and where the lift is recorded. There is no --force.
ESCAPES = [
    {"what": "a DERIVED threshold", "how": "decisions.yml with value, class, approved_by, verbatim",
     "recorded": "reports/payload.json parameters[].decided_by and .verbatim"},
    {"what": "a REVIEW finding", "how": "proceed; the finding stays in gates[]",
     "recorded": "reports/payload.json gates[] — no decision is recorded against it yet"},
    {"what": "the apply step", "how": "--mode evidence runs no removal",
     "recorded": "INPUTS.json described.mode"},
]

CANNOT_SHOW = [
    "An observation passing QC is not an observation that is intact; it is one whose summary "
    "statistics fell inside thresholds this run derived or was given.",
    "A threshold is a finding about ONE cohort's distributions. Nothing here establishes that it "
    "transfers to another cohort, and the calibration cohort is a single tissue on a single platform.",
    "A removal rate that is equal across arms of the design was not differential; it was not "
    "therefore harmless, and when the overall removal exceeds one third the 3x differential test "
    "is inert and says so.",
    "A called doublet is a barcode the detector scored high; without hashing or genotypes no "
    "call here is a verified doublet, and a barcode below the light floor was never scored.",
    "The denoised matrix is a checkpoint: every check on it is cohort-relative and none can tell "
    "a well-removed ambient transcript from a badly-removed real one.",
]


def cannot_establish_by_step() -> dict:
    """The per-step 'cannot establish' text, read from engine/steps.py, so the two never drift."""
    try:
        from . import steps
        return {k: v[1] for k, v in steps.STEP_TEXT.items()}
    except Exception:  # noqa: BLE001 — a host without the step modules still gets the declaration
        return {}


def apply_criteria() -> list:
    src = (ROOT / "adapters" / "scanpy_ops.py").read_text(encoding="utf-8", errors="replace")
    m = re.search(r"APPLY_CRITERIA\s*=\s*\((.*?)\)", src, re.S)
    return re.findall(r"\"([a-z_]+)\"", m.group(1)) if m else []


def describe() -> dict:
    from .provenance import git_provenance
    version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    g = git_provenance(ROOT) or {}
    return {
        "contract": "1.0", "profile": "single-cell/1.0", "name": "scqc", "version": version,
        "commit": g.get("commit"), "dirty": g.get("dirty"), "state_version": STATE_VERSION,
        "class": "method", "layer": "stack", "reversible": True,
        "summary": "quality-control decisions and gates for single-cell and single-nuclei RNA-seq; "
                   "thresholds are findings, not arguments",
        "needs": NEEDS, "provides": PROVIDES, "sees": SEES,
        "criteria": apply_criteria(), "gates": GATES, "verdicts": VERDICTS, "escapes": ESCAPES,
        "cannot_show": CANNOT_SHOW, "cannot_establish": cannot_establish_by_step(),
        "status_files": ["STATUS.json", "RUNNING.txt", "SEALED.txt", "FAILED.txt"],
        "products": {"masks": "tables/<sample>.percell.csv, one boolean column per criterion, every barcode",
                     "ledger": "tables/removal_ledger.csv, removed barcodes with the criteria that fired",
                     "report": "reports/report.json and reports/payload.json",
                     "objects": "objects/<sample>.filtered.h5ad with uns['scqc']"},
        "language": "python", "entry": "scqc_cli.py", "needs_env": True,
    }


def dumps() -> str:
    return json.dumps(describe(), indent=1, sort_keys=True, default=str)
