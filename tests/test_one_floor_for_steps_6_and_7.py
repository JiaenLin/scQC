# Builds the task graph and runs steps 6 and 7 against a stand-in for the analysis subprocess,
# then runs the real CLI in apply mode on a synthetic single-cell cohort with a decisions file; a
# stand-in replaces the doublet detector (tests/fixtures/fake_rscript.py), so no R is needed.
"""Steps 6 and 7 apply the same thresholds, whoever set them.

Step 6 clusters the population the deliverable will hold - quality-filtered, doublets attached and
not applied - and step 7 carries step 6's cluster flags onto every delivered cell. Step 6 read the
DERIVED floors and ceiling itself, while step 7 resolved them through `_apply_thresholds`, which
lets a decisions file adjudicate them. With a decisions file the two steps therefore applied
different thresholds: step 6 clustered and flagged cells the deliverable did not hold, and the
flags on the deliverable described a different population (single-cell-harness ADR-0027, open
from the audit; autodev gap GAP-024).

  A  the graph and the two steps, with the analysis subprocess stood in: the population step 6
     clusters and the thresholds step 7 measures are the same numbers - with an adjudicated
     floor and ceiling, with none, and with a decision missing its attestation
  B  step 6's task signature follows the decisions file, so a changed decision re-runs step 6
     rather than letting it be skipped as unchanged
  C  the CLI, apply mode, an adjudicated UMI floor above the declared one: no cell step 6
     clustered sits below it
"""
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
fails: list[str] = []

from engine import graph, steps  # noqa: E402

ATTESTED = {"approved_by": "a test operator",
            "verbatim": "apply this value to the synthetic cohort, read off its report"}
DECIDED = {"quality": {"umi_floor": {"value": 3000, **ATTESTED},
                       "gene_floor": {"value": 900, **ATTESTED},
                       "mito_ceiling_pct": {"value": 12.5, **ATTESTED}}}
# A value with no words is a number somebody typed: step 7 applies the DERIVED value instead.
UNATTESTED = {"quality": {"umi_floor": {"value": 3000, "approved_by": "a test operator"}}}


def _row(s):
    return {"sample": s, "platform": "10x", "species": "homo_sapiens",
            "reference": "homo_sapiens/GRCh38", "assay": "scrna", "matrix": f"/d/{s}/raw",
            "aligner_cells": f"/d/{s}/filtered", "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]",
            "genotype": "REF" if s.endswith("1") else "ALT"}


ROWS = [_row("A1"), _row("B2")]
CEILINGS = {"A1": 14.25, "B2": 16.5}


class _Done:
    ok = True

    def __init__(self, metrics):
        self.metrics = metrics


class _P:
    """The attributes the graph and steps 6-7 read, over a real tables directory."""

    def __init__(self, td: Path, decisions):
        self.results = td / "results"
        self.work = td / "work"
        self.scratch = td / "scratch"
        self.project = td / "proj"
        self.mode = "apply"
        self.samples = ROWS
        self.decisions = decisions
        self.results_by_key = {"05_quality": _Done({
            "umi_proposed": 500, "genes_proposed": 250,
            "umi_provenance": "declared", "genes_provenance": "declared"})}
        (self.results / "tables").mkdir(parents=True, exist_ok=True)
        with open(self.results / "tables" / "mito_ceiling_per_sample.csv", "w",
                  newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["sample", "ceiling"])
            for s, c in CEILINGS.items():
                w.writerow([s, f"{c:.6f}"])


def _steps_6_and_7(td: Path, decisions):
    """{sample: (what step 6 clusters on, what step 7 measures against)}."""
    p = _P(td, decisions)
    tasks = {t.key: t for t in graph.main_stage(p, "python", {}, {r["sample"]: {"mode": "accept"}
                                                                   for r in ROWS})}
    seen = []
    real = steps._scanpy
    steps._scanpy = lambda pipeline, op, h5, prefix, params, log, python_exe: (
        seen.append((op, params)) or {"outputs": [], "metrics": {}})
    try:
        out = {}
        for s in ("A1", "B2"):
            seen.clear()
            tasks[f"06_cluster/{s}"].fn(task=tasks[f"06_cluster/{s}"], pipeline=p, log=None)
            pop = next(pr["population"] for op, pr in seen if op == "cluster")
            seen.clear()
            tasks[f"07_measure/{s}"].fn(task=tasks[f"07_measure/{s}"], pipeline=p, log=None)
            meas = next(pr for op, pr in seen if op == "apply_measure")
            out[s] = ((float(pop["umi_floor"]), float(pop["gene_floor"]),
                       float(pop["mito_ceiling"])),
                      (float(meas["umi_floor"]), float(meas["gene_floor"]),
                       float(meas["mito_ceiling_pct"])))
        return out, tasks
    finally:
        steps._scanpy = real


# ---- A. one set of thresholds
for label, decisions, expect in (
        ("an adjudicated floor and ceiling", DECIDED,
         {s: (3000.0, 900.0, 12.5) for s in CEILINGS}),
        ("no decisions file", {}, {s: (500.0, 250.0, c) for s, c in CEILINGS.items()}),
        ("a decision without its words", UNATTESTED,
         {s: (500.0, 250.0, c) for s, c in CEILINGS.items()})):
    with tempfile.TemporaryDirectory() as td:
        got, _ = _steps_6_and_7(Path(td), decisions)
    for s, (six, seven) in got.items():
        if six != seven:
            fails.append(f"A: {label}, {s}: step 6 clusters on (umi, genes, mito) {six}, step 7 "
                         f"applies {seven}")
        elif six != expect[s]:
            fails.append(f"A: {label}, {s}: both steps apply {six}; expected {expect[s]}")
print("A. step 6 clusters on the thresholds step 7 applies: adjudicated, absent, unattested")

# ---- B. a changed decision re-runs step 6
#
# The same directory both times, so the decisions file is the only thing that differs: a task's
# signature hashes its parameters, and every path in them would differ between two directories.
with tempfile.TemporaryDirectory() as td:
    _, t1 = _steps_6_and_7(Path(td), DECIDED)
    changed = {"quality": {**DECIDED["quality"], "umi_floor": {"value": 2000, **ATTESTED}}}
    _, t2 = _steps_6_and_7(Path(td), changed)
    _, t3 = _steps_6_and_7(Path(td), DECIDED)
    if t1["06_cluster/A1"].signature() != t3["06_cluster/A1"].signature():
        fails.append("B: step 6's signature changes with nothing changed; this check cannot "
                     "tell a changed decision from noise")
    elif t1["06_cluster/A1"].signature() == t2["06_cluster/A1"].signature():
        fails.append("B: step 6's signature does not change with the decisions file, so a "
                     "changed floor would leave step 6 skipped on the old population")
print("B. step 6's signature follows the decisions file")

# ---- C. the CLI
try:
    import anndata  # noqa: F401
    from fixtures import cohorts as C
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")

FAKE_R = ROOT / "tests" / "fixtures" / "fake_rscript.py"
COLS = ["sample", "platform", "species", "reference", "assay", "matrix", "aligner_cells",
        "genotype", "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]
FLOOR = 3000

if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        lines = []
        for i, s in enumerate(("REF_1", "ALT_1")):
            raw, filt, _ = C.cellranger_library(td / "data", s, 30 + i, C.HUMAN)
            lines.append("\t".join(str(v) for v in (
                s, "10x", "homo_sapiens", "homo_sapiens/ref", "scrna", raw, filt,
                "REF" if i == 0 else "ALT", "MT-", "^RP[SL]", 0.05, 0.05)))
        sheet = td / "samplesheet.tsv"
        sheet.write_text("\t".join(COLS) + "\n" + "\n".join(lines) + "\n")
        reg = td / "registry.tsv"
        reg.write_text("species\tbuild\tpath\nhomo_sapiens\tref\t-\n")
        dec = td / "decisions.yml"
        dec.write_text("quality:\n  umi_floor:\n    value: %d\n    approved_by: a test operator\n"
                       "    verbatim: apply %d UMI as the floor for this synthetic cohort\n"
                       % (FLOOR, FLOOR))
        proj = td / "proj"
        p = subprocess.run(
            [sys.executable, str(ROOT / "scqc_cli.py"), "run", "--project", str(proj),
             "--samplesheet", str(sheet), "--registry", str(reg), "--mode", "apply",
             "--decisions", str(dec), "--executor", "local", "--python", sys.executable,
             "--rscript", str(FAKE_R), "--jobs", "4", "--seed", "0"],
            capture_output=True, text=True, timeout=1800)
        res = [d for d in (proj / "results").iterdir()
               if d.is_dir() and not d.is_symlink() and d.name != "latest"] \
            if (proj / "results").exists() else []
        if p.returncode != 0 or not res:
            fails.append(f"C: the apply run did not complete (exit {p.returncode}):\n"
                         + (p.stdout + p.stderr)[-1500:])
        else:
            T = res[0] / "tables"
            for s in ("REF_1", "ALT_1"):
                with open(T / f"{s}.percell.csv", encoding="utf-8") as fh:
                    depth = {r["barcode"]: float(r["total_counts"]) for r in csv.DictReader(fh)}
                with open(T / f"{s}.embedding.csv", encoding="utf-8") as fh:
                    clustered = [r["barcode"] for r in csv.DictReader(fh)
                                 if r["clustered"] == "True"]
                below = [b for b in clustered if depth.get(b, 0.0) < FLOOR]
                if not clustered:
                    fails.append(f"C: {s}: step 6 clustered nothing")
                elif not any(v < FLOOR for v in depth.values()):
                    fails.append(f"C: {s}: no called cell lies below {FLOOR} UMI, so this "
                                 f"check could not fail")
                elif below:
                    fails.append(f"C: {s}: step 6 clustered {len(below)} of {len(clustered)} "
                                 f"cells below the adjudicated {FLOOR}-UMI floor step 7 "
                                 f"applies")
        print(f"C. the CLI: no cell step 6 clustered sits below the adjudicated floor")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("steps 6 and 7 apply one set of thresholds, whoever set them")
