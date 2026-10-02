# Runs the real CLI end to end on small synthetic cohorts in temporary directories, with a stand-in
# for the doublet detector (tests/fixtures/fake_rscript.py); needs the analysis stack, not R.
"""The same pipeline under different conditions - the cohorts a reusable QC tool will be given.

Until 2026-10-02 scQC had run end to end on one cohort: one species, one platform, one assay, one
input route. Every defect the second cohort found hid because the first one's value equalled a
default or the first one never took the path (single-cell-harness ADR-0027). This runs the real
`scqc run`, evidence and apply, on synthetic cohorts that differ on the axes a cohort differs on,
and asserts each one completes - or is refused for its stated reason, never crashes:

  cells_human      single cells, human, 10x MatrixMarket, two libraries, evidence and apply
  cells_mouse      single cells, mouse gene naming (mt-, Rp[sl])
  cells_one_lib    a single library - no design, no sibling to compare
  cells_h5_adt     10x HDF5 input with antibody-capture features and a duplicated symbol
  cells_depths     one deep and one shallow library
  mixed_assays     a cohort that mixes cells and nuclei - refused, by name
  no_mito_genes    a reference whose symbols match no mitochondrial prefix - refused, by name

The detector is a stand-in, so nothing here says the doublet calls are right; it says that the
pipeline around them carries any honest cohort to a report, and to a deliverable where asked.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
fails: list[str] = []

try:
    import anndata as ad
    import numpy as np  # noqa: F401
    from fixtures import cohorts as C
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")

FAKE_R = ROOT / "tests" / "fixtures" / "fake_rscript.py"
COLS = ["sample", "platform", "species", "reference", "assay", "matrix", "aligner_cells",
        "genotype", "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]


def run(td: Path, rows, mode, registry):
    sheet = td / "samplesheet.tsv"
    sheet.write_text("\t".join(COLS) + "\n"
                     + "".join("\t".join(str(r.get(c, "")) for c in COLS) + "\n" for r in rows))
    reg = td / "registry.tsv"
    reg.write_text("species\tbuild\tpath\n" + "".join(f"{sp}\t{b}\t-\n" for sp, b in registry))
    proj = td / f"proj_{mode}"
    p = subprocess.run(
        [sys.executable, str(ROOT / "scqc_cli.py"), "run", "--project", str(proj),
         "--samplesheet", str(sheet), "--registry", str(reg), "--mode", mode,
         "--executor", "local", "--python", sys.executable, "--rscript", str(FAKE_R),
         "--jobs", "4", "--seed", "0"], capture_output=True, text=True, timeout=1800)
    res = [d for d in (proj / "results").iterdir()
           if d.is_dir() and not d.is_symlink() and d.name != "latest"] \
        if (proj / "results").exists() else []
    return p.returncode, p.stdout + p.stderr, (res[0] if res else None)


def cells(td, names, spec=C.HUMAN if HAVE else None, depths=None, **kw):
    rows = []
    for i, s in enumerate(names):
        extra = dict(kw)
        if depths:
            extra["cell_umi"] = depths[i]
        raw, filt, _ = C.cellranger_library(td / "data", s, 10 + i, spec, **extra)
        rows.append({"sample": s, "platform": "10x", "species": spec["species"],
                     "reference": f"{spec['species']}/ref", "assay": "scrna", "matrix": raw,
                     "aligner_cells": filt, "genotype": "WT" if i % 2 == 0 else "Mut",
                     "mt_prefix": spec["mt_prefix"], "ribo_pattern": spec["ribo_pattern"],
                     "dbr": 0.05, "dbr_sd": 0.05})
    return rows


def completes(label, td, rows, mode, registry, check=None):
    rc, log, R = run(td, rows, mode, registry)
    if rc != 0 or R is None:
        tail = "\n".join(l for l in log.splitlines() if l.strip())[-1200:]
        fails.append(f"{label}/{mode}: did not complete (exit {rc}):\n{tail}")
        return None
    st = json.loads((R / "STATUS.json").read_text())
    if mode == "apply" and not (R / "objects" / "cohort.deliverable.h5ad").exists():
        fails.append(f"{label}/apply: no deliverable written (status {st.get('status')})")
    if check:
        check(R)
    print(f"  {label:<14} {mode:<8} completed")
    return R


def refused(label, td, rows, registry, words):
    rc, log, _ = run(td, rows, "evidence", registry)
    if rc == 0:
        fails.append(f"{label}: completed; it should have been refused ({words!r})")
    elif "Traceback" in log and words not in log:
        fails.append(f"{label}: crashed rather than refusing:\n{log[-800:]}")
    elif words not in log:
        fails.append(f"{label}: refused, but not for the stated reason {words!r}:\n{log[-600:]}")
    else:
        print(f"  {label:<14} refused   ({words})")


HUMAN_REG = [("homo_sapiens", "ref")]
MOUSE_REG = [("mus_musculus", "ref")]

if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        rows = cells(td, ["WT_1", "MUT_1"])
        completes("cells_human", td, rows, "evidence", HUMAN_REG)

        def human_apply(R):
            d = ad.read_h5ad(R / "objects" / "cohort.deliverable.h5ad", backed="r")
            if d.obs_names.has_duplicates:
                fails.append("cells_human/apply: the deliverable has duplicate cell names")
            if "barcode" not in d.obs.columns:
                fails.append("cells_human/apply: the aligner's barcode was not kept")
            d.file.close()
            with open(R / "tables" / "removal_by_criterion.csv", encoding="utf-8") as fh:
                crit = {r["criterion"] for r in csv.DictReader(fh)}
            if "fail_not_aligner_cell" not in crit or "fail_not_cellbender_cell" in crit:
                fails.append(f"cells_human/apply: removal criteria {sorted(crit)}")
        completes("cells_human", td, rows, "apply", HUMAN_REG, human_apply)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        completes("cells_mouse", td, cells(td, ["WT_1", "MUT_1"], C.MOUSE), "apply", MOUSE_REG)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        completes("cells_one_lib", td, cells(td, ["ONLY_1"]), "apply", HUMAN_REG)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        def no_antibodies(R):
            # The QC object holds gene expression only: an antibody count is not a transcript,
            # and counted as one it inflates every cell's depth and dilutes its mito share.
            o = ad.read_h5ad(R / "objects" / "WT_1_ambient.h5", backed="r")
            ft = o.var["feature_types"].astype(str).tolist() if "feature_types" in o.var else []
            if any(t != "Gene Expression" for t in ft) or any("TotalSeq" in str(v)
                                                               for v in o.var_names):
                fails.append("cells_h5_adt: antibody-capture features reached the QC object")
            o.file.close()
        completes("cells_h5_adt", td, cells(td, ["WT_1", "MUT_1"], fmt="h5", n_adt=8,
                                            duplicate=True),
                  "apply", HUMAN_REG, no_antibodies)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        completes("cells_depths", td,
                  cells(td, ["DEEP_1", "SHALLOW_1"], depths=[(4000, 12000), (700, 1800)]),
                  "evidence", HUMAN_REG)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        rows = cells(td, ["WT_1", "MUT_1"])
        rows[1]["assay"] = "snrna"
        refused("mixed_assays", td, rows, HUMAN_REG, "mixes assays")

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        rows = cells(td, ["WT_1", "MUT_1"])
        for r in rows:
            r["mt_prefix"] = "MITO-"
        refused("no_mito_genes", td, rows, HUMAN_REG, "MITO-")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("every condition completes, or is refused for its stated reason")
