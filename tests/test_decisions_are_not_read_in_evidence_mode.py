# Runs the real CLI on a one-library synthetic single-cell cohort in a temporary directory, with the
# stand-in doublet detector; needs the analysis stack, not R.
"""An evidence run given a decisions file is not refused: the file is not read, and the run says so
(harness autodev W020, B1-051).

`scqc run --mode evidence --decisions d.yml` stopped before opening any file ("--decisions is only
read in apply mode"). Evidence mode applies nothing, so the refusal protected no result, and a host
passing --decisions to every run of a project was stopped by it. Checked here, through the real CLI:
the run exits 0, prints `decisions not read in evidence mode`, records it in STATUS.json `notes`,
and an override in the file - an absurd UMI floor - changes no threshold the run derived.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
fails: list[str] = []

try:
    import anndata  # noqa: F401
    from fixtures import cohorts as C
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")

COLS = ["sample", "platform", "species", "reference", "assay", "matrix", "aligner_cells",
        "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]

if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        raw, filt, _ = C.cellranger_library(td / "data", "ONLY_1", 5, C.HUMAN)
        row = {"sample": "ONLY_1", "platform": "10x", "species": "homo_sapiens",
               "reference": "homo_sapiens/ref", "assay": "scrna", "matrix": raw,
               "aligner_cells": filt, "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]",
               "dbr": 0.05, "dbr_sd": 0.05}
        sheet = td / "samplesheet.tsv"
        sheet.write_text("\t".join(COLS) + "\n" + "\t".join(str(row[c]) for c in COLS) + "\n")
        (td / "registry.tsv").write_text("species\tbuild\tpath\nhomo_sapiens\tref\t-\n")
        dec = td / "decisions.yml"
        dec.write_text("quality:\n  umi_floor:\n    value: 1999\n    approved_by: someone\n"
                       "    verbatim: a floor nobody would choose, to show it is not read\n")
        p = subprocess.run(
            [sys.executable, str(ROOT / "scqc_cli.py"), "run", "--project", str(td / "proj"),
             "--samplesheet", str(sheet), "--registry", str(td / "registry.tsv"),
             "--mode", "evidence", "--decisions", str(dec), "--executor", "local",
             "--python", sys.executable,
             "--rscript", str(ROOT / "tests" / "fixtures" / "fake_rscript.py"),
             "--jobs", "4", "--seed", "0"], capture_output=True, text=True, timeout=1800)
        out = p.stdout + p.stderr
        if p.returncode != 0:
            fails.append(f"evidence with --decisions exited {p.returncode}: {out[-400:]}")
        if "decisions not read in evidence mode" not in out:
            fails.append("the run did not say the decisions file was not read")
        res = [d for d in (td / "proj" / "results").iterdir()
               if d.is_dir() and not d.is_symlink()] if (td / "proj" / "results").exists() else []
        if res:
            st = json.loads((res[0] / "STATUS.json").read_text())
            if not any("decisions not read in evidence mode" in n for n in st.get("notes") or []):
                fails.append(f"STATUS.json does not record the note: {st.get('notes')}")
            pay = json.loads((res[0] / "reports" / "payload.json").read_text())
            if "1999" in json.dumps(pay.get("parameters") or []):
                fails.append("the decisions file's floor reached the run's parameters")
        print("an evidence run given --decisions completes, says the file was not read, and "
              "records it")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("decisions are not read in evidence mode, and the run says so instead of stopping")
