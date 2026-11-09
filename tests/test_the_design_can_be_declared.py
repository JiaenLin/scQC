# Runs the real CLI on a four-library synthetic single-cell cohort in a temporary directory, with the
# stand-in doublet detector calling each library's declared rate; needs the analysis stack, not R.
"""The design can be declared, and then only its factors can stop a run (harness autodev W020,
B1-104).

Design factors were discovered from the samplesheet and never declared: any column that varied -
a lane, a chemistry - fed the differential checks with the power to refuse the run. Checked here,
through the real CLI, on a cohort whose LANE carries a 3.5x doublet differential and whose genotype
carries none:

  1  with nothing declared, the lane differential REFUSES, as before (discovery is kept)
  2  with `--design-factors genotype`, the same differential is a REVIEW naming the lane as
     discovered and not declared, and the run completes
  3  a declared factor that is not a samplesheet column is refused by name
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
        "genotype", "lane", "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]


def run(td: Path, label: str, extra: list):
    proj = td / f"proj_{label}"
    p = subprocess.run(
        [sys.executable, str(ROOT / "scqc_cli.py"), "run", "--project", str(proj),
         "--samplesheet", str(td / "samplesheet.tsv"), "--registry", str(td / "registry.tsv"),
         "--mode", "evidence", "--executor", "local", "--python", sys.executable,
         "--rscript", str(ROOT / "tests" / "fixtures" / "fake_rscript.py"),
         "--jobs", "4", "--seed", "0", *extra], capture_output=True, text=True, timeout=1800)
    res = [x for x in (proj / "results").iterdir() if x.is_dir() and not x.is_symlink()] \
        if (proj / "results").exists() else []
    return p.returncode, p.stdout + p.stderr, (res[0] if res else None)


if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        rows = []
        for i, (s, g, lane, dbr) in enumerate((("WT_1", "WT", "L1", 0.02),
                                                ("MUT_1", "Mut", "L1", 0.02),
                                                ("WT_2", "WT", "L2", 0.07),
                                                ("MUT_2", "Mut", "L2", 0.07))):
            raw, filt, _ = C.cellranger_library(td / "data", s, 50 + i, C.HUMAN)
            rows.append({"sample": s, "platform": "10x", "species": "homo_sapiens",
                         "reference": "homo_sapiens/ref", "assay": "scrna", "matrix": raw,
                         "aligner_cells": filt, "genotype": g, "lane": lane, "mt_prefix": "MT-",
                         "ribo_pattern": "^RP[SL]", "dbr": dbr, "dbr_sd": dbr})
        (td / "samplesheet.tsv").write_text(
            "\t".join(COLS) + "\n" + "".join("\t".join(str(r[c]) for c in COLS) + "\n"
                                             for r in rows))
        (td / "registry.tsv").write_text("species\tbuild\tpath\nhomo_sapiens\tref\t-\n")

        rc, log, R = run(td, "discovered", [])
        if rc != 2:
            fails.append(f"1: with nothing declared the lane differential did not refuse (exit {rc})")
        print("1. nothing declared: the discovered lane's 3.5x differential refuses")

        rc, log, R = run(td, "declared", ["--design-factors", "genotype"])
        if rc != 0 or R is None:
            fails.append(f"2: with genotype declared the run did not complete (exit {rc}):\n"
                         + "\n".join(log.splitlines()[-10:]))
        else:
            text = (R / "reports" / "payload.json").read_text()
            if "was discovered from the samplesheet and is not a declared design factor" not in text:
                fails.append("2: the lane differential does not say it was discovered, not declared")
        print("2. genotype declared: the lane differential is REVIEW, and the run completes")

        rc, log, R = run(td, "absent", ["--design-factors", "genotype,diet"])
        if rc == 0 or "diet" not in log or "not a samplesheet column" not in log:
            fails.append(f"3: a declared factor with no column was not refused by name (exit {rc})")
        print("3. a declared factor that is not a column is refused by name")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the design can be declared; then only its factors refuse")
