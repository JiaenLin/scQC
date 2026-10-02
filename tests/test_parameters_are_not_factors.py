# Discovers a design from samplesheet rows and runs one installer function on a temp directory.
"""Two defects a second cohort found before its first run (single-cell-harness ADR-0027).

Q7  Design discovery admits any column with 2..n-1 distinct values. It skipped what the run IS
    (sample, platform, species, ...) but not this pipeline's own per-library PARAMETERS - so a
    doublet rate declared per library, with two libraries sharing a value, became a "factor",
    and every differential check compared libraries by the parameter they were run with.
Q8  The environment installer skipped any environment directory that existed. A build killed
    mid-transaction leaves one with no interpreter, and every rerun then printed "[skip] core
    already exists" and died on the missing bin/python with exit 127 and no reason.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
fails: list[str] = []

# ---- Q7
from engine import steps  # noqa: E402
rows = [{"sample": s, "platform": "10x", "species": "homo_sapiens", "reference": "r", "assay": "scrna",
         "matrix": f"/d/{s}", "genotype": g, "dbr": d, "dbr_sd": d, "mt_prefix": "MT-",
         "ribo_pattern": "^RP[SL]", "aligner_cells": f"/d/{s}_f"}
        for s, g, d in (("A1", "WT", "0.05"), ("A2", "WT", "0.06"), ("B1", "Mut", "0.05"), ("B2", "Mut", "0.06"))]
got = sorted(steps._design(rows))
if got != ["genotype"]:
    fails.append(f"Q7: the design discovered is {got}; a per-library parameter is not a factor")
print(f"Q7. design discovered from a sheet with per-library doublet rates: {got}")

# ---- Q8
fn = subprocess.run(["sed", "-n", "/^make_env() {/,/^}/p", str(ROOT / "setup" / "install_env.sh")],
                    capture_output=True, text=True).stdout
if "make_env" not in fn:
    fails.append("Q8: make_env() not found in setup/install_env.sh")
with tempfile.TemporaryDirectory() as td:
    (Path(td) / "core" / "conda-meta").mkdir(parents=True)        # what a killed build leaves
    script = f'set -e\n{fn}\nPREFIX="{td}"; CONDA=false; created=()\nmake_env core 3.12\necho MADE\n'
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                       env={"PATH": os.environ.get("PATH", "")})
    said = r.stdout + r.stderr
    if r.returncode == 0 or "interrupted build" not in said or "rm -rf" not in said:
        fails.append(f"Q8: a half-built environment was not refused with its remedy (exit "
                     f"{r.returncode}): {said.strip()[-160:]}")
print("Q8. an environment directory with no interpreter is refused, naming the remedy")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("parameters are not factors, and an interrupted build is not a built one")
