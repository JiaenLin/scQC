# Builds Pipeline objects in temporary directories, then runs the real CLI twice AT ONCE on two
# synthetic single-cell projects that share a folder name; a stand-in replaces the doublet detector
# (tests/fixtures/fake_rscript.py), so no R is needed.
"""Two runs never share a scratch directory.

The local executor put its scratch at `$TMPDIR/scqc_<project folder name>`. Two runs whose project
folders shared a NAME - two copies of the test suite, each building `proj_apply` in its own
temporary directory - therefore wrote the same intermediates at the same paths at the same time,
and failed each other (harness autodev W001's verifier: two suites at once, test_conditions failed;
alone 8 of 8). The name carried the folder and nothing that told two runs apart.

  A  the scratch of two projects with one folder name differs, under $TMPDIR and under
     SCQC_SCRATCH; two runs of one project with different run keys differ too; under a scheduler
     the scratch stays inside the run's own work directory, as before
  B  the CLI, twice at once, on two same-named projects: both complete, and each wrote its own
     deliverable
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
fails: list[str] = []

from engine.pipeline import Pipeline  # noqa: E402


class _Local:
    name = "local"


class _Pbs:
    name = "pbs"


ROWS = [{"sample": "S1", "platform": "10x", "species": "homo_sapiens",
         "reference": "homo_sapiens/ref", "assay": "scrna", "matrix": "/d/S1"}]


def scratch(project: Path, executor, rows=ROWS, mode="evidence") -> Path:
    return Pipeline(project=project, mode=mode, executor=executor, samples=rows).scratch


# ---- A. the scratch directory's name
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    saved = {k: os.environ.get(k) for k in ("TMPDIR", "SCQC_SCRATCH")}
    try:
        os.environ["TMPDIR"] = str(td / "tmp")
        (td / "tmp").mkdir()
        os.environ.pop("SCQC_SCRATCH", None)
        a = scratch(td / "one" / "proj_apply", _Local())
        b = scratch(td / "two" / "proj_apply", _Local())
        if a == b:
            fails.append(f"A1: two projects named proj_apply share the scratch {a}")
        if a.parent != td / "tmp" or not a.name.startswith("scqc_proj_apply"):
            fails.append(f"A1: the scratch left $TMPDIR or lost the project's name: {a}")
        c = scratch(td / "one" / "proj_apply", _Local(), mode="apply")
        if c == a:
            fails.append("A2: two runs of one project with different run keys share a scratch")
        os.environ["SCQC_SCRATCH"] = str(td / "shared")
        d = scratch(td / "one" / "proj_apply", _Local())
        e = scratch(td / "two" / "proj_apply", _Local())
        if d == e or d.parent != td / "shared":
            fails.append(f"A3: under SCQC_SCRATCH two same-named projects share {d}")
        os.environ.pop("SCQC_SCRATCH", None)
        p = Pipeline(project=td / "one" / "proj_apply", mode="evidence", executor=_Pbs(),
                     samples=ROWS)
        if p.scratch.parent != p.work:
            fails.append(f"A4: under a scheduler the scratch left the run's work directory: "
                         f"{p.scratch}")
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
print("A. scratch is named for the run: two same-named projects, and two run keys, never meet")

# ---- B. the CLI, twice at once
try:
    import anndata as ad  # noqa: F401
    from fixtures import cohorts as C
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")

FAKE_R = ROOT / "tests" / "fixtures" / "fake_rscript.py"
COLS = ["sample", "platform", "species", "reference", "assay", "matrix", "aligner_cells",
        "genotype", "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]

if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        shared_tmp = td / "tmp"
        shared_tmp.mkdir()
        procs = []
        for side in ("left", "right"):
            base = td / side
            rows = []
            for i, s in enumerate(("WT_1", "MUT_1")):
                raw, filt, _ = C.cellranger_library(base / "data", s, 20 + i, C.HUMAN)
                rows.append({"sample": s, "platform": "10x", "species": "homo_sapiens",
                             "reference": "homo_sapiens/ref", "assay": "scrna", "matrix": raw,
                             "aligner_cells": filt, "genotype": "WT" if i == 0 else "Mut",
                             "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]", "dbr": 0.05,
                             "dbr_sd": 0.05})
            sheet = base / "samplesheet.tsv"
            sheet.write_text("\t".join(COLS) + "\n" + "".join(
                "\t".join(str(r[c]) for c in COLS) + "\n" for r in rows))
            reg = base / "registry.tsv"
            reg.write_text("species\tbuild\tpath\nhomo_sapiens\tref\t-\n")
            env = dict(os.environ, TMPDIR=str(shared_tmp))
            env.pop("SCQC_SCRATCH", None)
            # The SAME project folder name on both sides - the condition that collided.
            procs.append((side, base / "proj_apply", subprocess.Popen(
                [sys.executable, str(ROOT / "scqc_cli.py"), "run",
                 "--project", str(base / "proj_apply"), "--samplesheet", str(sheet),
                 "--registry", str(reg), "--mode", "apply", "--executor", "local",
                 "--python", sys.executable, "--rscript", str(FAKE_R), "--jobs", "4",
                 "--seed", "0"],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)))
        for side, proj, pr in procs:
            out, _ = pr.communicate(timeout=1800)
            res = [d for d in (proj / "results").iterdir()
                   if d.is_dir() and not d.is_symlink() and d.name != "latest"] \
                if (proj / "results").exists() else []
            if pr.returncode != 0 or not res:
                fails.append(f"B: the {side} run did not complete (exit {pr.returncode}):\n"
                             + "\n".join(l for l in out.splitlines() if l.strip())[-1200:])
            elif not (res[0] / "objects" / "cohort.deliverable.h5ad").exists():
                fails.append(f"B: the {side} run wrote no deliverable")
        dirs = sorted(p.name for p in shared_tmp.iterdir() if p.name.startswith("scqc_"))
        if len(dirs) != 2:
            fails.append(f"B: two runs at once made {len(dirs)} scratch director(ies) under "
                         f"$TMPDIR: {dirs}")
    print("B. two runs at once, two same-named projects: both complete, each in its own scratch")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("two runs never share a scratch directory")
