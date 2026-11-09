# Runs the real CLI - `validate`, and `run` in evidence mode - on a one-library synthetic cohort in
# a temporary directory, with and without a project registry; needs the analysis stack, not R.
"""A project brings its own reference registry; the tool supplies none.

scQC shipped `references/_registry/registry.tsv` with one row - the first cohort's mouse
reference - and read it whenever a project gave no `--registry`, as `validate` did too. A project
that brought no registry was therefore resolved against that cohort's reference without a word:
a first-cohort value surviving as a default (harness autodev gap GAP-025). The samplesheet below
names exactly that reference, so on the old code every case here passed silently.

  A  `scqc validate` without `--registry` refuses, naming the flag
  B  `scqc validate --registry` with the project's own file accepts the row
  C  `scqc run` without `--registry` is refused at ingest - every library, naming the flag -
     and nothing after ingest runs; STATUS.json says refused
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
# The reference the tool used to ship, and fall back to.
REF_SPECIES, REF_BUILD = "mus_musculus", "ensembl_112_filtered"


def cli(*args):
    p = subprocess.run([sys.executable, str(ROOT / "scqc_cli.py"), *args],
                       capture_output=True, text=True, timeout=1800)
    return p.returncode, p.stdout + p.stderr


if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        raw, filt, _ = C.cellranger_library(td / "data", "ONLY_1", 7, C.MOUSE)
        sheet = td / "samplesheet.tsv"
        sheet.write_text("\t".join(COLS) + "\n" + "\t".join(str(v) for v in (
            "ONLY_1", "10x", REF_SPECIES, f"{REF_SPECIES}/{REF_BUILD}", "scrna", raw, filt,
            "mt-", "^Rp[sl]", 0.05, 0.05)) + "\n")
        reg = td / "registry.tsv"
        reg.write_text(f"species\tbuild\tpath\n{REF_SPECIES}\t{REF_BUILD}\t-\n")

        rc, out = cli("validate", "--samplesheet", str(sheet))
        if rc == 0 or "--registry" not in out:
            fails.append(f"A: validate without --registry: exit {rc}; it must refuse, naming "
                         f"the flag:\n{out[-500:]}")
        else:
            print("  A  validate without a project registry refuses, naming --registry")

        rc, out = cli("validate", "--samplesheet", str(sheet), "--registry", str(reg))
        if rc != 0:
            fails.append(f"B: validate with the project's registry: exit {rc}:\n{out[-500:]}")
        else:
            print("  B  validate with the project's own registry accepts the row")

        proj = td / "proj"
        rc, out = cli("run", "--project", str(proj), "--samplesheet", str(sheet),
                      "--mode", "evidence", "--executor", "local", "--python", sys.executable,
                      "--jobs", "1", "--seed", "0")
        res = [d for d in (proj / "results").iterdir()
               if d.is_dir() and not d.is_symlink() and d.name != "latest"] \
            if (proj / "results").exists() else []
        status = (json.loads((res[0] / "STATUS.json").read_text()).get("status")
                  if res and (res[0] / "STATUS.json").exists() else None)
        if rc != 2 or "--registry" not in out or "00_ingest" not in out:
            fails.append(f"C: run without --registry: exit {rc}; it must be refused at ingest, "
                         f"naming the flag:\n{out[-800:]}")
        elif "01_ambient" in out or "stage 2" in out:
            fails.append(f"C: run without --registry went past ingest:\n{out[-800:]}")
        elif status != "refused":
            fails.append(f"C: STATUS.json says {status!r}, not refused")
        else:
            print("  C  run without a project registry is refused at ingest; STATUS refused")

# ---- D. a new project is scaffolded with its own registry, and a template that declares nothing
#
# `setup/init_project.sh` wrote no registry, so a new project had none to pass. And its decisions
# template pre-filled step 6's three thresholds with one calibration cohort's numbers: a template
# copied to decisions.yml unedited replaced what step 6 proposes from the project's own clusters
# with that cohort's values. Every key of it that the code reads must be empty.
import shutil  # noqa: E402

try:
    import yaml  # noqa: F401 - the decisions template is YAML
    _YAML = True
except ImportError:
    _YAML = False

if shutil.which("bash") and _YAML:
    sys.path.insert(0, str(ROOT))
    from engine.decisions import load as _load_decisions  # noqa: E402

    with tempfile.TemporaryDirectory() as td:
        proj = Path(td) / "new"
        p = subprocess.run(["bash", str(ROOT / "setup" / "init_project.sh"), "--dir", str(proj),
                            "--assay", "scrna", "--samples", "2"],
                           capture_output=True, text=True, timeout=120)
        reg = proj / "registry.tsv"
        bad = len(fails)
        if p.returncode != 0:
            fails.append(f"D: init_project.sh failed (exit {p.returncode}): {p.stderr[-300:]}")
        else:
            if not reg.exists() or not reg.read_text().startswith("species\tbuild\tpath"):
                fails.append("D: init_project.sh wrote no registry.tsv with the registry's header")
            tpl = _load_decisions(proj / "decisions.template.yml")
            cc = tpl.get("cluster_check") or {}
            set_ = {k: cc.get(k) for k in ("a_umi_fraction", "b_mito_pct", "c_uninformative_pct")
                    if cc.get(k) is not None}
            if set_:
                fails.append(f"D: the decisions template declares step 6's thresholds before "
                             f"any evidence: {set_}")
        if len(fails) == bad:
            print("  D  a new project gets its registry's header and a template that declares "
                  "nothing")
else:
    print("  D  not checked: no bash, or no PyYAML, on this machine")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("a project brings its own registry: without one, validate and run both refuse")
