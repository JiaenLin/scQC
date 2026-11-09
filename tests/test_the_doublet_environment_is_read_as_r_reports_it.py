# Runs conf/env/install_rdoublet.sh under bash in temporary directories, against a stand-in package manager
# and a stand-in Rscript written there; installs nothing, downloads nothing, needs only the standard library.
"""The R doublet environment is checked as R reports it, and the conda build as conda recorded it.

conda-forge's r-xgboost 1.7.6 installs the CRAN release of the R package, which R reports as 1.7.6.1 (its
DESCRIPTION, in builds cpu_r43had0c348_6 and cpu_r43h6728c87_0), and 1.7.6.1 is what the calibration run's own
report records for the environment that scored its doublets. The recipe compared R's report with the conda
version, 1.7.6, so it refused every correct environment and `setup/install_env.sh --with-doublet` exited 1
(harness autodev W014's verifier, D1). An existing environment it refused was not told how to recover (D4).

R now only REPORTS what it holds; the recipe compares, so the comparison reads the same whoever reports - the
environment's R, or the stand-in here, which prints the versions a case gives it:

  1  --check names the conda specs - r-xgboost 1.7.6 as a cpu build, r-matrix 1.6_5 - and the versions R
     must report: xgboost 1.7.6.1, Matrix 1.6.5
  2  a fresh build: the manager is asked for exactly the pinned specs, and an environment whose R reports
     xgboost 1.7.6.1 and Matrix 1.6.5 is accepted
  3  an existing environment whose R reports xgboost 2.0.0 is refused, naming the shim and the fix (remove
     it and re-run, or build into a new prefix), and is not modified
  4  an existing environment whose R reports scDblFinder 1.14.0, or Matrix 1.6.4, is refused
  5  an existing environment whose conda record holds a GPU build of r-xgboost is refused
  6  an existing environment that matches is accepted, and not modified
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECIPE = ROOT / "conf" / "env" / "install_rdoublet.sh"
# What R reports in an environment built to the pins, and what conda records for it.
GOOD_R = {"R": "4.3.3", "scDblFinder": "1.16.0", "SingleCellExperiment": "1.24.0", "xgboost": "1.7.6.1",
          "Matrix": "1.6.5"}
GOOD_CONDA = [("r-base", "4.3.3", "he2d9a6e_8"), ("r-xgboost", "1.7.6", "cpu_r43h6728c87_0"),
              ("r-matrix", "1.6_5", "r43he966344_1"), ("bioconductor-scdblfinder", "1.16.0", "r43hdfd78af_0"),
              ("bioconductor-singlecellexperiment", "1.24.0", "r43hdfd78af_0")]
SPECS = ("r-base=4.3.3", "r-matrix=1.6_5", "r-xgboost=1.7.6=cpu_*", "bioconductor-scdblfinder=1.16.0",
         "bioconductor-singlecellexperiment=1.24.0")
fails: list[str] = []

# The stand-in Rscript ignores the code it is given and prints the versions written beside it, as the
# recipe's report asks R to print them: one `name<TAB>version` line each.
RSCRIPT = '#!/bin/sh\ncat "$(dirname "$0")/../standin_report.tsv"\n'

MANAGER = r"""#!/bin/bash
# Stand-in for a conda-family manager: records every call. `create -p P <specs>` makes P with a stand-in
# Rscript reporting $STANDIN_R_REPORT and one conda record per spec; nothing is downloaded or installed.
printf 'manager %s\n' "$*" >> "@LOG@"
case "${1:-}" in
    --version) echo "stand-in 0.0"; exit 0 ;;
    create) shift ;;
    *) exit 0 ;;
esac
p=""; specs=()
while [ $# -gt 0 ]; do
    case "$1" in
        -p|--prefix|-c|--channel) [ "$1" = "-p" ] || [ "$1" = "--prefix" ] && p="$2"; shift 2 ;;
        -*) shift ;;
        *) specs+=("$1"); shift ;;
    esac
done
mkdir -p "$p/bin" "$p/conda-meta"
cp "@RSCRIPT@" "$p/bin/Rscript"; chmod +x "$p/bin/Rscript"
printf '%s\n' "$STANDIN_R_REPORT" | tr ',' '\n' | tr '=' '\t' > "$p/standin_report.tsv"
for s in "${specs[@]}"; do
    name="${s%%=*}"; rest="${s#"$name"}"; rest="${rest#=}"
    version="${rest%%=*}"; build="${rest#"$version"}"; build="${build#=}"
    [ -n "$version" ] || version=0
    build="${build:-standin_0}"; build="${build//\*/r43standin_0}"
    echo "{\"name\": \"$name\", \"version\": \"$version\", \"build\": \"$build\"}" \
        > "$p/conda-meta/$name-$version-$build.json"
done
"""


def setup(td: Path):
    stand = td / "standin"
    stand.mkdir()
    log = td / "calls.log"
    (stand / "Rscript.in").write_text(RSCRIPT, encoding="utf-8")
    (stand / "micromamba").write_text(MANAGER.replace("@LOG@", str(log))
                                      .replace("@RSCRIPT@", str(stand / "Rscript.in")), encoding="utf-8")
    for name in ("curl", "wget"):
        (stand / name).write_text(f'#!/bin/sh\necho "{name} $*" >> "{log}"\nexit 22\n', encoding="utf-8")
    for p in stand.iterdir():
        p.chmod(0o755)
    return stand, log


def existing(prefix: Path, report: dict, conda: list):
    """An rdoublet environment as a build would have left it: an Rscript, conda records, what R reports."""
    env = prefix / "rdoublet"
    (env / "bin").mkdir(parents=True)
    (env / "conda-meta").mkdir()
    (env / "bin" / "Rscript").write_text(RSCRIPT, encoding="utf-8")
    (env / "bin" / "Rscript").chmod(0o755)
    (env / "standin_report.tsv").write_text("".join(f"{k}\t{v}\n" for k, v in report.items()), encoding="utf-8")
    for name, version, build in conda:
        (env / "conda-meta" / f"{name}-{version}-{build}.json").write_text(
            json.dumps({"name": name, "version": version, "build": build}), encoding="utf-8")


def run(td: Path, *args: str, report: dict | None = None):
    stand, log = td / "standin", td / "calls.log"
    env = {"PATH": os.pathsep.join([str(stand), "/usr/bin", "/bin"]), "HOME": str(td), "LC_ALL": "C",
           "TMPDIR": str(td), "STANDIN_R_REPORT": ",".join(f"{k}={v}" for k, v in (report or {}).items())}
    r = subprocess.run(["bash", str(RECIPE), *args], capture_output=True, text=True, env=env, cwd=str(td),
                       timeout=120)
    calls = log.read_text(encoding="utf-8").splitlines() if log.exists() else []
    return r.returncode, r.stdout + r.stderr, calls


print("conf/env/install_rdoublet.sh - the R doublet environment, read as R reports it")
print("=" * 74)

# ---- 1. --check
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    setup(td)
    rc, said, calls = run(td, "--check")
    named = ("r-xgboost=1.7.6=cpu_*", "r-matrix=1.6_5", "xgboost 1.7.6.1", "Matrix 1.6.5")
    if rc != 0 or calls or not all(n in said for n in named):
        fails.append(f"1: --check exit {rc}; named {[n for n in named if n in said]} of {list(named)}; "
                     f"calls {calls[:1]}: {said.strip()[-240:]!r}")
    else:
        print("  1  --check names the conda specs and the versions R must report (xgboost 1.7.6.1, Matrix 1.6.5)")

# ---- 2. a fresh build, R reporting what conda-forge's packages report
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    setup(td)
    rc, said, calls = run(td, str(td / "envs"), report=GOOD_R)
    create = next((c for c in calls if c.startswith("manager create")), "")
    missing = [s for s in SPECS if f" {s}" not in f" {create} "]
    if rc != 0 or missing:
        fails.append(f"2: a fresh build whose R reports xgboost 1.7.6.1: exit {rc}; specs not asked for "
                     f"{missing}: {said.strip()[-300:]!r}")
    else:
        print("  2  a fresh build asks for exactly the pinned specs; R's xgboost 1.7.6.1 and Matrix 1.6.5 accepted")

# ---- 3-6. existing environments: checked, never modified
CASES = [
    ("3 xgboost 2.0.0", dict(GOOD_R, xgboost="2.0.0"), GOOD_CONDA, ("xgboost 2.0.0", "shim"), True),
    ("4a scDblFinder 1.14.0", dict(GOOD_R, scDblFinder="1.14.0"), GOOD_CONDA, ("scDblFinder 1.14.0",), True),
    ("4b Matrix 1.6.4", dict(GOOD_R, Matrix="1.6.4"), GOOD_CONDA, ("Matrix 1.6.4",), True),
    ("5 a GPU build of r-xgboost", GOOD_R,
     [c if c[0] != "r-xgboost" else ("r-xgboost", "1.7.6", "cuda120_r43ha4afa4e_0") for c in GOOD_CONDA],
     ("cuda120_r43ha4afa4e_0",), True),
    ("6 an environment that matches", GOOD_R, GOOD_CONDA, (), False),
]
for case, report, conda, named, refuse in CASES:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        setup(td)
        existing(td / "envs", report, conda)
        rc, said, calls = run(td, str(td / "envs"))
        modified = [c for c in calls if c.startswith("manager") and " --version" not in f" {c}"]
        fix = "Remove it and re-run, or build into a new prefix" in " ".join(said.split())
        if refuse and (rc == 0 or not all(n in said for n in named) or not fix or modified):
            fails.append(f"{case}: exit {rc}; named {[n for n in named if n in said]}; fix named {fix}; "
                         f"calls that modify it {modified[:1]}: {said.strip()[-260:]!r}")
        elif not refuse and (rc != 0 or modified):
            fails.append(f"{case}: exit {rc}; calls that modify it {modified[:1]}: {said.strip()[-260:]!r}")
        else:
            print(f"  {case}: {'refused, the fix named' if refuse else 'accepted'} (exit {rc}); not modified")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the R doublet environment is checked as R reports it, and an existing one is never modified")
