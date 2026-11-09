# Runs setup/install_env.sh and conf/env/install_rdoublet.sh under bash in temporary directories, against a
# stand-in package manager, pip and interpreter written there; installs nothing, downloads nothing.
"""An interrupted build is repaired with --repair, not by hand (harness autodev W020, B1-108).

A build killed mid-transaction leaves an environment directory holding conda-meta/ and no interpreter.
The installer refuses it - correctly, an installer does not delete what it did not verify it made - and
the refusal told every reader to type `rm -rf` before the environment could build again. The refusal
stays; `--repair` is the other half: the one shape an interrupted build leaves is removed, rebuilt and
printed.

  1  core half-built (conda-meta/, no bin/python), no --repair: refused, the refusal names --repair, and
     nothing is removed
  2  the same with --repair: removed, rebuilt from the lock, printed as a repair, and core verifies ok
  3  a directory with no conda-meta/ is not known to be a build: --repair refuses it and removes nothing
  4  the R doublet environment half-built (conda-meta/, no bin/Rscript): refused without SCQC_REPAIR,
     naming --repair; with it, removed and created again

The stand-ins are the ones tests/test_the_denoiser_environment_is_pinned.py uses, copied: a package
manager that makes an environment with a stand-in interpreter, and a pip that writes a record per pin.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "setup" / "install_env.sh"
RDOUBLET = ROOT / "conf" / "env" / "install_rdoublet.sh"
REAL = sys.executable
fails: list[str] = []

MANAGER = r"""#!/bin/bash
# Stand-in for a conda-family package manager: records every call; `create -p P ... python=V` makes P with
# a stand-in interpreter and pip. Nothing is downloaded or installed.
printf 'manager %s\n' "$*" >> "@LOG@"
case "${1:-}" in
    --version) echo "stand-in 0.0"; exit 0 ;;
    create) shift ;;
    *) exit 0 ;;
esac
p=""; v=""
while [ $# -gt 0 ]; do
    case "$1" in
        -p|--prefix) p="$2"; shift 2 ;;
        python=*) v="${1#python=}"; shift ;;
        *) shift ;;
    esac
done
[ -n "$p" ] || { echo "stand-in: create without -p" >&2; exit 2; }
mkdir -p "$p/bin" "$p/site"
cat > "$p/bin/python" <<EOF
#!/bin/bash
case "\$*" in *version_info*) echo "$v"; exit 0 ;; esac
[ "\${1:-}" = "--version" ] && { echo "Python $v"; exit 0; }
if [ "\${1:-}" = "-m" ] && [ "\${2:-}" = "pip" ]; then shift 2; exec "$p/bin/pip" "\$@"; fi
if [ "\${1:-}" = "-c" ]; then case "\${2:-}" in import\ *) exit 0 ;; esac; fi
PYTHONPATH="$p/site" exec "@REAL@" -S "\$@"
EOF
cat > "$p/bin/pip" <<EOF
#!/bin/bash
exec "@REAL@" -S "@PIP@" "$p" "\$@"
EOF
chmod +x "$p/bin/python" "$p/bin/pip"
"""

PIP = r'''"""Stand-in pip: records the call; `install` writes one distribution record per requirement given."""
import json, os, re, sys

env, argv = sys.argv[1], sys.argv[2:]
with open(os.environ["STANDIN_LOG"], "a", encoding="utf-8") as fh:
    fh.write("pip " + json.dumps({"env": env, "argv": argv}) + "\n")
if not argv or argv[0] != "install":
    sys.exit(0)
reqs, i = [], 1
while i < len(argv):
    a = argv[i]
    if a in ("-r", "--requirement"):
        with open(argv[i + 1], encoding="utf-8") as fh:
            reqs += [ln.strip() for ln in fh if ln.strip() and not ln.strip().startswith(("#", "-"))]
        i += 2
    elif a in ("-i", "--index-url", "--extra-index-url"):
        i += 2
    elif a.startswith("-"):
        i += 1
    else:
        reqs.append(a)
        i += 1
drift = dict(x.split("=", 1) for x in os.environ.get("STANDIN_DRIFT", "").split(",") if "=" in x)
site = os.path.join(env, "site")
for r in reqs:
    spec = r.split(";")[0].strip()
    m = re.match(r"^([A-Za-z0-9][A-Za-z0-9._-]*)\s*(?:==\s*(\S+))?", spec)
    name, version = m.group(1), m.group(2) or "0+any"
    key = re.sub(r"[-_.]+", "_", name).lower()
    version = drift.get(key, version)
    for old in os.listdir(site):
        if old.lower().startswith(key + "-") and old.endswith(".dist-info"):
            for f in os.listdir(os.path.join(site, old)):
                os.remove(os.path.join(site, old, f))
            os.rmdir(os.path.join(site, old))
    d = os.path.join(site, f"{key}-{version}.dist-info")
    os.makedirs(d)
    with open(os.path.join(d, "METADATA"), "w", encoding="utf-8") as fh:
        fh.write(f"Metadata-Version: 2.1\nName: {name}\nVersion: {version}\n")
'''



def stand_ins(td: Path) -> tuple[Path, Path]:
    """A directory of stand-ins (manager, and the downloaders the installer could reach for), and its log."""
    bindir, log = td / "standin", td / "standin.log"
    bindir.mkdir()
    (bindir / "standin_pip.py").write_text(PIP, encoding="utf-8")
    manager = MANAGER.replace("@LOG@", str(log)).replace("@REAL@", REAL)
    manager = manager.replace("@PIP@", str(bindir / "standin_pip.py"))
    for name in ("micromamba", "mamba", "conda"):
        (bindir / name).write_text(manager, encoding="utf-8")
    for name in ("curl", "wget"):
        (bindir / name).write_text(f'#!/bin/sh\necho "{name} $*" >> "{log}"\nexit 22\n', encoding="utf-8")
    for p in bindir.iterdir():
        p.chmod(0o755)
    return bindir, log


def install(td: Path, installer: Path, *args: str, drift: str = ""):
    bindir, log = td / "standin", td / "standin.log"
    env = {"PATH": os.pathsep.join([str(bindir), "/usr/bin", "/bin"]), "HOME": str(td), "LC_ALL": "C",
           "TMPDIR": str(td), "STANDIN_LOG": str(log), "STANDIN_DRIFT": drift}
    r = subprocess.run(["bash", str(installer), *args], capture_output=True, text=True, env=env,
                       cwd=str(td), timeout=600)
    return r.returncode, r.stdout + r.stderr


def calls(log: Path) -> list[str]:
    return log.read_text(encoding="utf-8").splitlines() if log.exists() else []




def half_built(prefix: Path, name: str, conda_meta: bool = True) -> Path:
    """What a killed build leaves: the directory, its conda-meta/ (or not), and no interpreter."""
    d = prefix / name
    (d / "conda-meta").mkdir(parents=True) if conda_meta else d.mkdir(parents=True)
    (d / "left_by_the_killed_build.txt").write_text("x", encoding="utf-8")
    return d


def run_rdoublet(td: Path, prefix: Path, repair: bool):
    bindir, log = td / "standin", td / "standin.log"
    env = {"PATH": os.pathsep.join([str(bindir), "/usr/bin", "/bin"]), "HOME": str(td), "LC_ALL": "C",
           "TMPDIR": str(td), "STANDIN_LOG": str(log), "SCQC_CONDA": str(bindir / "micromamba"),
           **({"SCQC_REPAIR": "1"} if repair else {})}
    r = subprocess.run(["bash", str(RDOUBLET), str(prefix)], capture_output=True, text=True, env=env,
                       cwd=str(td), timeout=600)
    return r.returncode, r.stdout + r.stderr


# ---- 1-2. core, refused without --repair and rebuilt with it
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    stand_ins(td)
    prefix = td / "envs"
    core = half_built(prefix, "core")
    rc, out = install(td, INSTALLER, "--prefix", str(prefix))
    if rc == 0:
        fails.append("1: a half-built core environment was not refused")
    if "--repair" not in out:
        fails.append(f"1: the refusal does not name --repair: {out[-300:]!r}")
    if not (core / "left_by_the_killed_build.txt").exists():
        fails.append("1: the half-built environment was removed without --repair")
    print("1. a half-built core is refused without --repair, naming it, and left as it was")

    rc, out = install(td, INSTALLER, "--prefix", str(prefix), "--repair")
    if rc != 0:
        fails.append(f"2: --repair did not complete (exit {rc}): {out[-400:]!r}")
    if "[repair]" not in out:
        fails.append("2: the repair was not printed")
    if (core / "left_by_the_killed_build.txt").exists():
        fails.append("2: --repair did not remove what the interrupted build left")
    if not (core / "bin" / "python").exists():
        fails.append("2: --repair did not build core again")
    if not any(ln.strip().startswith("core") and " ok " in ln + " " for ln in out.splitlines()):
        fails.append(f"2: core did not verify ok after the repair: {out[-400:]!r}")
    print("2. --repair removes it, builds core again from the lock, and core verifies ok")

# ---- 3. a directory that is not known to be a build is never removed
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    stand_ins(td)
    prefix = td / "envs"
    other = half_built(prefix, "core", conda_meta=False)
    rc, out = install(td, INSTALLER, "--prefix", str(prefix), "--repair")
    if rc == 0:
        fails.append("3: a directory with no conda-meta/ was accepted")
    if not (other / "left_by_the_killed_build.txt").exists():
        fails.append("3: --repair removed a directory it could not identify as a build")
    print("3. a directory without conda-meta/ is refused by --repair and not removed")

# ---- 4. the R doublet environment
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    _bindir, log = stand_ins(td)
    prefix = td / "envs"
    rd = half_built(prefix, "rdoublet")
    rc, out = run_rdoublet(td, prefix, repair=False)
    if rc == 0 or "--repair" not in out:
        fails.append(f"4: a half-built rdoublet was not refused naming --repair (exit {rc}): {out[-300:]!r}")
    if not (rd / "left_by_the_killed_build.txt").exists():
        fails.append("4: the half-built rdoublet was removed without --repair")
    rc, out = run_rdoublet(td, prefix, repair=True)
    if "[repair]" not in out or (rd / "left_by_the_killed_build.txt").exists():
        fails.append(f"4: SCQC_REPAIR=1 did not remove the half-built rdoublet: {out[-300:]!r}")
    if not any(ln.startswith("manager create") and str(rd) in ln for ln in calls(log)):
        fails.append("4: after removing it, rdoublet was not created again")
    print("4. a half-built rdoublet is refused naming --repair; with it, removed and created again")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("an interrupted build is repaired by --repair, and only an interrupted build")
