# Runs setup/install_env.sh under bash in temporary directories, against a stand-in package manager, pip and
# interpreter written there; installs nothing, downloads nothing, and needs only the standard library.
"""The denoiser's environment is built from a lock that pins torch and numpy, and no other is accepted.

The installer built the CellBender environment with `pip install cellbender==0.3.2` and named no torch, so
pip took the newest: a cohort's first build got torch 2.14.1, and every denoiser task then died inside
CellBender's checkpoint (harness ADR-0027, finding Q9). `conf/env/install_cs_cb.sh` installed torch with no
version at all. The environment is now one lock, applied without resolving, and three versions are required
of it whatever the lock says: cellbender 0.3.2, torch 2.5.1, numpy 1.26.4.

  1  the shipped lock: every line exact, its interpreter exact, and the three versions in it
  2  `install_env.sh --check` accepts the shipped locks, names the R doublet pins, and builds nothing
  3  a lock naming torch unpinned, as a range, at another version, or not at all - or another numpy - is
     refused before any package manager is called, and nothing is created
  4  a build against the stand-ins applies the core and denoiser locks with `--no-deps`, and the denoiser
     environment holds torch 2.5.1 and numpy 1.26.4
  5  a pip that installs another torch than the lock's is caught: the installer exits non-zero naming it
  6  an existing denoiser environment holding another torch is refused, and not modified

The stand-in manager makes an environment whose interpreter answers the installer's version query and runs
anything else in this interpreter over that environment's own site directory; the stand-in pip records its
call and writes a distribution record for each requirement it is given. What the installer checks is
therefore read off records the installer did not write.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "setup" / "install_env.sh"
LOCK = ROOT / "conf" / "env" / "cellbender.lock.txt"
REQUIRED = {"cellbender": "0.3.2", "torch": "2.5.1", "numpy": "1.26.4"}
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


def installed(env: Path) -> dict:
    """{normalised name: version} from the distribution records in a stand-in environment."""
    out = {}
    for d in (env / "site").glob("*.dist-info"):
        meta = (d / "METADATA").read_text(encoding="utf-8")
        name = re.search(r"^Name: (.+)$", meta, re.M).group(1)
        out[re.sub(r"[-_.]+", "-", name).lower()] = re.search(r"^Version: (.+)$", meta, re.M).group(1)
    return out


print("The denoiser environment: one lock, applied without resolving, torch and numpy required")
print("=" * 74)

# ---- 1. the shipped lock
if not LOCK.is_file():
    fails.append(f"1: no denoiser lock at {LOCK.relative_to(ROOT)}: nothing pins torch or numpy")
else:
    text = LOCK.read_text(encoding="utf-8")
    py = re.search(r"^#\s*scqc-python:\s*([0-9][0-9.]*)", text, re.M)
    pins, loose = {}, []
    for n, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#") or re.match(r"^--(extra-)?index-url[\s=]+\S+$", line):
            continue
        m = re.match(r"^([A-Za-z0-9][A-Za-z0-9._-]*)==([A-Za-z0-9][A-Za-z0-9.+!_-]*)\s*(;.*)?$", line)
        if not m:
            loose.append(f"line {n}: {line}")
            continue
        pins[re.sub(r"[-_.]+", "-", m.group(1)).lower()] = m.group(2)
    wrong = {k: pins.get(k) for k, v in REQUIRED.items() if (pins.get(k) or "").split("+")[0] != v}
    if loose or wrong or not py or py.group(1).count(".") != 2:
        fails.append(f"1: the shipped lock: interpreter {py.group(1) if py else None}; not exact {loose[:3]}; "
                     f"required versions not pinned {wrong}")
    else:
        print(f"  1  {len(pins)} exact pins, python {py.group(1)}; torch {pins['torch']}, numpy "
              f"{pins['numpy']}, cellbender {pins['cellbender']}")

# ---- 2. --check: the shipped locks accepted, nothing built
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    _, log = stand_ins(td)
    rc, said = install(td, INSTALLER, "--check", "--with-cellbender", "--with-doublet")
    rpins = ("r-base=4.3.3", "bioconductor-scdblfinder=1.16.0", "r-xgboost=1.7.6")
    if rc != 0 or calls(log) or not all(p in said for p in rpins):
        fails.append(f"2: --check exit {rc}; calls {calls(log)[:2]}; R pins named "
                     f"{[p for p in rpins if p in said]}: {said.strip()[-300:]!r}")
    else:
        print(f"  2  --check accepts the shipped locks, names {', '.join(rpins)}, and calls nothing")

# ---- 3. refusals, before anything is looked for or built
MUTATIONS = [
    ("torch unpinned", "torch", lambda name, line: "torch" if name == "torch" else line, "2.5.1"),
    ("torch as a range", "torch", lambda name, line: "torch>=2.5" if name == "torch" else line, "2.5.1"),
    ("torch 2.14.1", "torch", lambda name, line: "torch==2.14.1" if name == "torch" else line, "2.14.1"),
    ("torch not named", "torch", lambda name, line: None if name == "torch" else line, "2.5.1"),
    ("numpy 2.1.3", "numpy", lambda name, line: "numpy==2.1.3" if name == "numpy" else line, "1.26.4"),
]
for case, pkg, edit, must in MUTATIONS:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        _, log = stand_ins(td)
        repo = td / "repo"
        shutil.copytree(ROOT / "setup", repo / "setup")
        shutil.copytree(ROOT / "conf" / "env", repo / "conf" / "env")
        lock = repo / "conf" / "env" / "cellbender.lock.txt"
        if not lock.is_file():
            fails.append(f"3 {case}: no denoiser lock to edit - nothing pins {pkg}, so nothing can refuse it")
            continue
        out = []
        for line in lock.read_text(encoding="utf-8").splitlines():
            m = re.match(r"^([A-Za-z0-9][A-Za-z0-9._-]*)", line.strip())
            name = re.sub(r"[-_.]+", "-", m.group(1)).lower() if m else ""
            new = edit(name, line)
            if new is not None:
                out.append(new)
        lock.write_text("\n".join(out) + "\n", encoding="utf-8")
        prefix = td / "envs"
        rc, said = install(td, repo / "setup" / "install_env.sh", "--prefix", str(prefix), "--with-cellbender")
        if rc == 0 or pkg not in said or must not in said or calls(log) or prefix.exists():
            fails.append(f"3 {case}: exit {rc}; '{pkg}' named {pkg in said}, '{must}' named {must in said}; "
                         f"calls before refusing {calls(log)[:2]}; prefix created {prefix.exists()}: "
                         f"{said.strip()[-240:]!r}")
        else:
            print(f"  3  {case}: refused (exit {rc}) before any package manager was called; nothing created")

# ---- 4-6. a build against the stand-ins
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    _, log = stand_ins(td)
    prefix = td / "envs"
    rc, said = install(td, INSTALLER, "--prefix", str(prefix), "--with-cellbender")
    pips = [json.loads(c[4:]) for c in calls(log) if c.startswith("pip ")]
    by_env = {Path(p["env"]).name: p["argv"] for p in pips if p["argv"][:1] == ["install"]}
    held = installed(prefix / "cellbender") if (prefix / "cellbender" / "site").is_dir() else {}
    torch_ok = held.get("torch", "").split("+")[0] == REQUIRED["torch"]
    if (rc != 0 or set(by_env) != {"core", "cellbender"}
            or not all("--no-deps" in a for a in by_env.values())
            or not torch_ok or held.get("numpy") != REQUIRED["numpy"]):
        fails.append(f"4: exit {rc}; pip installs {by_env}; the denoiser environment holds torch "
                     f"{held.get('torch')}, numpy {held.get('numpy')}: {said.strip()[-300:]!r}")
    else:
        print(f"  4  core and denoiser locks applied with --no-deps; the denoiser holds torch "
              f"{held['torch']}, numpy {held['numpy']} ({len(held)} distributions)")

    # 6: the environment now exists; make it hold another torch, and run the installer again
    before = len(calls(log))
    rec = next((prefix / "cellbender" / "site").glob("torch-*.dist-info"), None)
    if rec is None:
        fails.append("6: the stand-in build left no torch in the denoiser environment to change")
    else:
        meta = rec / "METADATA"
        meta.write_text(re.sub(r"^Version: .*$", "Version: 2.14.1", meta.read_text(encoding="utf-8"),
                               flags=re.M), encoding="utf-8")
        rc, said = install(td, INSTALLER, "--prefix", str(prefix), "--with-cellbender")
        again = calls(log)[before:]
        touched = [c for c in again if "cellbender" in c and (c.startswith("pip ") or " create " in c)]
        if rc == 0 or "2.14.1" not in said or "torch" not in said or touched:
            fails.append(f"6: an existing environment holding torch 2.14.1: exit {rc}; named "
                         f"{'2.14.1' in said}; calls that modify it {touched[:2]}: {said.strip()[-240:]!r}")
        else:
            print(f"  6  an existing environment holding torch 2.14.1: refused (exit {rc}), not modified")

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    stand_ins(td)
    rc, said = install(td, INSTALLER, "--prefix", str(td / "envs"), "--with-cellbender",
                       drift="torch=2.14.1")
    if rc == 0 or "2.14.1" not in said or "torch" not in said:
        fails.append(f"5: pip installed torch 2.14.1 for the lock's 2.5.1: exit {rc}, named {'2.14.1' in said}: "
                     f"{said.strip()[-240:]!r}")
    else:
        print(f"  5  a pip that installed torch 2.14.1 for the lock's 2.5.1 is caught (exit {rc})")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the denoiser environment is pinned, applied without resolving, and nothing else is accepted")
