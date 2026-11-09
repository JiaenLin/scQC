# Reads conf/env/requirements.lock.txt and the import statements of this repository's Python files; runs
# nothing, and needs only the standard library.
"""The analysis environment's lock carries what this repository's own code and suite need.

`setup/install_env.sh` builds `core` from `conf/env/requirements.lock.txt` and from nothing else. That lock
held no pyarrow, so the environment the installer built failed this repository's own suite:
`tests/test_native_io.py` sets pandas' string storage to "pyarrow", which pandas refuses where pyarrow is not
installed (harness autodev W001's verifier, gap C05). The cluster's lock beside it held pyarrow; no script
read it.

  1  every requirement is one package at one exact version, and the interpreter line names
     major.minor.patch: the installer applies the lock without resolving and compares the interpreter's
     full version, which every analysis call record writes down (`versions.python`)
  2  pyarrow is pinned
  3  every third-party module imported anywhere in this repository is provided by a pinned distribution
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "conf" / "env" / "requirements.lock.txt"
EXACT = re.compile(r"^([A-Za-z0-9][A-Za-z0-9._-]*)\s*==\s*([A-Za-z0-9][A-Za-z0-9.+!_-]*)\s*(;.*)?$")
INDEX = re.compile(r"^--(extra-)?index-url[\s=]+\S+$")
# An import name that is not its distribution's name.
DISTRIBUTION = {"yaml": "pyyaml", "sklearn": "scikit-learn", "PIL": "pillow", "skimage": "scikit-image"}
# Imported, and deliberately absent from the lock, with the reason.
NOT_LOCKED = {"pytest": "imported only once a pytest run has already imported it "
                        "(tests/test_audit_ambient.py); the suite runs as scripts"}
fails: list[str] = []


def norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


print("conf/env/requirements.lock.txt - the analysis environment's lock")
print("=" * 74)
if not LOCK.is_file():
    raise SystemExit(f"FAILED - no lock at {LOCK}")
text = LOCK.read_text(encoding="utf-8")

# ---- 1. exact, interpreter included
py = re.search(r"^#\s*scqc-python:\s*([0-9][0-9.]*)", text, re.M)
pins: dict[str, str] = {}
loose = []
for n, raw in enumerate(text.splitlines(), 1):
    line = raw.strip()
    if not line or line.startswith("#") or INDEX.match(line):
        continue
    m = EXACT.match(line)
    if not m:
        loose.append(f"line {n}: {line}")
        continue
    pins[norm(m.group(1))] = m.group(2)
if not py or py.group(1).count(".") != 2:
    fails.append(f"1: the interpreter line names {py.group(1) if py else 'nothing'}, not major.minor.patch")
if loose:
    fails.append(f"1: {len(loose)} line(s) are not one package at one exact version: {loose[:3]}")
if py and not loose and py.group(1).count(".") == 2:
    print(f"  1  {len(pins)} exact pins; interpreter {py.group(1)}")

# ---- 2. pyarrow
if "pyarrow" not in pins:
    fails.append("2: pyarrow is not pinned, so tests/test_native_io.py fails in the environment the "
                 "installer builds from this lock")
else:
    print(f"  2  pyarrow {pins['pyarrow']}")

# ---- 3. every third-party import is pinned
if not hasattr(sys, "stdlib_module_names"):
    print("  3  NOT CHECKED: this interpreter has no sys.stdlib_module_names (Python < 3.10)")
else:
    local = {p.stem for p in ROOT.rglob("*.py")} | {p.name for p in ROOT.rglob("*") if p.is_dir()}
    seen: dict[str, str] = {}
    for p in sorted(ROOT.rglob("*.py")):
        if ".git" in p.parts:
            continue
        try:
            tree = ast.parse(p.read_text(encoding="utf-8"))
        except SyntaxError as e:
            fails.append(f"3: {p.relative_to(ROOT)} does not parse: {e}")
            continue
        for node in ast.walk(tree):
            names = ([a.name for a in node.names] if isinstance(node, ast.Import) else
                     [node.module] if isinstance(node, ast.ImportFrom) and node.module and not node.level
                     else [])
            for name in names:
                top = name.split(".")[0]
                if top not in sys.stdlib_module_names and top not in local:
                    seen.setdefault(top, str(p.relative_to(ROOT)))
    unpinned = sorted(m for m in seen if m not in NOT_LOCKED
                      and norm(DISTRIBUTION.get(m, m)) not in pins)
    if unpinned:
        fails.append("3: imported and not pinned: "
                     + ", ".join(f"{m} ({seen[m]})" for m in unpinned))
    else:
        print(f"  3  every third-party module imported is pinned: "
              f"{', '.join(sorted(m for m in seen if m not in NOT_LOCKED))}"
              f"; not locked, by design: {', '.join(sorted(m for m in seen if m in NOT_LOCKED))}")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the core lock is exact and carries what the code and the suite import")
