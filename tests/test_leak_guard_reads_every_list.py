# Runs the leak guard (tests/test_portability.py) as a subprocess against term lists written to a
# temporary directory; reads no cohort's list and needs only the standard library.
"""The leak guard reads every term list it is given, and refuses one that gives it nothing.

`tests/test_portability.py` read ONE path from $SCQC_FORBIDDEN_TERMS. Handed two lists joined as
PATH joins them - the form the harness's own conform check takes - it found no file of that name,
scanned 0 terms and printed PASS (harness autodev W001's verifier, gap C06). A guard that proves
nothing in the words of one that proved everything is worse than no guard, because it is the one
nobody re-checks. Each case below fails on the guard that read one path:

  1  two lists joined by os.pathsep: a term that is ONLY in the second list and IS in this
     repository must be found - so the second list was read
  2  a named list that holds no term: refused, not a pass over 0 terms
  3  a named list that does not exist: refused, naming it
  4  one ordinary list whose terms are absent from the tree: still a pass (the guard is not
     refusing everything)
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GUARD = ROOT / "tests" / "test_portability.py"
fails: list[str] = []

# A word certain to be in this repository (its own CLI's name) and one certain not to be: a list
# with only the second must pass, and a list carrying the first must be caught.
PRESENT = "selftest"
ABSENT = "zq" + "xv" + "wj" + "absentterm"


def guard(spec: str):
    env = dict(os.environ)
    env["SCQC_FORBIDDEN_TERMS"] = spec
    p = subprocess.run([sys.executable, str(GUARD)], capture_output=True, text=True, env=env,
                       cwd=str(ROOT), timeout=300)
    return p.returncode, p.stdout + p.stderr


with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    first, second = td / "first.txt", td / "second.txt"
    first.write_text(f"# a cohort's list\n{ABSENT}\n", encoding="utf-8")
    second.write_text(f"# another cohort's list\n{PRESENT}\n", encoding="utf-8")
    empty = td / "empty.txt"
    empty.write_text("# a list with a header and no term\n\n", encoding="utf-8")

    rc, out = guard(f"{first}{os.pathsep}{second}")
    if rc == 0 or PRESENT not in out:
        fails.append(f"1: two lists joined by {os.pathsep!r}: exit {rc}; the second list's term "
                     f"{PRESENT!r} was not searched for:\n{out[-400:]}")
    else:
        print(f"  1  both lists read: the second list's term was found (exit {rc})")

    rc, out = guard(str(empty))
    if rc == 0 or "0 terms" not in out:
        fails.append(f"2: a list with no term: exit {rc}; it must refuse, naming 0 terms:\n"
                     f"{out[-400:]}")
    else:
        print("  2  a named list with no term is refused")

    missing = td / "no_such_list.txt"
    rc, out = guard(f"{first}{os.pathsep}{missing}")
    if rc == 0 or "no such file" not in out:
        fails.append(f"3: a list that does not exist: exit {rc}; it must refuse, naming it:\n"
                     f"{out[-400:]}")
    else:
        print("  3  a named list that does not exist is refused")

    rc, out = guard(str(first))
    if rc != 0:
        fails.append(f"4: one list whose term is absent from the tree: exit {rc}; it must pass:\n"
                     f"{out[-400:]}")
    else:
        print("  4  one ordinary list, nothing found: PASS")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the leak guard reads every list it is given and refuses one that gives it nothing")
