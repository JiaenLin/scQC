# Copies the leak guard (tests/test_portability.py) into a temporary tree of made-up files and runs
# it there as a subprocess, against made-up lists; reads no cohort's list and needs only the stdlib.
"""The leak guard reads a generated word list the way the harness's `sch conform` reads it.

The harness generates each cohort's list from the cohort's records (`sch conform terms build`), and
nearly every line of it carries a prefix that says how it is matched: `word:T` (a whole word, any
case, a number never inside a longer number), `exact:T` (the same, in T's case), `T*` (a word that
begins with T). The site's shapes come as a second file of `<regex> :: <what it is> :: <literal>`.
Until harness autodev W048 this guard read every line as literal text and read no shapes file, so
the generated lists could match almost nothing: on one tree it printed PASS where the harness's own
scan, given the same lists, found every hit (harness GAP-C78). Each case below fails on that guard:

  1  `word:T` finds T as a whole word in any case - in a file with no extension too, which the
     old guard never opened - and skips (and counts) a binary file
  2  `word:T` does not find T inside a longer word
  3  `word:T*` finds a word that begins with T
  4  `exact:T` finds T in its own case and not in another
  5  `word:N` for a number does not find it inside a longer number, and does find it alone
  6  spaces in `word:T` match any run of whitespace
  7  a line with no prefix is still the text anywhere, in any case
  8  a site shape from $SCQC_SITE_SHAPES is found; a malformed shape line, or a shapes file that
     does not exist, fails the guard rather than passing it
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GUARD = ROOT / "tests" / "test_portability.py"
fails: list[str] = []

# Made-up words, built so this file never spells them whole.
W = "zor" + "vat"
LOW_EXACT, UP_EXACT = "qe" + "lp", ("qe" + "lp").upper()
NUM = "7," + "654"
HOST = "vq" + "1234" + "x"
NOTES = "\n".join([
    f"the arm {W.capitalize()} and the arm {W}ine",          # 1
    f"{UP_EXACT} is the label",                               # 2
    f"{LOW_EXACT} is a variable",                             # 3
    f"1{NUM} and {NUM}0 and {NUM}.5 cells",                  # 4
    f"{NUM} cells",                                           # 5
    "alpha    beta",                                          # 6
    f"{HOST} ran it",                                         # 7
]) + "\n"


def tree(td: Path) -> Path:
    """A tree with the guard at tests/, so its ROOT is the tree, and three made-up files."""
    root = td / "repo"
    (root / "tests").mkdir(parents=True)
    shutil.copyfile(GUARD, root / "tests" / "test_portability.py")
    (root / "notes.txt").write_text(NOTES, encoding="utf-8")
    (root / "bin").mkdir()
    (root / "bin" / "tool").write_text(f"#!/bin/sh\necho {W}\n", encoding="utf-8")
    (root / "blob.bin").write_bytes(b"\x00\x01" + W.encode() + b"\x00")
    return root


def run(root: Path, td: Path, words: list, shapes: list | None = None):
    lst = td / "terms.txt"
    lst.write_text("# made-up terms\n" + "".join(w + "\n" for w in words), encoding="utf-8")
    env = {k: v for k, v in os.environ.items() if not k.startswith("SCQC_")}
    env["SCQC_FORBIDDEN_TERMS"] = str(lst)
    if shapes is not None:
        sh = td / "shapes.txt"
        sh.write_text("# made-up shapes\n" + "".join(s + "\n" for s in shapes), encoding="utf-8")
        env["SCQC_SITE_SHAPES"] = str(sh)
    p = subprocess.run([sys.executable, str(root / "tests" / "test_portability.py")],
                       capture_output=True, text=True, env=env, cwd=str(root), timeout=120)
    out = p.stdout + p.stderr
    found = set(re.findall(r"^  LEAK (\S+?):(\d+): ", out, re.M))
    return p.returncode, out, {(f, int(n)) for f, n in found}


def expect(case: str, got, want_rc: int, want: set, also: str = ""):
    rc, out, found = got
    if rc != want_rc or found != want or (also and also not in out):
        fails.append(f"{case}: exit {rc} (want {want_rc}), found {sorted(found)} (want "
                     f"{sorted(want)}){'; output lacks ' + repr(also) if also and also not in out else ''}"
                     f":\n{out[-500:]}")
    else:
        print(f"  {case}: as the harness reads it")


with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    root = tree(td)
    NOTE = "notes.txt"
    expect("1 word:, every file by content, binary skipped",
           run(root, td, [f"word:{W}"]), 1, {(NOTE, 1), ("bin/tool", 2)}, "1 binary skipped")
    expect("2 word: is not inside a longer word", run(root, td, [f"word:{W[:-1]}"]), 0, set())
    expect("3 word:T* is a prefix", run(root, td, [f"word:{W[:-1]}*"]), 1,
           {(NOTE, 1), ("bin/tool", 2)})
    expect("4 exact: keeps its case", run(root, td, [f"exact:{UP_EXACT}"]), 1, {(NOTE, 2)})
    expect("5 a number is never inside a longer one", run(root, td, [f"word:{NUM}"]), 1,
           {(NOTE, 5)})
    expect("6 spaces match a run of whitespace", run(root, td, ["word:alpha beta"]), 1, {(NOTE, 6)})
    expect("7 no prefix is the text anywhere", run(root, td, [W[1:-1]]), 1,
           {(NOTE, 1), ("bin/tool", 2)})
    absent = "zq" + "xv" + "wj" + "absentterm"
    expect("8a a site shape is read", run(root, td, [absent], [r"\bvq\d{4}x\b :: a made-up host :: vq"]),
           1, {(NOTE, 7)}, "a made-up host")
    expect("8b a malformed shape fails the guard", run(root, td, [absent], [r"(unclosed :: a broken shape"]),
           1, set(), "not a regular expression")
    rc, out, found = run(root, td, [absent])
    env_missing = dict(os.environ, SCQC_FORBIDDEN_TERMS=str(td / "terms.txt"),
                       SCQC_SITE_SHAPES=str(td / "no_such_shapes.txt"))
    p = subprocess.run([sys.executable, str(root / "tests" / "test_portability.py")],
                       capture_output=True, text=True, env=env_missing, cwd=str(root), timeout=120)
    if p.returncode == 0 or "no such file" not in p.stdout + p.stderr:
        fails.append(f"8c: a shapes file that does not exist: exit {p.returncode}; it must refuse, "
                     f"naming it:\n{(p.stdout + p.stderr)[-400:]}")
    else:
        print("  8c a named shapes file that does not exist is refused")
    if rc != 0:
        fails.append(f"control: a list whose only term is absent must pass, exit {rc}:\n{out[-400:]}")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the leak guard reads the generated lists' syntax and the site's shapes as sch conform does")
