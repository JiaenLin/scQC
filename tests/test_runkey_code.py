# The pipeline's own version is part of the run key. Without it, a code change republishes.
"""Does a change to the CODE send a run to a different directory?

WHY THIS TEST EXISTS

engine/runkey.py opens by promising that "nothing is ever overwritten by a run that would have
produced something different". The key covered the samplesheet, the declared parameters and the
mode - everything the operator supplies - and not the code.

A DERIVED threshold is a function of the inputs AND of the code that derives it. So the promise
had a hole in exactly the case that matters most, and the 2026-08-13 mitochondrial change walked
straight into it: same samplesheet, same parameters, same mode, a deliverable a few per cent smaller. A
re-run would have resolved to the previous directory, found every task complete, skipped all of
them, and republished the old numbers under the new code - exit zero, report renders, nothing
anywhere says the deliverable is not what the current code produces.

WHAT IS CHECKED

  1. same inputs + same code  -> same digest        (resume still works)
  2. same inputs + DIFFERENT code -> different digest  (the hole, closed)
  3. an unidentifiable checkout does not key as an identified one
  4. a dirty tree does not key as a clean one, and unchecked is neither
  5. the description carries the code, so a reader can see WHY two runs diverged
  6. the index gains a `code` column and an older index is migrated, not appended to crookedly
  7. THE CODE IS KEYED BY ITS CONTENT, NOT ITS COMMIT (harness autodev W020, B1-106): in a copy of
     this tree, a docs-only or tests-only edit keeps the key - a commit that changed only a
     docstring no longer costs every project its completed work - while an edit to a module that
     runs moves it, committed or not; the commit and dirty flag are kept beside the key in
     INPUTS.json, one record per run that used the directory
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import runkey  # noqa: E402

fails: list[str] = []
print("Run key - the pipeline's version is part of it")
print("=" * 74)


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'ok    ' if ok else 'FAILED'} {name}{('   ' + detail) if detail else ''}")
    if not ok:
        fails.append(f"{name}: {detail}")


ROWS = [{"sample": "A", "matrix": "/data/A", "assay": "snrna"},
        {"sample": "B", "matrix": "/data/B", "assay": "snrna"}]
TOOLS = {"light_floor": 200, "dbr": 0.06, "seed": 0}
CLEAN_A = {"commit": "a" * 40, "dirty": False}
CLEAN_B = {"commit": "b" * 40, "dirty": False}


def key(code, rows=ROWS, tools=TOOLS, mode="apply"):
    return runkey.compute(samplesheet_rows=rows, tools=tools, mode=mode, code=code)[0]


# 1 - resume must still work. If this fails the change has cost every project its cached work.
check("same inputs and same code give the same digest", key(CLEAN_A) == key(CLEAN_A),
      "resume depends on this")

# 2 - the hole. This is the whole point of the file.
check("a different commit gives a different digest", key(CLEAN_A) != key(CLEAN_B),
      f"{key(CLEAN_A)} vs {key(CLEAN_B)}")

# 3 - an unidentifiable checkout is its own state, not a free pass into an identified one.
check("an unidentified checkout keys apart from an identified one",
      key({"commit": "unidentified", "dirty": None}) != key(CLEAN_A))

# 4 - dirty, clean and unchecked are three states, not two. `dirty=None` must not read as False.
check("a dirty tree keys apart from a clean one at the same commit",
      key({"commit": "a" * 40, "dirty": True}) != key(CLEAN_A))
check("an unchecked tree keys apart from a verified-clean one",
      key({"commit": "a" * 40, "dirty": None}) != key(CLEAN_A),
      "unknown is not a value")

# ...and the parameters still matter, or the code component has swallowed them.
check("a parameter change still changes the digest",
      key(CLEAN_A, tools={**TOOLS, "light_floor": 500}) != key(CLEAN_A))
check("a samplesheet change still changes the digest",
      key(CLEAN_A, rows=[{**ROWS[0], "matrix": "/data/OTHER"}, ROWS[1]]) != key(CLEAN_A))
check("mode still changes the digest", key(CLEAN_A, mode="evidence") != key(CLEAN_A))

# 5 - explainable, not just different.
_, described = runkey.compute(samplesheet_rows=ROWS, tools=TOOLS, mode="apply", code=CLEAN_A)
check("the description carries the code identity",
      described.get("code") == CLEAN_A, f"{described.get('code')}")

# The real accessors must return something usable on this checkout, or the default is a fiction.
import json  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402

ident = runkey.code_identity()
check("code_identity() is a digest of this checkout's code, and holds no commit",
      isinstance(ident, dict) and re.fullmatch(r"[0-9a-f]{64}", str(ident.get("digest") or ""))
      is not None and "commit" not in ident, f"{ident}")
rec = getattr(runkey, "code_record", lambda: {})()
check("code_record() reads this checkout's commit, for the manifest",
      isinstance(rec.get("commit"), str) and rec["commit"] != "", f"{rec}")

# 7 - content, not commit. A copy of the tree, edited one class of file at a time.
with tempfile.TemporaryDirectory() as td:
    tree = Path(td) / "tree"
    shutil.copytree(ROOT, tree, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.pyc"))

    def ident_of():
        return runkey.code_identity(tree)

    first = ident_of()
    check("the same tree keys the same twice", ident_of() == first, f"{first}")
    (tree / "docs" / "FILTERS.md").write_text(
        (tree / "docs" / "FILTERS.md").read_text(encoding="utf-8") + "\nan edited sentence\n",
        encoding="utf-8")
    (tree / "README.md").write_text("a README rewritten\n", encoding="utf-8")
    check("a docs-only edit keeps the key", ident_of() == first,
          "a commit that changed only a document must not cost a project its completed work")
    (tree / "tests" / "test_runkey_code.py").write_text("# a test rewritten\n", encoding="utf-8")
    check("a tests-only edit keeps the key", ident_of() == first)
    (tree / "engine" / "__pycache__").mkdir(exist_ok=True)
    (tree / "engine" / "__pycache__" / "x.cpython-312.pyc").write_bytes(b"\x00compiled")
    check("a compiled cache beside the code keeps the key", ident_of() == first)
    q = tree / "modules" / "05_quality" / "quality.py"
    q.write_text(q.read_text(encoding="utf-8") + "\n# a modules/ edit, uncommitted\n",
                 encoding="utf-8")
    check("an edit to a module that runs moves the key, committed or not", ident_of() != first,
          "the code that derives a threshold changed, and the key must say so")
    after_module = ident_of()
    a = tree / "adapters" / "scanpy_ops.py"
    a.write_text(a.read_text(encoding="utf-8") + "\n# edited\n", encoding="utf-8")
    check("an edit to an adapter moves it again", ident_of() != after_module)

# The record goes BESIDE the description, never into it: reusing the directory from a second
# commit with an equal description is allowed, and both commits are listed.
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    try:
        d1 = runkey.claim(root, "samekey00000", described,
                          record={"commit": "a" * 40, "dirty": False})
        d2 = runkey.claim(root, "samekey00000", described,
                          record={"commit": "b" * 40, "dirty": False})
        man = json.loads((d2 / runkey.MANIFEST).read_text(encoding="utf-8"))
        commits = [r.get("commit") for r in man.get("code_records") or []]
        check("a second commit with an equal description reuses the directory", d1 == d2)
        check("...and the manifest lists both commits beside the description",
              commits == ["a" * 40, "b" * 40], f"{commits}")
    except TypeError as e:
        check("claim() takes the commit record beside the description", False, str(e))

# 6 - the index. A pre-existing index written without the column must be migrated, not appended
#     to with rows of a different width.
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    old = root / "INDEX.tsv"
    old.write_text("digest\tmode\tsamples\tparameters\tfirst_seen\tnote\n"
                   "oldkey000000\tapply\t2\tlight_floor=200\t2026-08-11T00:00:00+0800\tapply\n",
                   encoding="utf-8")
    runkey.index(root, "newkey000000", described, note="apply",
                 **({"record": CLEAN_A} if "record" in runkey.index.__code__.co_varnames else {}))
    lines = old.read_text(encoding="utf-8").splitlines()
    widths = {len(ln.split("\t")) for ln in lines}
    check("the index has one row width after migration", len(widths) == 1, f"widths {widths}")
    check("...and gained a code column", lines[0].split("\t")[3] == "code", lines[0])
    check("...the old row is marked as predating the column",
          "pre-2026-08-13" in lines[1], lines[1])
    check("...and the new row carries the commit",
          CLEAN_A["commit"][:12] in lines[2], lines[2])

    dirty_desc = dict(described, code={"commit": "c" * 40, "dirty": True})
    runkey.index(root, "dirtykey0000", dirty_desc, note="apply")
    check("a dirty run is marked +dirty in the index",
          "+dirty" in old.read_text(encoding="utf-8"), "a modified tree must be visible")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("run key OK - a code change writes beside the old run, not over its conclusions")
