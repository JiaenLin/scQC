# The pipeline's own version is part of the run key. Without it, a code change republishes.
"""Does a change to the CODE send a run to a different directory?

WHY THIS TEST EXISTS

engine/runkey.py opens by promising that "nothing is ever overwritten by a run that would have
produced something different". The key covered the samplesheet, the declared parameters and the
mode - everything the operator supplies - and not the code.

A DERIVED threshold is a function of the inputs AND of the code that derives it. So the promise
had a hole in exactly the case that matters most, and the 2026-08-13 mitochondrial change walked
straight into it: same samplesheet, same parameters, same mode, a deliverable 5.9% smaller. A
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

# The real accessor must return something usable on this checkout, or the default is a fiction.
ident = runkey.code_identity()
check("code_identity() reads this checkout",
      isinstance(ident, dict) and isinstance(ident.get("commit"), str)
      and ident["commit"] != "", f"{ident}")

# 6 - the index. A pre-existing index written without the column must be migrated, not appended
#     to with rows of a different width.
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    old = root / "INDEX.tsv"
    old.write_text("digest\tmode\tsamples\tparameters\tfirst_seen\tnote\n"
                   "oldkey000000\tapply\t2\tlight_floor=200\t2026-08-11T00:00:00+0800\tapply\n",
                   encoding="utf-8")
    runkey.index(root, "newkey000000", described, note="apply")
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
