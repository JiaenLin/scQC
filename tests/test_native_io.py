"""WHY THIS SUITE EXISTS

ONE DEFECT, and it takes down the whole run without a Python traceback.

Under `--executor pbs`, `pipeline.run`'s worker threads read `.h5ad` objects through
anndata/pandas while sibling threads fork for `qsub`/`qstat`. pandas 3.0 builds every categorical's
`categories` Index as an `ArrowStringArray`, and that combination segfaults - repeatedly, and on
commits that had run successfully for months. It is neither a regression nor a resource limit:
what decides it is whether pyarrow is importable in the interpreter, because `mode.string_storage`
defaults to `auto`, which means *use pyarrow if you can*.

`engine.native_io.harden()` selects pandas' python string storage, which removes pyarrow from the
string path so the crashing frame is never entered, and enables `faulthandler` so that if anything
else crashes natively it names a line instead of printing one shell message.

WHAT THIS SUITE CHECKS: that the setting is real, that it actually changes the array class the
crash occurred in, that it is idempotent, that it never raises, and that `main()` calls it before
doing any work. It does NOT try to reproduce a segfault - a test that races threads to trigger a
native crash passes on the machine you debug it on and fails at random elsewhere, and a flaky
suite gets switched off.
"""

from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name if cond else f"{name}  [{detail}]")
    print(f"  {'ok    ' if cond else 'FAILED'} {name}" + (f"   {detail}" if not cond else ""))


from engine import native_io  # noqa: E402

# ------------------------------------------------------------------------- A. it applies at all

first = native_io.harden()
check("harden() returns what it applied", isinstance(first, dict))
check("faulthandler is enabled", first.get("faulthandler") is True)

# ------------------------------------------------- B. it changes the class the crash occurred in

try:
    import pandas as pd
except ImportError:
    print("SKIP: pandas is not installed; the storage assertions need it")
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    sys.exit(1 if FAIL else 0)

check("pandas string storage is 'python' after harden()",
      pd.get_option("mode.string_storage") == "python",
      f"got {pd.get_option('mode.string_storage')!r}")

# THE ASSERTION THAT MATTERS. The crash was in ArrowStringArray._from_sequence, reached through
# Categorical.from_codes -> validate_categories -> Index.__new__. Build exactly that and check
# pyarrow is not in the path any more.
cat = pd.Categorical.from_codes([0, 1, 0], categories=pd.Index(["a", "b"]))
backend = type(cat.categories.array).__name__
check("Categorical.from_codes no longer builds an ArrowStringArray", backend != "ArrowStringArray",
      f"categories backend is {backend} - the crashing frame is still reachable")

# ...and that the option is what decides it, so the check above cannot pass by accident on a build
# where pyarrow is simply absent.
pd.set_option("mode.string_storage", "pyarrow")
arrow_backend = type(
    pd.Categorical.from_codes([0, 1], categories=pd.Index(["x", "y"])).categories.array).__name__
pd.set_option("mode.string_storage", "python")
check("...and the option is what decides it, not the absence of pyarrow",
      arrow_backend != backend,
      f"both storages give {backend} - this environment cannot distinguish them, so this suite "
      f"proves nothing about the fix here")

# ------------------------------------------------------------------------------ C. it is safe

again = native_io.harden()
check("harden() is idempotent", again.get("string_storage") == "python")
check("...and says it had already run", again.get("already") is True)
check("the setting survives a second call", pd.get_option("mode.string_storage") == "python")

# ---------------------------------------------------------------- D. main() actually calls it
#
# A setting nothing invokes is a setting that does not exist. Source grep rather than an import,
# because importing the CLI runs its parser construction.
_cli = (ROOT / "scqc_cli.py").read_text(encoding="utf-8")
_main = _cli.split("\ndef main(", 1)[1].split("\n\ndef ", 1)[0]
check("main() imports harden", "from engine.native_io import harden" in _main)
check("main() calls it BEFORE parsing arguments",
      _main.index("harden()") < _main.index("parse_args"),
      "a crash handler installed after the work has started is not installed")

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
for f in FAIL:
    print(f"  FAILED: {f}")
sys.exit(1 if FAIL else 0)
