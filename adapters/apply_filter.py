# Applies an agreed keep-mask to a matrix and writes the ledger of what left. This is the only
# adapter in the execution layer that removes anything. It writes NOTHING until
# modules/07_apply/apply.py has accepted the approval, and it writes the ledger BEFORE the
# filtered object, so that a crash between the two leaves a record of a removal that did not
# happen rather than a removal with no record.
"""Step 7 execution - remove what the decision layer approved, and name every observation that
left along with every criterion that removed it.

WHY THIS ADAPTER IS SHAPED AROUND THE LEDGER RATHER THAN AROUND THE WRITE

The removal itself is two lines of anndata. Everything else here exists because a removal that
cannot be enumerated afterwards is a removal nobody can question, and because the question that
decides whether a threshold was worth having is not "how many were removed" but "how many were
removed by THIS criterion ALONE". A ledger that records the first criterion that matched cannot
answer it: a barcode failing the count floor, the gene floor and the mitochondrial ceiling would
be filed under whichever criterion the loop happened to test first, and the mitochondrial
ceiling would then look load-bearing when removing it would have changed nothing. So every row
carries the FULL set of criteria that fired, and `criterion_summary()` reports both the total
per criterion and the count for which that criterion was the only one.

THREE INDEPENDENT ROUTES TO ONE NUMBER

How many observations left is derived three ways and the removal stops if any two disagree:

  1. from the keep-mask         n_in - sum(keep_mask)
  2. from the criteria          the number of observations with at least one criterion firing
  3. from the gate's own record `build_removal_record()` in modules/07_apply/apply.py

Route 1 is what is handed to `apply_removal()` as the mask sum and route 3 is handed to it as
the record, so the gate's existing cross-check does real work: handing it route 3 for both sides
would make it compare a number with itself. Route 2 is checked here as SET equality rather than
by count, because the counts can agree while the identities do not - a mask and a criterion table
that have come apart by one observation each in opposite directions produce the same total and a
ledger naming the wrong barcodes, and nothing downstream can detect that they did.

ORDER OF OPERATIONS, AND WHY IT IS THIS ORDER

  validate -> reconcile -> approval gate -> ledger -> object -> read both artifacts back

Nothing touches the filesystem before the gate. After it, the ledger is written first: if the
process dies between the two artifacts, what remains on disk is a ledger for a removal that never
took effect, which the next run overwrites harmlessly. The reverse order leaves a filtered object
whose removed observations are unnamed, and that state is indistinguishable from a correct one.

The final steps read both artifacts back off disk. `verify_recoverable()` checks that the
identifiers the ledger names, plus the identifiers that were kept, are exactly the identifiers
that went in; and the written object's own observation identifiers are read back out of the .h5ad
and compared, one by one and in order, against the identifiers this call meant to keep. A
recoverability claim that is never checked is the same class of defect as a regression test that
was never run against the bug.

WRITING THE OBJECT: A LOSSLESS CAST, AND A TOKEN RATHER THAN A TIMESTAMP

pandas gives a plain string index the `str`/`StringDtype` backing by default, and anndata refuses
to write those unless `anndata.settings.allow_write_nullable_strings` is turned on - so on a
current stack the deliverable could not be written at all. The setting is global and would change
how every object in the process is stored, which is not a decision this adapter makes on a
caller's behalf, so the labels are re-backed by object arrays instead, with
`adapters/matrix.py::object_backed_labels()`. The cast is lossless - the same `str` objects, a
different array behind them - and the claim is not taken on trust: the identifiers are captured
BEFORE the cast, and after the write they are read back out of the file and compared element by
element. A column carrying missing values is refused rather than given an invented label.

Freshness of that write is established by the token `adapters/matrix.py::write_h5ad()` documents,
stored in the object's `uns` and read back off disk - never by comparing mtimes. A time comparison
needs a tolerance for clock and filesystem resolution, and a tolerance window is exactly where a
restored artifact lands: a copy of a previous run's object back-dated by one second passed the
mtime check that refused the same construction at two hours.

WHAT RECOVERY MEANS HERE, AND WHAT IT DOES NOT

The ledger stores IDENTIFIERS, not counts and not data. Recovering an observation means
re-reading the input matrix with the identifier the ledger names, which is why writing the
filtered object over its own input is refused: it destroys the only copy of what left, and the
ledger then names observations that no longer exist anywhere. The ledger also cannot know about a
criterion applied upstream and not passed in - it then reads as though an observation left for
fewer reasons than it did, and that is the one failure mode a reader of this file should watch
for.

FORMAT

Gzipped CSV, columns `identifier`, `n_criteria`, `criteria` (names, `|`-separated) and one 0/1
column per criterion in the order the criteria were supplied. That is exactly the layout
`write_removal_record()` documents in modules/07_apply/apply.py, compressed - a per-cell ledger
for a cohort runs to millions of rows, and keeping the layout identical means the format
documented there stays the one format a reader has to know. The gzip member is written with
mtime 0 so that identical content hashes identically; a ledger whose bytes change on every write
cannot be compared between runs. No comment lines and no preamble, so `pandas.read_csv` and
`csv.reader(gzip.open(...))` both read it without arguments.

A ledger with a header and no rows is a valid and meaningful artifact: it says the criteria were
evaluated and none of them fired. It is not the same thing as a zero-byte file, which is a failed
write, and the two are told apart by reading the file back rather than by looking at its size.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

try:  # normal import when the repository root is on sys.path
    from engine.task import TaskFailure
except ImportError:  # loaded by file path, as scqc_cli.py loads the step modules
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from engine.task import TaskFailure

__all__ = [
    "Approval",
    "RESERVED_COLUMNS",
    "apply_filter",
    "apply_module",
    "as_bool",
    "coerce_approval",
    "criterion_summary",
    "h5ad_n_obs",
    "h5ad_obs_identifiers",
    "ledger_columns",
    "ledger_row",
    "mask_to_bools",
    "matrix_adapter",
    "normalise_identifier",
    "normalise_identifiers",
    "read_ledger",
    "reconcile_masks",
    "removal_ledger",
    "run_apply_filter",
    "verify_recoverable",
    "write_ledger",
    "written_n_obs",
]

#: The three columns every ledger starts with, in this order. A criterion may not be named after
#: one of them: the row is a flat mapping, so a criterion called `criteria` would overwrite the
#: column that lists which criteria fired and the row would silently lose its own explanation.
RESERVED_COLUMNS = ("identifier", "n_criteria", "criteria")

#: Loaded under the same key scqc_cli.py uses, so both routes share one module object. Loading a
#: second copy under a different name would give `ApplyRefusal` two distinct class objects, and a
#: caller's `except ApplyRefusal` would then fail to catch the refusal raised by this one.
_APPLY_MODULE_KEY = "scqc_apply"
_APPLY_ATTRS = ("apply_removal", "build_removal_record", "ApplyRefusal", "RemovalRecord")

#: `adapters/matrix.py`, loaded the same way, for the two things this adapter needs from it: the
#: lossless object-array cast that lets a current pandas index be written at all, and the write
#: token that says which call produced a file.
_MATRIX_MODULE_KEY = "scqc_matrix"
_MATRIX_ATTRS = ("object_backed_labels", "new_write_token", "read_write_token", "WRITE_TOKEN_KEY")

_EXAMPLES = 5  # how many identifiers a discrepancy message prints before it stops


# ------------------------------------------------------------------------------------------
# a previous run's output is not this run's result


def _clear_previous_output(path, context: str = "") -> float:
    """Delete a declared output BEFORE writing it, and return the time the write started.

    Checking that an output exists after a writer returns proves nothing when the same path
    already held the previous run's file: a writer that returns without writing - a subprocess
    that exits 0 having failed, a library that swallows an error - leaves the old artifact
    standing, it passes the existence check, and the previous parameters' numbers are then
    reported under the new ones. Nothing downstream can see that this happened, because a stale
    .h5ad opens, reads and counts exactly like a fresh one.

    So the path is removed first and its absence is confirmed. Deleting a previous run's output is
    the point: it is this step's own artifact from an earlier invocation, its ledger names what it
    contained, and the alternative to removing it is silently republishing it as the current
    result. Nothing else is touched - the input is never a candidate, `apply_filter()` having
    already refused to write over it.
    """
    path = Path(path)
    try:
        if path.is_symlink() or path.exists():
            path.unlink()
    except OSError as exc:
        raise TaskFailure(
            f"the previous {path} could not be removed before writing this run's: "
            f"{type(exc).__name__}: {exc}\n"
            f"  It is removed first so that a write which produces nothing cannot leave the "
            f"earlier run's file behind to be reported as this run's result.{context}") from exc
    if path.exists():
        raise TaskFailure(
            f"{path} still exists after being deleted. Another process is writing to it, or the "
            f"path is a directory or a mount point. This run will not write there, because its "
            f"output could not then be told from the one already present.{context}")
    return time.time()


def _confirm_fresh_output(path, started: float, token: str, context: str = "") -> None:
    """The output exists, is not empty, and carries THIS call's write token. Silence is the pass.

    `token` was invented by `adapters/matrix.py::new_write_token()` after the path was observed
    absent and handed to the writer to store inside the object; reading it back off disk is the
    second, independent statement that the bytes there are the ones just written.

    It replaces an mtime comparison, which could not make that statement. A clock allowance is
    unavoidable when a compute node and a shared filesystem disagree to the second, and the
    allowance is a window: a copy of a previous run's object back-dated by one second passed the
    check that refused the same construction at two hours, and `os.utime` places a file's mtime
    anywhere at all, so no width of window closes it.

      Detects: a writer that returned without writing and left an earlier artifact standing; a
      file restored, copied, renamed or hardlinked into this path with any mtime whatever; a file
      written by any other producer, which carries no token.

      Does not detect: anything about the contents beyond their provenance. That the object holds
      the right observations is established separately, by reading its identifiers back and
      comparing them with the ones this call meant to keep.

    `started` is used only to describe the file in a refusal. It decides nothing.
    """
    path = Path(path)
    if not path.exists():
        raise TaskFailure(f"{path} is absent after a write that reported success.{context}")
    st = path.stat()
    if st.st_size == 0:
        raise TaskFailure(f"{path} is zero bytes after the write.{context}")
    mx = matrix_adapter()
    stored = mx.read_write_token(path)
    age = started - st.st_mtime
    if stored is None:
        raise TaskFailure(
            f"{path} exists ({st.st_size:,} bytes, mtime {age:+.1f}s relative to the start of "
            f"this write) but carries no {mx.WRITE_TOKEN_KEY} in its /uns, so it is not the file "
            f"this run produced.\n"
            f"  The path was deleted and observed absent immediately before the writer ran and "
            f"the object it was handed carried the token. An .h5ad here without one was written "
            f"by something else - a restored artifact, or another tool - and it does not describe "
            f"this removal.{context}")
    if stored != token:
        raise TaskFailure(
            f"{path} carries the write token {stored!r}, not this run's {token!r} "
            f"({st.st_size:,} bytes, mtime {age:+.1f}s relative to the start of this write).\n"
            f"  It is a different write's output standing at this path: the file was deleted and "
            f"observed absent immediately before the writer ran, so what is here now was restored "
            f"or copied rather than written. It describes a different removal.{context}")


# ------------------------------------------------------------------------------------------
# loading the decision layer


def repo_root() -> Path:
    """The scQC checkout this adapter belongs to. `adapters/` sits directly under it."""
    return Path(__file__).resolve().parents[1]


def apply_module():
    """Load `modules/07_apply/apply.py` - the approval gate - by path.

    The step directories are numeric-prefixed and so are not importable as packages; the CLI
    loads them the same way. The attribute check is not defensiveness for its own sake: this
    adapter is only permitted to remove anything because that module agrees to it, and a version
    of it that no longer exposes the gate must stop the run rather than be worked around.
    """
    mod = sys.modules.get(_APPLY_MODULE_KEY)
    if mod is None:
        path = repo_root() / "modules" / "07_apply" / "apply.py"
        if not path.exists():
            raise TaskFailure(
                f"the approval gate is missing: {path}\n"
                f"  This adapter removes observations only with modules/07_apply/apply.py's "
                f"agreement. Check out the full repository, or run from a tree where that file "
                f"exists - there is no path through this adapter that removes anything without "
                f"it.")
        spec = importlib.util.spec_from_file_location(_APPLY_MODULE_KEY, path)
        if spec is None or spec.loader is None:
            raise TaskFailure(f"could not build an import spec for the approval gate at {path}")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        try:
            spec.loader.exec_module(mod)
        except Exception as exc:  # noqa: BLE001 - the reason is reported, not swallowed
            del sys.modules[spec.name]
            raise TaskFailure(f"the approval gate at {path} failed to import: {exc!r}") from exc
    missing = [a for a in _APPLY_ATTRS if not hasattr(mod, a)]
    if missing:
        raise TaskFailure(
            f"the module loaded as the approval gate does not expose {missing}.\n"
            f"  Expected modules/07_apply/apply.py. Removal stops here rather than proceeding "
            f"through whatever this module is.")
    return mod


def matrix_adapter():
    """Load `adapters/matrix.py`, by ordinary import if that works and by path if it does not.

    Two things are needed from it and neither is duplicated here. `object_backed_labels()` is the
    lossless cast that lets a current pandas string index be written at all - re-implementing it
    beside this one is how the two spellings drift until one of them invents a label for a missing
    value. `new_write_token()` and `read_write_token()` are the identity a written file carries,
    and an adapter that minted its own token in a different format could not read back a file the
    other adapter wrote.

    Loaded under a fixed key so both routes share one module object, and refused when the module
    found does not expose all four names: a checkout whose matrix adapter predates the token would
    otherwise be worked around silently, which is the freshness check quietly turning itself off.
    """
    mod = sys.modules.get(_MATRIX_MODULE_KEY)
    if mod is None:
        try:
            from adapters import matrix as mod  # normal import when the repo root is importable
        except ImportError:
            path = repo_root() / "adapters" / "matrix.py"
            if not path.exists():
                raise TaskFailure(
                    f"the matrix adapter is missing: {path}\n"
                    f"  This adapter writes the filtered object through its lossless label cast "
                    f"and verifies the result with its write token. Check out the full "
                    f"repository, or run from a tree where that file exists.") from None
            spec = importlib.util.spec_from_file_location(_MATRIX_MODULE_KEY, path)
            if spec is None or spec.loader is None:
                raise TaskFailure(f"could not build an import spec for the matrix adapter at "
                                  f"{path}") from None
            mod = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = mod
            try:
                spec.loader.exec_module(mod)
            except Exception as exc:  # noqa: BLE001 - the reason is reported, not swallowed
                del sys.modules[spec.name]
                raise TaskFailure(
                    f"the matrix adapter at {path} failed to import: {exc!r}") from exc
        sys.modules[_MATRIX_MODULE_KEY] = mod
    missing = [a for a in _MATRIX_ATTRS if not hasattr(mod, a)]
    if missing:
        raise TaskFailure(
            f"the module loaded as the matrix adapter does not expose {missing}.\n"
            f"  Expected adapters/matrix.py. The write stops here rather than falling back to a "
            f"weaker check or to a second copy of the label cast.")
    return mod


# ------------------------------------------------------------------------------------------
# unknown is not a value


#: The missing-value sentinels the scientific stack delivers, matched by TYPE NAME so that the
#: predicate below recognises them without importing anything. `pandas.NA` (`NAType`),
#: `pandas.NaT` (`NaTType`) and `numpy.ma.masked` (`MaskedConstant`) are three different objects
#: from two libraries, none of them is None, none of them is a float, and every one of them
#: reaches this module from an ordinary table: `NAType` from a nullable column, `NaTType` from a
#: datetime column, `MaskedConstant` from indexing any masked array at a masked position.
_UNKNOWN_TYPE_NAMES = frozenset(("NAType", "NaTType", "MaskedConstant"))

#: Module roots whose scalars need the library itself to answer the question. Everything outside
#: them is settled by the two cheap tests above it.
_UNKNOWN_LIBRARIES = frozenset(("numpy", "pandas"))


def _library_unknown(v) -> bool:
    """The unknown test for a value that came out of numpy or pandas. One import site, lazy.

    numpy and pandas are imported here rather than at module scope because `scqc_cli.py` imports
    this file on interpreters that have neither, and both imports are absorbed: a value whose type
    lives in numpy cannot exist in a process where numpy failed to import, so an absent library
    can only mean the answer is "not one of that library's sentinels".

    `numpy.isnan` is asked rather than `v != v` because the numpy scalars that are not float
    subclasses - `float32`, `datetime64('NaT')`, `timedelta64('NaT')` - answer the first correctly
    and the second inconsistently, and because `numpy.ma.masked != numpy.ma.masked` evaluates to
    `masked`, which is FALSY: the self-inequality spelling of this test reports the standard
    masked sentinel as a known value.
    """
    try:
        import numpy as np
    except ImportError:  # pragma: no cover - exercised only on a stdlib-only interpreter
        np = None
    if np is not None:
        masked = getattr(getattr(np, "ma", None), "masked", None)
        if masked is not None and v is masked:
            return True
        if isinstance(v, np.generic):
            try:
                return bool(np.isnan(v))
            except (TypeError, ValueError):
                return False  # a numpy string or object scalar: not a number, not a NaN
    try:
        import pandas as pd
    except ImportError:  # pragma: no cover - as above
        return False
    try:
        answer = pd.isna(v)
    except (TypeError, ValueError):
        return False
    if isinstance(answer, bool):
        return answer
    if np is not None and isinstance(answer, np.bool_):
        return bool(answer)
    return False  # an array answer means v was a container, which is not a missing scalar


def _is_unknown(v) -> bool:
    """True when a value carries no information. THE unknown test for this module - use no other.

    One predicate, because every place that asks "was this measured?" has to agree, and the
    sentinels are not one thing: `None`, `float('nan')`, `numpy.float64('nan')`,
    `numpy.float32('nan')`, `pandas.NA`, `pandas.NaT` and `numpy.ma.masked` all mean unmeasured
    and only the first is None, only some are floats, and no two of them are equal.

    Each one that slips through fails in the same direction, which is why the list is exhaustive
    rather than the two spellings that were met first: `nan is not None` is True and `nan >= x` is
    False, so a NaN passes an `is not None` guard and then reads as a value that FAILED the test -
    a pass. `numpy.ma.masked != numpy.ma.masked` is falsy, so the self-inequality test reports it
    as known. `bool(pandas.NA)` raises rather than returning False, so a missing cell reached
    through truthiness crashes somewhere unrelated to its cause.

    Truthiness of a possibly-numpy boolean belongs downstream of this call and is spelled
    `bool(x)`, never `x is True`: `numpy.bool_(True) is True` is False, so an identity test reads
    a genuinely flagged row from any numpy-backed table as not flagged.
    """
    if v is None:
        return True
    if type(v).__name__ in _UNKNOWN_TYPE_NAMES:
        return True
    if isinstance(v, float):
        # bool(): numpy.float64 is a float subclass and its comparison returns numpy.bool_, and a
        # predicate that answers numpy.bool_ propagates the same identity trap it exists to close.
        return bool(v != v)
    if getattr(type(v), "__module__", "").split(".")[0] in _UNKNOWN_LIBRARIES:
        return _library_unknown(v)
    try:
        return bool(v != v)  # any other library's NaN scalar, e.g. decimal.Decimal('NaN')
    except Exception:  # noqa: BLE001 - an exotic __ne__ is not evidence of missingness
        return False


#: Returned when a value is not a numpy scalar. A sentinel rather than None, because None is
#: itself one of the values this module has to be able to tell apart from "no answer".
_NO_UNWRAP = object()


def _numpy_scalar_value(v):
    """Unwrap a numpy scalar to its Python value, or return the sentinel `_NO_UNWRAP`.

    `numpy.bool_` is not a subclass of `bool`, so `isinstance(mask[i], bool)` is False for every
    entry of a perfectly ordinary boolean numpy array. Rejecting those would refuse the most
    common input this adapter is given, and a check that fires on correct behaviour gets removed.
    """
    if getattr(type(v), "__module__", "").split(".")[0] != "numpy":
        return _NO_UNWRAP
    item = getattr(v, "item", None)
    if not callable(item):
        return _NO_UNWRAP
    try:
        u = item()
    except (TypeError, ValueError):
        return _NO_UNWRAP
    return _NO_UNWRAP if u is v else u


def as_bool(v, what: str) -> bool:
    """A mask entry is True or False, or the run stops. Never a default, never a truthiness test.

    Refused deliberately, each for a failure this pipeline has met or is one keystroke away from:
    an unknown entry, because removing on a criterion that was never evaluated - or recording it
    as not having fired - states more than was measured; a number other than 0 or 1, because a
    count that arrived where a mask was expected would otherwise read as True for every non-zero
    observation; and a string, because `bool("False")` is True and a mask read from a CSV without
    conversion removes everything.
    """
    if _is_unknown(v):
        raise TaskFailure(
            f"{what} is unknown ({v!r}).\n"
            f"  An unknown is not a False. Removing an observation on a criterion that was never "
            f"evaluated, and recording it as not having fired, both state more than was "
            f"measured. Compute the criterion for every observation, or drop the criterion.")
    if isinstance(v, bool):
        return v
    if isinstance(v, (int, float)):
        if v == 0:
            return False
        if v == 1:
            return True
        raise TaskFailure(
            f"{what} is {v!r}, which is a number but not a mask entry.\n"
            f"  A criterion mask holds True/False or 0/1 per observation. A count or a score has "
            f"arrived where a mask was expected; every non-zero entry would read as 'remove'.")
    u = _numpy_scalar_value(v)
    if u is not _NO_UNWRAP:
        return as_bool(u, what)
    raise TaskFailure(
        f"{what} is {v!r} of type {type(v).__name__}, which is not a boolean.\n"
        f"  Strings are refused rather than coerced: bool('False') is True, so a mask read from "
        f"a text column without conversion removes every observation and reports it as intended. "
        f"Convert the column to booleans at the point it is read.")


# ------------------------------------------------------------------------------------------
# identifiers


def normalise_identifier(x, what: str = "identifier") -> str:
    """The single spelling of an identifier used in the ledger, on disk and in every comparison.

    Everything the ledger promises rests on one identifier meaning one observation in all three
    places, and a CSV round trip returns text. Fixing the spelling here - once, at the boundary -
    is what allows `verify_recoverable()` to compare what was written against what was passed in
    without the comparison quietly succeeding because both sides were coerced differently.

    Leading and trailing whitespace is refused rather than stripped, because an identifier that
    survives this function but not a spreadsheet, a `cut`, or another CSV reader is an identifier
    that recovers the wrong observation later, in a session where nobody is looking for it.
    """
    if _is_unknown(x):
        raise TaskFailure(
            f"{what} is unknown ({x!r}).\n"
            f"  An observation with no identifier cannot be named in the ledger, so its removal "
            f"could not be undone. Supply the barcode - or the sample-plus-barcode key - for "
            f"every observation.")
    if isinstance(x, (bytes, bytearray)):
        try:
            s = bytes(x).decode("utf-8")
        except UnicodeDecodeError as exc:
            raise TaskFailure(
                f"{what} is bytes that are not UTF-8: {x!r}. Decode the identifiers with the "
                f"encoding they were written in before passing them here.") from exc
    else:
        u = _numpy_scalar_value(x)
        s = str(u if u is not _NO_UNWRAP else x)
    if not s.strip():
        raise TaskFailure(
            f"{what} is blank ({x!r}).\n"
            f"  A blank key names no observation and cannot be told from a missing cell when the "
            f"ledger is read back.")
    if s != s.strip():
        raise TaskFailure(
            f"{what} has leading or trailing whitespace: {s!r}.\n"
            f"  It is refused rather than trimmed here, because other readers of the ledger trim "
            f"it and this one would not, and the two would then disagree about which observation "
            f"a row names. Trim the identifiers upstream.")
    if "\n" in s or "\r" in s:
        raise TaskFailure(
            f"{what} contains a line break: {s!r}.\n"
            f"  A ledger row must stay one line, or every line-oriented tool that reads it - "
            f"grep, head, wc - reports a different number of removals than the CSV parser does.")
    return s


def normalise_identifiers(obs_ids, what: str = "identifier") -> list:
    """Normalise a sequence of identifiers and refuse duplicates.

    Duplicates are refused rather than merged: a ledger keyed on an ambiguous identifier cannot
    answer "what if that one had not been dropped", which is the only question it exists to
    answer. Barcodes repeat across samples, so the key for a multi-sample object is the sample
    combined with the barcode - built by the caller, who knows which sample this is, and not
    invented here.
    """
    if isinstance(obs_ids, (str, bytes, dict)):
        raise TaskFailure(
            f"{what}s must be a sequence of identifiers, one per observation; got "
            f"{type(obs_ids).__name__}.")
    try:
        raw = list(obs_ids)
    except TypeError as exc:
        raise TaskFailure(
            f"{what}s is not iterable ({type(obs_ids).__name__}); expected one identifier per "
            f"observation.") from exc
    out = [normalise_identifier(x, f"{what}[{i}]") for i, x in enumerate(raw)]
    seen, dupe_set, dupes = set(), set(), []
    for s in out:
        if s in seen and s not in dupe_set:
            dupe_set.add(s)
            dupes.append(s)
        seen.add(s)
    if dupes:
        raise TaskFailure(
            f"{len(dupes)} duplicate {what}(s) among {len(out)} observations, e.g. "
            f"{dupes[:_EXAMPLES]!r}.\n"
            f"  A ledger keyed on an ambiguous identifier cannot say which of the observations "
            f"sharing that key was removed. Barcodes repeat across samples: build the key as "
            f"sample plus barcode before calling this adapter, rather than merging them here.")
    return out


# ------------------------------------------------------------------------------------------
# masks and criteria


def _labels_are_positional(labels, n: int) -> bool:
    """True when an index is the default 0..n-1 and therefore carries no observation identity.

    `pandas.Series([...])` with no index argument gets a RangeIndex, and a RangeIndex says nothing
    about which observation each entry belongs to - the position is all there is, so consuming it
    positionally is not an assumption, it is the only reading available.
    """
    if len(labels) != n:
        return False
    for i, lab in enumerate(labels):
        if isinstance(lab, bool) or not isinstance(lab, int):
            u = _numpy_scalar_value(lab)
            if u is _NO_UNWRAP or isinstance(u, bool) or not isinstance(u, int):
                return False
            lab = u
        if lab != i:
            return False
    return True


def _labelled_mask_values(mask, labels, n: int, obs_ids, what: str) -> list:
    """Values of a labelled mask (a pandas Series) in OBSERVATION order, or a refusal.

    `list(series)` reads a Series positionally and throws its index away. When the index is a set
    of barcodes and the caller's observations are in another order - which is what a `.groupby`,
    a `.reindex`, a merge or a scanpy round trip routinely produces - that silently records the
    wrong observations, which is precisely the failure the docstring above promises to catch.

    So the index is used when there is one to use:

      * index equal to the observation ids, in any order -> the values are taken BY LABEL;
      * a default 0..n-1 RangeIndex -> positional, because such an index holds no identity;
      * anything else -> refused, naming the labels that do not match.

    With `obs_ids=None` there is nothing to align against, so only the second case is accepted. A
    labelled mask and no identifiers is not something this function can check, and reading it
    positionally anyway is the defect, not the fallback.
    """
    values = list(mask)
    if len(values) != n:
        return values  # the caller's length check reports this; alignment cannot fix a length
    if obs_ids is None:
        if _labels_are_positional(labels, n):
            return values
        raise TaskFailure(
            f"{what} is a labelled {type(mask).__name__} and no observation identifiers were "
            f"supplied to check its index against.\n"
            f"  Reading it positionally would discard the labels, and a mask whose order differs "
            f"from the observations records the wrong ones while every count still agrees. Pass "
            f"obs_ids=<the identifiers, in observation order>, or hand in a plain sequence.")
    ids = list(obs_ids)
    try:
        keys = [normalise_identifier(lab, f"{what} index[{i}]") for i, lab in enumerate(labels)]
    except TaskFailure:
        if _labels_are_positional(labels, n):
            return values
        raise
    pos = {}
    dupes = []
    for i, k in enumerate(keys):
        if k in pos:
            dupes.append(k)
        else:
            pos[k] = i
    if dupes:
        raise TaskFailure(
            f"{what} has a repeated index label, e.g. {sorted(set(dupes))[:_EXAMPLES]!r}.\n"
            f"  Which entry belongs to that observation cannot be decided, so the mask cannot be "
            f"aligned. De-duplicate it upstream.")
    key_set, id_set = set(keys), set(ids)
    if key_set == id_set:
        return [values[pos[i]] for i in ids]
    if _labels_are_positional(labels, n):
        return values
    only_mask = sorted(key_set.difference(id_set))
    only_obs = sorted(id_set.difference(key_set))
    raise TaskFailure(
        f"{what} is indexed by labels that are not the observation identifiers.\n"
        f"  {len(only_mask):,} label(s) appear only in the mask, e.g. "
        f"{only_mask[:_EXAMPLES]!r}\n"
        f"  {len(only_obs):,} observation(s) appear only in the matrix, e.g. "
        f"{only_obs[:_EXAMPLES]!r}\n"
        f"  It is refused rather than consumed positionally: the two are not the same removal, "
        f"and a positional read of a mismatched index names the wrong observations in a ledger "
        f"that agrees with itself. Rebuild the mask against this object's observations.")


def mask_to_bools(mask, n: int, what: str, obs_ids=None) -> list:
    """Turn a per-observation mask into exactly `n` real booleans, IN OBSERVATION ORDER.

    Length is checked before content. A mask that does not line up with the identifiers records
    the wrong observations, and nothing downstream can detect that it did - the ledger is
    internally consistent, every count agrees, and the named barcodes are the wrong ones.

    That protection is why `obs_ids` exists. A pandas Series carries its own index, and consuming
    it positionally is exactly the misalignment being guarded against, so a labelled mask is
    aligned on its index when the identifiers are available and REFUSED when it cannot be checked.
    Every caller inside this module passes the identifiers; the argument is optional only so that
    a plain sequence - which has no index and can only be positional - still works unchanged.
    """
    if mask is None:
        raise TaskFailure(
            f"{what} is missing.\n"
            f"  There is no default mask. A removal decided by an absent criterion is not a "
            f"decision.")
    if isinstance(mask, (str, bytes, dict, set, frozenset)):
        raise TaskFailure(
            f"{what} is a {type(mask).__name__}, which is not a per-observation mask.\n"
            f"  Iterating a dict yields its keys and iterating a string yields characters; "
            f"either would produce a mask of the right length made of the wrong things. Pass a "
            f"sequence of booleans in observation order.")
    labels = getattr(mask, "index", None)
    if labels is not None and hasattr(mask, "columns"):
        raise TaskFailure(
            f"{what} is a {type(mask).__name__} with columns {list(mask.columns)[:5]!r}, not a "
            f"per-observation mask.\n"
            f"  Iterating a DataFrame yields its column names. Select the column you mean and "
            f"pass that.")
    if labels is not None:
        try:
            labels = list(labels)
        except TypeError:
            labels = None
    if labels is not None:
        seq = _labelled_mask_values(mask, labels, n, obs_ids, what)
    else:
        try:
            seq = list(mask)
        except TypeError as exc:
            raise TaskFailure(
                f"{what} is not a sequence of booleans ({type(mask).__name__}); expected one "
                f"entry per observation.") from exc
    if len(seq) != n:
        raise TaskFailure(
            f"{what} has {len(seq):,} entries for {n:,} observations.\n"
            f"  A mask that does not line up with the identifiers names the wrong observations, "
            f"and every count computed from it still agrees with itself. Check that the mask and "
            f"the matrix are in the same order and from the same run.")
    return [as_bool(v, f"{what}[{i}]") for i, v in enumerate(seq)]


def _criteria_items(criteria) -> list:
    """Normalise `criteria` to an ordered list of (name, mask) and validate the names.

    Order is the caller's insertion order, not sorted: it becomes the column order of the ledger,
    and a caller comparing two ledgers should see the columns they supplied in the order they
    supplied them.
    """
    if criteria is None:
        raise TaskFailure(
            "no criteria were supplied.\n"
            "  A ledger that cannot say WHY an observation left is a list of casualties, not a "
            "record. Pass at least one named criterion, including any that removed nothing - a "
            "criterion with zero rows is evidence it was evaluated.")
    if hasattr(criteria, "items"):
        items = list(criteria.items())
    else:
        try:
            items = [(k, v) for k, v in criteria]
        except (TypeError, ValueError) as exc:
            raise TaskFailure(
                "criteria must be a mapping of criterion name to boolean mask, or a sequence of "
                "(name, mask) pairs.") from exc
    if not items:
        raise TaskFailure(
            "the criteria mapping is empty.\n"
            "  A removal with no named criterion cannot be explained afterwards, and this "
            "adapter does not invent a name for it.")
    names = []
    for name, _ in items:
        if not isinstance(name, str) or not name.strip():
            raise TaskFailure(
                f"criterion name {name!r} is not a non-empty string. The name is what a reader "
                f"of the ledger sees as the reason an observation left.")
        if name != name.strip():
            raise TaskFailure(
                f"criterion name {name!r} has leading or trailing whitespace, which does not "
                f"survive a round trip through every CSV reader. Trim it.")
        if name in RESERVED_COLUMNS:
            raise TaskFailure(
                f"criterion name {name!r} collides with a reserved ledger column "
                f"{RESERVED_COLUMNS}.\n"
                f"  A row is a flat mapping: this criterion would overwrite the column that "
                f"lists which criteria fired, and the row would lose its own explanation. Rename "
                f"the criterion.")
        if "|" in name:
            raise TaskFailure(
                f"criterion name {name!r} contains '|', the separator used by the `criteria` "
                f"column. The set of criteria could not be parsed back out of a written ledger. "
                f"Rename the criterion.")
        if any(c in name for c in ("\n", "\r")):
            raise TaskFailure(
                f"criterion name {name!r} contains a line break, which would split the header "
                f"row across two lines. Rename the criterion.")
        if name in names:
            raise TaskFailure(
                f"criterion {name!r} is supplied twice. Two different masks under one name "
                f"cannot both be recorded, and silently keeping the last one records a removal "
                f"under a reason that did not cause it.")
        names.append(name)
    return items


def _normalise_criteria(criteria, ids):
    """(ordered names, {name: [bool, ...]}) - the one place criteria are validated.

    Takes the identifiers rather than a count so that every criterion mask is aligned against the
    observations it claims to describe; a labelled mask read positionally is the misalignment
    `mask_to_bools()` exists to refuse.
    """
    items = _criteria_items(criteria)
    names = [name for name, _ in items]
    masks = {name: mask_to_bools(mask, len(ids), f"criterion {name!r}", obs_ids=ids)
             for name, mask in items}
    return names, masks


# ------------------------------------------------------------------------------------------
# the ledger


def ledger_columns(criteria_names) -> list:
    """Column order for a ledger over these criteria. Pure, so it can be checked without a file."""
    return list(RESERVED_COLUMNS) + list(criteria_names)


def ledger_row(identifier: str, fired, criteria_names) -> dict:
    """One ledger row: the identifier, and EVERY criterion that removed it.

    Both spellings of the same fact are written - the `|`-joined names and a 0/1 column per
    criterion - because they are what makes the file usable from two directions: a reader
    grepping for one criterion wants the column, and a reader asking why a particular barcode
    went wants the list. `read_ledger()` checks the two against each other on the way back in, so
    the redundancy is a check rather than a duplication.
    """
    fired_set = set(fired)
    unknown = sorted(fired_set.difference(criteria_names))
    if unknown:
        raise TaskFailure(
            f"row for {identifier!r} cites criteria that are not in the criteria list: "
            f"{unknown!r}. A ledger column that does not exist cannot be read back.")
    row = {"identifier": identifier,
           "n_criteria": len(fired_set),
           "criteria": "|".join(n for n in criteria_names if n in fired_set)}
    for name in criteria_names:
        row[name] = 1 if name in fired_set else 0
    return row


def _rows_from_masks(ids, names, masks) -> list:
    """Build the rows from already-validated identifiers and masks.

    Shared by `removal_ledger()` and `apply_filter()` so that the rows written are built from the
    same validated masks the reconciliation and the gate's record were built from, and not from a
    second pass that could see something different.
    """
    rows = []
    for pos, ident in enumerate(ids):
        fired = [name for name in names if masks[name][pos]]
        if fired:
            rows.append(ledger_row(ident, fired, names))
    return rows


def removal_ledger(obs_ids, criteria) -> list:
    """One row per REMOVED observation, naming every criterion that removed it.

    An observation is removed when ANY criterion fires, and it carries ALL of them. Recording
    only the first match is the shape of this that looks correct and is not: it makes the
    per-criterion tallies depend on the order the criteria happen to be evaluated in, and it
    makes the question that decides whether a threshold was worth keeping - how many observations
    that criterion removed ON ITS OWN - unanswerable from the record.

    Returns a list of plain dicts in observation order, ready for `write_ledger()`; it touches no
    files, so a caller can inspect exactly what would be recorded before anything is approved.
    """
    ids = normalise_identifiers(obs_ids)
    names, masks = _normalise_criteria(criteria, ids)
    return _rows_from_masks(ids, names, masks)


def _flag(v, what: str) -> int:
    """A 0/1 ledger cell. A blank cell is unknown and is refused, not read as zero."""
    if _is_unknown(v):
        raise TaskFailure(
            f"{what} is blank.\n"
            f"  A blank cell in a criterion column is not a zero: it cannot be told from a "
            f"criterion that was never evaluated, and reading it as 'did not fire' understates "
            f"why the observation left.")
    if isinstance(v, str):
        s = v.strip()
        if s not in ("0", "1"):
            raise TaskFailure(f"{what} is {v!r}; a criterion column holds 0 or 1.")
        return int(s)
    return 1 if as_bool(v, what) else 0


def _fired_from_row(row, criteria_names, where: str) -> list:
    """The criteria that fired for a row, checked against both spellings of the answer."""
    for key in ("criteria", "n_criteria"):
        if key not in row:
            raise TaskFailure(f"{where} has no {key!r} column; this is not a scQC removal ledger.")
    absent = [name for name in criteria_names if name not in row]
    if absent:
        raise TaskFailure(
            f"{where} is missing the criterion column(s) {absent!r}. An absent column is not a "
            f"zero - it says the criterion was never recorded for this observation, which is a "
            f"different statement from 'it did not fire'.")
    from_cols = [name for name in criteria_names if _flag(row[name], f"{where}[{name!r}]")]
    listed_raw = row["criteria"]
    if _is_unknown(listed_raw):
        raise TaskFailure(
            f"{where} has a blank `criteria` column. Every removed observation carries the "
            f"reason it left; a row without one is not recoverable as a decision.")
    listed = [c for c in str(listed_raw).split("|") if c]
    if sorted(listed) != sorted(from_cols):
        raise TaskFailure(
            f"{where} disagrees with itself: the `criteria` column says {listed!r} and the 0/1 "
            f"columns say {from_cols!r}.\n"
            f"  The ledger has been edited or truncated since it was written. Re-run the "
            f"removal rather than trusting either half.")
    try:
        declared = int(str(row["n_criteria"]).strip())
    except (TypeError, ValueError) as exc:
        raise TaskFailure(
            f"{where} has n_criteria={row['n_criteria']!r}, which is not a count.") from exc
    if declared != len(from_cols):
        raise TaskFailure(
            f"{where} declares n_criteria={declared} but {len(from_cols)} criteria fired. The "
            f"row's own summary and its columns disagree; the file has been altered.")
    if not from_cols:
        raise TaskFailure(
            f"{where} names no criterion at all.\n"
            f"  The ledger holds removed observations only, so a row with nothing firing means "
            f"either an observation was removed for a reason that was never recorded, or a kept "
            f"observation was written into the ledger. Both are unrecoverable as they stand.")
    return from_cols


def criterion_summary(rows, criteria_names) -> dict:
    """How much each criterion actually did - including how often it acted alone.

    `by_criterion` counts every observation a criterion fired on, and those counts overlap, so
    they sum to more than the total removed. `sole_criterion` counts the observations a criterion
    removed that NO other criterion would have removed, and that is the number that says whether
    a threshold mattered: a mitochondrial ceiling with a large total and a sole count of zero
    removed nothing the count floors were not already removing, and moving it would change
    nothing at all.
    """
    names = list(criteria_names)
    by = {n: 0 for n in names}
    sole = {n: 0 for n in names}
    hist = {}
    for i, row in enumerate(rows):
        fired = _fired_from_row(row, names, f"ledger row {i}")
        for n in fired:
            by[n] += 1
        if len(fired) == 1:
            sole[fired[0]] += 1
        hist[len(fired)] = hist.get(len(fired), 0) + 1
    return {"n_removed": len(rows),
            "by_criterion": by,
            "sole_criterion": sole,
            "by_n_criteria": {k: hist[k] for k in sorted(hist)}}


def write_ledger(rows, path, columns=None) -> Path:
    """Write the ledger as gzipped CSV and read it back before claiming it exists.

    Written to a temporary file and renamed into place, so an interrupted write cannot leave a
    truncated ledger under the final name where the next reader takes it for the whole record.
    The verification is a full re-read rather than a size check: decompressing the file end to
    end also checks the gzip CRC, and a header-only ledger - a real outcome, meaning every
    criterion was evaluated and none fired - is thereby told apart from a zero-byte file, which
    is a failed write. The two look identical to any test that only asks whether the path exists.

    The temporary file is also what keeps a previous run's ledger from being accepted as this
    one's: `os.replace` puts bytes written by THIS call at the final name or raises, so there is
    no path through here where the file that is then read back and reported is the earlier run's.
    That is the same requirement `_clear_previous_output()` meets for the filtered object, met by
    a different mechanism because this writer, unlike anndata's, is ours.

    `columns` is required when `rows` is empty, because the header is then the only evidence left
    of which criteria were evaluated, and it cannot be inferred from no rows.
    """
    path = Path(path)
    if path.suffix != ".gz":
        raise TaskFailure(
            f"the ledger path {path} does not end in .gz but the file is gzip-compressed.\n"
            f"  Every reader that opens it by name would get binary. Name it '*.csv.gz'.")
    if columns is None:
        if not rows:
            raise TaskFailure(
                "an empty ledger needs its columns passed explicitly.\n"
                "  A header naming every criterion is the only record that the criteria were "
                "evaluated and none of them fired; inferred from zero rows, it would say nothing "
                "was checked.")
        columns = list(rows[0].keys())
    columns = list(columns)
    if columns[:len(RESERVED_COLUMNS)] != list(RESERVED_COLUMNS):
        raise TaskFailure(
            f"ledger columns must begin {list(RESERVED_COLUMNS)}; got {columns[:3]}.")
    for i, row in enumerate(rows):
        extra = sorted(set(row).difference(columns))
        missing = sorted(set(columns).difference(row))
        if extra or missing:
            raise TaskFailure(
                f"ledger row {i} does not match the columns: missing {missing}, unexpected "
                f"{extra}. Every row must carry every criterion, including the ones that did not "
                f"fire for it, or a 0 cannot be told from an absence.")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    try:
        with open(tmp, "wb") as raw:
            # filename="" and mtime=0: the gzip header otherwise records the write time and the
            # source name, so two ledgers with identical content would hash differently and no
            # provenance record could show that a rerun produced the same removal.
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as fh:
                    w = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n")
                    w.writeheader()
                    for row in rows:
                        w.writerow(row)
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(tmp, path)
    except OSError as exc:
        if Path(tmp).exists():
            try:
                os.unlink(tmp)
            except OSError:
                pass
        raise TaskFailure(f"could not write the removal ledger to {path}: {exc}") from exc
    if not path.exists():
        raise TaskFailure(f"the removal ledger is absent after a write that reported success: "
                          f"{path}")
    size = path.stat().st_size
    if size == 0:
        raise TaskFailure(f"the removal ledger {path} is zero bytes; the write produced nothing.")
    back = read_ledger(path)
    if back["columns"] != columns:
        raise TaskFailure(
            f"the ledger read back from {path} has columns {back['columns']}, not {columns}.")
    if len(back["rows"]) != len(rows):
        raise TaskFailure(
            f"the ledger read back from {path} holds {len(back['rows']):,} rows, not the "
            f"{len(rows):,} that were written. The file on disk is not the record that was built.")
    return path


def read_ledger(path) -> dict:
    """Read a ledger back and check it against itself.

    Returns `{"columns", "criteria", "rows", "identifiers"}`. Every row is validated on the way
    through - the `|`-joined names against the 0/1 columns, the declared count against both - so
    a caller never has to decide which half of a row to believe.
    """
    path = Path(path)
    if not path.exists():
        raise TaskFailure(
            f"no removal ledger at {path}.\n"
            f"  Without it the removal cannot be enumerated, so it cannot be undone or "
            f"questioned. Re-run the apply step; do not proceed from the filtered object alone.")
    try:
        with gzip.open(str(path), "rt", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh)
            try:
                header = next(reader)
            except StopIteration:
                raise TaskFailure(
                    f"the removal ledger {path} is empty - not even a header row. An empty file "
                    f"cannot be told from 'nothing was removed'; the two are different outcomes "
                    f"and this one is a failed write.") from None
            body = [r for r in reader]
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise TaskFailure(f"the removal ledger {path} could not be read as gzip: {exc}") from exc
    if header[:len(RESERVED_COLUMNS)] != list(RESERVED_COLUMNS):
        raise TaskFailure(
            f"{path} does not look like a scQC removal ledger: its first columns are "
            f"{header[:3]}, expected {list(RESERVED_COLUMNS)}.")
    criteria = header[len(RESERVED_COLUMNS):]
    if not criteria:
        raise TaskFailure(
            f"{path} names no criteria at all. A ledger without criterion columns records that "
            f"observations left without recording why.")
    rows, ids = [], []
    for i, raw in enumerate(body):
        if len(raw) != len(header):
            raise TaskFailure(
                f"{path} row {i + 2} has {len(raw)} fields for {len(header)} columns. The file "
                f"is truncated or was edited; a short row silently shifts every criterion by one.")
        row = dict(zip(header, raw))
        row["identifier"] = normalise_identifier(row["identifier"], f"{path} row {i + 2}")
        _fired_from_row(row, criteria, f"{path} row {i + 2}")
        rows.append(row)
        ids.append(row["identifier"])
    # Counted in one pass. `ids.count(x)` inside a comprehension is quadratic, and a cohort
    # ledger runs to millions of rows - the check would appear to hang rather than to fail.
    seen, dupe_set, dupes = set(), set(), []
    for x in ids:
        if x in seen and x not in dupe_set:
            dupe_set.add(x)
            dupes.append(x)
        seen.add(x)
    if dupes:
        raise TaskFailure(
            f"{path} names {len(dupes)} identifier(s) more than once, e.g. "
            f"{dupes[:_EXAMPLES]!r}. A removed observation appears once; a repeated key means "
            f"two ledgers were concatenated or the identifier is not unique.")
    return {"columns": header, "criteria": criteria, "rows": rows, "identifiers": ids}


# ------------------------------------------------------------------------------------------
# reconciliation and the recoverability check


def reconcile_masks(obs_ids, keep_mask, criteria) -> dict:
    """Check that the keep-mask and the criteria describe the SAME removal, by identity.

    Two routes to what leaves - the mask that will subset the matrix, and the criteria that will
    be recorded as the reasons - and they are compared as sets of identifiers rather than as
    counts. Equal counts over different observations is a real state: one observation dropped by
    the mask and not by the criteria, another by the criteria and not by the mask, every total in
    agreement, and a ledger naming barcodes that were not the ones removed.

    Returns the two identifier sets and the kept set; raises on any disagreement.
    """
    ids = normalise_identifiers(obs_ids)
    n = len(ids)
    keep = mask_to_bools(keep_mask, n, "keep_mask", obs_ids=ids)
    names, masks = _normalise_criteria(criteria, ids)

    dropped_by_mask = {ids[i] for i in range(n) if not keep[i]}
    kept_by_mask = {ids[i] for i in range(n) if keep[i]}
    dropped_by_criteria = {ids[i] for i in range(n)
                           if any(masks[name][i] for name in names)}

    unexplained = sorted(dropped_by_mask.difference(dropped_by_criteria))
    unapplied = sorted(dropped_by_criteria.difference(dropped_by_mask))
    if unexplained or unapplied:
        hint = ""
        if dropped_by_criteria and not dropped_by_criteria.intersection(dropped_by_mask):
            # Zero overlap, not equality: criteria supplied in pass/keep semantics fire on the
            # kept observations and on nothing the mask drops, but their union is usually a
            # superset of the kept set rather than equal to it, so an equality test would miss
            # the very case it was written for. Requiring zero overlap keeps it from firing on
            # masks that are merely a little out of step, which always overlap.
            hint = ("\n  Every observation a criterion fired on is one the keep-mask KEEPS, and "
                    "no criterion fired on anything the mask drops. That is the signature of "
                    "criteria supplied in pass/keep semantics: in `criteria`, True must mean "
                    "'this criterion REMOVES this observation', the opposite of `keep_mask`.")
        raise TaskFailure(
            f"the keep-mask and the criteria describe different removals.\n"
            f"  {len(unexplained):,} observation(s) are dropped by the mask with no criterion "
            f"recorded, e.g. {unexplained[:_EXAMPLES]!r} - those removals could never be "
            f"explained afterwards.\n"
            f"  {len(unapplied):,} observation(s) have a criterion firing but are kept by the "
            f"mask, e.g. {unapplied[:_EXAMPLES]!r} - the ledger would record a removal that did "
            f"not happen.\n"
            f"  Both come from a mask and a criterion table built at different times or in "
            f"different orders. Rebuild them together.{hint}")
    return {"identifiers": ids,
            "keep": keep,
            "kept": kept_by_mask,
            "removed": dropped_by_mask,
            "criteria_names": names,
            "criteria_masks": masks}


def verify_recoverable(ledger_path, original_ids, kept_ids) -> None:
    """Read the ledger off disk and prove that removed + kept accounts for every input.

    This is the claim the pipeline makes in its README, checked against the artifact rather than
    against the variables that produced it: the file is re-read, so a ledger written to the wrong
    path, truncated, or written from a stale row list fails here rather than being discovered by
    whoever later tries to recover an observation and finds it named nowhere.

    Returns None and raises `TaskFailure` on any discrepancy. Silence is the pass.
    """
    ledger = read_ledger(ledger_path)
    removed = ledger["identifiers"]
    original = normalise_identifiers(original_ids, "original identifier")
    kept = normalise_identifiers(kept_ids, "kept identifier")

    removed_set, original_set, kept_set = set(removed), set(original), set(kept)
    problems = []

    both = sorted(removed_set.intersection(kept_set))
    if both:
        problems.append(
            f"{len(both):,} identifier(s) are recorded as removed AND kept, e.g. "
            f"{both[:_EXAMPLES]!r}")
    strangers = sorted(removed_set.difference(original_set))
    if strangers:
        problems.append(
            f"{len(strangers):,} ledger identifier(s) were never in the input, e.g. "
            f"{strangers[:_EXAMPLES]!r} - the ledger describes a different object")
    kept_strangers = sorted(kept_set.difference(original_set))
    if kept_strangers:
        problems.append(
            f"{len(kept_strangers):,} kept identifier(s) were never in the input, e.g. "
            f"{kept_strangers[:_EXAMPLES]!r}")
    lost = sorted(original_set.difference(removed_set).difference(kept_set))
    if lost:
        problems.append(
            f"{len(lost):,} input identifier(s) are neither kept nor named in the ledger, e.g. "
            f"{lost[:_EXAMPLES]!r} - those observations left with no record of why")
    if len(removed) + len(kept) != len(original):
        problems.append(
            f"{len(removed):,} removed + {len(kept):,} kept = {len(removed) + len(kept):,}, "
            f"but {len(original):,} went in")
    if problems:
        raise TaskFailure(
            "the removal is not recoverable from the ledger at "
            f"{ledger_path}:\n  - " + "\n  - ".join(problems) +
            "\n  Every removed observation must be recoverable by re-reading the input with the "
            "identifier the ledger names. Do not use the filtered object until this reconciles.")


# ------------------------------------------------------------------------------------------
# the approval


@dataclass(frozen=True)
class Approval:
    """The two things the gate compares, kept apart on purpose.

    `user_verbatim` is what the operator is asserting they said about THIS action; `approvals` is
    the record of what was actually recorded, keyed by action text. The gate's whole function is
    to compare the two, so they must arrive from separate sources. Deriving one from the other -
    `user_verbatim = approvals[action]` - makes the comparison self-satisfying and turns the only
    gate in this pipeline that stops a removal into a no-op that still reads like a check in
    every log it writes.
    """

    user_verbatim: str
    approvals: dict = field(default_factory=dict)


def _as_mapping(value, action: str, refusal):
    """`dict(value)`, or the refusal - never the TypeError/ValueError `dict()` raises itself.

    `dict("CONFIRM")` raises ValueError and `dict(7)` raises TypeError, and both are the ordinary
    way this function is reached: someone passes the approved words where the record of approvals
    belongs. That is a rejected approval, so it leaves as one.
    """
    try:
        return dict(value)
    except (TypeError, ValueError) as exc:
        raise refusal(
            f"the recorded approvals are not a mapping - refused: {action!r}\n"
            f"  Expected {{action text: verbatim words}}, got {type(value).__name__} "
            f"({value!r}), which {type(exc).__name__} says cannot be read as one. The approvals "
            f"are the record the operator's words are CHECKED AGAINST; a value that is not a "
            f"mapping cannot hold an entry for this action.") from exc


def coerce_approval(approval, action: str) -> Approval:
    """Accept an `Approval`, or a mapping carrying both of its fields, and nothing else.

    Raises `ApplyRefusal` - not `TaskFailure` - for everything it rejects: an approval that
    cannot be established is a consent problem, and it must reach the caller as the same class
    of outcome as the gate's own refusals rather than as something to debug. That promise is kept
    by reading the approval inside a handler, because the reading itself can fail: `dict(...)` on
    a non-mapping raises TypeError or ValueError, `"user_verbatim" in approval` raises TypeError
    for an object whose `__contains__` does not take a string, and a property can raise anything
    at all. Every one of those escaped this function as a bare exception, which is a caller's
    `except ApplyRefusal` not firing on an approval that was in fact rejected.

    The one exception, deliberately: `apply_module()` raises `TaskFailure` when the approval gate
    itself is missing or unloadable. That is not a rejection of this approval - there is no
    ApplyRefusal class to raise it with - and it must not be reported as one.
    """
    refusal = apply_module().ApplyRefusal
    if approval is None:
        raise refusal(
            f"no approval supplied - refused: {action!r}\n"
            f"  A removal requires the operator's own words recorded against this exact action. "
            f"There is no force flag and no default approval.")
    unreadable = (
        f"the approval is not in a form the gate can check - refused: {action!r}\n"
        f"  Pass Approval(user_verbatim=..., approvals=...): the words being offered for "
        f"this action, and the recorded approvals they are checked against, from separate "
        f"sources. A bare string or a bare {{action: words}} mapping would have to serve as "
        f"both sides of the comparison, which always matches - the check would pass on every "
        f"input including no approval at all.")
    try:
        if isinstance(approval, Approval):
            appr = approval
        elif hasattr(approval, "user_verbatim") and hasattr(approval, "approvals"):
            appr = Approval(user_verbatim=approval.user_verbatim,
                            approvals=_as_mapping(approval.approvals, action, refusal))
        elif hasattr(approval, "get") and "user_verbatim" in approval and "approvals" in approval:
            appr = Approval(user_verbatim=approval["user_verbatim"],
                            approvals=_as_mapping(approval["approvals"], action, refusal))
        else:
            raise refusal(unreadable)
    except refusal:
        raise
    except Exception as exc:  # noqa: BLE001 - re-raised as the refusal this function promises
        raise refusal(
            f"{unreadable}\n"
            f"  Reading the approval that was supplied ({type(approval).__name__}) raised "
            f"{type(exc).__name__}: {exc}") from exc
    if not hasattr(appr.approvals, "get"):
        raise refusal(
            f"the recorded approvals are not a mapping - refused: {action!r}\n"
            f"  Expected {{action text: verbatim words}}, got {type(appr.approvals).__name__}.")
    # No check is made that the two sides are distinct OBJECTS. It is tempting - passing
    # `approvals[action]` back as `user_verbatim` does defeat the comparison - but Python interns
    # short string literals, so `{action: "CONFIRM"}` with a separately written "CONFIRM" is the
    # same object, and the check would refuse the ordinary correct call. Keeping the two fields
    # separate in the signature is the part that can be enforced without firing on correct use.
    return appr


# ------------------------------------------------------------------------------------------
# the object


def h5ad_n_obs(path) -> int:
    """Number of observations in an .h5ad, read from the file structure with h5py.

    Two layouts are handled because both are in circulation: modern files store `obs` as a group
    whose `_index` attribute names the index dataset, and files written by anndata before 0.7
    store `obs` as a single compound dataset. An unrecognised layout raises rather than returning
    a count, because a wrong observation count here would confirm the write and hide a truncation.
    """
    import h5py

    path = Path(path)
    with h5py.File(str(path), "r") as f:
        if "obs" not in f:
            raise TaskFailure(f"{path} has no /obs; it is not an .h5ad this adapter can verify.")
        obs = f["obs"]
        if isinstance(obs, h5py.Dataset):  # anndata < 0.7 compound layout
            return int(obs.shape[0])
        key = obs.attrs.get("_index")
        if key is None:
            raise TaskFailure(
                f"{path} has an /obs group with no `_index` attribute, so which dataset holds "
                f"the observation names cannot be established. Guessing one risks confirming a "
                f"truncated write.")
        if isinstance(key, bytes):
            key = key.decode("utf-8")
        key = str(key)
        if key not in obs:
            raise TaskFailure(
                f"{path} declares its observation index as {key!r} but /obs holds no such "
                f"dataset.")
        return int(obs[key].shape[0])


def written_n_obs(path) -> int:
    """How many observations the written object actually holds.

    h5py first because it is cheap and reads no matrix; anndata's own reader as a fallback,
    because it necessarily understands any layout anndata just wrote and a verification that
    fires on a correct write is a verification somebody will delete. Both failing is reported
    with both reasons - never as a count.
    """
    try:
        return h5ad_n_obs(path)
    except Exception as first:  # noqa: BLE001 - re-raised below with the second attempt's reason
        try:
            import anndata

            ad = anndata.read_h5ad(str(path), backed="r")
            try:
                return int(ad.n_obs)
            finally:
                fh = getattr(ad, "file", None)
                if fh is not None:
                    try:
                        fh.close()
                    except Exception:  # noqa: BLE001 - closing a read handle is best effort
                        pass
        except Exception as second:  # noqa: BLE001
            raise TaskFailure(
                f"the filtered object at {path} was written but its observation count could not "
                f"be read back.\n"
                f"  h5py: {first!r}\n"
                f"  anndata: {second!r}\n"
                f"  The file exists and the ledger beside it is complete; check the object by "
                f"hand before using it. The count is reported as unreadable rather than assumed."
            ) from second


def _decode_h5_strings(values) -> list:
    """An HDF5 string dataset as a list of `str`. bytes on some builds, str on others."""
    out = []
    for v in values:
        out.append(v.decode("utf-8", "replace") if isinstance(v, (bytes, bytearray)) else str(v))
    return out


def _h5py_obs_identifiers(path, obs_id_key=None) -> list:
    """The observation identifiers in a written .h5ad, read with h5py and no matrix loaded.

    Handles the modern group layout only - `/obs` a group whose `_index` attribute names the index
    dataset, and a column stored either as a plain dataset or as a categorical `categories`/`codes`
    pair, which is what anndata's writer produces for a string column. Anything else raises, and
    `h5ad_obs_identifiers()` falls back to anndata's own reader rather than guessing.

    A categorical code of -1 is a missing value and stops this: an observation whose identifier is
    absent cannot be compared with the one it was supposed to be written under, and reading the
    absence as a label would make the comparison pass on a file that lost an identifier.
    """
    import h5py

    path = Path(path)
    with h5py.File(str(path), "r") as f:
        if "obs" not in f:
            raise TaskFailure(f"{path} has no /obs; it is not an .h5ad this adapter can verify.")
        obs = f["obs"]
        if not isinstance(obs, h5py.Group):
            raise TaskFailure(f"{path}: /obs is a {type(obs).__name__}, not a group; this adapter "
                              f"reads identifiers only from the group layout.")
        if obs_id_key is None:
            key = obs.attrs.get("_index")
            if key is None:
                raise TaskFailure(
                    f"{path} has an /obs group with no `_index` attribute, so which dataset holds "
                    f"the observation names cannot be established.")
            key = key.decode("utf-8") if isinstance(key, bytes) else str(key)
        else:
            key = str(obs_id_key)
        if key not in obs:
            raise TaskFailure(f"{path}: /obs holds no {key!r}; present: "
                              f"{', '.join(sorted(obs.keys()))}")
        node = obs[key]
        if isinstance(node, h5py.Group):
            if "categories" not in node or "codes" not in node:
                raise TaskFailure(f"{path}: /obs/{key} is a group without categories and codes; "
                                  f"this adapter cannot read it as identifiers.")
            cats = _decode_h5_strings(node["categories"][:])
            codes = [int(c) for c in node["codes"][:]]
            out = []
            for i, c in enumerate(codes):
                if c < 0 or c >= len(cats):
                    raise TaskFailure(
                        f"{path}: observation {i} has no value in /obs/{key} (category code "
                        f"{c}).\n"
                        f"  A missing identifier is not a label. The written object cannot be "
                        f"checked against the identifiers this removal kept.")
                out.append(cats[c])
            return out
        return _decode_h5_strings(node[:])


def h5ad_obs_identifiers(path, obs_id_key=None) -> list:
    """The observation identifiers a written .h5ad actually holds, in stored order.

    `obs_id_key=None` reads the observation index; a key reads that column of `.obs`. This is the
    read half of the round trip: what went to the writer is compared against what came back off
    disk, so the lossless label cast the write needs is a checked claim rather than an asserted
    one.

    h5py first because it reads no matrix, anndata second because it necessarily understands any
    layout anndata just wrote - a verification that fires on a correct write is a verification
    somebody will delete. Both failing is reported with both reasons and never as a list.
    """
    path = Path(path)
    try:
        return _h5py_obs_identifiers(path, obs_id_key)
    except Exception as first:  # noqa: BLE001 - re-raised below with the second attempt's reason
        try:
            import anndata

            ad = anndata.read_h5ad(str(path), backed="r")
            try:
                if obs_id_key is None:
                    return [str(x) for x in ad.obs_names]
                if obs_id_key not in ad.obs.columns:
                    raise TaskFailure(
                        f"{path}: .obs has no column {obs_id_key!r}; it holds "
                        f"{list(ad.obs.columns)[:10]}...")
                col = ad.obs[obs_id_key]
                if bool(col.isna().any()):
                    raise TaskFailure(
                        f"{path}: .obs[{obs_id_key!r}] contains missing values, so the written "
                        f"object cannot be checked against the identifiers this removal kept.")
                return [str(x) for x in col]
            finally:
                fh = getattr(ad, "file", None)
                if fh is not None:
                    try:
                        fh.close()
                    except Exception:  # noqa: BLE001 - closing a read handle is best effort
                        pass
        except Exception as second:  # noqa: BLE001
            raise TaskFailure(
                f"the identifiers of the object at {path} could not be read back.\n"
                f"  h5py: {first!r}\n"
                f"  anndata: {second!r}\n"
                f"  The write is therefore unverified: nothing here can say the object names the "
                f"observations this removal kept. Do not use it until it has been checked by "
                f"hand.") from second


def _require_same_identifiers(got, expected, path, what: str, context: str = "") -> None:
    """Two lists of identifiers are equal, in order, or the run stops naming the first difference.

    Compared element by element rather than as sets or by length: a written object holding the
    right identifiers in the wrong order is a different object, and every count taken from it
    still agrees with itself.
    """
    if len(got) != len(expected):
        raise TaskFailure(
            f"{path} holds {len(got):,} identifiers in {what} where {len(expected):,} were "
            f"written.{context}")
    for i, (a, b) in enumerate(zip(got, expected)):
        if a != b:
            raise TaskFailure(
                f"{path} does not name the observations that were written to it: {what} entry "
                f"{i} read back as {a!r} where {b!r} was written.\n"
                f"  The identifiers were captured before the object was prepared for writing and "
                f"compared against the file afterwards, so this says the write changed them. The "
                f"ledger names what left by identifier; an object whose identifiers differ from "
                f"it cannot be reconciled with it.{context}")


def _module_version(mod, name: str) -> str:
    """The version of a library this run actually imported.

    Distribution metadata is asked first because `anndata.__version__` is deprecated and warns,
    and will eventually be absent; the module attribute is the fallback for anything installed
    without metadata. A library that answers neither is recorded as having answered neither -
    never as an empty string, which would read in a report as a version nobody looked up.
    """
    from importlib.metadata import PackageNotFoundError, version as _dist_version

    try:
        return str(_dist_version(name))
    except (PackageNotFoundError, ValueError, OSError):
        pass
    v = getattr(mod, "__version__", None)
    if v is None:
        return f"{name} imported from {getattr(mod, '__file__', 'an unknown path')}, " \
               f"but reports no version"
    return str(v)


def _resolve_adata(adata, obs_id_key):
    """(AnnData, source path or None, identifiers). Accepts an object or a path to one."""
    src = None
    if isinstance(adata, (str, Path)):
        src = Path(adata)
        if not src.exists():
            raise TaskFailure(f"no such input matrix: {src}")
        import anndata

        obj = anndata.read_h5ad(str(src))
    else:
        obj = adata
        for attr in ("obs", "obs_names", "n_obs"):
            if not hasattr(obj, attr):
                raise TaskFailure(
                    f"`adata` is a {type(obj).__name__}, which is neither a path to an .h5ad nor "
                    f"an AnnData (no .{attr}).")
    if obs_id_key is None:
        raw = list(obj.obs_names)
        label = "obs_names"
    else:
        if obs_id_key not in obj.obs.columns:
            raise TaskFailure(
                f"obs_id_key={obs_id_key!r} is not a column of .obs; it holds "
                f"{list(obj.obs.columns)[:10]}...")
        raw = list(obj.obs[obs_id_key])
        label = f"obs[{obs_id_key!r}]"
    ids = normalise_identifiers(raw, label)
    if len(ids) != int(obj.n_obs):
        raise TaskFailure(
            f"{label} yielded {len(ids):,} identifiers for {int(obj.n_obs):,} observations.")
    return obj, src, ids


def apply_filter(adata, keep_mask, criteria, out_h5ad, ledger_path, action, approval,
                 obs_id_key=None, allow_empty=False, compression="gzip") -> dict:
    """Remove what was approved, record what left, and prove afterwards that it can be recovered.

    `keep_mask` is True where an observation STAYS. `criteria` maps a criterion name to a mask
    that is True where that criterion REMOVES the observation. The two conventions are opposite
    and both are checked against each other by identity before anything happens, because an
    inverted criteria mask produces a plausible object and a ledger describing the complement of
    the removal.

    Nothing is written until `modules/07_apply/apply.py::apply_removal()` has accepted the
    approval, and its refusal propagates unchanged. The gate is handed the mask-derived count and
    the record built from the criteria, which are the two independent routes it exists to
    compare; handing it one number twice would leave that comparison passing on anything.

    Both declared outputs are produced fresh, never inherited. The filtered object's path is
    deleted and observed absent before anndata is called, and the object it is handed carries a
    token invented after that moment, which is read back off disk afterwards; the ledger is
    written to a temporary file, renamed into place and read back. Either way a previous run's
    artifact cannot satisfy the check that this run wrote one - which it would, silently, if the
    writer returned without writing and only existence were asked, and which it also would under
    the mtime comparison this replaced, since a restored file back-dated by a second lands inside
    any tolerance a clock allowance makes necessary. The deletion happens after the gate and after
    the ledger, so the ordering guarantee is unchanged: nothing on disk is touched until the
    removal has been approved and recorded.

    The object's own observation identifiers are then read back out of the .h5ad and compared, in
    order, with the ones captured before it was prepared for writing - which is what makes the
    lossless label cast that write needs a checked claim. `.obs` and `.var` are re-backed by object
    arrays first, because pandas gives string labels a dtype anndata refuses to write and the
    alternative, `anndata.settings.allow_write_nullable_strings`, is a global that would change how
    every other object in the process is stored.

    Returns `{"outputs": [...], "metrics": {...}, "versions": {...}}`. The counts - `n_in`,
    `n_kept`, `n_removed`, the per-criterion totals and the per-criterion SOLE totals - are in
    `metrics`, and `outputs` lists only files whose existence was checked after the write.
    """
    gate = apply_module()
    # Loaded here, with the gate, rather than at the point of use: it supplies the cast without
    # which the object cannot be written and the token by which the write is verified, so a tree
    # missing it must stop before the approval is spent and the ledger is on disk.
    mx = matrix_adapter()
    out_h5ad = Path(out_h5ad)
    ledger_path = Path(ledger_path)
    if not isinstance(action, str) or not action.strip():
        raise TaskFailure(
            "`action` must be the exact text the approval was recorded against; a removal with "
            "no action text cannot be matched to any approval.")
    if out_h5ad.resolve() == ledger_path.resolve():
        raise TaskFailure(
            f"the filtered object and the ledger are the same path ({out_h5ad}); one would "
            f"overwrite the other and the removal would end up with no record.")

    appr = coerce_approval(approval, action)
    obj, src, ids = _resolve_adata(adata, obs_id_key)
    if src is not None and src.resolve() == out_h5ad.resolve():
        raise TaskFailure(
            f"the filtered object would be written over its own input ({src}).\n"
            f"  The ledger stores identifiers, not data: recovering a removed observation means "
            f"re-reading the input with the identifier the ledger names. Overwriting the input "
            f"deletes the only copy of what left and makes the ledger a list of names that "
            f"resolve to nothing. Write to a new path.")

    n_in = len(ids)
    rec = reconcile_masks(ids, keep_mask, criteria)
    names, masks, keep = rec["criteria_names"], rec["criteria_masks"], rec["keep"]
    n_removed_mask = sum(1 for k in keep if not k)
    n_kept_mask = n_in - n_removed_mask
    if n_kept_mask == 0 and not allow_empty:
        raise TaskFailure(
            f"the keep-mask keeps 0 of {n_in:,} observations.\n"
            f"  An empty matrix is written and read as a successful result, and every count "
            f"downstream of it is zero for a reason nobody can see. If removing everything is "
            f"genuinely the finding, pass allow_empty=True and say so in the run log; if it is "
            f"not, check whether the mask is inverted.")

    rows = _rows_from_masks(ids, names, masks)
    record = gate.build_removal_record(ids, masks)  # the gate's own route to the same removal
    if record.n_removed != len(rows) or record.n_in != n_in:
        raise TaskFailure(
            f"two independent constructions of the removal disagree: this adapter's ledger names "
            f"{len(rows):,} of {n_in:,} observations, the gate's record names "
            f"{record.n_removed:,} of {record.n_in:,}. Neither is used until they agree.")

    # ---- the gate. Nothing above this line has written to disk; nothing below it runs if the
    # approval is absent, empty, or recorded against different words.
    kept_by_gate = gate.apply_removal(n_in, n_removed_mask, action, appr.user_verbatim,
                                      appr.approvals, record=record, record_path=None)
    if int(kept_by_gate) != n_kept_mask:
        raise TaskFailure(
            f"the gate returned {int(kept_by_gate):,} kept where the mask keeps {n_kept_mask:,}.")

    # ---- ledger first: a ledger for a removal that did not happen is recoverable, a removal
    # with no ledger is not.
    written_ledger = write_ledger(rows, ledger_path, columns=ledger_columns(names))

    # Everything from here to the write happens AFTER the ledger, so every failure in it has to
    # carry the same explanation: the ledger is already on disk and nothing has been lost. This
    # sentence is built once and reused rather than being attached to one of the failures.
    after_ledger = (
        f"\n  The approval was accepted and the ledger at {written_ledger} was written first, so "
        f"what would have left is still named and nothing has been lost. The input is untouched. "
        f"Fix this and re-run the step; do not hand-edit the ledger to match a partial object.")

    try:
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        # This mkdir used to sit outside any handler, one line below the ledger write, so the
        # failure a full or read-only output directory produces - the commonest failure at this
        # point in the step - escaped as a bare OSError and the explanation below it could not be
        # reached from the line that raises it.
        raise TaskFailure(
            f"the output directory {out_h5ad.parent} could not be created: "
            f"{type(exc).__name__}: {exc}{after_ledger}") from exc

    import numpy as np

    keep_arr = np.asarray(keep, dtype=bool)
    filtered = obj[keep_arr].copy()
    if int(filtered.n_obs) != n_kept_mask:
        raise TaskFailure(
            f"subsetting produced {int(filtered.n_obs):,} observations where the mask keeps "
            f"{n_kept_mask:,}; the mask and the object are not aligned.")
    kept_ids = [ids[i] for i in range(n_in) if keep[i]]

    # Captured from the object as it stands, BEFORE it is prepared for writing. These are what the
    # read-back is compared against, so the preparation below is a checked claim and not an
    # asserted one: if the cast altered an identifier, or the writer did, the comparison says so.
    expected_index = [str(x) for x in filtered.obs_names]
    expected_key = None if obs_id_key is None else [str(x) for x in filtered.obs[obs_id_key]]

    try:
        # pandas gives a plain string index the `str`/StringDtype backing by default and anndata
        # declines to write it, so on a current stack the deliverable could not be written at all.
        # Re-backing the labels with object arrays is lossless - the same `str` objects, a
        # different array behind them - and is preferred over
        # `anndata.settings.allow_write_nullable_strings`, which is global: flipping it would
        # change how every other object written anywhere in this process is stored, and would
        # produce files anndata < 0.11 cannot read.
        new_obs, new_var = mx.object_backed_labels(filtered, where="apply_filter")
    except TaskFailure as exc:
        raise TaskFailure(f"{exc}{after_ledger}") from exc
    if new_obs is not None:
        filtered.obs = new_obs
    if new_var is not None:
        filtered.var = new_var

    started = _clear_previous_output(out_h5ad, after_ledger)
    # Invented after the path was observed absent, so nothing already on disk can be carrying it.
    # `filtered` is this call's own copy, so stamping it changes nothing the caller handed in.
    token = mx.new_write_token()
    filtered.uns[mx.WRITE_TOKEN_KEY] = token
    try:
        filtered.write_h5ad(str(out_h5ad), compression=compression)
    except Exception as exc:  # noqa: BLE001 - reported with the state it leaves behind
        hint = ""
        if "allow_write_nullable_strings" in str(exc):
            # The lossless cast above covers the labels and the string columns anndata declines.
            # Reaching here anyway means something it does not cover, and the remedy is still not
            # to flip the global setting on the caller's behalf.
            hint = ("\n  Cause: a column or index anndata declines to write as a nullable/arrow "
                    "string array survived `adapters/matrix.py::object_backed_labels()`, which "
                    "this adapter applies to .obs and .var before writing. This adapter does not "
                    "flip `anndata.settings.allow_write_nullable_strings` for you - it is a "
                    "global setting that changes how every object in the run is written and "
                    "produces files anndata < 0.11 cannot read. Cast the offending field to "
                    "object dtype upstream, where the decision is recorded, or set the setting "
                    "deliberately in the orchestrator and say that you did.")
        raise TaskFailure(
            f"the filtered object could not be written to {out_h5ad}: {type(exc).__name__}: "
            f"{exc}{after_ledger}{hint}") from exc
    _confirm_fresh_output(out_h5ad, started, token, after_ledger)
    n_written = written_n_obs(out_h5ad)
    if n_written != n_kept_mask:
        raise TaskFailure(
            f"{out_h5ad} holds {n_written:,} observations but {n_kept_mask:,} were kept.\n"
            f"  The ledger at {written_ledger} describes the intended removal, so what left is "
            f"still named; the object does not match it and must not be used.")

    # ---- the round trip. The count above says how many observations survived; these say WHICH,
    # by reading them back out of the file and comparing them with what went to the writer. A
    # count cannot tell a correct removal from one that kept the right number of the wrong rows.
    written_index = h5ad_obs_identifiers(out_h5ad)
    _require_same_identifiers(written_index, expected_index, out_h5ad, ".obs_names", after_ledger)
    if obs_id_key is None:
        written_raw = written_index
    else:
        written_raw = h5ad_obs_identifiers(out_h5ad, obs_id_key)
        _require_same_identifiers(written_raw, expected_key, out_h5ad, f".obs[{obs_id_key!r}]",
                                  after_ledger)
    written_ids = normalise_identifiers(written_raw, "written identifier")
    if written_ids != kept_ids:
        differ = next((i for i, (a, b) in enumerate(zip(written_ids, kept_ids)) if a != b), None)
        where = ("" if differ is None else
                 f" First difference at position {differ}: the file names "
                 f"{written_ids[differ]!r}, the removal kept {kept_ids[differ]!r}.")
        raise TaskFailure(
            f"{out_h5ad} does not name the observations this removal kept: {len(written_ids):,} "
            f"identifiers on disk against {len(kept_ids):,} kept.{where}\n"
            f"  The ledger at {written_ledger} names what left, by identifier, and it cannot be "
            f"reconciled with an object naming anything else.{after_ledger}")

    verify_recoverable(written_ledger, ids, kept_ids)

    summary = criterion_summary(rows, names)
    import anndata

    versions = {"anndata": _module_version(anndata, "anndata"),
                "numpy": _module_version(np, "numpy")}
    try:
        import h5py

        versions["h5py"] = _module_version(h5py, "h5py")
    except ImportError:
        versions["h5py"] = "not invoked"
    versions["scqc.modules.07_apply"] = _gate_fingerprint()

    return {
        "outputs": [out_h5ad, written_ledger],
        "metrics": {
            "action": action,
            "n_in": n_in,
            "n_kept": n_kept_mask,
            "n_removed": n_removed_mask,
            "pct_removed": round(100.0 * n_removed_mask / n_in, 4) if n_in else None,
            "n_criteria": len(names),
            "criteria": list(names),
            "by_criterion": summary["by_criterion"],
            "sole_criterion": summary["sole_criterion"],
            "by_n_criteria": summary["by_n_criteria"],
            "ledger": str(written_ledger),
            "ledger_bytes": written_ledger.stat().st_size,
            "input_path": str(src) if src is not None else "in-memory AnnData",
            "obs_id_key": "obs_names" if obs_id_key is None else obs_id_key,
            "recoverability_verified": True,
            # Both read back off disk, not asserted: the identifiers in the written object were
            # compared one by one with the ones handed to the writer, and the token stored in it
            # is this call's.
            "identifiers_verified": True,
            "write_token": token,
        },
        "versions": versions,
    }


def _gate_fingerprint() -> str:
    """Hash of the approval-gate source that authorised this removal.

    The gate has no version string, and which text a removal was approved under is exactly the
    thing a later reader needs. Recorded as an observation of the file on disk rather than as a
    claim about which release it came from.
    """
    from engine.provenance import file_hash

    return file_hash(repo_root() / "modules" / "07_apply" / "apply.py", cap_mb=None)


def run_apply_filter(adata, keep_mask, criteria, out_h5ad, ledger_path, action, approval,
                     obs_id_key=None, allow_empty=False, compression="gzip",
                     executor=None) -> dict:
    """Task-shaped wrapper over `apply_filter()`, for an orchestrator that calls every adapter
    the same way.

    The executor is accepted and deliberately not used: this is the one step that removes
    observations, and the approval check, the ledger and the write must happen in one process.
    Shipping them to a scheduler puts a process boundary between the gate and the write, where
    the refusal becomes an exit code that a wrapper script can drop and the masks that were
    approved are no longer demonstrably the masks that were applied. Which executor was offered
    is recorded, because a run that expected to be scheduled and was not is worth seeing.
    """
    out = apply_filter(adata, keep_mask, criteria, out_h5ad, ledger_path, action, approval,
                       obs_id_key=obs_id_key, allow_empty=allow_empty, compression=compression)
    out["metrics"]["executor"] = "in-process (removal is never delegated)"
    out["metrics"]["executor_offered"] = getattr(executor, "name", type(executor).__name__
                                                 if executor is not None else "none")
    return out
