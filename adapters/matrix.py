# Matrix input, matrix output, and the summary statistics the gates downstream consume.
# This adapter reads matrices and writes objects. It removes no barcode, no gene and no cluster,
# it subsets nothing, and it renames no identifier. Step 7 is the only code path in this pipeline
# permitted to drop an observation, and it needs a recorded approval to do it.
"""Step 0's measurement layer: produce the numbers `lib/verify_raw.py` returns a verdict on.

WHY THIS EXISTS

`verify_raw.verify()` decides whether a matrix is raw from `n_barcodes`, `n_genes`, `min_counts`,
`max_counts`, `p98_counts` and `integer_counts`, and `modules/00_ingest/ingest.plan_one()` obtains
them by calling a `stats_fn` its caller supplies. No code in this repository could produce that
function's return value, so the verdict at step 0 rested on statistics computed somewhere else and
typed in by hand. A verdict is worth exactly what the measurement under it is worth, and a
hand-typed measurement carries no provenance at all: it cannot be recomputed, it cannot be
attributed to a file, and nothing detects it drifting away from the matrix it claims to describe.

`summary_stats()` is that function. The keys it returns are the keyword arguments `verify()`
accepts, plus `n_nonzero_barcodes`.

UNKNOWN IS NOT A VALUE, AND THIS IS WHERE IT WOULD ENTER

Every number `verify()` reads is a float, and every one of its tests is a comparison. A NaN
compares False against everything, so a single NaN reaching `verify()` turns
`min_counts >= 100` into "this matrix has raw droplets" - the exact shape of the bug
docs/PRINCIPLES.md section 4 records three instances of, and the one that is hardest to see
because the resulting verdict reads as a clean pass. So:

  * per-barcode totals are checked with `isfinite` and a non-finite total is a TaskFailure,
    never a value that flows on;
  * `integer_check()` answers False for NaN, for infinity and for negative values, because none
    of them is a count;
  * the JSON written by the subprocess entry point is dumped with `allow_nan=False`, so a NaN
    fails at the point of writing rather than being read back as a number;
  * `parse_stats_json()` rejects the JSON `NaN` and `Infinity` literals rather than turning them
    into floats, and refuses a payload with a missing or null required key rather than defaulting
    it.

NEVER DENSIFY

A raw droplet matrix is hundreds of thousands of barcodes by tens of thousands of genes, and its
dense form does not fit in memory on any machine this pipeline runs on. Per-barcode totals come
from a sparse row sum; the percentile is taken over that 1-D vector of totals, which is one
float64 per barcode and small. `integer_check()` walks the stored values in bounded blocks - for
a sparse matrix the stored values are the non-zeros, and the implicit zeros need no examination
because zero is integral and non-negative.

`p98_counts` IS THE 98TH PERCENTILE OF PER-BARCODE TOTAL COUNTS OVER EVERY BARCODE IN THE FILE

Not over the non-empty ones, and not over a called subset. The consequence is worth stating
because it looks like a defect and is not: on a genuinely raw droplet matrix the 98th percentile
sits deep among the empties, so `verify()`'s ceiling test - max sitting on its own p98 - cannot
fire there. That test exists for a delivered matrix that has already been through cell calling,
where p98 is computed over cells and a ceiling is visible. A matrix that has been both capped and
cell-called fails both of `verify()`'s tests, which is the case the test was written for.

MEASURE THE STRUCTURE, DO NOT ASSUME IT

The HDF5 reader opens the file and inspects the layout before reading it, and builds the AnnData
itself. Delegating to a reader written for one producer's layout means a file from another
producer either fails with an error about a missing dataset or, worse, is read under an
assumption that happens to hold for the shape and not for the contents. CellRanger v2, CellRanger
v3 and CellBender all write a CSC block of features x barcodes with the feature and barcode
tables beside it, under three different group names; that structure is what is checked for, and
an unrecognised layout is refused with the group names the file actually contains.

A DECLARED OUTPUT IS DELETED BEFORE IT IS WRITTEN, AND CARRIES A TOKEN THIS CALL INVENTED

`write_h5ad()` and `run_summary_stats()` both remove their output path, confirm it is gone, and
only then invoke the writer. Asking only whether the output exists is satisfied by the previous
run's file, so a tool that exits 0 having written nothing is recorded as a success and the earlier
parameters' numbers are reported under the new ones - invisibly, because a stale .h5ad or
statistics JSON opens, parses and reads exactly like a fresh one. Deleting first is what makes the
existence check mean something. It removes an artifact of an earlier run of the same step and
nothing else; no input is ever a candidate.

The second statement - that the bytes now at that path are this call's - is made with a token and
NOT with an mtime. An mtime was the previous spelling and it was blind by construction: clocks and
filesystems disagree, so the comparison needed a tolerance window, and a window is precisely where
a restored artifact lands - a copy of a previous object back-dated by one second was accepted,
while the same construction at two hours was refused. `os.utime` will set a file's mtime to any
value at all, so no width of window fixes that. Instead `new_write_token()` invents a uuid4 AFTER
the path has been observed absent, the writer stores it INSIDE the artifact - under `uns` for an
.h5ad, as a payload field for the statistics JSON - and it is read back off disk and compared.
A file that does not carry this call's token was not produced by this call, however new it looks.

WHAT THIS ADAPTER DOES NOT DO

It does not make variable names unique, does not drop duplicate features, does not filter, does
not normalise, does not transpose anything it has not identified by comparing the stored shape
against the lengths of the barcode and feature tables, and does not choose between two matrices
found inside one archive. Each of those is a decision, decisions belong to the orchestrator and
the modules, and one of them - choosing between a raw and a filtered matrix - is precisely the
choice step 0 exists to make explicit.

The single transformation it does apply is on the way out and is lossless: `write_h5ad()` re-backs
`pandas` nullable/arrow string columns and indices with object arrays, because anndata refuses to
write them and pandas 3.0 makes them the default - without it, `write_h5ad(read_matrix(x))` raises
on the stack this repository is pinned against. Same `str` objects, different array behind them;
`object_backed_labels()` states the reasoning, and a column that actually contains missing values
is refused rather than given an invented label.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import tarfile
import time
import uuid
from pathlib import Path
from typing import Any, Iterator, Optional, Sequence

try:  # The repo root is not necessarily on sys.path: `pip install -e .` exposes only scqc_cli,
    from engine.task import TaskFailure  # and this file is also loaded by path by the subprocess
    from engine.provenance import NOT_INVOKED  # entry point below.
except ImportError:  # pragma: no cover - exercised only outside an in-tree import
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from engine.task import TaskFailure
    from engine.provenance import NOT_INVOKED

#: The keyword arguments `lib/verify_raw.py::verify` accepts, in its own order.
VERIFY_KEYS = ("name", "n_barcodes", "n_genes", "min_counts", "max_counts", "p98_counts",
               "expected_genes", "integer_counts")
#: Measured here and deliberately not passed to `verify`: it is what tells a reader whether a
#: minimum of zero means "this file contains genuine empty droplets" or "a barcode list was
#: padded", and it is the denominator figure F1 needs.
EXTRA_KEYS = ("n_nonzero_barcodes",)
#: Every key `summary_stats` returns. The set is fixed; a key is never omitted when its value is
#: unavailable, because an absent key and a key whose value is None are read differently by
#: everything downstream.
STATS_KEYS = VERIFY_KEYS + EXTRA_KEYS
#: Keys whose value must be a number or a bool. `name` is a string and `expected_genes` may be
#: None - a reference size that was not declared is not a measurement that failed.
STATS_REQUIRED_VALUES = ("n_barcodes", "n_genes", "min_counts", "max_counts", "p98_counts",
                         "integer_counts", "n_nonzero_barcodes")
#: Of those, the ones that count things: whole numbers, never negative.
STATS_COUNT_KEYS = ("n_barcodes", "n_genes", "n_nonzero_barcodes")
#: And the ones that are sums over a row: any finite real. Not required to be positive, because a
#: scaled or centred matrix has negative row sums and whether that disqualifies it as raw input is
#: `verify()`'s verdict to render, not a parser's.
STATS_REAL_KEYS = ("min_counts", "max_counts", "p98_counts")

#: The payload `main()` writes and `parse_stats_json()` reads. Versioned so that a file written by
#: an older checkout is refused rather than parsed under the wrong assumptions.
STATS_SCHEMA = "scqc.matrix.stats/1"

#: Bounds on how much is held at once while scanning values. Neither is a threshold on the data.
_VALUE_BLOCK = 1 << 22          # stored values examined per block by integer_check
_ROW_BLOCK = 8192               # rows sliced per block from a backed or on-disk matrix

#: Where a write records, INSIDE the artifact it produced, WHICH CALL produced it. An mtime cannot
#: carry that: `os.utime` sets it to any value, a restored file lands wherever the restorer chose,
#: and the clock allowance a cross-host comparison needs is itself a window a restored file fits
#: inside. A uuid4 invented after the path was observed absent exists nowhere else until this
#: call's writer stores it, so an artifact that does not carry it was not written by this call.
#: In an .h5ad it is an entry in `uns`; in the statistics JSON it is a top-level field.
WRITE_TOKEN_KEY = "scqc_write_token"

_TAR_SUFFIXES = (".tar", ".tar.gz", ".tgz", ".tar.bz2", ".tbz2", ".tbz", ".tar.xz", ".txz")
_MTX_NAMES = ("matrix.mtx.gz", "matrix.mtx")
_BARCODE_NAMES = ("barcodes.tsv.gz", "barcodes.tsv")
_FEATURE_NAMES = ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
#: HDF5 datasets that together make a CSC block, whatever group they are found under.
_CSC_PARTS = ("data", "indices", "indptr", "shape")


# ---------------------------------------------------------------------------------------------
# small pure helpers


#: The missing-value sentinels the scientific stack delivers, matched by TYPE NAME so that
#: `_is_missing()` recognises them without importing anything - this module is imported by
#: `scqc_cli.py` on interpreters that have neither pandas nor numpy. `pandas.NA` (`NAType`),
#: `pandas.NaT` (`NaTType`) and `numpy.ma.masked` (`MaskedConstant`) are three distinct objects
#: from two libraries; none is None, none is a float, and each arrives here from an ordinary
#: table - a nullable column, a datetime column, any masked array indexed at a masked position.
_MISSING_TYPE_NAMES = frozenset(("NAType", "NaTType", "MaskedConstant"))
_MISSING_LIBRARIES = frozenset(("numpy", "pandas"))


def _library_missing(value: Any) -> bool:
    """The missing test for a value out of numpy or pandas. Lazy imports, absorbed if absent.

    `numpy.isnan` is asked rather than `value != value` because the numpy scalars that are not
    float subclasses - `float32`, `datetime64('NaT')` - answer the first correctly and the second
    inconsistently, and because `numpy.ma.masked != numpy.ma.masked` evaluates to `masked`, which
    is FALSY: the self-inequality spelling reports the standard masked sentinel as a measurement.
    """
    try:
        import numpy as np
    except ImportError:  # pragma: no cover - exercised only on a stdlib-only interpreter
        np = None
    if np is not None:
        masked = getattr(getattr(np, "ma", None), "masked", None)
        if masked is not None and value is masked:
            return True
        if isinstance(value, np.generic):
            try:
                return bool(np.isnan(value))
            except (TypeError, ValueError):
                return False  # a numpy string or object scalar: not a number, not a NaN
    try:
        import pandas as pd
    except ImportError:  # pragma: no cover - as above
        return False
    try:
        answer = pd.isna(value)
    except (TypeError, ValueError):
        return False
    if isinstance(answer, bool):
        return answer
    if np is not None and isinstance(answer, np.bool_):
        return bool(answer)
    return False  # an array answer means value was a container, not a missing scalar


def _is_missing(value: Any) -> bool:
    """True when a value is not a measurement. THE missing test for this module - use no other.

    One predicate, because every place that asks "was this measured?" has to agree, and unmeasured
    has more than the two spellings this function used to enumerate: `None`, `float('nan')`,
    `numpy.float64('nan')`, `numpy.float32('nan')`, `pandas.NA`, `pandas.NaT` and
    `numpy.ma.masked` all mean the same thing, only one is None, only some are floats, and none of
    them is equal to itself in the way the naive test assumes.

    Each miss fails in the same direction and that direction is a PASS: every test in
    `verify()` is a comparison, a sentinel that survives this predicate is compared with `>=`, and
    a comparison a sentinel loses reads as a matrix that failed the threshold - which is what
    `verify()` was asked to detect. `numpy.ma.masked` is worse still: `masked != masked` is falsy,
    so the old body answered "measured" for a value that is definitionally absent, and
    `bool(pandas.NA)` raises instead of answering at all.

    Returns a real `bool`, never `numpy.bool_`: a caller writing `x is True` around a numpy
    boolean gets False for a true answer, so this predicate does not hand one out.
    """
    if value is None:
        return True
    if type(value).__name__ in _MISSING_TYPE_NAMES:
        return True
    if isinstance(value, float):
        return bool(value != value)  # float nan, and numpy.float64 nan, a float subclass
    if getattr(type(value), "__module__", "").split(".")[0] in _MISSING_LIBRARIES:
        return _library_missing(value)
    try:
        return bool(value != value)  # any other library's NaN scalar, e.g. Decimal('NaN')
    except Exception:  # noqa: BLE001 - an exotic __ne__ is not evidence of missingness
        return False


def _plain_scalar(value: Any) -> Any:
    """A numpy scalar as its Python value; anything else unchanged.

    `numpy.int64` is not an `int` and `numpy.bool_` is not a `bool`, so the type tests below would
    reject an ordinary measurement taken with numpy. Unwrapping once, here, is what lets those
    tests be strict about the type without firing on correct input.
    """
    if getattr(type(value), "__module__", "").split(".")[0] != "numpy":
        return value
    item = getattr(value, "item", None)
    if not callable(item):
        return value
    try:
        out = item()
    except (TypeError, ValueError):
        return value
    return value if out is value else out


def _require_positive_int_or_none(value: Any, what: str) -> Optional[int]:
    """A declared count is an int above zero, or it was not declared. Nothing in between."""
    if value is None:
        return None
    if _is_missing(value):
        raise TaskFailure(
            f"{what} is missing ({value!r}). A blank cell read from a table arrives as float nan, "
            f"pandas.NA or numpy.ma.masked, and every one of them compares False against every "
            f"threshold - pass None if it was not declared, so it is absent rather than silently "
            f"passing every test.")
    value = _plain_scalar(value)
    if isinstance(value, str) or isinstance(value, (bytes, bytearray)):
        raise TaskFailure(
            f"{what} must be an integer or None, got the {type(value).__name__} {value!r}. Text "
            f"is refused rather than converted: a number that arrived as text arrived from "
            f"somewhere that did not parse it, and the next value from there may not convert.")
    try:
        out = int(value)
        exact = out == float(value)
    except (TypeError, ValueError, OverflowError):
        # OverflowError is the infinity case: JSON's 1e999 parses to float('inf') without ever
        # touching the NaN/Infinity literals, and int(inf) raises a class this used to let
        # escape - past every caller's `except TaskFailure`.
        raise TaskFailure(f"{what} must be a whole number above zero or None, got "
                          f"{value!r}") from None
    if not exact:
        raise TaskFailure(f"{what} must be a whole number, got {value!r}")
    if out <= 0:
        raise TaskFailure(f"{what} must be greater than zero, got {out}")
    return out


def _clear_previous_output(path) -> float:
    """Delete a declared output BEFORE it is written, and return the time the write starts.

    Checking that an output exists after a writer returns proves nothing when the same path
    already held the previous run's file. A tool that exits 0 having written nothing - and every
    tool this pipeline drives can - leaves the earlier artifact standing, it passes the existence
    check, and the earlier parameters' numbers are then reported under the new ones. Nothing
    downstream can see it happen: a stale .h5ad or a stale statistics JSON opens, parses and reads
    exactly like a fresh one.

    So the path is removed first and its absence is confirmed, and the returned timestamp lets the
    caller make the second, independent statement afterwards - that the file now there was
    modified after this write began.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        if p.is_symlink() or p.exists():
            p.unlink()
    except OSError as exc:
        raise TaskFailure(
            f"the previous {p} could not be removed before writing this run's: "
            f"{type(exc).__name__}: {exc}\n"
            f"  It is removed first so that a run which produces nothing cannot leave the earlier "
            f"file behind to be read as this run's result.") from exc
    if p.exists():
        raise TaskFailure(
            f"{p} still exists after being deleted. Another process is writing to it, or the path "
            f"is a directory or a mount point; either way this run's output could not be told "
            f"from the one already there.")
    return time.time()


def _finite_number(value: Any, key: str, where: str, *, whole: bool, non_negative: bool) -> float:
    """One required statistic, as a finite Python float. Anything else names the key and stops.

    Presence and non-missingness are not enough, and this is where that was learned:
    `parse_stats_json` checked only that a required key was there and not None or NaN, so
    `"abc"` and `1e999` reached
    `verify()` untouched. JSON has no infinity literal but `1e999` parses to `float('inf')` without
    ever passing through `parse_constant`, and text passes every missing-value test there is. Both
    then meet `>=` in a gate: the string raises TypeError from inside a comparison three frames
    away from anything that names the file, and the infinity silently satisfies every floor.
    """
    raw = _plain_scalar(value)
    if isinstance(raw, bool):
        raise TaskFailure(
            f"{where}: {key} is the boolean {raw!r}, not a measurement. True compares equal to 1 "
            f"against every threshold in verify(), so a flag that arrived where a count belongs "
            f"would pass or fail on nothing at all.")
    if not isinstance(raw, (int, float)):
        raise TaskFailure(
            f"{where}: {key} is {value!r} ({type(value).__name__}), which is not a number.\n"
            f"  Every test in verify() is a comparison. Text is refused rather than converted: "
            f"'{raw}' >= 100 raises TypeError from inside a gate that cannot say which file it "
            f"was reading, and a value that arrived as text arrived from something that did not "
            f"measure it.")
    out = float(raw)
    if out != out or out in (float("inf"), float("-inf")):
        raise TaskFailure(
            f"{where}: {key} is {value!r}, which is not finite.\n"
            f"  NaN loses every comparison and infinity wins every one of them, so either turns "
            f"verify() into a verdict about nothing. Note that no JSON literal is needed for "
            f"this: 1e999 parses to inf through the ordinary number path.")
    if whole and out != int(out):
        raise TaskFailure(f"{where}: {key} is {value!r}; a count of barcodes or genes is a whole "
                          f"number.")
    if non_negative and out < 0:
        raise TaskFailure(f"{where}: {key} is {value!r}; a count cannot be negative.")
    return out


def require_measured_stats(stats: dict, where: str) -> None:
    """Every statistic `verify()` compares is present, measured, numeric and finite. Or it stops.

    The one place that judgement is made, so that the JSON route and the in-process route cannot
    diverge: `parse_stats_json()` reads a payload written by another interpreter and
    `verify_kwargs()` takes a dict built in this one, and a value that one of them refuses must
    not be a value the other passes through to the same `verify()`.

    Silence is the pass. Every refusal names the offending key, because a report that the
    statistics are unusable without saying which one sends the reader to re-measure all of them.
    """
    covered = set(STATS_COUNT_KEYS) | set(STATS_REAL_KEYS) | {"integer_counts"}
    uncovered = [k for k in STATS_REQUIRED_VALUES if k not in covered]
    if uncovered:
        # A required statistic with no rule would be checked for presence and nothing else, which
        # is the state this function was written to end. It fails loudly here rather than passing
        # the new key through unexamined the way the old code passed every key through.
        raise TaskFailure(
            f"{where}: the required statistics {uncovered} have no validation rule in this "
            f"module. Add them to STATS_COUNT_KEYS or STATS_REAL_KEYS - a required value that "
            f"nothing checks reaches verify() and is compared with >=.")
    absent = [k for k in STATS_REQUIRED_VALUES if k not in stats]
    if absent:
        raise TaskFailure(f"{where}: the statistics {absent} are absent; there is no default for "
                          f"any of them and verify() cannot be called without them")
    unmeasured = [k for k in STATS_REQUIRED_VALUES if _is_missing(stats[k])]
    if unmeasured:
        raise TaskFailure(
            f"{where}: the statistics {unmeasured} are None, NaN or a missing-value sentinel. "
            f"Every test in verify() is a comparison and a sentinel is False against all of them, "
            f"so passing these through would return a verdict of USABLE on a matrix nothing was "
            f"measured from.")
    for key in STATS_COUNT_KEYS:
        _finite_number(stats[key], key, where, whole=True, non_negative=True)
    for key in STATS_REAL_KEYS:
        # Not required to be positive: these are sums over a row, and a matrix that has been
        # scaled or centred has negative ones. Whether that disqualifies it as raw is verify()'s
        # verdict to render; what is refused here is a value it cannot render one on at all.
        _finite_number(stats[key], key, where, whole=False, non_negative=False)
    flag = _plain_scalar(stats["integer_counts"])
    if not isinstance(flag, bool) and not (isinstance(flag, int) and flag in (0, 1)):
        raise TaskFailure(
            f"{where}: integer_counts is {stats['integer_counts']!r}, which is not a yes or a no. "
            f"It answers whether the values are counts; a string or a number here would be read "
            f"through truthiness, where 'False' and 0.5 are both True.")
    if "name" in stats:
        name = stats["name"]
        if not isinstance(name, str) or not name.strip():
            raise TaskFailure(
                f"{where}: name is {name!r}. verify() prints the name on every reason it returns, "
                f"and an unattributable verdict in a cohort report is not usable evidence.")
    if "expected_genes" in stats:
        _require_positive_int_or_none(stats["expected_genes"], f"{where}: expected_genes")


def detect_format(path) -> str:
    """Which reader this path needs: 'mtx_dir', 'h5ad', 'h5' or 'tar'.

    Decided from the path, so it is testable without a filesystem, except for the directory case
    which can only be decided by asking. An unrecognised extension is refused rather than guessed:
    a matrix read by the wrong reader either fails loudly or - the case worth refusing for -
    succeeds with the axes the other way round.
    """
    p = Path(path)
    if p.is_dir():
        return "mtx_dir"
    name = p.name.lower()
    if name.endswith(".h5ad"):
        return "h5ad"
    if name.endswith(".h5") or name.endswith(".hdf5"):
        return "h5"
    for suf in _TAR_SUFFIXES:
        if name.endswith(suf):
            return "tar"
    raise TaskFailure(
        f"cannot tell what kind of matrix {p} is from its name.\n"
        f"  recognised: a directory holding matrix.mtx[.gz] with barcodes.tsv[.gz] and "
        f"features/genes.tsv[.gz]; a .h5ad; a .h5 or .hdf5 (10x or CellBender); "
        f"a tar archive ({', '.join(_TAR_SUFFIXES)}) containing an MTX directory.")


def _first_existing(directory: Path, names: Sequence[str]) -> Optional[Path]:
    for n in names:
        c = directory / n
        if c.is_file():
            return c
    return None


def _open_text(path: Path):
    """Open a .tsv or .tsv.gz as text. Nothing else in this module knows about compression."""
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def called_barcodes(path) -> list:
    """The barcodes a cell caller kept, from a directory or a barcode file.

    Step 2 compares the aligner's cell call with the denoiser's, and the number that decides the
    gate is how many the aligner called that the denoiser did NOT - which needs the SETS, not two
    counts. Two callers can agree on 21,000 cells and disagree about which 21,000.

    Accepts what the two callers actually emit:
      * a directory holding `barcodes.tsv[.gz]` - CeleScope `outs/filtered`, CellRanger
        `filtered_feature_bc_matrix`;
      * that file directly;
      * a one-column CSV, which is CellBender's `_cell_barcodes.csv` shape.

    Deliberately stdlib-only: this runs in the orchestrator, a bare interpreter on a cluster, and
    reading a list of strings is not a reason to require the analysis stack.
    """
    p = Path(path)
    if p.is_dir():
        for n in _BARCODE_NAMES:
            if (p / n).is_file():
                p = p / n
                break
        else:
            raise TaskFailure(
                f"{path} holds no {' or '.join(_BARCODE_NAMES)}, so the cells it called cannot be "
                f"read. A directory that is not a called-cell matrix is not an empty one.")
    if not p.is_file():
        raise TaskFailure(f"{p} does not exist, so its cell call cannot be read.")

    out, seen = [], set()
    with _open_text(p) as fh:
        for line in fh:
            # The first field, whether the separator is a tab (10x/CeleScope) or a comma
            # (CellBender). Suffixes such as `-1` are NOT stripped: they distinguish barcodes
            # within a run, and trimming them here would merge distinct cells into one.
            bc = line.strip().split("\t")[0].split(",")[0].strip().strip('"')
            if not bc:
                continue
            if bc.lower() in ("barcode", "barcodes", "cell_barcode", "cell_barcodes"):
                continue
            if bc in seen:
                raise TaskFailure(
                    f"{p} lists barcode {bc!r} more than once. Counting cells from it would "
                    f"over-count, and the overlap with the other caller would be wrong in the "
                    f"direction that flatters the comparison.")
            seen.add(bc)
            out.append(bc)
    if not out:
        raise TaskFailure(
            f"{p} lists no barcodes. That is not a cell call of zero to carry forward; a caller "
            f"that produced nothing has failed, and the step comparing it must say so.")
    return out


def _read_table_rows(path: Path) -> list:
    """Every non-empty line of a TSV, split on tabs. Blank trailing lines are not rows.

    Nothing is deduplicated and nothing is reordered: the row order IS the matrix's column or row
    order, and a reader that tidies it silently mislabels every gene after the first duplicate.
    """
    rows = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    if not rows:
        raise TaskFailure(f"{path} contains no rows - the matrix cannot be labelled from it")
    return rows


def _decode(values) -> list:
    """HDF5 string datasets come back as bytes on some builds and str on others."""
    out = []
    for v in values:
        out.append(v.decode("utf-8", "replace") if isinstance(v, bytes) else str(v))
    return out


# ---------------------------------------------------------------------------------------------
# reading


def _read_mtx_dir(directory: Path):
    """Read a 10x/CeleScope MTX triple into an AnnData of barcodes x features.

    Orientation is DETERMINED, not assumed: the stored shape is compared against the lengths of
    the barcode and feature tables. Market Matrix files from these producers are features x
    barcodes, but a matrix silently read the wrong way round produces a plausible object with
    every gene named as a cell, and that survives a long way downstream.
    """
    import numpy as np
    import pandas as pd
    import anndata as ad
    from scipy import sparse as sp
    from scipy import io as sio

    directory = Path(directory)
    mtx = _first_existing(directory, _MTX_NAMES)
    bcs = _first_existing(directory, _BARCODE_NAMES)
    fts = _first_existing(directory, _FEATURE_NAMES)
    missing = []
    if mtx is None:
        missing.append("matrix.mtx[.gz]")
    if bcs is None:
        missing.append("barcodes.tsv[.gz]")
    if fts is None:
        missing.append("features.tsv[.gz] or genes.tsv[.gz]")
    if missing:
        present = sorted(p.name for p in directory.iterdir()) if directory.is_dir() else []
        raise TaskFailure(
            f"{directory} is not an MTX directory: missing {', '.join(missing)}.\n"
            f"  files present: {', '.join(present) if present else '(none)'}")

    m = sio.mmread(str(mtx))
    X = sp.csr_matrix(m) if sp.issparse(m) else sp.csr_matrix(np.asarray(m))

    barcodes = [r[0] for r in _read_table_rows(bcs)]
    feat_rows = _read_table_rows(fts)
    ids = [r[0] for r in feat_rows]
    names = [r[1] if len(r) > 1 and r[1] else r[0] for r in feat_rows]
    types = [r[2] for r in feat_rows] if all(len(r) > 2 for r in feat_rows) else None

    n_b, n_f = len(barcodes), len(ids)
    if X.shape == (n_f, n_b):
        X = X.T.tocsr()
        orientation = "features x barcodes, transposed on read"
    elif X.shape == (n_b, n_f):
        orientation = "barcodes x features as stored"
    else:
        raise TaskFailure(
            f"{mtx} has shape {X.shape}, which matches neither "
            f"{n_f} features x {n_b} barcodes nor its transpose.\n"
            f"  barcodes read from {bcs.name}: {n_b}\n"
            f"  features read from {fts.name}: {n_f}\n"
            f"  The three files do not describe one matrix; nothing here can decide which is "
            f"wrong.")

    var = pd.DataFrame({"gene_ids": ids}, index=pd.Index(names, name=None))
    if types is not None:
        var["feature_types"] = types
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(barcodes, name=None)), var=var)
    adata.uns["scqc_source"] = str(directory)
    adata.uns["scqc_format"] = "mtx_dir"
    adata.uns["scqc_orientation"] = orientation
    adata.uns["scqc_duplicate_var_names"] = int(len(names) - len(set(names)))
    return adata


def _read_hdf5(path: Path):
    """Read a 10x or CellBender HDF5 into an AnnData of barcodes x features.

    The layout is read from the file rather than assumed from the extension. Three producers are
    recognised by structure and not by name: CellRanger v3 and CellBender 0.3 write the CSC block
    under `matrix/`, CellRanger v2 writes it under a group named for the genome, and CellBender
    0.2 writes it under `background_removed/`. All three are the same shape, so what is looked
    for is a single group carrying data/indices/indptr/shape. An unrecognised file is refused
    with the group names it actually contains, because a guess here reads the axes backwards.
    """
    import numpy as np
    import pandas as pd
    import anndata as ad
    import h5py
    from scipy import sparse as sp

    path = Path(path)
    with h5py.File(str(path), "r") as fh:
        groups = [k for k in fh.keys()
                  if isinstance(fh[k], h5py.Group) and set(_CSC_PARTS) <= set(fh[k].keys())]
        if len(groups) != 1:
            raise TaskFailure(
                f"{path} does not hold exactly one count matrix.\n"
                f"  groups carrying {', '.join(_CSC_PARTS)}: "
                f"{', '.join(groups) if groups else '(none)'}\n"
                f"  top-level keys: {', '.join(sorted(fh.keys())) or '(none)'}\n"
                f"  Recognised layouts: CellRanger v3 and CellBender 0.3 ('matrix'), "
                f"CellRanger v2 (one group per genome), CellBender 0.2 "
                f"('background_removed').")
        g = fh[groups[0]]
        data = np.asarray(g["data"][:])
        indices = np.asarray(g["indices"][:])
        indptr = np.asarray(g["indptr"][:])
        shape = tuple(int(x) for x in np.asarray(g["shape"][:]).ravel())
        if len(shape) != 2:
            raise TaskFailure(f"{path}: {groups[0]}/shape is {shape}, not a 2-D shape")
        if "barcodes" not in g:
            raise TaskFailure(f"{path}: {groups[0]} has no 'barcodes' dataset; "
                              f"datasets present: {', '.join(sorted(g.keys()))}")
        barcodes = _decode(g["barcodes"][:])

        ids = names = types = genome = None
        feats = g["features"] if "features" in g else None
        if isinstance(feats, h5py.Group):
            if "id" in feats:
                ids = _decode(feats["id"][:])
            if "name" in feats:
                names = _decode(feats["name"][:])
            if "feature_type" in feats:
                types = _decode(feats["feature_type"][:])
            if "genome" in feats:
                genome = _decode(feats["genome"][:])
        else:
            if "genes" in g:
                ids = _decode(g["genes"][:])
            if "gene_names" in g:
                names = _decode(g["gene_names"][:])
        if ids is None and names is None:
            raise TaskFailure(
                f"{path}: no feature table under {groups[0]} - looked for features/id, "
                f"features/name, genes and gene_names; datasets present: "
                f"{', '.join(sorted(g.keys()))}")

    index = names if names is not None else ids
    n_f, n_b = shape
    if len(barcodes) != n_b or len(index) != n_f:
        raise TaskFailure(
            f"{path}: stored shape {shape} says {n_f} features x {n_b} barcodes, but the file "
            f"carries {len(index)} feature labels and {len(barcodes)} barcodes. The matrix and "
            f"its labels do not describe the same object.")

    X = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()

    var_cols = {}
    if ids is not None:
        var_cols["gene_ids"] = ids
    if types is not None:
        var_cols["feature_types"] = types
    if genome is not None:
        var_cols["genome"] = genome
    var = pd.DataFrame(var_cols, index=pd.Index(index, name=None)) if var_cols \
        else pd.DataFrame(index=pd.Index(index, name=None))
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(barcodes, name=None)), var=var)
    adata.uns["scqc_source"] = str(path)
    adata.uns["scqc_format"] = "h5"
    adata.uns["scqc_h5_group"] = groups[0]
    adata.uns["scqc_orientation"] = "features x barcodes CSC, transposed on read"
    adata.uns["scqc_duplicate_var_names"] = int(len(index) - len(set(index)))
    return adata


def _read_h5ad(path: Path, backed: Optional[str] = None):
    import anndata as ad

    adata = ad.read_h5ad(str(path), backed=backed) if backed else ad.read_h5ad(str(path))
    if adata.X is None:
        raise TaskFailure(
            f"{path} has no X. An object carrying only layers or .raw cannot be verified as raw "
            f"input; name the layer to promote and do it explicitly, upstream of this adapter.")
    return adata


def _tar_destination(archive: Path, tmp_dir) -> Path:
    """Where an archive is unpacked. Never beside the archive."""
    archive = Path(archive)
    if tmp_dir is None:
        raise TaskFailure(
            f"{archive} is an archive and no tmp_dir was given. It is not unpacked beside itself: "
            f"vendor deliveries are read-only and stay unmodified, and an extraction written into "
            f"them destroys the boundary between what was delivered and what was computed. Pass "
            f"tmp_dir=<a scratch directory>.")
    dest_root = Path(tmp_dir).resolve()
    parent = archive.resolve().parent
    if dest_root == parent or parent in dest_root.parents:
        raise TaskFailure(
            f"tmp_dir {dest_root} is inside the archive's own directory {parent}. Unpacking there "
            f"writes into the delivery this pipeline must leave untouched; give a scratch "
            f"directory outside it.")
    return dest_root / (archive.name + ".extracted")


def _safe_extract(tar: tarfile.TarFile, dest: Path) -> None:
    """Unpack, refusing any member that would write outside `dest`.

    An archive is untrusted input like any other. A member named `../../etc/x`, or a symlink
    pointing out of the tree, writes wherever it likes during extraction, and the extraction of a
    vendor delivery is exactly where this pipeline promises not to write.
    """
    dest = Path(dest).resolve()

    def inside(target: Path) -> bool:
        target = Path(target)
        return target == dest or dest in target.parents

    for m in tar.getmembers():
        target = (dest / m.name).resolve()
        if not inside(target):
            raise TaskFailure(f"archive member {m.name!r} would be written to {target}, outside "
                              f"the extraction directory {dest}. Refusing to unpack it.")
        if m.issym() or m.islnk():
            link = (target.parent / m.linkname).resolve()
            if not inside(link):
                raise TaskFailure(f"archive member {m.name!r} links to {link}, outside {dest}")
        elif not (m.isfile() or m.isdir()):
            raise TaskFailure(f"archive member {m.name!r} is neither a file, a directory nor a "
                              f"link; refusing to unpack an archive containing it")
    try:
        tar.extractall(str(dest), filter="data")            # Python 3.12+
    except TypeError:                                       # pragma: no cover - older Pythons
        tar.extractall(str(dest))


def _extract_tar_mtx(archive: Path, tmp_dir, reuse_extracted: bool = False) -> Path:
    """Unpack an archive to a caller-supplied directory and return the one MTX directory in it.

    Two archives in this pipeline's world contain two matrices - a raw one and a filtered one -
    and choosing between them is the decision step 0 exists to force into the open. So an archive
    containing more than one MTX directory is refused with both paths named, rather than resolved
    by an ordering rule nobody would ever read.
    """
    archive = Path(archive)
    if not archive.is_file():
        raise TaskFailure(f"archive does not exist: {archive}")
    if not tarfile.is_tarfile(str(archive)):
        raise TaskFailure(f"{archive} is named like a tar archive but is not one")

    dest = _tar_destination(archive, tmp_dir)
    if dest.exists() and any(dest.iterdir()):
        if not reuse_extracted:
            raise TaskFailure(
                f"{dest} already exists and is not empty. Reading it would use a tree this run "
                f"did not produce, and an extraction older than its archive reads exactly like a "
                f"current one. Remove it, give a different tmp_dir, or pass reuse_extracted=True "
                f"to take responsibility for its contents.")
    else:
        dest.mkdir(parents=True, exist_ok=True)
        with tarfile.open(str(archive), "r:*") as tar:
            _safe_extract(tar, dest)

    found = set()
    for pattern in _MTX_NAMES:
        for hit in dest.rglob(pattern):
            found.add(hit.parent)
    if not found:
        top = sorted(p.name for p in dest.iterdir())[:20]
        raise TaskFailure(
            f"{archive} contains no matrix.mtx or matrix.mtx.gz.\n"
            f"  unpacked into {dest}; top level: {', '.join(top) if top else '(empty)'}")
    if len(found) > 1:
        listing = "\n".join(f"    {p}" for p in sorted(found))
        raise TaskFailure(
            f"{archive} contains {len(found)} MTX directories:\n{listing}\n"
            f"  Which one is the input is a decision, not a lookup - a raw and a filtered matrix "
            f"are the two candidates and they are not interchangeable. Pass the directory you "
            f"mean instead of the archive.")
    return sorted(found)[0]


def read_matrix(path, tmp_dir=None, *, backed: Optional[str] = None,
                reuse_extracted: bool = False):
    """Read a count matrix into an AnnData of barcodes x features.

    Accepts a 10x/CeleScope MTX directory, a .h5ad, a 10x or CellBender .h5, or a tar archive
    containing an MTX directory. The format is decided by `detect_format()`.

    `tmp_dir` is required for an archive and is where it is unpacked; nothing is ever written
    beside the input, which may be a read-only vendor delivery this pipeline must leave
    unmodified. `backed` applies only to .h5ad, where anndata supports it - requesting it for a
    format that cannot honour it is refused rather than quietly ignored, because a caller asking
    for backed mode is asking because the object does not fit in memory.

    Whatever goes wrong underneath surfaces as a TaskFailure naming the file and the format it
    was read as. A reader raising an OSError about an HDF5 object header, from three frames inside
    a library, does not tell the operator which of ten inputs it was.
    """
    p = Path(path)
    kind = detect_format(p)
    if backed is not None and kind != "h5ad":
        raise TaskFailure(f"backed={backed!r} was requested for a {kind} input. Only .h5ad can be "
                          f"opened backed; reading {p} loads it into memory and pretending "
                          f"otherwise would hide that.")
    if kind != "mtx_dir" and not p.is_file():
        raise TaskFailure(f"matrix does not exist: {p}")
    if kind == "mtx_dir" and not p.is_dir():
        raise TaskFailure(f"matrix directory does not exist: {p}")
    try:
        if kind == "mtx_dir":
            return _read_mtx_dir(p)
        if kind == "h5ad":
            return _read_h5ad(p, backed=backed)
        if kind == "h5":
            return _read_hdf5(p)
        return _read_mtx_dir(_extract_tar_mtx(p, tmp_dir, reuse_extracted=reuse_extracted))
    except TaskFailure:
        raise
    except ImportError as e:
        raise TaskFailure(
            f"reading {p} as {kind} needs a library this interpreter does not have: {e}. The "
            f"measurement runs where anndata is installed - see build_stats_argv() and "
            f"run_summary_stats(), which hand the work to another interpreter.") from e
    except Exception as e:
        raise TaskFailure(f"failed to read {p} as {kind}: "
                          f"{type(e).__name__}: {e}") from e


def _as_adata(path_or_adata, tmp_dir=None, backed: Optional[str] = None):
    """An AnnData is used as given; anything else is a path and is read."""
    if hasattr(path_or_adata, "n_obs") and hasattr(path_or_adata, "X"):
        return path_or_adata, None
    p = Path(path_or_adata)
    return read_matrix(p, tmp_dir=tmp_dir, backed=backed), p


# ---------------------------------------------------------------------------------------------
# measuring


def _iter_values(X, block: int = _VALUE_BLOCK) -> Iterator:
    """Yield the STORED values of X in blocks, never materialising it densely.

    For a sparse matrix the stored values are the non-zeros; the implicit zeros are not examined
    because zero is integral and non-negative, so it cannot change any answer this module asks of
    the values. Everything else - a dense array, an h5py dataset, an anndata backed sparse
    dataset - is sliced by rows, so the block bounds what is held at once.
    """
    import numpy as np
    from scipy import sparse as sp

    if sp.issparse(X):
        fmt = getattr(X, "format", "")
        data = X.data if fmt in ("csr", "csc", "coo", "bsr") else X.tocoo().data
        data = np.asarray(data).ravel()
        for i in range(0, data.size, block):
            yield data[i:i + block]
        return

    shape = getattr(X, "shape", None)
    if shape is None or len(shape) != 2:
        raise TaskFailure(f"cannot read values from a matrix of type {type(X).__name__} with "
                          f"shape {shape!r}")
    n_rows, n_cols = int(shape[0]), int(shape[1])
    rows_per = max(1, block // max(1, n_cols))
    for i in range(0, n_rows, rows_per):
        chunk = X[i:i + rows_per]
        if sp.issparse(chunk):
            yield np.asarray(chunk.data).ravel()
        else:
            yield np.asarray(chunk).reshape(-1)


def integer_check(X, block: int = _VALUE_BLOCK) -> bool:
    """Are these values integer counts?

    False - not an exception - for a normalised, scaled or log-transformed matrix, because that
    is a property of the input `verify()` is entitled to render a verdict on. False also for
    negative and for non-finite values, and the reason is that the question asked is whether these
    are COUNTS, not whether they are integers: -3 is an integer and is not a count, and a NaN
    reported as integral would travel into `verify()` as a pass. An unsigned or boolean dtype is
    integral by construction and is not scanned value by value.

    Nothing is densified. Sparse stored values are scanned in blocks; any other array is sliced by
    rows in blocks. A non-numeric dtype is a TaskFailure rather than a False, because it is not an
    answer about counts - it means the object is not a count matrix at all.
    """
    if X is None:
        raise TaskFailure("integer_check received no matrix (X is None)")
    if block < 1:
        raise TaskFailure(f"block must be at least 1, got {block}")
    # Imported after the argument checks: a bad argument is a bad argument whether or not numpy
    # is installed, and reporting a missing dependency instead sends the reader to the wrong fix.
    import numpy as np

    for chunk in _iter_values(X, block=block):
        a = np.asarray(chunk)
        if a.size == 0:
            continue
        kind = a.dtype.kind
        if kind == "b" or kind == "u":
            continue
        if kind == "i":
            if bool((a < 0).any()):
                return False
            continue
        if kind != "f":
            raise TaskFailure(
                f"matrix values have dtype {a.dtype} (kind {kind!r}), which is not numeric. "
                f"This is not a count matrix; nothing here can say whether it is raw.")
        if not bool(np.isfinite(a).all()):
            return False
        if bool((a < 0).any()):
            return False
        if not bool(np.all(a == np.rint(a))):
            return False
    return True


def _row_totals(X, block: int = _ROW_BLOCK):
    """Per-barcode total counts as float64, from sparse row sums. Never densifies.

    Summed in float64 rather than the stored dtype: an int32 matrix whose deepest barcode carries
    more than 2^31 counts overflows a same-dtype sum silently and negatively, and a negative total
    would then be read by `verify()` as a matrix with raw droplets. float64 is exact for integers
    up to 2^53, which is far past any real library.

    A non-finite total raises. It cannot be returned: every test in `verify()` is a comparison,
    and NaN compares False against all of them, so a NaN minimum reads as "raw droplets, pass".
    """
    import numpy as np
    from scipy import sparse as sp

    if X is None:
        raise TaskFailure("matrix has no X - there is nothing to total")
    shape = getattr(X, "shape", None)
    if shape is None or len(shape) != 2:
        raise TaskFailure(f"cannot total a matrix of type {type(X).__name__} with shape "
                          f"{shape!r}")
    n_rows = int(shape[0])
    if n_rows == 0:
        raise TaskFailure("matrix has zero barcodes; there is no distribution to summarise")

    if sp.issparse(X) or isinstance(X, np.ndarray):
        totals = np.asarray(X.sum(axis=1, dtype=np.float64)).ravel()
    else:
        totals = np.empty(n_rows, dtype=np.float64)
        for i in range(0, n_rows, block):
            chunk = X[i:i + block]
            if sp.issparse(chunk):
                totals[i:i + block] = np.asarray(chunk.sum(axis=1, dtype=np.float64)).ravel()
            else:
                totals[i:i + block] = np.asarray(chunk).sum(axis=1, dtype=np.float64)

    if totals.shape[0] != n_rows:
        raise TaskFailure(f"row totals came back with {totals.shape[0]} entries for {n_rows} "
                          f"barcodes")
    bad = ~np.isfinite(totals)
    if bool(bad.any()):
        where = np.flatnonzero(bad)[:5].tolist()
        raise TaskFailure(
            f"{int(bad.sum())} of {n_rows} barcodes have a non-finite total (first at row "
            f"index {where}). Returning them would put NaN into every threshold comparison "
            f"downstream, where NaN is False against all of them and therefore reads as a pass.")
    return totals


def summary_stats(path_or_adata, expected_genes=None, *, name=None, tmp_dir=None,
                  backed: Optional[str] = None) -> dict:
    """The statistics `lib/verify_raw.py::verify` renders a verdict on, plus n_nonzero_barcodes.

    Returns exactly `STATS_KEYS`, always all of them. `p98_counts` is the 98th percentile of
    per-barcode total counts over EVERY barcode in the file, taken with numpy's default linear
    interpolation; see this module's docstring for what that means on a raw droplet matrix.

    `name` labels the verdict. For a path it defaults to the path as given, which is what
    `modules/00_ingest/ingest.plan_one` passes; for an AnnData already in memory there is nothing
    to default to and it must be supplied, because `verify()` prints the name on every reason it
    returns and an unattributable verdict in a cohort report is not usable evidence.

    The returned dict cannot be splatted into `verify()` - it carries `n_nonzero_barcodes`, which
    `verify()` does not accept. Use `verify_kwargs()`.
    """
    expected = _require_positive_int_or_none(expected_genes, "expected_genes")
    import numpy as np

    adata, source = _as_adata(path_or_adata, tmp_dir=tmp_dir, backed=backed)
    # `_is_missing`, not `is None`: a name read from a samplesheet cell that was blank arrives as
    # float nan or pandas.NA, and `str()` of either is a label - "nan", "<NA>" - that looks like a
    # sample and identifies nothing. An absent name is absent however it was spelled.
    if _is_missing(name):
        if source is None:
            raise TaskFailure(
                "summary_stats was given an AnnData and no name. An in-memory object carries no "
                "path, and verify() prints the name on every reason it returns - pass "
                "name='<what this matrix is>'.")
        label = str(path_or_adata)
    else:
        label = str(name)

    n_barcodes = int(adata.n_obs)
    n_genes = int(adata.n_vars)
    if n_barcodes == 0:
        raise TaskFailure(f"{label}: the matrix has zero barcodes")
    if n_genes == 0:
        raise TaskFailure(f"{label}: the matrix has zero features")

    totals = _row_totals(adata.X)
    stats = {
        "name": label,
        "n_barcodes": n_barcodes,
        "n_genes": n_genes,
        "min_counts": float(totals.min()),
        "max_counts": float(totals.max()),
        "p98_counts": float(np.percentile(totals, 98.0)),
        "expected_genes": expected,
        "integer_counts": bool(integer_check(adata.X)),
        "n_nonzero_barcodes": int((totals > 0).sum()),
    }
    missing = [k for k in STATS_KEYS if k not in stats]
    if missing:
        raise TaskFailure(f"summary_stats built an incomplete result, missing {missing}")
    # Checked on the way out as well as on the way in. This is the constructor of the payload
    # every gate downstream reads, so a statistic that is not a finite number must fail here -
    # beside the matrix that produced it and the code that computed it - rather than in a parser
    # two processes away that can only say the file is unusable.
    require_measured_stats(stats, f"summary_stats({label!r})")
    return stats


def verify_kwargs(stats: dict, include_name: bool = True) -> dict:
    """The subset of `stats` that `verify()` accepts, validated before it gets there.

    `summary_stats()` returns one key `verify()` does not take, so splatting its result raises
    TypeError. This selects, and it also runs `require_measured_stats()` - the point of that check
    being that `verify()` cannot make it itself: a NaN reaching it produces a clean PASS rather
    than an error, a string raises TypeError from inside a comparison, and an infinity satisfies
    every floor it is tested against.

    `include_name=False` drops `name`, which is what `modules/00_ingest/ingest.plan_one` needs: it
    passes the name itself, so a `stats_fn` returning one raises TypeError for a duplicate
    keyword.
    """
    if not isinstance(stats, dict):
        raise TaskFailure(f"verify_kwargs expects a dict of statistics, got "
                          f"{type(stats).__name__}")
    missing = [k for k in VERIFY_KEYS if k not in stats]
    if missing:
        raise TaskFailure(f"statistics are missing {missing}; verify() cannot be called without "
                          f"them and there is no default for any of them")
    require_measured_stats(stats, "the statistics handed to verify_kwargs")
    out = {k: stats[k] for k in VERIFY_KEYS}
    if not include_name:
        out.pop("name")
    return out


def ingest_stats_fn(expected_genes=None, *, tmp_dir=None, backed: Optional[str] = None):
    """A `stats_fn` for `modules/00_ingest/ingest.plan_one`, which calls `stats_fn(matrix_path)`.

    plan_one does `verify(name=matrix, **stats_fn(matrix))`, so the callable must return verify's
    keyword arguments WITHOUT `name` and without anything verify does not accept. Constructing it
    here keeps that contract in one place rather than in every orchestrator that wires the two
    together.
    """
    def stats_fn(matrix_path) -> dict:
        return verify_kwargs(
            summary_stats(matrix_path, expected_genes=expected_genes, tmp_dir=tmp_dir,
                          backed=backed),
            include_name=False)
    return stats_fn


def barcode_rank(path_or_adata, n: int = 5000, *, tmp_dir=None,
                 backed: Optional[str] = None) -> list:
    """The barcode-rank curve for figure F1: [(rank, total_counts), ...], rank 1 the deepest.

    Downsampled log-uniformly to at most `n` points, because the curve is read on log-log axes
    where hundreds of thousands of points at the tail overplot into a solid band and the knee -
    the only part anyone reads - is decided by the first few hundred. The spacing is deterministic
    and there is no sampling, so no seed is threaded through: log-spaced ranks collide at the head
    and are deduplicated, which is why the result can be shorter than `n`.

    Barcodes with a zero total are KEPT in the curve. They are not plottable on a log axis and the
    figure may drop them, but they are the evidence that a matrix still contains its empty
    droplets, and dropping them here would make a cell-called matrix and a raw one produce the
    same picture. Removing observations is step 7's job, not a plotting helper's.
    """
    if n < 2:
        raise TaskFailure(f"barcode_rank needs at least 2 points, got n={n}")
    import numpy as np

    adata, _ = _as_adata(path_or_adata, tmp_dir=tmp_dir, backed=backed)
    totals = _row_totals(adata.X)
    order = np.argsort(-totals, kind="stable")
    ordered = totals[order]
    size = int(ordered.size)
    if size == 0:
        raise TaskFailure("the matrix has no barcodes; there is no rank curve to draw")

    if size <= n:
        ranks = np.arange(1, size + 1, dtype=np.int64)
    else:
        spaced = np.rint(np.logspace(0.0, np.log10(float(size)), n)).astype(np.int64)
        ranks = np.unique(np.concatenate(([1], spaced, [size])))
        ranks = ranks[(ranks >= 1) & (ranks <= size)]
    return [(int(r), float(ordered[r - 1])) for r in ranks]


# ---------------------------------------------------------------------------------------------
# writing


def _nullable_string_dtype(dtype) -> bool:
    """True for the pandas string dtypes anndata's writer declines: StringDtype and arrow strings.

    Detected by asking pandas, then by name, because the arrow-backed spelling is an `ArrowDtype`
    rather than a `StringDtype` and the two are not related by inheritance.
    """
    import pandas as pd

    string_dtype = getattr(pd, "StringDtype", None)
    if string_dtype is not None and isinstance(dtype, string_dtype):
        return True
    return str(dtype).lower() in ("string", "str", "string[python]", "string[pyarrow]",
                                  "string[pyarrow_numpy]", "large_string[pyarrow]")


def _object_backed(labels, what: str, where: str):
    """One pandas Index or Series, re-backed by a plain object array. Values are not touched.

    A missing value stops it. `astype(object)` turns `pandas.NA` into a sentinel sitting in an
    object array, h5py then refuses the array with `Can't implicitly convert non-string objects to
    strings`, and the alternatives - writing "nan" as a label, or dropping it - are both a change
    to the data that this adapter is not entitled to make on the caller's behalf.
    """
    if not _nullable_string_dtype(labels.dtype):
        return labels
    hasnans = getattr(labels, "hasnans", None)
    if hasnans is None:
        hasnans = bool(labels.isna().any())
    if bool(hasnans):
        raise TaskFailure(
            f"{where}: {what} is a nullable string column that contains missing values, and "
            f"anndata will not write it.\n"
            f"  Casting it to plain strings here would have to invent a spelling for the missing "
            f"entries - 'nan', or an empty label - and a label this adapter invented is "
            f"indistinguishable downstream from one the data carried. Decide what the missing "
            f"entries mean and fill or drop them upstream, where the decision is recorded.")
    return labels.astype(object)


def object_backed_labels(adata, where: str = "write_h5ad"):
    """`(obs, var)` re-backed by object arrays where anndata declines to write them.

    pandas 3.0 gives every string column and every string index a `StringDtype` by default, and
    anndata refuses to write those unless `anndata.settings.allow_write_nullable_strings` is
    turned on - so on that stack `write_h5ad(read_matrix(...))` raises, which is the round trip
    this adapter exists to perform. Both readers here produce exactly such frames: the barcode
    index, the feature index, `gene_ids`, `feature_types` and `genome` are all built from lists of
    Python strings.

    The fix is a dtype cast rather than the setting, and the choice is deliberate:

      * the global setting changes how EVERY object written in the process is stored, which is not
        a decision an adapter makes on the caller's behalf - the same reasoning is recorded in
        `adapters/apply_filter.py`;
      * a file written with nullable string arrays cannot be read by anndata < 0.11, so turning
        the setting on trades a write failure here for a read failure on someone else's machine;
      * the cast is lossless. The values are the same Python `str` objects either way; only the
        array backing them differs, and object-backed strings are what pandas 2 and every anndata
        version can write and read.

    Returns `(None, None)` when nothing needs changing, so the common case touches nothing.

    `where` names the caller in any refusal. It is not cosmetic: `adapters/apply_filter.py` also
    applies this cast, and a message attributing its refusal to `write_h5ad` sends the reader to
    the wrong function on the one path where the removal has already been approved and recorded.
    """
    # pandas is imported by `_nullable_string_dtype()` on the first dtype it is asked about, which
    # happens below for every frame - there is no path through here that reaches anndata's writer
    # without pandas having been imported first.
    obs, var = adata.obs, adata.var
    new = {}
    for what, frame in (("obs", obs), ("var", var)):
        cols = [c for c in frame.columns if _nullable_string_dtype(frame[c].dtype)]
        index_needs = _nullable_string_dtype(frame.index.dtype)
        if not cols and not index_needs:
            continue
        out = frame.copy(deep=False)
        for c in cols:
            out[c] = _object_backed(frame[c], f".{what}[{c!r}]", where)
        if index_needs:
            out.index = _object_backed(frame.index, f".{what}_names", where)
        new[what] = out
    if not new:
        return None, None
    return new.get("obs"), new.get("var")


def new_write_token() -> str:
    """A value that exists nowhere - on disk, in any other process - until this call invents it.

    A uuid4 rather than a timestamp, a counter or a digest of the object being written. Each of
    those can be reproduced by something that is not this call: a clock by reading the clock, a
    counter by counting, and a content digest by the previous run's own output, which is exactly
    the artifact the token has to be able to tell apart from this one.
    """
    return "scqc-write-" + uuid.uuid4().hex


def _decode_scalar_string(value) -> Optional[str]:
    """One HDF5 or numpy scalar as `str`; None when it is not a string at all.

    h5py returns variable-length strings as `bytes` on some builds and `str` on others, and a
    scalar dataset read with `[()]` can arrive as a 0-d numpy array wrapping either.
    """
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, str):
        return value
    item = getattr(value, "item", None)
    if callable(item):
        try:
            inner = item()
        except (TypeError, ValueError):
            return None
        if isinstance(inner, (bytes, str)):
            return _decode_scalar_string(inner)
    return None


def read_write_token(path) -> Optional[str]:
    """The token `write_h5ad()` stored inside an .h5ad, or None when the file carries none.

    None means "this file has no token" and is a real answer about a real file - a .h5ad written
    by anything other than `write_h5ad()`, including an earlier checkout of this adapter. A file
    that cannot be opened at all is NOT that answer and raises instead, because a read failure
    reported as "no token" would be indistinguishable from a file that simply predates the check.

    h5py first, because it reads one scalar and no matrix; anndata second, because it necessarily
    understands any layout anndata just wrote and a verification that fires on a correct write is
    a verification somebody will delete. Both failing is reported with both reasons.
    """
    p = Path(path)
    try:
        import h5py

        with h5py.File(str(p), "r") as fh:
            uns = fh.get("uns")
            if uns is None or WRITE_TOKEN_KEY not in uns:
                return None
            return _decode_scalar_string(uns[WRITE_TOKEN_KEY][()])
    except Exception as first:  # noqa: BLE001 - re-raised below with the second attempt's reason
        try:
            import anndata

            obj = anndata.read_h5ad(str(p), backed="r")
            try:
                uns = getattr(obj, "uns", None)
                if uns is None or WRITE_TOKEN_KEY not in uns:
                    return None
                return _decode_scalar_string(uns[WRITE_TOKEN_KEY])
            finally:
                fh = getattr(obj, "file", None)
                if fh is not None:
                    try:
                        fh.close()
                    except Exception:  # noqa: BLE001 - closing a read handle is best effort
                        pass
        except Exception as second:  # noqa: BLE001
            raise TaskFailure(
                f"{p} could not be opened to read back the write token that says which call "
                f"produced it.\n"
                f"  h5py: {first!r}\n"
                f"  anndata: {second!r}\n"
                f"  The file is reported as unverifiable rather than assumed to be this run's; "
                f"nothing downstream may treat it as an output.") from second


def write_h5ad(adata, path, compression: str = "gzip"):
    """Write an AnnData and confirm THIS call produced the file. Returns the path that was checked.

    The check is the point. A write that raised nothing and produced nothing is discovered
    otherwise three steps later, as a confusing error about a different file, and `engine/task.py`
    exists to stop exactly that: an output is a file that was verified to exist, never a path that
    was passed to a writer.

    Existence alone is not that check, which is why the path is deleted first: an .h5ad left at the
    same path by an earlier run satisfies every test a writer's silence would otherwise escape - it
    opens, it reads, it counts - and the previous parameters' object is then reported as this run's
    output.

    WHAT THE FRESHNESS CHECK IS, AND WHAT IT CAN AND CANNOT DETECT

    After the write the file must exist, be non-empty, start with the HDF5 signature, and carry in
    `uns[WRITE_TOKEN_KEY]` the uuid4 this call invented after observing the path absent. That token
    is the evidence. No time is compared: the previous spelling compared mtimes with a two-second
    allowance for clock and filesystem resolution, and a restored artifact back-dated by one second
    was accepted while the same construction at two hours was refused - a tolerance window is
    exactly where a restored file is placed, and `os.utime` puts a file's mtime anywhere at all, so
    no width of window closes it.

      Detects: a writer that returned without writing and left the previous run's file standing; an
      artifact restored from a backup, a copy, a rename or a hardlink, whatever its mtime; a file
      written by an earlier checkout of this adapter, which carries no token; a file written by any
      other producer; a write that landed somewhere other than this path.

      Does not detect: anything about the CONTENT beyond its provenance - the token says these
      bytes came from this call's writer, not that X, obs and var hold what the caller intended,
      and not that the file is complete beyond the token being readable. It says nothing about
      what happens to the path after this function returns, and it cannot distinguish two writes
      of the same object from one. A process that could read this process's memory could copy the
      token; nothing short of that can produce it.

    Anything anndata declines to store as it stands is normalised losslessly first, by
    `object_backed_labels()`; on pandas 3.0 that is what makes `write_h5ad(read_matrix(x))` work at
    all. Nothing about the caller's object is left changed: the re-backed frames and the token are
    attached for the duration of the write and the originals are put back afterwards.
    """
    if not hasattr(adata, "write_h5ad"):
        raise TaskFailure(f"write_h5ad expects an AnnData, got {type(adata).__name__}")
    if compression not in ("gzip", "lzf", None):
        raise TaskFailure(f"compression={compression!r} is not one h5py provides; use 'gzip', "
                          f"'lzf' or None")
    p = Path(path)
    if p.is_dir():
        raise TaskFailure(f"{p} is a directory; refusing to write an .h5ad over it")
    p.parent.mkdir(parents=True, exist_ok=True)
    new_obs, new_var = object_backed_labels(adata)
    old_obs, old_var = adata.obs, adata.var
    uns = getattr(adata, "uns", None)
    if uns is None or not hasattr(uns, "__setitem__"):
        raise TaskFailure(
            f"write_h5ad was given a {type(adata).__name__} whose .uns cannot hold an entry, so "
            f"the write token that identifies this call's output could not be stored in the file. "
            f"An output whose provenance cannot be recorded is not written.")
    had_token = WRITE_TOKEN_KEY in uns
    old_token = uns[WRITE_TOKEN_KEY] if had_token else None
    # Invented AFTER the path is observed absent, so it cannot have reached any file already there.
    started = _clear_previous_output(p)
    token = new_write_token()
    try:
        if new_obs is not None:
            adata.obs = new_obs
        if new_var is not None:
            adata.var = new_var
        adata.uns[WRITE_TOKEN_KEY] = token
        adata.write_h5ad(str(p), compression=compression)
    except TaskFailure:
        raise
    except Exception as e:
        # anndata's own message carries the remedy for the cases that actually occur here - a
        # dtype its writer declines, a full filesystem, a path that cannot be created - so it is
        # quoted rather than reinterpreted. What is added is which file was being written.
        raise TaskFailure(f"failed to write {p}: {type(e).__name__}: {e}") from e
    finally:
        # Restored unconditionally. anndata's writer also converts string columns to categoricals
        # in place as it goes, so without this the caller's object comes back from a write with
        # dtypes it did not have - and a function that quietly rewrites its argument is a worse
        # defect than the one being fixed. The token is removed again for the same reason: it
        # describes one write, not the object.
        if new_obs is not None:
            adata.obs = old_obs
        if new_var is not None:
            adata.var = old_var
        if had_token:
            adata.uns[WRITE_TOKEN_KEY] = old_token
        else:
            try:
                del adata.uns[WRITE_TOKEN_KEY]
            except KeyError:
                pass

    if not p.is_file():
        raise TaskFailure(f"write_h5ad returned but {p} does not exist. The object was not "
                          f"written; nothing downstream may treat this as an output.")
    size = p.stat().st_size
    if size == 0:
        raise TaskFailure(f"{p} was created but is empty (0 bytes) - the write did not complete")
    with p.open("rb") as fh:
        magic = fh.read(8)
    if magic != b"\x89HDF\r\n\x1a\n":
        raise TaskFailure(f"{p} exists ({size} bytes) but does not start with the HDF5 signature; "
                          f"it is not a readable .h5ad")
    stored = read_write_token(p)
    # Reported for the reader, never compared: the age of a file is not evidence about who wrote
    # it, and this line exists so that a refusal can be recognised in a log, not to decide one.
    age = started - p.stat().st_mtime
    if stored is None:
        raise TaskFailure(
            f"{p} exists ({size:,} bytes, mtime {age:+.1f}s relative to the start of this write) "
            f"but carries no {WRITE_TOKEN_KEY} in its /uns, so it is not the file this call "
            f"produced.\n"
            f"  The path was deleted and observed absent immediately beforehand and the writer was "
            f"handed a token to store. An .h5ad here without one was written by something else - "
            f"a restored artifact, or an earlier checkout of this adapter - and its contents "
            f"describe different parameters.")
    if stored != token:
        raise TaskFailure(
            f"{p} carries the write token {stored!r}, not this call's {token!r} (size {size:,} "
            f"bytes, mtime {age:+.1f}s relative to the start of this write).\n"
            f"  It is a different write's output standing at this path: the file was deleted and "
            f"observed absent immediately before the writer ran, so what is here now was restored "
            f"or copied rather than written. Its numbers describe different parameters.")
    return p


# ---------------------------------------------------------------------------------------------
# argument construction and output parsing - pure, and testable without either tool present


def build_stats_argv(python_exe, matrix, out_json, *, expected_genes=None, name=None,
                     tmp_dir=None, rank_points: int = 0, module_path=None,
                     reuse_extracted: bool = False, write_token=None) -> list:
    """The argv that makes a Python interpreter measure one matrix and write the JSON payload.

    Separated from running it because the interpreter that measures the matrix is usually not the
    one running the orchestrator: the aligner, the denoiser and the analysis stack have
    incompatible pins and live in different environments, so scanpy and anndata are frequently
    absent from the process making the decisions. This builds a command; `run_summary_stats()`
    hands it to an executor.

    The module is invoked by FILE PATH rather than as `-m adapters.matrix`, so the measuring
    environment needs nothing installed and nothing on PYTHONPATH.

    `write_token` is echoed into the payload the subprocess writes, so the caller can establish
    that the JSON now at `out_json` is the one this invocation produced rather than an earlier
    one left there. It is passed on the command line and never through the environment, because an
    inherited environment is the one thing a restored artifact's producer also had.
    """
    mod = Path(module_path) if module_path is not None else Path(__file__).resolve()
    argv = [str(python_exe), str(mod), "stats", "--matrix", str(matrix), "--out", str(out_json)]
    if write_token is not None:
        tok = str(write_token)
        if not tok.strip() or any(c.isspace() for c in tok):
            raise TaskFailure(
                f"write_token={write_token!r} is blank or contains whitespace. It is compared "
                f"verbatim against the value read back out of the payload; a token that does not "
                f"survive an argv round trip would refuse every correct run.")
        argv += ["--write-token", tok]
    if expected_genes is not None:
        argv += ["--expected-genes", str(_require_positive_int_or_none(expected_genes,
                                                                      "expected_genes"))]
    if name is not None:
        argv += ["--name", str(name)]
    if tmp_dir is not None:
        argv += ["--tmp-dir", str(tmp_dir)]
    if rank_points:
        if rank_points < 2:
            raise TaskFailure(f"rank_points must be 0 or at least 2, got {rank_points}")
        argv += ["--rank-points", str(int(rank_points))]
    if reuse_extracted:
        argv += ["--reuse-extracted"]
    return argv


def _reject_json_constant(token: str):
    raise TaskFailure(
        f"the statistics JSON contains the literal {token}. A non-finite number here would be "
        f"compared against every threshold downstream and lose silently; the file is not a "
        f"measurement and is refused rather than parsed.")


def parse_stats_json(text: str) -> dict:
    """Validate and return the payload written by this module's `stats` entry point.

    Every required key must be present, measured, numeric and FINITE. A missing key is refused
    rather than defaulted, because the whole value of this payload is that the numbers in it were
    measured from a file somebody can open.

    `parse_constant` catches only the bare JSON literals `NaN`, `Infinity` and `-Infinity`, and
    that is not the same as catching non-finite values: `1e999` is an ordinary JSON number that
    `json.loads` turns into `float('inf')` without consulting `parse_constant` at all, and
    `"abc"` is an ordinary JSON string that no missing-value test rejects. Both used to reach
    `verify()` and meet `>=` there - the infinity satisfying every floor silently, the string
    raising TypeError from a frame that cannot name the file. `require_measured_stats()` is
    therefore run over the parsed payload as well, and it names the key it refuses.
    """
    try:
        payload = json.loads(text, parse_constant=_reject_json_constant)
    except json.JSONDecodeError as e:
        raise TaskFailure(f"the statistics JSON could not be parsed: {e}") from None
    if not isinstance(payload, dict):
        raise TaskFailure(f"the statistics JSON is a {type(payload).__name__}, not an object")
    if payload.get("schema") != STATS_SCHEMA:
        raise TaskFailure(f"the statistics JSON declares schema {payload.get('schema')!r}, "
                          f"expected {STATS_SCHEMA!r}. It was written by a different version of "
                          f"this adapter and its keys cannot be assumed.")
    for key in ("source", "format"):
        if not isinstance(payload.get(key), str) or not payload[key]:
            raise TaskFailure(f"the statistics JSON has no '{key}'; a measurement that cannot say "
                              f"which file it came from is not traceable to anything")
    stats = payload.get("stats")
    if not isinstance(stats, dict):
        raise TaskFailure("the statistics JSON has no 'stats' object")
    missing = [k for k in STATS_KEYS if k not in stats]
    if missing:
        raise TaskFailure(f"the statistics JSON is missing {missing}")
    require_measured_stats(stats, "the statistics JSON")
    versions = payload.get("versions")
    if not isinstance(versions, dict) or not versions:
        raise TaskFailure("the statistics JSON records no tool versions; a result whose "
                          "provenance was not observed cannot be told from one that was")
    if "write_token" in payload:
        # Checked for shape here and for VALUE in run_summary_stats(), which is the only place
        # that knows which token was issued. A payload carrying a blank one is refused rather than
        # read as "no token", because the two are different claims and only one of them is honest.
        tok = payload["write_token"]
        if not isinstance(tok, str) or not tok.strip():
            raise TaskFailure(
                f"the statistics JSON carries write_token={tok!r}, which is not a token. It is "
                f"what says which invocation produced this file; a blank one identifies nothing.")
    rank = payload.get("barcode_rank")
    if rank is not None:
        if not isinstance(rank, list):
            raise TaskFailure("'barcode_rank' is present but is not a list of [rank, total] pairs")
        for row in rank:
            if not (isinstance(row, (list, tuple)) and len(row) == 2):
                raise TaskFailure(f"'barcode_rank' contains a malformed row: {row!r}")
    return payload


def run_summary_stats(matrix, out_json, log, python_exe, *, expected_genes=None, name=None,
                      tmp_dir=None, rank_points: int = 0, module_path=None,
                      reuse_extracted: bool = False, env=None, cwd=None, timeout_s=None,
                      executor=None) -> dict:
    """Measure one matrix in another interpreter and return {'outputs', 'metrics', 'versions'}.

    `executor` is required; it is keyword-only because a default would have to be a local one, and
    silently running on the submitting host a step that was meant for a compute node is the kind
    of substitution this repository refuses everywhere else.

    `out_json` is deleted and observed absent BEFORE the command runs, and the payload it writes
    must carry the token this call issued on the command line. Existence alone is not evidence: a
    command that exits 0 and writes nothing leaves the previous run's numbers in place, they
    parse, they are returned as metrics, and they describe a different matrix measured under
    different parameters. Nothing downstream can detect that, which is why it is settled here -
    and settled with a token rather than an mtime, for the reason `write_h5ad()` sets out at
    length: a time comparison needs a tolerance window and a restored file lands inside it.

    That check refuses a payload written by a copy of this module old enough not to understand
    `--write-token`, and does so by name rather than by passing the run: an unverifiable
    measurement and a verified one must not look alike.

    The versions returned are the ones the MEASURING process reported for the libraries it
    actually imported - not the ones the orchestrator's own interpreter has. On a scheduler those
    are different machines and frequently different stacks.
    """
    if executor is None:
        raise TaskFailure("run_summary_stats needs an executor; pass executor=LocalExecutor() or "
                          "a PBSExecutor explicitly")
    matrix = Path(matrix)
    out = Path(out_json)
    # Invented after the path is observed absent, so no file already present can be carrying it.
    _clear_previous_output(out)
    token = new_write_token()

    argv = build_stats_argv(python_exe, matrix, out, expected_genes=expected_genes, name=name,
                            tmp_dir=tmp_dir, rank_points=rank_points, module_path=module_path,
                            reuse_extracted=reuse_extracted, write_token=token)
    captured = executor.shell(argv, log=Path(log), env=env, cwd=cwd, timeout_s=timeout_s)

    if not out.is_file() or out.stat().st_size == 0:
        tail = "\n".join((captured or "").strip().splitlines()[-15:])
        raise TaskFailure(
            f"{argv[0]} exited 0 but wrote no statistics to {out}.\n  log: {log}\n"
            f"  last lines:\n{tail}")
    payload = parse_stats_json(out.read_text(encoding="utf-8"))
    stamped = payload.get("write_token")
    if stamped is None:
        raise TaskFailure(
            f"{out} carries no write_token, so it cannot be shown to be this run's measurement.\n"
            f"  {argv[0]} was given --write-token {token}. A payload without one was written "
            f"either by an older copy of this module - check module_path={module_path!r} - or by "
            f"an earlier run whose file was restored to this path after it was deleted. The "
            f"numbers are refused rather than reported under parameters they may not describe.\n"
            f"  log: {log}")
    if stamped != token:
        raise TaskFailure(
            f"{out} carries write_token {stamped!r}, not this run's {token!r}, so it is another "
            f"invocation's measurement standing at this path.\n"
            f"  The path was deleted and observed absent immediately before the command ran; what "
            f"is here now was restored or copied rather than written. Its numbers describe a "
            f"different matrix or different parameters.\n"
            f"  log: {log}")

    metrics = dict(payload["stats"])
    metrics["source_format"] = payload["format"]
    if rank_points:
        rank = payload.get("barcode_rank")
        if rank is None:
            raise TaskFailure(f"{rank_points} rank points were requested but {out} carries no "
                              f"'barcode_rank'")
        metrics["n_rank_points"] = len(rank)
    return {"outputs": [out], "metrics": metrics, "versions": dict(payload["versions"])}


# ---------------------------------------------------------------------------------------------
# subprocess entry point


def _observed_versions() -> dict:
    """Versions of the libraries this process imported. Nothing is inferred from an environment.

    What is recorded is IMPORTED, which is a weaker claim than USED and is stated as the weaker
    one: anndata pulls h5py and pandas in transitively, so an MTX directory read without ever
    touching HDF5 still reports an h5py version. `NOT_INVOKED` therefore means the library was
    never imported at all in this process, not that its code did not run.

    A library that was not imported is recorded as not invoked rather than omitted or blanked, so
    a reader can tell that from "it was imported and its version was not captured". Every string
    here was obtained by asking the module in the process that read the matrix - never from a
    lockfile, an environment name, or the orchestrator's own interpreter, which on a scheduler is
    a different machine with a different stack.
    """
    out = {"python": sys.version.split()[0]}
    for mod in ("anndata", "scanpy", "numpy", "scipy", "h5py", "pandas"):
        if mod in sys.modules:
            v = getattr(sys.modules[mod], "__version__", None)
            out[mod] = str(v) if v is not None else "imported, version not exposed"
        else:
            out[mod] = NOT_INVOKED
    return out


def main(argv=None) -> int:
    """Measure one matrix and write the JSON payload `parse_stats_json()` reads.

    This exists so the measurement can happen in the environment that has anndata, which is
    usually not the environment running the pipeline.
    """
    ap = argparse.ArgumentParser(
        prog="adapters/matrix.py",
        description="Measure a count matrix: the statistics step 0 renders a verdict on.")
    sub = ap.add_subparsers(dest="cmd")
    s = sub.add_parser("stats", help="write summary statistics for one matrix as JSON")
    s.add_argument("--matrix", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--expected-genes", dest="expected_genes", type=int)
    s.add_argument("--name")
    s.add_argument("--tmp-dir", dest="tmp_dir")
    s.add_argument("--rank-points", dest="rank_points", type=int, default=0)
    s.add_argument("--reuse-extracted", dest="reuse_extracted", action="store_true")
    # Echoed into the payload so the caller can tell this invocation's file from one restored to
    # the same path. It is not a measurement and nothing here reads it.
    s.add_argument("--write-token", dest="write_token")
    a = ap.parse_args(argv)
    if a.cmd != "stats":
        ap.print_help()
        return 2

    try:
        kind = detect_format(a.matrix)
        adata = read_matrix(a.matrix, tmp_dir=a.tmp_dir, reuse_extracted=a.reuse_extracted)
        stats = summary_stats(adata, expected_genes=a.expected_genes,
                              name=a.name if a.name else str(a.matrix))
        payload = {
            "schema": STATS_SCHEMA,
            "source": str(a.matrix),
            "format": kind,
            "stats": stats,
            "versions": _observed_versions(),
        }
        if a.write_token is not None:
            payload["write_token"] = str(a.write_token)
        if a.rank_points:
            payload["barcode_rank"] = [[r, t] for r, t in barcode_rank(adata, n=a.rank_points)]
        out = Path(a.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        # allow_nan=False: a NaN that reached here would be written as the JSON literal NaN,
        # read back as a float, and compared against a threshold it can only lose against.
        out.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False),
                       encoding="utf-8")
    except TaskFailure as e:
        sys.stderr.write(f"matrix adapter: {e}\n")
        return 1
    except ValueError as e:
        sys.stderr.write(f"matrix adapter: refusing to write a non-finite statistic: {e}\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
