# Execution adapter for step 4: it exports a matrix, runs scDblFinder, and reports what came
# back. It removes nothing, it chooses no threshold, and it never converts a barcode that was
# not scored into a barcode that was scored and cleared.
"""Step 4's execution layer - scDblFinder through Rscript, and the bookkeeping that must survive
the round trip.

WHAT THIS FILE DECIDES: NOTHING

`dbr` is DECLARED by the platform, `dbr.sd` is DERIVED from a sweep, and the light floor is
step 3's. All three arrive here as arguments and none of them has a default. A module that
supplies a default for an adjudicated parameter is reproducible only in the sense that it does
the same wrong thing every time, so the failure when one is missing is a refusal rather than a
fallback.

THE ONE THING THIS ADAPTER MUST NOT LOSE

A barcode below the light floor was never handed to the detector. Its doublet status is UNKNOWN.
It is not a singlet, and the difference is not cosmetic: in the calibration cohort the
never-scored set was 59,322 nuclei, 24% of the total, and folding it into the negatives moves
the reported rate from 7.52% of nuclei scored to 5.70% of all cells. Both numbers are true and
only one of them answers any given question.

So the calls are not a plain dict. `DoubletCalls` holds the scored barcodes as
`{barcode: (score, is_doublet)}` and carries the never-scored set beside them, and a lookup of a
never-scored barcode raises `UnknownDoubletStatus` instead of returning anything at all -
including instead of returning a caller-supplied default. That is deliberate and it is the whole
point of the class: `calls.get(bc, False)` is precisely the defect docs/PRINCIPLES.md section 4
records, where an unknown doublet fraction read through `or 0` counted every cell in a cluster
as surviving. `UnknownDoubletStatus` inherits from `LookupError` rather than `KeyError`, so a
caller wrapping the lookup in `except KeyError` does not quietly re-acquire the bug.

WHY THE MATRIX GOES OUT AS MatrixMarket

R is not asked to read .h5ad. A triple of matrix.mtx, barcodes.tsv and features.tsv is a format
both runtimes read with what they already have - Matrix::readMM on one side, scipy.io on the
other - and it makes the orientation explicit rather than implicit. The export is genes x cells,
which is AnnData's transpose, and `adapters/scdblfinder.R` checks the dimensions against the two
label files before it computes anything: a transposed matrix is still a valid matrix and still
returns scores.

Genes are never dropped on export, all-zero ones included. Dropping a gene is a removal under
the checklist, scDblFinder performs its own feature selection regardless, and an empty gene
costs one row of nothing.

A PREVIOUS RUN'S CALLS ARE NOT THIS RUN'S CALLS

`run_scdblfinder` used to check that `out_csv` existed after Rscript returned. A file left by an
earlier run satisfies that, and satisfies every cross-check built on top of it - the barcodes
match, because the export is the same; the doublet count matches, because it is compared against
the R metrics printed by a run that DID happen; and the whole sweep then reports the previous
setting's calls under the new one, which is the one comparison the sweep exists to make. The
output is therefore DELETED before Rscript is launched, and the R adapter deletes it too, so that
what is at the path afterwards can only have come from this invocation.

UNKNOWN IS NOT A VALUE, AND IT HAS MORE THAN TWO SHAPES

`is_missing()` is the single predicate this module asks that question with, and it covers None,
float nan, numpy nan, `pandas.NA`, `pandas.NaT` and `numpy.ma.masked`. A blank `dbr` that reaches
`_require_number` as `pandas.NA` is not None and is not a float; unguarded it survives to a
comparison, compares False against every bound, and becomes a threshold nobody chose.

WHAT IS SEPARATELY TESTABLE

The two halves that can be wrong without the tool being present are argument construction and
output parsing, so both are pure functions that touch no file: `resolve_dbr_sd`,
`build_rscript_cmd`, `parse_versions`, `parse_r_metrics`, `parse_calls_rows`,
`parse_calls_csv_text`, `parse_mtx_header` and `deep_decile_rate`.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import re
import sys
import time
from collections import abc as _abc
from dataclasses import dataclass, field as _dc_field
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

# engine.task is standard-library only, so importing it here keeps this module importable on a
# bare clone. The path bootstrap exists because an adapter is sometimes loaded by file path - by
# a test, or by a runner that has not put the repository root on sys.path - and a locally
# defined stand-in for TaskFailure would break `except TaskFailure` in the orchestrator.
try:  # pragma: no cover - whichever import path is in use is the one exercised
    from engine.task import TaskFailure
except ImportError:  # pragma: no cover
    _ROOT = str(Path(__file__).resolve().parents[1])
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    from engine.task import TaskFailure

__all__ = [
    "UNKNOWN",
    "is_missing",
    "is_true",
    "SCDBLFINDER_MIN_UMI",
    "XGBOOST_SHIM_FROM",
    "VERSION_PREFIX",
    "METRIC_PREFIX",
    "CALLS_HEADER",
    "REQUIRED_VERSIONS",
    "UnknownDoubletStatus",
    "DoubletCalls",
    "ExportedMatrix",
    "ScDblFinderDetector",
    "resolve_dbr_sd",
    "build_rscript_cmd",
    "parse_versions",
    "parse_r_metrics",
    "parse_calls_rows",
    "parse_calls_csv_text",
    "parse_mtx_header",
    "read_calls",
    "read_labels_file",
    "deep_decile_rate",
    "export_matrix",
    "run_scdblfinder",
    "sweep",
]

#: The doublet status of a barcode that was never handed to the detector. Not a call.
UNKNOWN = "UNKNOWN"

#: scDblFinder's own documented minimum coverage: "it might be necessary to remove cells with a
#: very low coverage (e.g. <200 reads) to avoid errors." Recorded so this adapter's declared
#: contract matches modules/03_light_floor and modules/04_doublets. It is applied as a default
#: nowhere in this file; the floor arrives as an argument or not at all.
SCDBLFINDER_MIN_UMI = 200

#: The xgboost version from which scDblFinder's top-level training arguments are accepted only
#: through a deprecation shim. See conf/env/install_rdoublet.sh: scds crashes on the same
#: change, scDblFinder does not, and the scores it then returns look entirely normal.
XGBOOST_SHIM_FROM = "2.0.0"

VERSION_PREFIX = "##scqc-version"
METRIC_PREFIX = "##scqc-metric"

#: Column order written by adapters/scdblfinder.R.
CALLS_HEADER = ("barcode", "doublet_score", "doublet_class")

#: Versions the R script must have printed. A tool that ran and whose version was not captured
#: is not the same as a tool that did not run, and only the second has a recorded form
#: (engine.provenance.NOT_INVOKED). The first is a failure.
REQUIRED_VERSIONS = ("R", "scDblFinder")

_NOT_A_NUMBER = frozenset({"", "na", "n/a", "nan", "none", "null", "nil", "-", "."})
_SAFE_LABEL = re.compile(r"[^A-Za-z0-9._-]+")
_BOM = "﻿"
_MISSING = object()


# --------------------------------------------------------------------- unknown is not a value

#: Type NAMES of the missing-value scalars that are neither None nor a float. Matched by name so
#: that the check costs nothing, and works, when pandas and numpy are not installed - this module
#: must stay importable with no third-party package at module scope.
_MISSING_TYPE_NAMES = frozenset({"NAType", "NaTType", "MaskedConstant"})


def is_missing(value) -> bool:
    """True when a value carries no information, in every shape one actually arrives in.

    The ONE predicate this module asks that question with. `None` is the shape everyone remembers
    and the rarest in practice: a blank cell in a samplesheet or a parameter table is
    `float('nan')` when pandas read it through numpy, `pandas.NA` through the nullable or
    pyarrow-backed dtypes, `pandas.NaT` for a parsed date, and `numpy.ma.masked` through a masked
    array. None of those is `None`, only the first is a `float`, and every one of them survives
    an `is not None` guard and then compares False against every bound - which downstream is
    indistinguishable from a value that was measured and did not exceed it. For `dbr` or `dbr.sd`
    that is a threshold nobody chose, applied to a whole cohort and reported as though declared.

    Four routes, cheapest first:

      * identity against `None`;
      * blank and whitespace-only text, in `str` and in `bytes`;
      * the type NAME, which catches `pandas.NA`, `pandas.NaT` and `numpy.ma.masked` without
        importing anything;
      * `value != value`, which catches `float('nan')`, `numpy.float64('nan')`,
        `numpy.float32('nan')` - not a `float` subclass - and `numpy.datetime64('NaT')`.

    pandas is then consulted, but only if it is ALREADY in `sys.modules`: a value cannot be a
    pandas scalar in a process that never imported pandas, so looking there is both sufficient
    and free, whereas importing pandas inside a predicate that runs per barcode would cost a
    second of start-up to a CLI that does not otherwise need it.

    For a value that may be a numpy boolean the rule here is `bool(x)` AFTER `is_missing(x)` has
    been checked, never `x is True`: `numpy.bool_(True) is True` is False, so identity reads a
    genuinely flagged row as unflagged. See `is_true`.
    """
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip()
    if isinstance(value, bytes):
        return not value.decode("utf-8", "replace").strip()
    if type(value).__name__ in _MISSING_TYPE_NAMES:
        return True
    try:
        if bool(value != value):
            return True
    except (TypeError, ValueError):
        pass                      # an object whose __ne__ refuses is not thereby missing
    pandas = sys.modules.get("pandas")
    if pandas is not None and not hasattr(value, "__len__"):
        try:
            verdict = pandas.isna(value)
        except (TypeError, ValueError):
            return False
        if not hasattr(verdict, "__len__"):
            try:
                return bool(verdict)
            except (TypeError, ValueError):
                return False
    return False


def is_true(value, name: str = "flag") -> bool:
    """Truthiness of a flag that may be a numpy boolean, with unknown refused rather than read.

    `numpy.bool_(True) is True` is False and `bool(pandas.NA)` raises, so the rule is: refuse
    unknown, then `bool()`.
    """
    if is_missing(value):
        raise TaskFailure(
            f"{name} is {value!r}, which is not True or False. A flag that was never set is not "
            f"the same as one set to False, and reading it as False here would record a decision "
            f"nobody made.")
    return bool(value)


# --------------------------------------------------------------------------- argument guards


def _require_number(value, name: str, *, minimum: float | None = None,
                    maximum: float | None = None, exclusive: bool = False) -> float:
    """Coerce to a finite float, refusing every shape a missing value arrives in.

    Everything `is_missing()` names stops here rather than becoming a parameter nobody chose:
    None, float nan, `pandas.NA`, `pandas.NaT`, `numpy.ma.masked`, and the blank or `'NA'` string
    a standard-library CSV reader produces. Which of those a blank table cell becomes depends
    only on the dtype backend the caller read it with, and only the first two were caught before.

    Numeric scalars that are NOT Python numbers are accepted, once past the unknown check:
    `numpy.int64` and `numpy.float32` are not subclasses of `int` or `float`, so an `isinstance`
    gate refuses an ordinary depth read out of a numpy array - a check that fires on correct
    behaviour, which is how checks end up removed.
    """
    if is_missing(value):
        raise TaskFailure(
            f"{name} is {value!r}, which carries no value. There is no default for it in this "
            f"adapter: a missing value must stop the run rather than become a threshold nobody "
            f"chose. None, NaN, pandas.NA, pandas.NaT, numpy.ma.masked and a blank cell all "
            f"arrive here, and not one of them is a number.")
    if isinstance(value, bool):
        raise TaskFailure(f"{name} is a boolean ({value!r}); a number was required.")
    if isinstance(value, str):
        text = value.strip()
        if text.lower() in _NOT_A_NUMBER:
            raise TaskFailure(
                f"{name} is {value!r}, which is a missing value written down rather than a "
                f"number. Unknown is not a value; supply {name} or do not run this step.")
        try:
            number = float(text)
        except ValueError:
            raise TaskFailure(f"{name} is not a number: {value!r}") from None
    elif isinstance(value, (int, float)):
        number = float(value)
    else:
        try:
            number = float(value)      # numpy scalars, and anything else carrying __float__
        except (TypeError, ValueError):
            raise TaskFailure(
                f"{name} is a {type(value).__name__}; a number was required.") from None
    if math.isnan(number):
        raise TaskFailure(
            f"{name} is NaN. NaN is not None and it is not a number - it passes an "
            f"`is not None` guard and then compares False against every threshold, which reads "
            f"downstream as permissive.")
    if not math.isfinite(number):
        raise TaskFailure(f"{name} is {number}, which is not finite.")
    if minimum is not None:
        if (number <= minimum) if exclusive else (number < minimum):
            raise TaskFailure(
                f"{name} must be {'greater than' if exclusive else 'at least'} {minimum}, "
                f"got {number}.")
    if maximum is not None:
        if (number >= maximum) if exclusive else (number > maximum):
            raise TaskFailure(
                f"{name} must be {'less than' if exclusive else 'at most'} {maximum}, "
                f"got {number}.")
    return number


def _require_int(value, name: str, *, minimum: int | None = None) -> int:
    """`_require_number`, then refuse anything carrying a fractional part."""
    number = _require_number(value, name, minimum=minimum)
    if number != int(number):
        raise TaskFailure(f"{name} must be a whole number, got {value!r}.")
    return int(number)


def _require_label(value, name: str) -> str:
    """A non-empty string that is not a missing value spelled out.

    `is_missing` runs first so that the refusal names the real problem. A barcode column read
    through pandas' nullable dtypes hands blanks over as `pandas.NA`, not as `''`, and
    `str(pandas.NA)` is the entirely plausible `'<NA>'` - which would have become a key.
    """
    if is_missing(value):
        raise TaskFailure(
            f"{name} is {value!r}, which carries no value. None, NaN, pandas.NA, pandas.NaT, "
            f"numpy.ma.masked and a blank cell all arrive here; each stringifies into something "
            f"that looks like a real label ('nan', '<NA>', 'NaT') and would become a key.")
    if not isinstance(value, str):
        raise TaskFailure(f"{name} is a {type(value).__name__}; a string was required.")
    text = value.strip()
    if not text or text.lower() in _NOT_A_NUMBER:
        raise TaskFailure(
            f"{name} is {value!r}, which is a blank or a missing value written down. A barcode "
            f"read out of an empty table cell arrives as 'nan' and looks like a real key.")
    return text


def _require_barcodes(values, name: str) -> tuple:
    """A sequence of distinct, non-blank barcode strings, order preserved."""
    if values is None:
        raise TaskFailure(f"{name} was not supplied.")
    if isinstance(values, (str, bytes)):
        raise TaskFailure(f"{name} is a single string; a sequence of barcodes was required.")
    out = []
    seen: dict = {}
    for i, v in enumerate(values):
        barcode = _require_label(v, f"{name}[{i}]")
        if barcode in seen:
            raise TaskFailure(
                f"{name} contains {barcode!r} at positions {seen[barcode]} and {i}. A duplicated "
                f"barcode merges two nuclei into one key, and the second call silently replaces "
                f"the first.")
        seen[barcode] = i
        out.append(barcode)
    if not out:
        raise TaskFailure(f"{name} is empty; there is nothing to score.")
    return tuple(out)


def _safe_label(text) -> str:
    """Filename-safe form of a label, so `dbr.sd = 1` cannot collide with `default`."""
    return _SAFE_LABEL.sub("_", str(text).strip()) or "unnamed"


def _package_version(module, distribution: str) -> str:
    """The version of an imported package, asked of the installed distribution first.

    `module.__version__` is not universal and is actively being withdrawn - anndata 0.12 raises
    a FutureWarning on it - so the installed metadata is preferred and the attribute is the
    fallback. A version that could not be obtained is reported as `unreported`, which is
    distinguishable from a version string; it is never guessed from a requirements file, because
    a recorded version that was not observed cannot be told apart from a real one.
    """
    from importlib import metadata

    try:
        return str(metadata.version(distribution))
    except Exception:  # noqa: BLE001 - any metadata failure falls back to the attribute
        attr = getattr(module, "__version__", None)
        return "unreported" if attr is None else str(attr)


# --------------------------------------------------------------------------- the calls object


class UnknownDoubletStatus(LookupError):
    """Asked for the doublet status of a barcode that was never scored.

    Deliberately a `LookupError` and not a `KeyError`. A caller writing
    `try: calls[bc] except KeyError: singlet` has re-created the exact defect this class exists
    to prevent, and inheriting from `KeyError` would let that code swallow the refusal without
    anyone noticing it had.
    """


class DoubletCalls(dict):
    """`{barcode: (score, is_doublet)}` for the barcodes that were SCORED, and nothing else.

    The never-scored barcodes are carried in `.unscored` and are absent from the mapping, so
    iterating the calls can never yield an invented one. Looking one up raises rather than
    returning a value, and `get()` raises rather than returning a default, because supplying a
    default for a barcode known never to have been examined is the failure in
    docs/PRINCIPLES.md section 4 rather than a defence against it.

    `.unscored` is `None` when the never-scored population was not supplied. That third state is
    kept distinct from an empty set: "nothing was left unscored" and "nobody said what was left
    unscored" support different claims, and only the first supports a rate against the whole
    deliverable.
    """

    def __init__(self, scored: Mapping, unscored=None, *, sample: str | None = None,
                 source=None):
        super().__init__(scored)
        if unscored is None:
            self.unscored = None
        else:
            self.unscored = frozenset(
                _require_label(b, "unscored barcode") for b in unscored)
            overlap = self.unscored.intersection(self.keys())
            if overlap:
                raise TaskFailure(
                    f"{len(overlap)} barcode(s) are recorded both as scored and as never "
                    f"scored, e.g. {sorted(overlap)[:5]}. One of the two records is wrong and "
                    f"there is no safe way to pick which.")
        self.sample = sample
        self.source = None if source is None else str(source)

    # -- lookups that refuse rather than invent ----------------------------------------------

    def __missing__(self, barcode):
        raise UnknownDoubletStatus(self._explain(barcode))

    def get(self, barcode, default=_MISSING):
        """Refuses a default. `calls.get(bc, False)` is the bug, not the guard against it."""
        if dict.__contains__(self, barcode):
            return dict.__getitem__(self, barcode)
        raise UnknownDoubletStatus(
            self._explain(barcode)
            + ("" if default is _MISSING else
               f" A default of {default!r} was supplied and is ignored: a default for a barcode "
               f"nobody examined records an answer that was never obtained."))

    def _explain(self, barcode) -> str:
        where = f" in sample {self.sample}" if self.sample else ""
        if self.unscored is None:
            return (f"barcode {barcode!r} was not scored{where}, and the never-scored population "
                    f"was not supplied to these calls, so it cannot be said whether it sat below "
                    f"the light floor or was absent from the sample. Either way its doublet "
                    f"status is UNKNOWN, which is not a singlet.")
        if barcode in self.unscored:
            return (f"barcode {barcode!r} was never scored{where} - it sat below the light floor "
                    f"and was never handed to the detector. Its doublet status is UNKNOWN, which "
                    f"is not the same as a singlet, and folding it into the negatives moves "
                    f"every rate computed from these calls.")
        return (f"barcode {barcode!r} is not part of these calls{where}: it is neither among the "
                f"{len(self)} scored nor among the {len(self.unscored)} recorded as never "
                f"scored.")

    # -- accessors ---------------------------------------------------------------------------

    def status(self, barcode) -> str:
        """`"doublet"`, `"singlet"` or `"UNKNOWN"`; raises for a barcode from another sample."""
        if dict.__contains__(self, barcode):
            return "doublet" if dict.__getitem__(self, barcode)[1] else "singlet"
        if self.unscored is not None and barcode in self.unscored:
            return UNKNOWN
        raise UnknownDoubletStatus(self._explain(barcode))

    def score_of(self, barcode) -> float:
        return self[barcode][0]

    def is_doublet(self, barcode) -> bool:
        return self[barcode][1]

    def known(self, barcode) -> bool:
        """True if these calls have anything at all to say about this barcode."""
        return (dict.__contains__(self, barcode)
                or (self.unscored is not None and barcode in self.unscored))

    @property
    def n_scored(self) -> int:
        return len(self)

    @property
    def n_unscored(self) -> int | None:
        """None when the never-scored population was not supplied. Never 0 as a stand-in."""
        return None if self.unscored is None else len(self.unscored)

    @property
    def n_called(self) -> int:
        return sum(1 for v in self.values() if v[1])

    def called(self) -> tuple:
        return tuple(bc for bc, v in self.items() if v[1])

    def fraction_called(self, denominator: str) -> float:
        """A doublet rate always carries its denominator; there is no bare-number form.

        `"scored"` divides by the nuclei actually examined and is what
        `modules/04_doublets.SweepResult.per_sample_rate` is documented to hold. `"all"` divides
        by scored plus never-scored and needs that population to have been supplied.
        """
        if denominator == "scored":
            if not self:
                raise TaskFailure(
                    "no nuclei were scored, so a rate over the scored set has no denominator. "
                    "Reported as undefined rather than as zero: a zero denominator is a missing "
                    "measurement, not a rate of nothing.")
            return self.n_called / len(self)
        if denominator == "all":
            if self.unscored is None:
                raise TaskFailure(
                    "a rate against all cells needs the never-scored population, which was not "
                    "supplied. Pass `unscored=` from the light-floor step, or ask for the rate "
                    "over the scored set instead and say which it is wherever it is printed.")
            total = len(self) + len(self.unscored)
            if total == 0:
                raise TaskFailure("no nuclei at all; the rate has no denominator.")
            return self.n_called / total
        raise TaskFailure(
            f"denominator must be 'scored' or 'all', got {denominator!r} - a rate without one "
            f"is not comparable to anything, including a published figure.")

    def __repr__(self) -> str:
        un = "unsupplied" if self.unscored is None else f"{len(self.unscored):,}"
        return (f"<DoubletCalls sample={self.sample!r} scored={len(self):,} "
                f"called={self.n_called:,} unscored={un}>")


# --------------------------------------------------------------------------- pure: arguments


def resolve_dbr_sd(setting, dbr) -> tuple:
    """Turn one sweep token into `(label, value_or_None)`. Pure.

    `modules/04_doublets.SWEEP` is `("default", "dbr", "1")`, and those tokens are what
    `SweepResult.setting` and `recommend()` match on, so the label is returned verbatim and
    round-trips unchanged. Only the value is resolved:

      "default"   the argument is not passed at all and the installed package applies its own
      "dbr"       dbr.sd is set equal to the declared dbr
      otherwise   a number, or a string holding one

    `dbr` is validated only where it is used, so a sweep can be planned before the platform's
    rate is known and will refuse at the point the missing value would have mattered.

    `None` means "do not pass the argument" and is the only value that means it. A blank cell -
    NaN, `pandas.NA`, `pandas.NaT`, `numpy.ma.masked` - is refused rather than read as `None`,
    because a sweep setting nobody filled in is not a decision to use the package default; it is
    a row of the sweep plan that was never written, and silently running the default under a
    label that says something else is how a sweep stops measuring what it claims to.
    """
    if setting is None:
        return ("default", None)
    if is_missing(setting):
        raise TaskFailure(
            f"dbr.sd setting is {setting!r}, which carries no value. Pass the token 'default' to "
            f"omit the argument deliberately - a blank is not a decision to use the package "
            f"default, and this sweep's whole purpose is to distinguish the two.")
    if isinstance(setting, bool):
        raise TaskFailure(
            f"dbr.sd setting is a boolean ({setting!r}); a token or a number was required.")
    if isinstance(setting, str):
        token = setting.strip()
        if token.lower() == "default":
            return ("default", None)
        if token.lower() == "dbr":
            return ("dbr", _require_number(dbr, "dbr", minimum=0.0, maximum=1.0, exclusive=True))
        label = token
    elif isinstance(setting, float):
        label = repr(setting)
    elif isinstance(setting, int):
        label = str(setting)
    else:
        raise TaskFailure(
            f"dbr.sd setting is a {type(setting).__name__}; expected 'default', 'dbr' or a "
            f"number.")
    return (label, _require_number(setting, "dbr.sd", minimum=0.0))


def build_rscript_cmd(rscript, script, mtx_dir, out_csv, dbr, dbr_sd, seed, *,
                      threads: int = 1, refuse_xgboost_ge: str | None = XGBOOST_SHIM_FROM,
                      features_col: int = 1) -> list:
    """Build the exact argv for `adapters/scdblfinder.R`. Pure: no file is touched.

    Separated from `run_scdblfinder` because the two failures this call can carry - a swapped
    argument and a silently defaulted one - are both visible in the argv, and neither needs R
    installed to test for. `dbr_sd` of None becomes the token `default`, which tells the R
    script to omit the argument entirely rather than to pass a number chosen here.
    """
    rscript_s = _require_label(str(rscript), "rscript")
    script_s = _require_label(str(script), "script")
    mtx_s = _require_label(str(mtx_dir), "mtx_dir")
    out_s = _require_label(str(out_csv), "out_csv")
    dbr_v = _require_number(dbr, "dbr", minimum=0.0, maximum=1.0, exclusive=True)
    seed_v = _require_int(seed, "seed")
    threads_v = _require_int(threads, "threads", minimum=1)
    features_v = _require_int(features_col, "features_col", minimum=1)
    if dbr_sd is None:
        sd_token = "default"
    else:
        sd_token = "{0:.10g}".format(_require_number(dbr_sd, "dbr.sd", minimum=0.0))
    cmd = [rscript_s, script_s, mtx_s, out_s, "{0:.10g}".format(dbr_v), sd_token, str(seed_v),
           "threads={0}".format(threads_v), "features_col={0}".format(features_v)]
    if refuse_xgboost_ge is None:
        # Written into the argv, and so into the log, rather than left implicit: scoring through
        # the xgboost deprecation shim is permitted and has to be visible afterwards.
        cmd.append("xgboost_max=")
    else:
        cmd.append("xgboost_max=" + _require_label(refuse_xgboost_ge, "refuse_xgboost_ge"))
    return cmd


# --------------------------------------------------------------------------- pure: parsing


def _tagged_lines(text, prefix: str, what: str) -> dict:
    if text is None:
        raise TaskFailure(f"no output was captured, so no {what} could be read from it.")
    out: dict = {}
    for line in str(text).splitlines():
        if not line.startswith(prefix):
            continue
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) != 3:
            raise TaskFailure(
                f"malformed {what} line from the R adapter: {line!r}. Expected "
                f"{prefix}<TAB>name<TAB>value.")
        name, value = parts[1].strip(), parts[2].strip()
        if not name:
            raise TaskFailure(f"a {what} line has an empty name: {line!r}")
        if name in out and out[name] != value:
            raise TaskFailure(
                f"the R adapter reported {name} twice with different values, {out[name]!r} then "
                f"{value!r}. One run cannot have used both.")
        out[name] = value
    return out


def parse_versions(text) -> dict:
    """Read the `##scqc-version` lines the R script printed. Pure.

    Every version here was obtained by asking the installed library on this run. Nothing is
    filled in from a lock file or an environment name: a recorded version that was not observed
    is a fabricated provenance record and cannot be told apart from a real one. A tool that ran
    and whose version could not be captured is therefore a failure, not `not invoked`.
    """
    versions = _tagged_lines(text, VERSION_PREFIX, "version")
    missing = [k for k in REQUIRED_VERSIONS if k not in versions]
    if missing:
        raise TaskFailure(
            f"the R adapter did not report {', '.join(missing)}. It prints "
            f"'{VERSION_PREFIX}<TAB>name<TAB>value' before it scores anything, so an absent line "
            f"means the output was truncated, redirected, or produced by a different script - "
            f"and a version that was not observed must not be written into a report.")
    return versions


def parse_r_metrics(text) -> dict:
    """Read the `##scqc-metric` lines the R script printed, as raw strings. Pure."""
    return _tagged_lines(text, METRIC_PREFIX, "metric")


def _r_metric(metrics: Mapping, key: str, where: str) -> str:
    """Fetch a metric the R script must have printed, refusing an absent one.

    Written out rather than reached with `.get(key, 0)`, because the number this returns is
    cross-checked against one derived independently, and a zero standing in for an absent value
    would make the two agree exactly when they should not.
    """
    if key not in metrics:
        raise TaskFailure(
            f"the R adapter did not report '{key}'. It is printed as "
            f"'{METRIC_PREFIX}<TAB>{key}<TAB>value' on every successful run, so its absence "
            f"means the output was truncated or came from a different script. Log: {where}")
    return metrics[key]


def parse_calls_rows(rows: Iterable, *, source: str = "<rows>") -> dict:
    """Turn CSV rows into `{barcode: (score, is_doublet)}`. Pure, and strict on purpose.

    Nothing is coerced or repaired. An empty score is not zero, an `NA` class is not a singlet,
    and an unrecognised class is not mapped onto the nearest known one - all three raise,
    because each silently converts a nucleus nobody could classify into one the detector
    cleared.
    """
    iterator = iter(rows)
    try:
        header = next(iterator)
    except StopIteration:
        raise TaskFailure(
            f"{source} is empty; the R adapter writes a header even for an empty result.") \
            from None
    got = tuple(str(c).lstrip(_BOM).strip().lower() for c in header)
    if got != CALLS_HEADER:
        raise TaskFailure(
            f"{source} has header {got} but {CALLS_HEADER} was expected. Column order is not "
            f"guessed: a barcode column read as a score is silently numeric-looking nonsense.")
    calls: dict = {}
    for lineno, row in enumerate(iterator, start=2):
        if not row or (len(row) == 1 and not str(row[0]).strip()):
            continue
        if len(row) != 3:
            raise TaskFailure(
                f"{source} line {lineno} has {len(row)} field(s), expected 3: {row!r}. A barcode "
                f"containing a comma produces exactly this, and shifts every column after it.")
        barcode = _require_label(row[0], f"{source} line {lineno} barcode")
        if barcode in calls:
            raise TaskFailure(
                f"{source} line {lineno} repeats barcode {barcode!r}. The second row would "
                f"replace the first and one nucleus's call would vanish without a count.")
        score = _require_number(row[1], f"{source} line {lineno} doublet_score")
        klass = str(row[2]).strip().lower()
        if klass not in ("singlet", "doublet"):
            raise TaskFailure(
                f"{source} line {lineno} has doublet_class {row[2]!r}; only 'singlet' and "
                f"'doublet' are accepted. An unrecognised label is UNKNOWN and is not mapped "
                f"onto either of them.")
        calls[barcode] = (score, klass == "doublet")
    if not calls:
        raise TaskFailure(
            f"{source} holds a header and no calls. An empty result is not an absence of "
            f"doublets; it is a detector that scored nothing.")
    return calls


def parse_calls_csv_text(text, *, source: str = "<text>") -> dict:
    """`parse_calls_rows` over CSV text. Pure."""
    if text is None:
        raise TaskFailure(f"{source}: no text to parse.")
    return parse_calls_rows(csv.reader(io.StringIO(str(text))), source=source)


def parse_mtx_header(lines: Sequence) -> dict:
    """Read a MatrixMarket banner and size line. Pure.

    Used to verify an export after writing it. Orientation is the one property of this hand-off
    that cannot be recovered downstream - a transposed matrix scores without complaint - so the
    dimensions are read back out of the file rather than assumed from the object that wrote it.
    """
    banner = None
    for raw in lines:
        text = str(raw).strip()
        if not text:
            continue
        if banner is None:
            banner = text
            continue
        if text.startswith("%"):
            continue
        parts = text.split()
        if len(parts) != 3:
            raise TaskFailure(f"MatrixMarket size line is {text!r}; expected three integers.")
        try:
            n_rows, n_cols, n_nonzero = (int(p) for p in parts)
        except ValueError:
            raise TaskFailure(
                f"MatrixMarket size line is not three integers: {text!r}") from None
        fields = banner.split()
        if not banner.startswith("%%MatrixMarket") or len(fields) < 4:
            raise TaskFailure(f"not a MatrixMarket file; its first line is {banner!r}")
        if fields[2].lower() != "coordinate":
            raise TaskFailure(
                f"MatrixMarket format is {fields[2]!r}; the R adapter reads 'coordinate' "
                f"(sparse) only.")
        return {"n_rows": n_rows, "n_cols": n_cols, "n_nonzero": n_nonzero,
                "field": fields[3].lower(), "banner": banner}
    raise TaskFailure("no MatrixMarket size line found; the file is empty or truncated.")


# --------------------------------------------------------------------------- small file I/O


def _open_maybe_gz(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(str(path), "rt", encoding="utf-8", newline="")
    return open(str(path), "rt", encoding="utf-8", newline="")


def _resolve_optional_gz(directory, stem: str) -> Path:
    plain = Path(directory) / stem
    gz = Path(directory) / (stem + ".gz")
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise TaskFailure(
        f"neither {plain} nor {gz} exists. The export written by export_matrix() holds "
        f"matrix.mtx, barcodes.tsv and features.tsv, each optionally gzipped.")


def _clear_paths(paths, *, what: str) -> list:
    """Delete each path that exists, before the thing that writes it runs. Returns what went.

    A file left by an earlier run satisfies every check made AFTER a command returns - it exists,
    it parses, its dimensions agree with themselves - so "the output is present" proves only that
    an output is present. Emptying the paths first is what turns that check into a statement
    about THIS run. A file that cannot be removed stops the caller here rather than being worked
    around, because working around it means writing beside it and reading whichever one wins.
    """
    removed = []
    for p in paths:
        path = Path(p)
        if not path.exists():
            continue
        try:
            path.unlink()
        except OSError as exc:
            raise TaskFailure(
                f"could not remove {path} ({what}): {exc}. It is removed BEFORE the run rather "
                f"than checked afterwards, because a file that survives the command cannot be "
                f"told apart from one the command wrote, and the previous run's numbers would be "
                f"recorded under this run's parameters.") from None
        removed.append(str(path))
    return removed


def read_labels_file(path, *, column: int = 0, name: str = "labels") -> tuple:
    """Read one column of a TSV of labels, refusing blanks and duplicates."""
    path = Path(path)
    if not path.exists():
        raise TaskFailure(f"{name} file does not exist: {path}")
    out = []
    with _open_maybe_gz(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            text = line.rstrip("\r\n")
            if not text:
                continue
            parts = text.split("\t")
            if len(parts) <= column:
                raise TaskFailure(
                    f"{path} line {lineno} has {len(parts)} column(s); column {column + 1} was "
                    f"requested.")
            out.append(parts[column])
    return _require_barcodes(out, f"{name} in {path}")


def read_calls(csv_path, *, unscored=None, sample: str | None = None) -> DoubletCalls:
    """Read the R adapter's CSV back as `DoubletCalls`.

    `unscored` is the population that was never handed to the detector - normally the barcodes
    the light floor left out. Leaving it None is permitted and is recorded as "not supplied",
    which is a different claim from "nothing was left out" and supports fewer statements.
    """
    path = Path(csv_path)
    if not path.exists():
        raise TaskFailure(f"calls file does not exist: {path}")
    if path.stat().st_size == 0:
        raise TaskFailure(
            f"calls file is empty: {path}. The R adapter writes a .partial file and renames it, "
            f"so a zero-byte file at this path was not written by a completed run.")
    with _open_maybe_gz(path) as handle:
        calls = parse_calls_rows(csv.reader(handle), source=str(path))
    return DoubletCalls(calls, unscored, sample=sample, source=path)


# --------------------------------------------------------------------------- pure: deep decile


def deep_decile_rate(umi_by_key: Mapping, called_keys: Iterable, scored_keys: Iterable,
                     *, decile: float = 0.1) -> dict:
    """Share of the deepest UMI decile that was called, over the SCORED set only. Pure.

    `modules/04_doublets.recommend` compares this against `DEEP_DECILE_ALARM`, so the number
    that feeds `SweepResult.deep_decile_rate` is `rate_over_scored`. The key names its
    denominator because a doublet rate quoted without one is not comparable to anything,
    including a published figure.

    Never-scored nuclei are excluded from both the numerator and the denominator rather than
    counted as not called, which is the arithmetic that would turn a coverage gap into evidence
    of cleanliness. They are not merely dropped, though: `n_unscored_at_or_above_cut` reports how
    many of them sit inside the depth range the decile covers. A light floor removes the
    shallowest nuclei, so a non-zero count there means the scoring set was chosen by something
    other than depth and the decile is not the population it appears to be.

    Keys may be anything hashable. Pooling a cohort means keying by `(sample, barcode)`, because
    barcode sequences repeat across libraries and a bare barcode would silently merge them.
    """
    fraction = _require_number(decile, "decile", minimum=0.0, maximum=1.0, exclusive=True)
    scored = list(scored_keys)
    scored_set = set(scored)
    if len(scored_set) != len(scored):
        raise TaskFailure(
            "scored_keys contains duplicates; pool a cohort by (sample, barcode) rather than by "
            "barcode, which repeats across libraries.")
    if not scored_set:
        raise TaskFailure("no scored nuclei, so the deepest decile has no members.")
    called = set(called_keys)
    stray = called - scored_set
    if stray:
        raise TaskFailure(
            f"{len(stray)} barcode(s) are called doublets but are not in the scored set, e.g. "
            f"{sorted(map(str, stray))[:5]}. A call for a nucleus nobody scored means the two "
            f"records describe different populations.")
    missing = [k for k in scored_set if k not in umi_by_key]
    if missing:
        raise TaskFailure(
            f"{len(missing)} scored nucleus/nuclei have no UMI total, e.g. "
            f"{sorted(map(str, missing))[:5]}. A missing depth must not be read as zero: it "
            f"would place the nucleus in the shallowest decile and so out of this one.")
    depths = {k: _require_number(umi_by_key[k], f"UMI total for {k}", minimum=0.0)
              for k in scored_set}
    ordered = sorted(depths.items(), key=lambda kv: (-kv[1], str(kv[0])))
    n_top = max(1, int(math.ceil(fraction * len(ordered))))
    cut = ordered[n_top - 1][1]
    in_decile = [k for k, v in ordered if v >= cut]
    n_called_in = sum(1 for k in in_decile if k in called)
    above = 0
    for key, value in umi_by_key.items():
        if key in scored_set:
            continue
        if _require_number(value, f"UMI total for {key}", minimum=0.0) >= cut:
            above += 1
    return {
        "rate_over_scored": n_called_in / len(in_decile),
        "cut_umi": cut,
        "decile": fraction,
        "n_scored": len(scored_set),
        "n_in_decile": len(in_decile),
        "n_called_in_decile": n_called_in,
        "n_unscored_at_or_above_cut": above,
    }


# --------------------------------------------------------------------------- export


@dataclass
class ExportedMatrix:
    """Where the MatrixMarket triple went, and who is in it and who is not.

    `below_floor` and `not_selected` are two different reasons a barcode was never scored, and
    they are kept apart: one is the light floor doing its documented job, the other is the caller
    having handed in a smaller population than the object holds. Both end up UNKNOWN, and only
    the first is explained by a threshold.
    """

    mtx_dir: Path
    matrix_path: Path
    barcodes_path: Path
    features_path: Path
    exported: tuple
    below_floor: tuple
    not_selected: tuple
    umi_by_barcode: dict
    n_genes: int
    mtx_field: str
    min_umi: Optional[int]
    versions: dict = _dc_field(default_factory=dict)
    checks: dict = _dc_field(default_factory=dict)

    @property
    def unscored(self) -> tuple:
        """Every barcode in the source that was not exported, and so cannot have been scored."""
        return tuple(self.below_floor) + tuple(self.not_selected)

    def result(self) -> dict:
        """The adapter result shape, for an orchestrator that records this as its own step."""
        metrics = {
            "mtx_dir": str(self.mtx_dir),
            "n_exported": len(self.exported),
            "n_below_floor": len(self.below_floor),
            "n_not_selected": len(self.not_selected),
            "n_genes": self.n_genes,
            "mtx_field": self.mtx_field,
            "min_umi": self.min_umi,
        }
        metrics.update(self.checks)
        return {
            "outputs": [self.matrix_path, self.barcodes_path, self.features_path],
            "metrics": metrics,
            "versions": dict(self.versions),
        }


def export_matrix(h5ad_path, out_dir, min_umi, *, barcodes=None, layer=None,
                  gzip_mtx: bool = False, require_integer_counts: bool = True) -> ExportedMatrix:
    """Write a genes x cells MatrixMarket triple for `adapters/scdblfinder.R`.

    `out_dir` is normally a scratch directory: the triple is an inter-process format, not a
    deliverable, and nothing downstream reads it once the calls exist. It is an explicit
    argument rather than a temporary directory chosen here, because a path this adapter invents
    is a path the run log cannot show.

    `min_umi` has no default and must be passed explicitly, including as None to apply no floor
    at all. The light floor is step 3's DERIVED parameter, and a detector's scoring set is not
    something an adapter is entitled to choose; passing it in writing is what makes the choice
    reviewable.

    Every gene is exported, all-zero ones included. Dropping a gene is a removal under the
    checklist in docs/PRINCIPLES.md, scDblFinder performs its own feature selection regardless,
    and an empty gene costs one row of nothing.

    The whole matrix is read into memory. That is a real limit, and it is stated rather than
    worked around: the backed-mode alternative changes which values are read, and this adapter
    cannot verify that without the file in front of it.

    EVERY FILE THIS FUNCTION COULD WRITE IS REMOVED FIRST, in both its plain and its gzipped
    spelling. Writing `matrix.mtx.gz` beside a `matrix.mtx` left by an earlier export does not
    overwrite it, and both the R adapter and `_resolve_optional_gz` prefer the plain name - so
    the next run would score the PREVIOUS export's counts against this run's barcode list, with
    every dimension check passing because it is a complete and self-consistent matrix. It is
    simply the wrong one.
    """
    import anndata           # heavy; the CLI must stay importable without any of these
    import numpy as np
    import scipy.io
    import scipy.sparse as sp

    gzip_mtx = is_true(gzip_mtx, "gzip_mtx")
    require_integer_counts = is_true(require_integer_counts, "require_integer_counts")

    src = Path(h5ad_path)
    if not src.exists():
        raise TaskFailure(f"input matrix does not exist: {src}")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stale = _clear_paths(
        [out / n for n in ("matrix.mtx", "matrix.mtx.gz", "barcodes.tsv", "barcodes.tsv.gz",
                           "features.tsv", "features.tsv.gz")],
        what=f"a previous export in {out}")

    floor = None if min_umi is None else _require_int(min_umi, "min_umi", minimum=0)

    adata = anndata.read_h5ad(str(src))
    if layer is None:
        matrix = adata.X
        layer_name = "X"
    else:
        if layer not in adata.layers:
            raise TaskFailure(
                f"layer {layer!r} is not in {src}. Available: "
                f"{sorted(adata.layers.keys()) or 'none'}. scQC does not fall back to X - a "
                f"different layer is different data, and the report would describe the wrong "
                f"one.")
        matrix = adata.layers[layer]
        layer_name = f"layers[{layer!r}]"
    if matrix is None:
        raise TaskFailure(f"{src} has no counts in {layer_name}.")

    obs_names = _require_barcodes([str(b) for b in adata.obs_names], f"obs_names of {src}")
    var_names = [str(g) for g in adata.var_names]
    if not var_names:
        raise TaskFailure(f"{src} has no genes.")

    csr = matrix if sp.issparse(matrix) else sp.csr_matrix(np.asarray(matrix))
    csr = csr.tocsr()
    if tuple(csr.shape) != (len(obs_names), len(var_names)):
        raise TaskFailure(
            f"{layer_name} of {src} is {tuple(csr.shape)} but the object has {len(obs_names)} "
            f"obs and {len(var_names)} var. It is internally inconsistent and the export "
            f"orientation cannot be trusted.")

    data = csr.data
    if data.size and not bool(np.all(np.isfinite(data))):
        raise TaskFailure(
            f"{layer_name} of {src} holds non-finite values. NaN in a count matrix is a missing "
            f"measurement and it must not be written out as a number.")
    integral = bool(data.size == 0 or np.all(np.mod(data, 1) == 0))
    if require_integer_counts and not integral:
        raise TaskFailure(
            f"{layer_name} of {src} holds fractional values, so it is not raw counts. "
            f"scDblFinder builds its null by summing pairs of observed transcriptomes, and its "
            f"behaviour on non-count input is undocumented rather than loud. Point `layer=` at "
            f"the counts, or pass require_integer_counts=False to export a real-valued matrix "
            f"and record in the decisions file that the scores were computed on one.")

    totals = np.asarray(csr.sum(axis=1)).ravel()
    if totals.size and not bool(np.all(np.isfinite(totals))):
        raise TaskFailure(f"non-finite per-cell totals in {layer_name} of {src}.")
    umi_by_barcode = {bc: float(totals[i]) for i, bc in enumerate(obs_names)}
    index = {bc: i for i, bc in enumerate(obs_names)}

    if barcodes is None:
        wanted = list(obs_names)
        not_selected: list = []
    else:
        wanted = list(_require_barcodes(barcodes, "barcodes"))
        absent = [b for b in wanted if b not in index]
        if absent:
            raise TaskFailure(
                f"{len(absent)} requested barcode(s) are not in {src}, e.g. {absent[:5]}. "
                f"Dropping them here would remove nuclei the caller believes were scored; the "
                f"two populations have to be reconciled before scoring, not after.")
        chosen = set(wanted)
        not_selected = [b for b in obs_names if b not in chosen]

    if floor is None:
        exported = list(wanted)
        below: list = []
    else:
        exported = [b for b in wanted if umi_by_barcode[b] >= floor]
        below = [b for b in wanted if umi_by_barcode[b] < floor]
    if not exported:
        deepest = max((umi_by_barcode[b] for b in wanted), default=0.0)
        raise TaskFailure(
            f"no barcode in {src} reaches the floor of {floor} UMI ({len(wanted)} considered, "
            f"deepest {deepest:.0f}). An empty scoring set is a collapsed threshold, not a "
            f"sample without doublets.")
    empty = [b for b in exported if umi_by_barcode[b] <= 0]
    if empty:
        raise TaskFailure(
            f"{len(empty)} exported barcode(s) have zero total counts, e.g. {empty[:5]}. "
            f"scDblFinder's vignette records the consequence - 'Size factors should be "
            f"positive' - and that error arrives from inside normalisation, where it reads as a "
            f"bug in the tool. Apply the light floor (step 3) before scoring; it selects the "
            f"scoring set and removes nothing from the analysis.")

    rows = np.fromiter((index[b] for b in exported), dtype=np.int64, count=len(exported))
    sub = csr[rows, :]
    if integral:
        sub = sub.astype(np.int64)
        mm_field = "integer"
    else:
        mm_field = "real"
    # AnnData is cells x genes and the export is its transpose. Checked here, and again by the
    # R adapter against the two label files, because a transposed matrix scores without error.
    genes_by_cells = sub.T.tocsc()
    if tuple(genes_by_cells.shape) != (len(var_names), len(exported)):
        raise TaskFailure(
            f"transposed export is {tuple(genes_by_cells.shape)}, expected "
            f"({len(var_names)}, {len(exported)}).")

    matrix_path = out / ("matrix.mtx.gz" if gzip_mtx else "matrix.mtx")
    if gzip_mtx:
        with gzip.open(str(matrix_path), "wb") as handle:
            scipy.io.mmwrite(handle, genes_by_cells, field=mm_field, symmetry="general")
    else:
        with open(str(matrix_path), "wb") as handle:
            scipy.io.mmwrite(handle, genes_by_cells, field=mm_field, symmetry="general")

    barcodes_path = out / "barcodes.tsv"
    barcodes_path.write_text("".join(b + "\n" for b in exported), encoding="utf-8")
    features_path = out / "features.tsv"
    features_path.write_text("".join("{0}\t{0}\tGene Expression\n".format(g) for g in var_names),
                             encoding="utf-8")

    for path in (matrix_path, barcodes_path, features_path):
        if not path.exists() or path.stat().st_size == 0:
            raise TaskFailure(f"the export did not produce {path}; nothing downstream can run.")

    if gzip_mtx:
        with gzip.open(str(matrix_path), "rt", encoding="utf-8") as handle:
            head = [handle.readline() for _ in range(6)]
    else:
        with open(str(matrix_path), "rt", encoding="utf-8") as handle:
            head = [handle.readline() for _ in range(6)]
    header = parse_mtx_header(head)
    if (header["n_rows"], header["n_cols"]) != (len(var_names), len(exported)):
        raise TaskFailure(
            f"{matrix_path} declares {header['n_rows']} x {header['n_cols']} but "
            f"{len(var_names)} genes x {len(exported)} cells were written.")
    if header["n_nonzero"] != int(genes_by_cells.nnz):
        raise TaskFailure(
            f"{matrix_path} declares {header['n_nonzero']} non-zeros but the matrix holds "
            f"{int(genes_by_cells.nnz)}.")

    unscored_all = list(below) + list(not_selected)
    checks = {
        "source": str(src),
        "layer": layer_name,
        "counts_are_integral": integral,
        "n_nonzero": int(genes_by_cells.nnz),
        "n_duplicate_gene_names": len(var_names) - len(set(var_names)),
        "deepest_unscored_umi": (max(umi_by_barcode[b] for b in unscored_all)
                                 if unscored_all else None),
        "shallowest_exported_umi": min(umi_by_barcode[b] for b in exported),
        "gzip_mtx": gzip_mtx,
        "mtx_banner": header["banner"],
        "stale_export_files_removed": stale,
    }
    versions = {"anndata": _package_version(anndata, "anndata"),
                "scipy": _package_version(scipy, "scipy"),
                "numpy": _package_version(np, "numpy")}
    return ExportedMatrix(
        mtx_dir=out, matrix_path=matrix_path, barcodes_path=barcodes_path,
        features_path=features_path, exported=tuple(exported), below_floor=tuple(below),
        not_selected=tuple(not_selected), umi_by_barcode=umi_by_barcode,
        n_genes=len(var_names), mtx_field=mm_field, min_umi=floor, versions=versions,
        checks=checks)


# --------------------------------------------------------------------------- run


def run_scdblfinder(rscript, mtx_dir, out_csv, dbr, dbr_sd, seed, log, executor, *,
                    script=None, unscored=None, sample: str | None = None, threads: int = 1,
                    r_libs=None, refuse_xgboost_ge: str | None = XGBOOST_SHIM_FROM,
                    features_col: int = 1, timeout_s: int | None = None,
                    env: Mapping | None = None) -> dict:
    """Score one already-exported sample, verify what came back, and report it.

    Returns `{"outputs": [calls csv], "metrics": {...}, "versions": {...}}`. The calls
    themselves are read back with `read_calls(out_csv, unscored=...)`, which is a separate call
    so that the never-scored population travels with them explicitly rather than being carried
    implicitly in a result dict.

    `dbr_sd` of None passes the token `default`, which makes the R script omit the argument and
    lets the installed package apply its own. What that default IS is reported as
    `metrics["dbr_sd_formal_default"]`, together with `metrics["dbr_sd_formal_default_route"]`
    saying whether that string is the value the package uses or only the placeholder in its
    signature - see `adapters/scdblfinder.R`, which observes the installed function's formals and
    refuses to present a `NULL` placeholder as an observed default.

    THE OUTPUT IS DELETED BEFORE Rscript RUNS. Checking that `out_csv` exists afterwards proves a
    calls file is there and nothing else: a file left by an earlier setting of the sweep passes
    that check, then passes the barcode cross-check (the export is the same barcodes), then
    passes the doublet- and cell-count cross-checks (those compare the CSV against the R metrics
    that this run printed - and a stale CSV from a DIFFERENT dbr.sd will disagree, but only by
    luck of the counts differing). With the path emptied first, the file that is there afterwards
    was written by this invocation, which is what every check below then describes.

    The result is checked five ways before it is returned, because what this step can produce is
    a plausible-looking CSV rather than an error:

      * the CSV's barcodes must be exactly the exported barcodes, so a detector that silently
        dropped nuclei is caught here rather than becoming a coverage gap nobody counted;
      * the doublet count parsed from the CSV must equal the count R printed - one number
        reached by two independent routes;
      * so must the cell count;
      * the seed and dbr R reports having used must equal the ones asked for, which catches an
        argv that arrived shifted by one;
      * and so must dbr.sd. It was the one parameter never cross-checked, which is the wrong one
        to leave out: it is what the sweep VARIES and what step 4 DERIVES, so an argv that lost
        it - or an R release that ignored it - would produce a sweep in which every setting ran
        under the same value and the flat rate that follows would read as a property of the data.
        `default` and a number are also checked against each other in both directions, so a run
        that silently fell back to the package default cannot be recorded as an explicit value.
    """
    rscript_s = _require_label(str(rscript), "rscript")
    mtx = Path(mtx_dir)
    if not mtx.is_dir():
        raise TaskFailure(f"mtx_dir is not a directory: {mtx}. Run export_matrix() first.")
    r_script = (Path(script) if script is not None
                else Path(__file__).resolve().with_name("scdblfinder.R"))
    if not r_script.exists():
        raise TaskFailure(
            f"the scDblFinder R adapter is missing: {r_script}. It ships beside this module; "
            f"pass script= if it lives elsewhere.")

    out_path = Path(out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    log_path = Path(log)

    exported = read_labels_file(_resolve_optional_gz(mtx, "barcodes.tsv"), name="barcodes")
    seed_v = _require_int(seed, "seed")
    dbr_v = _require_number(dbr, "dbr", minimum=0.0, maximum=1.0, exclusive=True)
    threads_v = _require_int(threads, "threads", minimum=1)
    sd_label, sd_value = ("default", None) if dbr_sd is None else resolve_dbr_sd(dbr_sd, dbr_v)

    cmd = build_rscript_cmd(rscript_s, r_script, mtx, out_path, dbr_v, sd_value, seed_v,
                            threads=threads_v, refuse_xgboost_ge=refuse_xgboost_ge,
                            features_col=features_col)

    # Thread counts are pinned rather than left to the machine because xgboost sums in an order
    # that depends on how many threads it was given: unpinned, the scores depend on which node
    # the job landed on, and that is not something a recorded seed can repair.
    run_env = {"OMP_NUM_THREADS": str(threads_v),
               "OPENBLAS_NUM_THREADS": str(threads_v),
               "MKL_NUM_THREADS": str(threads_v)}
    if r_libs is not None:
        run_env["R_LIBS_USER"] = str(r_libs)
    if env:
        run_env.update({str(k): str(v) for k, v in env.items()})

    # A previous run's calls must never be accepted as this run's. Both the target and the
    # .partial the R adapter renames from are removed BEFORE Rscript is launched, so the file
    # found afterwards can only have been written by this invocation.
    stale = _clear_paths([out_path, Path(str(out_path) + ".partial")],
                         what=f"a previous scDblFinder run for {sample or out_path.name}")
    started = time.time()

    output = executor.shell(cmd, log=log_path, env=run_env, timeout_s=timeout_s)

    versions = parse_versions(output)
    r_metrics = parse_r_metrics(output)

    if not out_path.exists():
        raise TaskFailure(
            f"{r_script.name} exited 0 but {out_path} does not exist. The path was emptied "
            f"before the run, so this is not a stale-file question: the adapter writes a "
            f".partial file and renames it, and a missing output means it never reached the "
            f"write. Log: {log_path}")
    if out_path.stat().st_size == 0:
        raise TaskFailure(f"{out_path} is empty. Log: {log_path}")
    # Corroboration only. The removal above is the proof of authorship; an mtime is stamped by
    # whichever host ran Rscript and comparing it against this host's clock would fail runs on
    # skew alone.
    try:
        out_age_s = round(started - out_path.stat().st_mtime, 3)
    except OSError:
        out_age_s = None

    calls = read_calls(out_path, unscored=unscored, sample=sample)

    scored_set = set(calls.keys())
    exported_set = set(exported)
    if scored_set != exported_set:
        lost = sorted(exported_set - scored_set)
        extra = sorted(scored_set - exported_set)
        raise TaskFailure(
            f"the calls do not cover the export: {len(lost)} exported barcode(s) have no call "
            f"(e.g. {lost[:5]}), and {len(extra)} call(s) name a barcode that was not exported "
            f"(e.g. {extra[:5]}). A nucleus that quietly loses its call is recorded downstream "
            f"as not a doublet, which is not the same as never having been examined. "
            f"Log: {log_path}")

    where = str(log_path)
    reported_doublets = _require_int(_r_metric(r_metrics, "n_doublets", where), "n_doublets")
    if reported_doublets != calls.n_called:
        raise TaskFailure(
            f"the R adapter reported {reported_doublets} doublet(s) and the CSV holds "
            f"{calls.n_called}. One number reached by two routes must agree; it does not, so "
            f"one of the two records describes a different run. Log: {log_path}")
    reported_cells = _require_int(_r_metric(r_metrics, "n_cells", where), "n_cells")
    if reported_cells != len(calls):
        raise TaskFailure(
            f"the R adapter reported {reported_cells} cell(s) and the CSV holds {len(calls)}. "
            f"Log: {log_path}")
    reported_seed = _require_int(_r_metric(r_metrics, "seed_used", where), "seed_used")
    if reported_seed != seed_v:
        raise TaskFailure(
            f"the R adapter used seed {reported_seed} and {seed_v} was asked for. An argv that "
            f"arrived shifted produces exactly this, and every other argument is shifted too. "
            f"Log: {log_path}")
    reported_dbr = _require_number(_r_metric(r_metrics, "dbr_used", where), "dbr_used")
    if not math.isclose(reported_dbr, dbr_v, rel_tol=1e-9, abs_tol=1e-12):
        raise TaskFailure(
            f"the R adapter used dbr = {reported_dbr} and {dbr_v} was asked for. "
            f"Log: {log_path}")

    # dbr.sd, checked in both directions. This is the parameter the sweep varies and step 4
    # derives from the sweep, so a run that quietly used something else does not produce a wrong
    # number - it produces a sweep that measured one setting several times, and the resulting
    # flatness reads as a property of the data rather than of the prior.
    reported_sd = _r_metric(r_metrics, "dbr_sd_used", where)
    if sd_value is None:
        if reported_sd != "package-default":
            raise TaskFailure(
                f"dbr.sd was deliberately NOT passed - the token 'default' was sent so the "
                f"installed package would apply its own - but the R adapter reports having used "
                f"{reported_sd!r}. Whatever ran was not the run that was asked for, and the "
                f"sweep setting labelled 'default' would be recording a number instead. "
                f"Log: {log_path}")
    else:
        if reported_sd == "package-default":
            raise TaskFailure(
                f"dbr.sd = {sd_value} was passed and the R adapter reports having used the "
                f"package default instead. The argument was dropped somewhere between the argv "
                f"and scDblFinder(), so this setting of the sweep did not happen; recording it "
                f"would make two settings look identical because they were. Log: {log_path}")
        used_sd = _require_number(reported_sd, "dbr_sd_used")
        if not math.isclose(used_sd, sd_value, rel_tol=1e-9, abs_tol=1e-12):
            raise TaskFailure(
                f"the R adapter used dbr.sd = {used_sd} and {sd_value} was asked for. One number "
                f"reached by two routes must agree. dbr.sd is what this sweep varies, so a run "
                f"at the wrong one is not a small error: it is a point of the sweep that "
                f"describes a setting nobody chose. Log: {log_path}")

    metrics = {
        "sample": sample,
        "mtx_dir": str(mtx),
        "calls_csv": str(out_path),
        "log": str(log_path),
        "n_scored": len(calls),
        "n_called": calls.n_called,
        "rate_over_scored": calls.fraction_called("scored"),
        "n_unscored": calls.n_unscored,
        "unscored_supplied": calls.unscored is not None,
        "dbr": dbr_v,
        "dbr_sd_label": sd_label,
        "dbr_sd_value": sd_value,
        "dbr_sd_used_by_r": reported_sd,
        # What the installed scDblFinder DECLARES as its dbr.sd default, and whether that string
        # is the value it would use or only the placeholder in its signature. The R adapter
        # observes both and refuses to conflate them; see its header.
        "dbr_sd_formal_default": _r_metric(r_metrics, "dbr_sd_formal_default", where),
        "dbr_sd_formal_default_route": _r_metric(
            r_metrics, "dbr_sd_formal_default_route", where),
        "seed": seed_v,
        "threads": threads_v,
        "n_genes": _require_int(_r_metric(r_metrics, "n_genes", where), "n_genes"),
        "xgboost_ceiling": refuse_xgboost_ge,
        "stale_outputs_removed": stale,
        "output_age_s_at_start": out_age_s,
        "started_epoch_s": started,
        "freshness_proof": ("the calls CSV and its .partial were deleted before Rscript ran, so "
                            "the file read back was written by this invocation; the age is "
                            "corroboration from whichever host ran R, not the proof"),
        "cross_check": ("n_cells, n_doublets, seed, dbr and dbr.sd agree between the CSV and the "
                        "R adapter's own report"),
    }
    if calls.unscored is not None:
        metrics["rate_over_all"] = calls.fraction_called("all")
    return {"outputs": [out_path], "metrics": metrics, "versions": versions}


def sweep(dbr_sd_values, rscript, mtx_dirs, out_dir, dbr, seed, log_dir, executor, *,
          script=None, umi_by_barcode=None, unscored=None, threads: int = 1, r_libs=None,
          refuse_xgboost_ge: str | None = XGBOOST_SHIM_FROM, features_col: int = 1,
          timeout_s: int | None = None, env: Mapping | None = None) -> dict:
    """Run the dbr.sd sweep across every sample and return the numbers step 4 judges.

    The sweep exists because a flat rate is suspect rather than reassuring: at the package
    default scDblFinder returned 10.04-10.79% across libraries differing 2.5-fold in size, and
    that flatness belonged to the prior rather than to the data. Only a sweep separates the two,
    which is why `dbr.sd` is DERIVED and not a command-line argument.

    Nothing here judges the result - `modules/04_doublets` does - and the shapes are returned to
    fit it without translation::

        s = metrics["per_setting"][label]
        SweepResult(setting=label,
                    per_sample_rate=s["per_sample_rate_over_scored"],
                    deep_decile_rate=s["deep_decile"]["rate_over_scored"])

    `per_sample_rate_over_scored` names its denominator because the same calls give 7.52% of
    nuclei scored and 5.70% of all cells, and a rate compared against a published band without
    saying which it is, is not a comparison. `deep_decile` is None unless `umi_by_barcode` was
    supplied, rather than 0.0: not computed is not the same as computed and small, and
    `recommend()` reads a missing decile as no evidence rather than as reassurance.

    `mtx_dirs`, `unscored` and `umi_by_barcode` are keyed by sample, and the last two must cover
    every sample in the first if they are given at all. A partial map is refused rather than
    filled in, because a cohort rate computed from a subset of the cohort still prints.

    Every sample is scored with the same seed at every setting, so a difference between two
    settings is the setting.
    """
    if not isinstance(mtx_dirs, _abc.Mapping):
        raise TaskFailure("mtx_dirs must be a mapping of {sample: exported mtx directory}.")
    if not mtx_dirs:
        raise TaskFailure("mtx_dirs is empty; there is nothing to sweep.")
    settings = list(dbr_sd_values)
    if not settings:
        raise TaskFailure(
            "no dbr.sd values were given. A single setting is not a sweep, and the point of the "
            "sweep is to show whether the rate is a measurement or the prior.")

    dbr_v = _require_number(dbr, "dbr", minimum=0.0, maximum=1.0, exclusive=True)
    seed_v = _require_int(seed, "seed")
    out_root = Path(out_dir)
    log_root = Path(log_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    for label, mapping in (("umi_by_barcode", umi_by_barcode), ("unscored", unscored)):
        if mapping is None:
            continue
        if not isinstance(mapping, _abc.Mapping):
            raise TaskFailure(f"{label} must be a mapping keyed by sample.")
        absent = [s for s in mtx_dirs if s not in mapping]
        if absent:
            raise TaskFailure(
                f"{label} is missing {len(absent)} sample(s): {sorted(map(str, absent))[:5]}. A "
                f"partial map would silently change what the cohort numbers are computed over, "
                f"and the result would still print.")

    resolved = []
    seen: dict = {}
    for token in settings:
        label, value = resolve_dbr_sd(token, dbr_v)
        if label in seen:
            raise TaskFailure(
                f"dbr.sd setting {label!r} appears twice in the sweep. Two runs under one label "
                f"overwrite each other's calls, and the surviving one is whichever ran last.")
        seen[label] = value
        resolved.append((label, value))

    outputs: list = []
    versions: dict = {}
    per_setting: dict = {}

    for label, value in resolved:
        slug = _safe_label(label)
        rates: dict = {}
        n_scored: dict = {}
        n_called: dict = {}
        csvs: dict = {}
        pooled_umi: dict = {}
        pooled_scored: list = []
        pooled_called: list = []

        for name in sorted(mtx_dirs):
            sample_unscored = None if unscored is None else unscored[name]
            out_csv = out_root / f"{_safe_label(name)}.dbrsd_{slug}.calls.csv"
            log_path = log_root / f"{_safe_label(name)}.dbrsd_{slug}.log"
            result = run_scdblfinder(
                rscript, mtx_dirs[name], out_csv, dbr_v, value, seed_v, log_path, executor,
                script=script, unscored=sample_unscored, sample=name, threads=threads,
                r_libs=r_libs, refuse_xgboost_ge=refuse_xgboost_ge, features_col=features_col,
                timeout_s=timeout_s, env=env)
            outputs.extend(result["outputs"])
            csvs[name] = str(out_csv)
            rates[name] = result["metrics"]["rate_over_scored"]
            n_scored[name] = result["metrics"]["n_scored"]
            n_called[name] = result["metrics"]["n_called"]
            for tool, version in result["versions"].items():
                if tool in versions and versions[tool] != version:
                    raise TaskFailure(
                        f"{tool} reported {versions[tool]!r} and then {version!r} within one "
                        f"sweep. A cohort processed by two versions of the same tool is not one "
                        f"cohort, and the sweep would be comparing the versions rather than the "
                        f"settings.")
                versions[tool] = version
            if umi_by_barcode is not None:
                calls = read_calls(out_csv, unscored=sample_unscored, sample=name)
                pooled_scored.extend((name, bc) for bc in calls)
                pooled_called.extend((name, bc) for bc in calls.called())
                for barcode, depth in umi_by_barcode[name].items():
                    pooled_umi[(name, barcode)] = depth

        entry = {
            "dbr_sd_value": value,
            "per_sample_rate_over_scored": rates,
            "per_sample_n_scored": n_scored,
            "per_sample_n_called": n_called,
            "calls_csv": csvs,
            "deep_decile": None,
            "deep_decile_supplied": umi_by_barcode is not None,
        }
        if umi_by_barcode is not None:
            entry["deep_decile"] = deep_decile_rate(pooled_umi, pooled_called, pooled_scored)
        per_setting[label] = entry

    metrics = {
        "dbr": dbr_v,
        "seed": seed_v,
        "samples": sorted(mtx_dirs),
        "settings": [label for label, _ in resolved],
        "per_setting": per_setting,
        "note": ("per_sample_rate_over_scored is the fraction of the nuclei actually scored, "
                 "which is what SweepResult.per_sample_rate is documented to hold. Nuclei below "
                 "the light floor are UNKNOWN and appear in no denominator here."),
    }
    return {"outputs": outputs, "metrics": metrics, "versions": versions}


# --------------------------------------------------------------------------- detector contract


class ScDblFinderDetector:
    """scDblFinder, as `modules/04_doublets.Detector` requires a detector to declare itself.

    Every field below is read out of the tool's own documentation rather than inferred from its
    behaviour, because that is what `check_detector()` is checking - whether the properties that
    decide whether a configuration is safe were ever stated:

      reproducible                a seed fixes the result; contrast DoubletDetection, whose own
                                  documentation says PhenoGraph does not support random seeds
      needs_empty_drops_removed   yes; the null is built by summing observed transcriptomes, so
                                  a pool holding empty droplets calibrates it on debris
      min_umi_floor               200, "to avoid errors" - the only one of the four common tools
                                  to put a number on it
      imports_rate_prior          yes, and its default is the 10x loading formula, which on a
                                  Singleron cohort tracked library size at r = 0.872

    `score()` applies no floor of its own. The protocol forbids subsetting the input, so a
    scoring set that still contains empty droplets is refused rather than quietly trimmed:
    choosing it belongs to step 3, in writing.
    """

    name = "scDblFinder"
    reproducible = True
    needs_empty_drops_removed = True
    min_umi_floor = SCDBLFINDER_MIN_UMI
    imports_rate_prior = True

    def __init__(self, rscript, work_dir, executor, *, script=None, r_libs=None,
                 refuse_xgboost_ge: str | None = XGBOOST_SHIM_FROM, threads: int = 1,
                 timeout_s: int | None = None):
        self.rscript = _require_label(str(rscript), "rscript")
        self.work_dir = Path(work_dir)
        self.executor = executor
        self.script = script
        self.r_libs = r_libs
        self.refuse_xgboost_ge = refuse_xgboost_ge
        self.threads = _require_int(threads, "threads", minimum=1)
        self.timeout_s = timeout_s
        #: The full adapter result of the most recent `score()`, or None if it has not run.
        #: Kept because `score()` returns the protocol's `(scores, calls)` and the outputs,
        #: metrics and observed versions must not be discarded to fit that shape.
        self.last_result = None

    def score(self, matrix, sample: str, seed: int, *, dbr=None, dbr_sd=None, unscored=None,
              layer=None, gzip_mtx: bool = False):
        """Return `(scores, calls)` for exactly the nuclei given, and score nothing else.

        `matrix` is either a directory already holding a MatrixMarket triple or a path to an
        .h5ad, in which case every barcode in it is exported and no floor is applied here. A
        detector that trims its own input turns a threshold nobody recorded into a coverage gap
        nobody counted, so an input still holding zero-count droplets is refused by
        `export_matrix`, with the step that owns the decision named in the message.
        """
        sample_name = _require_label(sample, "sample")
        src = Path(matrix)
        work = self.work_dir / _safe_label(sample_name)
        if src.is_dir():
            mtx_dir = src
            export_versions: dict = {}
        elif src.is_file():
            exported = export_matrix(src, work / "mtx", None, layer=layer, gzip_mtx=gzip_mtx)
            mtx_dir = exported.mtx_dir
            export_versions = exported.versions
        else:
            raise TaskFailure(
                f"matrix is neither an exported MatrixMarket directory nor a file: {src}")
        slug = _safe_label(sample_name)
        out_csv = work / f"{slug}.calls.csv"
        result = run_scdblfinder(
            self.rscript, mtx_dir, out_csv, dbr, dbr_sd, seed, work / f"{slug}.log",
            self.executor, script=self.script, unscored=unscored, sample=sample_name,
            threads=self.threads, r_libs=self.r_libs,
            refuse_xgboost_ge=self.refuse_xgboost_ge, timeout_s=self.timeout_s)
        calls = read_calls(out_csv, unscored=unscored, sample=sample_name)
        self.last_result = {
            "outputs": result["outputs"],
            "metrics": result["metrics"],
            "versions": {**export_versions, **result["versions"]},
        }
        scores = {barcode: value[0] for barcode, value in calls.items()}
        return scores, calls
