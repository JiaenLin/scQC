# Execution layer: this module computes numbers and writes files. It decides nothing.
# Nothing here removes an observation. The valley finder reports where a density minimum sits
# and whether there are two modes to put one between; modules/05_quality/quality.py decides
# whether that is a threshold. The cluster profiler describes clusters;
# modules/06_cluster_check/cluster_flags.py decides which are flagged. The removal itself
# happens only in modules/07_apply/apply.py, under a recorded approval.
"""scanpy operations: per-cell QC metrics, the KDE valley, clustering, and the cluster profile.

WHAT THIS FILE IS, AND WHAT IT IS NOT ALLOWED TO BE

The decision layer consumes numbers it cannot produce. `quality.derive()` takes `Valley(sample,
metric, value, bimodal)` objects; `cluster_flags.apply_flags()` takes rows carrying
`umi_frac_of_sample`, `median_pct_mt`, `pct_uninformative` and `pct_doublet`. This adapter
produces exactly those and nothing more opinionated: no policy, no samplesheet, no default that
stands in for a declaration. Every function takes explicit paths and parameters and the caller
decides what they mean.

THE VALLEY, WHICH IS THE PART TO READ CAREFULLY

`quality.py` records the measurement that motivates the whole procedure: across ten libraries of
one cohort the density valley ranged 274-473 UMI and 184-352 genes. The PROCEDURE transfers; the
number does not. `find_valley()` is that procedure.

A density minimum exists in any smooth curve. Returning the minimum of a unimodal density means
returning the flank of the only mode, which dresses an arbitrary cut as a measurement - and it
does so in a form nobody downstream can question, because a number on a page carries no record of
how it was obtained. So bimodality is TESTED here, the test is named in the docstring of
`find_valley()` together with what it cannot detect, and `bimodal=False` is returned rather than a
value with a caveat attached. `quality.derive()` then refuses the cohort rather than proposing a
floor, which is the correct outcome: the operator chooses the cut explicitly and records it as a
judgement, the way the mitochondrial ceiling already is.

NO GENE CLASS IS EXCLUDED FROM VARIABLE-GENE SELECTION

`cluster()` computes highly-variable genes over every gene in the object and passes the result to
PCA untouched. There is no exclusion list and no `& ~excluded` term, and one must not be added.
The reason is in docs/PRINCIPLES.md section 1: a regex written to match ribosomal genes matched a
ribosomal protein KINASE - an mTOR signalling gene - in a metabolic study, and it was called
"ribosomal genes" until somebody printed the list. What this module does instead is REPORT: the
symbols of every flagged gene that selection chose are recorded in
`adata.uns["scqc_cluster"]["flagged_hvg"]`, so the influence of a gene class can be read without
anything having been decided about it. Reporting influence is permitted; preventing it is not.

Where a gene class is used at all - the mitochondrial percentage, the locked uninformative set
behind criterion C - the actual matched symbols are recorded beside the count, for the same
reason. A gene list nobody printed is a hypothesis about the reference.

ORDER MATTERS, AND GETTING IT WRONG IS SILENT

`qc_metrics()` must run BEFORE `cluster()`. `cluster()` normalises `adata.X` in place, after which
every cell's counts sum to the same number; `total_counts` computed at that point is a constant
and `median_umi` becomes meaningless while still printing as a plausible figure. `cluster()`
therefore refuses to start unless `obs["total_counts"]` is already present.

UNITS, STATED BECAUSE AN IMPLICIT SCALE IS HOW TWO CORRECT FUNCTIONS PRODUCE A WRONG NUMBER

  find_valley           returns the valley in the metric's OWN units (UMI, or detected genes),
                        not in log10, because that is the scale `quality.UMI_BOUNDS` is in.
  median_umi            raw UMI.
  umi_frac_of_sample    a FRACTION (cluster median / that sample's median). Criterion A is
                        `< 0.5`, so a percentage here would silently never fire.
  median_pct_mt         PERCENT, 0-100, as scanpy's `pct_counts_mt` already is.
  pct_uninformative     PERCENT of the examined top markers, 0-100.
  pct_mt_markers        PERCENT, 0-100. Same denominator.
  pct_ribo_markers      PERCENT, 0-100. Same denominator.
  pct_doublet           PERCENT of the SCORED cells in the cluster, 0-100.

UNKNOWN IS NEVER A VALUE, AND UNKNOWN IS NOT ONLY None AND float NaN

Every quantity that could not be computed is emitted as None: never 0, never NaN, never False.
`cluster_flags` is three-valued and depends on receiving the gap rather than a floor
(docs/PRINCIPLES.md section 4). The two places this bites hardest are handled explicitly:
`pct_doublet` is computed over the cells that were SCORED, so a never-scored cell is absent from
both the numerator and the denominator rather than counted as a singlet; and a cluster whose
denominator is zero gets None rather than 0%.

Recognising an unknown is one predicate, `_unknown()`, used everywhere in this module, and it is
wider than `v is None or v != v` on purpose. A table delivers pandas.NA, pandas.NaT,
numpy.ma.masked and numpy scalars that are not float subclasses; every one of those is
`is not None`, most are not `float`, and each then meets `>=` (False), `int()` or `bool()`
somewhere downstream, where a False reads as "measured and failed the test" - which is a PASS at
every gate here. `_is_true()` is its companion for flags: `x is True` is False for
numpy.bool_(True), so an identity test drops every genuinely flagged row coming from a
numpy-backed table. Both import numpy and pandas lazily; this file stays importable with neither.

RUNNING IT SOMEWHERE ELSE, AND NOT ACCEPTING THE LAST RUN'S ANSWER

The analysis stack is usually not installed where the orchestrator runs. `run_scanpy_op()` sends
one operation through an `Executor` - locally or to a scheduler - by invoking this same file as a
script under the interpreter of the analysis environment, then reads back a metrics JSON and
checks every file it claims. Argument construction (`build_scanpy_cmd`) and output parsing
(`parse_metrics_json`) are separate pure functions so both can be tested without the tool.

That a claimed output EXISTS afterwards is not evidence this run produced it - every run of a
resolution sweep writes to the same prefix, and a tool that exits 0 having written nothing leaves
the previous run's files exactly where the check looks. So `run_scanpy_op()` deletes every output
the previous run declared BEFORE launching, passes a token this invocation alone knows, and
afterwards requires the metrics file to carry that token and every output to match the size and
mtime the writer recorded as it wrote it. Read its docstring before changing any of the three.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

# This file is also executed as a script by the analysis environment's interpreter, in which case
# the repository root is not necessarily on sys.path and `engine.task` would not import. Inserting
# it here rather than catching ImportError is deliberate: a fallback definition of TaskFailure
# would be a DIFFERENT class from the one the orchestrator catches, so `except TaskFailure` would
# stop working in exactly the situation it exists for.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from adapters.cellbender import ID_COLUMNS  # noqa: E402  (see the sys.path note above)
from engine.fs import VISIBILITY_TIMEOUT_S, await_visible  # noqa: E402  (see the note above)
from engine.task import TaskFailure  # noqa: E402  (see the sys.path note above)

#: Must equal `TOPN` in modules/06_cluster_check/cluster_flags.py. That module's criterion C is
#: "share of the top-N markers in the locked mt+ribo set"; if the two numbers drift apart the
#: proposed C threshold is derived against one denominator and applied against another.
MARKER_TOPN = 20

#: The keys `cluster_flags.apply_flags()` reads. Every row this module emits carries all of them,
#: with None where the value could not be computed. Read that module before changing this tuple.
PROFILE_KEYS = (
    "sample",
    "cluster",
    "n",
    "median_umi",
    "umi_frac_of_sample",
    "median_pct_mt",
    "pct_uninformative",
    "pct_mt_markers",
    "pct_ribo_markers",
    "pct_doublet",
)

#: Which per-cell column each valley metric is measured on. There is no default metric and no
#: fallback column: a valley found on the wrong column is a number with no way to be wrong out
#: loud. `quality.derive()` treats "umi" specially (UMI_BOUNDS) and everything else as genes.
VALLEY_METRIC_COLUMNS = {
    "umi": "total_counts",
    "genes": "n_genes_by_counts",
}

#: DECLARED. The mitochondrial percentage above which a droplet is excluded from the population
#: the mitochondrial CEILING is derived over. It is not a filter and removes nothing: it says
#: which droplets are eligible to describe where the healthy population ends. A droplet more than
#: half mitochondrial by count is ambient RNA rather than a cell with high mitochondrial content,
#: and letting it set the MAD widens the fence that decides which real cells survive. Override
#: per run with the `mito_derivation_max` parameter. The reasoning, and the measurement behind
#: the value, are in `_op_valley` beside the code that applies it.
MITO_DERIVATION_MAX = 50.0

#: A library that is present but was not used has a version; one that is not installed does not.
#: Kept distinct from `engine.provenance.NOT_INVOKED` ("not invoked") because they are different
#: facts and a reader has to be able to tell them apart.
ABSENT = "not installed"

#: Operations `run_scanpy_op()` and the script entry point understand.
OPS = ("qc", "valley", "cluster", "apply_measure", "apply_write")

#: Ops that do not read an input object. `apply_write` takes a LIST of objects in its params and
#: `main()` would otherwise load one of them for nothing - a hundred megabytes to be discarded.
OPS_WITHOUT_INPUT = ("apply_write",)

# Values below this maximum, in a matrix that is not integral, are the signature of data that has
# already been log1p'd. Running normalize_total on it produces an embedding that looks entirely
# ordinary and describes data transformed twice.
_LOG1P_MAX = 30.0

# Relative spread of per-cell totals below which the matrix has already been normalised to a
# constant sum. Real libraries do not produce this.
_CONSTANT_TOTAL_CV = 1e-6


# --------------------------------------------------------------------------------------------
# small shared helpers
# --------------------------------------------------------------------------------------------

def _try_import(name):
    """Import a third-party module if it is there, else None.

    This module is imported by the stdlib-only decision layer as well as run inside the analysis
    environment, so no unknown-value predicate may make numpy or pandas a hard requirement.
    """
    try:
        return __import__(name)
    except Exception:                                        # noqa: BLE001 - absence is the fact
        return None


def _unknown(v) -> bool:
    """THE unknown predicate for this module. True when a value carries no information.

    Same three-valued contract as `cluster_flags._unknown`, and deliberately wider than
    `v is None or v != v`, because a table does not deliver only those two. Every sentinel below
    is `is not None`, and most of them are not `float`, so a predicate that enumerates None and
    float NaN lets them through - and the value then meets `>=` (False), `int()` (TypeError or a
    silent 0) or `bool()` (pd.NA raises), where a False reads as "measured and failed the test",
    which is a PASS in every gate this pipeline has.

    Covered, all of them tested:

      None                     the caller's own "not computed"
      float('nan')             arithmetic that could not be done
      numpy.float64('nan')     a float SUBCLASS, so it is caught by the isinstance branch
      numpy.float32('nan')     NOT a float subclass - caught by the numpy.generic branch
      numpy.datetime64('NaT')  likewise
      pandas.NA                pandas' own missing scalar; `bool(pd.NA)` RAISES
      pandas.NaT               missing timestamp
      numpy.ma.masked          the MaskedConstant singleton; an ndarray subclass, not a float

    NOT unknown: 0, 0.0, False, '' and numpy.bool_(False). An absent measurement and a measured
    zero are different findings and this predicate must not merge them.

    numpy and pandas are imported lazily and their absence is not an error: with neither
    installed this degrades to None-and-float-NaN, which is all that can arrive in that case.
    An array is never "unknown" - only scalars are - so a vectorised result is answered False
    and the caller is expected to have gone elementwise.
    """
    if v is None:
        return True
    if isinstance(v, float):                # includes numpy.float64, which subclasses float
        # bool(), because numpy.float64('nan') != itself yields numpy.bool_(True), and this
        # predicate returning a numpy boolean would recreate the identity trap it exists to
        # close: `_unknown(x) is True` would be False for a value it had just called unknown.
        return bool(v != v)
    if isinstance(v, (bool, int, str, bytes, bytearray)):
        return False                        # a real value; none of these can be a sentinel

    np = _try_import("numpy")
    if np is not None:
        masked = getattr(getattr(np, "ma", None), "masked", None)
        if masked is not None and v is masked:
            return True
        if isinstance(v, np.generic):       # np.float32('nan'), np.datetime64('NaT'), np.bool_
            try:
                return bool(np.isnan(v))
            except (TypeError, ValueError): # np.str_, np.bool_: real values, isnan undefined
                return False

    pd = _try_import("pandas")
    if pd is not None:
        try:
            r = pd.isna(v)
        except (TypeError, ValueError):     # an object pandas cannot judge is not thereby missing
            return False
        if r is True or r is False:
            return bool(r)
        if np is not None and isinstance(r, np.bool_):
            return bool(r)
        return False                        # array-valued: v was not a scalar
    return False


def _is_true(v) -> bool:
    """Truthiness of a possibly-numpy, possibly-missing flag. Unknown is False, and never raises.

    `x is True` is False for numpy.bool_(True), so an identity test silently drops every genuinely
    flagged row that arrives from a numpy-backed table. `bool(x)` is right - but `bool(pd.NA)`
    raises and `bool(nan)` is True, so the unknown check has to come first and cannot be skipped.
    """
    if _unknown(v):
        return False
    try:
        return bool(v)
    except (TypeError, ValueError):         # a value that refuses to be judged is not a True
        return False


def _clean(v):
    """Map any unknown to None so it is emitted as an empty cell rather than 'nan' or '<NA>'.

    A blank cell is read back by pandas as NaN and caught by `cluster_flags._unknown`; the literal
    text 'nan' or '<NA>' in a CSV is read back as a string and is not.
    """
    return None if _unknown(v) else v


def _float_array(values, *, unknown_as_nan: bool = True):
    """`values` as a float64 array, with every unknown sentinel mapped to NaN.

    `np.asarray(col, dtype=float)` raises on an object column holding pd.NA and silently coerces
    numpy.ma.masked to its fill value, so the conversion itself is where a sentinel becomes a
    number. Object columns are converted elementwise through `_unknown`; numeric columns are
    already float and take the fast path.
    """
    import numpy as np

    arr = np.asarray(values)
    if arr.dtype.kind in "fiub":
        return arr.astype(float, copy=False)
    out = np.empty(arr.size, dtype=float)
    flat = arr.ravel()
    for i in range(flat.size):
        v = flat[i]
        if _unknown(v):
            out[i] = float("nan") if unknown_as_nan else 0.0
            continue
        try:
            out[i] = float(v)
        except (TypeError, ValueError):
            out[i] = float("nan")
    return out.reshape(arr.shape)


def _bool_array(values) -> tuple:
    """`(mask, n_unknown)`: elementwise truthiness with every unknown mapped to False.

    Returns the count of unknowns as well as the mask because "flagged False" and "never
    established" are different findings and the caller has to be able to record which it had.
    Mapping unknown to False is the conservative direction everywhere this is used - an unknown
    highly-variable flag does not put a gene into the embedding, an unknown mt flag does not put a
    gene into a reported class - and it is reversible: nothing is deleted from the object.
    """
    import numpy as np

    arr = np.asarray(values)
    if arr.dtype.kind == "b":
        return arr.astype(bool, copy=False), 0
    flat = arr.ravel()
    out = np.zeros(flat.size, dtype=bool)
    n_unknown = 0
    for i in range(flat.size):
        v = flat[i]
        if _unknown(v):
            n_unknown += 1
            continue
        out[i] = _is_true(v)
    return out.reshape(arr.shape), n_unknown


#: The key `run_scanpy_op()` puts a per-invocation token under in the params file, and that
#: `main()` copies into the metrics file. A metrics file left over from an earlier run carries an
#: earlier token, which is how "this file exists" is turned into "this invocation wrote it".
RUN_TOKEN_KEY = "_scqc_run_token"


def _require_str(params: dict, key: str) -> str:
    """Fetch a DECLARED parameter, refusing rather than substituting anything for its absence."""
    if key not in params or _unknown(params[key]) or str(params[key]).strip() == "":
        raise TaskFailure(
            f"required parameter {key!r} was not supplied. There is no default for it - it is a "
            f"property of the species, platform or design, and a pipeline that guesses it "
            f"produces a result nobody can attribute. Supplied keys: "
            f"{', '.join(sorted(params)) or '(none)'}")
    return str(params[key])


def observed_versions() -> dict:
    """Versions of the analysis stack, as reported by the libraries themselves in this process.

    Nothing is read from a lockfile or an environment name. A library that is not installed is
    recorded as `not installed` rather than omitted, so a reader can tell a stack that lacked
    leidenalg from one where nobody looked.
    """
    out = {}
    for name in ("scanpy", "anndata", "numpy", "scipy", "pandas", "leidenalg", "igraph",
                 "sklearn"):
        try:
            mod = __import__(name)
        except Exception:                                    # noqa: BLE001 - absence is the fact
            out[name] = ABSENT
            continue
        v = getattr(mod, "__version__", None)
        out[name] = str(v) if v else "installed, no __version__"
    out["python"] = sys.version.split()[0]
    return out


def _matrix_value_sample(matrix, max_values: int = 100_000):
    """A strided sample of the stored values of a dense or sparse matrix.

    Strided rather than head-of-array: the first N entries of a CSR matrix are the first few
    cells, and a check that only ever looks at the first few cells is a check that a per-cell
    transform can walk straight past.
    """
    import numpy as np

    data = getattr(matrix, "data", None)
    arr = np.asarray(data if data is not None else matrix).ravel()
    if arr.size == 0:
        return arr
    step = max(1, int(arr.size // max(1, max_values)))
    return arr[::step][:max_values]


def _refuse_if_transformed(matrix, where: str, allow_transformed: bool = False) -> list:
    """Refuse a matrix that has already been normalised or logged; report non-integrality.

    Two signatures are decisive and are refused, because in both cases the arithmetic downstream
    is applied twice and the result reads as ordinary:

      negative values          scaled or z-scored data - counts have no negative values
      non-integral and small   the range of log1p'd data

    Non-integrality ALONE is reported, not refused. CellBender's denoised output is this
    pipeline's own step-1 product and is not guaranteed integral; a gate that refuses it fires on
    correct behaviour, and a gate that fires on correct behaviour gets switched off.
    """
    import numpy as np

    notes = []
    s = _matrix_value_sample(matrix)
    if s.size == 0:
        return ["the matrix has no stored values - every entry is zero"]
    s = s.astype(float, copy=False)
    if not np.all(np.isfinite(s)):
        raise TaskFailure(f"{where}: the matrix contains non-finite values (NaN or inf). "
                          f"Every downstream sum, median and percentage inherits them silently.")
    smin, smax = float(s.min()), float(s.max())
    integral = bool(np.all(np.equal(np.mod(s, 1), 0)))
    if smin < 0 and not allow_transformed:
        raise TaskFailure(
            f"{where}: the matrix contains negative values (minimum {smin:.4g}), so it has been "
            f"scaled or z-scored. Counts are needed here. Pass the counts layer, or "
            f"allow_transformed=True if this is deliberate and you have recorded why.")
    if not integral and smax < _LOG1P_MAX and not allow_transformed:
        raise TaskFailure(
            f"{where}: the matrix is non-integral with a maximum of {smax:.4g}, which is the "
            f"range of log1p'd data rather than of counts. Normalising it again would produce a "
            f"perfectly ordinary-looking embedding of data transformed twice. Pass the counts "
            f"layer, or allow_transformed=True if this is deliberate.")
    if not integral:
        notes.append(f"{where}: values are not integral (sampled range {smin:.4g}-{smax:.4g}). "
                     f"Expected for a denoised matrix; recorded rather than refused")
    return notes


def _refuse_if_totals_constant(totals, where: str) -> None:
    """Refuse when per-cell totals are already equalised - the fingerprint of normalize_total."""
    import numpy as np

    t = _float_array(totals).ravel()      # every unknown sentinel becomes NaN, then drops out
    t = t[np.isfinite(t)]
    if t.size == 0:
        raise TaskFailure(f"{where}: no finite per-cell totals; nothing can be measured from this")
    m = float(t.mean())
    if m > 0 and float(t.std()) / m < _CONSTANT_TOTAL_CV:
        raise TaskFailure(
            f"{where}: every cell has the same total ({m:.6g}), which is what normalize_total "
            f"leaves behind. Counts computed from this are a constant, and a median UMI derived "
            f"from a constant still prints as a number. Supply the counts.")


# --------------------------------------------------------------------------------------------
# step 5 - the valley
# --------------------------------------------------------------------------------------------

def find_modes(density) -> tuple:
    """Indices of the interior local maxima and minima of a sampled curve.

    Pure, and deliberately free of numpy so the mode logic can be tested on a hand-written list.

    Plateaus are handled by comparing RUNS of equal value rather than neighbouring points: a
    density that is flat across three grid points at its peak has no strictly-greater-than-both
    -neighbours sample, and a naive scan reports zero modes on a curve that plainly has one. The
    index reported for a run is its midpoint.

    Endpoints are excluded. A maximum at the first or last grid point is not a mode inside the
    observed range - it is the curve still rising at the edge - and treating it as one invents a
    mode wherever the data was truncated. `find_valley()` reports that condition separately, as
    `mode_at_left_edge` / `mode_at_right_edge`, because a lower mode pinned against the left edge
    is the signature of an input that was cell-called before it arrived.
    """
    d = list(density)
    n = len(d)
    if n < 3:
        return [], []
    runs = []                                     # (start, end_inclusive, value)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and d[j + 1] == d[i]:
            j += 1
        runs.append((i, j, d[i]))
        i = j + 1
    peaks, troughs = [], []
    for k in range(1, len(runs) - 1):
        start, end, val = runs[k]
        prev, nxt = runs[k - 1][2], runs[k + 1][2]
        mid = (start + end) // 2
        if val > prev and val > nxt:
            peaks.append(mid)
        elif val < prev and val < nxt:
            troughs.append(mid)
    return peaks, troughs


def valley_between(grid, density, peaks) -> dict:
    """The minimum between the two TALLEST modes, with the two quantities that judge it.

    Pure, stdlib-only, and separated from the density estimation so it can be tested against a
    curve whose answer is known by construction.

    Returns a dict always - never a bare number - because the caller must be able to see WHY
    there is no valley as easily as where one is:

      mode_lo_index, mode_hi_index   the two tallest modes, ordered by POSITION not by height
      valley_index                   argmin of the density strictly between them, or None
      depth                          1 - d(valley) / min(d(mode_lo), d(mode_hi)); None when that
                                     minimum is zero, because the ratio is then undefined rather
                                     than large. Guarding the denominator with an epsilon would
                                     report a depth of ~1.0 for a curve with no lower mode at all
      separation_log10               distance between the two modes on the grid's own scale
      reason                         why there is no valley, when there is none
    """
    g, d, pk = list(grid), list(density), list(peaks)
    out = {"mode_lo_index": None, "mode_hi_index": None, "valley_index": None,
           "depth": None, "separation_log10": None, "reason": ""}
    if len(pk) < 2:
        out["reason"] = (f"{len(pk)} interior mode(s) found at this bandwidth - a valley needs "
                         f"two modes to sit between")
        return out
    tallest = sorted(pk, key=lambda i: d[i], reverse=True)[:2]
    lo, hi = sorted(tallest)
    out["mode_lo_index"], out["mode_hi_index"] = lo, hi
    out["separation_log10"] = g[hi] - g[lo]
    inner = range(lo + 1, hi)
    if not inner:
        out["reason"] = "the two tallest modes are adjacent grid points with nothing between them"
        return out
    v = min(inner, key=lambda i: d[i])
    out["valley_index"] = v
    floor_ = min(d[lo], d[hi])
    if floor_ > 0:
        out["depth"] = 1.0 - (d[v] / floor_)
    else:
        out["reason"] = ("the shorter of the two modes has zero density, so relative valley "
                         "depth is UNDEFINED rather than complete")
    return out


def mode_masses(values_log10, cut) -> dict:
    """How many OBSERVATIONS sit on each side of the valley, and what share of the total.

    Pure and stdlib-only, so it can be checked against a hand-counted list.

    This is the empirical mass, counted from the data, NOT the area under the kernel density
    estimate. That distinction is the whole point. A KDE places a kernel of fixed width on every
    point, so two or three barcodes sitting together far out in a tail produce a local maximum
    whose HEIGHT is a property of the bandwidth and whose supporting evidence is three
    observations. Density height cannot tell those apart from a population; a count can.

    Returns `{"n_total", "n_below", "n_above", "mass_below", "mass_above"}` with the shares as
    fractions of `n_total`, or all-None when there is nothing to count.
    """
    out = {"n_total": 0, "n_below": 0, "n_above": 0, "mass_below": None, "mass_above": None}
    if _unknown(cut):
        return out
    c = float(cut)
    total = below = 0
    for v in values_log10:
        if _unknown(v):
            continue
        total += 1
        if float(v) <= c:
            below += 1
    out["n_total"] = total
    out["n_below"] = below
    out["n_above"] = total - below
    if total > 0:
        out["mass_below"] = below / total
        out["mass_above"] = (total - below) / total
    return out


def mode_spread(values_log10, cut, side: str = "smaller") -> float:
    """Robust log10 dispersion of the observations on one side of the valley.

    THE DISCRIMINATOR MASS COULD NOT PROVIDE. Counting alone cannot separate a small real
    population from a kernel bump: 30 genuine cells in a 3,000-barcode object are 1.0% of it,
    while 25 co-located barcodes in a 20,000-barcode object are 0.125% - and any threshold placed
    between those two numbers is a number chosen to fit two examples.

    Dispersion separates them by their nature rather than their size. A population of cells is
    drawn from a distribution and is SPREAD OUT: a lognormal with sigma 0.5 has a log10 spread
    near 0.22 whatever its size. A kernel bump is the kernel - the barcodes under it sit at
    essentially one value, so their spread is near zero, and it stays near zero however many of
    them there are.

    Returns the interquartile range in log10 units, which ignores the single most extreme point
    rather than being defined by it. Returns None when the side holds too few observations to
    have a dispersion at all.
    """
    vals = sorted(float(v) for v in values_log10 if not _unknown(v))
    if not vals or _unknown(cut):
        return None
    c = float(cut)
    below = [v for v in vals if v <= c]
    above = [v for v in vals if v > c]
    group = (below if len(below) <= len(above) else above) if side == "smaller" else (
        above if side == "above" else below)
    if len(group) < 4:
        # Fewer than four points have no interquartile range worth the name. That is not a
        # spread of zero - it is an absence of evidence, and it is returned as one.
        return None
    g = sorted(group)
    q1 = g[int(0.25 * (len(g) - 1))]
    q3 = g[int(0.75 * (len(g) - 1))]
    return float(q3 - q1)


def modes_persist(valley_log10, mode_lo_log10, mode_hi_log10, wide_mode_positions_log10) -> dict:
    """Whether the SEPARATION survives a wider kernel: a mode still on each side of the valley.

    Pure and stdlib-only.

    A different question from the one this file used to ask, and from the first thing one reaches
    for instead of it. Both of those were tried and both were wrong:

      counting modes at the wider bandwidth (what the code did, and what `find_valley`'s docstring
      described as the criterion doing most of the work) passes whenever two bumps of ANY kind
      remain. On a unimodal log-normal with two extreme barcodes appended, the main mode and the
      kernel bump over those two barcodes are two modes, so it passed - on the exact input the
      test exists to reject;

      matching each selected mode to a nearby wide-bandwidth mode within a tolerance is stricter
      but is a knife-edge. Widening a kernel MERGES neighbouring bumps, and the merged mode sits
      between them, which is further than one kernel width from either. Measured: a genuine
      mixture with a 2% cell population, whose cell mode is resolved into two sub-bumps at the
      reference bandwidth, was refused because the merged wide mode landed 0.057 log10 from the
      selected one against a 0.054 tolerance. A criterion that turns on the third decimal place
      of a smoothing parameter is not measuring the data.

    So what is asked is the structural question, which needs no tolerance: at the wider bandwidth
    is there still at least one mode BELOW the valley and at least one ABOVE it? Modes may move,
    merge and be renumbered; what must survive is that the curve is still separated there. A
    single population whose flank was dressed as a second mode has all its wide-bandwidth modes on
    one side.

    Returns `{"assessed", "lo_side_has_wide_mode", "hi_side_has_wide_mode", "both",
    "n_wide_below", "n_wide_above", "lo_nearest", "hi_nearest", "lo_shift", "hi_shift"}`. The
    nearest wide mode to each selected mode and how far it moved are reported but NOT decided on -
    they are what a reader needs to see whether a pass was comfortable or bare.
    """
    out = {"assessed": False, "lo_side_has_wide_mode": False, "hi_side_has_wide_mode": False,
           "both": False, "n_wide_below": 0, "n_wide_above": 0,
           "lo_nearest": None, "hi_nearest": None, "lo_shift": None, "hi_shift": None}
    if _unknown(valley_log10):
        return out
    wide = [float(w) for w in wide_mode_positions_log10 if not _unknown(w)]
    out["assessed"] = True
    cut = float(valley_log10)
    below = [w for w in wide if w < cut]
    above = [w for w in wide if w > cut]
    out["n_wide_below"], out["n_wide_above"] = len(below), len(above)
    out["lo_side_has_wide_mode"] = bool(below)
    out["hi_side_has_wide_mode"] = bool(above)
    out["both"] = bool(below and above)
    if wide:
        for key, pos in (("lo", mode_lo_log10), ("hi", mode_hi_log10)):
            if _unknown(pos):
                continue
            nearest = min(wide, key=lambda w: abs(w - float(pos)))
            out[key + "_nearest"] = nearest
            out[key + "_shift"] = abs(nearest - float(pos))
    return out


def assess_bimodality(selection: dict, n_modes_reference: int, n_modes_wide: int, *,
                      min_valley_depth: float, min_mode_separation_log10: float,
                      mass: dict = None, persistence: dict = None,
                      min_mode_mass: float = 0.0, min_mode_observations: int = 0,
                      min_mode_spread: float = 0.0,
                      smaller_mode_iqr=None) -> tuple:
    """Reduce the measured shape to bimodal yes/no, with every criterion reported separately.

    Pure. Returns `(bimodal, criteria, reasons)`; `bimodal` is the conjunction of the criteria.

    A criterion whose input was UNDEFINED - a depth that could not be computed, a mass nobody
    passed - is False here, never True and never dropped: "not established" and "measured and
    passed" must not collapse into the same answer. Which of the two it was is stated in
    `reasons`, so a reader is never left to infer it from a bare False. That is why `mass` and
    `persistence` default to None and a None fails: a caller who does not measure them does not
    thereby get a pass on them.

    The six criteria, and which defect each one exists for:

      two_modes                    there is more than one interior maximum at all
      valley_between_modes         a grid point exists between the two tallest, so there is a
                                   minimum to report
      both_modes_carry_mass        each side of the valley holds at least `min_mode_mass` of the
                                   observations AND at least `min_mode_observations` of them.
                                   THIS IS THE ONE THAT REJECTS A KERNEL BUMP OVER A HANDFUL OF
                                   EXTREME BARCODES, which passed every other criterion here
      valley_survives_wider_bandwidth
                                   at the wider kernel there is still a mode on EACH SIDE of the
                                   valley. The previous version of this criterion counted modes
                                   at the wider bandwidth and passed whenever two of anything
                                   remained - including two bumps both sitting on the same side
                                   of the valley it had just reported
      valley_deep_enough           the dip is a dip and not a shoulder
      modes_separated_enough       the modes are further apart than one population's ripple
    """
    c = {
        "two_modes": bool(n_modes_reference >= 2),
        "valley_between_modes": not _unknown(selection.get("valley_index")),
        "both_modes_carry_mass": False,
        "valley_survives_wider_bandwidth": False,
        "valley_deep_enough": False,
        "modes_separated_enough": False,
    }
    reasons = []
    if not c["two_modes"]:
        reasons.append(f"{n_modes_reference} interior mode(s) at the reference bandwidth")
    if not c["valley_between_modes"] and selection.get("reason"):
        reasons.append(selection["reason"])

    if not isinstance(mass, dict) or _unknown(mass.get("mass_below")):
        reasons.append("mode mass NOT MEASURED - the share of barcodes under each mode was not "
                       "supplied, so it is not established rather than established and passed")
    else:
        n_below, n_above = int(mass["n_below"]), int(mass["n_above"])
        m_below, m_above = float(mass["mass_below"]), float(mass["mass_above"])
        smaller = min(m_below, m_above)
        smaller_n = min(n_below, n_above)
        c["smaller_mode_is_dispersed"] = bool(
            smaller_mode_iqr is not None and smaller_mode_iqr >= float(min_mode_spread))
        c["both_modes_carry_mass"] = bool(
            smaller >= min_mode_mass
            and smaller_n >= int(min_mode_observations)
            and c["smaller_mode_is_dispersed"])
        if not c["both_modes_carry_mass"]:
            reasons.append(
                f"the smaller mode holds {smaller_n:,} barcode(s), {100 * smaller:.4f}% of the "
                f"{int(mass['n_total']):,} used, under the {100 * min_mode_mass:.4f}% / "
                f"{int(min_mode_observations)}-barcode minimum - that is a kernel placed over a "
                f"few extreme points, not a population (below the valley {n_below:,}, above "
                f"{n_above:,})")

    if not isinstance(persistence, dict) or not persistence.get("assessed"):
        reasons.append("bandwidth persistence NOT MEASURED - not established rather than passed")
    else:
        c["valley_survives_wider_bandwidth"] = _is_true(persistence.get("both"))
        if not c["valley_survives_wider_bandwidth"]:
            empty = [n for n, k in (("below", "lo_side_has_wide_mode"),
                                    ("above", "hi_side_has_wide_mode"))
                     if not persistence.get(k)]
            reasons.append(
                f"widening the kernel leaves all {n_modes_wide} of its mode(s) on one side of the "
                f"valley - nothing {' or '.join(empty)} it "
                f"({persistence.get('n_wide_below')} below, {persistence.get('n_wide_above')} "
                f"above) - so the separation is a feature of the smoothing, not of the data")

    depth = selection.get("depth")
    if _unknown(depth):
        reasons.append("valley depth UNDEFINED - not measured and failed, but not measurable")
    else:
        depth = float(depth)
        c["valley_deep_enough"] = bool(depth >= min_valley_depth)
        if not c["valley_deep_enough"]:
            reasons.append(f"valley is {100 * depth:.1f}% below the shorter mode, under the "
                           f"{100 * min_valley_depth:.1f}% required - a shoulder, not a dip")

    sep = selection.get("separation_log10")
    if _unknown(sep):
        reasons.append("mode separation UNDEFINED - fewer than two modes to separate")
    else:
        sep = float(sep)
        c["modes_separated_enough"] = bool(abs(sep) >= min_mode_separation_log10)
        if not c["modes_separated_enough"]:
            reasons.append(f"modes are {abs(sep):.3f} log10 apart ({10 ** abs(sep):.2f}x), under "
                           f"the {min_mode_separation_log10:.3f} required - one population's "
                           f"ripple, not two populations")
    return all(c.values()), c, reasons


def find_valley(values, metric, *, bw_method="scott", grid_size=512, min_valley_depth=0.10,
                min_mode_separation_log10=0.30, bw_stability_factor=1.5, min_positive=200,
                min_mode_mass=0.001, min_mode_observations=8, min_mode_spread=0.10, max_points=None, seed=0) -> tuple:
    """The density valley of a per-barcode metric, and whether there are two modes to justify it.

    Returns `(value, bimodal, diagnostics)`. `value` is in the metric's OWN units - UMI, or
    detected genes - because that is the scale `quality.UMI_BOUNDS` and `quality.GENE_BOUNDS` are
    written in. It is the grid point of the minimum, so its precision is the grid spacing, which
    is reported as `grid_step_umi_at_valley`.

    THE PROCEDURE

    A Gaussian kernel density estimate (scipy.stats.gaussian_kde) is fitted to log10 of the
    positive values and evaluated on a linear grid of `grid_size` points spanning the observed
    range. Log10 because the two populations this separates - near-empty droplets and cells - are
    separated by orders of magnitude, and on a linear axis the lower one occupies a handful of
    pixels. The reported density is therefore a density with respect to log10(metric); a figure
    drawing it must label the axis that way or the area under the curve means nothing.

    Modes are the interior local maxima of that curve, the valley is the minimum between the two
    TALLEST modes, and the value returned is 10 ** grid[valley].

    THE BIMODALITY TEST, EXACTLY

    SIX conditions, all of which must hold. Each is reported individually in
    `diagnostics["criteria"]` so a refusal can be read rather than guessed at. This is not a
    hypothesis test and none of the numbers below is a p-value; they are shape thresholds, and
    the two marked MEASURED HERE were set against the synthetic cases in the calibration note at
    the end of this docstring:

      1. `two_modes` - two interior modes exist at the reference bandwidth;
      2. `valley_between_modes` - a grid point lies strictly between them, so a minimum exists;
      3. `both_modes_carry_mass` - counting the actual observations either side of the valley,
         the SMALLER side holds at least `min_mode_mass` of them (default 0.001, i.e. 0.1%) AND
         at least `min_mode_observations` of them (default 8). MEASURED HERE, and this is the
         criterion that rejects a mode that is a kernel bump over a handful of extreme barcodes.
         It is a count of data, not a density height: a Gaussian KDE puts a kernel of the same
         width on every point, so three barcodes together in a tail make a local maximum whose
         height says nothing about how many observations are under it;
      4. `valley_survives_wider_bandwidth` - re-estimating the density with the bandwidth
         multiplied by `bw_stability_factor`, there is still at least one mode BELOW the reported
         valley and at least one ABOVE it. Read `modes_persist()` for why this is not the same
         question as "are there still two modes" (which is what this criterion used to ask, and
         which passed on the unimodal case below), and for the tolerance-matching version that
         was tried and rejected for refusing a genuine 2% cell population;
      5. `valley_deep_enough` - relative valley depth `1 - d(valley) / min(d(mode_lo), d(mode_hi))`
         is at least `min_valley_depth`. This separates a dip from a shoulder: the flank of a
         single mode has a depth near zero against its own peak;
      6. `modes_separated_enough` - the modes are at least `min_mode_separation_log10` apart.

    WHY THE MASS CRITERION EXISTS, AND WHAT WAS MEASURED TO SET IT

    Without it this function returned `bimodal=True` on a STRICTLY UNIMODAL input - the exact
    input it exists to refuse. On 200,000 draws from a log-normal (mu=2.6, sigma=0.55 in log10,
    rounded to integers) with as few as TWO extra barcodes appended at 300,000 UMI, criteria 1,
    2, 5 and 6 all passed and the old form of criterion 4 passed as well, and a valley of 3.4 UMI
    was returned as a measurement. Counting observations shows what the density hid: 18 barcodes,
    0.009% of the sample, lay below that valley. With the mass criterion the same input is
    refused, on that count, in `reasons`.

    The 0.1% default is a judgement and it is a trade. Measured on synthetic mixtures of a debris
    mode and a cell mode 1.6 log10 apart, a cell population at 1% of barcodes carries ~0.96% mass
    and passes with an order of magnitude to spare; the spurious tail modes above carry 0.001% to
    0.01% and fail by two orders. 0.1% sits between them, nearer the spurious end so that a
    genuine but small population is not thrown away. The absolute floor of 30 observations covers
    the other regime: on a 2,000-barcode object 0.1% is two barcodes, which is precisely the
    artefact being excluded, and 30 is the point below which the share itself has a Poisson
    relative standard error above ~18% and is not a measurement of anything either. Both are
    parameters; a caller who moves them is moving the criterion that does the real work and
    should say so where the run is recorded.

    What the six criteria together do and do not achieve, measured rather than asserted. Over 216
    synthetic unimodal libraries (n in 2,000-100,000, sigma 0.35-1.1, six seeds, integer-rounded
    and not, alone and with 3 or 25 SCATTERED extreme barcodes appended) none is reported bimodal.
    Over 96 synthetic two-population mixtures none is refused by the mass or persistence criteria
    alone - every mixture that is refused is refused by the depth or mode-count criteria that were
    already here, on gaps of 1.2 log10 where the two populations genuinely overlap. The one shape
    that still passes is a tight cluster of CO-LOCATED extreme barcodes big enough to clear the
    thresholds; see the last entry under WHAT THIS TEST CANNOT DETECT

  A SMALL POPULATION AND A KERNEL BUMP CONVERGE AS THE OBJECT GETS SMALLER. Twenty-five barcodes
  sitting far above the mode are 0.0125% of a 200,000-barcode object and 1.25% of a 2,000-barcode
  one. At 1.25% they are indistinguishable, by mass or by dispersion, from twenty-five genuinely
  very deep nuclei - because on the evidence available here they ARE indistinguishable. This
  function reports bimodal=True in that case, and it should: the density does have two modes.

  That is not the end of the check. The valley it returns then sits far out in the tail - tens of
  thousands of UMI - and `modules/05_quality/quality.py` refuses any valley outside its sanity
  bounds. The two checks are independent and neither subsumes the other: this one asks whether
  the distribution has two modes, and the bound asks whether the resulting threshold could
  plausibly be a debris/nucleus boundary. A caller that uses find_valley WITHOUT that bound has
  removed the check that catches this case.

  Also undetectable here: a genuinely tight real population, whose members sit close together by
  biology rather than by artifact, is refused by the dispersion criterion. That is the cost of
  the criterion and it is paid deliberately - a false refusal is visible and recoverable, a false
  threshold is neither.

    WHAT THIS TEST CANNOT DETECT - read this before quoting a valley

    - **It is not a hypothesis test.** There is no null distribution, no p-value and no false
      positive rate. It reports that a curve has the SHAPE of two modes, not that a unimodal
      density would be unlikely to produce that shape. Hartigan's dip test and Silverman's
      critical-bandwidth bootstrap do give a p-value; both need machinery this file deliberately
      does not carry, and neither has been run on this pipeline's data.
    - **Anything finer than the bandwidth is invisible.** Two populations closer together than
      the kernel width are smoothed into one mode and reported unimodal. The bandwidth is
      reported as `bandwidth` so this bound is at least visible.
    - **No single criterion above is sufficient, and DEPTH AND SEPARATION ARE THE WEAKEST.**
      Relative depth is measured against the SHORTER of the two modes, so a bump far out in a
      tail has a tiny reference height and the dip in front of it computes as deep; separation is
      then automatically satisfied because the bump is far away. On the synthetic unimodal
      log-normal described above, both passed on an artefact.
    - **Bandwidth persistence is a weak criterion and is not what previous versions of this
      docstring claimed.** It was described as the criterion doing most of the work; it was not
      doing that work. What it tested was that TWO modes of any kind survived a wider kernel, and
      on the unimodal case above two did - the main mode and the tail bump - so it passed, which
      is how a strictly unimodal input came back bimodal. It now asks whether the wider kernel
      still has a mode on each side of the valley. That is a real question, and on the cases
      measured here it does reject the unimodal artefacts; it is still weak, because widening a
      Gaussian kernel does not remove points, so a bump over a handful of co-located barcodes
      survives widening perfectly well and only needs one other bump on the far side to pass.
      Treat criterion 3, the mass, as the one doing the work, and look at it first when this
      function is asked why it refused - or why it did not.
    - **The mass criterion bounds how small a population can be seen, and nothing more.** A real
      subpopulation below `min_mode_mass` of the barcodes is reported unimodal. That is a refusal
      to measure, not a finding of absence, and `mass_below` / `mass_above` are in the
      diagnostics so the reader can see how close it was.
    - **The mass criterion counts barcodes; it cannot ask whether they are a population.** A
      cluster of co-located barcodes that clears both thresholds is reported as a mode whatever it
      is - a spike-in, a clipped value, a plate of debris. Measured, so the boundary is not
      guessed at: 25 barcodes at exactly 300,000 UMI added to 20,000 log-normal draws come back
      `bimodal=True` with the valley at 182,000 UMI (31 barcodes above it, 0.155% - past both
      thresholds); the same 25 added to 200,000 draws come back False (0.0125%). Nothing here can
      tell those two apart, and nothing here is meant to: what the valley IS remains an
      interpretation. The number is bounds-checked downstream by `quality.UMI_BOUNDS`, which is
      where a valley of 182,000 UMI is caught, and this is why that check exists as well as this
      one.
    - **Overlapping populations with one mode cannot be found by any density method.** Two
      populations sharing a mean and differing in spread produce a unimodal mixture. Unimodal
      here means "one mode", never "one population".
    - **It cannot say what the two modes ARE.** Debris against nuclei, two cell types of very
      different RNA content, or two libraries of different depth pooled into one object all give
      two modes. That the lower mode is empty droplets is an interpretation, not this output.
    - **It cannot see what was removed before it ran.** A matrix that was already cell-called has
      had its lower mode deleted; the remaining curve can still be bimodal for an unrelated
      reason, and the valley then means something else entirely. `lib/verify_raw.py` is the check
      for that, and it runs at step 0 for this reason. A mode pinned against the edge of the
      observed range is reported here as `mode_at_left_edge` / `mode_at_right_edge` - suggestive,
      not conclusive, and deliberately NOT folded into `bimodal`.
    - **Non-positive and non-finite values are outside the log transform** and are excluded from
      the density, counted, and reported as `n_nonpositive_excluded` / `n_nonfinite_excluded`.
      They are not removed from anything; this function writes no data.
    - **The value moves with the bandwidth rule.** `bw_method` is Scott's rule by default;
      Silverman's gives a different valley on the same data. The choice is recorded, not hidden.

    DETERMINISM

    Nothing here is stochastic unless `max_points` is set, which subsamples with a seeded
    generator before FITTING; the subsample size and `seed` are then recorded in the diagnostics.
    The mass in criterion 3 is counted over every positive value regardless, because the question
    it asks is how many barcodes there are, and a subsample would answer a different one.
    Leaving `max_points` at None uses every value everywhere and is exact.
    """
    import numpy as np
    from scipy.stats import gaussian_kde

    if grid_size < 8:
        raise TaskFailure(f"grid_size={grid_size} is too coarse to locate a minimum on")

    # Not np.asarray(..., dtype=float): an object column carrying pd.NA raises there, and
    # numpy.ma.masked is coerced to its fill value, which is how a sentinel becomes a barcode with
    # a count. _float_array maps every unknown to NaN, where the non-finite branch below counts it.
    v = _float_array(values).ravel()
    n_total = int(v.size)
    finite = np.isfinite(v)
    n_nonfinite = int((~finite).sum())
    positive = finite & (v > 0)
    n_nonpositive = int((finite & ~positive).sum())
    x = np.log10(v[positive])
    n_used = int(x.size)

    if n_used < min_positive:
        raise TaskFailure(
            f"{metric}: only {n_used} of {n_total} barcodes carry a positive value; a kernel "
            f"density estimate on that many points is not a measurement of anything. "
            f"({n_nonpositive} non-positive, {n_nonfinite} non-finite.) Check that the right "
            f"column was passed and that the matrix is the unfiltered droplet matrix.")
    if float(x.max() - x.min()) <= 0:
        raise TaskFailure(f"{metric}: every barcode has the same value ({10 ** float(x[0]):.6g}); "
                          f"there is no distribution to find a valley in")

    x_all = x                       # every positive value; the mass criterion counts on this
    subsampled = None
    if max_points is not None and n_used > max_points:
        rng = np.random.default_rng(seed)
        keep = rng.choice(n_used, size=int(max_points), replace=False)
        keep.sort()
        x = x[keep]
        subsampled = int(max_points)

    kde = gaussian_kde(x, bw_method=bw_method)
    grid = np.linspace(float(x.min()), float(x.max()), int(grid_size))
    dens = np.asarray(kde(grid), dtype=float)
    if not np.all(np.isfinite(dens)):
        raise TaskFailure(f"{metric}: the density estimate returned non-finite values; the valley "
                          f"cannot be located and must not be guessed at")

    peaks, troughs = find_modes(dens.tolist())
    sel = valley_between(grid.tolist(), dens.tolist(), peaks)

    wide = gaussian_kde(x, bw_method=float(kde.factor) * float(bw_stability_factor))
    dens_wide = np.asarray(wide(grid), dtype=float)
    peaks_wide, _ = find_modes(dens_wide.tolist())

    vi = sel["valley_index"]
    value = float(10 ** grid[vi]) if vi is not None else None
    step = float(grid[1] - grid[0])

    # Criterion 3, counted on the observations rather than read off the density. x_all, not x:
    # a subsample would answer "how many of the points I fitted", which is not the question.
    mass = mode_masses(x_all.tolist(), (float(grid[vi]) if vi is not None else None))

    # Criterion 4. The wide kernel's own sd on the log10 axis is recorded beside it because it is
    # the scale at which the wide curve can resolve anything, but the criterion does not compare
    # against it - see modes_persist() for why a tolerance was tried and rejected.
    bw_wide_log10 = float(wide.factor) * float(np.std(x, ddof=1))
    persistence = modes_persist(
        (float(grid[vi]) if vi is not None else None),
        (float(grid[sel["mode_lo_index"]]) if sel["mode_lo_index"] is not None else None),
        (float(grid[sel["mode_hi_index"]]) if sel["mode_hi_index"] is not None else None),
        [float(grid[i]) for i in peaks_wide])

    bimodal, criteria, reasons = assess_bimodality(
        sel, len(peaks), len(peaks_wide),
        min_valley_depth=min_valley_depth,
        min_mode_separation_log10=min_mode_separation_log10,
        mass=mass, persistence=persistence,
        min_mode_mass=float(min_mode_mass),
        min_mode_observations=int(min_mode_observations),
        min_mode_spread=float(min_mode_spread),
        # Criterion 5, on the same observations and the same cut criterion 3 counts.
        smaller_mode_iqr=mode_spread(x_all.tolist(),
                                     (float(grid[vi]) if vi is not None else None),
                                     "smaller"))

    diagnostics = {
        "metric": metric,
        "n_input": n_total,
        "n_used": int(x.size),
        "n_nonpositive_excluded": n_nonpositive,
        "n_nonfinite_excluded": n_nonfinite,
        "subsampled_to": subsampled,
        "seed": seed if subsampled is not None else None,
        "bw_method": bw_method if isinstance(bw_method, str) else float(bw_method),
        "bandwidth": float(kde.factor),
        "bandwidth_wide": float(wide.factor),
        "bw_stability_factor": float(bw_stability_factor),
        "grid_size": int(grid_size),
        "grid_step_log10": step,
        # The grid is linear in log10, so its spacing in the metric's own units depends on where
        # you are on it. Reported AT THE VALLEY, which is the only place the number is used.
        "grid_step_umi_at_valley": (float(value * (10 ** step - 1)) if value is not None
                                    else None),
        "grid_log10": grid.tolist(),
        "grid": (10 ** grid).tolist(),
        "density": dens.tolist(),
        "density_wide": dens_wide.tolist(),
        "density_is_wrt": "log10(" + str(metric) + ")",
        "mode_indices": list(peaks),
        "mode_values": [float(10 ** grid[i]) for i in peaks],
        "trough_indices": list(troughs),
        "mode_indices_wide": list(peaks_wide),
        "mode_lo_index": sel["mode_lo_index"],
        "mode_hi_index": sel["mode_hi_index"],
        "valley_index": vi,
        "valley_log10": (float(grid[vi]) if vi is not None else None),
        "depth": sel["depth"],
        "separation_log10": sel["separation_log10"],
        # Counts, not densities. n_below + n_above == n_used, and the pair is what criterion 3
        # is decided on, so a reader can recompute the decision from the table.
        "n_below_valley": mass["n_below"],
        "n_above_valley": mass["n_above"],
        "mass_below": mass["mass_below"],
        "mass_above": mass["mass_above"],
        "mass_counted_on": "all positive values, not the KDE subsample",
        "persistence": persistence,
        "bandwidth_wide_log10": bw_wide_log10,
        "mode_at_left_edge": bool(dens[0] > dens[1]),
        "mode_at_right_edge": bool(dens[-1] > dens[-2]),
        "criteria": criteria,
        "reasons": reasons,
        "thresholds": {"min_valley_depth": float(min_valley_depth),
                       "min_mode_separation_log10": float(min_mode_separation_log10),
                       "min_mode_mass": float(min_mode_mass),
                       "min_mode_observations": int(min_mode_observations),
                       "bw_stability_factor": float(bw_stability_factor),
                       "min_positive": int(min_positive)},
        "range_log10": [float(x.min()), float(x.max())],
    }
    return value, bool(bimodal), diagnostics


def valley_note(diagnostics: dict) -> str:
    """One line for `quality.Valley.note`: what was measured and, if it failed, which criterion.

    Written from the diagnostics rather than from the caller's memory of them, so the note and
    the number cannot disagree.
    """
    bits = [f"KDE bw={diagnostics['bandwidth']:.4g} ({diagnostics['bw_method']}) on "
            f"{diagnostics['n_used']:,} of {diagnostics['n_input']:,} barcodes",
            f"{len(diagnostics['mode_indices'])} mode(s)"]
    d = diagnostics.get("depth")
    bits.append(f"depth {100 * d:.1f}%" if not _unknown(d) else "depth UNDEFINED")
    mb, ma = diagnostics.get("mass_below"), diagnostics.get("mass_above")
    if _unknown(mb) or _unknown(ma):
        bits.append("mode mass UNMEASURED")
    else:
        # The counts as well as the shares: 0.01% and 18 barcodes are the same fact, and the one
        # a reader can act on is the count.
        bits.append(f"mass {100 * mb:.3f}%/{100 * ma:.3f}% "
                    f"({diagnostics.get('n_below_valley'):,}/"
                    f"{diagnostics.get('n_above_valley'):,} barcodes) either side")
    if diagnostics["n_nonpositive_excluded"]:
        bits.append(f"{diagnostics['n_nonpositive_excluded']:,} non-positive excluded from the "
                    f"log transform")
    if diagnostics["mode_at_left_edge"]:
        bits.append("a mode sits against the LEFT edge of the observed range - check the input "
                    "was not already cell-called")
    if diagnostics["reasons"]:
        bits.append("NOT bimodal: " + "; ".join(diagnostics["reasons"]))
    return "; ".join(bits)


# --------------------------------------------------------------------------------------------
# per-cell QC metrics
# --------------------------------------------------------------------------------------------

def qc_metrics(adata, mt_prefix, ribo_pattern, *, allow_empty_mt=False, allow_empty_ribo=False,
               layer=None, add_log1p=False, allow_transformed=False):
    """Per-cell QC via `sc.pp.calculate_qc_metrics`, with the matched gene lists recorded.

    Adds `total_counts`, `n_genes_by_counts`, `pct_counts_mt` and `pct_counts_ribo` to `obs`, and
    the boolean `mt` / `ribo` columns to `var`. Returns the same object, modified in place.

    GENE-CLASS PATTERNS ARE PARAMETERS WITH NO DEFAULT

    `mt_prefix` and `ribo_pattern` are species-specific and this module will not choose them.
    Mouse and human differ in case alone for the mitochondrial prefix, which is precisely the kind
    of difference that produces a plausible zero rather than an error. Matching is CASE-SENSITIVE
    for that reason: a case-insensitive `MT` prefix also matches MTOR.

    `mt_prefix` is a literal prefix. `ribo_pattern` is a regular expression matched with re.search
    semantics, so anchor it with `^` if you mean the start of the symbol - and print the list it
    produces before trusting it. A pattern intended for ribosomal genes has matched a ribosomal
    protein KINASE before now (docs/PRINCIPLES.md section 1). The symbols matched by each pattern
    are recorded in full in `adata.uns["scqc_qc_metrics"]`, not summarised as a count, because a
    count cannot be audited.

    A PATTERN THAT MATCHES NOTHING IS REFUSED

    Zero matched mitochondrial genes gives every cell `pct_counts_mt == 0`, which is
    indistinguishable on the page from a clean library and passes every mitochondrial gate in this
    pipeline. That is the failure in docs/PRINCIPLES.md section 4 with a different variable name,
    so it raises. `allow_empty_mt` / `allow_empty_ribo` exist for a reference that genuinely lacks
    the class; they are recorded in `uns` when used, because a bypass that leaves no trace is a
    bypass nobody can review.
    """
    import pandas as pd
    import scanpy as sc

    if not isinstance(mt_prefix, str) or not mt_prefix:
        raise TaskFailure(
            "mt_prefix must be a non-empty string and has no default: it is species-specific "
            "(mouse 'mt-', human 'MT-'), and a wrong or absent prefix reports 0% mitochondrial "
            "content for every cell, which passes every gate downstream.")
    if not isinstance(ribo_pattern, str) or not ribo_pattern:
        raise TaskFailure(
            "ribo_pattern must be a non-empty regular expression and has no default: it is "
            "species-specific, and an empty pattern would report every gene or no gene as "
            "ribosomal without either being visible.")
    try:
        re.compile(ribo_pattern)
    except re.error as e:
        raise TaskFailure(f"ribo_pattern {ribo_pattern!r} is not a valid regular expression: {e}")

    if layer is not None and layer not in adata.layers:
        raise TaskFailure(f"layer {layer!r} is not in this object. Layers present: "
                          f"{', '.join(map(str, adata.layers.keys())) or '(none)'}")
    matrix = adata.layers[layer] if layer is not None else adata.X
    notes = _refuse_if_transformed(matrix, "qc_metrics", allow_transformed=allow_transformed)

    names = pd.Index(adata.var_names).astype(str)
    # `.str.startswith` returns NaN (or pd.NA on a nullable dtype) for a missing symbol, and a
    # mask carrying either is not a mask: numpy indexing with it either raises or, worse, treats
    # it as True and puts an unnamed gene into a reported gene class. Unknown means NOT matched,
    # and the count of unknowns is recorded below rather than left to be inferred.
    mt_mask, n_mt_unknown = _bool_array(names.str.startswith(mt_prefix))
    ribo_mask, n_ribo_unknown = _bool_array(
        names.str.contains(ribo_pattern, regex=True, na=False))
    mt_genes = [str(g) for g in names[mt_mask]]
    ribo_genes = [str(g) for g in names[ribo_mask]]

    if not mt_genes and not allow_empty_mt:
        raise TaskFailure(
            f"mt_prefix {mt_prefix!r} matched 0 of {len(names):,} genes. Every cell would be "
            f"reported at 0% mitochondrial content, which is what a clean library looks like and "
            f"what a wrong prefix looks like. First few symbols in this reference: "
            f"{', '.join(map(str, names[:8]))}. Pass allow_empty_mt=True only if this reference "
            f"genuinely has no mitochondrial genes.")
    if not ribo_genes and not allow_empty_ribo:
        raise TaskFailure(
            f"ribo_pattern {ribo_pattern!r} matched 0 of {len(names):,} genes; matching uses "
            f"re.search, so a pattern written for a different species' capitalisation matches "
            f"nothing and reports 0% rather than failing. First few symbols in this reference: "
            f"{', '.join(map(str, names[:8]))}. Pass allow_empty_ribo=True to declare the "
            f"absence deliberate.")

    # A RIBOSOMAL PATTERN THAT CATCHES KINASES IS THE ONE MISTAKE THIS METRIC HAS. `^Rp[sl]`
    # also matches Rps6ka1-6, Rps6kb1-2, Rps6kc1 and Rps6kl1 - the S6 kinases, which are mTOR
    # SIGNALLING genes and not ribosomal proteins. This tool's own documentation suggested that
    # pattern, so every samplesheet copied from it counted ten signalling genes as ribosome.
    #
    # It is not caught by the empty-match check above, and it never fails: it returns a
    # plausible percentage that is slightly wrong in a direction nobody can see. Named here
    # rather than left to the reader, because a count cannot be audited and a percentage cannot
    # be un-published.
    _kin = [g for g in ribo_genes if re.match(r"^(?:Rps6k|RPS6K)", str(g))]
    if _kin:
        notes.append(
            f"ribo_pattern {ribo_pattern!r} matched {len(_kin)} S6 KINASE(S) - "
            f"{', '.join(map(str, _kin[:6]))}{' ...' if len(_kin) > 6 else ''} - which are mTOR "
            f"signalling genes, not ribosomal proteins. pct_counts_ribo therefore includes them. "
            f"Use '^Rp[sl](?!6k)' (mouse) or '^RP[SL](?!6K)' (human) to exclude them.")

    adata.var["mt"] = mt_mask
    adata.var["ribo"] = ribo_mask

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None,
                               log1p=bool(add_log1p), inplace=True, layer=layer)

    # Measure the structure rather than assume it: the column names depend on `qc_vars`, and a
    # renamed column would otherwise surface three steps later as a missing threshold input.
    required = ("total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo")
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        raise TaskFailure(
            f"calculate_qc_metrics did not produce {', '.join(missing)}. It produced: "
            f"{', '.join(map(str, adata.obs.columns))}. The threshold modules read these names "
            f"directly, so a rename has to be handled here rather than absorbed downstream.")

    _refuse_if_totals_constant(_float_array(adata.obs["total_counts"]), "qc_metrics")

    if n_mt_unknown or n_ribo_unknown:
        notes.append(f"{n_mt_unknown} gene(s) gave an unknown mt match and {n_ribo_unknown} an "
                     f"unknown ribo match (a missing symbol); each was counted as NOT matched")
    adata.uns["scqc_qc_metrics"] = {
        "mt_prefix": mt_prefix,
        "ribo_pattern": ribo_pattern,
        "match_semantics": "mt: literal prefix, case-sensitive. ribo: re.search, case-sensitive",
        "n_mt_genes": len(mt_genes),
        "n_ribo_genes": len(ribo_genes),
        "n_mt_match_unknown": int(n_mt_unknown),
        "n_ribo_match_unknown": int(n_ribo_unknown),
        "mt_genes": mt_genes,          # the LIST, not a count: a count cannot be audited
        "ribo_genes": ribo_genes,
        "layer": layer,
        "allow_empty_mt": bool(allow_empty_mt),
        "allow_empty_ribo": bool(allow_empty_ribo),
        "allow_transformed": bool(allow_transformed),
        "notes": notes,
        "versions": observed_versions(),
    }
    return adata


# --------------------------------------------------------------------------------------------
# clustering
# --------------------------------------------------------------------------------------------

def cluster(adata, resolution, seed, n_hvg=2000, n_pcs=50, n_neighbors=15, *,
            key_added=None, hvg_flavor="seurat", target_sum=1e4, require_qc_metrics=True,
            allow_transformed=False):
    """normalize_total -> log1p -> HVG -> PCA -> neighbours -> leiden, with nothing excluded.

    Returns the same object, modified IN PLACE: `adata.X` is normalised and logged, which is what
    the reference chain does. Keep the counts elsewhere before calling this if you need them - a
    layer, or the file on disk.

    `key_added` defaults to `leiden_res<resolution>` rather than to `leiden`, so a resolution
    sweep does not overwrite its own earlier answer under one column name. The resolved key is
    recorded in `adata.uns["scqc_cluster"]["key"]` and is what `cluster_profile()` should be given.

    NO GENE CLASS IS EXCLUDED FROM SELECTION, AND NONE MAY BE ADDED

    `highly_variable_genes` is computed over every gene in the object and the resulting flag is
    used as scanpy produced it. There is no exclusion term here. What is produced instead is a
    report: the symbols of the flagged mitochondrial and ribosomal genes that selection CHOSE, in
    `uns["scqc_cluster"]["flagged_hvg"]`, so their influence can be read and adjudicated. If
    `var["mt"]` / `var["ribo"]` are absent the report is None - never an empty list, which would
    read as "examined, none selected".

    Selection is also the reversible form of the removal question (docs/PRINCIPLES.md section 1,
    question 4): no gene is dropped from the object, so every gene remains testable no matter what
    the embedding used.

    ORDERING

    `qc_metrics()` must have run first. After normalize_total every cell's counts sum to
    `target_sum`, so a `total_counts` computed afterwards is a constant and every depth statistic
    derived from it is meaningless while still printing. That is checked, not trusted.
    """
    import numpy as np
    import scanpy as sc

    if require_qc_metrics and "total_counts" not in adata.obs.columns:
        raise TaskFailure(
            "obs['total_counts'] is absent, so qc_metrics() has not run. Clustering normalises X "
            "in place; computing depth afterwards yields a constant that still prints as a "
            "number. Run qc_metrics() first, or pass require_qc_metrics=False and be explicit "
            "about where the depth statistics are coming from.")
    _refuse_if_transformed(adata.X, "cluster", allow_transformed=allow_transformed)
    if "total_counts" in adata.obs.columns and not allow_transformed:
        # The same check qc_metrics makes, from the other side: an object that arrives already
        # normalised carries a constant total_counts, and this is the last point at which that is
        # still visible - after normalize_total below it is true of every object.
        _refuse_if_totals_constant(_float_array(adata.obs["total_counts"]), "cluster")

    key = key_added if key_added else f"leiden_res{float(resolution):g}"

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=int(n_hvg), flavor=hvg_flavor)

    if "highly_variable" not in adata.var.columns:
        raise TaskFailure(
            f"highly_variable_genes(flavor={hvg_flavor!r}) produced no 'highly_variable' column. "
            f"var columns present: {', '.join(map(str, adata.var.columns))}")
    # `np.asarray(col, dtype=bool)` on a nullable boolean column raises on pd.NA, and on an object
    # column it makes every non-empty sentinel True - which would put a gene into the embedding
    # because its flag was missing. Unknown is not selected, and the count is recorded.
    hv, n_hv_unknown = _bool_array(adata.var["highly_variable"])
    n_hv = int(hv.sum())
    if n_hv == 0:
        raise TaskFailure("no genes were selected as highly variable; PCA has nothing to run on. "
                          "Check that the object holds more than one cell type's worth of cells "
                          "and that the matrix is counts.")

    flagged = {}
    for cls in ("mt", "ribo"):
        if cls in adata.var.columns:
            m, n_cls_unknown = _bool_array(adata.var[cls])
            chosen = [str(g) for g in adata.var_names[m & hv]]
            flagged[cls] = {"n_in_reference": int(m.sum()),
                            "n_selected": len(chosen),
                            "n_flag_unknown": int(n_cls_unknown),
                            "selected": chosen}
        else:
            # NOT an empty list. "Never examined" and "examined, none selected" are different
            # findings and must not print the same way.
            flagged[cls] = None

    max_comps = min(int(adata.n_obs), n_hv) - 1
    if int(n_pcs) > max_comps:
        raise TaskFailure(
            f"n_pcs={n_pcs} but arpack can return at most {max_comps} components here "
            f"({adata.n_obs:,} cells, {n_hv:,} highly-variable genes). Lower n_pcs or raise "
            f"n_hvg; silently reducing it would make two runs of the same pipeline "
            f"incomparable without saying so.")

    sc.tl.pca(adata, n_comps=int(n_pcs), svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=int(n_neighbors), n_pcs=int(n_pcs), random_state=seed)
    try:
        sc.tl.leiden(adata, resolution=float(resolution), random_state=seed, key_added=key)
    except ImportError as e:
        raise TaskFailure(
            f"leiden clustering is unavailable in this environment: {e}. scanpy needs leidenalg "
            f"or python-igraph installed; scQC does not substitute a different community "
            f"detection algorithm, because the cluster labels would not be the ones the report "
            f"describes.")

    if key not in adata.obs.columns:
        raise TaskFailure(f"leiden did not write obs[{key!r}]. obs columns: "
                          f"{', '.join(map(str, adata.obs.columns))}")
    labels = adata.obs[key].astype(str)
    n_clusters = int(labels.nunique())

    # Which leiden backend scanpy resolved to. scanpy has changed this default across versions
    # and the two backends do not return identical partitions, so the resolved value belongs in
    # the record beside the resolution rather than in a reader's assumptions.
    try:
        import inspect
        leiden_flavor = inspect.signature(sc.tl.leiden).parameters["flavor"].default
    except Exception:                                       # noqa: BLE001 - absence is the fact
        leiden_flavor = "not exposed by this scanpy version"

    pca_var = None
    try:
        pca_var = float(np.asarray(adata.uns["pca"]["variance_ratio"]).sum())
    except Exception:                                       # noqa: BLE001 - optional diagnostic
        pca_var = None
    # Whether PCA restricted itself to the highly-variable genes is scanpy's DEFAULT, and the
    # argument that controls it has been renamed across versions. Recording what scanpy resolved
    # to is cheaper than asserting it and being wrong quietly: if this says the mask was not
    # applied, the embedding used every gene and n_hvg meant nothing.
    try:
        pca_params = {str(k): str(v) for k, v in dict(adata.uns["pca"]["params"]).items()}
    except Exception:                                       # noqa: BLE001 - optional diagnostic
        pca_params = None

    adata.uns["scqc_cluster"] = {
        "key": key,
        "algorithm": "leiden",
        "resolution": float(resolution),
        "seed": int(seed),
        "n_hvg_requested": int(n_hvg),
        "n_hvg_selected": n_hv,
        "n_hvg_flag_unknown": int(n_hv_unknown),
        "hvg_flavor": hvg_flavor,
        "n_pcs": int(n_pcs),
        "n_neighbors": int(n_neighbors),
        "target_sum": float(target_sum),
        "n_clusters": n_clusters,
        "pca_variance_ratio_sum": pca_var,
        "pca_params_as_scanpy_resolved_them": pca_params,
        "leiden_flavor_default": str(leiden_flavor),
        "flagged_hvg": flagged,
        "exclusions": "none - every gene entered variable-gene selection",
        "versions": observed_versions(),
    }
    return adata


def embed(adata, seed, *, n_hvg=2000, n_pcs=50, n_neighbors=15, min_dist=0.5,
          hvg_flavor="seurat", target_sum=1e4, allow_transformed=False) -> tuple:
    """normalize_total -> log1p -> HVG -> PCA -> neighbours -> UMAP. Returns `(coords, record)`.

    `coords` is one `(x, y)` pair per observation IN THE ORDER `adata.obs_names` holds them, so a
    caller zips it against the barcodes rather than trusting two sorted lists to agree.

    THE COORDINATES ARE WRITTEN DOWN BECAUSE NOTHING ELSE CAN RECOVER THEM

    Figures F10 and F11 are the only two that can show whether a removal took a COHERENT REGION or
    scattered points, and both need the SAME coordinates: re-embedding the survivors gives a
    different layout, and a reader comparing before with after would be reading the projection
    rather than the data. UMAP is also not stable across versions or across a change in the input
    set, so an embedding that is not stored is not reproducible - it is regenerated, which is a
    different picture under the same name.

    THE SAME PROCEDURE AS `cluster()`, DELIBERATELY, AND NOT THE SAME CALL

    Every step up to the neighbour graph is what `cluster()` does, because a figure drawn on a
    differently-built graph would not be showing the clustering it sits beside. It is not
    `cluster()` itself because the two run over DIFFERENT POPULATIONS on purpose - see
    `_op_cluster` - and because leiden over a population nothing will use is wasted work.

    `min_dist` is a FIXED procedure parameter of the layout, not a threshold: it changes how tight
    the picture looks and changes no number anywhere. It is recorded so two runs' figures can be
    told apart.
    """
    import numpy as np
    import scanpy as sc

    _refuse_if_transformed(adata.X, "embed", allow_transformed=allow_transformed)

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=int(n_hvg), flavor=hvg_flavor)
    if "highly_variable" not in adata.var.columns:
        raise TaskFailure(
            f"highly_variable_genes(flavor={hvg_flavor!r}) produced no 'highly_variable' column, "
            f"so there is nothing to embed on. var columns present: "
            f"{', '.join(map(str, adata.var.columns))}")
    hv, _hv_unknown = _bool_array(adata.var["highly_variable"])
    n_hv = int(hv.sum())
    if n_hv == 0:
        raise TaskFailure("no genes were selected as highly variable, so PCA has nothing to run "
                          "on and there is no embedding to draw.")

    # The same arithmetic `cluster()` refuses on, and refused here too rather than quietly
    # lowering n_pcs: an embedding built on fewer components than the clustering is a different
    # projection of the same cells, and nothing on the figure would say so.
    max_comps = min(int(adata.n_obs), n_hv) - 1
    if int(n_pcs) > max_comps:
        raise TaskFailure(
            f"n_pcs={n_pcs} but arpack can return at most {max_comps} components over the "
            f"population being embedded ({adata.n_obs:,} cells, {n_hv:,} highly-variable "
            f"genes). Lower n_pcs or raise n_hvg.")

    sc.tl.pca(adata, n_comps=int(n_pcs), svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=int(n_neighbors), n_pcs=int(n_pcs), random_state=seed)
    try:
        sc.tl.umap(adata, min_dist=float(min_dist), random_state=seed)
    except ImportError as e:
        raise TaskFailure(
            f"UMAP is unavailable in this environment: {e}. scanpy needs umap-learn installed. "
            f"scQC does not substitute a different layout, because the figure's caption names "
            f"the one it describes.") from None

    if "X_umap" not in adata.obsm:
        raise TaskFailure(f"sc.tl.umap wrote no obsm['X_umap']. obsm keys: "
                          f"{', '.join(map(str, adata.obsm.keys()))}")
    xy = np.asarray(adata.obsm["X_umap"], dtype=float)
    if xy.ndim != 2 or xy.shape[0] != int(adata.n_obs) or xy.shape[1] < 2:
        raise TaskFailure(
            f"obsm['X_umap'] has shape {xy.shape}, which cannot be matched one-to-one to "
            f"{adata.n_obs:,} observations. Coordinates that do not line up with their barcodes "
            f"colour the wrong points, and there is no symptom of that on the page.")

    record = {
        "layout": "umap",
        "n_obs": int(adata.n_obs),
        "seed": int(seed),
        "n_hvg_requested": int(n_hvg),
        "n_hvg_selected": n_hv,
        "hvg_flavor": hvg_flavor,
        "n_pcs": int(n_pcs),
        "n_neighbors": int(n_neighbors),
        "min_dist": float(min_dist),
        "target_sum": float(target_sum),
        "versions": observed_versions(),
    }
    return [(float(a), float(b)) for a, b in xy[:, :2]], record


# --------------------------------------------------------------------------------------------
# the cluster profile - the rows step 6 consumes
# --------------------------------------------------------------------------------------------

def _cluster_sort_key(label: str):
    """Numeric where the label is numeric, textual otherwise, so 2 sorts before 10."""
    s = str(label)
    return (0, int(s), "") if s.isdigit() else (1, 0, s)


def _resolve_uninformative(adata, uninformative_genes, mt_key, ribo_key,
                           allow_empty_uninformative):
    """The locked uninformative set, split into its mt and ribo halves where that is knowable.

    Accepts either a mapping of class -> symbols (the split is then explicit and authoritative)
    or a flat iterable (the split is derived from `var[mt_key]` / `var[ribo_key]`, and is None if
    those columns are absent - never an empty set, which would report 0% for a quantity nobody
    computed).

    criterion C is reported SPLIT by `cluster_flags`, which says why: in the calibration cohort
    the ribosomal half was empty and folding the two together would have hidden which one fired.
    """
    import numpy as np

    mt_set = ribo_set = None
    if hasattr(uninformative_genes, "items"):
        classes = {str(k).lower(): {str(g) for g in v if not _unknown(g)}
                   for k, v in uninformative_genes.items()}
        all_set = set()
        for v in classes.values():
            all_set |= v
        mt_set = classes.get("mt")
        ribo_set = classes.get("ribo")
        source = ("mapping supplied by the caller: "
                  + ", ".join(f"{k}={len(v)}" for k, v in sorted(classes.items())))
    else:
        all_set = {str(g) for g in uninformative_genes if not _unknown(g)}
        if mt_key in adata.var.columns and ribo_key in adata.var.columns:
            names = np.asarray([str(g) for g in adata.var_names])
            mt_flag, _ = _bool_array(adata.var[mt_key])
            ribo_flag, _ = _bool_array(adata.var[ribo_key])
            mt_set = {g for g in names[mt_flag] if g in all_set}
            ribo_set = {g for g in names[ribo_flag] if g in all_set}
            source = (f"flat list of {len(all_set)}, split by var[{mt_key!r}]/var[{ribo_key!r}]")
        else:
            source = (f"flat list of {len(all_set)}; var[{mt_key!r}]/var[{ribo_key!r}] absent, so "
                      f"the mt/ribo split is UNKNOWN and reported as such")

    if not all_set and not allow_empty_uninformative:
        raise TaskFailure(
            "the locked uninformative gene set is empty. Criterion C would then read 0% for every "
            "cluster, which is indistinguishable from a measured clean result and passes the gate "
            "for the wrong reason. Pass the set, or allow_empty_uninformative=True to declare the "
            "absence deliberate.")

    in_object = {str(g) for g in adata.var_names}
    present = {g for g in all_set if g in in_object}
    if all_set and not present and not allow_empty_uninformative:
        raise TaskFailure(
            f"none of the {len(all_set)} uninformative genes appear in this object's var_names "
            f"(first few of the set: {', '.join(sorted(all_set)[:6])}; first few in the object: "
            f"{', '.join(map(str, list(adata.var_names)[:6]))}). A set locked against a different "
            f"reference reports 0% for every cluster rather than failing.")
    return all_set, mt_set, ribo_set, present, source


def attach_doublet_calls(adata, csv_path, *, key, barcode_column="barcode",
                         class_column="doublet_class", score_column="doublet_score"):
    """Put a per-barcode doublet call into `obs[key]`. Attached, never applied.

    Criterion D is computable at exactly one point in this pipeline: on an object carrying the
    doublet calls ATTACHED and NOT applied (modules/06_cluster_check/cluster_flags.py, first
    section). After a removal every cluster is 0% doublet by construction and D passes for the
    wrong reason. The detector runs in a different environment and leaves a CSV beside the object,
    so something has to bring the two together, and this is it.

    A barcode absent from the CSV was never scored - it sat below the light floor - and is left
    None rather than filled with "singlet". `_doublet_masks` then drops it from BOTH sides of D's
    ratio, instead of counting a cell nobody examined as a cell that passed.

    ZERO OVERLAP IS REFUSED, and that is the point of the function existing rather than a one-line
    merge. In this project's own cohort the object's barcodes carried a `<sample>_` prefix and a
    detector's CSV did not: the same run, described two ways, intersecting in NO barcodes. The
    merge would have succeeded and written None into every cell, and criterion D would then have
    read "not evaluated" for the entire cohort - which on the page is very hard to tell from a
    cohort that was evaluated and had no doublet-driven clusters.
    """
    import pandas as pd

    p = Path(csv_path)
    if not p.exists():
        raise TaskFailure(
            f"doublet_csv {p} does not exist. Criterion D is measured from it; running without it "
            f"reports D as unknown for every cluster, so the file is required once it is named.")
    with p.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        cols = list(reader.fieldnames or [])
        for need in (barcode_column, class_column):
            if need not in cols:
                raise TaskFailure(
                    f"{p} has no {need!r} column. Its columns are: {', '.join(cols) or '(none)'}. "
                    f"Measure the structure rather than assume it - a renamed column here becomes "
                    f"a criterion silently not evaluated.")
        # The SCORE travels with the class where the file carries one. It is not a criterion -
        # the class is what the filter reads - but a deliverable that says a nucleus was called a
        # doublet and cannot say how narrowly leaves its reader unable to judge the call at all.
        score_column = str(score_column or "")
        has_score = bool(score_column) and score_column in cols
        calls, scores = {}, {}
        for row in reader:
            b = str(row[barcode_column])
            v = row[class_column]
            calls[b] = None if _unknown(v) or not str(v).strip() else str(v).strip()
            if has_score:
                try:
                    scores[b] = float(row[score_column])
                except (TypeError, ValueError):
                    scores[b] = None

    names = [str(b) for b in adata.obs_names]
    attached = sum(1 for b in names if b in calls)
    if calls and not attached:
        raise TaskFailure(
            f"none of the {len(calls):,} barcodes in {p.name} appear in this object's obs_names. "
            f"First few in the file: {', '.join(list(calls)[:3])}. First few in the object: "
            f"{', '.join(names[:3])}. The two are describing the same library under two naming "
            f"conventions; joining them anyway would report criterion D as unevaluated for every "
            f"cluster, which reads on the page like a cohort with no doublet-driven clusters.")
    adata.obs[key] = pd.Series([calls.get(b) for b in names],
                               index=adata.obs_names, dtype=object)
    if has_score:
        adata.obs[f"{key}_score"] = pd.Series([scores.get(b) for b in names],
                                              index=adata.obs_names, dtype=object)
    adata.uns["scqc_doublet_attach"] = {
        "csv": str(p), "key": key,
        "n_in_csv": len(calls), "n_in_object": len(names),
        "n_attached": attached, "n_never_scored": len(names) - attached,
        "classes_seen": sorted({v for v in calls.values() if v is not None}),
        "score_attached": bool(has_score),
        "applied": False,
    }
    return adata


def _doublet_masks(adata, doublet_key, doublet_positive):
    """(scored, positive) boolean arrays, with never-scored cells excluded from BOTH.

    A cell with no doublet score was never examined. Reading it through `or 0` counts it as a
    surviving singlet and inflates every cluster's denominator - one of the three shapes in
    docs/PRINCIPLES.md section 4. The two masks are built separately so the exclusion is visible
    rather than implied by a fillna.

    "Never examined" is decided by `_unknown`, not by `.notna()` alone: an object column can carry
    numpy.ma.masked or a numpy NaN that pandas does not always answer for, and a sentinel that
    survives into `scored` becomes a cell in criterion D's denominator. Truthiness goes through
    `_is_true` for the same reason in reverse - `v is True` is False for numpy.bool_(True), so an
    identity test would silently score every genuinely flagged cell as a singlet.
    """
    import numpy as np

    col = adata.obs[doublet_key]
    vals = np.asarray(col.to_numpy(), dtype=object).ravel()
    scored = np.asarray([not _unknown(v) for v in vals], dtype=bool)
    positive = np.zeros(vals.size, dtype=bool)
    if not _unknown(doublet_positive):
        positive[scored] = np.asarray([_is_true(v == doublet_positive) for v in vals[scored]],
                                      dtype=bool)
        return scored, positive

    # Set membership, not identity: numpy.bool_(True) hashes and compares equal to True, and
    # numpy.int64(1) to 1, so a numpy-backed column is recognised here rather than refused.
    seen = {v for v in vals[scored]}
    if seen <= {True, False, 0, 1, 0.0, 1.0}:
        positive[scored] = np.asarray([_is_true(v) for v in vals[scored]], dtype=bool)
        return scored, positive
    raise TaskFailure(
        f"obs[{doublet_key!r}] holds {sorted(map(str, seen))[:6]} and is neither boolean nor 0/1, "
        f"so which value means 'doublet' is a guess. Pass doublet_positive explicitly - guessing "
        f"it wrong inverts criterion D without changing how the table reads.")


def cluster_profile(adata, cluster_key, uninformative_genes, *, sample=None, sample_key=None,
                    doublet_key=None, doublet_positive=None, markers=True,
                    marker_topn=MARKER_TOPN, marker_method="wilcoxon", min_cells_for_markers=3,
                    mt_key="mt", ribo_key="ribo", extra_columns=True,
                    allow_empty_uninformative=False) -> list:
    """One row per (sample, cluster), carrying exactly the keys `apply_flags()` reads.

    Every row has all of `PROFILE_KEYS`. A value that could not be computed is None - never 0,
    never NaN, never False, and never pd.NA or numpy.ma.masked either, all of which read
    downstream as a value that was measured (see `_unknown`) - because `cluster_flags` is
    three-valued and its whole design depends on receiving the gap. Writing 0 for an uncomputed
    marker share turns criterion C into "evaluated and passed" for a quantity no file contains.

    `n` is every cell in the cluster. The medians (`median_umi`, `median_pct_mt`) are taken over
    the cells that have that quantity, so their denominator can be smaller than `n`; the number
    of cells missing a depth is recorded in `uns["scqc_cluster_profile"]["notes"]`. The
    alternative - letting one missing depth make the median NaN and withhold criterion A for the
    entire library - hides a one-barcode problem behind a whole sample's worth of blanks.

    Exactly one of `sample` (a literal name for a single-library object) or `sample_key` (an obs
    column, for an object holding several) must be given. There is no default and no inference
    from the file name: `umi_frac_of_sample` is a ratio against THAT sample's own median, and
    getting the grouping wrong changes criterion A everywhere at once.

    With `sample_key`, depth and mitochondrial content are per (sample, cluster) but the MARKERS
    are not: `rank_genes_groups` ranks a cluster against the rest of the object, so every sample's
    row for one cluster carries the same `pct_uninformative`. Stated here because the table does
    not show it, and a reader comparing criterion C across samples would otherwise be comparing a
    number with itself.

    `markers=False`, or an absent doublet column, produce None in those fields rather than an
    omitted key. That is the case `cluster_flags` documents as C/FLAG/WATCH being withheld: the
    flag counts must differ across a resolution sweep because of the DATA, never because of what
    was calculated.

    When `extra_columns` is on, the row also carries the denominators (`n_markers_examined`,
    `n_doublet_scored`, `n_doublet_called`) and the clustering's identity (`resolution`,
    `algorithm`, `cluster_key`) where `cluster()` recorded it. They are additive - `apply_flags()`
    copies the row and ignores what it does not read - and they exist because a rate without its
    denominator cannot be checked, and because a profile table that does not say which resolution
    produced it is unusable in a sweep.
    """
    # scanpy is imported inside the markers branch rather than here: with markers=False this
    # function is arithmetic over obs, and requiring the analysis stack to compute a table that
    # does not need it turns an available profile into an unavailable one.
    import numpy as np
    import pandas as pd

    # _unknown on both, not `is None`: a caller reading these out of a table hands over pd.NA for
    # "not set", and `pd.NA is None` is False - so both would look supplied and the object would
    # be grouped by a column literally named '<NA>'.
    sample_given, key_given = not _unknown(sample), not _unknown(sample_key)
    if sample_given == key_given:
        raise TaskFailure(
            "give exactly one of sample= (a literal library name) or sample_key= (an obs column). "
            "umi_frac_of_sample is a ratio against that sample's own median depth, so the "
            "grouping cannot be inferred.")
    if cluster_key not in adata.obs.columns:
        raise TaskFailure(f"obs[{cluster_key!r}] is absent. obs columns: "
                          f"{', '.join(map(str, adata.obs.columns))}")
    if "total_counts" not in adata.obs.columns:
        raise TaskFailure("obs['total_counts'] is absent; median_umi and umi_frac_of_sample "
                          "cannot be computed. Run qc_metrics() on the counts before clustering.")
    if key_given and sample_key not in adata.obs.columns:
        raise TaskFailure(f"obs[{sample_key!r}] is absent. obs columns: "
                          f"{', '.join(map(str, adata.obs.columns))}")

    notes = []
    obs = adata.obs
    labels = np.asarray(obs[cluster_key].astype(str))
    samples = (np.asarray(obs[sample_key].astype(str)) if key_given
               else np.asarray([str(sample)] * adata.n_obs))
    # _float_array, not asarray(dtype=float): an object column carrying pd.NA raises there and
    # numpy.ma.masked silently becomes its fill value, i.e. a cell with a made-up depth.
    totals = _float_array(obs["total_counts"])

    if "pct_counts_mt" in obs.columns:
        pct_mt = _float_array(obs["pct_counts_mt"])
    else:
        pct_mt = None
        notes.append("obs['pct_counts_mt'] absent - median_pct_mt is UNKNOWN on every row, so "
                     "criterion B is not evaluated rather than passed")

    if _unknown(doublet_key):
        scored = positive = None
        notes.append("no doublet column given - pct_doublet is UNKNOWN on every row. Criterion D "
                     "is withheld, not cleared")
    elif doublet_key not in obs.columns:
        raise TaskFailure(f"obs[{doublet_key!r}] is absent. obs columns: "
                          f"{', '.join(map(str, obs.columns))}. Pass doublet_key=None to declare "
                          f"the doublet flags genuinely unavailable.")
    else:
        scored, positive = _doublet_masks(adata, doublet_key, doublet_positive)

    all_set, mt_set, ribo_set, present, uninf_source = _resolve_uninformative(
        adata, uninformative_genes, mt_key, ribo_key, allow_empty_uninformative)
    notes.append(f"uninformative set: {uninf_source}; {len(present)} of {len(all_set)} present in "
                 f"this object")

    # ---- markers, computed once over the whole object at this clustering
    top_by_cluster = None
    if markers:
        import scanpy as sc

        sizes = pd.Series(labels).value_counts()
        too_small = sorted(sizes[sizes < int(min_cells_for_markers)].index.tolist())
        if too_small:
            raise TaskFailure(
                f"cluster(s) {', '.join(map(str, too_small[:8]))} have fewer than "
                f"{min_cells_for_markers} cells, and rank_genes_groups cannot rank a group that "
                f"small. Lower the resolution, or pass markers=False - which reports criterion C "
                f"as UNKNOWN rather than fabricating a marker share for a cluster of one.")
        mkey = f"scqc_rank_{cluster_key}"
        try:
            sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=marker_method,
                                    n_genes=int(marker_topn), key_added=mkey)
        except Exception as e:                              # noqa: BLE001 - reported, not hidden
            raise TaskFailure(
                f"rank_genes_groups(method={marker_method!r}) failed on obs[{cluster_key!r}]: "
                f"{type(e).__name__}: {e}. Criterion C depends on it; pass markers=False to "
                f"report C as UNKNOWN rather than continuing without saying so.")
        try:
            top_by_cluster = pd.DataFrame(adata.uns[mkey]["names"])
        except Exception as e:                              # noqa: BLE001 - structure, measured
            raise TaskFailure(f"could not read uns[{mkey!r}]['names'] as a table of ranked "
                              f"symbols: {type(e).__name__}: {e}")
    else:
        notes.append("markers=False - pct_uninformative, pct_mt_markers and pct_ribo_markers are "
                     "UNKNOWN on every row, and criteria C, FLAG and WATCH are withheld")

    # ---- per-sample denominators
    # Medians are taken over the barcodes that HAVE a depth. np.median of an array holding one
    # NaN is NaN, so without this a single unknown depth anywhere in the library would make that
    # library's median unknown, and with it umi_frac_of_sample on EVERY row - criterion A
    # withheld for the whole sample because of one barcode. That is the same treatment
    # pct_counts_mt already gets a few lines below, and the count excluded is recorded in `notes`
    # rather than absorbed: `n` still counts every cell in the cluster, so `n` and the number the
    # median was taken over can differ, and the note is how a reader sees that they did.
    n_total_unknown = int((~np.isfinite(totals)).sum())
    if n_total_unknown:
        notes.append(f"{n_total_unknown} of {totals.size} barcode(s) carry no usable "
                     f"total_counts (missing or non-finite); they are counted in n but excluded "
                     f"from every median depth and from the per-sample denominator")
    sample_median = {}
    for s in sorted(set(samples.tolist())):
        v = totals[samples == s]
        v = v[np.isfinite(v)]
        m = float(np.median(v)) if v.size else None
        sample_median[s] = m

    uns_cluster = adata.uns.get("scqc_cluster")
    if not isinstance(uns_cluster, dict) or uns_cluster.get("key") != cluster_key:
        uns_cluster = None

    groups = pd.DataFrame({"_s": samples, "_c": labels}).groupby(["_s", "_c"]).indices

    rows = []
    for (s, c) in sorted(groups.keys(), key=lambda k: (str(k[0]), _cluster_sort_key(k[1]))):
        idx = np.asarray(groups[(s, c)], dtype=int)
        n = int(idx.size)
        t = totals[idx]
        t = t[np.isfinite(t)]
        # None, not 0, when a cluster's every barcode lacks a depth: nothing was measured.
        med_umi = float(np.median(t)) if t.size else None

        smed = sample_median.get(str(s))
        if _unknown(med_umi) or _unknown(smed) or smed <= 0:
            # A ratio against zero or against a missing median is UNDEFINED, not large. An
            # epsilon in the denominator would print a number nobody could question.
            frac = None
        else:
            frac = float(med_umi / smed)

        if pct_mt is None:
            med_mt = None
        else:
            v = pct_mt[idx]
            v = v[np.isfinite(v)]
            med_mt = float(np.median(v)) if v.size else None

        if scored is None:
            n_scored = n_called = None
            pct_doublet = None
        else:
            n_scored = int(scored[idx].sum())
            n_called = int((scored[idx] & positive[idx]).sum())
            pct_doublet = (100.0 * n_called / n_scored) if n_scored > 0 else None

        pct_uninf = pct_mt_mk = pct_ribo_mk = None
        n_examined = None
        if top_by_cluster is not None:
            col = str(c)
            if col not in top_by_cluster.columns:
                raise TaskFailure(
                    f"cluster {col!r} has no column in the ranked-marker table (columns: "
                    f"{', '.join(map(str, top_by_cluster.columns))[:200]}). The label used for "
                    f"grouping and the label stored by rank_genes_groups have diverged.")
            # A short group is padded with NaN rather than truncated. Those are dropped and the
            # denominator shrinks with them, so the share is over what was actually examined.
            top = [str(g) for g in top_by_cluster[col].tolist()[:int(marker_topn)]
                   if not _unknown(g)]
            n_examined = len(top)
            if n_examined > 0:
                pct_uninf = 100.0 * sum(1 for g in top if g in all_set) / n_examined
                if mt_set is not None:
                    pct_mt_mk = 100.0 * sum(1 for g in top if g in mt_set) / n_examined
                if ribo_set is not None:
                    pct_ribo_mk = 100.0 * sum(1 for g in top if g in ribo_set) / n_examined

        row = {
            "sample": str(s),
            "cluster": str(c),
            "n": n,
            "median_umi": _clean(med_umi),
            "umi_frac_of_sample": _clean(frac),
            "median_pct_mt": _clean(med_mt),
            "pct_uninformative": _clean(pct_uninf),
            "pct_mt_markers": _clean(pct_mt_mk),
            "pct_ribo_markers": _clean(pct_ribo_mk),
            "pct_doublet": _clean(pct_doublet),
        }
        if extra_columns:
            row["n_markers_examined"] = n_examined
            row["n_doublet_scored"] = n_scored
            row["n_doublet_called"] = n_called
            row["cluster_key"] = str(cluster_key)
            row["algorithm"] = _clean(uns_cluster["algorithm"] if uns_cluster else None)
            row["resolution"] = _clean(uns_cluster["resolution"] if uns_cluster else None)
        rows.append(row)

    check_profile_rows(rows)
    adata.uns["scqc_cluster_profile"] = {
        "cluster_key": str(cluster_key),
        "n_rows": len(rows),
        "marker_method": (marker_method if markers else None),
        "marker_topn": (int(marker_topn) if markers else None),
        "doublet_key": doublet_key,
        "uninformative_n": len(all_set),
        "uninformative_present": len(present),
        "notes": notes,
    }
    return rows


def check_profile_rows(rows) -> None:
    """Every row carries every key `apply_flags()` reads. A missing key is not a missing value.

    `apply_flags()` calls `r["umi_frac_of_sample"]` and `r["pct_doublet"]` directly, so an absent
    key is a KeyError three steps from here, while an absent VALUE is the three-valued unknown the
    module was written to handle. The two must not be confused, so this refuses the first and
    permits the second.
    """
    for i, r in enumerate(rows):
        missing = [k for k in PROFILE_KEYS if k not in r]
        if missing:
            raise TaskFailure(f"profile row {i} is missing {', '.join(missing)}. Every row must "
                              f"carry every key in PROFILE_KEYS, with None where the value could "
                              f"not be computed.")
        for k, v in r.items():
            # Not `isinstance(v, float) and v != v`: pd.NA, pd.NaT, numpy.ma.masked and a
            # numpy.float32 NaN are none of them a float, and each of them reads downstream as a
            # value that was measured and failed - which is a pass at every gate in step 6.
            if v is not None and _unknown(v):
                raise TaskFailure(f"profile row {i} holds {v!r} ({type(v).__name__}) in {k!r}. A "
                                  f"missing-value sentinel passed as a value reads as 'measured "
                                  f"and failed the test' downstream; unknown must be None.")


# --------------------------------------------------------------------------------------------
# writing what the decision layer and the report read
# --------------------------------------------------------------------------------------------

def write_profile_csv(rows, path) -> Path:
    """The per-cluster profile as CSV. None is written as an empty cell, deliberately.

    A blank cell reads back through pandas as NaN, which `cluster_flags._unknown` catches. The
    string "nan" does not, and neither does 0.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    fields = list(PROFILE_KEYS)
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        for r in rows:
            # _unknown, not `is None`: writing pd.NA or a numpy NaN here puts the text '<NA>' or
            # 'nan' in the cell, which pandas reads back as a STRING and no unknown check catches.
            w.writerow({k: ("" if _unknown(r.get(k)) else r.get(k)) for k in fields})
    return p


def write_valleys_csv(records, path) -> Path:
    """Per-library valleys, in the columns `quality.Valley` is built from: sample, metric, valley,
    bimodal, note. A valley that could not be located is an empty cell, not a zero.

    `bimodal` may NOT be blank, and an unknown one raises rather than being written. The two
    readers of this column disagree about what an empty cell means - `scqc_cli` treats anything
    that is not false/0/no as bimodal, the tier-1 acceptance runner treats a blank as not bimodal
    - so a blank here is a claim whose meaning depends on who opens the file, and one of the two
    readings hands `quality.derive()` a threshold from a distribution nobody assessed.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "metric", "valley", "bimodal", "note"])
        for r in records:
            if _unknown(r.get("bimodal")):
                raise TaskFailure(
                    f"{r.get('sample')}/{r.get('metric')}: bimodality is {r.get('bimodal')!r}, "
                    f"neither True nor False. find_valley() always returns a bool; a record that "
                    f"reaches here without one was assembled somewhere else, and writing it "
                    f"blank would be read as 'bimodal' by scqc_cli and as 'not bimodal' by the "
                    f"acceptance runner.")
            w.writerow([r["sample"], r["metric"],
                        "" if _unknown(r.get("valley")) else r["valley"],
                        _is_true(r["bimodal"]), r.get("note", "")])
    return p


def write_density_csv(records, path) -> Path:
    """The KDE curve behind every valley, so figure F6 draws the measurement and not a redrawing.

    One row per grid point per (sample, metric): the metric's own units, the log10 grid the
    density is with respect to, the density, and whether that grid point is the valley. F6 has to
    be able to show the cut sitting where the data separates - or not sitting there.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "metric", "grid", "grid_log10", "density", "is_valley", "is_mode"])
        for r in records:
            d = r["diagnostics"]
            modes = set(d["mode_indices"])
            vi = d["valley_index"]
            has_valley = not _unknown(vi)
            for i, (g, gl, dv) in enumerate(zip(d["grid"], d["grid_log10"], d["density"])):
                # bool(), not the bare comparison: a numpy index compares to numpy.bool_, which
                # csv would write as 'True' here and as something else in another numpy version.
                w.writerow([r["sample"], r["metric"], g, gl, dv,
                            bool(has_valley and i == vi), bool(i in modes)])
    return p


# --------------------------------------------------------------------------------------------
# running one operation somewhere else
# --------------------------------------------------------------------------------------------

def build_scanpy_cmd(python_exe, op, h5ad_in, out_prefix, params_path, script=None) -> list:
    """The argv for one operation. Pure: it touches no filesystem and decides nothing.

    Separated from `run_scanpy_op` so the command can be asserted on in a test that has no
    scanpy, no cluster and no data - which is the only kind of test this half of the pipeline can
    have before it meets a real environment.

    `python_exe` is the interpreter of the ANALYSIS environment, not this one: the decision layer
    is stdlib-only by design and the analysis stack usually lives in a different prefix.
    """
    if op not in OPS:
        raise TaskFailure(f"unknown operation {op!r} (known: {', '.join(OPS)})")
    return [str(python_exe), str(script if script else Path(__file__).resolve()), str(op),
            "--h5ad", str(h5ad_in), "--out-prefix", str(out_prefix),
            "--params", str(params_path)]


def write_params_file(path, params: dict) -> Path:
    """Parameters as JSON beside the outputs. Sorted, so two identical runs produce one file."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(params, indent=2, sort_keys=True, default=str) + "\n",
                 encoding="utf-8")
    return p


def parse_metrics_json(text: str, expect_token: str = None) -> dict:
    """Parse and validate what one operation reported. Pure - takes the text, not a path.

    Refuses a partial record rather than filling the gaps: a metrics file missing its `outputs`
    list would let `run_scanpy_op` report success having checked nothing.

    `expect_token` is the per-invocation token `run_scanpy_op` put in the params file. When it is
    given, a record carrying a different token - or none - is refused, because a metrics file at
    the expected path is otherwise indistinguishable from the one an earlier run left there.
    """
    try:
        d = json.loads(text)
    except (ValueError, TypeError) as e:
        raise TaskFailure(f"the operation's metrics file is not valid JSON: {e}. The step exited "
                          f"0, so this is a bug in the writer, not a failed run.")
    if not isinstance(d, dict):
        raise TaskFailure(f"the metrics file holds {type(d).__name__}, not an object")
    for k, t in (("outputs", list), ("metrics", dict), ("versions", dict)):
        if k not in d:
            raise TaskFailure(f"the metrics file has no {k!r} key (keys: "
                              f"{', '.join(sorted(map(str, d))) or '(none)'})")
        if not isinstance(d[k], t):
            raise TaskFailure(f"metrics[{k!r}] is {type(d[k]).__name__}, expected {t.__name__}")
    if expect_token is not None:
        got = d.get("run_token")
        if _unknown(got) or str(got) != str(expect_token):
            raise TaskFailure(
                f"the metrics file carries run_token {got!r}, not {expect_token!r} - it is not "
                f"the record of this invocation. Either the command exited 0 without writing "
                f"anything and this is an earlier run's file, or the script that ran is not the "
                f"one in this repository. Numbers from a previous run reported under this run's "
                f"parameters is the failure this check exists for; it is not weakened by "
                f"deleting the file by hand.")
    return d


def _stat_signature(path) -> dict:
    """`{"size", "mtime_ns"}` for a file, as the identity of one particular written version."""
    st = Path(path).stat()
    return {"size": int(st.st_size), "mtime_ns": int(st.st_mtime_ns)}


def run_scanpy_op(op, h5ad_in, out_prefix, params, log, executor, python_exe=None,
                  script=None, timeout_s=None) -> dict:
    """Run one scanpy operation through an Executor and report only what THIS RUN wrote.

    Returns `{"outputs": [Path, ...], "metrics": {...}, "versions": {...}}`.

    EXISTENCE AFTER THE COMMAND IS NOT EVIDENCE THE COMMAND PRODUCED IT

    This function used to check that the metrics file existed once the command returned. A file
    left by an earlier run at the same prefix satisfies that check exactly as well as a fresh one,
    so a tool that exited 0 having written nothing was recorded as a success and the previous
    run's numbers were carried forward and reported under the new parameters - undetectably,
    because a stale metrics file opens and parses like a current one. Every run of a resolution
    sweep writes to the same prefix, which is where this bites.

    Three things happen here instead, and all three are required:

      BEFORE   every output the previous run DECLARED at this prefix, read out of the metrics file
               it left behind, is unlinked, and then the metrics file itself is unlinked. Nothing
               that could be mistaken for this run's work survives the start of it. Anything that
               cannot be unlinked is a refusal now, not a wrong number later.
      DURING   the params file carries a token unique to this invocation, which `main()` copies
               into the metrics record it writes.
      AFTER    the metrics file must exist, must carry THIS token, and each output it declares
               must exist with the size and mtime the writer recorded as it wrote it. The token
               proves the record is this invocation's; the size/mtime pair proves each output on
               disk is the file that invocation wrote and not a survivor beside it. Neither
               depends on comparing clocks across two machines, which is why it is done this way
               and not by timestamping the launch.

    `python_exe` defaults to the interpreter running this call, which is right only when the
    analysis stack is installed in it. On a cluster it is the analysis environment's python and
    the caller must pass it - there is no search of PATH, because silently finding a different
    scanpy produces a result the report does not describe.
    """
    import uuid

    h5ad_in = Path(h5ad_in)
    if not h5ad_in.exists():
        raise TaskFailure(f"input object does not exist: {h5ad_in}")
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    metrics_path = Path(str(out_prefix) + f".{op}.metrics.json")

    # --- BEFORE: clear this run's own landing site, including everything the last run declared.
    stale = []
    if metrics_path.exists():
        try:
            previous = json.loads(metrics_path.read_text(encoding="utf-8"))
            declared = previous.get("outputs") if isinstance(previous, dict) else None
            stale = [Path(o) for o in declared] if isinstance(declared, list) else []
        except (OSError, ValueError, TypeError):
            # An unreadable leftover tells us nothing about what it declared. The metrics file
            # still goes, so it cannot be mistaken for this run's record.
            stale = []
    removed = []
    for p in stale + [metrics_path]:
        try:
            if p.exists():
                p.unlink()
                removed.append(str(p))
        except OSError as e:
            raise TaskFailure(
                f"could not remove the previous run's output {p}: {e}. It is removed BEFORE the "
                f"command runs because a leftover that survives is indistinguishable afterwards "
                f"from a file this run wrote, and the numbers in it would be reported under this "
                f"run's parameters.")

    token = uuid.uuid4().hex
    params_path = write_params_file(str(out_prefix) + f".{op}.params.json",
                                    dict(params, **{RUN_TOKEN_KEY: token}))
    cmd = build_scanpy_cmd(python_exe if python_exe else sys.executable, op, h5ad_in,
                           out_prefix, params_path, script=script)

    # A PRIVATE NUMBA CACHE PER INVOCATION, and it is a correctness fix rather than a tuning one.
    #
    # UMAP is numba-compiled, and numba caches its compiled functions into the INSTALLED PACKAGE
    # directory by default. On a cluster that directory is one NFS path shared by every concurrent
    # job, and ten libraries embedding at once raced on it: one job died with
    # `OSError: [Errno 116] Stale file handle` inside `numba/core/caching.py` while reading the
    # cache index another job was rewriting. The traceback names numba, so it reads as an
    # environment fault rather than as a pipeline one - it is neither; it is concurrency over a
    # shared cache that nothing was serialising.
    #
    # The path is unique per (prefix, op), so no two jobs can meet in it. The cost is that numba
    # recompiles each time instead of reusing a warm cache; the alternative is a race that fails
    # one library in ten, and a failed library is a different cohort.
    cache = Path(str(out_prefix) + f".{op}.numba")
    cache.mkdir(parents=True, exist_ok=True)
    executor.shell(cmd, log=Path(log), env={"NUMBA_CACHE_DIR": str(cache)}, timeout_s=timeout_s)

    # --- AFTER: it exists, it is ours, and every file it names is the one it wrote.
    #
    # Waited for rather than checked. The op ran on another machine and these paths are on shared
    # storage; a file written a moment ago can be invisible here while the client's directory
    # cache lasts, and "it wrote nothing" is then said about a run that wrote everything.
    if await_visible([metrics_path]):
        raise TaskFailure(
            f"{op} exited 0 but wrote no metrics file at {metrics_path}, and none appeared within "
            f"{VISIBILITY_TIMEOUT_S}s (any earlier one was "
            f"removed before the run: {', '.join(removed) or 'none present'}). A step that "
            f"succeeds and produces nothing is the failure this layer refuses to pass on. "
            f"Log: {log}")
    record = parse_metrics_json(metrics_path.read_text(encoding="utf-8"), expect_token=token)

    signatures = record.get("output_stat")
    if not isinstance(signatures, dict):
        raise TaskFailure(
            f"{op}: the metrics file has no 'output_stat' map, so there is no way to tell an "
            f"output this invocation wrote from one that was already there. Log: {log}")

    outputs, missing, mismatched = [], [], []
    # One wait for the whole set rather than one per file, so a slow filesystem costs the step
    # its timeout once instead of once per output.
    invisible = set(await_visible(record["outputs"]))
    for o in record["outputs"]:
        p = Path(o)
        if str(o) in invisible or not p.exists():
            missing.append(p)
            continue
        want = signatures.get(str(o))
        if not isinstance(want, dict):
            mismatched.append(f"{p} (the run recorded no size/mtime for it)")
            continue
        got = _stat_signature(p)
        if (int(want.get("size", -1)) != got["size"]
                or int(want.get("mtime_ns", -1)) != got["mtime_ns"]):
            mismatched.append(f"{p} (on disk {got}, as written {want})")
            continue
        outputs.append(p)
    if missing:
        raise TaskFailure(f"{op} reported {len(missing)} output(s) that do not exist: "
                          f"{', '.join(map(str, missing))}. Log: {log}")
    if mismatched:
        raise TaskFailure(
            f"{op} reported {len(mismatched)} output(s) that are not the files this invocation "
            f"wrote: {'; '.join(mismatched)}. Either something else is writing into this prefix "
            f"or the file predates the run. Log: {log}")
    outputs.append(metrics_path)
    outputs.append(params_path)
    return {"outputs": outputs, "metrics": record["metrics"], "versions": record["versions"]}


# --------------------------------------------------------------------------------------------
# the script entry point - what runs inside the analysis environment
# --------------------------------------------------------------------------------------------

def _load(h5ad_in):
    import anndata as ad

    p = Path(h5ad_in)
    if not p.exists():
        raise TaskFailure(f"input object does not exist: {p}")
    return ad.read_h5ad(str(p))


def _op_qc(adata, params, out_prefix) -> tuple:
    """Compute per-cell QC and write the per-cell table the later steps read."""
    mt_prefix = _require_str(params, "mt_prefix")
    ribo_pattern = _require_str(params, "ribo_pattern")
    qc_metrics(adata, mt_prefix, ribo_pattern,
               allow_empty_mt=bool(params.get("allow_empty_mt", False)),
               allow_empty_ribo=bool(params.get("allow_empty_ribo", False)),
               layer=params.get("layer"),
               allow_transformed=bool(params.get("allow_transformed", False)))
    per_cell = Path(str(out_prefix) + ".qc_per_cell.csv.gz")
    per_cell.parent.mkdir(parents=True, exist_ok=True)
    adata.obs.to_csv(per_cell)
    outputs = [per_cell]
    if params.get("write_h5ad"):
        h5 = Path(str(out_prefix) + ".qc.h5ad")
        adata.write_h5ad(str(h5))
        outputs.append(h5)
    u = adata.uns["scqc_qc_metrics"]
    metrics = {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
               "n_mt_genes": u["n_mt_genes"], "n_ribo_genes": u["n_ribo_genes"],
               "mt_genes": u["mt_genes"], "ribo_genes": u["ribo_genes"],
               "notes": u["notes"]}
    return outputs, metrics


def _op_valley(adata, params, out_prefix) -> tuple:
    """Find the valley for each requested metric, and write the curve behind it."""
    sample = _require_str(params, "sample")
    wanted = params.get("metrics")
    if not isinstance(wanted, list) or not wanted:
        raise TaskFailure(
            f"required parameter 'metrics' must be a non-empty list naming which valleys to "
            f"measure (known: {', '.join(sorted(VALLEY_METRIC_COLUMNS))}). It is not defaulted: "
            f"step 5 derives a floor per metric, and a metric silently not measured produces a "
            f"cohort whose floor for it was never proposed rather than one that was refused.")
    unknown = [m for m in wanted if m not in VALLEY_METRIC_COLUMNS]
    if unknown:
        raise TaskFailure(f"unknown valley metric(s) {', '.join(map(str, unknown))}; known: "
                          f"{', '.join(sorted(VALLEY_METRIC_COLUMNS))}. There is no fallback "
                          f"column - a valley found on the wrong column cannot announce itself.")
    if "total_counts" not in adata.obs.columns:
        qc_metrics(adata, _require_str(params, "mt_prefix"), _require_str(params, "ribo_pattern"),
                   allow_empty_mt=bool(params.get("allow_empty_mt", False)),
                   allow_empty_ribo=bool(params.get("allow_empty_ribo", False)),
                   layer=params.get("layer"),
                   allow_transformed=bool(params.get("allow_transformed", False)))

    kw = {}
    for k in ("bw_method", "grid_size", "min_valley_depth", "min_mode_separation_log10",
              "bw_stability_factor", "min_positive", "min_mode_mass", "min_mode_observations",
              "max_points", "seed"):
        if k in params and not _unknown(params[k]):
            kw[k] = params[k]

    records = []
    for m in wanted:
        col = VALLEY_METRIC_COLUMNS[m]
        if col not in adata.obs.columns:
            raise TaskFailure(f"obs[{col!r}] is absent, so metric {m!r} cannot be measured. "
                              f"obs columns: {', '.join(map(str, adata.obs.columns))}")
        value, bimodal, diag = find_valley(adata.obs[col].to_numpy(), m, **kw)
        records.append({"sample": sample, "metric": m, "valley": value, "bimodal": bimodal,
                        "note": valley_note(diag), "diagnostics": diag})

    # The mitochondrial ceiling is derived from QUARTILES, not from a valley - that distribution
    # is unimodal, so the valley route does not apply to it (see modules/05_quality). Four numbers
    # are returned rather than the per-nucleus array, so nothing crossing the task boundary scales
    # with cell count.
    #
    # THE CEILING AND THE FLOORS ARE DERIVED ON DIFFERENT POPULATIONS, AND THAT IS THE POINT.
    #
    # It is the same pass over the same object, which is what keeps them comparable - but the
    # ceiling takes only the barcodes at or above `mito_floor_umi`, and the valleys take all of
    # them. The two quantities need opposite things:
    #
    #   a valley is a boundary BETWEEN TWO MODES, so it needs both modes present. Pre-cutting at
    #   a count floor deletes most of the debris mode, and the minimum then lands inside the
    #   nucleus mode or disappears - measured on the calibration cohort: one library's UMI valley
    #   moves from 279 to 1,018 and another's stops existing.
    #
    #   a PERCENTAGE needs a denominator big enough to mean something. A 30-UMI droplet with 10
    #   mitochondrial counts reads 33%, and a fence built on Q3 is set by exactly those. Over
    #   every called barcode the calibration cohort's ceilings come out 1.14x-1.71x higher than
    #   over the floored population, in 7 of 10 libraries, and always higher.
    #
    # This ran without the floor until 2026-08-10 and produced ceilings for a population the
    # ceiling is never applied to: every barcode it added is removed by the count floor first.
    #
    # THE SECOND RESTRICTION: `MITO_DERIVATION_MAX`, ADDED 2026-08-13.
    #
    # The floor is a statement about the DENOMINATOR. This one is a statement about what is being
    # measured, and it is the same argument the floor makes, from the other end of the axis:
    #
    #   the fence estimates WHERE THE HEALTHY POPULATION ENDS. A droplet that is more than half
    #   mitochondrial by count is not a cell whose mitochondrial content is high; on a nuclear
    #   preparation it is ambient RNA, because the mitochondrial genome is not in the nucleus and
    #   a nucleus therefore cannot be mostly mitochondrial. Those droplets are not the tail of the
    #   distribution being described - they are a different distribution, and letting them set the
    #   MAD widens the fence that decides which real cells survive.
    #
    # Measured on the calibration cohort, over its 2,462 such droplets: median 88-122 genes
    # against a population median of 339-734, median 248-351 UMI, and NOT ONE of them survived the
    # full filter. They could not influence the deliverable through any route except this one.
    #
    # THE ASYMMETRY WITH THE VALLEYS IS THE SAME ASYMMETRY, AND IT IS WHY THIS IS NOT GENERAL.
    # A derivation population may exclude what is definitionally not the thing being estimated. It
    # must never exclude what DEFINES the boundary being estimated - which is exactly why the count
    # valleys take no pre-filter at all: they are a boundary between two modes and need the debris
    # mode this would delete. Before adding a third restriction here, ask which of those two a
    # candidate is; most are the second wearing the clothes of the first.
    #
    # It is DECLARED, not derived. No statistic in the data says where a cell stops being a cell,
    # and a value chosen by sweeping until the fences looked right would be a threshold pretending
    # to be a population definition. 50% is the line at which the mitochondrial fraction exceeds
    # everything else in the droplet combined.
    mito_floor = params.get("mito_floor_umi")
    if "mito_floor_umi" not in params or _unknown(mito_floor):
        raise TaskFailure(
            "required parameter 'mito_floor_umi' was not supplied - the count floor the "
            "mitochondrial quartiles are taken above. It has no default: a ceiling derived over "
            "every called barcode is set by droplets too shallow for a percentage to mean "
            "anything, and it is looser than the one derived over the population it is applied "
            "to. Pass the light floor. It does NOT apply to the valleys, which need the debris "
            "mode this would delete.")
    mito_floor = float(mito_floor)
    mito_max = params.get("mito_derivation_max", MITO_DERIVATION_MAX)
    if _unknown(mito_max):
        mito_max = MITO_DERIVATION_MAX
    mito_max = float(mito_max)
    if not 0.0 < mito_max <= 100.0:
        raise TaskFailure(
            f"mito_derivation_max is {mito_max:g}; it must lie in (0, 100]. It is the line above "
            f"which a droplet is not a cell at all, not a quality threshold, so a value outside "
            f"the range a percentage can take is a mistake rather than a strict setting.")
    mito, mito_pop = (None, {"floor_umi": mito_floor, "derivation_max_pct": mito_max,
                             "n_above_floor": None, "n_excluded_high_mito": None,
                             "max_above_floor": None, "n_at_or_above": None,
                             "n_all_with_a_value": None})
    if "pct_counts_mt" in adata.obs.columns:
        if "total_counts" not in adata.obs.columns:
            raise TaskFailure(
                "obs['total_counts'] is absent, so the mitochondrial quartiles cannot be "
                "restricted to the barcodes above the floor, and taking them over everything is "
                "the defect this parameter exists to prevent.")
        mito, mito_pop = mito_quartiles(
            _float_array(adata.obs["pct_counts_mt"]),
            _float_array(adata.obs["total_counts"]),
            floor_umi=mito_floor, derivation_max=mito_max)

    # RIBOSOMAL CONTENT, SUMMARISED HERE BECAUSE THE MATRIX IS ALREADY OPEN. It is reported and
    # never filtered on - no threshold is derived from it. It exists in `obs` on every run and
    # was summarised only per CLUSTER, so a per-library figure could not be read at all, and the
    # confounding panel needs one: ribosomal fraction tracks the preparation, so a difference in
    # it across a confounded split is a technical difference far more readily than a biological
    # one.
    ribo = None
    if "pct_counts_ribo" in adata.obs.columns and "total_counts" in adata.obs.columns:
        ribo = ribo_summary(_float_array(adata.obs["pct_counts_ribo"]),
                            _float_array(adata.obs["total_counts"]),
                            floor_umi=mito_floor)

    # THE DENOISER'S CELL CALL, read from THIS object rather than from a file beside it.
    #
    # Step 2 compares the aligner's call with the denoiser's. It took the denoiser's from a
    # `_cell_barcodes.csv` named in the samplesheet - a SECOND artefact, from a run nothing
    # verified was this one. Two CellBender runs of the same library at different --fpr produce
    # different calls and identically-shaped files, so the comparison could be against the wrong
    # denoising with no symptom at all.
    #
    # It is also how a naming difference hides: on the calibration cohort the object carried a
    # `<sample>_` prefix and the CSV did not, so the two sets intersected in ZERO barcodes while
    # describing the same run. Reading the call from the object removes both failures, because
    # every step after step 1 then names its cells the same way.
    #
    # After denoising an empty droplet holds no counts, so the called cells are exactly the
    # barcodes with counts remaining. That is a property of the output rather than an inference
    # about it.
    import csv as _csv
    _tot = adata.X.sum(axis=1)
    _tot = _tot.A1 if hasattr(_tot, "A1") else _tot.ravel()
    called = [str(b) for b, keep in zip(adata.obs_names, _tot > 0) if keep]
    cells_path = Path(str(out_prefix) + ".called_barcodes.csv")
    with cells_path.open("w", encoding="utf-8", newline="") as fh:
        _w = _csv.writer(fh)
        _w.writerow(["barcode"])
        for _b in called:
            _w.writerow([_b])

    # ------------------------------------------------------------------ the nuclear fraction
    #
    # MEASURED HERE AND APPLIED NOWHERE IN THIS STEP. It is computed in the pass that already
    # holds the object and its barcodes, so the 131 MB per-droplet file is parsed ONCE per
    # library, in a task that already runs per library and in parallel. Step 7 then reads the
    # small per-barcode CSV written below rather than the aligner's file again.
    #
    # Absent `cellreads_stats` nothing here runs, nothing is written, and the metrics carry no
    # nuclear fraction - so a run that did not ask for it is unchanged.
    nf_path = None
    nf_metrics = {}
    if not _unknown(params.get("cellreads_stats")):
        from adapters import nuclear_fraction as _nf

        nf_by_bc, nf_stats = _nf.read_cellreads(
            params["cellreads_stats"], list(adata.obs_names), sample=sample,
            antisense=bool(params.get("nf_antisense", False)))
        nf_path = Path(str(out_prefix) + ".nuclear_fraction.csv")
        with nf_path.open("w", encoding="utf-8", newline="") as fh:
            _w = _csv.writer(fh)
            _w.writerow(["barcode", "nuclear_fraction"])
            for _b in adata.obs_names:
                _v = nf_by_bc.get(str(_b))
                # BLANK, not 0, for an undefined or unjoined barcode. A zero here is a claim that
                # the droplet held no intronic signal, which is a measurement nobody made.
                _w.writerow([str(_b), "" if _v is None else f"{_v:.6f}"])

        # THE MEDIAN IS TAKEN OVER THE POPULATION THE CRITERION CAN ACT ON: barcodes at or above
        # the light floor with a DEFINED fraction. Below that floor the barcodes are debris the
        # count floor removes anyway, and including them drags the median toward whatever the
        # ambient looks like in this library rather than toward what its nuclei look like. It is
        # the same population restriction the mitochondrial quartiles already use, for the same
        # reason, and it is recorded beside the number.
        _tot_all = adata.obs["total_counts"] if "total_counts" in adata.obs.columns else None
        _floor = float(params.get("mito_floor_umi") or 0)
        _pop = []
        for _b, _c in zip(adata.obs_names, (list(_tot_all) if _tot_all is not None
                                            else [None] * adata.n_obs)):
            _v = nf_by_bc.get(str(_b))
            if _v is None:
                continue
            if _c is not None and float(_c) < _floor:
                continue
            _pop.append(_v)
        nf_metrics = {
            "nf_median": _nf.median(_pop),
            "nf_n_in_median": len(_pop),
            "nf_median_floor_umi": _floor,
            **{f"nf_{k}": v for k, v in nf_stats.items()},
        }

    valleys = write_valleys_csv(records, str(out_prefix) + ".valleys.csv")
    density = write_density_csv(records, str(out_prefix) + ".valley_density.csv")
    jpath = Path(str(out_prefix) + ".valleys.json")
    jpath.write_text(json.dumps(records, indent=2, default=str) + "\n", encoding="utf-8")
    metrics = {"sample": sample,
               "valleys": {r["metric"]: r["valley"] for r in records},
               "bimodal": {r["metric"]: r["bimodal"] for r in records},
               "notes": {r["metric"]: r["note"] for r in records},
               # None when obs carries no pct_counts_mt, or the library is too small to place a
               # quartile. Absent is reported as absent; step 5 refuses rather than defaulting.
               "mito_quartiles": mito,
               "mito_population": mito_pop,
               # None when obs carries no pct_counts_ribo - which is every run made before this
               # was summarised per library. Absent is reported as absent, never as zero.
               "ribo_quartiles": ribo}
    metrics["n_called_by_denoiser"] = len(called)
    metrics.update(nf_metrics)
    outs = [valleys, density, jpath, cells_path]
    if nf_path is not None:
        outs.append(nf_path)
    return outs, metrics


def mito_quartiles(pct_mt, total_counts, *, floor_umi, derivation_max=MITO_DERIVATION_MAX):
    """`(quartiles, population)` for the mitochondrial ceiling. Selects; removes nothing.

    # rule-one: no-removal - this reads two per-barcode arrays and returns summary statistics.

    Lifted out of `_op_valley` so that the population the ceiling is derived over can be tested
    without an AnnData and a KDE. It is the whole of the derivation-population definition, in one
    place: a barcode contributes when it has both values, sits at or above `floor_umi`, and lies
    strictly below `derivation_max`. The reasoning for each restriction, and the measurements
    behind them, are in `_op_valley` beside the call.

    `quartiles` is None when fewer than four barcodes qualify - too few to place a quartile, which
    step 5 refuses on rather than defaulting. `population` is returned either way, so a reader can
    never see a ceiling without seeing what it was taken over, and can never see its absence
    without seeing why.

    Both are returned from ONE pass so they cannot disagree. Computing the population record
    separately is how `n_at_or_above` comes to describe a different set from the one the quartiles
    were placed on, and nothing downstream could detect it.
    """
    above = [float(m) for m, c in zip(pct_mt, total_counts)
             if m == m and c == c and float(c) >= floor_umi]
    v = sorted(m for m in above if m < derivation_max)
    # The true observed maximum, taken BEFORE the derivation cut. Keeping it is what makes the cut
    # safe: a run reporting only the derivation population's own maximum would say the worst
    # droplet in the library was 49.9% however bad it really was, and the one number that shows a
    # library is full of ambient would be the number the restriction erased. Never derived from.
    max_above = max(above) if above else None
    pop = {"floor_umi": float(floor_umi),
           "derivation_max_pct": float(derivation_max),
           "n_above_floor": len(above),
           "n_excluded_high_mito": len(above) - len(v),
           "max_above_floor": max_above,
           "n_at_or_above": len(v) if len(v) >= 4 else None,
           "n_all_with_a_value": int(sum(1 for x in pct_mt if x == x))}
    if len(v) < 4:
        return None, pop

    def q(p):
        h = (len(v) - 1) * p
        lo_i = int(h)
        hi_i = min(lo_i + 1, len(v) - 1)
        return v[lo_i] + (h - lo_i) * (v[hi_i] - v[lo_i])

    # MAD travels with the quartiles because both are needed and the matrix is open exactly once.
    # The ceiling is derived from median + k*1.4826*MAD; Tukey's Q3 + 1.5*IQR is carried alongside
    # as an independent second derivation, and the two disagreeing is a finding rather than
    # something for one of them to absorb. Computing it here rather than in the decision layer
    # keeps the per-nucleus array on this side of the task boundary - nothing that scales with
    # cell count crosses.
    med = q(0.5)
    dev = sorted(abs(x - med) for x in v)
    h = (len(dev) - 1) * 0.5
    lo_i = int(h)
    hi_i = min(lo_i + 1, len(dev) - 1)
    return ({"n": len(v), "median": med, "q1": q(0.25), "q3": q(0.75),
             "mad": dev[lo_i] + (h - lo_i) * (dev[hi_i] - dev[lo_i]),
             "max": v[-1], "max_above_floor": max_above}, pop)


def ribo_summary(pct_ribo, total_counts, *, floor_umi):
    """Median and quartiles of `pct_counts_ribo`, over the barcodes above the UMI floor.

    DELIBERATELY NOT `mito_quartiles`. That function derives a CEILING - it excludes high-mito
    barcodes so the derivation population is cells rather than debris, and carries a MAD because
    a threshold is computed from it. Nothing here decides a threshold: ribosomal content is
    reported, never filtered on, and reusing the mito machinery would import an exclusion rule
    that has no meaning for this metric.

    The UMI floor is shared, and only that: a summary taken over every droplet describes the
    ambient soup as much as the cells.

    Returns None when the population is too small to place quartiles in - the same rule the mito
    derivation uses - rather than a number computed from three points.
    """
    v = sorted(float(r) for r, t in zip(pct_ribo, total_counts)
               if r == r and t == t and float(t) >= float(floor_umi))
    if len(v) < 4:
        return None

    def q(pr):
        h = (len(v) - 1) * pr
        lo_i = int(h)
        hi_i = min(lo_i + 1, len(v) - 1)
        return v[lo_i] + (h - lo_i) * (v[hi_i] - v[lo_i])

    return {"n": len(v), "median": q(0.5), "q1": q(0.25), "q3": q(0.75), "max": v[-1]}


#: The declared population's threshold criteria, as (declared key, obs column, comparison).
#: The cell call is handled separately because it is a flag rather than a threshold.
POPULATION_CRITERIA = (("umi_floor", "total_counts", "ge"),
                       ("gene_floor", "n_genes_by_counts", "ge"),
                       ("mito_ceiling", "pct_counts_mt", "le"))


def _population_mask(adata, pop, *, honour=("cell_call_key", "umi_floor", "gene_floor",
                                            "mito_ceiling")) -> tuple:
    """`(keep, applied)` for a declared population - a boolean mask and what produced it.

    `honour` names which of the declared criteria are applied, and it is the whole reason this is
    a function rather than a block inside `_op_cluster`: the clustering takes all four, the
    embedding takes the cell call alone, and both read the SAME declaration. Two copies of this
    arithmetic would eventually disagree about what a cell call is, and the figure and the flags
    would then describe populations that differ by an amount nothing measures.
    """
    import numpy as _np

    keep = _np.ones(int(adata.n_obs), dtype=bool)
    applied: list = []
    cc = pop.get("cell_call_key")
    if "cell_call_key" in honour and not _unknown(cc):
        if cc not in adata.obs.columns:
            raise TaskFailure(
                f"population declares cell_call_key={cc!r} and obs has no such column. "
                f"Clustering the unfiltered object instead would silently reproduce the "
                f"defect this parameter exists to fix.")
        keep &= _float_array(adata.obs[cc].astype(float)) > 0.5
        applied.append(f"{cc}")
    for key, col, op in POPULATION_CRITERIA:
        if key not in honour:
            continue
        v = pop.get(key)
        if _unknown(v):
            continue
        if col not in adata.obs.columns:
            raise TaskFailure(
                f"population declares {key}={v} and obs has no {col!r} to apply it to.")
        arr = _float_array(adata.obs[col])
        # An unmeasured value is NOT a pass. NaN fails both comparisons in numpy, which is the
        # behaviour wanted here: a cell whose depth was never measured is not evidence of a
        # cell that passed.
        keep &= (arr >= float(v)) if op == "ge" else (arr <= float(v))
        applied.append(f"{col}{'>=' if op == 'ge' else '<='}{v}")
    return keep, applied


def _op_cluster(adata, params, out_prefix) -> tuple:
    """Cluster at one resolution and profile the result."""
    if "resolution" not in params or _unknown(params["resolution"]):
        raise TaskFailure("required parameter 'resolution' was not supplied")
    if "seed" not in params or _unknown(params["seed"]):
        raise TaskFailure("required parameter 'seed' was not supplied; a clustering with an "
                          "unrecorded seed cannot be reproduced or compared")
    # An omitted doublet column and a deliberately absent one are the same JSON. They are not the
    # same fact, so the key must be PRESENT - null is how the caller declares "no doublet flags",
    # which reports criterion D as unknown rather than clear.
    if "doublet_key" not in params:
        raise TaskFailure(
            "parameter 'doublet_key' must be present, and may be null. Omitting it and declaring "
            "no doublet flags produce the same profile, and step 6's criterion D is the one that "
            "becomes a tautology if it is evaluated at the wrong point in the pipeline.")

    # ORDER MATTERS AND GETTING IT WRONG IS SILENT (see this module's header). `cluster()`
    # normalises adata.X in place, after which every cell's counts sum to the same number and
    # `median_umi` - criterion A's whole basis - is a constant that still prints plausibly. So
    # cluster() refuses to start without obs["total_counts"], and this is where it comes from:
    # the object this step reads is the denoised one, which nothing has measured yet.
    #
    # Measured HERE rather than carried from step 5, because step 5 writes no object. Same
    # patterns, so the mitochondrial percentage criterion B thresholds and the one step 5's
    # ceiling came from are the same quantity.
    if "total_counts" not in adata.obs.columns:
        qc_metrics(adata, _require_str(params, "mt_prefix"), _require_str(params, "ribo_pattern"),
                   allow_empty_mt=bool(params.get("allow_empty_mt", False)),
                   allow_empty_ribo=bool(params.get("allow_empty_ribo", False)),
                   layer=params.get("layer"),
                   allow_transformed=bool(params.get("allow_transformed", False)))

    # ATTACHED, not applied. Nothing is removed here and nothing downstream of this op reads the
    # column except criterion D.
    if not _unknown(params.get("doublet_csv")):
        if _unknown(params.get("doublet_key")):
            raise TaskFailure(
                "doublet_csv was supplied but doublet_key is null, so the calls would be written "
                "to a column nothing reads and criterion D would report unknown with the data for "
                "it sitting in the object.")
        attach_doublet_calls(
            adata, params["doublet_csv"], key=str(params["doublet_key"]),
            barcode_column=str(params.get("doublet_barcode_column", "barcode")),
            class_column=str(params.get("doublet_class_column", "doublet_class")))

    # THE POPULATION STEP 6 IS SPECIFIED TO CLUSTER.
    #
    # `cluster_flags.py` requires "the step-5 object: quality-filtered, with the doublet flags
    # ATTACHED and NOT applied". The object this op reads is the DENOISED FULL DROPLET MATRIX, so
    # without this mask the clustering runs over empty droplets. Measured on the calibration
    # cohort, that cost: 1,398 of 1,531 clusters held no cell that reached the deliverable; leiden
    # spends a fixed resolution over what it is given, so the real cells were left UNDER-RESOLVED
    # at 245 clusters for a library with 8 real ones; and because every criterion is a cluster
    # MEDIAN, the libraries with the most droplet noise raised the FEWEST flags - backwards, and
    # on that cohort the bias ran along a design factor.
    #
    # THE MASK IS THE QUALITY FILTER WITHOUT THE DOUBLET CRITERION. Doublets stay attached and
    # unapplied because criterion D is only computable before a removal: afterwards every cluster
    # is 0% doublet by construction, which is the tautology the doublet_key check above guards.
    #
    # NOTHING IS REMOVED FROM A DELIVERABLE HERE. This op writes cluster labels for step 7 to
    # carry; step 7 remains the only place an observation is dropped. `population` is required to
    # be PRESENT and may be null - an omitted key and a deliberate "cluster everything" are the
    # same JSON and are not the same fact.
    if "population" not in params:
        raise TaskFailure(
            "parameter 'population' must be present, and may be null. Omitting it and declaring "
            "an unfiltered population produce the same clustering, and which cells were clustered "
            "is what decides whether step 6's flags describe nuclei or empty droplets.")
    # An omitted embedding and a declared "do not embed" are the same JSON and are not the same
    # fact, for the reason every required-present-may-be-null parameter here exists: a run that
    # produced no coordinates because nobody asked and one that produced none because the request
    # was lost look identical in the report, and only the second is a defect.
    if "embedding" not in params:
        raise TaskFailure(
            "parameter 'embedding' must be present, and may be null. It declares the population "
            "whose 2-D coordinates are written out for figures F10 and F11; null means no "
            "embedding is computed and the report says the coordinates were never asked for.")
    pop = params.get("population")
    n_before = int(adata.n_obs)

    # --- THE EMBEDDING, COMPUTED FIRST AND OVER A WIDER POPULATION THAN THE CLUSTERING.
    #
    # First, because `cluster()` normalises adata.X IN PLACE and the embedding needs counts.
    #
    # Wider, because of what the two are for. The clustering must run over the cells that reach
    # the deliverable, or its flags describe empty droplets. The embedding must contain the cells
    # a criterion REMOVED, or F11 - "did the removed nuclei leave as a population, or scattered?"
    # - is drawn over a population every one of them has already left, and it can only ever answer
    # "there were none". So the embedding takes the cell call and stops there: every nucleus the
    # count floors, the ceiling or the doublet call removed is still in the picture, marked.
    embedding = params.get("embedding")
    emb_barcodes, emb_coords, emb_record = None, None, None
    if not _unknown(embedding):
        if not isinstance(embedding, dict):
            raise TaskFailure(
                f"parameter 'embedding' is a {type(embedding).__name__}; it must be null or a "
                f"mapping declaring at least a population. `true` would have to mean a default "
                f"population chosen here, and which population is embedded is exactly what "
                f"decides whether F11 can see the nuclei a criterion removed.")
        emb_pop = embedding.get("population", "cell_called")
        if emb_pop not in ("cell_called", "clustered"):
            raise TaskFailure(
                f"embedding declares population={emb_pop!r}; it must be 'cell_called' (every "
                f"barcode the denoiser called, so the removed nuclei are in the picture) or "
                f"'clustered' (exactly what step 6 clusters, which contains no nucleus the "
                f"floors or the ceiling removed and cannot show where they sat).")
        honour = (("cell_call_key",) if emb_pop == "cell_called"
                  else ("cell_call_key", "umi_floor", "gene_floor", "mito_ceiling"))
        emb_keep, emb_applied = ((_population_mask(adata, pop, honour=honour))
                                 if not _unknown(pop) else (None, []))
        emb_adata = adata[emb_keep].copy() if emb_keep is not None else adata.copy()
        if int(emb_adata.n_obs) < 50:
            raise TaskFailure(
                f"the embedding population leaves {emb_adata.n_obs} cell(s) of {n_before} - too "
                f"few to embed. A projection of that many points is not a manifold.")
        emb_barcodes = [str(b) for b in emb_adata.obs_names]
        emb_coords, emb_record = embed(
            emb_adata, seed=int(params["seed"]),
            n_hvg=int(embedding.get("n_hvg", params.get("n_hvg", 2000))),
            n_pcs=int(embedding.get("n_pcs", params.get("n_pcs", 50))),
            n_neighbors=int(embedding.get("n_neighbors", params.get("n_neighbors", 15))),
            min_dist=float(embedding.get("min_dist", 0.5)),
            allow_transformed=bool(params.get("allow_transformed", False)))
        emb_record["population"] = emb_pop
        emb_record["population_applied"] = emb_applied or ["nothing - every droplet embedded"]
        emb_record["n_before"] = n_before
        # Freed before the clustering copy is taken. Two copies of a droplet matrix alive at once
        # is the difference between this op fitting in its memory request and being killed by the
        # scheduler, and the kill arrives as an exit code with no message about why.
        del emb_adata

    pop_note = "UNFILTERED - every droplet clustered, including empty ones"
    if not _unknown(pop):
        keep, applied = _population_mask(adata, pop)
        n_keep = int(keep.sum())
        if n_keep < 50:
            raise TaskFailure(
                f"the declared population leaves {n_keep} cell(s) of {n_before} - too few to "
                f"cluster. A threshold this severe is a mistake in the declaration, and "
                f"clustering the remainder would produce flags about nothing.")
        adata = adata[keep].copy()
        pop_note = (f"quality-filtered before clustering: {n_keep} of {n_before} "
                    f"({100.0 * n_keep / n_before:.1f}%) by {', '.join(applied)}; "
                    f"doublets ATTACHED and NOT applied")

    uninformative = params.get("uninformative_genes")
    if _unknown(uninformative):
        # The locked set is the mt+ribo symbols qc_metrics MATCHED in this object, taken as a
        # mapping so criterion C's split into C_mt and C_ribo is authoritative rather than
        # re-derived from var flags. It is resolved here rather than passed in because the
        # orchestrator never loads a matrix (tests/test_wiring.py check E) and so cannot know
        # which symbols a pattern matched against THIS reference - and a set locked against a
        # different one reports 0% for every cluster instead of failing.
        q = adata.uns.get("scqc_qc_metrics") or {}
        if "mt_genes" in q or "ribo_genes" in q:
            uninformative = {"mt": list(q.get("mt_genes") or []),
                             "ribo": list(q.get("ribo_genes") or [])}
    if _unknown(uninformative):
        raise TaskFailure(
            "the locked mt+ribo set criterion C is measured against could not be established. "
            "Supply 'uninformative_genes', or supply 'mt_prefix' and 'ribo_pattern' on an object "
            "that has not already been measured, so the symbols are matched against this "
            "reference. Without it C would read 0% for every cluster.")
    # The remaining .get() defaults below mirror the documented defaults of cluster() and
    # cluster_profile(); each is a FIXED procedure parameter, not a stand-in for a declaration.
    cluster(adata,
            resolution=float(params["resolution"]),
            seed=int(params["seed"]),
            n_hvg=int(params.get("n_hvg", 2000)),
            n_pcs=int(params.get("n_pcs", 50)),
            n_neighbors=int(params.get("n_neighbors", 15)),
            key_added=params.get("key_added"),
            allow_transformed=bool(params.get("allow_transformed", False)))
    key = adata.uns["scqc_cluster"]["key"]
    rows = cluster_profile(adata, key, uninformative,
                           sample=params.get("sample"),
                           sample_key=params.get("sample_key"),
                           doublet_key=params.get("doublet_key"),
                           doublet_positive=params.get("doublet_positive"),
                           markers=bool(params.get("markers", True)),
                           marker_topn=int(params.get("marker_topn", MARKER_TOPN)),
                           marker_method=str(params.get("marker_method", "wilcoxon")))
    profile = write_profile_csv(rows, str(out_prefix) + ".cluster_profile.csv")

    # WHICH BARCODE IS IN WHICH CLUSTER, written out rather than discarded.
    #
    # The profile is one row per (sample, cluster) and cannot be joined back to a barcode, so
    # without this the per-cell assignment leiden just computed dies with the process - and the
    # next stage re-clusters to ask a question this step has already answered. It is a two-column
    # table rather than a second copy of the matrix: the labels are what is needed downstream, and
    # the matrix is already on disk.
    labels = Path(str(out_prefix) + ".cluster_labels.csv")
    with labels.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["barcode", "sample", "cluster"])
        col = adata.obs[key]
        for b, c in zip(adata.obs_names, col):
            w.writerow([str(b), params.get("sample") or "", "" if _unknown(c) else str(c)])

    outputs = [profile, labels]

    # THE COORDINATES, WITH THE CLUSTER LABEL JOINED ON WHERE THERE IS ONE.
    #
    # Written after the clustering rather than beside the embedding, so the file carries both:
    # every embedded barcode, its position, whether it was one of the cells step 6 clustered, and
    # which cluster if so. A reader can then put F8's flagged clusters on F10's manifold without
    # this op having to draw anything.
    #
    # `clustered` is FALSE, not blank, for a cell the mask dropped. It was embedded and it was
    # not clustered; both are facts about it, and a blank would read as "not recorded".
    if emb_coords is not None:
        clustered_at = {str(b): ("" if _unknown(c) else str(c))
                        for b, c in zip(adata.obs_names, adata.obs[key])}
        emb_path = Path(str(out_prefix) + ".embedding.csv")
        with emb_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["barcode", "sample", "x", "y", "clustered", "cluster"])
            for b, (x, y) in zip(emb_barcodes, emb_coords):
                lab = clustered_at.get(b)
                w.writerow([b, params.get("sample") or "", f"{x:.6g}", f"{y:.6g}",
                            lab is not None, "" if lab is None else lab])
        outputs.append(emb_path)
        emb_record["n_also_clustered"] = sum(1 for b in emb_barcodes if b in clustered_at)

    # --- THE EXTRA RESOLUTIONS, PROFILED LAST AND INTO SIBLING FILES.
    #
    # LAST, because everything the deliverable depends on is already on disk by this point. An
    # extra resolution is a diagnostic; it must not be able to cost a run its cluster flags, so a
    # failure here is CAUGHT, named and counted rather than raised.
    #
    # `cluster()` is NOT called again - it normalises adata.X in place, and a second call would
    # normalise an already-normalised matrix (and be refused for it). Only leiden is re-run,
    # against the neighbour graph already built, which is resolution-independent. That also makes
    # the extras honest: they differ from the default in resolution and in nothing else.
    #
    # The file name carries the resolution and therefore does NOT end in `.cluster_profile.csv`,
    # which is what step 6's flag task selects on. Nothing downstream can pick one up by accident.
    extra_failed, extra_done = {}, []
    extras = params.get("extra_resolutions")
    if not _unknown(extras):
        import scanpy as sc

        primary_uns = dict(adata.uns["scqc_cluster"])
        for r in sorted({float(x) for x in extras}):
            if r == float(params["resolution"]):
                continue                      # the default is not also an extra
            ekey = f"leiden_res{r:g}"
            try:
                sc.tl.leiden(adata, resolution=r, random_state=int(params["seed"]),
                             key_added=ekey)
                # cluster_profile reads uns["scqc_cluster"] for the clustering's identity and
                # drops those columns when the key does not match. Swapped per resolution and
                # restored below, so the primary record is exactly what it was.
                adata.uns["scqc_cluster"] = {
                    **primary_uns, "key": ekey, "resolution": r,
                    "n_clusters": int(adata.obs[ekey].astype(str).nunique())}
                erows = cluster_profile(adata, ekey, uninformative,
                                        sample=params.get("sample"),
                                        sample_key=params.get("sample_key"),
                                        doublet_key=params.get("doublet_key"),
                                        doublet_positive=params.get("doublet_positive"),
                                        markers=bool(params.get("markers", True)),
                                        marker_topn=int(params.get("marker_topn", MARKER_TOPN)),
                                        marker_method=str(params.get("marker_method",
                                                                     "wilcoxon")))
                ep = write_profile_csv(erows, f"{out_prefix}.cluster_profile.res{r:g}.csv")
                el = Path(f"{out_prefix}.cluster_labels.res{r:g}.csv")
                with el.open("w", newline="", encoding="utf-8") as fh:
                    w = csv.writer(fh)
                    w.writerow(["barcode", "sample", "cluster"])
                    for b, c in zip(adata.obs_names, adata.obs[ekey]):
                        w.writerow([str(b), params.get("sample") or "",
                                    "" if _unknown(c) else str(c)])
                outputs.extend([ep, el])
                extra_done.append({"resolution": r, "n_clusters": len(erows)})
            except Exception as e:                          # noqa: BLE001 - named, not hidden
                # The commonest cause is a cluster too small for rank_genes_groups, which a finer
                # resolution reaches sooner. Recorded as a failure of THAT resolution; the run
                # and the default's flags are unaffected.
                extra_failed[f"{r:g}"] = f"{type(e).__name__}: {e}"
        adata.uns["scqc_cluster"] = primary_uns

    if params.get("write_h5ad"):
        h5 = Path(str(out_prefix) + ".clustered.h5ad")
        adata.write_h5ad(str(h5))
        outputs.append(h5)
    u = adata.uns["scqc_cluster"]
    metrics = {"key": key, "resolution": u["resolution"], "seed": u["seed"],
               "n_clusters": u["n_clusters"], "n_hvg_selected": u["n_hvg_selected"],
               "pca_variance_ratio_sum": u["pca_variance_ratio_sum"],
               "leiden_flavor_default": u["leiden_flavor_default"],
               "flagged_hvg": u["flagged_hvg"],
               "profile_notes": adata.uns["scqc_cluster_profile"]["notes"],
               "n_profile_rows": len(rows),
               # Present and null when no embedding was asked for, so the report can tell "not
               # requested" from "requested and it produced nothing".
               "embedding": emb_record,
               # Both, always. "which extras were profiled" cannot be read off the successes
               # alone: a resolution that was asked for and failed must not look like one that
               # was never asked for.
               "extra_resolutions": (extra_done if not _unknown(extras) else None),
               "extra_resolutions_failed": (extra_failed or None)}
    return outputs, metrics


#: The applied criteria, in the order they are written to the per-cell table. The names are the
#: column names, not a convention the reader has to reconstruct: `audit_removal.audit()` is told
#: which columns are criteria and checks that `removed` decomposes into exactly them.
APPLY_CRITERIA = ("fail_not_cellbender_cell", "fail_umi_floor", "fail_gene_floor",
                  "fail_mito_ceiling", "fail_doublet", "fail_mito_nf")


def _op_apply_measure(adata, params, out_prefix) -> tuple:
    """Every applied criterion, per barcode, for ONE library. Removes NOTHING and writes no object.

    This is the measuring half of step 7 and it exists separately for one reason: the removal
    gate has to run before anything is materialised. An op that filtered and wrote in one pass
    would have performed the removal before the operator's approval was checked, and an approval
    checked afterwards is a record, not a gate.

    So this writes a table covering EVERY barcode in the object - kept and removed alike - with
    one boolean per criterion. `build_removal_record()` turns that into the ledger, the gate
    checks the ledger against the mask, and only then does `apply_write` touch a matrix.

    Every criterion is TRUE or FALSE for every barcode, never unknown: `build_removal_record()`
    refuses an unknown, because an observation removed on a criterion nobody evaluated has not
    been judged. Where that cannot honestly be done the op refuses instead - see the light-floor
    check below.
    """
    import csv as _csv

    import numpy as np

    sample = _require_str(params, "sample")
    for k in ("umi_floor", "gene_floor", "mito_ceiling_pct", "light_floor"):
        if k not in params or _unknown(params[k]):
            raise TaskFailure(
                f"required parameter {k!r} was not supplied. Step 7 applies it; there is no "
                f"default for a threshold that removes observations.")
    if "doublet_csv" not in params:
        raise TaskFailure(
            "parameter 'doublet_csv' must be present, and may be null. Omitting it and declaring "
            "no doublet calls produce the same table, and the difference is whether the doublet "
            "criterion was evaluated or merely absent.")
    umi_floor = float(params["umi_floor"])
    gene_floor = float(params["gene_floor"])
    ceiling = float(params["mito_ceiling_pct"])
    light_floor = float(params["light_floor"])

    if "total_counts" not in adata.obs.columns:
        qc_metrics(adata, _require_str(params, "mt_prefix"), _require_str(params, "ribo_pattern"),
                   allow_empty_mt=bool(params.get("allow_empty_mt", False)),
                   allow_empty_ribo=bool(params.get("allow_empty_ribo", False)),
                   layer=params.get("layer"),
                   allow_transformed=bool(params.get("allow_transformed", False)))

    scored = np.zeros(adata.n_obs, dtype=bool)
    is_doublet = np.zeros(adata.n_obs, dtype=bool)
    dbl_class = [None] * adata.n_obs
    dbl_score = [None] * adata.n_obs
    if not _unknown(params.get("doublet_csv")):
        key = str(params.get("doublet_key") or "doublet_class")
        attach_doublet_calls(
            adata, params["doublet_csv"], key=key,
            barcode_column=str(params.get("doublet_barcode_column", "barcode")),
            class_column=str(params.get("doublet_class_column", "doublet_class")))
        scored, is_doublet = _doublet_masks(adata, key, params.get("doublet_positive"))
        dbl_class = list(adata.obs[key])
        if f"{key}_score" in adata.obs.columns:
            dbl_score = list(adata.obs[f"{key}_score"])

    counts = _float_array(adata.obs["total_counts"])
    genes = _float_array(adata.obs["n_genes_by_counts"])
    mt = _float_array(adata.obs["pct_counts_mt"])

    # The object holds what the denoiser produced; a barcode left with no counts is not a cell it
    # called. Recorded as its own criterion rather than folded into the UMI floor, or the ledger
    # would say a droplet was removed for being shallow when it was never a cell.
    is_cell = np.asarray([c == c and c > 0 for c in counts], dtype=bool)
    fail_cell = ~is_cell
    fail_umi = np.asarray([not (c == c and c >= umi_floor) for c in counts], dtype=bool)
    fail_gene = np.asarray([not (g == g and g >= gene_floor) for g in genes], dtype=bool)
    # A barcode with no mitochondrial value fails NOTHING here. It is already removed by the cell
    # criterion (no counts means no percentage), and inventing a failure for it would put a
    # criterion in the ledger that was never evaluated.
    fail_mito = np.asarray([m == m and m > ceiling for m in mt], dtype=bool)
    fail_doublet = np.asarray(is_doublet, dtype=bool)

    # ------------------------------------- the joint mitochondrial x nuclear-fraction criterion
    #
    #     fail_mito_nf = (mt > TRIGGER) AND (nf < FLOOR)
    #
    # ADDITIVE, never a replacement. `fail_mito_ceiling` is untouched and removes exactly what it
    # removed with this off, so the removal is a strict SUPERSET of the run without it and every
    # additional barcode comes from the band between the trigger and the applied ceiling.
    #
    # THE TRIGGER REMOVES NOTHING. It selects which barcodes the nuclear fraction is consulted on.
    # Below it the fraction is never read, because a droplet whose mitochondrial percentage is
    # already inside the declared bound is not one this criterion has anything to say about.
    #
    # All three parameters must be PRESENT and may be null, the same contract the other optional
    # parameters here carry: a run with no joint criterion is a supported state and is recorded as
    # one, while a missing parameter is a wiring defect - and the two produce identical tables.
    for _k in ("nf_csv", "nf_floor", "nf_trigger_pct"):
        if _k not in params:
            raise TaskFailure(
                f"parameter {_k!r} must be present, and may be null. A run with no joint "
                f"mitochondrial x nuclear-fraction criterion is a supported state; a missing "
                f"parameter is a wiring defect, and both produce the same table.")
    nf_values = [None] * adata.n_obs
    joint_armed = not any(_unknown(params[k]) for k in ("nf_csv", "nf_floor", "nf_trigger_pct"))
    fail_mito_nf = np.zeros(adata.n_obs, dtype=bool)
    if joint_armed:
        nf_floor = float(params["nf_floor"])
        nf_trigger = float(params["nf_trigger_pct"])
        if not 0 < nf_floor < 1:
            raise TaskFailure(
                f"{sample}: nf_floor is {nf_floor:g}; it must lie strictly in (0, 1). The nuclear "
                f"fraction is intronic / (intronic + exonic), a RATIO and not a percentage.")
        by_bc = {}
        with open(params["nf_csv"], encoding="utf-8", newline="") as fh:
            for row in _csv.DictReader(fh):
                raw = (row.get("nuclear_fraction") or "").strip()
                by_bc[row["barcode"]] = float(raw) if raw else None
        nf_values = [by_bc.get(str(b)) for b in adata.obs_names]
        n_missing = sum(1 for b in adata.obs_names if str(b) not in by_bc)
        if n_missing:
            raise TaskFailure(
                f"{sample}: {n_missing:,} of {adata.n_obs:,} barcodes are absent from "
                f"{params['nf_csv']}. The table is written per barcode of this same object in "
                f"step 5, so an absent row means the two describe different objects - which is "
                f"not a filter that can be applied to some of a library.")
        # AN UNKNOWN NUCLEAR FRACTION FAILS NOTHING, exactly as an unknown mitochondrial value
        # does. It is not evidence of an intact nucleus and it is not evidence of debris, so it
        # cannot be allowed to decide either way: a barcode over the trigger whose second axis is
        # unmeasured is KEPT, and it is counted separately from the ones the fraction spared.
        fail_mito_nf = np.asarray(
            [(m == m and m > nf_trigger) and (v is not None and float(v) < nf_floor)
             for m, v in zip(mt, nf_values)], dtype=bool)

    fails = {"fail_not_cellbender_cell": fail_cell, "fail_umi_floor": fail_umi,
             "fail_gene_floor": fail_gene, "fail_mito_ceiling": fail_mito,
             "fail_doublet": fail_doublet, "fail_mito_nf": fail_mito_nf}
    removed = np.zeros(adata.n_obs, dtype=bool)
    for v in fails.values():
        removed |= v
    keep = ~removed

    # A KEPT BARCODE THAT WAS NEVER DOUBLET-SCORED IS A HOLE IN THE FILTER, and it is silent:
    # `fail_doublet` reads False for it, which in the ledger means "this criterion did not remove
    # it" and on the page is indistinguishable from "it was examined and found to be a singlet".
    # It cannot arise while the UMI floor sits above the light floor, which is the usual case -
    # so it would appear only when someone lowered one of them, which is exactly when nobody is
    # looking at the interaction between the two.
    if not _unknown(params.get("doublet_csv")):
        holes = int(np.sum(keep & ~scored))
        if holes:
            raise TaskFailure(
                f"{sample}: {holes:,} barcodes would be KEPT that were never scored for "
                f"doublets. The doublet detector saw only barcodes at or above the light floor "
                f"({light_floor:g} UMI) and the applied UMI floor is {umi_floor:g}, so the gap "
                f"between them reaches the deliverable unexamined. In the ledger those barcodes "
                f"read as having passed the doublet criterion. Raise the UMI floor to the light "
                f"floor or above, or score the detector over the lower population.")

    # None when this object carries no pct_counts_ribo - an object written before the metric
    # existed, or a reference whose symbols matched none. The column is then blank on every row,
    # which is what "not measured" looks like and what 0 would not.
    ribo = (_float_array(adata.obs["pct_counts_ribo"])
            if "pct_counts_ribo" in adata.obs.columns else None)

    path = Path(str(out_prefix) + ".percell.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh)
        # `nuclear_fraction` sits beside `pct_counts_mt` because the two are the criterion's two
        # axes and a reader checking `fail_mito_nf` needs both. BLANK where undefined - never 0.
        # `pct_counts_ribo` sits beside `pct_counts_mt` because they are read together: both are
        # properties of the PREPARATION rather than of the condition under test, which is what
        # makes them readable across a design whose factors cannot be separated. It is reported
        # and never filtered on - no criterion below uses it. BLANK where absent, never 0: a
        # library whose reference matched no ribosomal gene must not read as one with none.
        w.writerow(["barcode", "sample", "cellbender_cell", "total_counts", "n_genes",
                    "pct_counts_mt", "pct_counts_ribo", "nuclear_fraction",
                    "doublet_score", "doublet_class", "doublet_scored",
                    *APPLY_CRITERIA, "removed", "keep"])
        for i, b in enumerate(adata.obs_names):
            w.writerow([str(b), sample, bool(is_cell[i]),
                        "" if counts[i] != counts[i] else f"{counts[i]:.0f}",
                        "" if genes[i] != genes[i] else f"{genes[i]:.0f}",
                        "" if mt[i] != mt[i] else f"{mt[i]:.6f}",
                        "" if ribo is None or ribo[i] != ribo[i] else f"{ribo[i]:.6f}",
                        "" if nf_values[i] is None else f"{float(nf_values[i]):.6f}",
                        "" if _unknown(dbl_score[i]) else f"{float(dbl_score[i]):.6f}",
                        "" if _unknown(dbl_class[i]) else str(dbl_class[i]),
                        bool(scored[i]),
                        *[bool(fails[c][i]) for c in APPLY_CRITERIA],
                        bool(removed[i]), bool(keep[i])])

    metrics = {"sample": sample, "n_in": int(adata.n_obs), "n_keep": int(keep.sum()),
               "n_removed": int(removed.sum()),
               "thresholds": {"umi_floor": umi_floor, "gene_floor": gene_floor,
                              "mito_ceiling_pct": ceiling, "light_floor": light_floor,
                              "nf_floor": (None if not joint_armed else float(params["nf_floor"])),
                              "nf_trigger_pct": (None if not joint_armed
                                                 else float(params["nf_trigger_pct"]))},
               "joint_armed": bool(joint_armed),
               # What the second axis SPARED and what it could not judge - the two numbers the
               # criterion cannot be read without. `n_mito_over_trigger` is its denominator.
               "n_mito_over_trigger": (0 if not joint_armed else int(sum(
                   1 for m in mt if m == m and m > float(params["nf_trigger_pct"])))),
               "n_nf_unknown_over_trigger": (0 if not joint_armed else int(sum(
                   1 for m, v in zip(mt, nf_values)
                   if m == m and m > float(params["nf_trigger_pct"]) and v is None))),
               "n_doublet_scored": int(scored.sum()),
               **{f"n_{c}": int(fails[c].sum()) for c in APPLY_CRITERIA}}
    return [path], metrics


#: Annotation columns that are LABELS, never quantities. A cluster id reads as a number and is
#: not one: coercing "10" to 10.0 turns an identifier into an arithmetic type, and anything that
#: then sorts or averages it produces a result with no meaning and no error.
_ANNOTATION_LABELS = ("cluster", "sample")


def _annotation_value(column: str, raw):
    """One annotation cell, with its type back and its unknowns intact.

    Read from a CSV, so everything arrives as text. An EMPTY cell is unknown and becomes None -
    never False and never 0.0. `x or None` was the first version of this and it is wrong in a way
    worth naming: a genuine 0.0 is falsy, so a cluster measured at 0% doublet came out as a
    cluster nobody measured.
    """
    if raw is None:
        return None
    text = str(raw).strip()
    if text == "":
        return None
    if column in _ANNOTATION_LABELS:
        return text
    if text == "True":
        return True
    if text == "False":
        return False
    try:
        return float(text)
    except ValueError:
        return text


def _unique_var_index(a, sample: str):
    """Give `var` a unique index, keeping the symbols, or refuse. Returns a note or None.

    Gene SYMBOLS are not unique in most references - this cohort's carries duplicates - and
    an object indexed by them cannot be merged: anndata.concat reindexes the alternate axis
    and pandas raises "cannot reindex on an axis with duplicate labels". It is also a hazard
    in the delivered object itself, where `adata[:, "Gm5773"]` silently returns two columns.

    So when, and only when, the index is not unique, it is replaced by the identifier axis -
    which is unique, and which the ambient audit already had to fall back on for the same
    reason - and the symbols are kept in `var["gene_symbol"]`. A reference whose symbols are
    unique is left exactly as it came in, so this changes nothing for most inputs.

    Matching on symbol - `mt-`, `^Rp[sl](?!6k)` - therefore reads var["gene_symbol"] on a delivered
    object whose index was replaced. Every threshold in this pipeline was computed before this
    point, on the per-library objects, so no measurement is affected by the choice.
    """
    import pandas as pd

    names = list(map(str, a.var_names))
    if len(set(names)) == len(names):
        return None
    for col in ID_COLUMNS:
        if col not in a.var.columns:
            continue
        ids = [str(v) for v in a.var[col]]
        if len(set(ids)) != len(ids):
            continue
        if "gene_symbol" not in a.var.columns:
            a.var["gene_symbol"] = pd.Categorical(names)
        a.var.index = pd.Index(ids, name=col)
        dup = len(names) - len(set(names))
        return (f"var was indexed by gene symbol, which is not unique ({dup:,} duplicate "
                f"label(s) over {len(names):,} genes); re-indexed by {col} and the symbols "
                f"kept in var['gene_symbol']")
    from collections import Counter
    worst = [g for g, _ in Counter(names).most_common(3)]
    raise TaskFailure(
        f"{sample}: var is indexed by a label that is not unique "
        f"({len(names) - len(set(names)):,} duplicates over {len(names):,} genes, e.g. "
        f"{', '.join(worst)}) and the object carries no unique identifier column to use "
        f"instead (looked for {', '.join(ID_COLUMNS)}). The libraries cannot be merged on "
        f"this axis, and making the labels unique by appending suffixes would invent gene "
        f"symbols that are not in the reference.")


def _annotation_column(values: list):
    """One annotation column, typed from its own values, with unknown still unknown.

    Every column used to be assigned as `dtype=object`, which is not writable: anndata
    encodes an object column as a STRING array, so writing the deliverable died on the
    first numeric one - `total_counts` - with "Can't implicitly convert non-string objects
    to strings", and the three-valued flags and any None-bearing label column were queued
    up behind it to fail the same way.

    The dtype is inferred from the values rather than from a list of column names, because
    the columns are whatever the annotation table carries and a hardcoded list silently
    mistypes any column added later. Each choice keeps missing distinct from present:

        bool   -> pandas nullable boolean   unknown is pd.NA, which is not False
        number -> float64                   unknown is NaN, which is not 0.0
        other  -> categorical               unknown is NaN, which is not ""

    All three round-trip through h5ad with the unknowns intact, so "this cluster was not
    flagged" and "this barcode was never examined" stay different facts in the delivered
    object and not only in the table it was built from.
    """
    import numpy as np
    import pandas as pd

    known = [v for v in values if v is not None]
    if known and all(isinstance(v, bool) for v in known):
        return pd.array([None if v is None else bool(v) for v in values], dtype="boolean")
    if known and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in known):
        return pd.array([np.nan if v is None else float(v) for v in values], dtype="float64")
    return pd.Categorical([None if v is None else str(v) for v in values])


def _op_apply_write(adata, params, out_prefix) -> tuple:
    """Materialise the deliverable: the kept barcodes of every library, in ONE object.

    Reached only after the gate has passed, which is why it takes a keep-list rather than a
    threshold - re-deriving the mask here would mean the object was filtered by something other
    than what was approved, and the two could differ without either being wrong on its own terms.
    """
    import anndata
    import numpy as np
    import pandas as pd

    # Absolute, like the imports at the top of this file and for the same reason: this module is
    # RUN AS A SCRIPT by the orchestrator, so it has no parent package and `from . import` raises
    # `attempted relative import with no known parent package` - at step 7, after every other
    # step has completed, which is the most expensive place in the pipeline to fail.
    from adapters import declaration as decl
    from adapters import matrix as mx

    libs = params.get("libraries")
    if not isinstance(libs, list) or not libs:
        raise TaskFailure("required parameter 'libraries' must be a non-empty list of "
                          "{sample, h5ad, keep_csv} entries.")
    # Run identity, from the orchestrator. Absent when this op is driven directly, in which case
    # the fields are written empty rather than invented - a declaration that names a run key the
    # object did not come from is worse than one that names none.
    prov = params.get("provenance") or {}
    out_h5 = Path(str(out_prefix) + ".deliverable.h5ad")
    # BOTH SHAPES, because they answer different questions and the per-library subsets exist in
    # memory a moment before the concatenation anyway. The combined object is what integration and
    # differential testing read; the per-library ones are what annotation reads, and this
    # pipeline's cluster check is per library precisely because identity is decided there.
    per_sample_dir = Path(str(out_prefix) + "_per_sample")
    per_sample_dir.mkdir(parents=True, exist_ok=True)
    per_lib, var_ref, var_from = [], None, None
    written_per_sample, var_notes = [], []
    for entry in libs:
        s = str(entry.get("sample"))
        keep_path = Path(str(entry.get("keep_csv")))
        wanted = [ln.strip() for ln in
                  keep_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
        a = _load(str(entry.get("h5ad")))
        note = _unique_var_index(a, s)
        if note and not var_notes:
            var_notes.append(note)
        # The var axis is not re-derived, it is CHECKED. anndata.concat with an outer join fills
        # a gene absent from one library with zeros, which is fabricated data wearing the shape
        # of a measurement; an inner join silently narrows the feature space instead.
        if var_ref is None:
            var_ref, var_from = list(map(str, a.var_names)), s
        elif list(map(str, a.var_names)) != var_ref:
            raise TaskFailure(
                f"{s} and {var_from} do not carry the same genes in the same order "
                f"({a.n_vars:,} against {len(var_ref):,}). Concatenating them would either "
                f"invent zeros for the genes one of them lacks or quietly drop them. Re-quantify "
                f"against one reference.")
        present = set(map(str, a.obs_names))
        missing = [b for b in wanted if b not in present]
        if missing:
            raise TaskFailure(
                f"{s}: {len(missing):,} approved barcodes are not in {entry.get('h5ad')} "
                f"(first: {', '.join(missing[:3])}). The keep-list and the object have come "
                f"apart; filtering anyway would deliver a different population from the one "
                f"that was approved.")
        sub = a[np.asarray([str(b) in set(wanted) for b in a.obs_names], dtype=bool)].copy()
        sub.obs["sample"] = pd.Categorical([s] * sub.n_obs)

        # WHAT STEP 6 FOUND, CARRIED ONTO THE NUCLEI IT FOUND IT ABOUT. Step 6 clusters each
        # library and flags each cluster, and until now none of that reached the deliverable: the
        # per-cell assignment was computed and discarded, so the next stage had to re-cluster to
        # ask a question already answered. A barcode with no annotation is left None, never a
        # blank string or a False - "this cluster was not flagged" and "this barcode was not in
        # the table" are different facts.
        ann_csv = entry.get("annotations_csv")
        if ann_csv:
            ann_path = Path(str(ann_csv))
            if not ann_path.exists():
                raise TaskFailure(
                    f"{s}: annotations_csv {ann_path} does not exist. It is the only route by "
                    f"which step 6's result reaches the object; writing the deliverable without "
                    f"it would silently produce one that has never been cluster-checked.")
            with ann_path.open(encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh)
                cols = [c for c in (reader.fieldnames or []) if c != "barcode"]
                table = {str(r["barcode"]): r for r in reader}
            for c in cols:
                sub.obs[c] = _annotation_column(
                    [_annotation_value(c, table.get(str(b), {}).get(c))
                     for b in sub.obs_names])

        per_lib.append({"sample": s, "n_in": int(a.n_obs), "n_kept": int(sub.n_obs)})
        # WHAT THE FLAG MEANS, on the object that carries it. Step 6's verdict reaches the
        # deliverable as `cluster_FLAG` and nothing on the object said what it was for, so a
        # downstream tool had to either ignore it - annotating nuclei this pipeline flagged - or
        # guess from the column name, which is that tool deciding what is technical on scQC's
        # behalf. The declaration is stamped on the PER-LIBRARY objects too, not only the merged
        # one, because the per-library objects are what annotation reads.
        decl.stamp(sub, sample=s, run_key=prov.get("run_key", ""),
                   commit=prov.get("commit", ""), version=prov.get("version", ""))
        one = per_sample_dir / f"{s}.filtered.h5ad"
        # Classic string encoding: a nullable-string `obs/_index` is unreadable by everything
        # except anndata, and the deliverable exists to be read by other things.
        with mx.classic_string_encoding():
            sub.write_h5ad(str(one))
        written_per_sample.append(one)
        del a, sub

    # THE COMBINED OBJECT IS THE MERGE OF THE PER-LIBRARY ONES, and is built by reading them back
    # rather than from the copies that were in memory when they were written. The two would
    # ordinarily be identical, which is the point: if they are not - a failed write, a truncated
    # file, an encoding that did not round-trip - the difference belongs in the combined object
    # rather than hidden by having kept a good copy in memory. What a reader can open is what was
    # merged.
    #
    # It also means the per-library objects are the primary artifact and the combined one is
    # derived, which is the honest ordering: each library is filtered on its own thresholds and
    # cluster-checked on its own clustering, and the cohort object is those results put together.
    parts = [_load(str(p)) for p in written_per_sample]
    combined = parts[0] if len(parts) == 1 else anndata.concat(
        parts, axis=0, join="inner", merge="same", label=None, index_unique=None)
    if combined.obs_names.has_duplicates:
        raise TaskFailure(
            "the combined object has duplicate barcodes, so a row in the removal ledger cannot "
            "be matched back to one observation. Prefix each library's barcodes before this "
            "point.")
    # COLUMNS, not a list of dicts. HDF5 has no record type for the latter and h5py fails on it
    # with "Can't implicitly convert non-string objects to strings", at the end of the only step
    # that removes anything - after the gate has passed and the filtering is done. Three parallel
    # arrays write natively and read back without a parser.
    combined.uns["scqc_apply"] = {
        "library": [str(x["sample"]) for x in per_lib],
        "library_n_in": [int(x["n_in"]) for x in per_lib],
        "library_n_kept": [int(x["n_kept"]) for x in per_lib],
        "n_delivered": int(combined.n_obs),
        "built_from": [str(p.name) for p in written_per_sample],
        "note": "the merge of the per-library filtered objects, read back from disk; the ledger "
                "names every barcode removed and the criteria that removed it",
        "var_index": var_notes or ["var is indexed as it arrived; gene labels are unique"],
    }
    # The merge must account for every nucleus the libraries kept. A concatenation that silently
    # dropped one - a var axis that did not align, a barcode colliding across libraries - would
    # produce a smaller cohort that still looks complete.
    expected = sum(int(x["n_kept"]) for x in per_lib)
    if int(combined.n_obs) != expected:
        raise TaskFailure(
            f"the merged object holds {combined.n_obs:,} nuclei and the per-library objects hold "
            f"{expected:,} between them. The merge lost or duplicated observations, so the cohort "
            f"object is not the sum of its libraries and no count taken from it describes them.")
    # Computed over the MERGED flag rather than carried from a part: the digest has to describe
    # the column in the object it is stamped on, or verifying it downstream proves nothing about
    # the file in hand. `sample` is empty here because this object is every library at once.
    decl.stamp(combined, sample="", run_key=prov.get("run_key", ""),
               commit=prov.get("commit", ""), version=prov.get("version", ""))
    with mx.classic_string_encoding():
        combined.write_h5ad(str(out_h5))
    return ([out_h5] + written_per_sample,
            {"n_delivered": int(combined.n_obs), "n_genes": int(combined.n_vars),
             "libraries": per_lib,
             "per_sample_objects": [str(p) for p in written_per_sample],
             "obs_columns": sorted(map(str, combined.obs.columns)),
             "var_index": str(combined.var.index.name or "gene label"),
             "var_index_note": var_notes[0] if var_notes else ""})


def main(argv=None) -> int:
    """Run one operation in this process. Invoked by `run_scanpy_op` under the analysis python.

    The metrics file is written LAST and only after every output has been stat'd here as well, so
    its existence is itself the evidence that the operation finished. A step that exits 0 and
    wrote nothing is caught on both sides of the executor boundary.

    Two things go into that file for the caller to check, and neither is optional:

      run_token     copied from the params file. It is what makes the caller's "the metrics file
                    is there" into "the metrics file is MINE" - a leftover from an earlier run at
                    the same prefix carries an earlier token and is refused.
      output_stat   the size and mtime of every output AS THIS PROCESS FINISHED WRITING IT, so
                    the caller can tell a file this run produced from one that was already there
                    under the same name. Two processes' clocks are never compared: the writer
                    reads its own stat and the reader compares against that.
    """
    ap = argparse.ArgumentParser(description="scQC scanpy operations (run by the adapter)")
    ap.add_argument("op", choices=OPS)
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--params", required=True)
    a = ap.parse_args(argv)

    params = json.loads(Path(a.params).read_text(encoding="utf-8"))
    if not isinstance(params, dict):
        raise TaskFailure(f"{a.params} does not hold a JSON object")

    adata = None if a.op in OPS_WITHOUT_INPUT else _load(a.h5ad)
    handler = {"qc": _op_qc, "valley": _op_valley, "cluster": _op_cluster,
               "apply_measure": _op_apply_measure, "apply_write": _op_apply_write}[a.op]
    outputs, metrics = handler(adata, params, a.out_prefix)

    missing = [str(p) for p in outputs if not Path(p).exists()]
    if missing:
        raise TaskFailure(f"{a.op} finished without writing: {', '.join(missing)}")

    record = {"op": a.op, "outputs": [str(p) for p in outputs], "metrics": metrics,
              "run_token": params.get(RUN_TOKEN_KEY),
              "output_stat": {str(p): _stat_signature(p) for p in outputs},
              "versions": observed_versions()}
    Path(str(a.out_prefix) + f".{a.op}.metrics.json").write_text(
        json.dumps(record, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"{a.op}: wrote {len(outputs)} output(s)")
    for p in outputs:
        print(f"  {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
