# Draws the nine report figures from PLAIN DATA. Nothing here reads a matrix, opens an AnnData,
# recomputes a statistic or decides a threshold: every number that appears on a panel arrives as
# an argument, so a figure cannot state something the pipeline did not compute.
"""The nine figures of docs/REPORT_DESIGN.md, one function each.

WHY THESE FUNCTIONS TAKE PLAIN LISTS AND DICTS

A figure that computes its own summary is a second implementation of a number the report already
has, and the two can disagree without anything noticing - the page then carries a value no file on
disk contains. So each function takes exactly what it draws. Where a quantity was not supplied it
is drawn as NOT SUPPLIED in the never-examined colour, never as zero and never omitted: an absent
denominator that silently becomes 1, or an absent count that silently becomes 0, is the failure
docs/PRINCIPLES.md section 4 exists to prevent, and on a figure it is invisible.

THE RULES EVERY FUNCTION HERE OBEYS, AND THE READING EACH ONE PREVENTS

  no invented legend     Every annotation is an argument. A legend entry stating a rate the
                         caller did not pass cannot be checked against anything.
  never-examined         Its own colour AND its own hatch, never merged into a negative. A
                         barcode nobody scored is not a barcode that passed.
  denominators           Every rate carries the population it is a rate OF, in the axis label.
                         5.70% of all cells and 7.52% of scored cells are both true and only one
                         of them is comparable with a published band.
  one scale              Panels meant for comparison share limits. A per-panel scale makes a 4x
                         difference look identical to a 4% one.
  n on every panel       A panel with no n reads as universal whether or not that was meant.
  thresholds are drawn   A cut at 350 is a line at 350, not a sentence under the figure.
  labels via label_text  An optional label is read with the unknown predicate, never with
                         `label or default`: a NaN is truthy, so `or` kept it and the axis read
                         literally `nan (log scale, original units)`.
  shape as well as hue   Every flag carries a marker or a hatch, so the figure survives being
                         printed in grey and being read by someone who cannot separate the hues.
  log axes in units      Ticks read 1,000, not 10^3.

WHAT THESE FIGURES CANNOT DO. None of them establishes that a threshold is correct, that a
removal was justified, or that a cluster is technical. They show where a cut sits against the
distribution it was derived from, and how evenly a loss falls across the design. Interpretation
is the reader's, and the report says so beside each one.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

# The engine package sits beside this one. Loading this module by path - which the report writer
# and the tests both do - leaves the repository root off sys.path, so it is put there explicitly
# rather than letting the import fail with a message about a package nobody was asked to install.
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from engine.task import TaskFailure  # noqa: E402

DEFAULT_DPI = 150

#: Okabe-Ito, which stays separable in the common forms of colour vision deficiency and in
#: greyscale print. `unknown` is reserved: no other category may use it.
PALETTE = {
    "ok": "#0072B2",
    "review": "#E69F00",
    "refuse": "#D55E00",
    "unknown": "#CC79A7",
    "neutral": "#56B4E9",
    "second": "#009E73",
    "threshold": "#000000",
    "muted": "#666666",
    # BEFORE and AFTER are a pair, and the pairing is the message: the prior state is present but
    # recessive, the retained state carries the ink. Deliberately not two hues - a reader
    # comparing two distributions of one quantity should not have to decode a colour legend to
    # know which is the result. It also survives greyscale, where two hues of equal value do not.
    "before": "#BFBFBF",
    "after": "#2F353C",
}

#: Never-examined carries a hatch as well as a colour, because colour is not the only encoding.
UNKNOWN_HATCH = "///"

NOT_SUPPLIED = "NOT SUPPLIED"


# ------------------------------------------------------------------------------------ helpers
# Everything above the first figure is pure and importable without matplotlib, so the argument
# construction and the parsing of a caller's data structures are testable on their own.


#: Modules consulted for their missing-value singletons, cached by name so a module that is not
#: installed is searched for once rather than on every value. A cached module object stays valid,
#: and a cached `None` stays valid too: a package that cannot be imported cannot be the source of
#: a sentinel later in the same process.
_OPTIONAL_MODULES: dict = {}


def _optional_module(name: str):
    """An already-imported module, or one imported now, or None. Never raises.

    `sys.modules` is consulted first because it answers without paying for an import - and it is
    sufficient on its own for the sentinels below, since an instance of `pandas.NA` cannot exist
    in a process that has not imported pandas. The real import is the fallback, so the predicate
    is still correct if this module is somehow handed a sentinel by a caller that imported it
    under another name, and it degrades to None rather than raising where the package is absent.
    """
    if name in _OPTIONAL_MODULES:
        return _OPTIONAL_MODULES[name]
    mod = sys.modules.get(name)
    if mod is None:
        try:
            import importlib
            mod = importlib.import_module(name)
        except Exception:                                                 # noqa: BLE001
            mod = None
    _OPTIONAL_MODULES[name] = mod
    return mod


def _missing_singletons() -> tuple:
    """`pandas.NA`, `pandas.NaT` and `numpy.ma.masked`, as far as they can be reached now.

    Compared by IDENTITY at the call site, never by equality: `pd.NA == pd.NA` is `pd.NA`, and
    truth-testing that raises TypeError, so an equality test would turn a missing value into a
    crash halfway through drawing a panel.
    """
    out = []
    pd = _optional_module("pandas")
    if pd is not None:
        out += [s for s in (getattr(pd, "NA", None), getattr(pd, "NaT", None)) if s is not None]
    np = _optional_module("numpy")
    if np is not None:
        masked = getattr(getattr(np, "ma", None), "masked", None)
        if masked is not None:
            out.append(masked)
    return tuple(out)


def _unknown(v) -> bool:
    """True when a value carries no information. THE unknown predicate for this module.

    Real tables deliver more spellings of "no value" than None and a Python NaN, and every one
    that slips through is then compared with `>=` or handed to `int()`, where it reads as a value
    that FAILED the test - which on every figure here is drawn as a pass. Covered:

      None                      the caller said nothing
      float('nan')              a blank cell from csv/float parsing
      numpy.float64('nan')      a numpy-backed column; it IS a float subclass, so `isinstance`
                                catches it, but `is None` never would
      numpy.float32('nan') etc  not a float subclass - caught by the self-inequality test below
      numpy.datetime64('NaT')   same
      pandas.NA, pandas.NaT     neither is None and neither is a float; identity only
      numpy.ma.masked           a masked element read out of a masked array

    Order is deliberate: the two fast paths (None, float) answer for almost every value without
    touching sys.modules, so this stays cheap when it is called once per plotted point.

    An object whose `!=` cannot be reduced to one boolean - a numpy ARRAY, most obviously - is
    NOT unknown. It is a container, and calling it unknown would silently drop a whole series.
    """
    if v is None:
        return True
    if isinstance(v, float):                     # Python float, and numpy.float64 (a subclass)
        # `bool(...)`, because `numpy.float64('nan') != itself` is `numpy.bool_(True)`, and a
        # predicate that returns a numpy bool has re-created the problem it exists to solve one
        # call up the stack.
        return bool(v != v)
    if isinstance(v, (int, str, bytes, bytearray)):     # bool is an int; none of these is NaN
        return False
    for sentinel in _missing_singletons():
        if v is sentinel:
            return True
    try:
        return bool(v != v)
    except Exception:                                                     # noqa: BLE001
        return False


def _tri(v):
    """Three-valued reading of a flag: True, False, or None for `never determined`.

    `v is True` is False for `numpy.bool_(True)`, so a genuinely flagged row arriving from any
    numpy-backed table read as not-flagged - the figure then shows the opposite of the finding.
    The unknown check comes first and the truth test second, so a sentinel can never be read as
    a determination, and a determination is read whatever type carries it.

    A value whose truth is ambiguous (an array) returns None: it was not a determination this
    figure can read, and NOT EVALUATED is the honest drawing of that.
    """
    if _unknown(v):
        return None
    if isinstance(v, str):
        # A string is read from its VOCABULARY, never by truthiness. `bool("False")` is True, so
        # a step-6 profile round-tripped through CSV - which is how a profile usually reaches a
        # figure - drew EVERY cluster as flagged, including the ones explicitly recorded False.
        # A word this function does not recognise is unknown rather than True: guessing at it is
        # how that defect arose.
        s = v.strip().lower()
        if s in ("true", "t", "yes", "y", "1"):
            return True
        if s in ("false", "f", "no", "n", "0"):
            return False
        return None
    try:
        return bool(v)
    except Exception:                                                     # noqa: BLE001
        return None


def _mapping(value) -> dict:
    """An optional mapping argument, normalised to `{}` when nothing was supplied.

    `thresholds or {}` was the old spelling and a NaN defeats it: `not float('nan')` is False, so
    a sentinel passed where a dict was expected sailed past the default and raised AttributeError
    on the first `.get` - a crash instead of the NOT SUPPLIED the figure is designed to draw.
    """
    if _unknown(value) or not value:
        return {}
    return value


def _median(values):
    """The median of the known values, or None if there are none.

    `sorted(v)[n // 2]` is the UPPER order statistic whenever n is even, not the median. It was
    labelled "median" on F2 for every even cohort - and a legend must not state a number the code
    did not compute, so the number is now the one the word means.
    """
    vals = sorted(float(v) for v in values if not _unknown(v))
    if not vals:
        return None
    mid = len(vals) // 2
    return vals[mid] if len(vals) % 2 else 0.5 * (vals[mid - 1] + vals[mid])


def n_label(n, what: str = "n") -> str:
    """`n = 1,234`, or a visible statement that it was not supplied.

    There is no third form. A panel whose n was not passed says so on its face rather than
    printing nothing, because a figure with no n reads as universal.

    A value that is present but cannot be read as a whole number is stated as such rather than
    being allowed to raise: it is a defect in what the caller passed, and a defect that kills the
    figure hides the eight panels that were fine along with the one that was not.
    """
    if _unknown(n):
        return f"{what}: {NOT_SUPPLIED}"
    try:
        return f"{what} = {int(n):,}"
    except (TypeError, ValueError, OverflowError):
        return f"{what}: NOT NUMERIC ({n!r})"


def denominator_label(text) -> str:
    """`denominator: <population>`, or a visible statement that it was not supplied."""
    if _unknown(text) or not str(text).strip():
        return f"denominator: {NOT_SUPPLIED}"
    return f"denominator: {text}"


def label_text(value, default: str) -> str:
    """An axis label, a title or a point's name from an optional argument. THE label reader.

    Every label on every figure below is built through this. `value or default` is not an
    equivalent spelling of it, and the difference is visible on the page:

      float('nan') or 'metric'   is float('nan'), because a NaN is TRUTHY. F7's y axis then read
                                 literally `nan (log scale, original units)` - a caption that
                                 names the metric `nan` and reads, to anyone who did not draw it,
                                 like a metric somebody called that.
      numpy.float32('nan'),
      numpy.datetime64('NaT'),
      pandas.NaT                 the same, printing `nan` or `NaT` into the axis label.
      pandas.NA                  never even reached the axis: `pd.NA or default` raises TypeError
                                 out of the middle of drawing the panel, so one blank cell in a
                                 metric name destroyed the whole figure.

    So the unknown predicate decides, not truthiness. An unknown label falls back to the default;
    so does a string that is empty or only whitespace, because a blank axis is indistinguishable
    from one nobody labelled. Anything else is rendered exactly as the caller wrote it - including
    `"0"` and `"False"`, which `or` would have silently replaced with the default.
    """
    if _unknown(value):
        return default
    text = str(value)
    return text if text.strip() else default


def grid_shape(n: int) -> tuple:
    """Rows and columns for `n` panels, wide rather than tall so the shared axis stays readable."""
    if n <= 0:
        raise TaskFailure("a figure was asked for with no libraries to draw; pass at least one")
    cols = 1 if n == 1 else (2 if n <= 4 else (3 if n <= 9 else 4))
    rows = -(-n // cols)
    return rows, cols


def _series(sample, entry, xkey: str, ykey: str) -> tuple:
    """Pull an (x, y) pair out of one library's entry, refusing anything that cannot be drawn."""
    if not isinstance(entry, dict):
        raise TaskFailure(f"library {sample!r}: expected a dict with {xkey!r} and {ykey!r}, "
                          f"got {type(entry).__name__}")
    x, y = entry.get(xkey), entry.get(ykey)
    # `_unknown`, not `is None`: a series handed over as a blank cell (NaN, pandas.NA) is an
    # absent series, and letting it through would fail later inside `list()` with a message that
    # names a type rather than the library and the key that were not supplied.
    if _unknown(x) or _unknown(y):
        raise TaskFailure(f"library {sample!r}: {xkey!r} and {ykey!r} are both required and "
                          f"{'both are' if _unknown(x) and _unknown(y) else 'one is'} absent. A "
                          f"partial series cannot be drawn, and drawing the half that exists "
                          f"would understate the curve without saying so")
    x, y = list(x), list(y)
    if len(x) != len(y):
        raise TaskFailure(f"library {sample!r}: {len(x):,} {xkey} values against {len(y):,} "
                          f"{ykey} values. A series whose axes disagree in length is plotted at "
                          f"the wrong positions and nothing downstream can detect that it was")
    if not x:
        raise TaskFailure(f"library {sample!r}: the series is empty")
    return x, y


def _positive_pairs(x, y) -> tuple:
    """The pairs a log axis can show, and how many it cannot.

    Points at or below zero are dropped by a log axis silently. They are counted and reported on
    the panel instead, because a curve that lost its left-hand end without saying so reads as a
    curve that never had one.
    """
    keep_x, keep_y, dropped = [], [], 0
    for a, b in zip(x, y):
        if _unknown(a) or _unknown(b) or a <= 0 or b <= 0:
            dropped += 1
        else:
            keep_x.append(float(a))
            keep_y.append(float(b))
    return keep_x, keep_y, dropped


def criterion_contributions(rows, criteria=None) -> dict:
    """Unique and shared removals per criterion, from a removal record's rows.

    `rows` is the `(identifier, [criterion, ...])` sequence `modules/07_apply` builds, one entry
    per REMOVED observation. Unique means that observation was removed by this criterion and no
    other; shared means at least one other criterion also fired on it. The two are separated
    because a criterion whose every removal is shared removes nothing on its own, and a total
    alone cannot show that.

    `criteria` names the full set the record was built over, so a criterion that fired on nothing
    still appears - evidence it was evaluated, which a table of only what fired cannot give.
    """
    names = [] if _unknown(criteria) else list(criteria)
    per = {k: {"unique": 0, "shared": 0, "total": 0} for k in names}
    n_removed = 0
    n_multi = 0
    for entry in rows:
        try:
            ident, why = entry
        except (TypeError, ValueError):
            raise TaskFailure(
                f"a removal-record row is not an (identifier, criteria) pair: {entry!r}") from None
        why = list(why)
        if not why:
            raise TaskFailure(
                f"observation {ident!r} appears in the removal record with no criterion. A "
                f"removed observation whose reason is empty cannot be questioned afterwards, "
                f"which is the whole purpose of the record")
        n_removed += 1
        if len(why) > 1:
            n_multi += 1
        for k in why:
            if k not in per:
                per[k] = {"unique": 0, "shared": 0, "total": 0}
                names.append(k)
            per[k]["total"] += 1
            per[k]["unique" if len(why) == 1 else "shared"] += 1
    return {"per_criterion": {k: per[k] for k in names},
            "n_removed": n_removed,
            "n_multi_criterion": n_multi}


# --------------------------------------------------------------------------------- matplotlib


def _mpl():
    """Import matplotlib with the Agg backend, inside the call that needs it.

    `force=False` because every figure below attaches its own Agg canvas: the global backend does
    not decide how these render, so switching a caller's interactive backend out from under them
    would be a side effect with no benefit. Older matplotlib has no `force` keyword and is
    handled rather than required.
    """
    import matplotlib
    try:
        matplotlib.use("Agg", force=False)
    except TypeError:                                                     # pragma: no cover
        matplotlib.use("Agg")
    return matplotlib


def _new_fig(figsize, dpi):
    """A Figure with an Agg canvas and no pyplot state.

    The object API rather than pyplot, so nothing accumulates in a global figure registry that
    the caller then has to remember to close.
    """
    _mpl()
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=figsize, dpi=dpi)
    FigureCanvasAgg(fig)
    return fig


def _log_axis(ax, which: str = "x") -> None:
    """A log scale whose ticks read in the original units - 1,000 rather than 10^3."""
    from matplotlib.ticker import ScalarFormatter

    axis = ax.xaxis if which == "x" else ax.yaxis
    (ax.set_xscale if which == "x" else ax.set_yscale)("log")
    fmt = ScalarFormatter()
    fmt.set_scientific(False)
    axis.set_major_formatter(fmt)


def _label(ax, x, y, text: str, **kw):
    """A mark placed in DATA coordinates - a threshold's value, a point's name.

    `set_in_layout(False)` on every one of them, because tight_layout measures an axes by the
    bounding box of its children: a long annotation makes the panel shrink to accommodate the
    text, and a set of panels that shrink by different amounts no longer share a scale. The one
    figure rule these labels exist to serve would then be broken by the labels themselves.
    """
    t = ax.text(x, y, text, **kw)
    t.set_in_layout(False)
    return t


def _annot(ax, text: str, xy, **kw):
    """A point's name, attached to the point and excluded from the layout for the same reason."""
    a = ax.annotate(text, xy, **kw)
    a.set_in_layout(False)
    return a


def _panel_note(ax, text: str, *, loc: str = "upper left", colour=None, wrap: int = 44):
    """A short annotation in AXES coordinates: n, a dropped-point count, a NOT SUPPLIED."""
    xy = {"upper left": (0.02, 0.97, "left", "top"),
          "upper right": (0.98, 0.97, "right", "top"),
          "lower left": (0.02, 0.03, "left", "bottom"),
          "lower right": (0.98, 0.03, "right", "bottom")}[loc]
    lines = []
    for para in str(text).split("\n"):
        lines += textwrap.wrap(para, wrap) or [""]
    # `_unknown`, not `colour or PALETTE["muted"]`: a NaN colour is truthy and would be handed
    # straight to matplotlib, which raises on it - the same defect as the labels, one layer down.
    t = ax.text(xy[0], xy[1], "\n".join(lines), transform=ax.transAxes, ha=xy[2], va=xy[3],
                fontsize=7, color=PALETTE["muted"] if _unknown(colour) else colour,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5})
    t.set_in_layout(False)
    return t


def _tidy(ax) -> None:
    ax.grid(True, which="major", linewidth=0.4, alpha=0.35)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _finish(fig, title: str, subtitle: str = "") -> None:
    head = title if not subtitle else f"{title}\n{subtitle}"
    fig.suptitle(head, fontsize=10 if not subtitle else 9, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.94 if subtitle else 0.96))


def _blank_panel(ax, message: str) -> None:
    """A panel for something that was never computed. It says so; it does not draw zero."""
    ax.text(0.5, 0.5, message, transform=ax.transAxes, ha="center", va="center", fontsize=8,
            color=PALETTE["unknown"], wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_color(PALETTE["unknown"])
        ax.spines[side].set_linestyle("--")


# ------------------------------------------------------------------------------- F1 · ingest


def fig_f1_barcode_rank(libraries, *, called_cells=None, dpi: int = DEFAULT_DPI):
    """F1 - is this really raw, unfiltered input?

    `libraries` maps a sample to `{"ranks": [...], "counts": [...], "n_barcodes": int|None}`;
    `called_cells` maps a sample to the number of barcodes the aligner called, drawn as a
    vertical line at that rank. A raw droplet matrix runs down into the single digits on the
    right of the curve; a matrix that has already been cell-called stops at a floor, and the
    floor is usually a round number somebody chose. Both are visible here and neither is visible
    in the file's shape or its name.

    Every panel shares limits: an unfiltered library and a pre-filtered one drawn to their own
    axes look identical.
    """
    if not libraries:
        raise TaskFailure("F1 needs at least one library's barcode-rank curve")
    called_cells = _mapping(called_cells)
    rows, cols = grid_shape(len(libraries))
    fig = _new_fig((3.1 * cols, 2.6 * rows), dpi)

    prepared, xs, ys = {}, [], []
    for s, entry in libraries.items():
        x, y = _series(s, entry, "ranks", "counts")
        px, py, dropped = _positive_pairs(x, y)
        if not px:
            raise TaskFailure(f"library {s!r}: every point is at or below zero, so a log-log "
                              f"barcode rank cannot be drawn from it")
        prepared[s] = (px, py, dropped, entry.get("n_barcodes"))
        xs += [min(px), max(px)]
        ys += [min(py), max(py)]

    for i, (s, (px, py, dropped, n_bc)) in enumerate(sorted(prepared.items()), start=1):
        ax = fig.add_subplot(rows, cols, i)
        ax.plot(px, py, linewidth=1.1, color=PALETTE["ok"])
        cut = called_cells.get(s)
        if not _unknown(cut) and float(cut) > 0:
            ax.axvline(float(cut), color=PALETTE["threshold"], linestyle="--", linewidth=1.0)
            _label(ax, float(cut), max(ys), f"  aligner cut: {int(cut):,}", rotation=90,
                   fontsize=6.5, va="top", ha="left", color=PALETTE["threshold"])
        _log_axis(ax, "x")
        _log_axis(ax, "y")
        ax.set_xlim(min(xs) * 0.8, max(xs) * 1.2)
        ax.set_ylim(min(ys) * 0.8, max(ys) * 1.2)
        ax.set_title(str(s), fontsize=8)
        note = n_label(n_bc, "barcodes")
        if _unknown(n_bc):
            note = f"points plotted = {len(px):,}; total {NOT_SUPPLIED}"
        if dropped:
            note += f"\n{dropped:,} point(s) at zero omitted by the log axis"
        _panel_note(ax, note, loc="lower left")
        _tidy(ax)
        if i % cols == 1 or cols == 1:
            ax.set_ylabel("UMI per barcode", fontsize=8)
        ax.set_xlabel("barcode rank", fontsize=8)
        ax.tick_params(labelsize=7)

    _finish(fig, "F1 · barcode rank per library (log-log, shared scale)",
            "the aligner's cell call is drawn where it was supplied; a floor on the right of a "
            "curve means the empties are already gone")
    return fig


# ------------------------------------------------------------------------------ F2 · ambient


def fig_f2_ambient_removal(per_library, *, per_gene_fractions=None, by_design=None,
                           gene_gut_threshold=None, material_floor=None, ratio_threshold=None,
                           count_denominator=None, dpi: int = DEFAULT_DPI):
    """F2 - how much ambient was removed, from what, and evenly across the design?

    The only figure in the report that can show a technical removal becoming an apparent
    biological difference, which is why the design panel is drawn even when there is nothing to
    put in it: an absent design map is stated on the panel rather than leaving the reader with
    two panels and no idea the third was possible.

    `per_library` maps a sample to the fraction of its counts removed. `per_gene_fractions` is a
    flat list of per-gene fractions. `by_design` maps a factor to `{level: rate}` or to
    `{level: {"rate": r, "n_samples": k, "denominator": "..."}}`. Thresholds are drawn only when
    passed; none has a default, because the value that binds is the caller's, not this module's.
    """
    if not per_library:
        raise TaskFailure("F2 needs at least one library's removed fraction")
    fig = _new_fig((12.0, 3.6), dpi)

    # (a) per library, with the cohort median drawn - a computed line, not a chosen one.
    ax = fig.add_subplot(1, 3, 1)
    names = sorted(per_library)
    known = [(s, float(per_library[s])) for s in names if not _unknown(per_library[s])]
    absent = [s for s in names if _unknown(per_library[s])]
    if known:
        ax.bar(range(len(known)), [100 * v for _, v in known], color=PALETTE["ok"], width=0.7)
        # A real median - the mean of the two middle values at even n. The label says median, so
        # the line has to be one; see `_median`.
        med = _median(v for _, v in known)
        if med is not None:
            ax.axhline(100 * med, color=PALETTE["threshold"], linestyle="--", linewidth=1.0)
            _label(ax, len(known) - 0.5, 100 * med, f" median {100 * med:.1f}%", fontsize=7,
                   va="bottom", ha="right", color=PALETTE["threshold"])
    for j, s in enumerate(absent, start=len(known)):
        ax.bar([j], [0], color=PALETTE["unknown"], hatch=UNKNOWN_HATCH, width=0.7)
        _label(ax, j, 0, f" {NOT_SUPPLIED}", rotation=90, fontsize=6.5, color=PALETTE["unknown"],
               va="bottom", ha="center")
    ax.set_xticks(range(len(known) + len(absent)))
    ax.set_xticklabels([s for s, _ in known] + absent, rotation=90, fontsize=7)
    ax.set_ylabel(f"% of counts removed\n({denominator_label(count_denominator)})", fontsize=8)
    ax.set_title("per library", fontsize=9)
    _panel_note(ax, n_label(len(per_library), "libraries"), loc="upper right")
    _tidy(ax)

    # (b) per gene: a shave is not a deletion, and only the distribution separates them.
    ax = fig.add_subplot(1, 3, 2)
    if _unknown(per_gene_fractions):
        _blank_panel(ax, "per-gene table NOT SUPPLIED\ngene gutting was not examined")
    else:
        vals = [100 * float(v) for v in per_gene_fractions if not _unknown(v)]
        skipped = sum(1 for v in per_gene_fractions if _unknown(v))
        if not vals:
            _blank_panel(ax, "per-gene table supplied but every value is unknown")
        else:
            ax.hist(vals, bins=40, color=PALETTE["neutral"])
            _log_axis(ax, "y")
            note = n_label(len(vals), "genes")
            if skipped:
                note += f"\n{skipped:,} gene(s) with no value, not counted"
            if not _unknown(gene_gut_threshold):
                thr = 100 * float(gene_gut_threshold)
                over = sum(1 for v in vals if v >= thr)
                ax.axvline(thr, color=PALETTE["refuse"], linestyle="--", linewidth=1.1)
                note += f"\n{over:,} gene(s) at or above the {thr:.0f}% line"
            _panel_note(ax, note, loc="upper right")
            ax.set_xlabel("% of a gene's counts removed", fontsize=8)
            ax.set_ylabel("genes (log)", fontsize=8)
    ax.set_title("per gene", fontsize=9)
    _tidy(ax)

    # (c) the panel that matters: is the removal even across the design?
    ax = fig.add_subplot(1, 3, 3)
    if not _mapping(by_design):
        _blank_panel(ax, "no design map supplied\nthe check that matters most was NOT RUN")
    else:
        pos, labels, heights, colours, boundaries = 0, [], [], [], []
        notes = []
        for factor in sorted(by_design):
            levels = by_design[factor]
            rates = []
            for level in sorted(levels, key=str):
                rate, n_samples, denom = _level_entry(levels[level])
                labels.append(f"{factor}={level}"
                              + ("" if _unknown(n_samples) else f" (n={int(n_samples)})"))
                if _unknown(rate):
                    heights.append(0.0)
                    colours.append(PALETTE["unknown"])
                else:
                    heights.append(100 * float(rate))
                    colours.append(PALETTE["ok"])
                    rates.append(float(rate))
                pos += 1
            boundaries.append(pos - 0.5)
            if len(rates) >= 2 and min(rates) > 0:
                notes.append(f"{factor}: max/min = {max(rates) / min(rates):.2f}x")
            elif len(rates) >= 2:
                notes.append(f"{factor}: max/min UNDEFINED (a level removed nothing)")
            else:
                notes.append(f"{factor}: NOT CHECKED ({len(rates)} level(s) with a rate)")
        bars = ax.bar(range(len(heights)), heights, color=colours, width=0.65)
        # Set per-bar, not with a list `hatch=` kwarg: older matplotlib takes only a scalar
        # there and would raise on the one figure whose whole point is marking what is unknown.
        for bar, colour in zip(bars, colours):
            if colour == PALETTE["unknown"]:
                bar.set_hatch(UNKNOWN_HATCH)
        for b in boundaries[:-1]:
            ax.axvline(b, color=PALETTE["muted"], linewidth=0.6, linestyle=":")
        if not _unknown(material_floor):
            ax.axhline(100 * float(material_floor), color=PALETTE["threshold"], linestyle="-.",
                       linewidth=1.0)
            _label(ax, len(heights) - 0.5, 100 * float(material_floor),
                   f" materiality floor {100 * float(material_floor):.0f}%", fontsize=6.5,
                   ha="right", va="bottom", color=PALETTE["threshold"])
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
        ax.set_ylabel("% of counts removed", fontsize=8)
        if not _unknown(ratio_threshold):
            notes.append(f"refusal line {float(ratio_threshold):.0f}x, where the removal is "
                         f"also material")
        # Headroom for the per-factor ratios, so the note sits above the bars rather than over
        # them: a number printed across the data it describes is read against the wrong bar.
        tallest = max(heights) if heights else 1.0
        ax.set_ylim(0, (tallest if tallest > 0 else 1.0) * (1.25 + 0.13 * len(notes)))
        _panel_note(ax, "\n".join(notes), loc="upper right", wrap=40)
    ax.set_title("by design level", fontsize=9)
    _tidy(ax)

    _finish(fig, "F2 · ambient removal: how much, from what, and evenly?",
            "a bar in the never-examined colour is a value that was not supplied - it is not a "
            "zero")
    return fig


def _level_entry(entry) -> tuple:
    """(rate, n_samples, denominator) from a design level given as a number or as a dict."""
    if isinstance(entry, dict):
        return entry.get("rate"), entry.get("n_samples"), entry.get("denominator")
    return entry, None, None


# ---------------------------------------------------------------------------- F3 · cell call


def fig_f3_cell_calls(calls, *, review_lost=None, refuse_lost=None, dpi: int = DEFAULT_DPI):
    """F3 - did the denoiser drop cells the aligner kept?

    `calls` maps a sample to `{"aligner": n, "cellbender": n, "lost": n|None}`. The left panel
    draws the equality line: below it the denoiser is calling fewer cells than the aligner, which
    is the point at which a cell-selection decision has migrated into a tool chosen for
    denoising. The right panel draws the REVIEW and REFUSE lines at the values the caller passed
    - a library sitting on one of them is visible without reading a table.

    A library whose loss was not supplied gets a marked slot in the never-examined colour and no
    bar. Drawing it at zero would report the best possible outcome for a quantity nobody measured.
    """
    if not calls:
        raise TaskFailure("F3 needs at least one library's cell calls")
    fig = _new_fig((9.5, 3.8), dpi)
    names = sorted(calls)

    ax = fig.add_subplot(1, 2, 1)
    xs, ys = [], []
    for s in names:
        c = calls[s]
        a, b = c.get("aligner"), c.get("cellbender")
        if _unknown(a) or _unknown(b):
            continue
        xs.append(float(a))
        ys.append(float(b))
        _annot(ax, str(s), (float(a), float(b)), fontsize=6, xytext=(3, 3),
               textcoords="offset points", color=PALETTE["muted"])
    if not xs:
        _blank_panel(ax, "no library supplied both an aligner and a denoiser count")
    else:
        ax.scatter(xs, ys, s=26, color=PALETTE["ok"], marker="o", zorder=3)
        lim = (min(xs + ys) * 0.9, max(xs + ys) * 1.1)
        ax.plot(lim, lim, color=PALETTE["threshold"], linestyle="--", linewidth=1.0)
        _label(ax, lim[1], lim[1], " equal calls", fontsize=7, ha="right", va="bottom",
               color=PALETTE["threshold"])
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        below = sum(1 for a, b in zip(xs, ys) if b < a)
        _panel_note(ax, f"{n_label(len(xs), 'libraries')}\n{below:,} below the equality line",
                    loc="upper left")
    ax.set_xlabel("cells called by the aligner", fontsize=8)
    ax.set_ylabel("cells called by the denoiser", fontsize=8)
    ax.set_title("caller against caller", fontsize=9)
    _tidy(ax)

    ax = fig.add_subplot(1, 2, 2)
    heights, colours, missing = [], [], []
    for j, s in enumerate(names):
        c = calls[s]
        lost, aligner = c.get("lost"), c.get("aligner")
        if _unknown(lost) or _unknown(aligner) or float(aligner) <= 0:
            heights.append(0.0)
            colours.append(PALETTE["unknown"])
            missing.append(j)
        else:
            frac = 100 * float(lost) / float(aligner)
            heights.append(frac)
            colours.append(PALETTE["refuse"] if (not _unknown(refuse_lost)
                                                 and frac > 100 * float(refuse_lost))
                           else (PALETTE["review"] if (not _unknown(review_lost)
                                                       and frac >= 100 * float(review_lost))
                                 else PALETTE["ok"]))
    bars = ax.bar(range(len(names)), heights, color=colours, width=0.7)
    for j in missing:
        bars[j].set_hatch(UNKNOWN_HATCH)
    for j in missing:
        _label(ax, j, 0, f" {NOT_SUPPLIED}", rotation=90, fontsize=6.5, ha="center", va="bottom",
               color=PALETTE["unknown"])
    for value, style, tag in ((review_lost, ":", "REVIEW"), (refuse_lost, "--", "REFUSE")):
        if not _unknown(value):
            ax.axhline(100 * float(value), color=PALETTE["threshold"], linestyle=style,
                       linewidth=1.0)
            _label(ax, len(names) - 0.5, 100 * float(value), f" {tag} {100 * float(value):.0f}%",
                   fontsize=6.5, ha="right", va="bottom", color=PALETTE["threshold"])
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=7)
    ax.set_ylabel("% of aligner calls lost\n(denominator: that library's aligner calls)",
                  fontsize=8)
    ax.set_title("loss at the cell call", fontsize=9)
    _panel_note(ax, n_label(len(names), "libraries"), loc="upper right")
    _tidy(ax)

    _finish(fig, "F3 · aligner against denoiser, with the gate lines drawn",
            "a permissive denoiser is the safe error; a strict one has become a filter")
    return fig


# --------------------------------------------------------------------------- F4 · light floor


def fig_f4_scoring_coverage(coverage, *, floor=None, umi_hist=None, dpi: int = DEFAULT_DPI):
    """F4 - what was never examined?

    `coverage` maps a sample to `{"scored": n, "below_floor": n, "never_scored": n}`. The third
    category is the whole point of the figure and is never merged into the first two: a nucleus
    below the light floor was not examined for a stated reason, a nucleus in `never_scored` was
    not examined for some other reason, and neither is a nucleus that was examined and cleared.

    `umi_hist` optionally supplies `{"edges": [...], "counts": [...]}` so the floor itself can be
    DRAWN on the axis it applies to; without it the floor is only nameable, which is why the
    parameter exists.
    """
    if not coverage:
        raise TaskFailure("F4 needs at least one library's scoring coverage")
    has_hist = not _unknown(umi_hist) and bool(umi_hist)
    panels = 2 if has_hist else 1
    fig = _new_fig((6.0 * panels, max(2.6, 0.34 * len(coverage) + 1.6)), dpi)
    names = sorted(coverage)

    ax = fig.add_subplot(1, panels, 1)
    cats = (("scored", PALETTE["ok"], ""),
            ("below_floor", PALETTE["review"], ".."),
            ("never_scored", PALETTE["unknown"], UNKNOWN_HATCH))
    left = [0.0] * len(names)
    incomplete = []
    for key, colour, hatch in cats:
        widths = []
        for j, s in enumerate(names):
            v = coverage[s].get(key)
            if _unknown(v):
                widths.append(0.0)
                if j not in incomplete:
                    incomplete.append(j)
            else:
                widths.append(float(v))
        ax.barh(range(len(names)), widths, left=left, color=colour, hatch=hatch, height=0.7,
                label=key.replace("_", " "))
        left = [a + b for a, b in zip(left, widths)]
    for j, s in enumerate(names):
        tag = f" {left[j]:,.0f}"
        if j in incomplete:
            tag += f"  ({NOT_SUPPLIED} in at least one category - this total is a FLOOR)"
        _label(ax, left[j], j, tag, va="center", fontsize=6.5, color=PALETTE["muted"])
    # Room for the totals, which are excluded from the layout so they cannot shrink the bars.
    # A cohort whose every count is unknown has a widest bar of zero; the axis is given a unit
    # range rather than a singular one, so a legitimately empty figure does not emit a warning
    # that reads, in a task log, like a defect.
    widest = max(left) if left else 0.0
    ax.set_xlim(0, (widest if widest > 0 else 1.0) * 1.45)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("nuclei (denominator: every nucleus the denoiser called, per library)",
                  fontsize=8)
    ax.legend(fontsize=7, loc="lower right", frameon=False)
    ax.set_title("who was handed to the doublet detector"
                 + ("" if _unknown(floor) else f" (floor {int(floor):,} UMI)"), fontsize=9)
    _panel_note(ax, n_label(len(names), "libraries"), loc="upper right")
    _tidy(ax)

    if has_hist:
        ax = fig.add_subplot(1, 2, 2)
        edges, counts = umi_hist.get("edges"), umi_hist.get("counts")
        if _unknown(edges) or _unknown(counts):
            _blank_panel(ax, "umi_hist supplied without 'edges' and 'counts'")
        else:
            edges, counts = list(edges), list(counts)
            if len(edges) != len(counts) + 1:
                raise TaskFailure(f"umi_hist: {len(edges):,} edges for {len(counts):,} bins; a "
                                  f"histogram needs one more edge than it has bins")
            widths = [b - a for a, b in zip(edges[:-1], edges[1:])]
            ax.bar(edges[:-1], counts, width=widths, align="edge", color=PALETTE["neutral"])
            if not _unknown(floor):
                ax.axvline(float(floor), color=PALETTE["threshold"], linestyle="--",
                           linewidth=1.1)
                _label(ax, float(floor), max(counts), f"  floor {int(floor):,} UMI", fontsize=7,
                       rotation=90, va="top", color=PALETTE["threshold"])
            _log_axis(ax, "x")
            ax.set_xlabel("UMI per nucleus", fontsize=8)
            ax.set_ylabel("nuclei", fontsize=8)
            ax.set_title("where the floor sits", fontsize=9)
            _panel_note(ax, n_label(sum(counts), "nuclei in the histogram"), loc="upper right")
        _tidy(ax)

    _finish(fig, "F4 · scored, below the floor, and never scored",
            "never scored is its own category and its own hatch; it is not a negative result")
    return fig


# ------------------------------------------------------------------------------ F5 · doublets


def fig_f5_doublet_sweep(sweep, *, published_band=None, chosen=None, rate_denominator=None,
                         param_label=None, dpi: int = DEFAULT_DPI):
    """F5 - is the doublet rate a measurement or the detector's prior?

    `sweep` maps a sample to `{"x": [prior values swept], "y": [called rate at each]}`. A rate
    that tracks the swept prior across its whole range was set by the prior; one that stays put
    was measured from the library. Both readings need the same y scale across libraries, so
    there is one axes rather than a panel each.

    `published_band` is `{"low": f, "high": f, "label": "..."}` and the label is REQUIRED: a
    shaded band with no stated source is a number the figure asserts and the code never computed.
    """
    if not sweep:
        raise TaskFailure("F5 needs at least one library's sweep")
    fig = _new_fig((7.6, 4.2), dpi)
    ax = fig.add_subplot(1, 1, 1)

    if not _unknown(published_band):
        low, high = published_band.get("low"), published_band.get("high")
        label = published_band.get("label")
        if _unknown(low) or _unknown(high):
            raise TaskFailure("F5: published_band needs both 'low' and 'high'")
        if _unknown(label) or not str(label).strip():
            raise TaskFailure(
                "F5: published_band needs a 'label' naming where the band comes from. An "
                "unattributed shaded region states a range the pipeline did not compute and "
                "cannot be checked against anything")
        ax.axhspan(100 * float(low), 100 * float(high), color=PALETTE["muted"], alpha=0.15)
        _label(ax, 0.99, 100 * float(high), f"{label}  ", transform=ax.get_yaxis_transform(),
               ha="right", va="bottom", fontsize=7, color=PALETTE["muted"])

    n_points = 0
    for s in sorted(sweep):
        x, y = _series(s, sweep[s], "x", "y")
        n_points += len(x)
        ax.plot([float(v) for v in x], [100 * float(v) for v in y], linewidth=1.0, marker="o",
                markersize=2.5, label=str(s))
    if not _unknown(chosen):
        ax.axvline(float(chosen), color=PALETTE["threshold"], linestyle="--", linewidth=1.1)
        _label(ax, float(chosen), ax.get_ylim()[1], f"  applied: {chosen}", rotation=90,
               fontsize=7, va="top", color=PALETTE["threshold"])
    ax.set_xlabel(label_text(param_label, f"swept parameter ({NOT_SUPPLIED} - name it)"),
                  fontsize=8)
    ax.set_ylabel(f"% called a doublet\n({denominator_label(rate_denominator)})", fontsize=8)
    ax.legend(fontsize=6.5, ncol=2, frameon=False)
    _panel_note(ax, f"{n_label(len(sweep), 'libraries')}   {n_points:,} swept points",
                loc="upper left")
    _tidy(ax)
    _finish(fig, "F5 · doublet rate against the detector's expected-rate prior",
            "a curve that follows the prior is a rate the prior set")
    return fig


# ------------------------------------------------------------------------------- F6 · quality


def fig_f6_quality_density(densities, *, valleys=None, cut=None, bounds=None, metric_label=None,
                           dpi: int = DEFAULT_DPI):
    """F6 - where is the cut, and why there?

    `densities` maps a sample to `{"x": [...], "y": [...], "n": int|None}` - the density the
    caller estimated, not one recomputed here. `valleys` maps a sample to its own measured
    minimum; `cut` is the single cohort constant, drawn on EVERY panel at the same place so a
    library the constant fits poorly is visible as a line sitting away from that library's own
    valley. `bounds` shades the range outside which a derived floor is refused.

    Alongside F2 this is the figure that carries the report: it is the only one that shows a
    derived threshold sitting where the data actually separates - or not.
    """
    if not densities:
        raise TaskFailure("F6 needs at least one library's density")
    valleys = _mapping(valleys)
    rows, cols = grid_shape(len(densities))
    fig = _new_fig((3.2 * cols, 2.5 * rows), dpi)

    prepared, xs, ys = {}, [], []
    for s in densities:
        x, y = _series(s, densities[s], "x", "y")
        prepared[s] = ([float(v) for v in x], [float(v) for v in y], densities[s].get("n"))
        xs += [min(prepared[s][0]), max(prepared[s][0])]
        ys += [0.0, max(prepared[s][1])]

    for i, (s, (x, y, n)) in enumerate(sorted(prepared.items()), start=1):
        ax = fig.add_subplot(rows, cols, i)
        if (not _unknown(bounds) and bounds is not None and len(bounds) == 2
                and not any(_unknown(b) for b in bounds)):
            ax.axvspan(min(xs), float(bounds[0]), color=PALETTE["muted"], alpha=0.12)
            ax.axvspan(float(bounds[1]), max(xs), color=PALETTE["muted"], alpha=0.12)
        ax.plot(x, y, linewidth=1.1, color=PALETTE["ok"])
        v = valleys.get(s)
        if not _unknown(v):
            ax.axvline(float(v), color=PALETTE["second"], linestyle=":", linewidth=1.1)
            _label(ax, float(v), max(ys), f" valley {float(v):,.0f}", fontsize=6.5, rotation=90,
                   va="top", color=PALETTE["second"])
        if not _unknown(cut):
            ax.axvline(float(cut), color=PALETTE["threshold"], linestyle="--", linewidth=1.2)
            _label(ax, float(cut), max(ys), f" cut {float(cut):,.0f}", fontsize=6.5, rotation=90,
                   va="top", ha="right", color=PALETTE["threshold"])
        ax.set_xlim(min(xs), max(xs))
        # A density that is flat at zero everywhere is a real input - it is what a caller passes
        # when the estimate failed - and it must not produce a singular axis warning on top of it.
        ax.set_ylim(0, (max(ys) if max(ys) > 0 else 1.0) * 1.08)
        ax.set_title(str(s), fontsize=8)
        note = n_label(n, "nuclei")
        if _unknown(v):
            note += "\nvalley NOT SUPPLIED"
        _panel_note(ax, note, loc="upper right")
        ax.set_xlabel(label_text(metric_label, f"metric ({NOT_SUPPLIED} - name it)"),
                      fontsize=7.5)
        ax.tick_params(labelsize=7)
        _tidy(ax)

    _finish(fig, "F6 · per-library density with the measured valley and the applied cut",
            "the same cut is drawn on every panel at the same scale; a shaded margin is outside "
            "the bounds a derived floor is allowed to take")
    return fig


def _violin_row(ax, distributions, names, *, log, cut, cap, metric_name):
    """One panel of paired per-library violins. Returns (positions, n_above_cap).

    PAIRED, AND LABELLED ONCE PER LIBRARY. Each library gets a `before` violin and an `after`
    violin side by side; the axis is labelled with the library, not with four lines of state and
    n under every shape. A reader compares within a pair first and across libraries second, and
    the layout should make the first of those the easy one.
    """
    import math

    data, positions, colours, ticks, tick_at = [], [], [], [], []
    dropped = above = 0
    for j, s in enumerate(names):
        entry = distributions[s] or {}
        centre = j * 2.2
        for k, (state, colour) in enumerate((("before", PALETTE["before"]),
                                             ("after", PALETTE["after"]))):
            raw = entry.get(state)
            if _unknown(raw):
                continue
            vals = [float(v) for v in raw if not _unknown(v)]
            if cap is not None:
                kept = [v for v in vals if v <= cap]
                above += len(vals) - len(kept)
                vals = kept
            if log:
                positive = [math.log10(v) for v in vals if v > 0]
                dropped += len(vals) - len(positive)
                vals = positive
            if not vals:
                continue
            data.append(vals)
            positions.append(centre + k * 0.85)
            colours.append(colour)
        ticks.append(s)
        tick_at.append(centre + 0.42)

    if not data:
        _blank_panel(ax, "no usable values were supplied")
        return [], 0

    parts = ax.violinplot(data, positions=positions, widths=0.78, showextrema=False,
                          showmedians=False)
    for body, colour in zip(parts["bodies"], colours):
        body.set_facecolor(colour)
        body.set_alpha(1.0 if colour == PALETTE["after"] else 0.9)
        body.set_edgecolor("none")

    # THE MEDIAN IS DRAWN AND ALSO WRITTEN. A rule says where the centre is; the number says what
    # it is, and a reader comparing libraries should not have to read a value off a log axis.
    for vals, pos, colour in zip(data, positions, colours):
        ordered = sorted(vals)
        m = ordered[len(ordered) // 2]
        ax.hlines(m, pos - 0.36, pos + 0.36, color=PALETTE["ok"], linewidth=2.0, zorder=4)
        if colour == PALETTE["after"]:
            shown = (10 ** m) if log else m
            ax.annotate(f"{shown:,.0f}", (pos, m), textcoords="offset points", xytext=(0, 6),
                        ha="center", fontsize=6.5, fontweight="bold", zorder=5)

    if not _unknown(cut):
        value = float(cut)
        if not (log and value <= 0):
            level = math.log10(value) if log else value
            ax.axhline(level, color=PALETTE["refuse"], linestyle="--", linewidth=1.2, zorder=3)
            _label(ax, max(positions), level, f" applied floor = {value:,.0f}", fontsize=6.5,
                   ha="right", va="bottom", color=PALETTE["refuse"])

    ax.set_xticks(tick_at)
    ax.set_xticklabels(ticks, fontsize=7, rotation=32, ha="right")
    ax.margins(y=0.14)
    if log:
        from matplotlib.ticker import FuncFormatter, MultipleLocator
        lo, hi = ax.get_ylim()
        if hi - lo >= 1.2:
            ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{10 ** v:,.0f}"))
        ax.set_ylabel(f"{metric_name} per nucleus", fontsize=8)
        if dropped:
            _panel_note(ax, f"{dropped:,} value(s) at or below zero cannot be shown on a log "
                            f"axis", loc="upper right")
    else:
        ax.set_ylabel(f"{metric_name} per nucleus", fontsize=8)
    _tidy(ax)
    return positions, above


def fig_f7_before_after(distributions, *, cut=None, metric_label=None, scale: str = "both",
                        cap=None, fig_id: str = "F7", dpi: int = DEFAULT_DPI):
    """F7 - what did the cut change?

    `distributions` maps a sample to `{"before": [...], "after": [...]}` of the metric's values.
    One pair of violins per library: the state before the cut in grey, what survived it in black,
    the median drawn AND printed, and the applied floor across both.

    THE TWO AXES ANSWER DIFFERENT QUESTIONS AND BOTH ARE NEEDED.

    A log axis shows that two populations exist and that the cut fell between them - on a linear
    axis the population being removed is a smear against the bottom and invisible, so the figure
    that justifies the threshold has to be the log one. A linear axis shows the scale a reader
    actually works in, and how little of the retained distribution sits anywhere near the cut.

    `scale` selects one - "log", "linear" - or "both", which stacks them. The report draws them as
    two figures so each can carry its own caption; "both" is kept because a reader looking at one
    figure in isolation is the case this is guarding against.

    A CAPPED AXIS MUST SAY IT IS CAPPED. `cap` clips the linear panel so the bulk is legible
    rather than compressed into the bottom decile by a long tail. The proportion of nuclei that
    fall outside is then counted and printed on the panel. A truncated axis that does not say so
    reports a distribution nobody drew, and a reader has no way to know the difference.

    Values at or below zero cannot appear on a log axis; they are counted and the count printed,
    rather than dropped where nobody can see how many went.
    """
    if not distributions:
        raise TaskFailure("F7 needs at least one library's before/after values")
    if scale not in ("both", "log", "linear"):
        raise TaskFailure(f"F7: scale must be 'both', 'log' or 'linear', got {scale!r}")

    names = sorted(distributions)
    # Read once, so two panels cannot disagree about what the metric is called, and through
    # `label_text` so an unsupplied name is stated rather than printed as `nan`.
    metric_name = label_text(metric_label, f"metric ({NOT_SUPPLIED} - name it)")
    panels = (("log",), ("linear",), ("log", "linear"))[("log", "linear", "both").index(scale)]
    width = max(6.5, 0.78 * len(names) + 2.6)
    fig = _new_fig((width, 3.4 * len(panels) + 0.5), dpi)

    total = 0
    for vals in distributions.values():
        for state in ("before", "after"):
            raw = (vals or {}).get(state)
            if not _unknown(raw):
                total += sum(1 for v in raw if not _unknown(v))

    for i, which in enumerate(panels):
        ax = fig.add_subplot(len(panels), 1, i + 1)
        use_cap = cap if (which == "linear" and cap is not None) else None
        _positions, above = _violin_row(ax, distributions, names, log=(which == "log"), cut=cut,
                                        cap=use_cap, metric_name=metric_name)
        if use_cap is not None:
            share = (100.0 * above / total) if total else 0.0
            _panel_note(ax, f"y capped at {float(use_cap):,.0f} - {share:.2f}% of plotted nuclei "
                            f"lie above and are not drawn", loc="upper left")

    axis_words = {"log": "a log axis", "linear": "a linear axis",
                  "both": "log and linear axes"}[scale]
    _finish(fig, f"{fig_id} - {metric_name} per library, before and after the cut, on {axis_words}",
            "grey is every called cell, black is what was retained; the median is drawn and "
            "printed, and the applied floor crosses every library")
    return fig


# ------------------------------------------------------------------------------ F8 · clusters


def fig_f8_cluster_flags(clusters, *, thresholds=None, x_key: str = "umi_frac_of_sample",
                         x_label=None, dpi: int = DEFAULT_DPI):
    """F8 - are any clusters technical?

    `clusters` is the step-6 profile: dicts carrying `sample`, `cluster`, the depth key named by
    `x_key`, `median_pct_mt`, and the three-valued `FLAG` and `WATCH`. A cluster whose FLAG was
    never determined was NOT EVALUATED - markers were not computed for it - and it is drawn in
    the never-examined colour with its own marker, never as a cluster that passed.

    `thresholds` may carry `a_umi_frac` and `b_pct_mt`; each is drawn as a line so a cluster that
    fell one side of a cut by a hair is visible as such. The booleans are a summary of two
    thresholds at once, and a reader cannot see from a FLAG column which of them a cluster failed.

    HOW A ROW IS CLASSIFIED, AND THE TWO WAYS THIS WAS WRONG

    Both flags are read through `_tri`, never through `is True` and never through a bare truth
    test. `numpy.bool_(True) is True` is False, so a genuinely FLAGGED cluster arriving from any
    numpy-backed profile used to be drawn as `clear` - the figure showed the opposite of the
    finding. And a blank cell arrives from pandas as NaN, not as None, so a cluster whose markers
    were NEVER COMPUTED used to fall through every branch to `clear` and be drawn identically to
    one that was computed and passed. never-examined is its own category for exactly that reading
    (docs/REPORT_DESIGN.md).

    The precedence, which is a judgement call and is recorded here rather than left in the code:

      1. FLAG determined true      -> FLAG. An established flag is the strongest statement on the
                                     panel and is drawn even if WATCH was never determined.
      2. either flag undetermined  -> NOT EVALUATED, including when WATCH is true and FLAG is not
                                     known. Drawing WATCH there would assert the cluster had been
                                     examined against the flag criteria and had cleared them,
                                     which is precisely the assertion nobody made.
      3. WATCH determined true     -> WATCH.
      4. both determined false     -> clear. This is the ONLY route to `clear`, and it requires
                                     two determinations rather than the absence of one.
    """
    if not clusters:
        raise TaskFailure("F8 needs at least one cluster row")
    thresholds = _mapping(thresholds)
    fig = _new_fig((6.8, 5.0), dpi)
    ax = fig.add_subplot(1, 1, 1)

    groups = {
        "FLAG": {"colour": PALETTE["refuse"], "marker": "X", "size": 58, "rows": []},
        "WATCH": {"colour": PALETTE["review"], "marker": "s", "size": 34, "rows": []},
        "clear": {"colour": PALETTE["ok"], "marker": "o", "size": 26, "rows": []},
        "NOT EVALUATED": {"colour": PALETTE["unknown"], "marker": "D", "size": 34, "rows": []},
    }
    skipped = 0
    for row in clusters:
        x, y = row.get(x_key), row.get("median_pct_mt")
        if _unknown(x) or _unknown(y):
            skipped += 1
            continue
        flag, watch = _tri(row.get("FLAG")), _tri(row.get("WATCH"))
        if flag is True:
            key = "FLAG"
        elif flag is None or watch is None:
            key = "NOT EVALUATED"
        elif watch is True:
            key = "WATCH"
        else:
            key = "clear"
        groups[key]["rows"].append((float(x), float(y), row))

    for key, g in groups.items():
        if not g["rows"]:
            continue
        ax.scatter([r[0] for r in g["rows"]], [r[1] for r in g["rows"]], s=g["size"],
                   c=g["colour"], marker=g["marker"], zorder=3, edgecolors="white",
                   linewidths=0.4, label=f"{key} ({len(g['rows']):,})")
        if key in ("FLAG", "WATCH"):
            for x, y, row in g["rows"]:
                # `dict.get(k, "?")` returns the DEFAULT only when the key is absent: a row
                # carrying `sample = NaN` - which is what a blank cell in a profile read from CSV
                # is - sails past it and labels the point `nan c3`.
                _annot(ax, f"{label_text(row.get('sample'), '?')} "
                           f"c{label_text(row.get('cluster'), '?')}", (x, y),
                       fontsize=5.5, xytext=(4, 3), textcoords="offset points",
                       color=PALETTE["muted"])

    for key, orientation, tag in (("a_umi_frac", "v", "A"), ("b_pct_mt", "h", "B")):
        value = thresholds.get(key)
        if _unknown(value):
            continue
        line = ax.axvline if orientation == "v" else ax.axhline
        line(float(value), color=PALETTE["threshold"], linestyle="--", linewidth=1.0)
        if orientation == "v":
            _label(ax, float(value), ax.get_ylim()[1], f" {tag} {value}", fontsize=7, rotation=90,
                   va="top", color=PALETTE["threshold"])
        else:
            _label(ax, ax.get_xlim()[1], float(value), f"{tag} {value} ", fontsize=7, ha="right",
                   va="bottom", color=PALETTE["threshold"])

    note = n_label(len(clusters), "clusters")
    if skipped:
        note += f"\n{skipped:,} cluster(s) missing a coordinate and not drawn"
    _panel_note(ax, note, loc="upper right")
    # The fallback is the key the depth was read from, which is itself an argument and so is
    # itself read through `label_text` rather than being assumed to be a string.
    ax.set_xlabel(label_text(x_label, label_text(x_key, f"depth key ({NOT_SUPPLIED})")
                             .replace("_", " ")), fontsize=8)
    ax.set_ylabel("cluster median % mitochondrial", fontsize=8)
    # Only where something was drawn: a legend call over an empty axes warns, and a warning in a
    # task log about a cohort whose coordinates were all missing points at the wrong problem -
    # the panel note above already states how many clusters could not be placed.
    if any(g["rows"] for g in groups.values()):
        ax.legend(fontsize=7, frameon=False, loc="best")
    _tidy(ax)
    _finish(fig, "F8 · cluster depth against mitochondrial content, with the flag cuts drawn",
            "NOT EVALUATED is a cluster whose markers were never computed; it is not a cluster "
            "that passed")
    return fig


# --------------------------------------------------------------------------------- F9 · apply


def fig_f9_criterion_contributions(per_criterion, *, n_removed=None, n_in=None,
                                   dpi: int = DEFAULT_DPI):
    """F9 - what did each criterion remove UNIQUELY?

    `per_criterion` maps a criterion to `{"unique": n, "shared": n}` - which is the
    `"per_criterion"` ENTRY of what `criterion_contributions()` returns, not its return value.
    That distinction is written out because passing the return value straight in fails here as
    `'int' object has no attribute 'get'`, having iterated `n_removed` as a criterion name. A
    criterion whose removals are all shared with another
    removed nothing on its own, and a total-only bar chart cannot show that; it is the number
    that decides whether dropping a criterion would change the deliverable at all.

    A criterion that fired on nothing is drawn at zero WITH its name, because a criterion that
    was evaluated and removed nothing is a result, and one that was never evaluated is not - the
    caller separates them by which criteria it passes in.
    """
    if not per_criterion:
        raise TaskFailure("F9 needs at least one criterion; a removal with no named criterion "
                          "cannot be questioned afterwards")
    fig = _new_fig((7.4, max(2.4, 0.44 * len(per_criterion) + 1.6)), dpi)
    ax = fig.add_subplot(1, 1, 1)

    names = list(per_criterion)
    uniq, shared, unknown_rows = [], [], []
    for j, k in enumerate(names):
        entry = per_criterion[k]
        u, s = entry.get("unique"), entry.get("shared")
        if _unknown(u) or _unknown(s):
            unknown_rows.append(j)
            uniq.append(0.0)
            shared.append(0.0)
        else:
            uniq.append(float(u))
            shared.append(float(s))

    ax.barh(range(len(names)), uniq, color=PALETTE["refuse"], height=0.6, label="removed by "
            "this criterion ALONE")
    ax.barh(range(len(names)), shared, left=uniq, color=PALETTE["neutral"], height=0.6,
            label="also removed by another criterion")
    # The widest STACKED bar, not the widest segment: `max(uniq + shared)` concatenates the two
    # lists and returns the largest single segment, which leaves the axis too short for the row
    # whose two segments are both middling - and its total label then runs off the figure.
    span = max((u + s) for u, s in zip(uniq, shared))
    for j in unknown_rows:
        ax.barh([j], [span if span > 0 else 1.0], color=PALETTE["unknown"],
                hatch=UNKNOWN_HATCH, height=0.6)
        _label(ax, 0, j, f"  {NOT_SUPPLIED}", va="center", fontsize=7, color="white")
    for j, k in enumerate(names):
        if j in unknown_rows:
            continue
        _label(ax, uniq[j] + shared[j], j,
               f"  {int(uniq[j]):,} unique / {int(uniq[j] + shared[j]):,} total",
               va="center", fontsize=6.5, color=PALETTE["muted"])

    ax.set_xlim(0, (span if span > 0 else 1.0) * 1.6)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7.5)
    denom = None
    if not _unknown(n_removed) and not _unknown(n_in):
        denom = f"{int(n_removed):,} removed of {int(n_in):,} observations"
    elif not _unknown(n_removed):
        denom = f"{int(n_removed):,} removed; the input total was not supplied"
    ax.set_xlabel(f"observations removed ({denominator_label(denom)})", fontsize=8)
    # Above the axes rather than inside it: a legend that lands on the longest bar hides the one
    # criterion the figure exists to draw attention to.
    ax.legend(fontsize=7, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=2)
    _panel_note(ax, n_label(len(names), "criteria"), loc="upper right")
    _tidy(ax)
    _finish(fig, "F9 · unique against shared contribution per removal criterion",
            "a criterion whose removals are entirely shared changes nothing if it is dropped")
    return fig


#: The report writer resolves a figure id to the function that draws it. Kept here rather than in
#: the writer so a new figure is added in one place.

# -------------------------------------------------------------------- F10 - per-library UMAP


def fig_f10_umap_per_library(embeddings, *, colour_by: str = "doublet", subtitle=None,
                             fig_id: str = "F10", dpi: int = DEFAULT_DPI):
    """F10 - where in the manifold did the removed nuclei sit?

    `embeddings` maps a sample to `{"x": [...], "y": [...], "doublet": [...], "removed": [...]}`,
    one entry per barcode. `colour_by` selects `doublet` - which barcodes the detector called -
    or `removed`, which barcodes any criterion removed.

    ONE PANEL PER LIBRARY, AND NO COHORT PANEL. Libraries are embedded separately because nothing
    in this pipeline pools them: clustering, the mitochondrial ceiling and the cluster check are
    all per library, and a pooled embedding would be an object no step produces. Drawing one
    would also invite a reader to read cross-library structure off a projection that no batch
    correction has been applied to.

    THE SAME COORDINATES ARE USED FOR BOTH CALLS, and that is the whole point of taking an
    embedding as input rather than computing one here. Re-embedding the retained cells produces a
    different layout, and a reader comparing the two would be looking at a difference that may be
    the projection rather than the data.

    What the figure can establish is whether removed nuclei are SCATTERED or CONCENTRATED. Doublets
    that fall on the boundaries between clusters are behaving as doublets should. A coherent
    region removed in one piece is a population leaving, which may be correct and is never
    something to discover after the fact.
    """
    if not embeddings:
        raise TaskFailure("F10 needs at least one library's embedding")
    if colour_by not in ("doublet", "removed"):
        raise TaskFailure(f"F10: colour_by must be 'doublet' or 'removed', got {colour_by!r}")

    names = sorted(embeddings)
    rows, cols = grid_shape(len(names))
    fig = _new_fig((min(12.0, 2.15 * cols + 0.6), 2.05 * rows + 0.9), dpi)
    highlight = PALETTE["refuse"] if colour_by == "doublet" else PALETTE["muted"]
    rest = PALETTE["muted"] if colour_by == "doublet" else PALETTE["ok"]

    for i, s in enumerate(names):
        ax = fig.add_subplot(rows, cols, i + 1)
        entry = embeddings[s] or {}
        xs = [v for v in (entry.get("x") or [])]
        ys = [v for v in (entry.get("y") or [])]
        marks = [v for v in (entry.get(colour_by) or [])]
        if not xs or not ys or len(xs) != len(ys):
            _blank_panel(ax, f"{s}: no embedding")
            continue
        if len(marks) != len(xs):
            # A flag array of the wrong length cannot be matched to a point. Drawing the cloud
            # uncoloured says less than the data holds and says so, rather than aligning two
            # sequences by position and hoping.
            marks = [None] * len(xs)
            _panel_note(ax, "flags do not match the embedding", loc="lower left")
        on, off = [], []
        for x, y, m in zip(xs, ys, marks):
            if _unknown(x) or _unknown(y):
                continue
            (on if _tri(m) is True else off).append((float(x), float(y)))
        if off:
            ax.scatter([p[0] for p in off], [p[1] for p in off], s=1.4, c=rest,
                       alpha=0.45 if colour_by == "doublet" else 0.7, linewidths=0, rasterized=True)
        if on:
            ax.scatter([p[0] for p in on], [p[1] for p in on], s=1.8, c=highlight,
                       alpha=0.9 if colour_by == "doublet" else 0.3, linewidths=0, rasterized=True)
        total = len(on) + len(off)
        share = (100.0 * len(on) / total) if total else float("nan")
        ax.set_title(f"{s}   {share:.1f}%" if total else f"{s}   {NOT_SUPPLIED}",
                     fontsize=7.5, pad=3)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor(PALETTE["muted"])
            spine.set_linewidth(0.5)

    what = ("called a doublet" if colour_by == "doublet"
            else "removed by any criterion")
    _finish(fig, f"{fig_id} - per-library embedding, coloured by whether a nucleus was {what}",
            label_text(subtitle, "one embedding per library; the percentage on each panel is the "
                                 "share highlighted. The same coordinates are used before and "
                                 "after, so a difference between the two is the data"))
    return fig


#: Figure id -> the function that draws it.
#:
#: TWO IDS MAY SHARE A FUNCTION. The id is what the report positions, numbers and captions, so
#: the same data shown two ways is two figures rather than one figure with two panels - a reader
#: referring to "F12" means the linear axis and nothing else. Each entry passes its own arguments,
#: and each function is told which id it is drawing under so its title cannot disagree with the
#: block around it.
FIGURE_FUNCTIONS = {
    "F1": fig_f1_barcode_rank,
    "F2": fig_f2_ambient_removal,
    "F3": fig_f3_cell_calls,
    "F4": fig_f4_scoring_coverage,
    "F5": fig_f5_doublet_sweep,
    "F6": fig_f6_quality_density,
    "F7": fig_f7_before_after,
    "F8": fig_f8_cluster_flags,
    "F9": fig_f9_criterion_contributions,
    "F10": fig_f10_umap_per_library,
    "F11": fig_f10_umap_per_library,
    "F12": fig_f7_before_after,
}
