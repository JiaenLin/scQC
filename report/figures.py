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
    """A panel for something that was never computed. It says so; it does not draw zero.

    WRAPPED HERE rather than by matplotlib. `wrap=True` measures against the FIGURE width, not
    the axes, so on a multi-panel figure a long reason ran straight out of its own dashed box
    and across its neighbours - the one panel whose whole job is to be legible.
    """
    ax.text(0.5, 0.5, "\n".join(_wrapped(line, 34) for line in str(message).split("\n")),
            transform=ax.transAxes, ha="center", va="center", fontsize=7.6,
            color=PALETTE["unknown"])
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
                           fig_id: str = "F6", dpi: int = DEFAULT_DPI):
    """F6 - where is the cut, and why there?

    `densities` maps a sample to `{"x": [...], "y": [...], "n": int|None}` - the density the
    caller estimated, not one recomputed here. `valleys` maps a sample to its own measured
    minimum; `cut` is the single cohort constant, drawn on EVERY panel at the same place so a
    library the constant fits poorly is visible as a line sitting away from that library's own
    valley. `bounds` shades the range outside which a derived floor is refused.

    Alongside F2 this is the figure that carries the report: it is the only one that shows a
    derived threshold sitting where the data actually separates - or not.

    THE X AXIS IS LOG, AND THAT IS NOT A PREFERENCE. The valley is a minimum between two modes
    that are orders of magnitude apart - a few hundred UMI against tens of thousands - and step 5
    finds it by fitting the density in log10 space, which is why the table it writes carries a
    `grid_log10` column beside the grid. Drawn against linear UMI the entire structure this figure
    exists to show, both modes and the cut between them, occupies the first fraction of a percent
    of the axis and the rest is empty: on the calibration cohort the valley sat at 350 and the
    axis ran to 140,000. The ticks read in ORIGINAL UNITS, so the reader is never asked to
    interpret an exponent.

    A density containing a non-positive grid point cannot go on a log axis; that panel is drawn
    linear and SAYS SO, rather than silently dropping the points or the axis.
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
    # ONE DECISION FOR THE WHOLE FIGURE. Panels meant for comparison share limits, so they must
    # also share a scale: two panels of the same quantity on different axes look like a
    # difference in the data.
    log_x = min(xs) > 0 and (_unknown(cut) or float(cut) > 0)

    for i, (s, (x, y, n)) in enumerate(sorted(prepared.items()), start=1):
        ax = fig.add_subplot(rows, cols, i)
        if log_x:
            _log_axis(ax, "x")
        if (not _unknown(bounds) and bounds is not None and len(bounds) == 2
                and not any(_unknown(b) for b in bounds)):
            ax.axvspan(min(xs), float(bounds[0]), color=PALETTE["muted"], alpha=0.12)
            ax.axvspan(float(bounds[1]), max(xs), color=PALETTE["muted"], alpha=0.12)
        ax.plot(x, y, linewidth=1.1, color=PALETTE["ok"])
        # SEPARATED VERTICALLY, not horizontally. The two lines are a measured valley and the
        # cohort constant derived from it, so they sit within a few percent of each other on
        # every well-behaved library - which is the good case, and the case in which two labels
        # at the same height overprint into an unreadable column of letters.
        top = (max(ys) if max(ys) > 0 else 1.0)
        v = valleys.get(s)
        if not _unknown(v):
            ax.axvline(float(v), color=PALETTE["second"], linestyle=":", linewidth=1.1)
            _label(ax, float(v), top * 0.99, f" valley {float(v):,.0f}", fontsize=6.5,
                   rotation=90, va="top", color=PALETTE["second"])
        if not _unknown(cut):
            ax.axvline(float(cut), color=PALETTE["threshold"], linestyle="--", linewidth=1.2)
            _label(ax, float(cut), top * 0.55, f" cut {float(cut):,.0f}", fontsize=6.5,
                   rotation=90, va="top", color=PALETTE["threshold"])
        ax.set_xlim(min(xs), max(xs))
        # A density that is flat at zero everywhere is a real input - it is what a caller passes
        # when the estimate failed - and it must not produce a singular axis warning on top of it.
        ax.set_ylim(0, (max(ys) if max(ys) > 0 else 1.0) * 1.08)
        ax.set_title(str(s), fontsize=8)
        note = n_label(n, "nuclei")
        if _unknown(v):
            note += "\nvalley NOT SUPPLIED"
        if not log_x:
            note += "\nLINEAR axis: a grid point is at or below zero"
        _panel_note(ax, note, loc="upper right")
        ax.set_xlabel(label_text(metric_label, f"metric ({NOT_SUPPLIED} - name it)"),
                      fontsize=7.5)
        ax.tick_params(labelsize=7)
        _tidy(ax)

    _finish(fig, f"{fig_id} · per-library density with the measured valley and the applied cut",
            _wrapped(f"{'log' if log_x else 'linear'} axis, in original units - the density is "
                     f"fitted in log space and the two modes are orders of magnitude apart. The "
                     f"same cut is drawn on every panel at the same scale; a shaded margin is "
                     f"outside the bounds a derived floor is allowed to take", 118))
    return fig


def _violin_row(ax, distributions, names, *, log, cut, cap, metric_name, cut_word="floor"):
    """One panel of paired per-library violins. Returns (positions, n_above_cap).

    PAIRED, AND LABELLED ONCE PER LIBRARY. Each library gets a `before` violin and an `after`
    violin side by side; the axis is labelled with the library, not with four lines of state and
    n under every shape. A reader compares within a pair first and across libraries second, and
    the layout should make the first of those the easy one.

    `cut` IS EITHER ONE NUMBER OR ONE PER LIBRARY, and the drawing differs because the claim
    differs. A cohort constant is a single rule and is drawn as a line crossing every library, so
    a library it fits badly is visible as a line in the wrong place. A per-library threshold is
    ten different rules and CANNOT honestly be drawn that way: one line across all of them would
    assert a cohort constant that was never applied, and the reader has no way to see it is not
    one. It is drawn as a segment over each library's own pair instead.

    `cut_word` is "floor" or "ceiling". The mitochondrial criterion removes cells ABOVE its
    threshold and every count criterion removes cells below theirs, so a figure that says "floor"
    over a ceiling is telling the reader the filter runs the other way.
    """
    import math

    data, positions, colours, ticks, tick_at = [], [], [], [], []
    centre_of = {}
    dropped = above = 0
    for j, s in enumerate(names):
        entry = distributions[s] or {}
        centre = j * 2.2
        centre_of[s] = centre
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

    def _level(v):
        """The y position of a threshold, or None where the axis cannot show it."""
        if _unknown(v):
            return None
        value = float(v)
        if log and value <= 0:
            return None
        return math.log10(value) if log else value

    if isinstance(cut, dict):
        drawn = []
        for s, v in cut.items():
            level = _level(v)
            if level is None or s not in centre_of:
                continue
            c = centre_of[s]
            ax.hlines(level, c - 0.45, c + 1.30, color=PALETTE["refuse"], linestyle="--",
                      linewidth=1.4, zorder=3)
            drawn.append(float(v))
        if drawn:
            span = (f"{min(drawn):,.2f}" if min(drawn) == max(drawn)
                    else f"{min(drawn):,.2f}-{max(drawn):,.2f}")
            _panel_note(ax, f"applied {cut_word}: PER LIBRARY, {span} - each segment is that "
                            f"library's own, and there is no cohort value",
                        loc="lower right", colour=PALETTE["refuse"], wrap=38)
    elif not _unknown(cut):
        level = _level(cut)
        if level is not None:
            ax.axhline(level, color=PALETTE["refuse"], linestyle="--", linewidth=1.2, zorder=3)
            _label(ax, max(positions), level, f" applied {cut_word} = {float(cut):,.0f}",
                   fontsize=6.5, ha="right", va="bottom", color=PALETTE["refuse"])

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
                        cap=None, fig_id: str = "F7", cut_word: str = "floor",
                        dpi: int = DEFAULT_DPI):
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
                                        cap=use_cap, metric_name=metric_name,
                                        cut_word=cut_word)
        if use_cap is not None:
            share = (100.0 * above / total) if total else 0.0
            _panel_note(ax, f"y capped at {float(use_cap):,.0f} - {share:.2f}% of plotted nuclei "
                            f"lie above and are not drawn", loc="upper left")

    axis_words = {"log": "a log axis", "linear": "a linear axis",
                  "both": "log and linear axes"}[scale]
    where = ("each library's own segment - the threshold is PER LIBRARY and there is no cohort "
             "value" if isinstance(cut, dict) else
             f"the applied {cut_word} crosses every library")
    _finish(fig, f"{fig_id} - {metric_name} per library, before and after the cut, on {axis_words}",
            f"grey is every called cell, black is what was retained; the median is drawn and "
            f"printed, and {where}")
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
        # THREE CATEGORIES, NOT TWO. A barcode the detector never examined is not a barcode it
        # cleared, and folding the two together draws every nucleus below the light floor in the
        # same colour as a measured singlet. On a cohort where a quarter of the embedded barcodes
        # sat under that floor, the "not a doublet" cloud would have been a quarter unexamined
        # with nothing on the page saying so.
        on, off, unknown = [], [], []
        for x, y, m in zip(xs, ys, marks):
            if _unknown(x) or _unknown(y):
                continue
            t = _tri(m)
            (on if t is True else (off if t is False else unknown)).append((float(x), float(y)))
        if unknown:
            ax.scatter([p[0] for p in unknown], [p[1] for p in unknown], s=1.4,
                       c=PALETTE["unknown"], alpha=0.5, linewidths=0, rasterized=True)
        if off:
            ax.scatter([p[0] for p in off], [p[1] for p in off], s=1.4, c=rest,
                       alpha=0.45 if colour_by == "doublet" else 0.7, linewidths=0, rasterized=True)
        if on:
            ax.scatter([p[0] for p in on], [p[1] for p in on], s=1.8, c=highlight,
                       alpha=0.9 if colour_by == "doublet" else 0.3, linewidths=0, rasterized=True)
        # The denominator is what was DETERMINED. A share over everything embedded would fall as
        # the never-examined population grows, which reads as the rate improving.
        determined = len(on) + len(off)
        share = (100.0 * len(on) / determined) if determined else float("nan")
        ax.set_title(f"{s}   {share:.1f}%" if determined else f"{s}   {NOT_SUPPLIED}",
                     fontsize=7.5, pad=3)
        if unknown:
            _panel_note(ax, f"{len(unknown):,} never determined", loc="lower left",
                        colour=PALETTE["unknown"], wrap=28)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor(PALETTE["muted"])
            spine.set_linewidth(0.5)

    what = ("called a doublet" if colour_by == "doublet"
            else "removed by any criterion")
    _finish(fig, f"{fig_id} - per-library embedding, coloured by whether a nucleus was {what}",
            label_text(subtitle, "one embedding per library; the percentage on each panel is the "
                                 "share highlighted, over the nuclei for which the flag was "
                                 "DETERMINED. Points in the third colour were never determined "
                                 "and are in no percentage"))
    return fig


# ------------------------------------------------------- F16-F18 · the confounding estimation
#
# THESE THREE ARE NOT A STEP. Every other figure here illustrates something the pipeline DID;
# these illustrate something the EXPERIMENT IS, and no pipeline can change it. If two factors
# partition the libraries identically they stay that way however well the libraries were made,
# and no threshold, no correction and no quantity of data separates them afterwards.
#
# The order is the order the questions have to be asked in:
#
#   F16  WHAT is confounded with what - read off the samplesheet by exact comparison of the
#        partitions. Not a statistic and not a p-value: two factors either induce the same
#        partition of the libraries or they do not, and that is a fact about the design.
#   F17  HOW FAR APART the confounded arms sit on every QC metric, drawn against how far apart
#        the LIBRARIES INSIDE an arm sit. The second is what makes the first readable: an arm
#        difference no larger than the spread between libraries of the same arm is a library
#        effect wearing the design's name.
#   F18  WHETHER THE FILTER made it worse - per criterion, per arm. This is the one part of the
#        confounding the pipeline itself creates, and therefore the one it must show.

#: How two design factors sit against each other.
#:
#: ONLY THE PROBLEM IS INKED. `crossed` - the outcome a reader wants - is left white, so the eye
#: goes to the coloured cells and finds exactly the pairs that make a claim impossible. Filling
#: all three would also have collided with the level colours in the panel beside it, where the
#: same two hues mean two arbitrary levels of some factor and nothing about a verdict.
#:
#: The word is drawn in every cell as well as the colour, so the matrix survives greyscale and a
#: reader who has not read the legend.
RELATION_COLOURS = {"aliased": PALETTE["refuse"], "nested": PALETTE["review"],
                    "crossed": "#FFFFFF"}

#: What each relationship entails for a claim. Printed on the figure: the word on its own does
#: not tell a reader what they may and may not conclude, which is the only thing they need.
RELATION_MEANING = {
    "aliased": "identical partition of the libraries - no analysis of these data can attribute "
               "a difference to one of the two rather than the other",
    "nested": "one is fixed once the other is known - its effect sits inside the other's and "
              "cannot be taken back out of it",
    "crossed": "both levels of each occur with both levels of the other - separable",
}

#: Level colours, one per level within a factor. `unknown` is deliberately NOT among them:
#: PALETTE reserves it for never-examined, and a design level borrowing it would read as an
#: absence of information rather than as a level somebody chose.
LEVEL_COLOURS = ("#0072B2", "#E69F00", "#009E73", "#56B4E9", "#D55E00", "#F0E442")


def _count(v) -> str:
    """A count for an annotation, or a stated absence. Never a silent blank."""
    if _unknown(v):
        return NOT_SUPPLIED
    try:
        return f"{int(v):,}"
    except (TypeError, ValueError, OverflowError):
        return f"{v!r}"


def _ink_on(hex_colour: str) -> str:
    """Black or white text, whichever stays legible on the given fill.

    Relative luminance rather than a per-colour lookup, so editing the palette later cannot
    leave white text on yellow - which is invisible in print and reads as an empty cell.
    """
    h = str(hex_colour).lstrip("#")
    try:
        r, g, b = (int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    except (ValueError, IndexError):
        return "#000000"
    return "#000000" if (0.2126 * r + 0.7152 * g + 0.0722 * b) > 0.55 else "#FFFFFF"


def _level_colours(levels) -> dict:
    """A colour per level, stable for a given ordering of levels."""
    return {lv: LEVEL_COLOURS[i % len(LEVEL_COLOURS)] for i, lv in enumerate(levels)}


def _wrapped(text: str, width: int) -> str:
    lines = []
    for para in str(text).split("\n"):
        lines += textwrap.wrap(para, width) or [""]
    return "\n".join(lines)


def _rect():
    """`matplotlib.patches.Rectangle`, imported where it is used.

    `import matplotlib` does NOT bind `matplotlib.patches`; it is a submodule and reaching it
    through the package attribute works only once something else has imported it. Every figure
    here has to be drawable on its own.
    """
    from matplotlib.patches import Rectangle
    return Rectangle


def fig_f16_design_map(levels, *, samples=None, factors=None, relations=None, arms=None,
                       dpi: int = DEFAULT_DPI):
    """F16 - what is confounded with what, before any number is read?

    `levels` maps a factor to `{sample: level}`; `relations` is a list of
    `{"a", "b", "kind", "detail"}` with `kind` in RELATION_COLOURS; `arms` maps an arm label to
    the libraries in it - an arm being what is left once the aliased factors are collapsed, and
    the unit F17 and F18 compare.

    THE LEFT PANEL IS THE EVIDENCE, THE RIGHT PANEL IS THE READING OF IT. The grid is the
    samplesheet redrawn and nothing is computed to make it; if two columns carry the same block
    structure the two factors are aliased, and that is visible in the grid before the matrix
    beside it says so in a word. A reader who distrusts the matrix can check it by eye.
    """
    Rectangle = _rect()
    if not _mapping(levels):
        raise TaskFailure("F16 needs the design: a mapping of factor -> {sample: level}")
    facs = [f for f in (factors or sorted(levels)) if f in levels]
    if not facs:
        raise TaskFailure("F16 was given a design with no factors to draw")
    names = list(samples or sorted({s for m in levels.values() for s in m}))
    if not names:
        raise TaskFailure("F16 was given no libraries to draw")

    rels = [r for r in (relations or []) if _mapping(r)]
    fig = _new_fig((max(7.6, 1.5 * len(facs) + 5.0), max(3.4, 0.34 * len(names) + 2.4)), dpi)
    # NO `wspace` HERE. `tight_layout` overrides a spacing the gridspec already carries and
    # warns "Axes that are not compatible with tight_layout" on every render when it has to -
    # a warning printed once per report, about a layout that was then silently not what was
    # asked for. The spacing is left to the one thing that is going to decide it anyway.
    gs = fig.add_gridspec(1, 2, width_ratios=(max(1.0, 0.95 * len(facs)),
                                              max(1.6, 1.05 * len(facs))))

    # --- the samplesheet as a grid: one row per library, one column per discovered factor
    ax = fig.add_subplot(gs[0, 0])
    palettes = {f: _level_colours(sorted({str(v) for v in levels[f].values()
                                          if not _unknown(v) and str(v).strip()}))
                for f in facs}
    for j, f in enumerate(facs):
        pal = palettes[f]
        for i, s in enumerate(names):
            raw = levels[f].get(s)
            if _unknown(raw) or not str(raw).strip():
                colour, text, ink, hatch = "#FFFFFF", "not stated", PALETTE["unknown"], \
                    UNKNOWN_HATCH
            else:
                text = str(raw)
                colour, hatch = pal.get(text, PALETTE["unknown"]), None
                ink = _ink_on(colour)
            ax.add_patch(Rectangle((j, len(names) - 1 - i), 1, 1, facecolor=colour,
                                   edgecolor="white", linewidth=1.2, hatch=hatch))
            ax.text(j + 0.5, len(names) - 1 - i + 0.5, text, ha="center", va="center",
                    fontsize=6.4, color=ink, fontweight="bold")
    # The arm bracket lives INSIDE the axes, in a margin the axes reserves for it. Drawn outside
    # with `annotation_clip=False` its label ran under the panel beside it and was cut off
    # mid-word - a label that names the confound, unreadable, on the figure about the confound.
    ax.set_xlim(0, len(facs) + (1.15 if _mapping(arms) else 0.0))
    ax.set_ylim(0, len(names))
    ax.set_xticks([j + 0.5 for j in range(len(facs))])
    ax.set_xticklabels(facs, fontsize=8, rotation=18, ha="right")
    ax.set_yticks([len(names) - 1 - i + 0.5 for i in range(len(names))])
    ax.set_yticklabels(names, fontsize=7)
    ax.tick_params(length=0)
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.set_title(f"the samplesheet, redrawn · {n_label(len(names), 'libraries')}", fontsize=8)

    # --- the arm bracket, on the RIGHT so it cannot collide with the library names. Drawn only
    # where an arm's rows are CONTIGUOUS: a bracket spanning rows that are not adjacent would
    # claim a grouping the picture does not show.
    if _mapping(arms):
        index = {s: i for i, s in enumerate(names)}
        for label, members in arms.items():
            rows = sorted(index[s] for s in members if s in index)
            if not rows or rows != list(range(rows[0], rows[-1] + 1)):
                continue
            top, bottom = len(names) - rows[0], len(names) - 1 - rows[-1]
            ax.annotate("", xy=(len(facs) + 0.16, bottom + 0.06),
                        xytext=(len(facs) + 0.16, top - 0.06),
                        arrowprops={"arrowstyle": "-", "linewidth": 2.6,
                                    "color": PALETTE["refuse"]})
            _label(ax, len(facs) + 0.44, (top + bottom) / 2.0,
                   str(label).replace("\n", " / "), ha="center", va="center", fontsize=6.6,
                   rotation=90, color=PALETTE["refuse"], fontweight="bold")

    # --- the reading: every pair of factors, compared exactly
    ax2 = fig.add_subplot(gs[0, 1])
    kind_of = {}
    for r in rels:
        a, b, kind = str(r.get("a")), str(r.get("b")), str(r.get("kind"))
        kind_of[(a, b)] = kind_of[(b, a)] = kind
    n = len(facs)
    for i, a in enumerate(facs):
        for j, b in enumerate(facs):
            y = n - 1 - i
            if i == j:
                ax2.add_patch(Rectangle((j, y), 1, 1, facecolor="#F2F2F2", edgecolor="white",
                                        linewidth=1.2))
                continue
            kind = kind_of.get((a, b))
            known = kind in RELATION_COLOURS
            colour = RELATION_COLOURS.get(kind, PALETTE["unknown"])
            # A white cell still needs an outline, or the matrix loses its shape where the news
            # is good and a reader cannot tell an empty cell from an absent one.
            ax2.add_patch(Rectangle((j, y), 1, 1, facecolor=colour,
                                    edgecolor="#CCCCCC" if colour == "#FFFFFF" else "white",
                                    linewidth=1.2, hatch=None if known else UNKNOWN_HATCH))
            ax2.text(j + 0.5, y + 0.5, kind if known else NOT_SUPPLIED, ha="center",
                     va="center", fontsize=6.2, color=_ink_on(colour), fontweight="bold")
    ax2.set_xlim(0, n)
    ax2.set_ylim(0, n)
    ax2.set_xticks([j + 0.5 for j in range(n)])
    ax2.set_xticklabels(facs, fontsize=7.5, rotation=18, ha="right")
    ax2.set_yticks([n - 1 - i + 0.5 for i in range(n)])
    ax2.set_yticklabels(facs, fontsize=7.5)
    ax2.tick_params(length=0)
    for side in ("top", "right", "bottom", "left"):
        ax2.spines[side].set_visible(False)
    ax2.set_title("every pair of factors, compared exactly", fontsize=8)

    # THE LEGEND IS THE FIGURE'S POINT, so it is placed on the FIGURE and not inside an axes.
    # Drawn in data coordinates below the matrix it was clipped away entirely by `bbox_inches`,
    # leaving a colour-coded matrix with nothing on the page saying what the colours entail -
    # which is the one thing a reader cannot supply for themselves.
    present = [k for k in ("aliased", "nested", "crossed")
               if any(str(r.get("kind")) == k for r in rels)]
    n_aliased = sum(1 for r in rels if str(r.get("kind")) == "aliased")
    _finish(fig, "F16 - what is confounded with what",
            _wrapped((f"{n_aliased} pair{'' if n_aliased == 1 else 's'} of design factors "
                      f"partition{'s' if n_aliased == 1 else ''} the libraries identically"
                      if n_aliased else
                      "no pair of design factors partitions the libraries identically")
                     + " - computed from the samplesheet by exact comparison, not estimated", 110))
    if present:
        # ONE LINE PER RELATIONSHIP. Joined with a separator and wrapped, the wrap fell on the
        # separator and left a dangling "·" at the end of a line, which reads as a truncation.
        text = "\n".join(_wrapped(f"{k.upper()}: {RELATION_MEANING[k]}", 150) for k in present)
        # Reserved in the layout and measured in inches, for the reasons F17 records at length.
        band = (0.20 * (1 + text.count("\n")) + 0.10) / fig.get_figheight()
        fig.tight_layout(rect=(0, min(band, 0.3), 1, 0.93))
        fig.text(0.5, 0.008, text, ha="center", va="bottom", fontsize=6.6,
                 color=PALETTE["muted"])
    return fig


#: The marker a library's own median is drawn with. TEN LIBRARIES IS COMMON AND TEN
#: DISTINGUISHABLE HUES ARE NOT, so identity is carried by colour AND shape together: the colour
#: cycles every five libraries and the shape changes when it wraps, which keeps every pair unique
#: up to thirty. Yellow is excluded - a 4-point yellow marker on white is not a mark.
MEDIAN_COLOURS = ("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9")
MEDIAN_MARKERS = ("o", "s", "^", "D", "v", "P")


def _median_style(i: int) -> tuple:
    return MEDIAN_COLOURS[i % len(MEDIAN_COLOURS)], MEDIAN_MARKERS[
        (i // len(MEDIAN_COLOURS)) % len(MEDIAN_MARKERS)]


def fig_f17_metrics_by_arm(metrics, *, arms=None, population=None, cap=None, absent=None,
                           dpi: int = DEFAULT_DPI):
    """F17 - how far apart do the confounded arms sit, and is that further than their libraries?

    `metrics` is a list of `{"key", "label", "log", "by_arm": {arm: [values]},
    "library_medians": {arm: {sample: median}}}`. One panel per metric, one violin per arm, and
    every library's OWN median drawn on top of the arm it belongs to.

    WHY THE LIBRARY MEDIANS ARE THE POINT OF THIS FIGURE

    A violin per arm shows a difference. It cannot show whether that difference is the arm or the
    libraries that happen to be in it - and in a confounded design those are the same libraries
    every time, which is exactly why the design is confounded. Drawing each library's own median
    inside its arm puts the two spreads on one axis: if the libraries of one arm sit as far from
    each other as the arms sit from each other, the arm difference is not evidence about the
    factor, whatever a test of it would report.

    So every panel states two quantities in its own heading - the GAP between the arms' median
    library medians, and the widest WITHIN-ARM span of them. Their ratio is this panel's estimate
    of how much of the metric the design could even in principle account for. IT IS AN ESTIMATE
    AND THE FIGURE SAYS SO. There is no test here and no p-value: a confounded factor cannot be
    tested, which is the whole reason this section exists.

    THE NUMBERS ARE IN THE HEADING AND THE LIBRARIES ARE IN ONE LEGEND, and both were learned by
    drawing it the other way. Printed inside the panel the numbers sat on top of the violin;
    written beside each dot the library names overplotted into an unreadable smear the moment two
    libraries had similar medians - which, in the case this figure exists for, they do. Every
    library keeps the same colour, the same marker and the same horizontal offset in EVERY panel,
    so a library that is an outlier on one metric can be found on the others.

    `absent` maps a metric label to the reason it could not be drawn, and those panels are drawn
    as stated absences - a panel silently left out reads as a metric nobody thought worth showing.
    """
    import math

    Rectangle = _rect()
    specs = [m for m in (metrics or []) if _mapping(m)]
    gone = _mapping(absent)
    if not specs and not gone:
        raise TaskFailure("F17 needs at least one metric, or a stated reason for every absence")
    order = list(arms or sorted({a for m in specs for a in _mapping(m.get("by_arm"))}))
    if not order:
        raise TaskFailure("F17 needs the confounded arms to compare")
    arm_x = {a: float(i) for i, a in enumerate(order)}

    # Every library, in arm order then name order, with its style and its offset fixed ONCE for
    # the whole figure. A library that moved between panels could not be compared across them.
    libraries, offset_of, style_of = [], {}, {}
    for arm in order:
        here = sorted({s for m in specs
                       for s in _mapping(_mapping(m.get("library_medians")).get(arm))})
        for j, s in enumerate(here):
            offset_of[s] = (arm_x[arm] - 0.26 + (0.52 * j / max(len(here) - 1, 1))
                            if len(here) > 1 else arm_x[arm])
            style_of[s] = _median_style(len(libraries))
            libraries.append(s)

    # NOT `grid_shape`, which is tuned for per-library panels and puts seven into a 3x3 with two
    # dead cells. Seven QC metrics is the ordinary case here, so the width is chosen to leave at
    # most one empty cell while staying wide rather than tall.
    total = len(specs) + len(gone)
    cols = total if total <= 4 else (3 if total in (5, 6, 9) else 4)
    rows = -(-total // cols)
    fig = _new_fig((3.6 * cols + 0.5, 3.15 * rows + 1.5), dpi)

    for idx, spec in enumerate(specs):
        ax = fig.add_subplot(rows, cols, idx + 1)
        log = bool(spec.get("log"))
        label = label_text(spec.get("label"), f"metric ({NOT_SUPPLIED} - name it)")
        by_arm = _mapping(spec.get("by_arm"))
        med_by_arm = _mapping(spec.get("library_medians"))

        data, positions, dropped = [], [], 0
        for arm in order:
            vals = [float(v) for v in (by_arm.get(arm) or []) if not _unknown(v)]
            if log:
                positive = [math.log10(v) for v in vals if v > 0]
                dropped += len(vals) - len(positive)
                vals = positive
            if not vals:
                continue
            data.append(vals)
            positions.append(arm_x[arm])
        if not data:
            _blank_panel(ax, f"{label}\n\nno values in any arm")
            continue

        parts = ax.violinplot(data, positions=positions, widths=0.66, showextrema=False,
                              showmedians=False)
        for body in parts["bodies"]:
            body.set_facecolor(PALETTE["before"])
            body.set_alpha(0.95)
            body.set_edgecolor("none")

        centres, spans = {}, {}
        for arm in order:
            meds = {s: float(v) for s, v in _mapping(med_by_arm.get(arm)).items()
                    if not _unknown(v) and (not log or float(v) > 0)}
            if not meds:
                continue
            y = {s: (math.log10(v) if log else v) for s, v in meds.items()}
            ordered = sorted(y.values())
            centres[arm] = ordered[len(ordered) // 2]
            spans[arm] = (ordered[0], ordered[-1])
            ax.add_patch(Rectangle((arm_x[arm] - 0.33, spans[arm][0]), 0.66,
                                   max(spans[arm][1] - spans[arm][0], 0.0),
                                   facecolor=PALETTE["review"], alpha=0.30,
                                   edgecolor=PALETTE["review"], linewidth=0.9, zorder=3))
            for s, v in y.items():
                colour, marker = style_of.get(s, (PALETTE["after"], "o"))
                ax.plot([offset_of.get(s, arm_x[arm])], [v], marker=marker, markersize=4.6,
                        color=colour, markeredgecolor="white", markeredgewidth=0.6,
                        linestyle="none", zorder=5)
            ax.hlines(centres[arm], arm_x[arm] - 0.35, arm_x[arm] + 0.35, color=PALETTE["ok"],
                      linewidth=2.4, zorder=6)

        # --- the two quantities, in the heading where nothing overlaps them
        if len(centres) >= 2:
            cs = [centres[a] for a in order if a in centres]
            gap = max(cs) - min(cs)
            widest = max((spans[a][1] - spans[a][0]) for a in spans)
            if log:
                shown_gap, shown_span = f"{10 ** gap:,.2f}x", f"{10 ** widest:,.2f}x"
            else:
                shown_gap, shown_span = f"{gap:,.3g}", f"{widest:,.3g}"
            head = (f"between arms {shown_gap}   ·   widest within an arm {shown_span}\n"
                    + ("one library in every arm - ratio not computable" if widest <= 0
                       else f"ratio {gap / widest:,.2f}"))
        elif centres:
            head = "only one arm has library medians\nthere is nothing to compare it with"
        else:
            head = "no library median in any arm"
        n_drawn = sum(len(v) for v in data)
        ax.set_title(f"{head}   ·   {n_label(n_drawn, 'nuclei drawn')}", fontsize=6.6,
                     color=PALETTE["muted"])

        ax.set_xticks([arm_x[a] for a in order])
        ax.set_xticklabels([str(a) for a in order], fontsize=6.8)
        ax.set_xlim(min(arm_x.values()) - 0.55, max(arm_x.values()) + 0.55)
        ax.set_ylabel(label, fontsize=7.5)
        ax.margins(y=0.12)
        if log:
            from matplotlib.ticker import FuncFormatter
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{10 ** v:,.0f}"))
            if dropped:
                _panel_note(ax, f"{dropped:,} value(s) at or below zero are not on this log "
                                f"axis", loc="lower right", wrap=28)
        _tidy(ax)

    for k, (label, reason) in enumerate(sorted(gone.items())):
        ax = fig.add_subplot(rows, cols, len(specs) + k + 1)
        _blank_panel(ax, f"{label}\n\n{reason}")

    pop = label_text(population, f"population {NOT_SUPPLIED}")
    cap_words = ("" if _unknown(cap) else
                 f"; at most {int(cap):,} nuclei per arm, at an even stride")
    head = ("F17 - every QC metric across the confounded arms",
            _wrapped(f"grey violin: every nucleus of the arm. marker: one library's OWN median, "
                     f"in the same colour, shape and position in every panel. amber band: the "
                     f"span those medians cover. blue rule: the arm's median of them. Over "
                     f"{pop}{cap_words}", 132))

    # ONE LEGEND FOR THE WHOLE FIGURE, at the foot, because a marker means the same thing in
    # every panel; per-panel legends would repeat it up to nine times and cover the data.
    #
    # THE BAND IT SITS IN IS RESERVED IN THE LAYOUT, not taken back afterwards. Calling
    # `subplots_adjust` after `tight_layout` moves the bottom edge and leaves the rows to absorb
    # the difference unevenly - which on a grid whose last row holds one panel gave that panel a
    # third of the figure, and it was the STATED-ABSENCE panel, drawn as a dashed box the size of
    # the other six put together.
    handles = []
    if libraries:
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], linestyle="none", marker=style_of[s][1], markersize=4.6,
                          color=style_of[s][0], markeredgecolor="white", markeredgewidth=0.6,
                          label=str(s)) for s in libraries]
    # RESERVED IN INCHES, converted to the fraction `rect` wants. A fraction is height-dependent:
    # the same 0.028-per-line rule that fits a two-row grid leaves a hand's width of nothing above
    # a four-row one, because a title does not get taller when the figure does.
    height = fig.get_figheight()
    legend_rows = -(-len(handles) // 6) if handles else 0
    bottom = (0.20 * legend_rows + 0.10) / height if legend_rows else 0.0
    top = 1.0 - (0.20 * (1 + head[1].count("\n")) + 0.10) / height
    fig.suptitle(f"{head[0]}\n{head[1]}", fontsize=9, y=0.995)
    fig.tight_layout(rect=(0, min(bottom, 0.3), 1, max(top, 0.6)))
    if handles:
        fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 6),
                   fontsize=6.6, frameon=False, handletextpad=0.4, columnspacing=1.4)
    return fig


def fig_f18_removal_by_arm(rates, *, arms=None, overall=None, criteria=None,
                           dpi: int = DEFAULT_DPI):
    """F18 - did the filter remove the same share from each confounded arm?

    `rates` maps an arm to `{criterion: percent of that arm's barcodes the criterion fired on}`;
    `overall` maps an arm to `{"n_in", "n_removed", "pct"}`.

    A REMOVAL RATE THAT DIFFERS ACROSS THE ARMS OF A CONFOUNDED DESIGN IS THE ONE PART OF THE
    CONFOUNDING THE PIPELINE ITSELF CREATES. Whatever the libraries were on arrival, if the
    filter took twice as large a share out of one arm then the arms now differ for a reason with
    nothing to do with the factor they are named after - and downstream neither sees it nor can
    undo it.

    The ratio between the extreme arms is printed on the left panel, because it is the number a
    reader has come here for.
    """
    by_arm = _mapping(rates)
    if not by_arm:
        raise TaskFailure("F18 needs a per-arm removal rate for at least one arm")
    order = list(arms or sorted(by_arm))
    crits = list(criteria or sorted({c for m in by_arm.values() for c in _mapping(m)}))
    tot = _mapping(overall)

    fig = _new_fig((max(8.0, 0.62 * len(crits) * max(1, len(order)) + 4.6), 4.0), dpi)
    gs = fig.add_gridspec(1, 2, width_ratios=(1.0, 2.4))      # spacing: see F16

    # --- left: the whole filter, per arm
    ax = fig.add_subplot(gs[0, 0])
    drawn = []
    for i, a in enumerate(order):
        entry = _mapping(tot.get(a))
        if not _unknown(entry.get("pct")):
            drawn.append((i, float(entry["pct"])))
    if drawn:
        ax.bar([i for i, _ in drawn], [p for _, p in drawn], width=0.58,
               color=PALETTE["after"], zorder=3)
        for i, p in drawn:
            entry = _mapping(tot.get(order[i]))
            ax.annotate(f"{p:.2f}%\n{_count(entry.get('n_removed'))} of "
                        f"{_count(entry.get('n_in'))}", (i, p), textcoords="offset points",
                        xytext=(0, 4), ha="center", fontsize=6.4, fontweight="bold")
        vals = [p for _, p in drawn]
        lo, hi = min(vals), max(vals)
        ratio = ("not computable - an arm removed nothing" if lo <= 0 else f"{hi / lo:,.2f}x")
        # IN THE HEADING, not on the panel. Inside the axes it sat across the top of the tallest
        # bar, which is the arm the reader is being asked to look at.
        ax.set_title(f"the filter as a whole\nwidest ratio across the arms: {ratio}",
                     fontsize=7.4, color=PALETTE["refuse"])
    else:
        _blank_panel(ax, "no per-arm removal total was supplied")
        ax.set_title("the filter as a whole", fontsize=8)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([str(a) for a in order], fontsize=6.4)
    ax.set_ylabel("removed (% of the arm's barcodes)", fontsize=7.5)
    ax.margins(y=0.20)
    _tidy(ax)

    # --- right: each criterion on its own
    ax2 = fig.add_subplot(gs[0, 1])
    if not crits:
        _blank_panel(ax2, "no criterion rates were supplied")
    else:
        width = 0.8 / max(len(order), 1)
        pal = _level_colours(order)
        for k, a in enumerate(order):
            row = _mapping(by_arm.get(a))
            xs = [i - 0.4 + width * (k + 0.5) for i in range(len(crits))]
            ys = [(0.0 if _unknown(row.get(c)) else float(row[c])) for c in crits]
            ax2.bar(xs, ys, width=width * 0.9, label=str(a).replace("\n", " / "), color=pal[a],
                    zorder=3, edgecolor="white", linewidth=0.4)
            # A criterion the arm has no value for is NOT a criterion that fired on nothing.
            # It is named in the never-examined colour rather than drawn as a zero-height bar,
            # which is indistinguishable from a criterion that fired on nobody.
            for i, c in enumerate(crits):
                if _unknown(row.get(c)):
                    _annot(ax2, "no value", (xs[i], 0.0), textcoords="offset points",
                           xytext=(0, 3), ha="center", fontsize=5.4,
                           color=PALETTE["unknown"], rotation=90)
        # THE RATIO PER CRITERION, which is the question the grouped bars are asked. Eyeballing
        # two bar heights answers "which is taller"; the number answers "by how much", and it is
        # the same arithmetic the whole section turns on.
        top = max([0.0] + [float(v) for a in order for v in _mapping(by_arm.get(a)).values()
                           if not _unknown(v)])
        for i, c in enumerate(crits):
            vals = [float(_mapping(by_arm.get(a))[c]) for a in order
                    if not _unknown(_mapping(by_arm.get(a)).get(c))]
            if len(vals) < 2:
                continue
            lo, hi = min(vals), max(vals)
            txt = "-" if lo <= 0 else f"{hi / lo:,.2f}x"
            _annot(ax2, txt, (i, max(vals)), textcoords="offset points", xytext=(0, 5),
                   ha="center", fontsize=6.2, fontweight="bold",
                   color=PALETTE["muted"] if lo > 0 else PALETTE["unknown"])
        ax2.set_ylim(0, top * 1.16 if top > 0 else 1.0)
        ax2.set_xticks(range(len(crits)))
        ax2.set_xticklabels([str(c).replace("fail_", "").replace("_", " ") for c in crits],
                            fontsize=6.8, rotation=22, ha="right")
        ax2.set_ylabel("fired on (% of the arm's barcodes)", fontsize=7.5)
        ax2.legend(fontsize=6.0, frameon=False, ncol=min(len(order), 2), loc="upper right")
        ax2.set_title("each criterion on its own - criteria overlap, so these do not sum",
                      fontsize=8)
        _tidy(ax2)

    _finish(fig, "F18 - what the filter removed from each confounded arm",
            _wrapped("a share that differs across the arms is a technical removal that will read "
                     "as the design downstream; the number over each pair is that ratio, and the "
                     "denominator throughout is every barcode the arm held", 118))
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
    # The gene axis gets its own id for the same reason F12 does. Step 5 derives TWO count floors
    # from two densities and applies both, and one F6 can carry only one metric: `cut` is a single
    # cohort constant drawn on every panel, and the UMI floor drawn over a gene density would be a
    # line in the wrong units on ten panels. Folding both into one figure under one id would also
    # make "the density figure" ambiguous in a report whose whole point is that a reader can name
    # what they looked at.
    "F13": fig_f6_quality_density,
    # The other two applied axes, same function as F7 for the same reason F12 shares it: the
    # three criteria step 7 applies should be read in one form, so a difference between the
    # panels is the filter and not the chart.
    "F14": fig_f7_before_after,
    "F15": fig_f7_before_after,
    # THE CONFOUNDING BLOCK. Not attached to a step: these three describe the DESIGN and
    # what the filter did to it, neither of which belongs to one stage of the pipeline.
    "F16": fig_f16_design_map,
    "F17": fig_f17_metrics_by_arm,
    "F18": fig_f18_removal_by_arm,
}
