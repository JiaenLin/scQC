# Assembles the QC report. This module decides nothing about the data: it reads a payload the
# orchestrator hands it, renders exactly what is in it, and states - visibly, on the page - every
# place the payload said nothing. It removes nothing, recomputes no threshold and never fills a
# gap with a plausible value.
"""The report of docs/REPORT_DESIGN.md: one self-contained HTML file and its JSON companion.

WHAT THIS DOCUMENT IS FOR, AND WHY THAT CHANGES ITS SHAPE

A metrics report answers what happened. A QC report has to answer what was DECIDED, on what
evidence, and what is still unresolved - because quality control is not a measurement, it is a
sequence of choices about what to discard. So the verdict is section 1 rather than an appendix,
a threshold appears with who chose it, and every step carries what it CANNOT establish.

THREE PROPERTIES THAT ARE NOT COSMETIC

  A partial run still produces a report.  A run that stopped at step 2 is exactly the run whose
     report matters most, because the document is the only record of why it stopped. Refusing to
     write it, or writing it with the unreached steps quietly omitted, destroys that record. Steps
     with no entry in the payload are rendered as NO RECORD SUPPLIED, never as steps that passed.

  Absent and empty are different claims.  `open_items: []` says nothing is open. No `open_items`
     key at all says nobody looked. The first renders `none`; the second renders a defect. The
     same distinction is applied to gate findings, to the parameter table and to every section,
     and it is the whole reason this module carries a MISSING sentinel rather than using
     `dict.get(k, default)`: a default here would turn "not recorded" into "recorded as fine",
     which is the failure docs/PRINCIPLES.md section 4 is about.

  The report reports on itself.  Anything the payload should have carried and did not becomes a
     DEFECT, counted on the front page and listed in the JSON. A report that silently omits a
     required block reads exactly like a complete one.

THE PAYLOAD

Every key is optional; every absence is stated on the page. Nothing below has a default value.

    {
      "run":        {"project": str, "mode": "evidence"|"apply", "invocation": str,
                     "started": iso},
      "deliverable":{"text": str,                # preferred: the one line, already composed
                     "n_kept": int, "n_in": int, # or these two, from which it is composed
                     "unit": str,                # what is being counted, e.g. "called cells"
                     "stopped_after": step_key, "stopped_because": str},
      "gates":      [{"step": step_key, "check": str, "severity": "REFUSE"|"REVIEW"|"ok",
                      "message": str, "detail": [str, ...]}, ...],   # or the modules' Finding
                                                                     # objects, read by attribute
      "parameters": [{"name": str, "value": any, "class": "FIXED"|"DERIVED"|"DECLARED"|
                      "ADJUDICATED", "basis": str, "verbatim": str, "decided_by": str,
                      "decided_on": str}, ...],
      "steps":      [{"key": step_key, "status": str, "what_it_does": str,
                      "found": [{"label": str, "value": any, "source": path}, ...],
                      "cannot_establish": str, "figures": ["F3", ...],
                      "sources": [path, ...]}, ...],
      "figures":    {"F1": {"data": {...}      # keyword arguments for report.figures.fig_f1_*
                            | "png_base64": str,   # or an image the caller rendered itself
                            "caption": str, "source": path}, ...},
      "provenance": {"pipeline": {"version": str, "commit": str, "dirty": bool,
                                  "describe": str, "branch": str},
                     "decisions": {"path": path, "hash": str},
                     "reference": {"registry": str, "genome": str, "annotation": str,
                                   "introns": str},
                     "tools": {"cellbender": "0.3.2", ...},   # OBSERVED strings only
                     "tools_expected": ["cellbender", "celescope", ...],
                     "input_check": [{"name": str, "p1_raw_values": bool,
                                      "p2_raw_droplets": bool, "reasons": [str, ...]}, ...],
                     "environment": {...},
                     "generated": iso, "newest_input": iso,
                     "inputs": {label: {"path": path, "mtime": iso}}},
      "open_items": [{"item": str, "closes_when": str, "blocked_on": str}, ...],
      "per_sample": {"source": path,       # every threshold, one row per library
                     "columns": [{"key": str, "label": str,
                                  "scope": "per library"|"cohort constant", "step": str}, ...],
                     "rows":    [{"sample": str, <key>: any, ...}, ...]},
    }

FRESHNESS IS COMPUTED HERE AND ENFORCED ELSEWHERE

`generated` earlier than `newest input` means the document describes data that has since changed,
and staleness has no symptom - a stale report opens, renders and reads exactly like a current one.
So it is computed, printed on the face of the document and recorded in the JSON as `stale`. It
does not stop the write: refusing to write the artifact would be refusing the remedy, and the
rebuilt document is the thing that fixes the staleness. `refuse_if_stale()` is the enforcing
form, for a caller that is about to PUBLISH one.

A timestamp that cannot be read, or a pair that cannot be compared because only one of them
carries a UTC offset, is recorded as NOT CHECKED with the reason - never as fresh.
"""

from __future__ import annotations

import base64
import html
import io
import json
import sys
import time
from datetime import datetime
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from engine.task import Refusal, TaskFailure  # noqa: E402

REPORT_VERSION = 1


class _MissingType:
    """A key that was not in the payload. It is not None, not empty, and not falsy.

    Truth-testing it raises rather than answering, because every bug this sentinel exists to
    prevent has the shape `if payload.get("gates"):` - which reads an absent key as "nothing to
    report" and prints a clean page over an unexamined run.
    """

    __slots__ = ()

    def __repr__(self) -> str:                                            # pragma: no cover
        return "MISSING"

    def __bool__(self):
        raise TypeError("MISSING must be compared with `is`, never truth-tested: an absent key "
                        "is not an empty one and the two are different claims")


MISSING = _MissingType()

NOT_STATED = "NOT STATED"
NOT_INVOKED = "not invoked"

#: The pipeline's own step list, with the one-line purpose and the removal column taken verbatim
#: from the steps table in README.md. These are properties of the pipeline, not of a run, so they
#: are constants here; anything a RUN determined must come from the payload. Where a step section
#: falls back to this text the JSON records `purpose_source: "README steps table"`, so a reader
#: can always tell the pipeline's description of a step from this run's.
STEPS = (
    ("00_ingest", "0 · ingest",
     "validate samplesheet, resolve reference, verify input is raw", "nothing", ("F1",)),
    ("01_ambient", "1 · ambient",
     "denoise (CellBender); audit the result; halve the learning rate if degenerate",
     "nothing", ("F2",)),
    ("02_cells", "2 · cell call",
     "compare aligner and denoiser calls; gate the loss", "nothing", ("F3",)),
    ("03_light_floor", "3 · light floor",
     "technical floor for doublet scoring - not a quality filter", "nothing", ("F4",)),
    ("04_doublets", "4 · doublets",
     "score per sample, before quality filtering; flag only", "nothing", ("F5",)),
    ("05_quality", "5 · quality",
     "derive count floors and mitochondrial ceiling", "nothing", ("F6", "F7")),
    ("06_cluster_check", "6 · cluster check",
     "per-cluster flags: depth, mitochondrial, markers, doublet", "nothing", ("F8",)),
    ("07_apply", "7 · apply", "pre-flight, verify approval, remove",
     "YES - the only step in the pipeline that removes anything", ("F9",)),
)

#: The question each figure answers, from docs/REPORT_DESIGN.md. Printed with the figure so a
#: reader knows what it is for before deciding whether it answers them.
FIGURE_QUESTIONS = {
    "F1": "is this really raw, unfiltered input?",
    "F2": "how much was removed, from what, and evenly?",
    "F3": "did the denoiser drop cells the aligner kept?",
    "F4": "what was never examined?",
    "F5": "is the rate a measurement or the prior?",
    "F6": "where is the cut and why there?",
    "F7": "what did the cut change?",
    "F8": "are any clusters technical?",
    "F9": "what did each criterion remove uniquely?",
    "F10": "where in the manifold did the doublets sit?",
    "F11": "did the removed nuclei leave as a population, or scattered?",
    "F12": "the same count distributions, on the scale people work in",
}

PARAM_CLASSES = {
    "FIXED": "true regardless of dataset; changing it is a code change",
    "DERIVED": "computed from this dataset; reproducible from the data alone",
    "DECLARED": "supplied by the operator up front, before seeing the result",
    "ADJUDICATED": "a human decided, after seeing the result, in their own words",
}

CANNOT_ESTABLISH_MISSING = ("cannot establish: NOT STATED — this is a defect in the report, not "
                            "an absence of limits.")

_SEVERITY_ALIASES = {"REFUSE": "REFUSE", "REVIEW": "REVIEW", "OK": "ok", "PASS": "ok",
                     "NOTE": "ok", "INFO": "ok"}


# ------------------------------------------------------------------------------------ helpers


#: Modules consulted for their missing-value singletons, cached by name so a package that is not
#: installed is searched for once rather than on every value. A cached module object stays valid,
#: and a cached `None` stays valid too: a package that cannot be imported cannot be the source of
#: a sentinel later in the same process.
_OPTIONAL_MODULES: dict = {}


def _optional_module(name: str):
    """An already-imported module, or one imported now, or None. Never raises.

    `sys.modules` is consulted first because it answers without paying for an import - and it is
    sufficient on its own for the sentinels below, since an instance of `pandas.NA` cannot exist
    in a process that has not imported pandas. The import is the fallback and is tolerated as
    failing, so this module stays importable with no third-party package installed - which the
    CLI depends on.
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

    Compared by IDENTITY at the call site, never by equality: `pd.NA == pd.NA` is `pd.NA` and
    truth-testing that raises TypeError, so an equality test here would turn a blank cell into a
    crash in the middle of assembling the document.
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

    A payload assembled from real tables carries more spellings of "no value" than None and a
    Python NaN, and one that slips through is rendered into the document as its repr - `<NA>`,
    `NaT`, `--` - which reads as a value somebody recorded. Covered:

      MISSING                   the key was not in the payload at all (this module's sentinel)
      None                      the key was there and said nothing
      float('nan')              a blank cell
      numpy.float64('nan')      a numpy-backed column; a float subclass, so `isinstance` catches
                                it - `is None` never would
      numpy.float32('nan'),
      numpy.datetime64('NaT')   not float subclasses - caught by the self-inequality test
      pandas.NA, pandas.NaT     neither is None, neither is a float; identity only
      numpy.ma.masked           a masked element read out of a masked array

    An object whose `!=` cannot be reduced to one boolean - a numpy ARRAY, a list - is NOT
    unknown. It is a container, and calling it unknown would drop a whole block from the report.
    """
    if v is None or v is MISSING:
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


def _tristate(value):
    """True, False, or None for `not a determination`. The only way a flag is read here.

    Two failures, one predicate. `value is True` is False for `numpy.bool_(True)`, so a check
    that genuinely PASSED was printed as NOT CHECKED; and a bare truth test in its place reads
    any non-empty string - `"false"` included - as a pass. So: an unknown value is None, a
    boolean or a number is taken as it is, and a string is mapped only from the spellings that
    unambiguously mean yes or no. Anything else is None, because an unrecognised value is not a
    passing one and must never be rendered as PASS.
    """
    if _unknown(value):
        return None
    if isinstance(value, str):
        s = value.strip().lower()
        if s in ("true", "yes", "y", "pass", "passed", "ok", "1"):
            return True
        if s in ("false", "no", "n", "fail", "failed", "0"):
            return False
        return None
    try:
        return bool(value)
    except Exception:                                                     # noqa: BLE001
        return None


def _as_list(value) -> list:
    """A list-shaped payload field, whatever shape it actually arrived in.

    `[str(d) for d in detail]` over a payload that carried a single string instead of a list of
    strings iterates it CHARACTER BY CHARACTER and renders one bullet per letter into the
    document - which looks like a defect in the data rather than in the reader. Scalars are
    therefore normalised to a one-element list, and that includes a dict: a caller who supplied
    one finding instead of a list of one finding meant one finding, not its keys.

    MISSING and unknown become the EMPTY list. Callers that must distinguish "absent" from
    "empty" - and several here do - test for MISSING themselves before calling this.
    """
    if _unknown(value):
        return []
    if isinstance(value, (str, bytes, bytearray, dict)):
        return [value]
    try:
        return list(value)
    except TypeError:
        return [value]


def _as_number(value):
    """A finite float from whatever the payload carried, or None if it is not one.

    Returns None rather than raising, so a caller can state the defect on the page. `float()` and
    `int()` on an unparseable value raise a bare ValueError out of the assembly, which destroys
    the report - and the report is the only record of why the run stopped, so losing it to a
    malformed count is the worst possible trade.
    """
    if _unknown(value):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    import math
    return out if math.isfinite(out) else None


def _count_text(x: float) -> str:
    """A count for the deliverable line: thousands-separated, and never silently truncated.

    A non-integral count is printed as it arrived rather than through `int()`, which would turn
    12.7 into 12 on the one line of the document everybody reads.
    """
    return f"{int(x):,}" if float(x).is_integer() else f"{x:,}"


def _get(mapping, key):
    """The value at `key`, or MISSING. Never a default, never a coerced empty."""
    if not isinstance(mapping, dict) or key not in mapping:
        return MISSING
    return mapping[key]


def _stated(value) -> bool:
    """True when a value was supplied AND says something.

    MISSING, every spelling of unknown `_unknown` knows about, and a blank string are not. The
    MISSING test is kept explicit even though `_unknown` covers it, because this is the predicate
    every section calls and the sentinel it must never let through is worth naming here.
    """
    if value is MISSING or _unknown(value):
        return False
    return bool(str(value).strip())


def _text(value, absent: str = NOT_STATED) -> str:
    return str(value) if _stated(value) else absent


def _defect(defects: list, where: str, what: str, severity: str = "defect") -> None:
    defects.append({"where": where, "what": what, "severity": severity})


def _attr(obj, name):
    """A field from a dict OR from one of the modules' Finding dataclasses. MISSING if absent."""
    if isinstance(obj, dict):
        return _get(obj, name)
    return getattr(obj, name) if hasattr(obj, name) else MISSING


def _parse_iso(value):
    """(epoch_seconds, has_offset) from an ISO-8601 timestamp, or (None, None) if unreadable.

    Tolerant about the spellings that reach a report - a trailing Z, an offset written without a
    colon, a space where the T should be - because refusing a legible timestamp on punctuation
    would turn a freshness check into a nuisance, and a nuisance check gets switched off. What it
    will not do is guess: an unreadable value returns None and the comparison is recorded as NOT
    CHECKED rather than as fresh.
    """
    if not _stated(value):
        return None, None
    s = str(value).strip()
    if "T" not in s and " " in s:
        s = s.replace(" ", "T", 1)
    if s.endswith(("Z", "z")):
        s = s[:-1] + "+00:00"
    if len(s) >= 5 and s[-5] in "+-" and s[-3] != ":" and s[-5:].strip("+-").isdigit():
        s = s[:-2] + ":" + s[-2:]
    try:
        dt = datetime.fromisoformat(s)
    except (ValueError, TypeError):
        return None, None
    return dt.timestamp(), dt.tzinfo is not None


def freshness(generated, newest_input) -> dict:
    """Is this artifact older than the data it describes?

    Returns `stale` as True, False or None - never False by default. None is NOT CHECKED and
    carries the reason, because a freshness check that reports "fresh" when it could not run is
    worse than no check: it is the only kind of wrongness that looks exactly like correctness.
    """
    out = {"generated": _text(generated), "newest_input": _text(newest_input),
           "checked": False, "stale": None, "reason": "", "margin_seconds": None}
    if not _stated(generated) or not _stated(newest_input):
        out["reason"] = ("NOT CHECKED - "
                         + ("no generation time; " if not _stated(generated) else "")
                         + ("no newest-input time was supplied, so nothing was compared"
                            if not _stated(newest_input) else ""))
        return out
    g, g_tz = _parse_iso(generated)
    n, n_tz = _parse_iso(newest_input)
    if g is None or n is None:
        out["reason"] = ("NOT CHECKED - "
                         f"{'the generation time' if g is None else 'the newest-input time'} "
                         f"could not be read as ISO-8601")
        return out
    if g_tz != n_tz:
        out["reason"] = ("NOT CHECKED - one timestamp carries a UTC offset and the other does "
                         "not; comparing them would assume a zone neither of them states")
        return out
    out["checked"] = True
    out["margin_seconds"] = g - n
    out["stale"] = g < n
    out["reason"] = ("the artifact is OLDER than its newest input" if out["stale"]
                     else "the artifact is at least as new as its newest input")
    return out


def newest_input_time(paths) -> dict:
    """The newest modification time among files that exist, and the names of those that do not.

    Offered for a caller that has paths rather than timestamps. Absent files are LISTED, not
    skipped: a newest-input time computed over the three inputs that happened to exist describes
    a different run from the one that had five.

    The timestamp carries a UTC offset, matching what `engine.provenance.environment()` writes
    for `generated`. Both sides of the freshness comparison must be spelled the same way or it
    reports NOT CHECKED - correct, and needlessly so if the two halves of one pipeline disagree
    about punctuation.
    """
    given = [] if _unknown(paths) else list(paths)
    newest, newest_path, absent = None, None, []
    for p in given:
        f = Path(p)
        if not f.exists():
            absent.append(str(p))
            continue
        m = f.stat().st_mtime
        if newest is None or m > newest:
            newest, newest_path = m, str(f)
    return {"newest_input": (None if newest is None else
                             time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(newest))),
            "newest_path": newest_path, "absent": absent, "n_checked": len(given)}


# ----------------------------------------------------------------------------- section 1 · verdict


def normalise_findings(payload, defects: list) -> dict:
    """Flatten every gate's findings into one list, keeping the step that raised each.

    Accepts dicts or the Finding dataclasses the step modules return, because making the
    orchestrator convert them is one more place a message can be dropped. A severity this
    function does not recognise is kept under its own heading rather than being rounded down to
    ok: an unrecognised verdict is not a passing one.

    Both list-shaped fields - the gate list itself and each finding's `detail` - go through
    `_as_list`, so a single finding supplied where a list was expected is one finding and a
    single detail string is one bullet, rather than one entry per character.
    """
    raw = _get(payload, "gates")
    # `_unknown` as well as MISSING: `gates: null` is not `gates: []`. An empty list says nothing
    # was raised; a null says nobody wrote anything down, and reading the second as the first
    # would turn an unexamined run into a PASS.
    if raw is MISSING or _unknown(raw):
        _defect(defects, "section 1 · verdict",
                "no gate findings were supplied. This report therefore cannot say the run "
                "passed - only that nobody recorded whether it did.")
        return {"supplied": False, "findings": [], "refuse": [], "review": [],
                "ok": [], "unrecognised": []}
    out = []
    for i, f in enumerate(_as_list(raw)):
        step = _attr(f, "step")
        if not _stated(step):
            _defect(defects, "section 1 · verdict",
                    f"finding {i} carries no step. Every REFUSE and REVIEW must name the step "
                    f"that raised it or it cannot be acted on.")
        sev_raw = _attr(f, "severity")
        sev = _SEVERITY_ALIASES.get(str(sev_raw).strip().upper()) if _stated(sev_raw) else None
        if sev is None:
            sev = "unrecognised"
            _defect(defects, "section 1 · verdict",
                    f"finding {i} ({_text(_attr(f, 'check'))}) has severity {sev_raw!r}, which "
                    f"is not REFUSE, REVIEW or ok. It is listed separately rather than treated "
                    f"as a pass.")
        out.append({"step": _text(step, "step NOT STATED"),
                    "check": _text(_attr(f, "check")),
                    "severity": sev,
                    "message": _text(_attr(f, "message")),
                    "detail": [str(d) for d in _as_list(_attr(f, "detail"))]})
    return {"supplied": True, "findings": out,
            "refuse": [f for f in out if f["severity"] == "REFUSE"],
            "review": [f for f in out if f["severity"] == "REVIEW"],
            "ok": [f for f in out if f["severity"] == "ok"],
            "unrecognised": [f for f in out if f["severity"] == "unrecognised"]}


def verdict_of(gates: dict) -> str:
    """One severity for the run. NOT DETERMINED when no findings were supplied - never PASS."""
    if not gates.get("supplied"):
        return "NOT DETERMINED"
    if gates["refuse"]:
        return "REFUSE"
    if gates["review"] or gates["unrecognised"]:
        return "REVIEW"
    return "PASS"


def build_deliverable(payload, defects: list) -> dict:
    """The one line at the top. Composed from counts only when both counts are present.

    Every arithmetic step here is guarded, because this function runs inside `assemble()` and a
    bare ValueError out of `assemble()` destroys the whole report - the one artifact that records
    why the run stopped. A count that was supplied and cannot be read as a finite number is
    therefore rendered as a visible DEFECT carrying the value verbatim, and an input population
    of zero states that the percentage is undefined instead of dividing by it. Neither case is
    allowed to look like a deliverable nobody wrote down.
    """
    d = _get(payload, "deliverable")
    if d is MISSING:
        _defect(defects, "section 1 · verdict",
                "no deliverable block was supplied, so the document cannot state what this run "
                "produced or that it produced nothing.")
        return {"text": f"{NOT_STATED} — no deliverable block was supplied", "source": "absent"}
    text = _get(d, "text")
    if _stated(text):
        return {"text": str(text), "source": "supplied"}
    kept, total = _get(d, "n_kept"), _get(d, "n_in")
    unit = _text(_get(d, "unit"), "observations")
    stopped = _get(d, "stopped_after")
    unusable = None
    if _stated(kept) and _stated(total):
        n_kept, n_in = _as_number(kept), _as_number(total)
        if n_kept is not None and n_in is not None and n_in > 0:
            removed = 100.0 * (n_in - n_kept) / n_in
            return {"text": f"{_count_text(n_kept)} of {_count_text(n_in)} {unit}  "
                            f"({removed:.1f}% removed)",
                    "source": "composed from n_kept and n_in",
                    "n_kept": n_kept, "n_in": n_in, "pct_removed": round(removed, 4)}
        if n_kept is None or n_in is None:
            bad = ", ".join(f"{name}={value!r}" for name, value, num in
                            (("n_kept", kept, n_kept), ("n_in", total, n_in)) if num is None)
            _defect(defects, "section 1 · verdict",
                    f"the deliverable counts were supplied and cannot be read as finite numbers "
                    f"({bad}). No deliverable line is composed from them: a count this report "
                    f"cannot parse must not be printed as though it had been checked.")
            unusable = {"text": f"{NOT_STATED} — the deliverable counts could not be read as "
                                f"numbers ({bad})",
                        "source": "counts supplied but not numeric",
                        "n_kept_raw": str(kept), "n_in_raw": str(total)}
        else:
            _defect(defects, "section 1 · verdict",
                    f"the deliverable states an input population of {_count_text(n_in)}. No "
                    f"removed percentage is defined against it, so none is printed and the two "
                    f"counts are shown as they were supplied.")
            unusable = {"text": f"{_count_text(n_kept)} of {_count_text(n_in)} {unit}  "
                                f"(% removed UNDEFINED — the input population was "
                                f"{_count_text(n_in)})",
                        "source": "composed from n_kept and n_in; no percentage is defined",
                        "n_kept": n_kept, "n_in": n_in, "pct_removed": None}
    # A stop reason outranks a pair of counts that could not compose a line: it says what the run
    # DID, where the counts only say that something about them was wrong - and the defect above
    # has already recorded that on the page either way.
    if _stated(stopped):
        why = _text(_get(d, "stopped_because"), "no reason recorded")
        return {"text": f"NO DELIVERABLE — the run stopped after {stopped}: {why}",
                "source": "composed from stopped_after", "stopped_after": str(stopped)}
    if unusable is not None:
        return unusable
    _defect(defects, "section 1 · verdict",
            "the deliverable block supplied neither a text, nor n_kept with n_in, nor a "
            "stopped_after. The reader cannot tell a run that produced nothing from one whose "
            "result was not written down.")
    return {"text": f"{NOT_STATED} — the deliverable block is empty", "source": "empty"}


def build_verdict(payload, defects: list) -> dict:
    gates = normalise_findings(payload, defects)
    return {"overall": verdict_of(gates),
            "deliverable": build_deliverable(payload, defects),
            "gates_supplied": gates["supplied"],
            "refuse": gates["refuse"], "review": gates["review"],
            "unrecognised": gates["unrecognised"],
            "n_ok": len(gates["ok"]), "n_findings": len(gates["findings"]),
            "findings": gates["findings"]}


def verdict_lines(doc: dict) -> list:
    """The section-1 block, in the shape docs/REPORT_DESIGN.md specifies.

    Built as data so the JSON companion carries exactly the lines the HTML shows; a reader
    comparing the two never has to wonder which was assembled from what.
    """
    width = 14
    v = doc["verdict"]
    lines = [f"{'DELIVERABLE':<{width}}{v['deliverable']['text']}"]
    if v["deliverable"].get("stopped_after"):
        lines.append(f"{'STOPPED':<{width}}after {v['deliverable']['stopped_after']}")

    def block(label, entries, fmt):
        if not entries:
            lines.append(f"{label:<{width}}none")
            return
        for i, e in enumerate(entries):
            lines.append(f"{label if i == 0 else '':<{width}}{fmt(e)}")

    def gate_fmt(f):
        return f"{f['step']}  {f['check']} — {f['message']}"

    if not v["gates_supplied"]:
        lines.append(f"{'REFUSE':<{width}}NOT DETERMINED — no gate findings were supplied")
        lines.append(f"{'REVIEW':<{width}}NOT DETERMINED — no gate findings were supplied")
    else:
        block("REFUSE", v["refuse"], gate_fmt)
        block("REVIEW", v["review"], gate_fmt)
        if v["unrecognised"]:
            block("UNRECOGNISED", v["unrecognised"], gate_fmt)
    items = doc["open_items"]
    if not items["stated"]:
        lines.append(f"{'OPEN':<{width}}{NOT_STATED} — nobody recorded whether anything is open")
    else:
        block("OPEN", items["items"], lambda e: e["item"])
    n_def = sum(1 for d in doc["defects"] if d["severity"] == "defect")
    lines.append(f"{'DEFECTS':<{width}}"
                 + ("none" if not n_def else
                    f"{n_def} — this report is incomplete; see the marked blocks"))
    fresh = doc["provenance"]["freshness"]
    if fresh["stale"] is True:
        lines.append(f"{'STALE':<{width}}this artifact is OLDER than its newest input "
                     f"({fresh['generated']} < {fresh['newest_input']})")
    elif fresh["stale"] is None:
        lines.append(f"{'FRESHNESS':<{width}}{fresh['reason']}")
    return lines


# -------------------------------------------------------------------------- section 2 · parameters


def build_parameter_table(payload, defects: list) -> dict:
    """One row per parameter, with the class that says who was allowed to set it.

    An ADJUDICATED parameter with no verbatim operator text is not rendered as a row. It renders
    as a defect in the row's place, because the class is the claim - it asserts that a human
    decided this after seeing the evidence - and a row with nobody's words behind it makes that
    claim without any way to check it.

    The table goes through `_as_list`: a single parameter dict supplied where a list of them was
    expected is one row, not one row per key of that dict.
    """
    rows = _get(payload, "parameters")
    if rows is MISSING or _unknown(rows):        # `parameters: null` is not `parameters: []`
        _defect(defects, "section 2 · what was decided",
                "no parameter table was supplied. The class of each parameter - who was allowed "
                "to set it - is the point of this report, and this run recorded none of it.")
        return {"stated": False, "rows": []}
    out = []
    for i, p in enumerate(_as_list(rows)):
        name = _get(p, "name")
        cls = _get(p, "class")
        where = f"section 2 · parameter {i}" + (f" ({name})" if _stated(name) else "")
        if not _stated(name):
            _defect(defects, where, "a parameter row has no name.")
        cls_txt = str(cls).strip().upper() if _stated(cls) else None
        if cls_txt is None:
            _defect(defects, where, "the parameter has no class. Without one the report cannot "
                                    "say whether the data produced this value or a person chose "
                                    "it, and the two carry completely different weight.")
        elif cls_txt not in PARAM_CLASSES:
            _defect(defects, where, f"class {cls!r} is not one of "
                                    f"{', '.join(sorted(PARAM_CLASSES))}.")
            cls_txt = None
        value, basis = _get(p, "value"), _get(p, "basis")
        verbatim = _get(p, "verbatim")
        if cls_txt == "ADJUDICATED" and not _stated(verbatim):
            _defect(defects, where,
                    "ADJUDICATED with no verbatim operator text. The row is withheld: an "
                    "adjudicated value asserts that a person decided it after reading the "
                    "evidence, and without their own words that assertion cannot be checked "
                    "and must not be printed as though it could.")
            out.append({"name": _text(name), "defect": True,
                        "message": "ADJUDICATED without verbatim operator text — row withheld"})
            continue
        if not _stated(value):
            _defect(defects, where, "the parameter has no value.")
        if cls_txt in ("DERIVED", "ADJUDICATED") and not _stated(basis):
            _defect(defects, where, f"a {cls_txt} parameter with no basis. What it was derived "
                                    f"from, or read against, is the only thing that makes it "
                                    f"reviewable.")
        out.append({"name": _text(name), "value": _text(value), "defect": False,
                    # `is None`, not `cls_txt or ...`: the class has already been validated
                    # against PARAM_CLASSES above, so None is the one thing it can be that is
                    # not a class, and truthiness is not how this file reads a value.
                    "class": NOT_STATED if cls_txt is None else cls_txt,
                    "class_meaning": "" if cls_txt is None else PARAM_CLASSES.get(cls_txt, ""),
                    "basis": _text(basis), "verbatim": _text(verbatim, ""),
                    "decided_by": _text(_get(p, "decided_by"), ""),
                    "decided_on": _text(_get(p, "decided_on"), "")})
    return {"stated": True, "rows": out}


# ------------------------------------------------------------------------------- section 3 · steps


def _step_findings(step_key: str, verdict_block: dict) -> list:
    return [f for f in verdict_block["findings"] if f["step"] == step_key]


def build_step_sections(payload, verdict_block: dict, defects: list) -> list:
    """One section per step, in pipeline order, whether or not the run reached it.

    A step the payload says nothing about is rendered as NO RECORD SUPPLIED rather than omitted.
    Omitting it would let a run that stopped at step 2 produce a document that looks like a
    complete five-step one, and the reader's only clue would be a section count nobody counts.
    """
    supplied = {}
    order = []
    raw = _get(payload, "steps")
    if raw is MISSING or _unknown(raw):          # `steps: null` is not `steps: []`
        _defect(defects, "section 3 · the steps",
                "no step records were supplied, so this report can say nothing about what any "
                "step did.")
    else:
        for s in _as_list(raw):
            key = _get(s, "key")
            if not _stated(key):
                _defect(defects, "section 3 · the steps", "a step record has no key; it is "
                                                          "listed at the end under its index.")
                key = f"(unkeyed step {len(order)})"
            supplied[str(key)] = s
            order.append(str(key))

    stopped_after = None
    d = _get(payload, "deliverable")
    if d is not MISSING and _stated(_get(d, "stopped_after")):
        stopped_after = str(_get(d, "stopped_after"))

    known = [k for k, *_ in STEPS]
    extra = [k for k in order if k not in known]
    sections = []
    reached_stop = False
    for key, title, purpose, removes, figures in STEPS:
        entry = supplied.get(key, MISSING)
        sections.append(_one_step(key, title, purpose, removes, figures, entry, payload,
                                  verdict_block, defects,
                                  after_stop=reached_stop, stopped_after=stopped_after))
        if stopped_after is not None and key == stopped_after:
            reached_stop = True
    for key in extra:
        sections.append(_one_step(key, key, None, None, (), supplied[key], payload,
                                  verdict_block, defects,
                                  after_stop=False, stopped_after=stopped_after))
        _defect(defects, f"section 3 · {key}",
                "this step key is not one of the pipeline's eight. It is rendered at the end "
                "rather than dropped.", severity="notice")
    return sections


def _one_step(key, title, purpose, removes, figures, entry, payload, verdict_block, defects,
              *, after_stop: bool, stopped_after) -> dict:
    where = f"section 3 · {key}"
    if entry is MISSING:
        # Two different absences, and collapsing them would state more than is known: a step
        # AFTER a declared stop demonstrably did not run, while a step with no record at all
        # might have run and gone unrecorded. Only the first entails anything about the data.
        status = ("did not run — the run stopped after " + str(stopped_after)) if after_stop \
            else "NO RECORD SUPPLIED — this report cannot say whether it ran"
        entry_dict = {}
        ran = False
    else:
        entry_dict = entry if isinstance(entry, dict) else {}
        status = _text(_get(entry_dict, "status"), "status NOT STATED")
        ran = str(status).strip().lower() not in (
            "not run", "did not run", "skipped", "no record supplied", "status not stated")

    what = _get(entry_dict, "what_it_does")
    if _stated(what):
        purpose_text, purpose_source = str(what), "supplied by the run"
    elif purpose is not None:
        purpose_text, purpose_source = purpose, "README steps table"
    else:
        purpose_text, purpose_source = NOT_STATED, "absent"
        _defect(defects, where, "no description of what this step does was supplied and the "
                                "step is not one the pipeline documents.", severity="notice")

    found_raw = _get(entry_dict, "found")
    found_stated = found_raw is not MISSING and not _unknown(found_raw)
    found = []
    if found_stated:
        for j, item in enumerate(_as_list(found_raw)):
            if isinstance(item, dict):
                label, value, source = _get(item, "label"), _get(item, "value"), \
                    _get(item, "source")
                if _stated(value) and not _stated(source):
                    _defect(defects, where,
                            f"finding {j} ({_text(label)}) states a value with no source file. "
                            f"Every number in this report must be traceable to a file a reader "
                            f"can open.")
                found.append({"label": _text(label, ""), "value": _text(value),
                              "source": _text(source, f"source {NOT_STATED}")})
            else:
                found.append({"label": "", "value": str(item),
                              "source": f"source {NOT_STATED}"})
                _defect(defects, where, f"finding {j} was supplied as free text with no source "
                                        f"file.", severity="notice")

    cannot = _get(entry_dict, "cannot_establish")
    no_record = entry is MISSING and not after_stop
    if _stated(cannot):
        cannot_text, cannot_source = str(cannot), "supplied by the run"
    elif no_record:
        cannot_text = ("NO RECORD SUPPLIED — this report cannot say what this step established, "
                       "cannot say what it could not establish, and cannot say that it ran.")
        cannot_source = "absent"
        _defect(defects, where, "no record of this step at all. Whether it ran is unknown, and "
                                "unknown is not the same as did-not-run: nothing here may be "
                                "read as evidence about the data either way.")
    elif not ran:
        cannot_text = ("THIS STEP DID NOT RUN — nothing was established by it and nothing about "
                       "the data is claimed here.")
        cannot_source = "entailed by the step's status"
        _defect(defects, where, "no 'cannot establish' text; the step did not run, so the "
                                "statement entailed by that status is shown instead.",
                severity="notice")
    else:
        cannot_text, cannot_source = CANNOT_ESTABLISH_MISSING, "absent"
        _defect(defects, where, "the step ran and supplied no 'cannot establish' text. The "
                                "omission is what makes an honest report read like an "
                                "over-confident one.")

    # `_as_list`, not a bare comprehension: a step that declared its figures as the string "F3"
    # would otherwise be read as three figures called F, 3 and nothing.
    declared_figures = _get(entry_dict, "figures")
    fig_ids = list(figures) if (declared_figures is MISSING or _unknown(declared_figures)) \
        else [str(f) for f in _as_list(declared_figures)]
    fig_blocks = []
    payload_figs = _get(payload, "figures")
    for fid in fig_ids:
        spec = MISSING if payload_figs is MISSING else _get(payload_figs, fid)
        block = {"id": fid, "question": FIGURE_QUESTIONS.get(fid, ""),
                 "caption": "", "source": "", "status": ""}
        if spec is MISSING:
            # A reason, when the assembler gave one. "NOT PRODUCED" alone reads as an oversight;
            # "NOT PRODUCED — step 5 fits a KDE, takes the minimum and discards the curve" tells
            # the reader it is a gap in what the pipeline records, and tells the next person
            # exactly what to change to close it.
            fig_notes = _get(payload, "figure_notes")
            reason = "" if fig_notes is MISSING else _text(_get(fig_notes, fid), "")
            block["status"] = ("NOT PRODUCED — " + reason if reason
                               else "NOT PRODUCED — no data or image was supplied for this "
                                    "figure")
            _defect(defects, where, f"figure {fid} is expected at this step and nothing was "
                                    f"supplied for it.", severity="notice" if not ran
                    else "defect")
        else:
            block["caption"] = _text(_get(spec, "caption"), "")
            block["source"] = _text(_get(spec, "source"), "")
            if not block["caption"]:
                _defect(defects, where, f"figure {fid} has no caption.", severity="notice")
            if not block["source"]:
                _defect(defects, where, f"figure {fid} names no source file for the data it "
                                        f"draws.", severity="notice")
        fig_blocks.append(block)

    sources = _get(entry_dict, "sources")
    return {"key": key, "title": title, "status": status, "ran": ran,
            "purpose": purpose_text, "purpose_source": purpose_source,
            "removes": removes if removes is not None else NOT_STATED,
            "found": found, "found_stated": found_stated,
            "cannot_establish": cannot_text, "cannot_establish_source": cannot_source,
            "figures": fig_blocks,
            "findings": _step_findings(key, verdict_block),
            "sources": [str(s) for s in _as_list(sources)]}


# -------------------------------------------------------------------------- section 4 · provenance


def build_tool_rows(tools, expected) -> list:
    """Every tool, with the version this run OBSERVED or the words `not invoked`.

    Three states, never two. A tool absent from the observed set was not invoked; a tool recorded
    as `not invoked` says so itself; anything else is a string the run got by asking the
    executable. A blank is none of those and is rendered as not invoked rather than as an empty
    cell, because an empty cell reads as "no version needed".

    `expected` goes through `_as_list`: a single tool name supplied as a string used to become
    one expected tool per character, and the table then listed nine tools called c, e, l, l...
    """
    expected = tuple(_as_list(expected))
    tools = {} if not isinstance(tools, dict) else tools
    names = sorted(set(expected) | set(tools))
    rows = []
    for name in names:
        raw = _get(tools if isinstance(tools, dict) else {}, name)
        if raw is MISSING:
            rows.append({"name": name, "version": NOT_INVOKED,
                         "basis": "absent from the observed set"})
        elif not _stated(raw) or str(raw).strip().lower() == NOT_INVOKED:
            rows.append({"name": name, "version": NOT_INVOKED,
                         "basis": "recorded as not invoked"})
        else:
            rows.append({"name": name, "version": str(raw).strip(), "basis": "observed"})
    return rows


def build_provenance(payload, defects: list, *, now=None) -> dict:
    prov = _get(payload, "provenance")
    if prov is MISSING:
        _defect(defects, "section 4 · provenance",
                "no provenance was supplied: no commit, no tool versions, no input check. "
                "Nothing in this report can be tied to the code or the tools that produced it.")
        prov = {}
    pipe = _get(prov, "pipeline")
    pipe = pipe if isinstance(pipe, dict) else {}
    # `_tristate`, so a numpy bool is read as the determination it is and an unrecognised value
    # is read as no determination at all - never rounded down to "clean".
    dirty_state = _tristate(_get(pipe, "dirty"))
    if dirty_state is None:
        dirty = "UNKNOWN — cleanliness of the tree was not determined"
        _defect(defects, "section 4 · provenance",
                "whether the working tree was dirty is unknown. It is printed as unknown and "
                "not as clean.")
    else:
        dirty = "yes" if dirty_state else "no"
    commit = _get(pipe, "commit")
    if not _stated(commit):
        _defect(defects, "section 4 · provenance", "no commit was recorded.")

    generated = _get(prov, "generated")
    if _stated(generated):
        gen, gen_source = str(generated), "supplied by the run"
    elif _stated(now):
        gen, gen_source = str(now), "supplied to the report writer"
    else:
        gen, gen_source = time.strftime("%Y-%m-%dT%H:%M:%S%z"), "observed when the file was "\
                                                                "written"

    newest = _get(prov, "newest_input")
    if not _stated(newest):
        inputs = _get(prov, "inputs")
        best = None
        if isinstance(inputs, dict):
            for meta in inputs.values():
                m = _get(meta, "mtime") if isinstance(meta, dict) else MISSING
                if _stated(m):
                    ts, _ = _parse_iso(m)
                    if ts is not None and (best is None or ts > best[0]):
                        best = (ts, str(m))
        newest = best[1] if best is not None else MISSING
    fresh = freshness(gen, newest)
    if fresh["stale"] is True:
        _defect(defects, "section 4 · provenance",
                "this artifact is older than its newest input. A stale report opens, renders "
                "and reads exactly like a correct one, which is why it is stated here and in "
                "report.json rather than warned about in a log.")
    elif fresh["stale"] is None:
        _defect(defects, "section 4 · provenance",
                f"freshness was not checked: {fresh['reason']}", severity="notice")

    decisions = _get(prov, "decisions")
    decisions = decisions if isinstance(decisions, dict) else {}
    if not _stated(_get(decisions, "hash")):
        _defect(defects, "section 4 · provenance",
                "no hash for a decisions file. Two runs with the same data and the same "
                "decisions file are the same run; a command line does not have that property.",
                severity="notice")

    reference = _get(prov, "reference")
    reference = reference if isinstance(reference, dict) else {}
    tool_rows = build_tool_rows(_get(prov, "tools"), _get(prov, "tools_expected"))
    if not tool_rows:
        _defect(defects, "section 4 · provenance", "no tools were recorded, invoked or not.")

    checks = _get(prov, "input_check")
    input_rows = []
    if checks is MISSING or _unknown(checks):    # `input_check: null` is not `input_check: []`
        _defect(defects, "section 4 · provenance",
                "the raw-input verification (P1 raw VALUES, P2 raw DROPLETS) was not recorded. "
                "A pre-filtered matrix cannot be un-filtered and nothing downstream detects it.")
    else:
        for c in _as_list(checks):
            p1, p2 = _tristate(_get(c, "p1_raw_values")), _tristate(_get(c, "p2_raw_droplets"))
            input_rows.append({
                "name": _text(_get(c, "name")),
                # Three states, not two: a check that was not recorded is NOT CHECKED and it is
                # not a PASS. Read through `_tristate` rather than `p1 is True`, which is False
                # for `numpy.bool_(True)` - so a verification that genuinely passed, arriving
                # from any numpy-backed table, was printed as NOT CHECKED - and rather than a
                # bare truth test, which would read the string "false" as a pass.
                "p1": "NOT CHECKED" if p1 is None else ("PASS" if p1 else "FAIL"),
                "p2": "NOT CHECKED" if p2 is None else ("PASS" if p2 else "FAIL"),
                "reasons": [str(r) for r in _as_list(_get(c, "reasons"))]})

    run = _get(payload, "run")
    run = run if isinstance(run, dict) else {}
    invocation = _get(run, "invocation")
    if not _stated(invocation):
        _defect(defects, "section 4 · provenance", "the invocation was not recorded.",
                severity="notice")

    block = {
        "pipeline": {"version": _text(_get(pipe, "version")), "commit": _text(commit),
                     "dirty": dirty, "describe": _text(_get(pipe, "describe"), ""),
                     "branch": _text(_get(pipe, "branch"), "")},
        "invocation": _text(invocation),
        "decisions": {"path": _text(_get(decisions, "path")),
                      "hash": _text(_get(decisions, "hash"))},
        "reference": {k: _text(_get(reference, k)) for k in
                      ("registry", "genome", "annotation", "introns")},
        "tools": tool_rows,
        "input_check": input_rows,
        "environment": _get(prov, "environment") if isinstance(_get(prov, "environment"), dict)
        else {},
        "generated": gen, "generated_source": gen_source,
        "freshness": fresh,
    }
    block["lines"] = _provenance_lines(block, run)
    return block



def _runtime_line(run: dict) -> str:
    """Wall-clock, summed task time and the speed-up, or NOT STATED if the run did not record it."""
    el = _get(run, "elapsed_s")
    ts = _get(run, "task_seconds_total")
    jobs = _get(run, "jobs")
    if not _stated(el):
        return "NOT STATED - the run did not record its own duration"
    def _hms(s):
        s = float(s)
        h, rem = divmod(int(s), 3600)
        m, sec = divmod(rem, 60)
        return (f"{h}h {m:02d}m {sec:02d}s" if h else
                (f"{m}m {sec:02d}s" if m else f"{sec}s"))
    out = f"wall {_hms(el)}"
    if _stated(ts) and float(ts) > 0 and float(el) > 0:
        out += f"   task-time {_hms(ts)}   speed-up {float(ts) / float(el):.1f}x"
    if _stated(jobs):
        out += f"   ({jobs} concurrent)"
    return out


def _provenance_lines(block: dict, run: dict) -> list:
    w = 16
    p = block["pipeline"]
    lines = [f"{'pipeline':<{w}}scQC {p['version']}  commit {p['commit']}  (dirty: {p['dirty']})",
             f"{'invocation':<{w}}{block['invocation']}",
             f"{'project':<{w}}{_text(_get(run, 'project'))}"
             f"   mode: {_text(_get(run, 'mode'))}",
             # Wall-clock AND summed task time. The ratio is what concurrency bought; one number
             # alone hides either how long a person waited or whether the machine was used.
             f"{'runtime':<{w}}" + _runtime_line(run),
             f"{'decisions':<{w}}{block['decisions']['path']}  {block['decisions']['hash']}",
             f"{'reference':<{w}}" + "  ".join(f"{k}: {v}" for k, v in
                                               block["reference"].items())]
    if not block["tools"]:
        lines.append(f"{'tools':<{w}}{NOT_STATED} — no tool was recorded, invoked or not")
    for i, t in enumerate(block["tools"]):
        lines.append(f"{'tools' if i == 0 else '':<{w}}{t['name']:<16}{t['version']}"
                     f"   [{t['basis']}]")
    if not block["input_check"]:
        lines.append(f"{'input check':<{w}}NOT CHECKED — P1 raw VALUES and P2 raw DROPLETS were "
                     f"not recorded")
    for i, c in enumerate(block["input_check"]):
        lines.append(f"{'input check' if i == 0 else '':<{w}}{c['name']:<28}"
                     f"P1 raw VALUES: {c['p1']}   P2 raw DROPLETS: {c['p2']}")
    f = block["freshness"]
    lines.append(f"{'generated':<{w}}{block['generated']}   [{block['generated_source']}]")
    lines.append(f"{'newest input':<{w}}{f['newest_input']}")
    lines.append(f"{'freshness':<{w}}"
                 + ("STALE — " + f["reason"] if f["stale"] is True else f["reason"]))
    return lines


# -------------------------------------------------------------------------- section 5 · open items


def build_open_items(payload, defects: list) -> dict:
    items = _get(payload, "open_items")
    # `open_items: null` is not `open_items: []`. The first is nobody looking, the second is a
    # statement that nothing is unresolved - and this section exists to keep them apart.
    if items is MISSING or _unknown(items):
        _defect(defects, "section 5 · open items",
                "no open-items list was supplied. A missing section and an empty one are not "
                "the same claim: this report cannot say that nothing is unresolved.")
        return {"stated": False, "items": []}
    out = []
    for i, it in enumerate(_as_list(items)):
        if isinstance(it, dict):
            text = _get(it, "item")
            if not _stated(text):
                _defect(defects, "section 5 · open items", f"open item {i} has no text.")
            out.append({"item": _text(text), "closes_when": _text(_get(it, "closes_when")),
                        "blocked_on": _text(_get(it, "blocked_on"))})
        else:
            out.append({"item": str(it), "closes_when": NOT_STATED, "blocked_on": NOT_STATED})
            _defect(defects, "section 5 · open items",
                    f"open item {i} says what is open but not what would close it or who it is "
                    f"blocked on.", severity="notice")
    return {"stated": True, "items": out}


def build_per_sample(payload, defects: list) -> dict:
    """Every threshold this run derived, one row per library.

    A pipeline that derives some thresholds per library and others for the cohort produces a
    result nobody can check without knowing WHICH IS WHICH, and that distinction survives nowhere
    in a report organised by step: a per-library ceiling and a cohort floor both appear as "a
    number step 5 produced". So each column carries its scope, and a cohort constant is repeated
    down its column rather than shown once - a reader comparing two libraries must be able to see
    which numbers differ because the libraries differ.

    A cell nothing recorded stays NOT STATED. It must not read as a zero, and it must not read as
    a threshold that was applied.
    """
    block = _get(payload, "per_sample")
    if block is MISSING or _unknown(block):
        _defect(defects, "section 3 · per-library thresholds",
                "no per-library threshold table was supplied. Which thresholds vary by library "
                "and which are one cohort constant is not derivable from the rest of this "
                "document, and it decides how every number in it may be compared.")
        return {"stated": False, "columns": [], "rows": [], "source": NOT_STATED}

    cols = []
    for i, c in enumerate(_as_list(_get(block, "columns"))):
        key = _get(c, "key")
        if not _stated(key):
            _defect(defects, "section 3 · per-library thresholds",
                    f"column {i} has no key, so its cells cannot be read.")
            continue
        scope = _get(c, "scope")
        if not _stated(scope):
            _defect(defects, "section 3 · per-library thresholds",
                    f"column {key!r} does not say whether it is per library or a cohort "
                    f"constant. Without that the column cannot be interpreted.", severity="notice")
        cols.append({"key": str(key), "label": _text(_get(c, "label"), str(key)),
                     "scope": _text(scope), "step": _text(_get(c, "step"), "")})

    rows = []
    for r in _as_list(_get(block, "rows")):
        sample = _get(r, "sample")
        if not _stated(sample):
            _defect(defects, "section 3 · per-library thresholds",
                    "a row carries no sample name, so its numbers belong to no library.")
            continue
        cells = {}
        for c in cols:
            v = _get(r, c["key"])
            cells[c["key"]] = NOT_STATED if (v is MISSING or _unknown(v)) else v
        rows.append({"sample": str(sample), "cells": cells})

    if not rows:
        _defect(defects, "section 3 · per-library thresholds",
                "the per-library table has no rows. An empty table and an absent one read the "
                "same on the page and are not the same claim.")
    source = _get(block, "source")
    if not _stated(source):
        _defect(defects, "section 3 · per-library thresholds",
                "the per-library table names no file. Every number in this report must be "
                "traceable to something a reader can open.")
    # Which columns actually VARY is the question the table exists to answer, so it is answered
    # here rather than left to the reader to scan for.
    varying = []
    for c in cols:
        seen = {str(r["cells"][c["key"]]) for r in rows}
        if len(seen) > 1:
            varying.append(c["key"])
    return {"stated": True, "columns": cols, "rows": rows, "source": _text(source),
            "varying": varying,
            "n_per_library": sum(1 for c in cols if c["scope"] == "per library"),
            "n_cohort": sum(1 for c in cols if c["scope"] == "cohort constant")}


# ---------------------------------------------------------------------------------- assembly


def assemble(payload, *, now=None) -> dict:
    """Everything the document says, as data. The HTML is rendered from this and nothing else.

    Separated from rendering so the JSON companion cannot drift from the page: there is one
    computation of every number and two presentations of it.
    """
    if not isinstance(payload, dict):
        raise TaskFailure(f"the report payload must be a dict, got {type(payload).__name__}. "
                          f"See the schema in this module's docstring.")
    defects: list = []
    verdict = build_verdict(payload, defects)
    parameters = build_parameter_table(payload, defects)
    steps = build_step_sections(payload, verdict, defects)
    provenance = build_provenance(payload, defects, now=now)
    open_items = build_open_items(payload, defects)
    run = _get(payload, "run")
    doc = {
        "report_version": REPORT_VERSION,
        "run": {"project": _text(_get(run, "project") if run is not MISSING else MISSING),
                "mode": _text(_get(run, "mode") if run is not MISSING else MISSING),
                "started": _text(_get(run, "started") if run is not MISSING else MISSING, "")},
        "verdict": verdict,
        "parameters": parameters,
        "steps": steps,
        "per_sample": build_per_sample(payload, defects),
        "provenance": provenance,
        "open_items": open_items,
        "figures": {},
        "defects": defects,
    }
    doc["verdict"]["defect_count"] = sum(1 for d in defects if d["severity"] == "defect")
    doc["verdict"]["notice_count"] = sum(1 for d in defects if d["severity"] == "notice")
    doc["verdict"]["lines"] = verdict_lines(doc)
    return doc


# ----------------------------------------------------------------------------------- figures


def _load_figures():
    """The figure module, whether this file was imported as a package or loaded by path."""
    try:
        from . import figures  # noqa: PLC0415
        return figures
    except ImportError:
        import importlib.util
        path = Path(__file__).resolve().parent / "figures.py"
        spec = importlib.util.spec_from_file_location("scqc_report_figures", path)
        if spec is None or spec.loader is None:                            # pragma: no cover
            raise TaskFailure(f"could not load the figure module from {path}") from None
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        return mod


def figure_data_uri(fig) -> str:
    """A PNG data URI for a matplotlib Figure. No file, no external request, no network.

    The Software tag is dropped where matplotlib allows it, so two runs over the same numbers
    produce the same bytes and a diff of two reports shows what changed rather than which
    matplotlib drew them.
    """
    buf = io.BytesIO()
    try:
        fig.savefig(buf, format="png", bbox_inches="tight", metadata={"Software": None})
    except (TypeError, ValueError, KeyError):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
    data = buf.getvalue()
    if not data:
        raise TaskFailure("a figure rendered to zero bytes; the report would carry an empty "
                          "image that looks like a blank result")
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")


def render_figures(payload, defects: list) -> dict:
    """Render every figure the payload describes, and record every one that could not be.

    A figure that fails to draw does NOT stop the report. The document is the only record of the
    run, and a run whose figure code raised is exactly the run somebody needs to read about; the
    failure is printed in the figure's place, with the exception text, rather than leaving a gap
    that reads like a figure nobody thought was needed.
    """
    figs = _get(payload, "figures")
    out = {"uris": {}, "index": {}, "matplotlib": NOT_INVOKED}
    if figs is MISSING or _unknown(figs):
        _defect(defects, "figures", "no figures were supplied.", severity="notice")
        return out
    if not isinstance(figs, dict):
        # Recorded, not raised. Raising here would take the whole document down over a malformed
        # figure block, and the document is the only record of the run that produced it.
        _defect(defects, "figures",
                f"payload['figures'] is a {type(figs).__name__}, not a dict keyed by figure id "
                f"(F1..F9). No figure could be resolved from it and every figure block in this "
                f"report is therefore empty.")
        return out
    module = None
    for fid in sorted(figs):
        spec = figs[fid] if isinstance(figs[fid], dict) else {}
        pre = _get(spec, "png_base64")
        if _stated(pre):
            out["uris"][fid] = "data:image/png;base64," + str(pre).strip()
            out["index"][fid] = {"rendered": True, "source": "supplied pre-rendered"}
            continue
        data = _get(spec, "data")
        if data is MISSING:
            out["index"][fid] = {"rendered": False,
                                 "source": "neither data nor png_base64 was supplied"}
            _defect(defects, f"figure {fid}",
                    "neither drawing data nor a rendered image was supplied.")
            continue
        if not isinstance(data, dict):
            out["index"][fid] = {"rendered": False,
                                 "source": f"data is {type(data).__name__}, not a dict of "
                                           f"keyword arguments"}
            _defect(defects, f"figure {fid}", "the 'data' entry must be a dict of keyword "
                                              "arguments for the figure function.")
            continue
        try:
            if module is None:
                module = _load_figures()
                import matplotlib
                out["matplotlib"] = str(matplotlib.__version__)
            fn = module.FIGURE_FUNCTIONS.get(fid)
            if fn is None:
                raise TaskFailure(f"{fid} is not one of {', '.join(sorted(FIGURE_QUESTIONS))}")
            fig = fn(**data)
            out["uris"][fid] = figure_data_uri(fig)
            out["index"][fid] = {"rendered": True, "source": f"drawn by {fn.__name__}"}
        except Exception as exc:                                          # noqa: BLE001
            out["index"][fid] = {"rendered": False,
                                 "source": f"{type(exc).__name__}: {exc}"}
            _defect(defects, f"figure {fid}",
                    f"could not be drawn — {type(exc).__name__}: {exc}")
    return out


# ------------------------------------------------------------------------------------- render


def _e(value) -> str:
    return html.escape(str(value), quote=True)


_CSS = """
:root{--ink:#1b1b1b;--paper:#ffffff;--rule:#d8d8d8;--muted:#5f5f5f;--ok:#0072B2;
--review:#E69F00;--refuse:#D55E00;--unknown:#CC79A7;--wash:#f6f6f4;}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);
font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
.wrap{max-width:1080px;margin:0 auto;padding:28px 22px 80px;}
h1{font-size:24px;margin:0 0 2px;letter-spacing:-.01em}
h2{font-size:18px;margin:38px 0 10px;padding-bottom:6px;border-bottom:2px solid var(--ink)}
h3{font-size:15px;margin:26px 0 6px}
p{margin:8px 0}
.sub{color:var(--muted);font-size:13px;margin:0 0 18px}
pre{background:var(--wash);border:1px solid var(--rule);border-radius:4px;padding:12px 14px;
overflow-x:auto;font:12.5px/1.5 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
white-space:pre;margin:10px 0}
table{border-collapse:collapse;width:100%;margin:10px 0;font-size:13px}
th,td{border:1px solid var(--rule);padding:6px 8px;text-align:left;vertical-align:top}
th{background:var(--wash);font-weight:600}
/* The per-library table is wider than the page on purpose - twenty-odd thresholds is what
   there are. It scrolls inside its own box rather than shrinking the type or dropping columns. */
.scroll{overflow-x:auto;border:1px solid var(--rule);border-radius:4px;margin:10px 0}
table.wide{margin:0;font-size:12px;white-space:nowrap;width:auto;min-width:100%}
table.wide th,table.wide td{border-left:none;border-right:1px solid var(--rule)}
table.wide tr td:first-child,table.wide tr th:first-child{position:sticky;left:0;
background:var(--paper);border-right:2px solid var(--ink)}
table.wide tr th:first-child{background:var(--wash)}
td.ns{color:var(--muted);font-style:italic}
.badge{display:inline-block;padding:3px 10px;border-radius:3px;color:#fff;font-weight:700;
font-size:12px;letter-spacing:.04em}
.b-REFUSE{background:var(--refuse)}.b-REVIEW{background:var(--review)}
.b-PASS{background:var(--ok)}.b-NOTDET{background:var(--unknown)}
.chip{display:inline-block;padding:1px 7px;border-radius:3px;font-size:11px;font-weight:700;
color:#fff}
.c-REFUSE{background:var(--refuse)}.c-REVIEW{background:var(--review)}.c-ok{background:var(--ok)}
.c-unrecognised{background:var(--unknown)}
.defect{border-left:4px solid var(--refuse);background:#fdf1ea;padding:9px 12px;margin:10px 0;
font-size:13px}
.notice{border-left:4px solid var(--unknown);background:#fbf0f6;padding:9px 12px;margin:10px 0;
font-size:13px}
.cannot{border-left:4px solid var(--muted);background:var(--wash);padding:9px 12px;margin:10px 0;
font-size:13px}
figure{margin:14px 0;padding:0}
figure img{max-width:100%;height:auto;border:1px solid var(--rule);border-radius:3px;
display:block}
figcaption{font-size:12px;color:var(--muted);margin-top:6px}
.q{font-style:italic}
.step{border:1px solid var(--rule);border-radius:5px;padding:2px 16px 14px;margin:16px 0;
background:var(--paper)}
.meta{font-size:12px;color:var(--muted)}
.src{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:11.5px;
color:var(--muted)}
.det{margin:4px 0 0 14px;padding:0;font-size:12px;color:var(--muted)}
.klass{font-weight:700;letter-spacing:.03em}
.k-ADJUDICATED{color:var(--refuse)}.k-DERIVED{color:var(--ok)}.k-DECLARED{color:var(--review)}
.k-FIXED{color:var(--muted)}
footer{margin-top:44px;padding-top:12px;border-top:1px solid var(--rule);font-size:12px;
color:var(--muted)}
"""


def _findings_table(findings) -> str:
    if not findings:
        return "<p class='meta'>no findings recorded for this step.</p>"
    rows = []
    for f in findings:
        detail = ""
        if f["detail"]:
            detail = "<ul class='det'>" + "".join(f"<li>{_e(d)}</li>" for d in f["detail"]) \
                     + "</ul>"
        rows.append(f"<tr><td><span class='chip c-{_e(f['severity'])}'>{_e(f['severity'])}"
                    f"</span></td><td>{_e(f['check'])}</td>"
                    f"<td>{_e(f['message'])}{detail}</td></tr>")
    return ("<table><tr><th>verdict</th><th>check</th><th>finding</th></tr>"
            + "".join(rows) + "</table>")


#: Elements that exist to pull something in from outside the document. Their presence is refused
#: whatever they point at: `<script src=...>` and `<link rel=stylesheet>` are the obvious ones,
#: `<base>` because it silently re-points every relative reference in the page.
_EXTERNAL_TAGS = frozenset(("script", "iframe", "frame", "frameset", "object", "embed", "applet",
                            "base", "link"))

#: Attributes whose value is a URL. Anything here must be a `data:` URI or a same-document
#: fragment; anything else is a request that leaves the file.
_URL_ATTRIBUTES = frozenset(("src", "srcset", "href", "data", "poster", "action", "formaction",
                             "background", "cite", "manifest", "codebase", "longdesc",
                             "xlink:href", "profile"))


def _is_local_reference(value) -> bool:
    """True when a URL-valued attribute keeps the document inside itself."""
    # `_unknown`, not `value or ""`: a NaN is truthy, so `or` would keep it and `str()` would
    # then hand this the four characters `nan`, which is reported as an external reference to a
    # host called nan. An absent value is an absent reference and fetches nothing.
    v = "" if _unknown(value) else str(value).strip()
    if not v:
        return True
    return v.startswith("#") or v[:5].lower() == "data:"


def _srcset_urls(value: str) -> list:
    """The candidate URLs in a `srcset`, split the way the HTML parser splits them.

    NOT `value.split(",")`. Every inline image in this report is a `data:` URI and every base64
    `data:` URI CONTAINS a comma - `data:image/png;base64,iVBOR...` - so splitting on commas tore
    one self-contained candidate into two, and the tail (`iVBOR...`) does not begin with `data:`
    and was reported as a reference leaving the document. `assert_self_contained()` then refused a
    page that was entirely self-contained, which is the failure docs/PRINCIPLES.md section 3 is
    about: a gate that fires on correct behaviour is a gate somebody switches off.

    The spec's rule is that a candidate's URL runs to the next WHITESPACE. A comma separates
    candidates only after the descriptor, or where it trails the URL itself - so:

        data:image/png;base64,AAA= 1x, data:image/png;base64,BBB= 2x

    is two candidates, both inline, and

        data:image/png;base64,AAA= 1x, https://cdn.example/x.png 2x

    is two candidates of which the second still has to be found. Parsing rather than splitting is
    what keeps both of those true at once: a comma inside a URL no longer hides the entry after
    it, and no longer invents an entry that was never written.
    """
    s = "" if _unknown(value) else str(value)
    out: list = []
    i, n = 0, len(s)
    while i < n:
        while i < n and (s[i].isspace() or s[i] == ","):        # between candidates
            i += 1
        if i >= n:
            break
        start = i
        while i < n and not s[i].isspace():                     # the URL: up to whitespace
            i += 1
        url = s[start:i]
        if url.endswith(","):
            # `url,` with no descriptor - the trailing commas are separators, not part of the URL.
            url = url.rstrip(",")
        else:
            # Skip this candidate's descriptor, ending at the comma that starts the next one.
            # Parenthesis depth is tracked because a descriptor may legally carry a `(1,2)`.
            depth = 0
            while i < n:
                c = s[i]
                if c == "(":
                    depth += 1
                elif c == ")" and depth:
                    depth -= 1
                elif c == "," and depth == 0:
                    i += 1
                    break
                i += 1
        if url:
            out.append(url)
    return out


def _css_external(text: str) -> list:
    """External references in a CSS fragment: `@import`, and `url(...)` that is not a data URI."""
    out = []
    low = str(text).lower()
    if "@import" in low:
        out.append("@import")
    i = 0
    while True:
        j = low.find("url(", i)
        if j < 0:
            break
        k = low.find(")", j)
        ref = text[j + 4:(k if k >= 0 else len(text))].strip().strip("\"'")
        if not _is_local_reference(ref):
            out.append(f"url({ref[:60]})")
        i = j + 4
    return out


def external_references(html_text: str) -> list:
    """Every place the finished document would reach outside itself, as readable strings.

    WHY THIS PARSES THE MARKUP INSTEAD OF SEARCHING FOR SUBSTRINGS

    The previous form searched the whole document for fragments like `srcset=`, `@import`,
    `url(http` and `<script`. `html.escape()` does not escape any of those - it escapes `&`, `<`,
    `>` and the quotes - so `<script` in a caption arrives as `&lt;script` and never matches,
    while a caption that merely MENTIONS `url(http://...)` or `srcset=` matches immediately and
    refuses a document that is perfectly self-contained. It was checking inert text and missing
    the structure.

    So the markup is parsed and the question is asked of the structure: is there an element whose
    job is to fetch something, or a URL-valued attribute pointing anywhere but at a data URI or a
    fragment of this same file? Escaped payload text produces no elements and no attributes at
    all, so it cannot trigger this, and a real external reference cannot hide from it.

    Returns a list of descriptions - empty when the document is self-contained.
    """
    from html.parser import HTMLParser

    found: list = []

    class _Scan(HTMLParser):
        def __init__(self):
            super().__init__(convert_charrefs=True)
            self._style_depth = 0

        def handle_starttag(self, tag, attrs):
            self._check(tag, attrs)
            if tag == "style":
                self._style_depth += 1

        def handle_startendtag(self, tag, attrs):
            self._check(tag, attrs)

        def handle_endtag(self, tag):
            if tag == "style" and self._style_depth:
                self._style_depth -= 1

        def handle_data(self, data):
            # Only inside <style>. Everywhere else this is escaped payload text, and text cannot
            # fetch anything - that is the whole reason the substring form gave false positives.
            if self._style_depth:
                found.extend(f"<style> contains {m}" for m in _css_external(data))

        def _check(self, tag, attrs):
            if tag in _EXTERNAL_TAGS:
                found.append(f"<{tag}> element")
            for raw_name, value in attrs:
                # `_unknown`, not `x or ""`. The parser hands back None for a valueless
                # attribute, and every `or` default in this file has been the same defect:
                # it also swallows the values that ARE information.
                name = ("" if _unknown(raw_name) else str(raw_name)).lower()
                text = "" if _unknown(value) else str(value)
                if name == "style" and text.strip():
                    found.extend(f"style attribute on <{tag}> contains {m}"
                                 for m in _css_external(text))
                    continue
                if name == "http-equiv" and text.strip().lower() == "refresh":
                    found.append(f"<{tag} http-equiv=refresh> redirects the document")
                    continue
                if name not in _URL_ATTRIBUTES or _unknown(value):
                    continue
                # srcset is a candidate list; `_srcset_urls` is why it is parsed rather than
                # split on commas - a base64 data URI contains commas of its own.
                refs = _srcset_urls(text) if name == "srcset" else [text.strip()]
                for ref in refs:
                    if not _is_local_reference(ref):
                        found.append(f"{name}={ref[:60]!r} on <{tag}>")

    scan = _Scan()
    try:
        scan.feed(html_text)
        scan.close()
    except Exception as exc:                                              # noqa: BLE001
        # Fail closed. A document that cannot be parsed has not been shown to be self-contained,
        # and "could not check" must never be recorded as "checked and clean".
        raise TaskFailure(f"the report could not be parsed as HTML, so it cannot be shown to be "
                          f"self-contained ({type(exc).__name__}: {exc})") from None
    return found


def assert_self_contained(html_text: str) -> None:
    """Raise unless the document can be opened with no network and no server.

    Every payload string is HTML-escaped before it reaches the page, so a caption containing
    `<script` arrives as text; this check is what proves that rather than assuming it. See
    `external_references()` for why it inspects the parsed markup rather than searching for
    substrings, which both false-positived on inert text and could not see the structure.
    """
    found = external_references(html_text)
    if found:
        raise TaskFailure(
            f"the report reaches outside itself: {sorted(set(found))}. The document must open "
            f"from a filesystem with no network: embed the resource as a data URI instead of "
            f"linking it.")


def render_html(doc: dict, figure_uris: dict) -> str:
    """The single-file document. Inline CSS, inline images, no external request of any kind."""
    v = doc["verdict"]
    badge = {"REFUSE": "b-REFUSE", "REVIEW": "b-REVIEW", "PASS": "b-PASS"}.get(
        v["overall"], "b-NOTDET")
    parts = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>scQC report — {_e(doc['run']['project'])}</title>",
        f"<style>{_CSS}</style></head><body><div class='wrap'>",
        f"<h1>scQC quality-control report <span class='badge {badge}'>{_e(v['overall'])}"
        f"</span></h1>",
        f"<p class='sub'>project {_e(doc['run']['project'])} · mode {_e(doc['run']['mode'])} · "
        f"generated {_e(doc['provenance']['generated'])} "
        f"({_e(doc['provenance']['generated_source'])})</p>",
    ]

    # ---------------------------------------------------------------- 1 verdict
    parts.append("<h2>1 · Verdict</h2>")
    parts.append("<pre>" + "\n".join(_e(line) for line in v["lines"]) + "</pre>")
    if not v["gates_supplied"]:
        parts.append("<div class='defect'><strong>Defect.</strong> No gate findings were "
                     "supplied. This document cannot state that the run passed; it can only "
                     "state that nobody recorded whether it did.</div>")
    if v["defect_count"]:
        parts.append(f"<div class='defect'><strong>{v['defect_count']} defect(s) in this "
                     f"report.</strong> Each is marked where it occurs and listed at the foot "
                     f"of the document. A report that omits a required block silently reads "
                     f"exactly like a complete one.</div>")
    if v["refuse"] or v["review"] or v["unrecognised"]:
        parts.append("<h3>every REFUSE and REVIEW, with the step that raised it</h3>")
        rows = []
        for f in v["refuse"] + v["review"] + v["unrecognised"]:
            detail = ""
            if f["detail"]:
                detail = "<ul class='det'>" + "".join(f"<li>{_e(d)}</li>" for d in f["detail"]) \
                         + "</ul>"
            rows.append(f"<tr><td><span class='chip c-{_e(f['severity'])}'>{_e(f['severity'])}"
                        f"</span></td><td class='src'>{_e(f['step'])}</td>"
                        f"<td>{_e(f['check'])}</td><td>{_e(f['message'])}{detail}</td></tr>")
        parts.append("<table><tr><th>verdict</th><th>step</th><th>check</th><th>finding</th>"
                     "</tr>" + "".join(rows) + "</table>")
    parts.append(f"<p class='meta'>{v['n_ok']} further finding(s) recorded as ok, shown with "
                 f"their steps in section 3.</p>")

    # ---------------------------------------------------------------- 2 parameters
    parts.append("<h2>2 · What was decided, and by whom</h2>")
    if not doc["parameters"]["stated"]:
        parts.append("<div class='defect'><strong>Defect.</strong> No parameter table was "
                     "supplied. The class of a parameter — who was allowed to set it — is the "
                     "point of this report.</div>")
    else:
        rows = []
        for r in doc["parameters"]["rows"]:
            if r.get("defect"):
                rows.append(f"<tr><td>{_e(r['name'])}</td><td colspan='3'>"
                            f"<span class='defect' style='display:block;margin:0'>"
                            f"<strong>Defect.</strong> {_e(r['message'])}</span></td></tr>")
                continue
            words = ""
            if r["verbatim"]:
                who = " — ".join(x for x in (r["decided_by"], r["decided_on"]) if x)
                words = (f"<br><span class='src'>&ldquo;{_e(r['verbatim'])}&rdquo;"
                         + (f" — {_e(who)}" if who else "") + "</span>")
            rows.append(
                f"<tr><td>{_e(r['name'])}</td><td>{_e(r['value'])}</td>"
                f"<td><span class='klass k-{_e(r['class'])}'>{_e(r['class'])}</span>"
                f"<br><span class='meta'>{_e(r['class_meaning'])}</span></td>"
                f"<td>{_e(r['basis'])}{words}</td></tr>")
        parts.append("<table><tr><th>parameter</th><th>value</th><th>class</th><th>basis</th>"
                     "</tr>" + "".join(rows) + "</table>")
        parts.append("<p class='meta'>An ADJUDICATED row carries the operator's own words. A "
                     "row without them is withheld and shown as a defect: the class asserts "
                     "that a person decided after reading the evidence, and that assertion is "
                     "not checkable without the words.</p>")

    # ---------------------------------------------------------------- 3 steps
    parts.append("<h2>3 · The steps</h2>")
    for s in doc["steps"]:
        parts.append("<div class='step'>")
        parts.append(f"<h3>{_e(s['title'])} <span class='meta'>— {_e(s['status'])}</span></h3>")
        parts.append(f"<p>{_e(s['purpose'])} <span class='meta'>[{_e(s['purpose_source'])}]"
                     f"</span></p>")
        parts.append(f"<p class='meta'>removes: {_e(s['removes'])}</p>")
        for fig in s["figures"]:
            uri = figure_uris.get(fig["id"])
            parts.append("<figure>")
            if _stated(uri):
                parts.append(f"<img alt='{_e(fig['id'])}' src='{_e(uri)}'>")
            else:
                parts.append(f"<div class='defect'><strong>Figure {_e(fig['id'])} is not "
                             f"here.</strong> {_e(_text(fig['status'], 'not rendered'))}</div>")
            cap = (f"<strong>{_e(fig['id'])}</strong> <span class='q'>"
                   f"{_e(fig['question'])}</span>")
            if _stated(fig["caption"]):
                cap += f" — {_e(fig['caption'])}"
            if fig["source"]:
                cap += f"<br><span class='src'>source: {_e(fig['source'])}</span>"
            parts.append(f"<figcaption>{cap}</figcaption></figure>")
        parts.append("<p><strong>What it found</strong></p>")
        if not s["found_stated"]:
            parts.append("<p class='meta'>no findings were recorded for this step.</p>")
        elif not s["found"]:
            parts.append("<p class='meta'>none recorded.</p>")
        else:
            rows = "".join(f"<tr><td>{_e(f['label'])}</td><td>{_e(f['value'])}</td>"
                           f"<td class='src'>{_e(f['source'])}</td></tr>" for f in s["found"])
            parts.append("<table><tr><th>quantity</th><th>value</th><th>source file</th></tr>"
                         + rows + "</table>")
        parts.append(_findings_table(s["findings"]))
        css = "defect" if s["cannot_establish_source"] == "absent" else "cannot"
        parts.append(f"<div class='{css}'><strong>What this step cannot establish.</strong> "
                     f"{_e(s['cannot_establish'])} <span class='meta'>"
                     f"[{_e(s['cannot_establish_source'])}]</span></div>")
        if s["sources"]:
            parts.append("<p class='src'>files: " + ", ".join(_e(x) for x in s["sources"])
                         + "</p>")
        parts.append("</div>")

    # ---------------------------------------------------------------- 4 provenance
    parts.append("<h2>4 · Provenance</h2>")
    parts.append("<pre>" + "\n".join(_e(line) for line in doc["provenance"]["lines"])
                 + "</pre>")
    fresh = doc["provenance"]["freshness"]
    if fresh["stale"] is True:
        parts.append("<div class='defect'><strong>STALE.</strong> This artifact is older than "
                     "its newest input. It is not out of date, it is wrong — and wrong in the "
                     "one way nobody checks, because it opens, renders and reads exactly like a "
                     "correct one. Rebuild it before it is shown to anyone.</div>")
    elif fresh["stale"] is None:
        parts.append(f"<div class='notice'><strong>Freshness NOT CHECKED.</strong> "
                     f"{_e(fresh['reason'])}</div>")
    for c in doc["provenance"]["input_check"]:
        if c["reasons"]:
            parts.append(f"<div class='notice'><strong>{_e(c['name'])}</strong><ul class='det'>"
                         + "".join(f"<li>{_e(r)}</li>" for r in c["reasons"]) + "</ul></div>")

    # ------------------------------------------- every threshold, per library
    ps = doc["per_sample"]
    parts.append("<h2>Every threshold this run derived, per library</h2>")
    if not ps["stated"]:
        parts.append("<div class='defect'><strong>Defect.</strong> No per-library threshold "
                     "table was supplied. Which thresholds vary by library and which are one "
                     "cohort constant is not derivable from the rest of this document, and it "
                     "decides how every number in it may be compared.</div>")
    else:
        parts.append(
            f"<p>{len(ps['rows'])} librar{'y' if len(ps['rows']) == 1 else 'ies'}, "
            f"<strong>{ps['n_per_library']}</strong> quantities derived per library and "
            f"<strong>{ps['n_cohort']}</strong> proposed once for the cohort. A cohort constant "
            f"repeats down its column on purpose: two libraries can only be compared by someone "
            f"who can see which numbers differ because the libraries do.</p>")
        head = ("<tr><th>library</th>"
                + "".join(f"<th>{_e(c['label'])}</th>" for c in ps["columns"]) + "</tr>")
        scope = ("<tr><td class='src'>scope</td>"
                 + "".join(f"<td class='src'>{_e(c['scope'])}</td>" for c in ps["columns"])
                 + "</tr>")
        body = "".join(
            "<tr><td><strong>" + _e(r["sample"]) + "</strong></td>"
            + "".join(
                "<td class='"
                + ("ns" if r["cells"][c["key"]] == NOT_STATED else "")
                + f"'>{_e(r['cells'][c['key']])}</td>" for c in ps["columns"])
            + "</tr>" for r in ps["rows"])
        parts.append("<div class='scroll'><table class='wide'>" + head + scope + body
                     + "</table></div>")
        constant = [c for c in ps["columns"] if c["key"] not in ps["varying"]]
        if constant:
            parts.append("<p class='src'>identical in every library: "
                         + _e(", ".join(c["label"] for c in constant)) + "</p>")
        parts.append(f"<p class='src'>source: {_e(ps['source'])}</p>")

    # ---------------------------------------------------------------- 5 open items
    parts.append("<h2>5 · Open items</h2>")
    if not doc["open_items"]["stated"]:
        parts.append("<div class='defect'><strong>Defect.</strong> No open-items list was "
                     "supplied. A missing section and an empty one are not the same claim: this "
                     "report cannot say that nothing is unresolved.</div>")
    elif not doc["open_items"]["items"]:
        parts.append("<p>none</p>")
    else:
        rows = "".join(f"<tr><td>{_e(i['item'])}</td><td>{_e(i['closes_when'])}</td>"
                       f"<td>{_e(i['blocked_on'])}</td></tr>"
                       for i in doc["open_items"]["items"])
        parts.append("<table><tr><th>open item</th><th>what would close it</th>"
                     "<th>blocked on</th></tr>" + rows + "</table>")

    # ---------------------------------------------------------------- defects appendix
    if doc["defects"]:
        parts.append("<h2>Defects and notices in this report</h2>")
        rows = "".join(f"<tr><td>{_e(d['severity'])}</td><td class='src'>{_e(d['where'])}</td>"
                       f"<td>{_e(d['what'])}</td></tr>" for d in doc["defects"])
        parts.append("<table><tr><th>kind</th><th>where</th><th>what is missing</th></tr>"
                     + rows + "</table>")

    parts.append("<footer>Self-contained: every figure is embedded, no request leaves this "
                 "file. Every number here is also in report.json beside it.</footer>")
    parts.append("</div></body></html>")
    return "".join(parts)


# ------------------------------------------------------------------------------------- write


def build_report(payload: dict, out_html, out_json, *, now=None, render_figures_flag: bool = True
                 ) -> dict:
    """Write the report and its JSON companion, and report what was actually written.

    Returns `{"outputs": [...], "metrics": {...}, "versions": {...}}` where `outputs` lists only
    files that exist, are non-empty, AND carry the bytes this call produced. A partial payload
    does not stop this: a run that halted at step 2 still produces the document, because the
    document is the only record of why it halted.

    NO PREVIOUS RUN'S FILE IS EVER ACCEPTED AS THIS RUN'S OUTPUT

    "The output exists afterwards" is satisfied by a leftover from an earlier run, so a build
    that wrote nothing would be recorded as a success and the old document - old numbers, old
    parameters, old verdict - would be carried forward under the new ones. Two things prevent it
    here. Both declared outputs are UNLINKED before either is written, so a failure between the
    two cannot leave a fresh HTML beside a previous run's JSON; and after the write each file is
    read back and its SHA-256 compared with the bytes this call produced, which is a stronger
    statement than a modification time (no clock resolution, no timezone, no copy that preserves
    mtime can defeat it).

    The bytes are written with `write_bytes`, not `write_text`: on Windows the latter translates
    every `\\n` into `\\r\\n`, so the file on disk would not be the bytes that were hashed and two
    platforms would produce two different documents from one payload.

    `now` pins the generation timestamp, so a test can produce byte-identical output twice.
    Staleness is computed and printed; it does not raise here - see `refuse_if_stale()`, which is
    the enforcing form for a caller about to publish.
    """
    out_html, out_json = Path(out_html), Path(out_json)
    if out_html.resolve() == out_json.resolve():
        raise TaskFailure(f"the HTML and the JSON would be written to the same path "
                          f"({out_html}); one would overwrite the other and the run would keep "
                          f"whichever was written last")

    doc = assemble(payload, now=now)
    if render_figures_flag:
        figures = render_figures(payload, doc["defects"])
    else:
        figures = {"uris": {}, "index": {}, "matplotlib": NOT_INVOKED}
        _defect(doc["defects"], "figures",
                "figure rendering was switched off for this build; every figure block is empty.")
    doc["figures"] = figures["index"]
    for fid, block in figures["index"].items():
        for step in doc["steps"]:
            for f in step["figures"]:
                if f["id"] == fid and not block["rendered"]:
                    f["status"] = block["source"]
    # The defect count and the verdict block are recomputed after the figures, so the front page
    # counts the figure failures too rather than reporting a number that was true before them.
    doc["verdict"]["defect_count"] = sum(1 for d in doc["defects"] if d["severity"] == "defect")
    doc["verdict"]["notice_count"] = sum(1 for d in doc["defects"] if d["severity"] == "notice")
    doc["verdict"]["lines"] = verdict_lines(doc)

    html_text = render_html(doc, figures["uris"])
    assert_self_contained(html_text)

    import hashlib
    blobs = {out_html: html_text.encode("utf-8"),
             out_json: json.dumps(doc, indent=2, ensure_ascii=False,
                                  default=str).encode("utf-8")}
    # Clear BOTH declared outputs before writing EITHER, so nothing that survives to the checks
    # below can be a file this call did not write.
    for path in (out_html, out_json):
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            if path.exists() or path.is_symlink():
                path.unlink()
        except OSError as exc:
            raise TaskFailure(f"the previous {path} could not be removed before writing this "
                              f"one ({exc}); a report that cannot replace its predecessor would "
                              f"be read as this run's output while carrying the last run's "
                              f"numbers") from None
    for path, blob in blobs.items():
        path.write_bytes(blob)

    outputs = []
    for path in (out_html, out_json):
        if not path.exists():
            raise TaskFailure(f"{path} was written and is not on disk afterwards; the report "
                              f"cannot claim an output it did not produce")
        if path.stat().st_size == 0:
            raise TaskFailure(f"{path} is zero bytes. An empty report reads as a run with "
                              f"nothing to say rather than as a failed write")
        on_disk = hashlib.sha256(path.read_bytes()).hexdigest()
        if on_disk != hashlib.sha256(blobs[path]).hexdigest():
            raise TaskFailure(f"{path} is not the file this build wrote: its contents hash to "
                              f"{on_disk[:12]} and the bytes produced here hash to "
                              f"{hashlib.sha256(blobs[path]).hexdigest()[:12]}. Something else "
                              f"is writing to that path, and reporting it as this run's output "
                              f"would attribute another run's numbers to these parameters")
        outputs.append(path)

    fresh = doc["provenance"]["freshness"]
    return {
        "outputs": outputs,
        "metrics": {
            "verdict": doc["verdict"]["overall"],
            "n_refuse": len(doc["verdict"]["refuse"]),
            "n_review": len(doc["verdict"]["review"]),
            "n_findings": doc["verdict"]["n_findings"],
            "n_defects": doc["verdict"]["defect_count"],
            "n_notices": doc["verdict"]["notice_count"],
            "n_parameters": len(doc["parameters"]["rows"]),
            "n_steps": len(doc["steps"]),
            "n_figures_requested": len(figures["index"]),
            "n_figures_rendered": sum(1 for b in figures["index"].values() if b["rendered"]),
            "n_open_items": len(doc["open_items"]["items"]),
            "open_items_stated": doc["open_items"]["stated"],
            "stale": fresh["stale"],
            "freshness_checked": fresh["checked"],
            "html_bytes": out_html.stat().st_size,
            "json_bytes": out_json.stat().st_size,
        },
        "versions": {"matplotlib": figures["matplotlib"], "python": sys.version.split()[0]},
        "document": doc,
    }


def refuse_if_stale(result_or_doc: dict) -> None:
    """Raise Refusal when the report is older than its inputs. The enforcing form.

    Kept out of `build_report` deliberately. Building is the REMEDY for staleness, and a check
    that blocks the command which would fix the problem is worse than no check - it is the one
    that gets switched off. This is called before PUBLISHING an artifact, not before writing one.
    """
    doc = result_or_doc.get("document", result_or_doc)
    fresh = doc.get("provenance", {}).get("freshness", {})
    if fresh.get("stale") is True:
        raise Refusal(
            f"the report is older than its newest input "
            f"(generated {fresh.get('generated')}, newest input {fresh.get('newest_input')}). "
            f"Staleness has no symptom: it opens, renders and reads exactly like a correct "
            f"report. Rebuild it rather than publishing this one.")
    if fresh.get("stale") is None:
        raise Refusal(
            f"freshness could not be checked: {fresh.get('reason', 'no reason recorded')}. "
            f"NOT CHECKED is its own outcome and must not be published as though it were a "
            f"pass; supply the generation time and the newest input time and build again.")
