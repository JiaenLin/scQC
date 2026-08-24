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
import re
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
     "score per sample, before quality filtering; flag only", "nothing", ("F5", "F10")),
    ("05_quality", "5 · quality",
     "derive count floors and mitochondrial ceiling", "nothing",
     ("F6", "F13", "F7", "F12", "F14", "F15")),
    ("06_cluster_check", "6 · cluster check",
     "per-cluster flags: depth, mitochondrial, markers, doublet", "nothing", ("F8",)),
    ("07_apply", "7 · apply", "pre-flight, verify approval, remove",
     "YES - the only step in the pipeline that removes anything", ("F9", "F11")),
)

#: The question each figure answers, from docs/REPORT_DESIGN.md. Printed with the figure so a
#: reader knows what it is for before deciding whether it answers them.
FIGURE_QUESTIONS = {
    "F1": "is this really raw, unfiltered input?",
    "F2": "how much was removed, from what, and evenly?",
    "F3": "did the denoiser drop cells the aligner kept?",
    "F4": "what was never examined?",
    "F5": "is the rate a measurement or the prior?",
    "F6": "where is the UMI cut and why there?",
    "F7": "what did the cut change?",
    "F8": "are any clusters technical?",
    "F9": "what did each criterion remove uniquely?",
    "F10": "where in the manifold did the doublets sit?",
    "F11": "did the removed nuclei leave as a population, or scattered?",
    "F12": "the same count distributions, on the scale people work in",
    "F13": "where is the GENE cut and why there? - step 5 derives two floors and applies both",
    "F14": "what did the GENE floor change?",
    "F15": "what did the MITOCHONDRIAL ceiling change, library by library?",
    # THE CONFOUNDING BLOCK. These three belong to no step: they describe the DESIGN and what the
    # filter did to it, and a pipeline cannot change either by running better.
    "F16": "what is confounded with what, before any number is read?",
    "F17": "how far apart do the confounded arms sit on each QC metric - further than the "
           "libraries inside an arm?",
    "F18": "did the filter remove the same share from each confounded arm?",
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


#: `finding 7 (per_setting (lib1))` and `finding 12 (per_setting (lib7))` are ONE fact about the
#: pipeline, written once per library. The index and the library differ; neither is the cause.
_FINDING_INSTANCE = re.compile(r"finding \d+ \(([^()]+?)(?: \(([^()]+)\))?\)")


def group_defects(defects: list) -> list:
    """One entry per distinct CAUSE, carrying how many times it occurred and in which libraries.

    WHY THE COUNT CHANGED, AND WHY IT IS NOT A WEAKENING

    Every rule here fires per occurrence, and several occurrences are per library. One structural
    fact - "this task declares no output file" - therefore arrived as thirty defects on a
    ten-library cohort, one per (metric x library). A reader seeing 44 assumes 44 problems; there
    were 17, and the three largest were one cause each. A count nobody can act on is a count nobody
    reads, and this document's whole claim is that its own gaps are legible.

    NOTHING IS HIDDEN. Each group keeps every instance, so the per-library detail is still in
    `report.json`; what changes is that the headline counts causes and the table carries `x10`
    beside the one line describing them.

    Order is preserved - first occurrence wins - because defects are emitted in the order the
    document is assembled, and that order is itself information about where the run thinned out.
    """
    grouped: dict = {}
    for d in defects:
        what = str(d.get("what", ""))
        key = (str(d.get("severity", "")), str(d.get("where", "")),
               _FINDING_INSTANCE.sub(lambda m: f"finding <{m.group(1)}>", what))
        m = _FINDING_INSTANCE.search(what)
        who = m.group(2) if m and m.group(2) else None
        slot = grouped.get(key)
        if slot is None:
            grouped[key] = {"severity": key[0], "where": key[1], "what": key[2],
                            "count": 1, "instances": [who] if who else []}
        else:
            slot["count"] += 1
            if who:
                slot["instances"].append(who)
    return list(grouped.values())


def _regroup(entries: list) -> list:
    """Group a list that already holds grouped entries mixed with raw ones.

    `render_figures` runs after `assemble` and appends raw defects to the already-grouped list, so
    what reaches the front page is a mixture. Passing that back through `group_defects` alone
    would reset every carried count to 1 - the figure failures counted right and everything before
    them undercounted, which is worse than either mistake alone because the total still looks
    plausible. This carries existing counts forward.
    """
    out: dict = {}
    for e in entries:
        n = int(e.get("count", 1) or 1)
        one = group_defects([e])[0]
        key = (one["severity"], one["where"], one["what"])
        slot = out.get(key)
        if slot is None:
            out[key] = {**one, "count": n,
                        "instances": list(e.get("instances") or one["instances"])}
        else:
            slot["count"] += n
            slot["instances"].extend(e.get("instances") or one["instances"])
    return list(out.values())


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
        # NOT a defect for `report` itself. Every run appends a `report` step, so this fired
        # on every run that ever completed - the document noticing its own last section. A
        # notice that is raised by correct behaviour teaches a reader to skip the list it is
        # in, which costs more than the one it would ever catch. An UNEXPECTED extra key is
        # still worth saying, so only the known one is exempt.
        if key == "report":
            continue
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
    # ONLY WHEN A DECISIONS FILE WAS USED. The file is optional by design: a run with none
    # applies what the pipeline derived and records every threshold as DERIVED, which is a
    # complete and honest account. Raising a notice because the operator overrode nothing
    # fires on correct behaviour on most runs. What IS worth saying is a decisions file that
    # was read and not hashed, because then two different files are one run.
    if _stated(_get(decisions, "path")) and not _stated(_get(decisions, "hash")):
        _defect(defects, "section 4 · provenance",
                "a decisions file was used and not hashed. Two runs with the same data and "
                "the same decisions file are the same run; a command line does not have that "
                "property.",
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


def build_removal_breakdown(payload, defects: list) -> dict:
    """How many observations each criterion removed, per library. The numbers behind F9.

    F9 draws the cohort, and the cohort is the wrong unit for the only question this section is
    asked: did a criterion fall evenly across the libraries? A criterion taking 2% of one and 42%
    of another has put a technical gradient exactly where the biology is measured, and a single
    bar cannot show it - nor can a bar be read off to a number, which is what a methods section
    needs.

    `fired` and `sole` are both carried. The totals overlap and sum to more than the removal, so
    the first says how much a criterion touched and only the second says whether it did any work
    of its own. A criterion with a large total and no sole removals changed nothing.

    An absent block is a notice, not a defect: a measure-only run has no removal to break down,
    and a report that called that a defect would train its reader to ignore the section.
    """
    block = _get(payload, "removal_breakdown")
    if block is MISSING or _unknown(block):
        _defect(defects, "what each criterion removed",
                "no per-library removal breakdown was supplied, so how much each criterion "
                "removed from each library cannot be read off this report. Expected when the run "
                "measured and did not apply.", severity="notice")
        return {"stated": False, "criteria": [], "samples": [], "rows": [],
                "source": NOT_STATED}

    criteria = [str(c) for c in _as_list(_get(block, "criteria")) if _stated(c)]
    samples = [str(s) for s in _as_list(_get(block, "samples")) if _stated(s)]
    rows = []
    for r in _as_list(_get(block, "rows")):
        crit, samp = _get(r, "criterion"), _get(r, "sample")
        if not _stated(crit) or not _stated(samp):
            _defect(defects, "what each criterion removed",
                    "a breakdown row names no criterion or no library, so its counts belong to "
                    "nothing.")
            continue
        cells = {k: (NOT_STATED if (_get(r, k) is MISSING or _unknown(_get(r, k)))
                     else _get(r, k))
                 for k in ("n_in", "n_fired", "n_sole", "n_removed_any",
                           "pct_of_library", "pct_sole_of_library")}
        rows.append({"criterion": str(crit), "sample": str(samp), **cells})

    if not rows:
        _defect(defects, "what each criterion removed",
                "the breakdown has no rows. An empty table and an absent one read the same on "
                "the page and are not the same claim.")
    source = _get(block, "source")
    if not _stated(source):
        _defect(defects, "what each criterion removed",
                "the breakdown names no file. Every number in this report must be traceable to "
                "something a reader can open.")

    # The widest ratio any single criterion shows across the libraries, computed here rather than
    # left for a reader to scan sixty cells for. It is the number this table exists to surface.
    worst, worst_crit = None, None
    for c in criteria:
        vals = [r["pct_of_library"] for r in rows
                if r["criterion"] == c and r["sample"] != "ALL"
                and isinstance(r["pct_of_library"], (int, float))]
        vals = [v for v in vals if v > 0]
        if len(vals) < 2:
            continue
        ratio = max(vals) / min(vals)
        if worst is None or ratio > worst:
            worst, worst_crit = ratio, c
    return {"stated": True, "criteria": criteria, "samples": samples, "rows": rows,
            "source": _text(source),
            "widest_ratio": None if worst is None else round(worst, 2),
            "widest_ratio_criterion": worst_crit}


#: The three figures of the confounding section, in the order the questions must be asked in.
CONFOUNDING_FIGURES = ("F16", "F17", "F18")


def build_confounding(payload, defects: list) -> dict:
    """Which design factors cannot be told apart in these libraries, and the evidence for it.

    THIS SECTION IS NOT ABOUT THE RUN. Everything else in this document reports something the
    pipeline did and could have done differently. A design in which two factors partition the
    libraries identically was fixed before the first library was made, and no threshold, no
    correction and no rerun changes it - so the section states it as a property of the experiment
    and offers no remedy, because there is none at this end of the process.

    The relationships are read out of F16's OWN DATA rather than recomputed here, for the reason
    the removal breakdown gives: a section that derives its own numbers can disagree with the
    figure beside it and nothing on the page would say which was right. One computation - in
    `report/collect.design_relations` - and two presentations of it.

    An absent figure here is a NOTICE and not a defect. The commonest reason for F17 and F18 to
    be absent is that no two factors are aliased, which is a clean design and the best possible
    outcome; recording it as a defect would train a reader to ignore the count.
    """
    where = "section 2 · confounding"
    payload_figs = _get(payload, "figures")
    fig_notes = _get(payload, "figure_notes")

    blocks = []
    for fid in CONFOUNDING_FIGURES:
        spec = MISSING if payload_figs is MISSING else _get(payload_figs, fid)
        block = {"id": fid, "question": FIGURE_QUESTIONS.get(fid, ""), "caption": "",
                 "source": "", "status": ""}
        if spec is MISSING or not isinstance(spec, dict):
            reason = "" if fig_notes is MISSING else _text(_get(fig_notes, fid), "")
            block["status"] = ("NOT PRODUCED — " + reason if reason else
                               "NOT PRODUCED — no data was supplied for this figure")
            _defect(defects, where, f"figure {fid} was not assembled.", severity="notice")
        else:
            block["caption"] = _text(_get(spec, "caption"), "")
            block["source"] = _text(_get(spec, "source"), "")
            if not block["source"]:
                _defect(defects, where, f"figure {fid} names no source for the data it draws.",
                        severity="notice")
        blocks.append(block)

    relations, factors, arms = [], [], []
    f16 = MISSING if payload_figs is MISSING else _get(payload_figs, "F16")
    data = MISSING if f16 is MISSING or not isinstance(f16, dict) else _get(f16, "data")
    if isinstance(data, dict):
        relations = [r for r in _as_list(_get(data, "relations")) if isinstance(r, dict)]
        factors = [str(f) for f in _as_list(_get(data, "factors"))]
        arm_map = _get(data, "arms")
        arms = sorted(arm_map) if isinstance(arm_map, dict) else []

    aliased = [r for r in relations if str(r.get("kind")) == "aliased"]
    nested = [r for r in relations if str(r.get("kind")) == "nested"]
    rows = [{"a": str(r.get("a")), "b": str(r.get("b")), "kind": str(r.get("kind")),
             "detail": _text(_get(r, "detail"), "")}
            for r in sorted(relations, key=lambda r: (
                {"aliased": 0, "nested": 1}.get(str(r.get("kind")), 2),
                str(r.get("a")), str(r.get("b"))))]

    if not factors:
        headline = ("No design factor was discovered in the samplesheet, so nothing here can be "
                    "confounded. A column is not a factor if it is constant, if it has more "
                    "levels than the libraries can support, or if it has a single library per "
                    "level.")
        severity = "ok"
    elif aliased:
        pairs = ", ".join(f"{r['a']} = {r['b']}" for r in rows if r["kind"] == "aliased")
        headline = (f"{len(aliased)} pair(s) of design factors partition these libraries "
                    f"identically: {pairs}. No analysis of these data can attribute a difference "
                    f"to one of an aliased pair rather than the other. This is a property of the "
                    f"experiment; nothing downstream can separate them.")
        severity = "REFUSE"
    elif nested:
        headline = (f"No pair is fully aliased, but {len(nested)} pair(s) are nested: one factor "
                    f"is fixed once the other is known, so its effect sits inside the other's "
                    f"and cannot be taken back out.")
        severity = "REVIEW"
    else:
        headline = (f"None of the {len(factors)} design factors is aliased with or nested in "
                    f"another: every pair is crossed, and each factor's effect can be estimated "
                    f"separately from the others.")
        severity = "ok"

    return {"stated": bool(relations) or bool(factors),
            "factors": factors, "arms": arms, "rows": rows, "figures": blocks,
            "n_aliased": len(aliased), "n_nested": len(nested),
            "headline": headline, "severity": severity,
            "source": "the run's samplesheet + tables/<sample>.percell.csv"}


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
        # Before the removal breakdown, because it is the question that decides what the
        # breakdown MEANS: the same removal rate is unremarkable in a crossed design and is the
        # whole finding in an aliased one.
        "confounding": build_confounding(payload, defects),
        "removal_breakdown": build_removal_breakdown(payload, defects),
        "provenance": provenance,
        "open_items": open_items,
        "figures": {},
        # GROUPED BY CAUSE. `defects` is the raw emission order; what the document carries
        # is one entry per cause with its count, because a per-library rule fired ten times
        # is one thing to fix. Every instance is kept inside the group.
        "defects": group_defects(defects),
    }
    doc["verdict"]["defect_count"] = sum(1 for d in doc["defects"]
                                         if d["severity"] == "defect")
    doc["verdict"]["notice_count"] = sum(1 for d in doc["defects"]
                                         if d["severity"] == "notice")
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
                f"(F1..F13). No figure could be resolved from it and every figure block in this "
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
/* The accepted layout design (docs/REPORT_DESIGN.md). Three theme states, not two: an explicit
   choice stamps data-theme, and the default "system" setting stamps nothing - so the bare :root
   carries the complete light palette and both dark rules only redefine tokens. */
:root{
  --ink:#16191d; --paper:#fbfbfc; --panel:#ffffff; --wash:#f1f4f7;
  --rule:#dde2e8; --muted:#5c6570; --spine:#c3ccd6;
  --ok:#0072B2; --review:#E69F00; --refuse:#D55E00; --unknown:#CC79A7;
  --before:#bcc5cf; --after:#2f353c; --median:#0072B2; --cutline:#D55E00;
  --shadow:0 1px 2px rgba(22,25,29,.06), 0 8px 24px -16px rgba(22,25,29,.28);
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --ink:#e6eaef; --paper:#0f1216; --panel:#161a1f; --wash:#1a1f25;
    --rule:#2b323a; --muted:#95a1ad; --spine:#3a434d;
    --ok:#57ABDD; --review:#EDA93B; --refuse:#E8734A; --unknown:#D68FB5;
    --before:#4a545f; --after:#cbd4dd; --median:#57ABDD; --cutline:#E8734A;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
  }
}
:root[data-theme="dark"]{
  --ink:#e6eaef; --paper:#0f1216; --panel:#161a1f; --wash:#1a1f25;
  --rule:#2b323a; --muted:#95a1ad; --spine:#3a434d;
  --ok:#57ABDD; --review:#EDA93B; --refuse:#E8734A; --unknown:#D68FB5;
  --before:#4a545f; --after:#cbd4dd; --median:#57ABDD; --cutline:#E8734A;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
}
*{box-sizing:border-box}
body{margin:0; background:var(--paper); color:var(--ink);
  font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  font-variant-numeric:tabular-nums}
.wrap{max-width:1180px; margin:0 auto; padding:34px 26px 100px}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.lbl{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  font-size:11px; letter-spacing:.13em; text-transform:uppercase; color:var(--muted)}
h1{font-size:19px; margin:0; letter-spacing:-.01em; font-weight:650}
h2{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:12px;
  letter-spacing:.16em; text-transform:uppercase; margin:0 0 4px; font-weight:600}
h3{margin:0 0 6px; font-size:13px; font-weight:620}
p{margin:6px 0}
a{color:var(--ok)}
.sub{color:var(--muted); font-size:13.5px; max-width:70ch}
.masthead{display:flex; justify-content:space-between; align-items:flex-start; gap:20px;
  padding-bottom:16px; border-bottom:2px solid var(--ink); flex-wrap:wrap}
.chip{display:inline-block; padding:3px 10px; border-radius:3px; color:#fff;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  font-size:11px; letter-spacing:.1em; font-weight:700}
.c-REVIEW,.c-review{background:var(--review)}
.c-ok,.c-PASS{background:var(--ok)}
.c-REFUSE,.c-refuse{background:var(--refuse)}
.c-NOTDET,.c-unrecognised{background:var(--unknown)}
.decision{display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:1px;
  background:var(--rule); border:1px solid var(--rule); border-radius:8px; overflow:hidden;
  margin:22px 0 10px}
.cell{background:var(--panel); padding:15px 17px}
.big{font-size:26px; font-weight:640; letter-spacing:-.02em; line-height:1.15; margin:2px 0 0;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.cell .note{font-size:12px; color:var(--muted); margin-top:3px}
.evenness{border-left:3px solid var(--review)}
.spine{position:relative; margin:26px 0 10px}
.spine::before{content:""; position:absolute; left:calc(100% - 86px); top:0; bottom:0; width:2px;
  background:var(--spine); transform:translateX(-50%)}
.axis{display:grid; grid-template-columns:1fr 172px; align-items:center; gap:16px;
  padding:22px 0; border-top:1px solid var(--rule)}
.axis:first-child{border-top:0}
.panel{background:var(--panel); border:1px solid var(--rule); border-radius:8px;
  padding:13px 15px 11px; box-shadow:var(--shadow); min-width:0}
.panel .cap{font-size:11.5px; color:var(--muted); margin:8px 0 0; line-height:1.45}
.panel .src{font-size:10.5px; color:var(--muted); margin-top:4px;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.cut{justify-self:center; text-align:center; background:var(--paper); padding:8px 4px;
  z-index:2; position:relative}
.cut .val{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:15px;
  font-weight:700; background:var(--ink); color:var(--paper); padding:5px 11px;
  border-radius:4px; display:inline-block; white-space:nowrap}
.cut .how{font-size:10.5px; color:var(--muted); margin-top:6px; letter-spacing:.04em}
img{display:block; width:100%; height:auto; max-width:100%}
.missing{border:1px dashed var(--unknown); border-radius:6px; padding:11px 13px;
  background:color-mix(in srgb, var(--unknown) 7%, transparent); font-size:12.5px;
  color:var(--ink)}
.missing b{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:11px;
  letter-spacing:.1em}
section{margin:44px 0 0}
.head{border-bottom:1px solid var(--rule); padding-bottom:7px; margin-bottom:16px}
.card{background:var(--panel); border:1px solid var(--rule); border-radius:8px;
  padding:15px 17px; box-shadow:var(--shadow)}
.grid2{display:grid; grid-template-columns:repeat(auto-fit,minmax(400px,1fr)); gap:14px}
.scroll{overflow-x:auto; border:1px solid var(--rule); border-radius:8px}
table{border-collapse:collapse; width:100%; font-size:12.5px}
th,td{padding:7px 11px; text-align:right; white-space:nowrap; border-bottom:1px solid var(--rule)}
th{background:var(--wash); font-weight:600; font-size:11px; letter-spacing:.06em;
  text-transform:uppercase; color:var(--muted)}
th .scope{display:block; font-weight:400; letter-spacing:.02em; text-transform:none;
  font-size:10px; opacity:.85}
td:first-child,th:first-child{text-align:left;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
tbody tr:last-child td{border-bottom:0}
tbody tr:hover td{background:var(--wash)}
.flag{color:var(--review); font-weight:700}
details{border:1px solid var(--rule); border-radius:7px; background:var(--panel); margin-bottom:8px}
summary{cursor:pointer; padding:10px 14px; font-size:13px; display:flex; gap:10px;
  align-items:center; list-style:none}
summary::-webkit-details-marker{display:none}
summary::after{content:"\25B8"; margin-left:auto; color:var(--muted)}
details[open] summary::after{content:"\25BE"}
details .body{padding:0 14px 13px; font-size:13px; color:var(--muted); max-width:82ch}
details .body ul{margin:6px 0 0 16px; padding:0}
summary:focus-visible{outline:2px solid var(--ok); outline-offset:-2px}
.limits{display:grid; grid-template-columns:repeat(auto-fit,minmax(258px,1fr)); gap:12px}
.limit{border-left:3px solid var(--unknown); padding:2px 0 2px 13px}
.limit b{display:block; font-size:12.5px; margin-bottom:2px}
.limit span{font-size:12.5px; color:var(--muted)}
.prov{display:grid; grid-template-columns:auto 1fr; gap:5px 20px; font-size:12.5px}
.prov dt{color:var(--muted)}
.prov dd{margin:0; font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
  overflow-wrap:anywhere}
.klass{font-weight:700; letter-spacing:.03em; font-size:11px;
  font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.k-ADJUDICATED{color:var(--refuse)} .k-DERIVED{color:var(--ok)}
.k-DECLARED{color:var(--review)} .k-FIXED{color:var(--muted)}
.defect{border-left:3px solid var(--refuse); padding:8px 0 8px 13px; margin:10px 0;
  font-size:13px}
footer{margin-top:52px; padding-top:16px; border-top:1px solid var(--rule); font-size:12px;
  color:var(--muted)}
@media (max-width:860px){
  .axis{grid-template-columns:1fr; gap:12px}
  .spine::before{display:none}
  .cut{justify-self:start}
}
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


def _fig_panel(block, uri, *, title=None) -> str:
    """One figure in the layout, or a stated reason there is none.

    A figure that could not be produced is drawn as a bordered note carrying the reason the
    assembler gave, never as a gap. A gap in a document reads as a figure nobody thought was
    needed, which is the one reading that is never true here.
    """
    head = f"<h3>{_e(title)}</h3>" if title else ""
    if block is None:
        return (f"<div class='panel'>{head}<div class='missing'><b>NOT PRODUCED</b><br>"
                f"this figure is not declared by any step, so nothing assembled it.</div></div>")
    fid = block.get("id", "")
    question = block.get("question") or ""
    cap = block.get("caption") or ""
    src = block.get("source") or ""
    if not uri:
        status = block.get("status") or "NOT PRODUCED"
        return (f"<div class='panel'>{head or ''}<div class='missing'>"
                f"<b>{_e(fid)} · NOT PRODUCED</b><br>{_e(status)}</div>"
                f"<p class='cap'>{_e(question)}</p></div>")
    return (f"<div class='panel'>{head}"
            f"<p class='lbl'>{_e(fid)} · {_e(question)}</p>"
            f"<img src='{uri}' alt='{_e(fid)}: {_e(question)}'>"
            + (f"<p class='cap'>{_e(cap)}</p>" if cap else "")
            + (f"<p class='src'>source: {_e(src)}</p>" if src else "")
            + "</div>")


def _confounding_section(block, uri) -> str:
    """The confounding block: one sentence, the pairs that produced it, and three figures.

    Placed immediately after the decision strip, because it decides what every number below it
    means. The evenness cell of that strip already asks whether the filter fell evenly across the
    design; this section is the same question asked of the design itself, and the answer bounds
    what any of the rest can be used to claim.

    The sentence, the table and the figures all come from one computation - see
    `build_confounding`. Nothing is derived here.
    """
    if not block or not block.get("stated"):
        return ""
    chip = {"REFUSE": "refuse", "REVIEW": "review"}.get(str(block.get("severity")), "ok")
    rows = block.get("rows") or []
    table = ""
    if rows:
        # The aliased rows carry the `flag` class, so the pairs that make a claim impossible are
        # the ones a reader's eye lands on rather than the ones that happen to sort first.
        def _row(r):
            cls = " class='flag'" if r["kind"] == "aliased" else ""
            return ("<tr>"
                    f"<td><span class='mono'>{_e(r['a'])}</span></td>"
                    f"<td><span class='mono'>{_e(r['b'])}</span></td>"
                    f"<td{cls}>{_e(r['kind'])}</td>"
                    f"<td>{_e(r['detail'])}</td></tr>")

        trs = "".join(_row(r) for r in rows)
        table = ("<div class='scroll'><table><thead><tr><th>factor</th><th>factor</th>"
                 "<th>relationship</th><th>how it was determined</th></tr></thead>"
                 f"<tbody>{trs}</tbody></table></div>")
    arms = block.get("arms") or []
    arms_line = ""
    if arms:
        shown = " · ".join("<span class='mono'>"
                           + _e(str(a).replace("\n", " / ")) + "</span>" for a in arms)
        arms_line = ("<p class='sub' style='margin-top:9px'>The arms F17 and F18 compare, one "
                     "per distinct combination of the aliased factors: " + shown + ".</p>")
    panels = "".join(_fig_panel(b, uri(b["id"])) for b in (block.get("figures") or []))
    return ("<section><div class='head'><h2>Where the design is confounded</h2>"
            "<p class='sub'>Computed from the samplesheet by exact comparison of the partitions "
            "each factor induces over the libraries. Not a statistic and not a test: a "
            "confounded factor cannot be tested, which is why it is reported here instead.</p>"
            "</div>"
            f"<div class='card'><p><span class='chip c-{chip}'>{_e(chip)}</span> "
            f"{_e(block.get('headline', ''))}</p>{table}{arms_line}"
            f"<p class='src'>source: {_e(str(block.get('source', '')))}</p></div>"
            f"{panels}</section>")


def _cut(value: str, label: str, how: str) -> str:
    """The threshold between a before and an after, with how it was arrived at."""
    return (f"<div class='cut'><span class='val'>{_e(value)}</span>"
            f"<p class='how'>{_e(label)}<br>{_e(how)}</p></div>")


def _removal_breakdown_table(block) -> str:
    """The numbers behind F9: one row per library, one column pair per criterion.

    Laid out library-by-row rather than criterion-by-row because the comparison a reader makes is
    ACROSS libraries within one criterion, and that comparison has to be a column to be made by
    eye. The cohort total is the last row, marked, so it cannot be mistaken for an eleventh
    library.

    Each criterion gets two figures: how many it fired on, and how many it removed on its own.
    Printing only the first makes every criterion look load-bearing; only the second makes a
    criterion that agrees with its neighbours look inert.

    Returns "" when there is nothing to show. A measure-only run has no removal to break down and
    a table of zeroes would say every criterion removed nothing, which is a claim about a removal
    that did not happen.
    """
    if not block or not block.get("stated"):
        return ""
    criteria = block.get("criteria") or []
    rows = block.get("rows") or []
    if not criteria or not rows:
        return ""
    by = {(r.get("sample"), r.get("criterion")): r for r in rows}
    samples = list(block.get("samples") or [])
    order = samples + (["ALL"] if any(s == "ALL" for s, _ in by) else [])

    head = ("<tr><th rowspan='2'>library</th><th rowspan='2'>in</th>"
            "<th rowspan='2'>removed</th>"
            + "".join(f"<th colspan='2'>{_e(c)}</th>" for c in criteria) + "</tr>"
            + "<tr>" + "".join("<th>fired</th><th>alone</th>" for _ in criteria) + "</tr>")

    trs = []
    for s in order:
        first = next((by[(s, c)] for c in criteria if (s, c) in by), None)
        if first is None:
            continue
        cls = " class='flag'" if s == "ALL" else ""
        tds = [f"<td{cls}><b>{_e(s)}</b></td>" if s == "ALL" else f"<td>{_e(s)}</td>",
               f"<td{cls}>{_fmt(first.get('n_in'))}</td>",
               f"<td{cls}>{_fmt(first.get('n_removed_any'))}</td>"]
        for c in criteria:
            r = by.get((s, c)) or {}
            pct = r.get("pct_of_library")
            pct_txt = (f" <span class='lbl'>{float(pct):.1f}%</span>"
                       if isinstance(pct, (int, float)) else "")
            tds.append(f"<td{cls}>{_fmt(r.get('n_fired'))}{pct_txt}</td>")
            tds.append(f"<td{cls}>{_fmt(r.get('n_sole'))}</td>")
        trs.append("<tr>" + "".join(tds) + "</tr>")

    ratio = block.get("widest_ratio")
    note = ("one criterion's removal rate varies "
            f"<b>{ratio:.2f}×</b> across the libraries "
            f"(<span class='mono'>{_e(str(block.get('widest_ratio_criterion') or ''))}</span>) — "
            "a filter that falls harder on one arm of the design puts a technical gradient where "
            "the biology is measured"
            if isinstance(ratio, (int, float)) else
            "no criterion fired on two or more libraries, so no rate could be compared across "
            "them")
    return ("<div class='head' style='margin-top:22px'><h3>The same, counted</h3>"
            "<p class='sub'><b>fired</b> is every observation the criterion removed and the "
            "columns overlap, so they sum to more than the removal. <b>alone</b> is the ones no "
            "other criterion would have removed — a criterion with a large total and none of "
            "these changed nothing.</p></div>"
            f"<div class='scroll'><table><thead>{head}</thead><tbody>{''.join(trs)}"
            "</tbody></table></div>"
            f"<p class='sub' style='margin-top:9px'>{note}. Source: "
            f"<span class='mono'>{_e(str(block.get('source', '')))}</span></p>")


def _column(per_sample, key):
    """One column of the per-library table, as a list of values.

    A row is `{"sample": ..., "cells": {key: value}}` - the values are nested, and reading them
    off the row itself returns None for every column, which renders as an empty table and as a
    threshold of "not recorded" beside a figure that plainly shows one.
    """
    return [(r.get("cells") or {}).get(key)
            for r in (per_sample.get("rows") or []) if isinstance(r, dict)]


def _numbers(values) -> list:
    """Only the values that are numbers. NOT STATED is a sentinel string, not a quantity."""
    out = []
    for v in values:
        try:
            out.append(float(v))
        except (TypeError, ValueError):
            continue
    return out


def _one_value(values):
    """The single value of a cohort constant, or None when the column does not agree.

    A constant that is not constant is not a constant, and printing the first of several
    different numbers as though it applied to every library is the error this prevents.
    """
    seen = set(_numbers(values))
    return next(iter(seen)) if len(seen) == 1 else None


def _fmt(value, suffix="") -> str:
    if value is None:
        return "not recorded"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    return (f"{int(round(f)):,}" if abs(f - round(f)) < 1e-9 and abs(f) >= 10
            else f"{f:g}") + suffix


def _evenness(findings):
    """The worst design-differential ratio any gate measured, and the factors it compared.

    Read from the gates rather than recomputed, so the headline number and the finding that
    raised it cannot disagree. None when no gate measured one - which is stated, not filled in
    with a reassuring number.
    """
    worst, where = None, ""
    for f in findings or []:
        if "design differential" not in str(f.get("check", "")):
            continue
        m = re.search(r"([\d.]+)\s*x", str(f.get("message", "")))
        if not m:
            continue
        ratio = float(m.group(1))
        if worst is None or ratio > worst:
            worst, where = ratio, str(f.get("check", "")).split(":")[-1].strip()
    return worst, where


def render_html(doc: dict, figure_uris: dict) -> str:
    """The single-file document, in the layout of docs/REPORT_DESIGN.md.

    Inline CSS, inline images, no external request of any kind.

    THE SHAPE, AND WHY IT IS NOT A LOG

    The first version of this renderer emitted the steps in order, each with its metrics - a run
    log. It is the natural thing to produce and the wrong thing to read: it puts what the
    pipeline did in front of what a reader needs to decide, which is whether the filter took the
    same proportion out of every arm of the design, and where each cut fell relative to the
    distribution it was cutting. Those two questions now open the document, and the steps are
    below them.
    """
    v = doc["verdict"]
    run, prov = doc["run"], doc["provenance"]
    per_sample = doc.get("per_sample") or {}
    steps = doc.get("steps") or []
    figs = {f["id"]: f for s in steps for f in (s.get("figures") or [])}
    uri = lambda fid: figure_uris.get(fid)                                    # noqa: E731
    d = v.get("deliverable") or {}

    parts = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>scQC report — {_e(Path(str(run['project'])).name)}</title>",
        f"<style>{_CSS}</style></head><body><div class='wrap'>",
    ]

    # ---------------------------------------------------------------- masthead
    n_libs = len(per_sample.get("rows") or []) or len(_column(per_sample, "sample"))
    parts.append(
        "<div class='masthead'><div>"
        f"<h1>Quality control · <span class='mono'>{_e(Path(str(run['project'])).name)}</span>"
        f"</h1>"
        f"<p class='sub'>{n_libs} librar{'y' if n_libs == 1 else 'ies'} · "
        f"{_e(run['mode'])} mode</p></div>"
        f"<div style='text-align:right'><span class='chip c-{_e(v['overall'])}'>"
        f"{_e(v['overall'])}</span>"
        f"<p class='lbl' style='margin-top:8px'>{_e(prov.get('generated', ''))} · "
        f"scQC {_e(doc.get('report_version', ''))}</p></div></div>")

    # ---------------------------------------------------------------- decision strip
    n_in, n_kept = d.get("n_in"), d.get("n_kept")
    pct = d.get("pct_removed")
    ratio, factor = _evenness(v.get("findings"))
    cells = [
        ("Observations in", _fmt(n_in), "before any criterion was applied"),
        ("Kept", _fmt(n_kept),
         (f"{100 - float(pct):.1f}% retained" if pct is not None else "")),
        ("Removed", _fmt(None if n_in is None or n_kept is None else n_in - n_kept),
         "every one named in the ledger"),
    ]
    strip = "".join(
        f"<div class='cell'><p class='lbl'>{_e(lbl)}</p><p class='big'>{_e(val)}</p>"
        f"<p class='note'>{_e(note)}</p></div>" for lbl, val, note in cells)
    strip += ("<div class='cell evenness'><p class='lbl'>Evenness across the design</p>"
              + (f"<p class='big' style='color:var(--review)'>{ratio:.2f}×</p>"
                 f"<p class='note'>widest ratio measured, on {_e(factor)}</p>"
                 if ratio is not None else
                 "<p class='big'>not measured</p><p class='note'>no gate compared removal "
                 "across the design</p>")
              + "</div>")
    parts.append(f"<div class='decision'>{strip}</div>")
    if not _stated(n_in) or not _stated(n_kept):
        parts.append(f"<p class='sub'>{_e(d.get('text') or '')}</p>")
    parts.append("<p class='sub' style='margin-bottom:0'>A filter that falls harder on one arm "
                 "of the design puts a technical gradient where the biology is measured. This "
                 "is the number to read first.</p>")

    # ---------------------------------------------------------------- the design itself
    # BEFORE the pipeline's own work, because it bounds what any of it can be used to claim. The
    # strip above asks whether the filter fell evenly across the design; this asks whether the
    # design can carry the question at all, and the answer is not something a rerun can improve.
    parts.append(_confounding_section(doc.get("confounding") or {}, uri))

    # ---------------------------------------------------------------- the spine
    umi_cut = _one_value(_column(per_sample, "umi_floor_proposed"))
    gene_cut = _one_value(_column(per_sample, "gene_floor_proposed"))
    mito = _numbers(_column(per_sample, "mito_ceiling_pct"))
    mito_txt = f"{min(mito):.2f}–{max(mito):.2f}%" if mito else "not recorded"
    # The mitochondrial rule's CLASS is taken from the parameter table, which reads it off the
    # run, rather than written here as fixed text. It was fixed text - "DERIVED · per library" -
    # and a bound-dominated or declared-k ceiling would have been labelled derived on the one
    # line of the report a reader looks at to see how the threshold was arrived at.
    mito_class = next((str(r.get("class")) for r in ((doc.get("parameters") or {}).get("rows") or [])
                       if str(r.get("name", "")).startswith("mitochondrial")), None)
    rows = [
        ("Counts per nucleus", "F7", f"≥ {_fmt(umi_cut)}", "UMI floor",
         "DERIVED · cohort constant"),
        ("Genes detected", "F14", f"≥ {_fmt(gene_cut)}", "gene floor",
         "DERIVED · cohort constant"),
        ("Mitochondrial %", "F15", f"≤ {mito_txt}", "mito ceiling",
         f"{mito_class} · per library" if mito_class else "class not recorded"),
        ("Where the doublets sat", "F10", "scDblFinder", "doublet call", "DECLARED · per library"),
        ("What survived, same embedding", "F11", "", "", ""),
    ]
    body = []
    for title, fid, val, label, how in rows:
        block = figs.get(fid) if fid else None
        if fid is None:
            block = {"id": "—", "question": "",
                     "status": "no figure is produced for this axis; the cut is stated from "
                               "the per-library table beside it"}
        cut = _cut(val, label, how) if val else "<div class='cut'></div>"
        body.append(f"<div class='axis'>{_fig_panel(block, uri(fid) if fid else None, title=title)}"
                    f"{cut}</div>")
    parts.append(
        "<section><div class='head'><h2>What quality control did</h2>"
        "<p class='sub'>One row per criterion: every library's distribution before the cut and "
        "what survived it, with the threshold on the rule at the right.</p></div>"
        f"<div class='spine'>{''.join(body)}</div></section>")

    # ---------------------------------------------------------------- linear axis + criteria
    parts.append(
        "<section><div class='head'><h2>The same distributions on a linear axis</h2>"
        "<p class='sub'>The log axis shows that two populations exist and where the boundary "
        "is. The linear axis shows the scale people work in, and how little of the retained "
        "distribution sits near the cut.</p></div>"
        f"{_fig_panel(figs.get('F12'), uri('F12'))}</section>")
    parts.append(
        "<section><div class='head'><h2>What each criterion removed</h2>"
        "<p class='sub'>Solid is what only this criterion removed; a criterion that removes "
        "nothing uniquely is doing no work of its own.</p></div>"
        f"{_fig_panel(figs.get('F9'), uri('F9'))}"
        f"{_removal_breakdown_table(doc.get('removal_breakdown') or {})}</section>")

    # ---------------------------------------------------------------- inputs and clusters
    parts.append(
        "<section><div class='head'><h2>What went in, and what was examined</h2>"
        "<p class='sub'>Before any threshold: whether the input was raw, how much ambient "
        "signal was removed and whether evenly, whether the denoiser dropped cells the aligner "
        "kept, and what was never examined at all.</p></div><div class='grid2'>"
        + "".join(_fig_panel(figs.get(f), uri(f)) for f in ("F1", "F2", "F3", "F4"))
        + "</div></section>")
    parts.append(
        "<section><div class='head'><h2>Cluster check</h2>"
        "<p class='sub'>Flags are reported; nothing is removed on them.</p></div>"
        f"<div class='grid2'>{_fig_panel(figs.get('F5'), uri('F5'))}"
        f"{_fig_panel(figs.get('F8'), uri('F8'))}</div></section>")

    # ---------------------------------------------------------------- per-library thresholds
    cols = per_sample.get("columns") or []
    rws = per_sample.get("rows") or []
    if cols and rws:
        head = "<tr><th>library</th>" + "".join(
            f"<th>{_e(c.get('label', c.get('key', '')))}"
            f"<span class='scope'>{_e(c.get('scope', ''))}</span></th>" for c in cols) + "</tr>"
        trs = []
        for r in rws:
            cells_of = r.get("cells") or {}
            clamped = str(cells_of.get("mito_clamped", "")).strip().lower()
            tds = [f"<td>{_e(str(r.get('sample', '')))}</td>"]
            for c in cols:
                key = c.get("key")
                val = cells_of.get(key)
                # The ceiling was decided by the declared bound rather than by this library's
                # own distribution. Marked because the two are different kinds of number in the
                # same column, and nothing else on the row says which this one is.
                mark = (" ▲" if key == "mito_ceiling_pct"
                        and clamped in ("true", "upper", "lower", "1") else "")
                cls = " class='flag'" if mark else ""
                tds.append(f"<td{cls}>{_e('' if val is None else str(val))}{mark}</td>")
            trs.append("<tr>" + "".join(tds) + "</tr>")
        parts.append(
            "<section><div class='head'><h2>Every threshold, per library</h2>"
            "<p class='sub'>Which numbers differ because the libraries differ, and which are "
            "one constant repeated. A cohort constant repeats down its column on purpose.</p>"
            "</div><div class='scroll'><table><thead>" + head + "</thead><tbody>"
            + "".join(trs) + "</tbody></table></div>"
            + (f"<p class='sub' style='margin-top:9px'>▲ a declared bound decided this ceiling, "
               f"not the library's own distribution. Source: "
               f"<span class='mono'>{_e(str(per_sample.get('source', '')))}</span></p>")
            + "</section>")

    # ---------------------------------------------------------------- what was decided
    # build_parameter_table() returns {"stated": bool, "rows": [...]}; the rows are one key
    # INSIDE it, and iterating the block itself yields its key strings.
    params = (doc.get("parameters") or {}).get("rows") or []
    if params:
        rows_p = "".join(
            f"<tr><td>{_e(str(p.get('name', '')))}</td>"
            f"<td>{_e(str(p.get('value', '')))}</td>"
            f"<td><span class='klass k-{_e(str(p.get('class', '')))}'>"
            f"{_e(str(p.get('class', '')))}</span></td>"
            f"<td>{_e(str(p.get('basis', '') or p.get('verbatim', '')))}</td></tr>"
            for p in params)
        parts.append(
            "<section><div class='head'><h2>What was decided, and by whom</h2>"
            "<p class='sub'>DERIVED is computed from this dataset. DECLARED was supplied before "
            "seeing the result. ADJUDICATED means a person decided after seeing it, in their own "
            "words.</p></div><div class='scroll'><table><thead><tr><th>parameter</th>"
            "<th>value</th><th>class</th><th>basis</th></tr></thead><tbody>"
            + rows_p + "</tbody></table></div></section>")

    # ---------------------------------------------------------------- limits
    limits = "".join(
        f"<div class='limit'><b>{_e(s.get('title', s.get('key', '')))}</b>"
        f"<span>{_e(s.get('cannot_establish', ''))}</span></div>"
        for s in steps if s.get("cannot_establish"))
    if limits:
        parts.append(
            "<section><div class='head'><h2>What this run could not establish</h2>"
            "<p class='sub'>Each step states its own limit. An omitted limit reads as no "
            "limit.</p></div>"
            f"<div class='limits'>{limits}</div></section>")

    # ---------------------------------------------------------------- findings
    refuse, review = v.get("refuse") or [], v.get("review") or []
    unrec = v.get("unrecognised") or []
    parts.append(
        f"<section><div class='head'><h2>Findings · {len(review)} review, {len(refuse)} "
        f"refusal{'' if len(refuse) == 1 else 's'}, {v.get('n_ok', 0)} ok</h2></div>")
    if not v.get("gates_supplied"):
        parts.append("<div class='defect'><b>Defect.</b> No gate findings were supplied. This "
                     "document cannot state that the run passed; only that nobody recorded "
                     "whether it did.</div>")
    for f in list(refuse) + list(review) + list(unrec):
        detail = ""
        if f.get("detail"):
            detail = "<ul>" + "".join(f"<li>{_e(str(x))}</li>" for x in f["detail"]) + "</ul>"
        parts.append(
            f"<details open><summary><span class='chip c-{_e(str(f.get('severity', '')))}'>"
            f"{_e(str(f.get('severity', '')))}</span>{_e(str(f.get('check', '')))}</summary>"
            f"<div class='body'>{_e(str(f.get('message', '')))}{detail}"
            f"<p class='src mono'>{_e(str(f.get('step', '')))}</p></div></details>")
    ok_rows = [f for f in (v.get("findings") or [])
               if str(f.get("severity", "")).lower() == "ok"]
    if ok_rows:
        parts.append("<details><summary><span class='chip c-ok'>ok</span>"
                     f"{len(ok_rows)} check(s) that passed</summary><div class='body'><ul>"
                     + "".join(f"<li><b>{_e(str(f.get('check', '')))}</b> — "
                               f"{_e(str(f.get('message', '')))}</li>" for f in ok_rows)
                     + "</ul></div></details>")
    parts.append("</section>")

    # ---------------------------------------------------------------- provenance
    prov_rows = "".join(f"<dt>{_e(str(k))}</dt><dd>{_e(str(val))}</dd>"
                        for k, val in prov.items() if not isinstance(val, (dict, list)))
    parts.append("<section><div class='head'><h2>Provenance</h2></div>"
                 f"<dl class='prov'>{prov_rows}"
                 f"<dt>wall clock</dt><dd>{_e(str(run.get('elapsed_s', '')))} s over "
                 f"{_e(str(run.get('jobs', '')))} job(s)</dd></dl></section>")

    # ---------------------------------------------------------------- open items and defects
    open_items = doc.get("open_items") or []
    if open_items:
        parts.append("<section><div class='head'><h2>Open items</h2></div><div class='card'><ul>"
                     + "".join(f"<li>{_e(str(o))}</li>" for o in open_items)
                     + "</ul></div></section>")
    if doc.get("defects"):
        def _row(x):
            n = int(x.get("count", 1) or 1)
            # The libraries are named where they fit. A cause that fired on three of ten is a
            # different problem from one that fired on all ten, and the count alone cannot say
            # which - so the names are printed while they are short enough to read.
            who = [str(i) for i in (x.get("instances") or []) if i]
            where = _e(str(x.get("where", "")))
            if n > 1:
                where += (f"<br><span class='note'>×{n}"
                          + (f" · {_e(', '.join(who[:6]))}" if who else "")
                          + (" …" if len(who) > 6 else "") + "</span>")
            return (f"<tr><td>{_e(str(x.get('severity', '')))}</td><td>{where}</td>"
                    f"<td>{_e(str(x.get('what', '')))}</td></tr>")

        rows_d = "".join(_row(x) for x in doc["defects"])
        parts.append(
            "<section><div class='head'><h2>Defects and notices in this report</h2>"
            "<p class='sub'>What this document could not state, and why. A report that omits a "
            "required block reads exactly like a complete one. One row per CAUSE: a rule that "
            "fires per library is one thing to fix, and its multiplicity is printed beside it "
            "rather than filling the table.</p></div>"
            "<div class='scroll'><table><thead><tr><th>kind</th><th>where</th>"
            "<th>what is missing</th></tr></thead><tbody>" + rows_d
            + "</tbody></table></div></section>")

    parts.append("<footer>Self-contained: every figure is embedded and no request leaves this "
                 "file. Every number here is also in <span class='mono'>report.json</span> "
                 "beside it.</footer>")
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
    # The confounding figures belong to no step, so they need the same pass over them. Without
    # it a figure that ASSEMBLED and then failed to draw would appear as an empty panel saying
    # only "NOT PRODUCED", with the exception that explains it recorded nowhere on the page.
    fig_holders = [step["figures"] for step in doc["steps"]]
    fig_holders.append((doc.get("confounding") or {}).get("figures") or [])
    for fid, block in figures["index"].items():
        for holder in fig_holders:
            for f in holder:
                if f["id"] == fid and not block["rendered"]:
                    f["status"] = block["source"]
    # The defect count and the verdict block are recomputed after the figures, so the front page
    # counts the figure failures too rather than reporting a number that was true before them.
    # Regrouped as well: `render_figures` appends raw entries to an already-grouped list, and
    # a count taken over the mixture would be neither causes nor occurrences. `group_defects`
    # is idempotent on entries it has already grouped except for the count, which is why the
    # grouped ones carry it forward rather than being recounted from one.
    doc["defects"] = _regroup(doc["defects"])
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
