# Exercises the report assembler and the figure module. Draws and asserts; removes nothing.
"""Report tests: nothing unknown is rendered as a value, and the self-containment check does not
fire on a self-contained page.

Every case here is a defect that was PROVED by execution before it was fixed, and each is written
as the shortest input that produced the wrong output rather than as a description of it:

  A  `label or default` kept a NaN, because a NaN is truthy. F7's y axis read literally
     `nan (log scale, original units)`, and `pandas.NA` did not even get that far - it raised
     TypeError out of the middle of drawing the panel.
  B  `srcset.split(",")` tore a base64 `data:` URI in half at the comma inside it, so the tail
     was reported as an external reference and `assert_self_contained()` REFUSED a page that was
     entirely self-contained. A gate that fires on correct behaviour is one somebody switches off
     (docs/PRINCIPLES.md section 3) - so the case that must keep passing is pinned here beside the
     case that must keep failing.
  C  an unparseable count raised ValueError out of `assemble()` and destroyed the whole report -
     the one artifact that records why a run stopped.
  D  a single string supplied where a list was expected was iterated CHARACTER BY CHARACTER into
     the document, one bullet per letter.

The figure cases need matplotlib and the sentinel cases need numpy and pandas; each block states
what it skipped rather than passing quietly, because a suite that reports success over a check it
never ran is the failure this pipeline exists to prevent.
"""
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from report import build as B                                             # noqa: E402
from report import figures as F                                           # noqa: E402

fails = []
skipped = []


def _try(name):
    try:
        return __import__(name)
    except ImportError:
        skipped.append(name)
        return None


np = _try("numpy")
pd = _try("pandas")
mpl = _try("matplotlib")

print("=" * 74)
print("A. every spelling of unknown is unknown, in both modules")

UNKNOWN = [("None", None), ("float nan", float("nan"))]
if np is not None:
    UNKNOWN += [("numpy.float64 nan", np.float64("nan")),
                ("numpy.float32 nan", np.float32("nan")),
                ("numpy.datetime64 NaT", np.datetime64("NaT", "ns")),
                ("numpy.ma.masked", np.ma.masked)]
if pd is not None:
    UNKNOWN += [("pandas.NA", pd.NA), ("pandas.NaT", pd.NaT)]
# Values that carry information and must NEVER be read as unknown. `0`, `False` and `"False"` are
# here because the cheap spellings of this predicate - `not v`, `v or default` - all swallow them.
KNOWN = [("0", 0), ("False", False), ("empty string", ""), ("'False'", "False"), ("[]", [])]
if np is not None:
    KNOWN += [("numpy array", np.array([1.0, 2.0])), ("numpy.bool_(True)", np.bool_(True))]

for mod in (B, F):
    for label, value in UNKNOWN:
        if mod._unknown(value) is not True:
            fails.append(f"A: {mod.__name__}._unknown({label}) must be True")
    for label, value in KNOWN:
        if mod._unknown(value) is not False:
            fails.append(f"A: {mod.__name__}._unknown({label}) must be False")
print(f" {len(UNKNOWN)} unknown spellings and {len(KNOWN)} known values, "
      f"against report.build and report.figures")
print(f" MISSING is unknown to the assembler: {B._unknown(B.MISSING)}")

print("\n" + "-" * 74)
print("B. an unknown LABEL falls back to the default; it is not printed as the word 'nan'")

for label, value in UNKNOWN + [("whitespace", "   ")]:
    got = F.label_text(value, "DEFAULT")
    if got != "DEFAULT":
        fails.append(f"B: label_text({label}) -> {got!r}, must be the default")
for label, value in (("'0'", "0"), ("'False'", "False"), ("0", 0), ("False", False)):
    if F.label_text(value, "DEFAULT") != str(value):
        fails.append(f"B: label_text({label}) must be kept, not replaced by the default")
print(f" label_text: {len(UNKNOWN) + 1} unknown spellings -> the default; "
      f"'0', 'False', 0 and False are kept")

if mpl is None:
    print(" SKIPPED the drawn-axis cases: matplotlib is not installed")
else:
    LABELLED = (
        ("F5", F.fig_f5_doublet_sweep, {"sweep": {"s": {"x": [0.02, 0.06], "y": [0.01, 0.05]}}},
         "param_label", lambda fig: [a.get_xlabel() for a in fig.get_axes()]),
        ("F6", F.fig_f6_quality_density,
         {"densities": {"s": {"x": [1, 2, 3], "y": [0.1, 0.5, 0.2], "n": 10}}, "cut": 2},
         "metric_label", lambda fig: [a.get_xlabel() for a in fig.get_axes()]),
        ("F7", F.fig_f7_before_after,
         {"distributions": {"s": {"before": [1, 2, 3, 4], "after": [2, 3, 4, 5]}}, "cut": 2},
         "metric_label", lambda fig: [a.get_ylabel() for a in fig.get_axes()]),
        ("F8", F.fig_f8_cluster_flags,
         {"clusters": [{"sample": "s1", "cluster": 3, "umi_frac_of_sample": 0.1,
                        "median_pct_mt": 5.0, "FLAG": True, "WATCH": False}]},
         "x_label", lambda fig: [a.get_xlabel() for a in fig.get_axes()]),
    )
    for fid, fn, kwargs, key, read in LABELLED:
        for label, value in UNKNOWN:
            try:
                fig = fn(**dict(kwargs, **{key: value}))
            except Exception as exc:                                      # noqa: BLE001
                fails.append(f"B: {fid} with {key}={label} raised "
                             f"{type(exc).__name__}: {exc}")
                continue
            for text in read(fig):
                low = text.lower()
                if "nan" in low or "nat" in low or "<na>" in low or not text.strip():
                    fails.append(f"B: {fid} {key}={label} drew the axis label {text!r}")
        fig = fn(**dict(kwargs, **{key: "UMI per nucleus"}))
        if not any("UMI per nucleus" in t for t in read(fig)):
            fails.append(f"B: {fid} dropped a label that WAS supplied")
    print(f" F5, F6, F7 and F8 draw the default for {len(UNKNOWN)} unknown spellings "
          f"and keep a supplied one")
    if pd is not None:
        fig = F.fig_f8_cluster_flags([{"sample": pd.NA, "cluster": float("nan"),
                                       "umi_frac_of_sample": 0.1, "median_pct_mt": 5.0,
                                       "FLAG": True, "WATCH": False}])
        texts = [t.get_text() for t in fig.get_axes()[0].texts]
        if any("nan" in t.lower() or "<na>" in t.lower() for t in texts):
            fails.append(f"B: F8 annotated a point from a blank cell: {texts}")
        print(f" F8 point annotations from blank cells: {texts[0]!r}")

print("\n" + "-" * 74)
print("C. srcset: a base64 data URI is ONE candidate, and it is inline")

SELF_CONTAINED = (
    '<img src="data:image/png;base64,iVBORw0KGgo=" '
    'srcset="data:image/png;base64,iVBORw0KGgo= 1x">',
    '<img srcset="data:image/png;base64,AAA= 1x, data:image/png;base64,BBB= 2x">',
    '<img srcset="data:image/png;base64,AAA=,data:image/png;base64,BBB=">',
    '<img srcset="">',
    '<p>a caption that merely mentions srcset= and url(http://x) is inert text</p>',
)
EXTERNAL = (
    '<img srcset="data:image/png;base64,AAA= 1x, https://cdn.example/x.png 2x">',
    '<img srcset="https://cdn.example/a,b.png 1x">',
    '<img srcset="a.png 1x">',
    '<img src="https://cdn.example/x.png">',
    '<link rel="stylesheet" href="data:text/css,body{}">',
)
for page in SELF_CONTAINED:
    refs = B.external_references(page)
    if refs:
        fails.append(f"C: a self-contained page was refused: {page[:60]} -> {refs}")
    try:
        B.assert_self_contained(page)
    except Exception as exc:                                              # noqa: BLE001
        fails.append(f"C: assert_self_contained raised on {page[:60]}: {exc}")
for page in EXTERNAL:
    if not B.external_references(page):
        fails.append(f"C: an external reference was MISSED: {page[:70]}")
print(f" {len(SELF_CONTAINED)} self-contained pages pass, "
      f"{len(EXTERNAL)} pages that reach outside are still caught")
print(f" _srcset_urls('data:image/png;base64,AAA= 1x, data:image/png;base64,BBB= 2x') -> "
      f"{B._srcset_urls('data:image/png;base64,AAA= 1x, data:image/png;base64,BBB= 2x')}")

print("\n" + "-" * 74)
print("D. an unreadable count is a visible defect, never an exception out of assemble()")

BAD_COUNTS = [("a word", "twelve"), ("a list", [1, 2]), ("an object", object())]
if pd is not None:
    BAD_COUNTS.append(("pandas.NA", pd.NA))
for label, value in BAD_COUNTS:
    try:
        doc = B.assemble({"deliverable": {"n_kept": value, "n_in": 100, "unit": "cells"}})
    except Exception as exc:                                              # noqa: BLE001
        fails.append(f"D: assemble() raised on n_kept={label}: {type(exc).__name__}: {exc}")
        continue
    text = doc["verdict"]["deliverable"]["text"]
    if "% removed)" in text:
        fails.append(f"D: n_kept={label} composed a deliverable line anyway: {text}")
    if not doc["defects"]:
        fails.append(f"D: n_kept={label} recorded no defect")
    print(f" n_kept={label:<10} -> {text[:60]}")
doc = B.assemble({"deliverable": {"n_kept": 0, "n_in": 0, "unit": "cells"}})
if "UNDEFINED" not in doc["verdict"]["deliverable"]["text"]:
    fails.append("D: an input population of zero must state that the percentage is UNDEFINED")
print(f" n_in=0            -> {doc['verdict']['deliverable']['text'][:60]}")

print("\n" + "-" * 74)
print("E. a string where a list was expected is ONE item, not one per character")

if B._as_list("F3") != ["F3"]:
    fails.append(f"E: _as_list('F3') -> {B._as_list('F3')}")
if B._as_list({"a": 1}) != [{"a": 1}]:
    fails.append("E: _as_list(dict) must be one entry, not its keys")
doc = B.assemble({
    "gates": {"step": "05_quality", "check": "floor", "severity": "REVIEW",
              "message": "one finding supplied as a dict", "detail": "one detail, not eight"},
    "parameters": {"name": "min_umi", "value": 350, "class": "DERIVED", "basis": "the valley"},
    "steps": {"key": "05_quality", "status": "ok", "figures": "F6", "sources": "a/b.csv"},
    "open_items": {"item": "one open item", "closes_when": "x", "blocked_on": "y"},
})
detail = doc["verdict"]["review"][0]["detail"]
if detail != ["one detail, not eight"]:
    fails.append(f"E: a detail string was iterated character-wise: {detail[:6]}")
if len(doc["parameters"]["rows"]) != 1:
    fails.append(f"E: one parameter dict became {len(doc['parameters']['rows'])} rows")
if len(doc["open_items"]["items"]) != 1:
    fails.append(f"E: one open item became {len(doc['open_items']['items'])} items")
figs = [f["id"] for s in doc["steps"] if s["key"] == "05_quality" for f in s["figures"]]
if figs != ["F6"]:
    fails.append(f"E: figures='F6' was iterated character-wise: {figs}")
sources = [s["sources"] for s in doc["steps"] if s["key"] == "05_quality"][0]
if sources != ["a/b.csv"]:
    fails.append(f"E: sources='a/b.csv' was iterated character-wise: {sources}")
print(f" detail={detail}, figures={figs}, sources={sources}, "
      f"parameters={len(doc['parameters']['rows'])}, open={len(doc['open_items']['items'])}")

print("\n" + "-" * 74)
print("F. end to end: a whole report, drawn from unknown-laden data, stays self-contained")

if mpl is None:
    print(" SKIPPED: matplotlib is not installed")
else:
    payload = {
        "run": {"project": "regression", "mode": "evidence",
                "started": "2026-08-09T10:00:00+00:00"},
        "deliverable": {"n_kept": 127050, "n_in": 244968, "unit": "nuclei"},
        "gates": [{"step": "05_quality", "check": "floor in bounds", "severity": "ok",
                   "message": "350 UMI", "detail": ["derived from the valley"]}],
        "figures": {
            "F7": {"data": {"distributions": {"s1": {"before": [1, 5, 20, 100],
                                                     "after": [20, 100, 300]}},
                            "cut": 20, "metric_label": float("nan")},
                   "caption": "the metric before and after the cut"},
            "F8": {"data": {"clusters": [{"sample": "s1", "cluster": 3,
                                          "umi_frac_of_sample": 0.1, "median_pct_mt": 5.0,
                                          "FLAG": True, "WATCH": False}]}},
        },
        "provenance": {"generated": "2026-08-09T12:00:00+00:00",
                       "newest_input": "2026-08-09T09:00:00+00:00"},
        "open_items": [],
    }
    with tempfile.TemporaryDirectory() as td:
        res = B.build_report(payload, Path(td) / "r.html", Path(td) / "r.json")
        text = (Path(td) / "r.html").read_text(encoding="utf-8")
    print(f" verdict {res['metrics']['verdict']}, "
          f"{res['metrics']['n_figures_rendered']} of "
          f"{res['metrics']['n_figures_requested']} figures rendered, "
          f"{res['metrics']['html_bytes']:,} bytes")
    if res["metrics"]["n_figures_rendered"] != 2:
        fails.append(f"F: figures did not render: {res['document']['figures']}")
    B.assert_self_contained(text)
    if "nan (log" in text or "nan (linear" in text:
        fails.append("F: the rendered page carries an axis label reading 'nan'")

print("\n" + "=" * 74)
if skipped:
    print(f"NOT RUN against: {', '.join(skipped)} (not installed) - the cases needing them were "
          f"skipped, not passed")
if fails:
    print("FAILED:")
    for x in fails:
        print(" -", x)
    raise SystemExit(1)
print("proved: an unknown label falls back to the default rather than printing itself as the")
print("word 'nan', a base64 data URI inside a srcset no longer makes the self-containment check")
print("refuse a self-contained page while a genuinely external candidate is still caught, an")
print("unreadable count becomes a defect on the page instead of an exception that destroys the")
print("report, and a string supplied where a list was expected is one item and not one per letter")
