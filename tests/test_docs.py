# Checks that the published methods document still describes the code it documents.
"""Does docs/FILTERS.md still describe the implementation?

WHY THIS TEST EXISTS

A methods document that has drifted from the code is worse than no document at all. It reads
exactly like a correct one, and it is the thing a user consults *instead of* reading the source —
so a threshold quoted there and changed here is a wrong number with a citation. Every other check
in this repository guards a computation; this one guards the description of it, on the grounds
that a published number nobody can check is the same class of problem.

It also caught its first defect on the day it was written, and not in the document: `find_valley`
declared `min_mode_observations=8` in its signature and told the reader "default 30" in its own
docstring, four lines apart. Whichever a reader believed, one of them was wrong.

THE CODE SIDE IS THE IMPORTED VALUE, NOT ITS SPELLING (harness autodev W020, B1-107). This compared
the document with the LITERAL TEXT of each assignment - `UMI_BOUNDS = {"snrna": (200, 1000), ...}`
as written - and with a second copy of the value typed into this file. Changing a threshold then
took three edits (the module, the document, and this test), and a change of spelling with no
change of value (a reformatted dict, a comment) failed it. Now each value is read from the
imported module, rendered the way the document writes it, and looked for in the document: edit
the constant and the document together and this passes unedited; edit one alone and it fails.

WHAT IT CANNOT DO

It checks the NUMBERS, not the prose. A document can quote every constant correctly and still
describe the wrong procedure, and nothing here would notice. What it makes impossible is the
cheapest and commonest kind of drift: a value edited in one place and not the other.
"""

from __future__ import annotations

import importlib.util
import inspect
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

fails: list[str] = []
DOC = ROOT / "docs" / "FILTERS.md"

print("Documentation - docs/FILTERS.md against the modules it describes")
print("=" * 74)

if not DOC.exists():
    raise SystemExit(f"FAILED - {DOC} does not exist; the README links to it")
doc = DOC.read_text(encoding="utf-8")


def load(rel: str):
    """The module at `rel`, imported the way the pipeline imports its step modules."""
    lib = str(ROOT / "lib")
    if lib not in sys.path:
        sys.path.insert(0, lib)
    spec = importlib.util.spec_from_file_location(f"docs_check_{Path(rel).stem}", ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


Q = load("modules/05_quality/quality.py")
D = load("modules/04_doublets/doublet_health.py")
L = load("modules/03_light_floor/light_floor.py")
C = load("modules/06_cluster_check/cluster_flags.py")
from adapters import scanpy_ops as S  # noqa: E402
from engine import declare  # noqa: E402

#: find_valley's procedure parameters, as its signature declares them.
FV = {k: v.default for k, v in inspect.signature(S.find_valley).parameters.items()
      if v.default is not inspect.Parameter.empty}


def counts(lo, hi) -> str:
    """A count range as the document writes it: **200 – 1,000**."""
    return f"**{int(lo):,} – {int(hi):,}**"


def percents(lo, hi) -> str:
    """A percentage range as the document writes it: **5 – 10%**."""
    return f"**{lo:g} – {hi:g}%**"


def per_assay(table, render) -> bool:
    """Every assay's value is in the document - one assay's alone would bound the other unseen."""
    return bool(table) and all(
        (render(*v) if isinstance(v, tuple) else render(v)) in doc for v in table.values())


#: (what it is, the document states the code's value)
CHECKS = [
    # Per assay since 2026-10-02: the nuclear pair is measured, the whole-cell pair declared,
    # and the document has to carry both or a reader of it bounds cells as nuclei.
    ("UMI bounds", per_assay(Q.UMI_BOUNDS, counts)),
    ("gene bounds", per_assay(Q.GENE_BOUNDS, counts)),
    ("valley spread review", f"**{Q.SPREAD_REVIEW:.1f}×**" in doc),
    ("Tukey multiplier", f"{Q.IQR_MULT:g} x IQR" in doc),
    ("mitochondrial bounds", per_assay(Q.MITO_BOUNDS, percents)),
    # k is DERIVED per cohort where the table holds None and DECLARED where it holds a number, so
    # what must stay true is that the document says which is which, for every assay.
    ("MAD k is per assay, derived or declared",
     callable(getattr(Q, "select_mad_k", None)) and all(
         ("**derived** from the cohort" in doc) if k is None
         else (f'`MITO_MAD_K["{a}"] = {k}`' in doc)
         for a, k in Q.MITO_MAD_K.items())),
    ("MAD k sanity limits", f"`MAD_K_BOUNDS = {tuple(Q.MAD_K_BOUNDS)}`" in doc),
    ("MAD scale factor", f"{Q.MAD_SCALE:g}" in doc),
    ("mito derivation max", f"`MITO_DERIVATION_MAX = {S.MITO_DERIVATION_MAX!r}`" in doc),
    ("light floor default", f"defaulting to **{L.DEFAULT_FLOOR}**" in doc),
    ("doublet silence", per_assay(D.ZERO_RATE, lambda v: f"< {100 * v:g}%")),
    ("doublet rate imposed", f"{D.SPREAD_IMPOSED:g}×" in doc),
    ("doublet rate unstable", f"> {D.SPREAD_UNSTABLE:g}×" in doc),
    ("design differential", per_assay(D.DESIGN_REFUSE, lambda v: f"≥ {v:g}×")),
    ("materiality floor", per_assay(D.MATERIAL, lambda v: f"≥ {100 * v:g}%")),
    ("cluster marker top-N", f"top-{C.TOPN} markers" in doc),
    ("KDE grid size", f"**{FV['grid_size']}**" in doc),
    ("valley depth", f"`min_valley_depth = {FV['min_valley_depth']:.2f}`" in doc),
    ("mode separation",
     f"`min_mode_separation_log10 = {FV['min_mode_separation_log10']:.2f}`" in doc),
    ("bandwidth stability", f"`bw_stability_factor = {FV['bw_stability_factor']:g}`" in doc),
    ("KDE bandwidth rule", f"{str(FV['bw_method']).capitalize()}'s bandwidth" in doc),
]

#: The criteria step 7 applies on the denoised route, from the tool's own declaration.
APPLIED = [c for c in declare.CRITERIA if c in ("fail_not_cellbender_cell", "fail_umi_floor",
                                                "fail_gene_floor", "fail_mito_ceiling",
                                                "fail_doublet")]
CHECKS.append(("applied criteria", len(APPLIED) == 5 and all(c in doc for c in APPLIED)))

for name, ok in CHECKS:
    ok = bool(ok)
    print(f"  {'ok    ' if ok else 'DRIFT '} {name}")
    if not ok:
        fails.append(
            f"{name}: the document does not state the value the module holds. One of them has "
            f"been edited without the other, and a reader cannot tell which.")

# A signature and its own docstring must agree, which is the drift that needs no second file.
sig = FV.get("min_mode_observations")
says = re.search(r"`min_mode_observations` of them \(default (\d+)\)", S.find_valley.__doc__ or "")
same = bool(sig is not None and says and int(says.group(1)) == int(sig))
print(f"  {'ok    ' if same else 'DRIFT '} {'signature vs docstring':<22}"
      f"signature {sig if sig is not None else '?'}   docstring {says.group(1) if says else '?'}")
if not same:
    fails.append(f"find_valley declares min_mode_observations={sig} and its docstring says "
                 f"{says.group(1) if says else '?'}.")

# The README must point at each document, or nobody finds it. And a document the README
# advertises but the repository does not contain is a broken promise on the front page.
readme = (ROOT / "README.md").read_text(encoding="utf-8")
for name in ("docs/QUICKSTART.md", "docs/USER_GUIDE.md", "docs/FILTERS.md", "docs/OUTPUTS.md"):
    exists = (ROOT / name).exists()
    linked = name in readme
    ok = exists and linked
    print(f"  {'ok    ' if ok else 'DRIFT '} {name:<22}exists {exists}   linked {linked}")
    if not ok:
        fails.append(f"{name}: exists={exists}, linked from README={linked}.")

# ---- docs/REPORT_DESIGN.md assigns every figure to a step; the code must put it there.
#
# The design document assigned F10 to step 4, F11 to step 7 and F12 to step 5. `STEPS` in
# report/build.py declared F1..F9 and stopped, and a step displays only the ids it declares - so
# those three were specified, implemented, assembled, RENDERED, and shown to nobody. The report
# was not wrong about them; it never mentioned them, which is the failure a reader cannot detect,
# because an absent figure reads like a figure nobody designed.
DESIGN = ROOT / "docs" / "REPORT_DESIGN.md"
if not DESIGN.exists():
    fails.append(f"{DESIGN} does not exist; it defines the report's figures.")
else:
    from report.build import FIGURE_QUESTIONS, STEPS

    declared = {fid: key for key, _t, _p, _r, ids in STEPS for fid in ids}
    rows = re.findall(r"^\|\s*(F\d+)\s*\|\s*(\d+)\s+[^|]*\|\s*([^|]+?)\s*\|",
                      DESIGN.read_text(encoding="utf-8"), re.MULTILINE)
    if not rows:
        fails.append("REPORT_DESIGN.md: no figure rows matched; this check tests nothing.")
    for fid, step_no, question in rows:
        where = declared.get(fid)
        ok = where is not None and where.startswith(f"{int(step_no):02d}_")
        print(f"  {'ok    ' if ok else 'DRIFT '} {fid:<4} doc: step {step_no:<2} "
              f"code: {where or 'NO STEP DECLARES IT'}")
        if not ok:
            fails.append(f"{fid}: REPORT_DESIGN.md assigns it to step {step_no}; "
                         f"report.build.STEPS declares it at {where or 'no step at all'}. A "
                         f"figure no step declares is never displayed, however well it draws.")
        # Markdown emphasis is formatting, not wording: `*uniquely*` and `uniquely` are the same
        # question, and failing on that would train the reader of this check to ignore it.
        question = re.sub(r"\*+", "", question)
        if fid in FIGURE_QUESTIONS and FIGURE_QUESTIONS[fid].strip() != question.strip():
            fails.append(f"{fid}: the question differs between the document and the code.\n"
                         f"    doc:  {question.strip()}\n"
                         f"    code: {FIGURE_QUESTIONS[fid].strip()}")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print(f"documentation OK - {len(CHECKS) + 2} published facts match the code they describe")
