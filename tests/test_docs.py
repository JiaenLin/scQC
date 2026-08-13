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

WHAT IT CANNOT DO

It checks the NUMBERS, not the prose. A document can quote every constant correctly and still
describe the wrong procedure, and nothing here would notice. What it makes impossible is the
cheapest and commonest kind of drift: a value edited in one place and not the other.
"""

from __future__ import annotations

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


def src(rel: str) -> str:
    return (ROOT / rel).read_text(encoding="utf-8")


def const(rel: str, name: str) -> str | None:
    """A module-level constant's literal, as written, with any trailing comment removed."""
    m = re.search(rf"^{name}\s*=\s*(.+?)(?:\s*#.*)?$", src(rel), re.M)
    return m.group(1).strip() if m else None


Q = "modules/05_quality/quality.py"
D = "modules/04_doublets/doublet_health.py"
L = "modules/03_light_floor/light_floor.py"
C = "modules/06_cluster_check/cluster_flags.py"
S = "adapters/scanpy_ops.py"

#: (what it is, is it still that in the code, is it still that in the document)
CHECKS = [
    ("UMI bounds", const(Q, "UMI_BOUNDS") == "(200, 1000)", "**200 – 1,000**" in doc),
    ("gene bounds", const(Q, "GENE_BOUNDS") == "(100, 600)", "**100 – 600**" in doc),
    ("valley spread review", const(Q, "SPREAD_REVIEW") == "2.0", "**2.0×**" in doc),
    ("Tukey multiplier", const(Q, "IQR_MULT") == "1.5", "1.5 x IQR" in doc),
    ("mitochondrial bounds",
     const(Q, "MITO_BOUNDS") == '{"snrna": (10.0, 25.0), "scrna": (10.0, 30.0)}',
     "10 – 25%" in doc and "10 – 30%" in doc),
    # k is DERIVED per cohort, so there is no k constant to check - what must stay true is that
    # the document says so, and that the sanity limits on the derived value match the code.
    ("MAD k is derived, not a constant",
     "MAD_K =" not in src(Q) and "def select_mad_k" in src(Q),
     "**derived, not declared**" in doc),
    ("MAD k sanity limits", const(Q, "MAD_K_BOUNDS") == "(2, 10)",
     "`MAD_K_BOUNDS = (2, 10)`" in doc),
    ("MAD scale factor", const(Q, "MAD_SCALE") == "1.4826", "1.4826" in doc),
    ("mito derivation max", const(S, "MITO_DERIVATION_MAX") == "50.0",
     "`MITO_DERIVATION_MAX = 50.0`" in doc),
    ("light floor default", const(L, "DEFAULT_FLOOR") == "200", "defaulting to **200**" in doc),
    ("doublet silence", const(D, "ZERO_RATE") == "0.005", "< 0.5%" in doc),
    ("doublet rate imposed", const(D, "SPREAD_IMPOSED") == "1.05", "1.05×" in doc),
    ("doublet rate unstable", const(D, "SPREAD_UNSTABLE") == "5.0", "> 5×" in doc),
    ("design differential", const(D, "DESIGN_REFUSE") == "3.0", "≥ 3×" in doc),
    ("materiality floor", const(D, "MATERIAL") == "0.01", "≥ 1%" in doc),
    ("cluster marker top-N", const(C, "TOPN") == "20", "top-20 markers" in doc),
    ("KDE grid size", "grid_size=512" in src(S), "**512**" in doc),
    ("valley depth", "min_valley_depth=0.10" in src(S), "`min_valley_depth = 0.10`" in doc),
    ("mode separation", "min_mode_separation_log10=0.30" in src(S),
     "`min_mode_separation_log10 = 0.30`" in doc),
    ("bandwidth stability", "bw_stability_factor=1.5" in src(S),
     "`bw_stability_factor = 1.5`" in doc),
    ("KDE bandwidth rule", 'bw_method="scott"' in src(S), "Scott's bandwidth" in doc),
]

APPLIED = ("fail_not_cellbender_cell", "fail_umi_floor", "fail_gene_floor",
           "fail_mito_ceiling", "fail_doublet")
CHECKS.append(("applied criteria", all(c in src(S) for c in APPLIED),
               all(c in doc for c in APPLIED)))

for name, in_code, in_doc in CHECKS:
    ok = bool(in_code) and bool(in_doc)
    print(f"  {'ok    ' if ok else 'DRIFT '} {name:<22}"
          f"code {'yes' if in_code else 'NO '}   document {'yes' if in_doc else 'NO '}")
    if not ok:
        fails.append(
            f"{name}: the module says {'yes' if in_code else 'something else'} and the document "
            f"says {'yes' if in_doc else 'something else'}. One of them has been edited without "
            f"the other, and a reader cannot tell which.")

# A signature and its own docstring must agree, which is the drift that needs no second file.
sig = re.search(r"min_mode_observations=(\d+)", src(S))
says = re.search(r"`min_mode_observations` of them \(default (\d+)\)", src(S))
same = bool(sig and says and sig.group(1) == says.group(1))
print(f"  {'ok    ' if same else 'DRIFT '} {'signature vs docstring':<22}"
      f"signature {sig.group(1) if sig else '?'}   docstring {says.group(1) if says else '?'}")
if not same:
    fails.append(f"find_valley declares min_mode_observations="
                 f"{sig.group(1) if sig else '?'} and its docstring says "
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
