# The per-library, per-criterion removal counts: the arithmetic, and that it reaches the page.
"""Does the report say how many observations each criterion removed, from each library?

WHY THIS TEST EXISTS

F9 has drawn the per-criterion removal since the report existed, and it draws the COHORT. The
cohort is the wrong unit for the only question the section is asked - did this criterion fall
evenly across the libraries? - because a criterion taking 2% of one library and 42% of another
sums to an unremarkable cohort bar. That is a technical gradient sitting exactly where the biology
is measured, and it is invisible in the figure by construction.

It is also the failure mode this repository has met most often at the report layer: a section
computed correctly, printed as computed, and never rendered - so the check runs from the counting
function all the way to the HTML rather than stopping at either end.

WHAT IS CHECKED

  1. `n_fired` counts every observation a criterion removed, overlaps included
  2. `n_sole` counts only the ones no other criterion would have removed
  3. the two are different where criteria overlap, and the test data makes them differ
  4. the cohort row is the sum of the libraries and is labelled, not an eleventh library
  5. the widest per-criterion ratio across libraries is computed, since that is the finding
  6. all of it reaches the rendered HTML
  7. an absent breakdown is a NOTICE and renders no table - a measure-only run has none
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.steps import BREAKDOWN_ALL, _removal_breakdown  # noqa: E402
from report import build as B  # noqa: E402

fails: list[str] = []
print("Per-library removal breakdown - engine/steps + report/build")
print("=" * 74)


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'ok    ' if ok else 'FAILED'} {name}{('   ' + detail) if detail else ''}")
    if not ok:
        fails.append(f"{name}: {detail}")


# ---------------------------------------------------------------- 1-4: the arithmetic
#
# Two libraries, deliberately lopsided. In A the mitochondrial ceiling fires on one barcode and
# that barcode ALSO fails the UMI floor, so the ceiling did no work of its own - the case that
# looks load-bearing in a total and is exposed only by the sole count. In B the ceiling fires on
# three and two of them on nothing else.
CRITERIA = ["fail_umi_floor", "fail_mito_ceiling"]
percell = (
    # library A: 5 barcodes, 2 removed, the ceiling never acting alone
    [{"sample": "A"} for _ in range(5)]
    # library B: 5 barcodes, 3 removed, the ceiling acting alone on 2
    + [{"sample": "B"} for _ in range(5)]
)
masks = {
    "fail_umi_floor":    [True, True, False, False, False,   True, False, False, False, False],
    "fail_mito_ceiling": [True, False, False, False, False,  True, True, True, False, False],
}
removed = [a or b for a, b in zip(masks["fail_umi_floor"], masks["fail_mito_ceiling"])]

rows = _removal_breakdown(percell, CRITERIA, masks, removed)
by = {(r["sample"], r["criterion"]): r for r in rows}

check("n_fired counts every observation a criterion removed",
      by[("A", "fail_mito_ceiling")]["n_fired"] == 1
      and by[("B", "fail_mito_ceiling")]["n_fired"] == 3,
      f"A={by[('A', 'fail_mito_ceiling')]['n_fired']}, "
      f"B={by[('B', 'fail_mito_ceiling')]['n_fired']}")

check("n_sole excludes what another criterion also removed",
      by[("A", "fail_mito_ceiling")]["n_sole"] == 0
      and by[("B", "fail_mito_ceiling")]["n_sole"] == 2,
      f"A={by[('A', 'fail_mito_ceiling')]['n_sole']} (the one hit also fails the UMI floor), "
      f"B={by[('B', 'fail_mito_ceiling')]['n_sole']}")

check("fired and sole differ where the criteria overlap",
      by[("A", "fail_mito_ceiling")]["n_fired"] != by[("A", "fail_mito_ceiling")]["n_sole"],
      "if these were equal the test data would not exercise the distinction")

check("n_removed_any is the union, not the sum",
      by[("A", "fail_umi_floor")]["n_removed_any"] == 2
      and by[("B", "fail_umi_floor")]["n_removed_any"] == 3,
      f"A={by[('A', 'fail_umi_floor')]['n_removed_any']}, "
      f"B={by[('B', 'fail_umi_floor')]['n_removed_any']}; the totals per criterion sum to "
      f"{by[('A', 'fail_umi_floor')]['n_fired'] + by[('A', 'fail_mito_ceiling')]['n_fired']} "
      f"in A")

check("the cohort row sums the libraries",
      by[(BREAKDOWN_ALL, "fail_mito_ceiling")]["n_fired"] == 4
      and by[(BREAKDOWN_ALL, "fail_mito_ceiling")]["n_sole"] == 2
      and by[(BREAKDOWN_ALL, "fail_umi_floor")]["n_in"] == 10,
      f"{by[(BREAKDOWN_ALL, 'fail_mito_ceiling')]}")

check("the cohort row is labelled, not blank", BREAKDOWN_ALL == "ALL",
      "a blank sample reads as a library whose name was lost")

check("the rate is per library, not per cohort",
      by[("A", "fail_mito_ceiling")]["pct_of_library"] == 20.0
      and by[("B", "fail_mito_ceiling")]["pct_of_library"] == 60.0,
      f"A={by[('A', 'fail_mito_ceiling')]['pct_of_library']}%, "
      f"B={by[('B', 'fail_mito_ceiling')]['pct_of_library']}%")

# ---------------------------------------------------------------- 5-6: it reaches the page
block = {"source": "tables/removal_by_criterion.csv",
         "criteria": CRITERIA, "samples": ["A", "B"], "rows": rows}
defects: list = []
built = B.build_removal_breakdown({"removal_breakdown": block}, defects)

check("the block validates with no defect", built["stated"] and not defects,
      f"defects={[d.get('what') for d in defects]}")
check("the widest ratio across libraries is computed",
      built["widest_ratio"] == 3.0 and built["widest_ratio_criterion"] == "fail_mito_ceiling",
      f"{built['widest_ratio']}x on {built['widest_ratio_criterion']!r}; "
      f"the ceiling takes 20% of A and 60% of B")

html = B._removal_breakdown_table(built)
for want in ("fail_mito_ceiling", "fail_umi_floor", ">A<", ">B<", "ALL",
             "3.00&times;" if "&times;" in html else "3.00×",
             "removal_by_criterion.csv"):
    check(f"rendered: {want!r} appears", want in html,
          "" if want in html else "absent from the rendered table")
check("the sole count reaches the page", ">2<" in html,
      "library B's two sole mitochondrial removals")

# End to end, through the whole document builder, because a section can validate and render and
# still never be placed in the page - which has happened here before, twice.
payload = {"run": {"project": "t", "mode": "apply"},
           "deliverable": {"n_in": 10, "n_kept": 5, "unit": "observations"},
           "removal_breakdown": block}
with tempfile.TemporaryDirectory() as td:
    out = Path(td) / "r.html"
    B.build_report(payload, out, Path(td) / "r.json")
    page = out.read_text(encoding="utf-8")
check("the table is placed in the built document", "The same, counted" in page
      and "fail_mito_ceiling" in page,
      "build_report() must render it, not merely validate it")

# ---------------------------------------------------------------- 7: absent is a notice
defects2: list = []
absent = B.build_removal_breakdown({}, defects2)
check("an absent breakdown is a notice, not a defect",
      absent["stated"] is False and len(defects2) == 1
      and defects2[0].get("severity") == "notice",
      f"{[(d.get('severity')) for d in defects2]}")
check("and renders no table at all", B._removal_breakdown_table(absent) == "",
      "a table of zeroes would claim every criterion removed nothing")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("removal breakdown OK - counted per library, and it reaches the page")
