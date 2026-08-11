"""Section 2 and the freshness check: two validators that had no producer.

`report/build.py` has validated a parameter table and computed a freshness verdict since the
report existed, and nothing in the engine ever supplied `payload["parameters"]` or
`provenance["newest_input"]`. So the section the design document calls "the point of this report"
never rendered, and every report said `NOT CHECKED` about its own staleness.

THE PROBE THAT MATTERS MOST is not that the table appears. It is that the CLASS in it is the class
step 7 actually applied. `steps._apply_thresholds` treats a decisions entry as ADJUDICATED only
when it carries a value, an `approved_by` AND the operator's own words; anything less is a number
somebody typed and the DERIVED value is used instead. If the report decided that differently it
would print a class the run did not apply - the document describing a different filter from the
one that ran, which is the failure this whole section exists to prevent.

Run: python tests/test_parameters_and_freshness.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.pipeline import Pipeline  # noqa: E402
from report.build import build_parameter_table, freshness, newest_input_time  # noqa: E402

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}"
          + (f"  --  {detail}" if detail and not cond else ""))


def fake(decisions=None, tools=None, samples=None, project=None):
    """The attributes `_parameters_block` and `_newest_input_block` read, and nothing else.

    A real Pipeline opens a project directory, builds a run key and takes a lock. The two methods
    under test read four attributes; constructing the rest would test the constructor.
    """
    return SimpleNamespace(
        FIXED_PARAMETERS=Pipeline.FIXED_PARAMETERS,
        DECLARED_PARAMETERS=Pipeline.DECLARED_PARAMETERS,
        decisions=decisions or {},
        tools=tools or {},
        samples=samples or [],
        project=Path(project) if project else Path("."),
        results_by_key={
            "05_quality": SimpleNamespace(metrics={"umi_proposed": 350, "genes_proposed": 260}),
            "07_apply": SimpleNamespace(metrics={"ceilings": "10.00-20.51% across 10 libraries"}),
        })


def rows_by_name(block):
    return {r["name"]: r for r in (block or [])}


print("\nthe parameter table is produced at all")
block = Pipeline._parameters_block(fake(tools={"dbr": 0.06, "light_floor": 200, "seed": 0}))
check("a block is produced", bool(block), str(block))
by = rows_by_name(block)
check("FIXED rows are present", any(r["class"] == "FIXED" for r in block))
check("DECLARED rows come from the tools dict", by.get("light floor (UMI)", {}).get("class")
      == "DECLARED", str(by.get("light floor (UMI)")))
check("a tools key that was not given produces no row", "clustering resolution" not in by,
      str(sorted(by)))
check("the derived floors are DERIVED", by.get("UMI floor", {}).get("class") == "DERIVED")
check("...and carry a basis, which is what makes them reviewable",
      bool(by.get("UMI floor", {}).get("basis")))

print("\nthe class is the one step 7 applied - the probe that matters")
FULL = {"quality": {"umi_floor": {"value": 500, "approved_by": "PI",
                                  "verbatim": "apply 500 as the floor"}}}
HALF = {"quality": {"umi_floor": {"value": 500, "approved_by": "PI"}}}       # no verbatim
BARE = {"quality": {"umi_floor": {"value": 500}}}                            # no attestation

by_full = rows_by_name(Pipeline._parameters_block(fake(decisions=FULL)))
check("a fully attested decision is ADJUDICATED",
      by_full["UMI floor"]["class"] == "ADJUDICATED", str(by_full["UMI floor"]))
check("...and carries the operator's own words",
      by_full["UMI floor"].get("verbatim") == "apply 500 as the floor")
check("...and its value is the declared one, not the derived one",
      by_full["UMI floor"]["value"] == 500)

for label, dec in (("no verbatim", HALF), ("no attestation at all", BARE)):
    r = rows_by_name(Pipeline._parameters_block(fake(decisions=dec)))["UMI floor"]
    check(f"a decision with {label} is DERIVED, as step 7 treats it",
          r["class"] == "DERIVED", str(r))
    check(f"...and shows the DERIVED value ({label})", r["value"] == 350, str(r["value"]))

print("\nthe report's own validator accepts what the pipeline emits")
defects = []
out = build_parameter_table({"parameters": Pipeline._parameters_block(
    fake(tools={"dbr": 0.06, "seed": 0}, decisions=FULL))}, defects)
check("every row validates", out["stated"] and not defects,
      "; ".join(d["what"][:60] for d in defects))
check("...and every class is one the report recognises",
      all(r.get("class") in ("FIXED", "DECLARED", "DERIVED", "ADJUDICATED")
          for r in out["rows"] if not r.get("defect")))

print("\nfreshness is computed, and NOT CHECKED where it cannot be")
with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    (d / "samplesheet.csv").write_text("sample\n", encoding="utf-8")
    (d / "matrix").write_text("x", encoding="utf-8")
    got = Pipeline._newest_input_block(
        fake(project=d, samples=[{"matrix": str(d / "matrix")}]))
    check("a newest-input time is supplied", bool(got.get("newest_input")), str(got))
    check("the file it came from is named", bool(got.get("newest_input_path")), str(got))

    # An input the run declared and that is not there must be LISTED, never skipped: a time taken
    # over the inputs that happened to exist describes a different run.
    got2 = Pipeline._newest_input_block(
        fake(project=d, samples=[{"matrix": str(d / "gone")}]))
    check("an absent declared input is listed, not silently dropped",
          got2.get("inputs_absent"), str(got2))

check("a run with nothing to compare reports NOT CHECKED rather than fresh",
      freshness("2026-08-11T23:00:00+0800", None)["stale"] is None)
check("an artifact older than its input is STALE",
      freshness("2026-08-11T21:00:00+0800", "2026-08-11T22:00:00+0800")["stale"] is True)
check("an artifact newer than its input is not stale",
      freshness("2026-08-11T23:00:00+0800", "2026-08-11T22:00:00+0800")["stale"] is False)
# Comparing a stamp with an offset against one without would assume a zone neither states.
check("mixed timezone spellings are NOT CHECKED, not guessed",
      freshness("2026-08-11T23:00:00+0800", "2026-08-11T22:00:00")["stale"] is None)
check("newest_input_time stamps with an offset, so both sides are spelled the same way",
      (newest_input_time([__file__])["newest_input"] or "")[-5] in "+-",
      newest_input_time([__file__])["newest_input"])

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
sys.exit(1 if FAIL else 0)
