"""SUPPLIED ambient correction: a third state, and not a kind of skip.

Re-entering the pipeline at step 2 with an already-denoised matrix is a normal thing to want.
Modelling it as a skip gets it wrong in both directions: for snRNA a skip is refused outright, so
a valid re-entry becomes impossible; and if it were allowed, the report would claim this run
performed a correction it never performed and would attach this run's parameters to someone
else's output.

The tests below are mostly refusals, because the value of the state is what it declines to let
you do anonymously.

Run: python tests/test_ambient_supplied.py
"""
from __future__ import annotations

import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ambient", ROOT / "modules/01_ambient/ambient.py")
amb = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = amb
spec.loader.exec_module(amb)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   -- {detail}" if detail and not cond
                                                       else ""))


PROV = {"tool": "CellBender", "version": "0.3.2",
        "params": "--fpr 0 --learning-rate 5e-5",
        "produced_by": "SAMBO stage 1, HPC, 2026-08-08"}

print("\nA. the three states are distinct")
check("RUN", amb.plan_ambient("s", "snrna").state == "RUN")
check("SUPPLIED", amb.plan_ambient("s", "snrna", supplied=PROV).state == "SUPPLIED")
check("SKIP (scrna only)",
      amb.plan_ambient("s", "scrna", skip=True, skip_reason="pre-cleaned").state == "SKIP")

print("\nB. SUPPLIED does not weaken the snRNA rule")
# The rule is that nuclei are never analysed UNCORRECTED. Supplied data IS corrected, so it does
# not breach the rule - but a bare skip must still be refused, or SUPPLIED becomes a loophole.
try:
    amb.plan_ambient("s", "snrna", skip=True, skip_reason="I already did it elsewhere")
    check("a bare skip is still refused for snRNA", False, "it was allowed")
except amb.AmbientRefusal as e:
    check("a bare skip is still refused for snRNA", "cannot be skipped" in str(e))
p = amb.plan_ambient("s", "snrna", supplied=PROV)
check("supplied snRNA still records the correction as mandatory", p.mandatory)
check("supplied is not run", p.run is False)

print("\nC. provenance is required, key by key")
for k in amb.SUPPLIED_REQUIRED:
    bad = {x: v for x, v in PROV.items() if x != k}
    try:
        amb.plan_ambient("s", "snrna", supplied=bad)
        check(f"missing {k!r} refuses", False, "accepted")
    except amb.AmbientRefusal as e:
        check(f"missing {k!r} refuses", k in str(e))
for k in amb.SUPPLIED_REQUIRED:
    blank = dict(PROV, **{k: "   "})
    try:
        amb.plan_ambient("s", "snrna", supplied=blank)
        check(f"whitespace-only {k!r} refuses", False, "accepted")
    except amb.AmbientRefusal:
        check(f"whitespace-only {k!r} refuses", True)

print("\nD. contradictions and wrong shapes refuse")
try:
    amb.plan_ambient("s", "snrna", skip=True, skip_reason="x", supplied=PROV)
    check("skip + supplied refuses", False, "accepted")
except amb.AmbientRefusal as e:
    check("skip + supplied refuses", "contradictory" in str(e))
for bad in ("yes", 1, ["CellBender"], object()):
    try:
        amb.plan_ambient("s", "snrna", supplied=bad)
        check(f"non-dict {type(bad).__name__} refuses", False, "accepted")
    except amb.AmbientRefusal:
        check(f"non-dict {type(bad).__name__} refuses", True)

print("\nE. the recorded reason carries what a reader needs")
r = amb.plan_ambient("s", "snrna", supplied=PROV).reason
for needle in ("CellBender", "0.3.2", "5e-5", "SAMBO"):
    check(f"reason names {needle!r}", needle in r)
check("reason says this run did NOT perform it", "NOT performed by this run" in r)
check("reason says the fraction removed is not measurable here", "NOT MEASURABLE" in r)
check("the provenance survives on the plan",
      amb.plan_ambient("s", "snrna", supplied=PROV).supplied == PROV)

print("\nF. the assay cross-check still fires under SUPPLIED")
# A supplied object does not exempt the declaration from being checked against the data.
try:
    amb.plan_ambient("s", "scrna", supplied=PROV, intronic_fraction=0.85)
    check("mis-declared assay still refuses when supplied", False, "accepted")
except amb.AmbientRefusal:
    check("mis-declared assay still refuses when supplied", True)

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
if FAIL:
    for x in FAIL:
        print(f"  FAILED: {x}")
    sys.exit(1)
