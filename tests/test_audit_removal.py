"""The step-7 auditor: does it catch the failures that leave no gap?

An auditor is only worth its runtime if it FAILS on broken input. Every case below is a way a
removal can be wrong while every object still has a plausible size and every count still
reconciles - which is the whole reason steps 2, 3 and 4 need auditing at all: they restrict and
flag without dropping, so nothing looks missing when they go wrong.

Run: python tests/test_audit_removal.py
"""
from __future__ import annotations

import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "audit_removal", ROOT / "modules/07_apply/audit_removal.py")
ar = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = ar
spec.loader.exec_module(ar)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"   -- {detail}" if detail and not cond
                                                       else ""))


def sev(findings, needle):
    for f in findings:
        if needle in f.check:
            return f.severity
    return "MISSING"


# A clean 10-droplet cohort: 2 uncalled, 3 below the light floor, 1 doublet.
#   idx  0  1  2  3  4  5  6  7  8  9
CELL = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
UMI = [50, 80, 100, 150, 180, 400, 500, 600, 700, 800]
SCORED = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]          # cell & umi >= 200
F_UMI = [0, 0, 1, 1, 1, 0, 0, 0, 0, 0]           # cell & umi < 350
F_MITO = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
F_DBL = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
REMOVED = [1, 1, 1, 1, 1, 0, 1, 0, 0, 1]         # not-cell OR any criterion
KEEP = [0, 0, 0, 0, 0, 1, 0, 1, 1, 0]

BASE = {"cellbender_cell": CELL, "keep": KEEP, "removed": REMOVED, "scored": SCORED,
        "total_counts": UMI, "f_umi": F_UMI, "f_mito": F_MITO, "f_dbl": F_DBL}
CRIT = ("f_umi", "f_mito", "f_dbl")
KW = dict(criteria=CRIT, scored_col="scored", light_floor=200, quality_floor=350,
          doublet_criterion="f_dbl")


def run(rows, **over):
    kw = dict(KW)
    kw.update(over)
    return ar.audit(rows, **kw)


print("\nA. a correct removal passes")
f = run(BASE, predoublet_keep=[k or d for k, d in zip(KEEP, F_DBL)])
for c in ("complementary", "decomposes", "step 2", "no unexamined", "light floor is strictly",
          "superset", "doublets ALONE", "outside the examined"):
    check(f"{c}: ok", sev(f, c) == "ok", f"got {sev(f, c)}")
check("verdict is not FAIL", ar.verdict(f) != "FAIL")

print("\nB. keep and removed disagreeing is caught")
bad = dict(BASE, keep=[1] + KEEP[1:])            # row 0 both kept and removed
check("both kept and removed -> FAIL", sev(run(bad), "complementary") == "FAIL")

print("\nC. a removal with no recorded criterion is caught")
# Row 7 is a called, scored, criterion-free cell - flipped to removed. This is what a removal
# performed outside step 7 looks like from the table.
bad = dict(BASE, removed=[*REMOVED[:7], 1, *REMOVED[8:]], keep=[*KEEP[:7], 0, *KEEP[8:]])
check("unexplained removal -> FAIL", sev(run(bad), "decomposes") == "FAIL")

print("\nD. step 2 contamination is caught")
bad = dict(BASE, keep=[1, *KEEP[1:]], removed=[0, *REMOVED[1:]])   # uncalled droplet kept
check("uncalled droplet in deliverable -> FAIL", sev(run(bad), "step 2") == "FAIL")

print("\nE. step 3 contamination is caught - the failure with no visible symptom")
# Quality floor below the light floor: an unexamined nucleus survives carrying 'not a doublet'
# because nothing looked at it. Object size is plausible; nothing is missing.
bad = dict(BASE, keep=[0, 0, 0, 0, 1, 1, 0, 1, 1, 0],
           removed=[1, 1, 1, 1, 0, 0, 1, 0, 0, 1], f_umi=[0, 0, 1, 1, 0, 0, 0, 0, 0, 0])
f = run(bad, quality_floor=150)
check("unexamined nucleus kept -> FAIL", sev(f, "no unexamined") == "FAIL")
check("floor ordering violated -> FAIL", sev(f, "light floor is strictly") == "FAIL")
check("both fire together (one is the cause of the other)",
      sev(f, "no unexamined") == "FAIL" and sev(f, "light floor is strictly") == "FAIL")

print("\nF. a doublet call on a never-examined nucleus is caught")
bad = dict(BASE, f_dbl=[0, 0, 1, 0, 0, 0, 0, 0, 0, 1])
check("call outside the examined set -> FAIL", sev(run(bad), "outside the examined") == "FAIL")

print("\nG. the step-6 exception is bounded, not a free pass")
f = run(BASE, predoublet_keep=[0] * 10)
check("not a superset -> FAIL", sev(f, "superset") == "FAIL")
# A superset that retains something OTHER than doublets is a second, unaudited filter.
f = run(BASE, predoublet_keep=[k or d for k, d in zip(KEEP, F_DBL)][:6] + [1, 1, 1, 1])
check("superset differing by non-doublets -> FAIL", sev(f, "doublets ALONE") == "FAIL")

print("\nH. an absent column fails the audit rather than skipping a check")
try:
    ar.audit({"cellbender_cell": CELL, "keep": KEEP}, **KW)
    check("missing column -> AuditFailure", False, "no exception")
except ar.AuditFailure as e:
    check("missing column -> AuditFailure", "removed" in str(e) or "column" in str(e))

print("\nI. a criterion with no unique contribution is flagged, not hidden")
# f_mito duplicated into a criterion that never fires alone.
rows = dict(BASE, f_dup=F_MITO)
f = run(rows, criteria=("f_umi", "f_mito", "f_dup", "f_dbl"))
check("duplicate criterion -> REVIEW", sev(f, "removes something no other does") == "REVIEW")

print("\nJ. enforce() raises on FAIL and is silent on REVIEW")
try:
    ar.enforce(run(dict(BASE, keep=[1] + KEEP[1:])))
    check("enforce raises on FAIL", False, "did not raise")
except ar.AuditFailure:
    check("enforce raises on FAIL", True)
ar.enforce(run(BASE, predoublet_keep=[k or d for k, d in zip(KEEP, F_DBL)]))
check("enforce is silent on a clean audit", True)

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
if FAIL:
    for x in FAIL:
        print(f"  FAILED: {x}")
    sys.exit(1)
