# Audits a recorded ambient-correction run against the auditor's checks and prints the findings.
"""Step 1c test: a real ambient-correction run must PASS; three injected faults must be caught.

A known-good alone proves nothing - a checker that returns "fine" unconditionally also passes it.
So each fault is injected into the same real cohort, one at a time.

This suite needs a cohort you hold locally and pandas to read it. Both absences are a SKIP with
a message rather than a failure, so the suite stays runnable where neither is present.
"""
import os
import sys
from pathlib import Path

# The cohort this pipeline is being checked against. It is NOT part of the pipeline: point
# COHORT_DIR at a directory you hold, whose recorded tables this suite reads.
_COHORT_ENV = os.environ.get("COHORT_DIR", "").strip()
if not _COHORT_ENV:
    print("SKIP: set COHORT_DIR to a cohort directory to audit"); raise SystemExit(0)
_COHORT = Path(_COHORT_ENV).expanduser()

COHORT = _COHORT / "results/tables"
if not (COHORT / "cellbender_before_after_summary.csv").exists():
    print(f"SKIP: no recorded ambient-correction tables under {COHORT}"); raise SystemExit(0)

try:
    import pandas as pd
except ImportError as e: # noqa: BLE001
    print(f"SKIP: pandas is needed to read the recorded tables ({e})"); raise SystemExit(0)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "01_ambient"))
from audit_ambient import audit, verdict # noqa: E402

summ = pd.read_csv(COHORT / "cellbender_before_after_summary.csv")
genes = pd.read_csv(COHORT / "cellbender_before_after_genes.csv.gz")
cfg = pd.read_csv(COHORT / "celescope_run_config.csv")
# Design factors are whatever this cohort recorded, not a fixed list. A factor is any
# column with 2-6 distinct values across the libraries; identifiers and free text are not
# factors and neither is a column that is constant, which cannot separate anything.
DESIGN_MAX_LEVELS = 6
_cand = [c for c in cfg.columns if c != "sample"
         and 2 <= cfg[c].nunique() <= DESIGN_MAX_LEVELS]
design = {f: dict(zip(cfg["sample"], cfg[f].astype(str))) for f in _cand}
if not design:
    print("SKIP: no design factor with 2-%d levels in celescope_run_config.csv"
          % DESIGN_MAX_LEVELS); raise SystemExit(0)
print(f"design factors discovered: {list(design)}")

fails = []
print("Step 1c - CellBender removal auditor\n" + "=" * 74)
print("A. the cohort as delivered (known good)")
f = audit(summ, genes, design)
for x in f:
    print(x)
v = verdict(f)
print(f"\n verdict: {v}")
if v == "REFUSE":
    fails.append("A: the accepted run must not be refused")

print("\n" + "-" * 74)
print("B. injected: one sample's removal doubled (degenerate fit)")
s2 = summ.copy()
# Take the target FROM the data. A hardcoded name selects zero rows on any other cohort,
# and a zero-row injection is indistinguishable from a detector that failed to fire.
_tgt = summ["sample"].iloc[0]
_sel = s2["sample"] == _tgt
assert _sel.sum() == 1, f"injection must hit exactly one library, hit {_sel.sum()}"
s2.loc[_sel, "fraction_removed_overall"] *= 2.6
print(f" injected into {_tgt}: "
      f"{summ.loc[_sel, 'fraction_removed_overall'].iloc[0]:.3f} -> "
      f"{s2.loc[_sel, 'fraction_removed_overall'].iloc[0]:.3f}")
f2 = audit(s2, None, design)
hit = [x for x in f2 if "cohort outlier" in x.check and x.severity != "ok"]
print(hit[0] if hit else " NOT CAUGHT")
if not hit:
    fails.append("B: a doubled removal must be flagged as a cohort outlier")

print("\n" + "-" * 74)
print("C. injected: removal 3.4x higher in one arm of a design factor (the design-differential check)")
s3 = summ.copy()
# Inject into ONE level of the first two-level factor, whichever it is.
_f = next((f for f, m in design.items() if len(set(m.values())) == 2), None)
if _f is None:
    print("SKIP: no two-level factor to inject a one-sided difference into")
    raise SystemExit(0)
_lvl = sorted(set(design[_f].values()))[-1]
arm = [k for k, v in design[_f].items() if v == _lvl]
_sel3 = s3["sample"].isin(arm)
assert _sel3.sum() == len(arm) and 0 < _sel3.sum() < len(s3), (
    f"injection must hit one arm, hit {_sel3.sum()} of {len(s3)}")
s3.loc[_sel3, "fraction_removed_overall"] *= 3.4
f3 = audit(s3, None, design)
# The check is named after the factor this cohort actually recorded, not after a fixed
# vocabulary - so ask for the factor that was injected into.
hit = [x for x in f3 if x.check.endswith(_f)]
print(hit[0] if hit else " NOT CAUGHT")
if not (hit and hit[0].severity == "REFUSE"):
    fails.append("C: a 3.4x one-arm differential must REFUSE")

print("\n" + "-" * 74)
print("D. injected: three genes gutted")
g4 = genes.copy()
# The three most widely detected genes in THIS cohort. Naming symbols would tie the test
# to one tissue and would inject nothing wherever those symbols are absent - which reads
# as a detector that failed rather than a test that was never run.
_det = "raw_detection_frac" if "raw_detection_frac" in g4.columns else "fraction_removed"
present = (g4.groupby("symbol")[_det].mean().sort_values(ascending=False)
           .head(3).index.tolist())
assert len(present) == 3, f"need 3 genes to gut, found {len(present)}"
print(f" gutting the 3 most-detected genes: {present}")
g4.loc[g4.symbol.isin(present), "fraction_removed"] = 0.97
g4.loc[g4.symbol.isin(present), "denoised_detection_frac"] = 0.001
f4 = audit(summ, g4, design)
hit = [x for x in f4 if x.check == "genes gutted"]
print(hit[0] if hit else " NOT CAUGHT")
if not (hit and hit[0].severity == "REVIEW"):
    fails.append("D: gutted genes must be flagged AND listed")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: a real run passes and all three injected faults are caught")
print("note WHICH check found B: the total removal is the number the denoiser reports and the")
print("LEAST informative of the five - the calibration cohort's degenerate library sat at 8.3%,")
print("LOW rather than high, comfortably inside a plausible-looking range")
