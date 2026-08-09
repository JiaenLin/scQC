# Recomputes a held cohort's recorded QC counts through the pipeline modules and compares them.
"""Tier 1 acceptance: drive steps 2-7 on a cohort's real per-cell table and reproduce its counts.

WHAT THIS TESTS AND WHAT IT DOES NOT

Tier 1 is the deterministic half. Given the stored denoised output, every step from the cell call
to the final removal is reproducible, so every checkpoint must match EXACTLY. Nothing in tier 1
has a legitimate reason to move.

It does NOT test the aligner or the denoiser - those are tier 2, whose tolerance is still recorded
as UNMEASURED because the run-to-run spread of variational inference has not been measured.

It does not validate the cohort either. It proves the pipeline does what that cohort's recorded
result says was done.

WHY IT USES A REAL TABLE AND NOT FIXTURES

The unit suites assert against transcribed numbers. That catches logic errors and misses
transcription errors - a boundary value copied from a rounded display flips a strict inequality
and the module looks wrong when the fixture is. This runner reads the actual per-cell table, so
the fixtures cannot be wrong about the data.

WHAT IS COHORT-SPECIFIC, AND HOW TO POINT IT AT YOUR OWN

Two things vary between cohorts and neither is baked in:

  expected.local.tsv  the checkpoint VALUES. Gitignored, because they describe data the reader
                      does not have. Copy expected.tsv.template and fill it in.
  schema.local.tsv    the LAYOUT - table filenames, column names, the design factor. Optional:
                      without it the defaults below apply, which are the layout the calibration
                      cohort recorded. Every key is documented in README.md.

A checkpoint with no expectation is a FAILURE, not a pass. A check that cannot fail reports
success under every regression it exists to catch.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# The cohort this pipeline is being checked against. It is NOT part of the pipeline: point
# COHORT_DIR at a directory you hold, whose recorded tables the expectations describe.
_COHORT_ENV = os.environ.get("COHORT_DIR", "").strip()
if not _COHORT_ENV:
    print("SKIP: set COHORT_DIR to the cohort directory to check against")
    raise SystemExit(0)
_COHORT = Path(_COHORT_ENV).expanduser()
if not _COHORT.is_dir():
    print(f"SKIP: COHORT_DIR is not a directory: {_COHORT}")
    raise SystemExit(0)

# Expectations describe ONE cohort and are not tracked: `.local.tsv` is gitignored. A missing
# file here is NOT a skip - COHORT_DIR was set, so the check was asked for, and passing it with
# nothing to compare against would be worse than failing it.
_EXPECT_FILE = HERE / "expected.local.tsv"
if not _EXPECT_FILE.exists():
    raise SystemExit(f"no expectations at {_EXPECT_FILE} - copy expected.tsv.template, fill it "
                     f"from this cohort's own recorded tables, and run again")

try:
    import pandas as pd
except ImportError as e: # noqa: BLE001
    print(f"SKIP: pandas is needed to read the recorded tables ({e})")
    raise SystemExit(0)

# The repository root, resolved from this file. A hardcoded install path makes the tests
# runnable only on the machine they were written on.
P = HERE.parents[1]
COHORT = _COHORT
for m in ("02_cells", "03_light_floor", "04_doublets", "05_quality", "06_cluster_check",
          "07_apply"):
    sys.path.insert(0, str(P / "modules" / m))

from cellcall_gate import gate, verdict as cc_verdict # noqa: E402
from light_floor import assess as floor_assess # noqa: E402
from doublet_health import health, verdict as dh_verdict # noqa: E402
from quality import Valley, derive # noqa: E402
from apply import preflight, propose_cluster_removal # noqa: E402

# ---------------------------------------------------------------------------- the local schema
# Defaults are the layout of the calibration cohort's recorded tables. They are DEFAULTS, not
# requirements: schema.local.tsv overrides any of them, key and value separated by a tab. See
# README.md for what each key means.
SCHEMA = {
    "tables_dir": "results/tables",
    "objects_dir": "results/objects/cellbender_h5ad",
    "denoised_glob": "*_lr5e-5.h5ad",
    "denoised_sample_sep": "_cellbender",
    "per_cell_table": "qc_filter_per_cell.csv.gz",
    "doublet_sweep_table": "doublet_dbrsd_sweep_per_cell.csv.gz",
    "cluster_profile_table": "cluster_check_profile.csv.gz",
    "valley_table": "quality_valleys.csv",
    "col_sample": "sample",
    "col_aligner_cell": "celescope_cell",
    "col_ambient_cell": "cellbender_cell",
    "col_doublet_score": "doublet_score",
    "col_total_counts": "total_counts",
    "col_keep": "keep",
    "col_doublet_call": "scdbl_sd0.06_call",
    "fail_columns": "fail_umi_lt_350,fail_genes_lt_250,fail_mito_gt_40,fail_scdblfinder_sd006",
    "doublet_fail_column": "fail_scdblfinder_sd006",
    "design_factor": "condition",
    "design_levels": "treat=treat,ctrl=ctrl",
    "light_floor": "200",
    "quality_floor": "350",
    "cluster_resolution": "1.0",
    "cluster_algorithm": "leiden",
}
_SCHEMA_FILE = HERE / "schema.local.tsv"
if _SCHEMA_FILE.exists():
    for raw in _SCHEMA_FILE.read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        k, _, v = raw.partition("\t")
        k, v = k.strip(), v.strip()
        if k not in SCHEMA:
            raise SystemExit(f"{_SCHEMA_FILE.name}: unknown key {k!r}. Known keys: "
                             f"{', '.join(sorted(SCHEMA))}")
        SCHEMA[k] = v
    print(f"schema overrides read from {_SCHEMA_FILE.name}")

FAIL_COLS = [c.strip() for c in SCHEMA["fail_columns"].split(",") if c.strip()]
DOUBLET_FAIL = SCHEMA["doublet_fail_column"].strip()
if DOUBLET_FAIL not in FAIL_COLS:
    raise SystemExit(f"doublet_fail_column {DOUBLET_FAIL!r} is not one of fail_columns "
                     f"({', '.join(FAIL_COLS)}) - step 5's checkpoint counts what survives the "
                     f"QUALITY criteria, so it has to know which column is the doublet one")
# Step 5 is the count and mitochondrial filters. The doublet criterion is step 4's and is
# deliberately excluded here: folding it in would report one number for two steps.
QUALITY_FAILS = [c for c in FAIL_COLS if c != DOUBLET_FAIL]
SAMPLE = SCHEMA["col_sample"]
FACTOR = SCHEMA["design_factor"]

def design_level(sample: str) -> str:
    """The design level of one library, by the first `level=substring` rule that matches.

    A level that matches nothing is reported as `unassigned` rather than silently folded into
    the last one: a library that quietly joins the wrong arm is exactly the fault the
    design-differential checks exist to catch.
    """
    for pair in SCHEMA["design_levels"].split(","):
        level, _, needle = pair.partition("=")
        if needle.strip() and needle.strip() in sample:
            return level.strip()
    return "unassigned"

def design_of(samples) -> dict:
    return {FACTOR: {s: design_level(str(s)) for s in samples}}

# ------------------------------------------------------------------------------ expectations
def read_expectations(path: Path) -> dict:
    out = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        parts = raw.split("\t")
        if len(parts) < 3 or parts[0].strip() != "1":
            continue
        out[parts[1].strip()] = parts[2].strip()
    return out

EXPECT = read_expectations(_EXPECT_FILE)

T = COHORT / SCHEMA["tables_dir"]
fails, checks = [], []

def check(name, got):
    """Compare one checkpoint. A missing or blank expectation FAILS.

    Falling through to "no expectation" makes a checkpoint that can never fail, and a suite of
    those reports success against every regression it was written to catch.
    """
    exp = EXPECT.pop(name, None)
    if exp is None or exp == "":
        checks.append((name, got, "-", "NO EXP"))
        fails.append(f"{name}: no expectation in {_EXPECT_FILE.name}. A checkpoint with nothing "
                     f"to compare against can never fail, so it is counted as a failure")
        return
    try:
        want = int(float(exp))
    except ValueError:
        checks.append((name, got, exp, "BAD"))
        fails.append(f"{name}: expectation {exp!r} is not a number. Tier 1 checkpoints are "
                     f"counts and carry no tolerance, so there is nothing else it could be")
        return
    ok = int(got) == want
    checks.append((name, got, exp, "PASS" if ok else "FAIL"))
    if not ok:
        fails.append(f"{name}: got {got:,}, expected {want:,}")

print("Tier 1 acceptance - steps 2-7 against a cohort's recorded tables")
print(f"cohort: {COHORT}")
print("=" * 78)

cells = pd.read_csv(T / SCHEMA["per_cell_table"])
cb = cells[cells[SCHEMA["col_ambient_cell"]].astype(bool)]
check("droplets_analysed", len(cells))
check("ambient_cells", len(cb))

# ---- step 2: the cell-call gate, against the aligner's own call
print("\nstep 2 - cell-call gate")
try:
    import anndata as ad
    calls = {}
    for f in sorted((COHORT / SCHEMA["objects_dir"]).glob(SCHEMA["denoised_glob"])):
        s = f.name.split(SCHEMA["denoised_sample_sep"])[0]
        a = ad.read_h5ad(f, backed="r")
        cs = a.obs[SCHEMA["col_aligner_cell"]].astype(bool)
        cbx = a.obs[SCHEMA["col_ambient_cell"]].astype(bool)
        calls[s] = dict(aligner=int(cs.sum()), cellbender=int(cbx.sum()),
                        lost=int((cs & ~cbx).sum()))
        a.file.close()
    if not calls:
        raise FileNotFoundError(f"no {SCHEMA['denoised_glob']} under {SCHEMA['objects_dir']}")
    g = gate(calls, design_of(calls))
    print(f" verdict: {cc_verdict(g)}")
    for x in g:
        if x.severity != "ok":
            print(" " + str(x).replace("\n", "\n ")[:230])
except Exception as e: # noqa: BLE001
    print(f" SKIP ({type(e).__name__}: {e})")

# ---- step 3: the light floor and its coverage
print("\nstep 3 - light floor")
scored = cb[SCHEMA["col_doublet_score"]].notna()
keep = cb[SCHEMA["col_keep"]].astype(bool)
cov = floor_assess(n_total=len(cb), n_scored=int(scored.sum()),
                   floor=int(SCHEMA["light_floor"]),
                   quality_floor=int(SCHEMA["quality_floor"]),
                   max_unscored_umi=float(cb.loc[~scored, SCHEMA["col_total_counts"]].max()),
                   min_scored_umi=float(cb.loc[scored, SCHEMA["col_total_counts"]].min()),
                   n_kept_unscored=int((~scored & keep).sum()))
print(" " + str(cov).replace("\n", "\n "))
check("scored_for_doublets", int(scored.sum()))

# ---- step 4: the doublet calls and their health
print("\nstep 4 - doublet calls")
sw = pd.read_csv(T / SCHEMA["doublet_sweep_table"])
called = sw[SCHEMA["col_doublet_call"]].astype(bool)
check("doublets_called", int(called.sum()))
rate = sw.assign(_c=called).groupby(SAMPLE, observed=True)["_c"].mean().to_dict()
h = health(rate, design_of(rate), n_kept_unscored=int((~scored & keep).sum()),
           detector_name="scDblFinder", reproducible=True)
print(f" verdict: {dh_verdict(h)}")
for x in h:
    if x.severity != "ok":
        print(" " + str(x).replace("\n", "\n ")[:200])

# ---- step 5: the count floors
# The valleys are DERIVED per cohort, so tier 1 cannot carry a list of them: one cohort's
# valleys are not another's, which is the reason this pipeline derives them at all. Zipping a
# stored list against whatever samples a cohort happens to have is worse than not checking,
# because zip() truncates in silence and the mismatch never surfaces.
#
# So: re-derive if this cohort recorded its own valley table, and otherwise say plainly that the
# DERIVATION is not covered. The checkpoints below cover the other half - that the thresholds
# recorded were the thresholds APPLIED.
print("\nstep 5 - quality thresholds")
_vf = T / SCHEMA["valley_table"]
if _vf.exists():
    v = pd.read_csv(_vf)

    def _bimodal(x) -> bool:
        """A blank cell is NOT a claim of bimodality.

        pandas reads an empty cell as NaN and `bool(nan)` is True, so a library whose bimodality
        was never assessed asserted that it WAS bimodal - and derive() then proposed a threshold
        from a distribution nobody had looked at. An absent answer is a no: bimodality has to be
        stated to count.
        """
        if x is None or (isinstance(x, float) and x != x):
            return False
        return str(x).strip().lower() not in ("", "false", "0", "no", "nan")

    umi_valleys = [Valley(str(r[SAMPLE]), "umi", float(r["valley"]), _bimodal(r["bimodal"]))
                   for _, r in v.iterrows() if str(r["metric"]) == "umi"]
    if umi_valleys:
        # derive() REFUSES rather than returns when the valleys are unusable. That refusal is a
        # correct verdict about the valleys, not a failure of the acceptance run, so it is
        # reported and the remaining checkpoints still execute. Uncaught, it aborted the whole
        # run and reported nothing about the checkpoints that follow it.
        try:
            pu = derive(umi_valleys, "umi", light_floor=int(SCHEMA["light_floor"]))
            print(" " + str(pu).replace("\n", "\n "))
        except Exception as e:                                        # noqa: BLE001
            print(f" REFUSED: {type(e).__name__}: {e}")
            print("  A verdict about this cohort's valleys, not about the pipeline; the")
            print("  checkpoints below are unaffected and still run.")
    else:
        print(f" NOT COVERED: {_vf.name} records no umi valleys")
else:
    print(f" NOT COVERED: no valley table at {_vf.name}, so threshold DERIVATION is outside")
    print("  tier 1 here. The checkpoints below test that the recorded thresholds were")
    print("  applied, not that they were derived.")

survives = ~cb[QUALITY_FAILS].any(axis=1)
check("step5_after_quality", int(survives.sum()))

# ---- unique contributions: how much each criterion removes that nothing else would
for a in FAIL_COLS:
    others = cb[[x for x in FAIL_COLS if x != a]].any(axis=1)
    check("only_" + a.removeprefix("fail_"), int((cb[a].astype(bool) & ~others).sum()))

# ---- step 6 + 7
print("\nsteps 6-7 - cluster pre-flight and the deliverable")
kept = int(keep.sum())
check("deliverable", kept)
prof_f = T / SCHEMA["cluster_profile_table"]
if prof_f.exists():
    pr = pd.read_csv(prof_f)
    pr = pr[(pr.resolution == float(SCHEMA["cluster_resolution"]))
            & (pr.algorithm == SCHEMA["cluster_algorithm"])]
    rows = pr.to_dict("records")
    for x in preflight(rows, kept_total=kept):
        print(" " + str(x).replace("\n", "\n ")[:240])
    prop = propose_cluster_removal(rows, design=design_of(pr[SAMPLE].unique())[FACTOR])
    print(f" proposed cluster removal: {prop['n_removed']} nuclei - {prop['status'][:60]}...")
else:
    print(f" SKIP: no cluster profile at {prof_f.name}")

print("\n" + "=" * 78)
# check() consumes each name it uses, so whatever is left was never checked. A row nobody reads
# is how a file comes to describe coverage the runner does not have.
for leftover in sorted(EXPECT):
    fails.append(f"{leftover}: recorded in {_EXPECT_FILE.name} but no checkpoint of that name "
                 f"is checked. Align it with the names in expected.tsv.template, or remove it")
w = max(len(c[0]) for c in checks)
for n, got, exp, st in checks:
    print(f" {st:6s} {n:<{w}} got {int(got):>9,} expected {exp:>9}")
print("=" * 78)
compared = sum(1 for c in checks if c[3] in ("PASS", "FAIL"))
if fails:
    print(f"TIER 1 FAILED - {len(fails)} problem(s) over {len(checks)} checkpoint(s), of which "
          f"{compared} were compared against {_EXPECT_FILE.name}:")
    for f in fails:
        print(" -", f)
    raise SystemExit(1)
print(f"TIER 1 PASSED - {compared} of {len(checks)} checkpoints compared against "
      f"{_EXPECT_FILE.name} and reproduced exactly")
