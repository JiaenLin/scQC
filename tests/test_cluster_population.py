"""Step 6 must cluster the cells that reach the deliverable, not the droplet matrix.

The defect this guards: `_op_cluster` read the denoised FULL DROPLET object, while
`cluster_flags.py` specifies "the step-5 object: quality-filtered". On one cohort that left 1,398
of 1,531 clusters holding no cell that reached the deliverable, under-resolved the real cells, and
- because every criterion is a cluster median - made the libraries with the MOST droplet noise
raise the FEWEST flags.

Each probe below must be able to FAIL. A fixture that cannot reach the failing comparison reports
"tolerated" whether or not the bug is present.

Run: python tests/test_cluster_population.py
"""
from __future__ import annotations

import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    import numpy as np
    import anndata as ad
    import pandas as pd
except ModuleNotFoundError as e:  # noqa: BLE001
    print(f"SKIP: needs {e.name}")
    raise SystemExit(0)

spec = importlib.util.spec_from_file_location("scanpy_ops", ROOT / "adapters/scanpy_ops.py")
so = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = so
spec.loader.exec_module(so)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  --  {detail}" if detail and not cond
                                                       else ""))


def cohort(n_real=400, n_empty=1600):
    """Real cells and empty droplets, as step 6 actually meets them."""
    rng = np.random.default_rng(0)
    n = n_real + n_empty
    X = np.zeros((n, 30), dtype=np.float32)
    X[:n_real] = rng.poisson(6.0, size=(n_real, 30))
    X[n_real:] = rng.poisson(0.02, size=(n_empty, 30))
    a = ad.AnnData(X)
    a.obs_names = [f"bc{i}" for i in range(n)]
    a.var_names = [f"g{i}" for i in range(30)]
    a.var["mt"] = False
    a.obs["cellbender_cell"] = np.r_[np.ones(n_real), np.zeros(n_empty)].astype(bool)
    a.obs["total_counts"] = np.asarray(X.sum(1)).ravel()
    a.obs["n_genes_by_counts"] = np.asarray((X > 0).sum(1)).ravel()
    a.obs["pct_counts_mt"] = np.r_[rng.uniform(0, 5, n_real), rng.uniform(0, 5, n_empty)]
    return a


POP = {"cell_call_key": "cellbender_cell", "umi_floor": 100.0,
       "gene_floor": 5.0, "mito_ceiling": 25.0}

print("\nthe population parameter")
try:
    so._op_cluster(cohort(), {"resolution": 1.0, "seed": 0, "doublet_key": None}, "/tmp/x")
    check("an ABSENT population refuses", False, "it clustered anyway")
except so.TaskFailure as e:
    check("an ABSENT population refuses", "population" in str(e).lower(), str(e)[:80])
except Exception as e:  # noqa: BLE001
    check("an ABSENT population refuses", False, f"{type(e).__name__}: {e}")

print("\nthe mask")
a = cohort()
n_before = a.n_obs
keep = (a.obs["cellbender_cell"].to_numpy()
        & (a.obs["total_counts"].to_numpy() >= POP["umi_floor"])
        & (a.obs["n_genes_by_counts"].to_numpy() >= POP["gene_floor"])
        & (a.obs["pct_counts_mt"].to_numpy() <= POP["mito_ceiling"]))
check("the fixture actually contains droplets the mask must drop",
      int(keep.sum()) < n_before, f"{keep.sum()} of {n_before}")
check("...and it is most of them - otherwise this probe proves nothing",
      int(keep.sum()) < 0.5 * n_before, f"{keep.sum()} of {n_before}")

print("\nrefusals that must not be silent")
try:
    so._op_cluster(cohort(), {"resolution": 1.0, "seed": 0, "doublet_key": None,
                              "population": {"cell_call_key": "not_a_column"}}, "/tmp/x")
    check("a cell_call_key that is not in obs refuses", False, "it clustered the whole object")
except so.TaskFailure as e:
    check("a cell_call_key that is not in obs refuses", "no such column" in str(e).lower(),
          str(e)[:80])
except Exception as e:  # noqa: BLE001
    check("a cell_call_key that is not in obs refuses", False, f"{type(e).__name__}: {e}")

try:
    so._op_cluster(cohort(), {"resolution": 1.0, "seed": 0, "doublet_key": None,
                              "population": {"umi_floor": 1e9}}, "/tmp/x")
    check("a population that leaves almost nothing refuses", False, "it clustered the remainder")
except so.TaskFailure as e:
    check("a population that leaves almost nothing refuses", "too few" in str(e).lower(),
          str(e)[:80])
except Exception as e:  # noqa: BLE001
    check("a population that leaves almost nothing refuses", False, f"{type(e).__name__}: {e}")

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
sys.exit(1 if FAIL else 0)
