# Profiles small synthetic clusterings with adapters.scanpy_ops.cluster_profile; needs the analysis stack.
"""A cluster too small to rank is not evaluated, and step 6 goes on (harness autodev W020, B1-045).

`rank_genes_groups` cannot rank a cluster of under three cells, and that used to fail the step-6
task and with it the run. Since `0e6a5f6` the small clusters are left out of the ranking and carry
criterion C as UNKNOWN; what remained was the case where NO cluster is large enough, which still
raised. Checked here:

  1  the register's CHECK: a 2-cell cluster beside larger ones - rows for every cluster, that
     cluster's `pct_uninformative` None (criterion C unknown), the others measured
  2  every cluster under three cells - rows for every cluster, C unknown on all, a note, no raise
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
fails: list[str] = []

try:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from adapters import scanpy_ops as so
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")


def obj(sizes):
    rng = np.random.default_rng(3)
    n = sum(sizes)
    genes = ["MT-CO1", "MT-ND1", "RPL3", "RPS6"] + [f"G{i}" for i in range(60)]
    X = rng.poisson(3.0, size=(n, len(genes))).astype("float32")
    labels = np.concatenate([[str(k)] * s for k, s in enumerate(sizes)])
    a = ad.AnnData(X=X, obs=pd.DataFrame({"leiden": pd.Categorical(labels)},
                                         index=[f"C{i}" for i in range(n)]),
                   var=pd.DataFrame(index=genes))
    a.obs["total_counts"] = X.sum(axis=1)
    a.obs["pct_counts_mt"] = 100.0 * X[:, :2].sum(axis=1) / np.maximum(X.sum(axis=1), 1)
    import scanpy as sc
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    return a


UNINF = {"mt": ["MT-CO1", "MT-ND1"], "ribo": ["RPL3", "RPS6"]}

if HAVE:
    for label, sizes, unknown in (("1 a 2-cell cluster beside larger ones", (30, 25, 2), {"2"}),
                                  ("2 every cluster under three cells", (2, 1, 2), {"0", "1", "2"})):
        try:
            rows = so.cluster_profile(obj(sizes), "leiden", UNINF, sample="S")
        except Exception as e:                                            # noqa: BLE001
            fails.append(f"{label}: cluster_profile raised {type(e).__name__}: {e}")
            continue
        got = {r["cluster"]: r["pct_uninformative"] for r in rows}
        if set(got) != {str(k) for k in range(len(sizes))}:
            fails.append(f"{label}: rows for {sorted(got)}, not every cluster")
        if any(got[c] is not None for c in unknown):
            fails.append(f"{label}: criterion C was given a value for a cluster too small to rank")
        if any(got[c] is None for c in set(got) - unknown):
            fails.append(f"{label}: a cluster large enough to rank was left unmeasured")
        print(f"{label}: {len(rows)} rows, C unknown on {sorted(c for c in got if got[c] is None)}")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("a cluster too small to rank is not evaluated, and the profile is still written")
