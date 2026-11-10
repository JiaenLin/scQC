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

# ---- 3. a LIBRARY too small to cluster is profiled as one group, and the op does not fail
# (harness autodev W020: the development suite's `tiny_sample` hazard, 7 cells, failed step 6)
if HAVE:
    import tempfile

    def tiny(n=7):
        rng = np.random.default_rng(5)
        genes = ["MT-CO1", "MT-ND1", "RPL3"] + [f"G{i}" for i in range(40)]
        X = rng.poisson(4.0, size=(n, len(genes))).astype("float32")
        a = ad.AnnData(X=X, obs=pd.DataFrame(index=[f"C{i}" for i in range(n)]),
                       var=pd.DataFrame(index=genes))
        a.obs["cell"] = True
        return a

    with tempfile.TemporaryDirectory() as td:
        try:
            outs, met = so._op_cluster(
                tiny(), {"sample": "S", "resolution": 2.0, "seed": 0, "doublet_key": None,
                         "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]",
                         "uninformative_classes": ["mt"],
                         "population": {"cell_call_key": "cell"},
                         "embedding": {"population": "cell_called"},
                         "extra_resolutions": [1.0, 3.0]},
                str(Path(td) / "S"))
            prof = Path(td) / "S.cluster_profile.csv"
            if not prof.is_file() or met.get("n_clusters") != 1 or met.get("n_profile_rows") != 1:
                fails.append(f"3: a 7-cell library was not profiled as one group: {met}")
            if (met.get("embedding") or {}).get("skipped") is None or \
                    (Path(td) / "S.embedding.csv").exists():
                fails.append("3: a 7-cell library was embedded, or its skipped embedding not recorded")
            if not any("too few" in n or "under the" in n for n in met.get("profile_notes") or []):
                fails.append(f"3: the profile does not say why the library was not clustered: "
                             f"{met.get('profile_notes')}")
            if "ONE group" not in str(met.get("profiled_as_one_group") or ""):
                fails.append(f"3: the op's metrics do not carry profiled_as_one_group, which step 6 "
                             f"reports: {sorted(met)}")
        except Exception as e:                                            # noqa: BLE001
            fails.append(f"3: _op_cluster failed on a 7-cell library: {type(e).__name__}: {e}")
    print("3 a 7-cell library: profiled as one group, not clustered, not embedded, and said so")

    # ---- 4. a SEVERE DECLARATION on a library that called enough cells is still refused: the
    # library is not small, the floors are - one group would hide a declaration that emptied it.
    # And a population that leaves no cell at all is refused whatever the library's size.
    for label, pop in (("4a floors that leave 10 of 120 called cells",
                        {"cell_call_key": "cell", "umi_floor": None}),
                       ("4b a population that leaves no cell", {"cell_call_key": "cell",
                                                                "umi_floor": 1e9})):
        with tempfile.TemporaryDirectory() as td:
            a = tiny(120)
            if pop["umi_floor"] is None:
                # the depth of every cell is known after qc_metrics; pick a floor only 10 clear
                so.qc_metrics(a, "MT-", "^RP[SL]")
                pop = {**pop, "umi_floor": float(np.sort(a.obs["total_counts"].to_numpy())[-10])}
            try:
                so._op_cluster(a, {"sample": "S", "resolution": 2.0, "seed": 0, "doublet_key": None,
                                   "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]",
                                   "uninformative_classes": ["mt"], "population": pop,
                                   "embedding": {"population": "cell_called"}},
                               str(Path(td) / "S"))
                fails.append(f"{label}: profiled instead of refused")
            except so.TaskFailure as e:
                if "too few" not in str(e):
                    fails.append(f"{label}: refused for another reason: {str(e)[:120]}")
    print("4 a declaration that empties a library of enough called cells is refused, as before")

    # ---- 5. step 6 REPORTS the library profiled as one group: a REVIEW finding naming it.
    import types
    from engine import steps
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "tables").mkdir()
        recorded = []
        profiles = {}
        for s, n in (("BIG", 120), ("SMALL", 7)):
            outs, met = so._op_cluster(
                tiny(n), {"sample": s, "resolution": 2.0, "seed": 0, "doublet_key": None,
                          "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]",
                          "uninformative_classes": ["mt"], "n_pcs": 10,
                          "population": {"cell_call_key": "cell"}, "embedding": None},
                str(td / s))
            profiles[f"06_cluster/{s}"] = types.SimpleNamespace(outputs=[str(o) for o in outs],
                                                                metrics=met)
        pipe = types.SimpleNamespace(results=td, results_by_key=profiles,
                                     record_findings=lambda step, fs: recorded.extend(fs))
        task = types.SimpleNamespace(params={"samples": ["BIG", "SMALL"],
                                             "assay": {"BIG": "scrna", "SMALL": "scrna"}})
        try:
            steps._cluster_flags(task, pipe, None)
            got = [f for f in recorded if f.get("check") == "libraries too small to cluster"]
            if len(got) != 1 or got[0].get("severity") != "REVIEW" \
                    or "SMALL" not in got[0].get("message", "") \
                    or "BIG" in got[0].get("message", ""):
                fails.append(f"5: step 6 did not report the one-group library as REVIEW: {recorded}")
        except Exception as e:                                            # noqa: BLE001
            fails.append(f"5: step 6 failed over a one-group library: {type(e).__name__}: {e}")
    print("5 step 6 names the library profiled as one group, as a REVIEW finding")
    if fails:
        print("FAILED - " + "; ".join(fails))
        raise SystemExit(1)
