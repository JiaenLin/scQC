# Builds a tiny denoiser-format .h5 and converts it; runs nothing external and needs no GPU.
"""The run route of the denoiser, which no cohort took until a second one did.

Every cohort before it arrived with its ambient correction already done, so `01_ambient` ran on
the supplied route only. Read before the first run on the other route (single-cell-harness
ADR-0027), three things could not have worked:

  Q4  CellBender wrote its own .h5 under `<sample>_ambient.h5`, the name every consumer opens as
      an AnnData - step 5 would have failed on the first library ever denoised inside this
      pipeline, after the GPU hours were spent.
  Q5  its cell-barcode CSV would have been `<sample>_ambient_cell_barcodes.csv`, while step 2
      looks for `<sample>_cellbender_cell_barcodes.csv`.
  --  every child job went to one queue, with ngpus=1 added for the denoiser; a site whose GPUs
      are in GPU queues had no way to send it there.

Each is checked here, and each check fails on the code that preceded it.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
fails: list[str] = []

# ---- 1. the conversion: the denoiser's format in, the object every step opens out
try:
    import h5py
    import numpy as np
    import scipy.sparse as sp
    import anndata as ad
    HAVE = True
except ImportError:
    HAVE = False
    print("SKIP 1: h5py / numpy / scipy / anndata not importable here")

if HAVE:
    from adapters import scanpy_ops as so
    rng = np.random.default_rng(0)
    n_bc, n_g = 40, 25
    X = sp.random(n_bc, n_g, density=0.3, random_state=0, format="csr")
    X.data = np.round(X.data * 9) + 1
    X = X.tolil()
    for i in range(30, 40):          # the droplets the denoiser called empty: zero counts
        X[i, :] = 0
    X = X.tocsr()
    G = X.T.tocsc()                  # 10x stores features x barcodes, column-compressed
    with tempfile.TemporaryDirectory() as td:
        native = Path(td) / "S1_cellbender.h5"
        with h5py.File(native, "w") as f:
            m = f.create_group("matrix")
            m["barcodes"] = np.array([f"AAAC{i:04d}-1" for i in range(n_bc)], dtype="S")
            m["data"] = G.data.astype("int32")
            m["indices"] = G.indices.astype("int64")
            m["indptr"] = G.indptr.astype("int64")
            m["shape"] = np.array([n_g, n_bc], dtype="int64")
            ft = m.create_group("features")
            ft["id"] = np.array([f"ENSG{i:011d}" for i in range(n_g)], dtype="S")
            ft["name"] = np.array(["MT-CO1", "MT-ND1"] + [f"G{i}" for i in range(n_g - 2)], dtype="S")
            ft["feature_type"] = np.array([b"Gene Expression"] * n_g)
            ft["genome"] = np.array([b"GRCh38"] * n_g)
        dest = Path(td) / "S1_ambient.h5"
        try:
            outs, met = so._op_convert(None, {"source": str(native), "dest": str(dest)},
                                       str(Path(td) / "S1_convert"))
            obj = ad.read_h5ad(dest)                  # what `_load` does at step 5
            tot = np.asarray(obj.X.sum(axis=1)).ravel()
            if obj.shape != (n_bc, n_g):
                fails.append(f"1: converted shape {obj.shape}, expected {(n_bc, n_g)}")
            if int((tot > 0).sum()) != 30 or int((tot == 0).sum()) != 10:
                fails.append(f"1: the zero-count droplets did not survive: {int((tot > 0).sum())} with counts")
            if not str(obj.var_names[0]).startswith("MT-"):
                fails.append(f"1: var_names are not symbols ({obj.var_names[0]!r}); the mito prefix matches on them")
            if met.get("droplets") != n_bc:
                fails.append(f"1: metrics say {met.get('droplets')} droplets")
        except Exception as e:                                        # noqa: BLE001
            fails.append(f"1: conversion failed: {type(e).__name__}: {e}")
    print("1. the denoiser's .h5 converts to the object step 5 opens, empties kept at zero")

# ---- 2. the names line up: what CellBender writes is what step 2 reads
from adapters import cellbender as cbd  # noqa: E402
from engine import graph  # noqa: E402


class _P:
    results = Path("/tmp/r"); work = Path("/tmp/w"); project = Path("/tmp/p")
    decisions: dict = {}
    mode = "evidence"
    samples = [{"sample": "S1", "platform": "10x", "species": "homo_sapiens",
                "reference": "homo_sapiens/GRCh38", "assay": "scrna", "matrix": "/d/S1",
                "aligner_cells": "/d/S1_filtered", "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]"}]


tasks = graph.main_stage(_P, "python", {}, {"S1": {"mode": "accept"}})
amb = next(t for t in tasks if t.key == "01_ambient/S1")
cells = next(t for t in tasks if t.key == "02_cells")
src = (ROOT / "engine" / "steps.py").read_text(encoding="utf-8")
body = src.split("def _ambient(", 1)[1].split("\ndef ", 1)[0]
if 'f"{sample}_cellbender.h5"' not in body or "output_h5=native" not in body:
    fails.append("2: the denoiser is not handed `<sample>_cellbender.h5` as its own output")
want_csv = str(cbd.expected_products(_P.results / "objects" / "S1_cellbender.h5")["cell_barcodes_csv"])
got_csv = (cells.params.get("call_paths") or {}).get("S1", {}).get("cellbender")
if got_csv != want_csv:
    fails.append(f"2: step 2 reads {got_csv}, CellBender writes {want_csv}")
if amb.outputs != (str(_P.results / "objects" / "S1_ambient.h5"),):
    fails.append(f"2: step 1 declares {amb.outputs}, not the object every consumer opens")
if not amb.params.get("python_exe"):
    fails.append("2: step 1 is not given the analysis interpreter its conversion runs under")
print("2. CellBender's barcode CSV is the one step 2 reads; step 1 declares the converted object")

# ---- 3. a GPU task goes to the GPU queue, when one is named
from engine import executor as ex  # noqa: E402
e = ex.PBSExecutor.__new__(ex.PBSExecutor)
e.queue, e.gpu_queue = "cpuq", "gpuq"
if not hasattr(e, "queue_for"):
    fails.append("3: the executor has no notion of a queue per task - every job goes to one queue")
else:
    if (e.queue_for(True), e.queue_for(False)) != ("gpuq", "cpuq"):
        fails.append(f"3: routed {(e.queue_for(True), e.queue_for(False))}")
    e.gpu_queue = None
    if e.queue_for(True) != "cpuq":
        fails.append("3: with no GPU queue named, a GPU task must go where every task goes")
if "queue = self.queue_for(gpu)" not in (ROOT / "engine" / "executor.py").read_text(encoding="utf-8"):
    fails.append("3: the job script does not take its queue from queue_for()")
print("3. --gpu-queue routes the denoiser; absent, nothing changes")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("ambient run path OK")
