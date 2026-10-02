"""Synthetic cohorts in the shapes an aligner actually delivers, for tests that drive the CLI.

Every library is three populations - barcodes with a handful of UMI, droplets an aligner would
not call, and cells - because a raw matrix is all three, and a pipeline tested on the cells alone
has never met the matrix it will be given. Options cover the conditions a reusable QC tool meets:
species (gene naming), input format (MatrixMarket directory or 10x HDF5), extra feature types
(antibody capture on a CITE-seq run), duplicated gene symbols (real references carry some), and
depth. Nothing here is any cohort's data or any cohort's numbers.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HUMAN = {"mt": ["MT-CO1", "MT-ND1", "MT-ND2"], "ribo": ["RPL3", "RPS6", "RPL13"],
         "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]", "species": "homo_sapiens", "id": "ENSG"}
MOUSE = {"mt": ["mt-Co1", "mt-Nd1", "mt-Nd2"], "ribo": ["Rpl3", "Rps6", "Rpl13"],
         "mt_prefix": "mt-", "ribo_pattern": "^Rp[sl]", "species": "mus_musculus", "id": "ENSMUSG"}


def genes(spec, n_other=1200, duplicate=False):
    other = [f"Gene{i}" if spec is MOUSE else f"GENE{i}" for i in range(n_other)]
    if duplicate:
        other[-1] = other[-2]            # one symbol twice, as real references have
    return spec["mt"] + spec["ribo"] + other


def library(seed, n_genes, *, n_noise=400, n_cells=160, n_uncalled=40, cell_umi=(1500, 6000),
            mito=(0.02, 0.12)):
    """(counts, barcodes, called): droplets x genes, the first three genes mitochondrial."""
    rng = np.random.default_rng(seed)
    rows, bcs, called = [], [], []

    def droplet(total, mito_frac):
        p = np.full(n_genes, (1 - mito_frac) / (n_genes - 3))
        p[:3] = mito_frac / 3
        return rng.multinomial(int(total), p)

    for i in range(n_noise):
        rows.append(droplet(rng.integers(1, 16), 0.05)); bcs.append(f"AAN{i:05d}-1")
    for i in range(n_cells):
        rows.append(droplet(rng.integers(*cell_umi), rng.uniform(*mito)))
        bcs.append(f"CCC{i:05d}-1"); called.append(bcs[-1])
    for i in range(n_uncalled):
        rows.append(droplet(rng.integers(400, 900), 0.45)); bcs.append(f"GGU{i:05d}-1")
    return sp.csr_matrix(np.vstack(rows)), bcs, called


def write_mtx(d: Path, X, bcs, names, spec, types=None):
    d.mkdir(parents=True, exist_ok=True)
    G = sp.coo_matrix(X.T)
    with gzip.open(d / "matrix.mtx.gz", "wt") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{G.shape[0]} {G.shape[1]} {G.nnz}\n")
        for r, c, v in zip(G.row, G.col, G.data):
            fh.write(f"{r + 1} {c + 1} {int(v)}\n")
    with gzip.open(d / "barcodes.tsv.gz", "wt") as fh:
        fh.write("".join(f"{b}\n" for b in bcs))
    types = types or ["Gene Expression"] * len(names)
    with gzip.open(d / "features.tsv.gz", "wt") as fh:
        fh.write("".join(f"{spec['id']}{i:011d}\t{g}\t{t}\n"
                         for i, (g, t) in enumerate(zip(names, types))))


def write_h5(path: Path, X, bcs, names, spec, types=None):
    import h5py
    path.parent.mkdir(parents=True, exist_ok=True)
    G = sp.csc_matrix(X.T)                                   # features x barcodes
    types = types or ["Gene Expression"] * len(names)
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        m["barcodes"] = np.array(bcs, dtype="S")
        m["data"] = G.data.astype("int32")
        m["indices"] = G.indices.astype("int64")
        m["indptr"] = G.indptr.astype("int64")
        m["shape"] = np.array(G.shape, dtype="int64")
        ft = m.create_group("features")
        ft["id"] = np.array([f"{spec['id']}{i:011d}" for i in range(len(names))], dtype="S")
        ft["name"] = np.array(names, dtype="S")
        ft["feature_type"] = np.array(types, dtype="S")
        ft["genome"] = np.array([b"ref"] * len(names))
        ft["_all_tag_keys"] = np.array([b"genome"])


def with_antibodies(X, n_adt, seed):
    """Append `n_adt` antibody-capture features, heavily counted, as a CITE-seq run would."""
    rng = np.random.default_rng(seed + 1000)
    # Scaled to each droplet's RNA, so an empty droplet carries antibody background and a cell
    # carries hundreds - as a CITE-seq run does - rather than every droplet carrying a cell's worth.
    lam = np.maximum(np.asarray(X.sum(axis=1)).ravel() / 20.0, 0.3)
    adt = rng.poisson(lam[:, None], size=(X.shape[0], n_adt))
    return sp.hstack([X, sp.csr_matrix(adt)]).tocsr()


def cellranger_library(root: Path, name: str, seed: int, spec=HUMAN, *, fmt="mtx",
                       n_adt=0, duplicate=False, **kw):
    """A library laid out as Cell Ranger lays it out. Returns (raw path, filtered path, called)."""
    names = genes(spec, duplicate=duplicate)
    X, bcs, called = library(seed, len(names), **kw)
    types = None
    if n_adt:
        X = with_antibodies(X, n_adt, seed)
        names = names + [f"CD{i}_TotalSeqB" for i in range(n_adt)]
        types = ["Gene Expression"] * (len(names) - n_adt) + ["Antibody Capture"] * n_adt
    keep = [bcs.index(b) for b in called]
    lib = root / name
    if fmt == "h5":
        raw = lib / "raw_feature_bc_matrix.h5"
        write_h5(raw, X, bcs, names, spec, types)
        filt = lib / "filtered_feature_bc_matrix.h5"
        write_h5(filt, X[keep], called, names, spec, types)
    else:
        raw = lib / "raw_feature_bc_matrix"
        write_mtx(raw, X, bcs, names, spec, types)
        filt = lib / "filtered_feature_bc_matrix"
        write_mtx(filt, X[keep], called, names, spec, types)
    return raw, filt, called
