# Writes a scQC project - per-sample count matrices, a samplesheet, a registry - from one object of
# cells. It filters nothing and decides nothing: every cell of the object is written, and the run
# that follows applies the pipeline's own rules.
"""A project from the development suite's two-shape fixture, so the fixture reaches the pipeline.

WHY THIS EXISTS (harness autodev W020, B1-041). The harness's development suite hands every tool
ONE object - its synthetic cohort, written twice under different column names - and scQC consumes a
PROJECT: a samplesheet and one matrix per library. So the suite's fixture tier ran `scqc describe`,
which checks the declaration and nothing of the arithmetic, and the step point's only test was an
acceptance runner that read one calibration cohort's directory layout and skipped without it. The
fixture never reached a step.

This splits the object into the project scQC is built to read, and nothing more:

  * one 10x-layout MatrixMarket matrix per library, from the counts the caller NAMES
    (`counts_layer`), grouped by the sample column the caller NAMES (`sample_key`) - neither is
    assumed, so the two shapes, which spell both differently, are both told and both run;
  * that library's own barcodes as its called cells (`aligner_cells`): the object is an object of
    cells - a cell somebody called is what each row is - and saying so is the honest cell call. Its
    rows are not droplets, and the pipeline is never asked to treat them as droplets;
  * a samplesheet: `assay` is `scrna` (the object holds cells, not nuclei), `platform` is
    `synthetic` (no aligner made it, and none is needed - a verified matrix from any platform is
    accepted), the mitochondrial and ribosomal patterns the caller gives, the doublet rate the
    caller declares, and every column of the caller's design table, under its own name;
  * a registry naming the synthetic reference.

Library names are the sample column's values; cell names are the object's own `obs_names`. A design
row with no cells in the object is left out and said to be. Nothing here is any cohort's data.
"""
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

from engine.task import TaskFailure

#: The registry key every fixture library resolves against. Synthetic, and named as such.
SPECIES, BUILD = "synthetic", "fixture"
#: What a library made by no aligner declares as its platform. Any declared platform is accepted
#: with a verified matrix (modules/00_ingest, B1-048), and this one says what it is.
PLATFORM = "synthetic"


def _safe(name: str) -> str:
    """A directory name for a library: its own name, with only path-hostile characters replaced."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("._")
    if not s:
        raise TaskFailure(f"sample name {name!r} leaves no usable directory name")
    return s


def _write_mtx(d: Path, X, barcodes, genes) -> None:
    """A 10x-layout MatrixMarket directory: features x barcodes, integer counts."""
    import numpy as np
    import scipy.sparse as sp

    d.mkdir(parents=True, exist_ok=True)
    G = sp.coo_matrix(sp.csr_matrix(X).T)
    vals = np.asarray(G.data)
    if vals.size and (not np.all(np.isfinite(vals)) or np.any(vals != np.rint(vals))):
        raise TaskFailure("the counts layer holds non-integer values; a count matrix is written from "
                          "counts, and the layer named is not one")
    with gzip.open(d / "matrix.mtx.gz", "wt") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{G.shape[0]} {G.shape[1]} {G.nnz}\n")
        for r, c, v in zip(G.row, G.col, vals):
            fh.write(f"{r + 1} {c + 1} {int(round(float(v)))}\n")
    with gzip.open(d / "barcodes.tsv.gz", "wt") as fh:
        fh.write("".join(f"{b}\n" for b in barcodes))
    with gzip.open(d / "features.tsv.gz", "wt") as fh:
        fh.write("".join(f"{g}\t{g}\tGene Expression\n" for g in genes))


def build_project(observations, out, *, sample_key: str, counts_layer: str, mt_prefix: str,
                  ribo_pattern: str, dbr: float, design=None) -> dict:
    """Write the project under `out`. Returns the paths and what was written and left out."""
    import anndata as ad

    out = Path(out)
    A = ad.read_h5ad(str(observations))
    if sample_key not in A.obs.columns:
        raise TaskFailure(f"obs[{sample_key!r}] is absent; the sample column must be named, and the "
                          f"object carries {', '.join(map(str, A.obs.columns))}")
    if counts_layer not in A.layers:
        raise TaskFailure(f"layer {counts_layer!r} is absent; the counts must be named, and the "
                          f"object carries {', '.join(sorted(str(k) for k in A.layers.keys()))}")
    X = A.layers[counts_layer]
    genes = [str(g) for g in A.var_names]
    samples = A.obs[sample_key].astype(str)
    levels = [s for s in dict.fromkeys(samples.tolist()) if s and s.lower() != "nan"]

    design_rows, left_out = {}, []
    if design:
        with open(design, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                key = str(r.get(sample_key, "")).strip()
                if not key:
                    raise TaskFailure(f"{design} has no {sample_key!r} column to join the design on")
                if key not in levels:
                    left_out.append(key)
                    continue
                design_rows[key] = {k: v for k, v in r.items() if k != sample_key}
    factors = sorted({k for r in design_rows.values() for k in r})

    rows = []
    data = out / "data"
    for s in levels:
        idx = (samples == s).to_numpy()
        names = [str(b) for b in A.obs_names[idx]]
        lib = data / _safe(s)
        _write_mtx(lib / "matrix", X[idx], names, genes)
        called = lib / "called"
        called.mkdir(parents=True, exist_ok=True)
        with gzip.open(called / "barcodes.tsv.gz", "wt") as fh:
            fh.write("".join(f"{b}\n" for b in names))
        row = {"sample": s, "platform": PLATFORM, "species": SPECIES,
               "reference": f"{SPECIES}/{BUILD}", "assay": "scrna",
               "matrix": str(lib / "matrix"), "aligner_cells": str(called),
               "mt_prefix": mt_prefix, "ribo_pattern": ribo_pattern,
               "dbr": f"{float(dbr):g}", "dbr_sd": f"{float(dbr):g}"}
        for f in factors:
            row[f] = (design_rows.get(s) or {}).get(f, "")
        rows.append(row)

    cols = list(rows[0]) if rows else []
    sheet = out / "samplesheet.tsv"
    sheet.write_text("\t".join(cols) + "\n" + "".join(
        "\t".join(str(r.get(c, "")) for c in cols) + "\n" for r in rows), encoding="utf-8")
    reg = out / "registry.tsv"
    reg.write_text(f"species\tbuild\tpath\n{SPECIES}\t{BUILD}\t-\n", encoding="utf-8")
    return {"samplesheet": sheet, "registry": reg, "libraries": levels,
            "cells": int(A.n_obs), "design_rows_without_cells": left_out}
