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
  * a samplesheet: `assay` is what the caller DECLARES (`--assay`, required, no default: the
    suite's object holds cells, so its declaration says `scrna`; a nucleus declaration is refused
    by ingest, which requires the raw droplets nuclei need and an object of cells does not hold),
    `platform` is `synthetic` (no aligner made it, and none is needed - a verified matrix from any platform is
    accepted), the mitochondrial and ribosomal patterns the caller gives, the doublet rate the
    caller declares, and every column of the caller's design table, under its own name;
  * a registry naming the synthetic reference.

The project lives beside the suite's `--out` (default `<out>.project`), and `--out` receives a flat
copy of the run's text products, because the project's run directories are named by a digest of the
samplesheet and the two shapes' samplesheets differ by construction.

Library names are the sample column's values; cell names are the object's own `obs_names`. A design
row with no cells in the object is left out and said to be. Nothing here is any cohort's data.

THE VALUES MUST REACH THE STEPS THAT READ THEM (harness autodev W048, GAP-023). A run that completes
shows that the pipeline carried the fixture; it does not show that what the caller TOLD it - the
assay, the design's factors, the mitochondrial prefix, the ribosomal pattern - is what the steps
read. A step that ignored the told prefix would still complete, with no mitochondrial fraction at
all. So `axes_reached` reads the run's own records after it and names every told value that did
not arrive: the assay in the run's samplesheet, its report and the record of every step that
writes the assay it ran under (each gate's `assay_received`, step 5's `mito_assay`); each declared
factor in the run's recorded parameters, samplesheet and design AND in the factors every gate
recorded it was given (`design_factors_received`, harness autodev W048 fix round 1: the samplesheet
and the parameters still say what was told when no gate received it); and each gene pattern
matching genes of the object AND giving a non-zero fraction in every library's per-cell table.
`fixture-run` fails on any miss. A limit of the fixture, not of the tool: both of the suite's
shapes are cells, so only `scrna` is told, and a step wired to `scrna` cannot be told apart here.
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


#: The assays a samplesheet may declare. No default: the caller says which (rule: no silent assay).
ASSAYS = ("scrna", "snrna")


def build_project(observations, out, *, sample_key: str, counts_layer: str, mt_prefix: str,
                  ribo_pattern: str, dbr: float, assay: str, design=None) -> dict:
    """Write the project under `out`. Returns the paths, what was written and left out, and how
    many genes of the object the told mitochondrial prefix and ribosomal pattern match."""
    import anndata as ad

    if assay not in ASSAYS:
        raise TaskFailure(f"assay {assay!r} is not one of {', '.join(ASSAYS)}; it is declared, never "
                          f"defaulted")

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
               "reference": f"{SPECIES}/{BUILD}", "assay": assay,
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
    rx = re.compile(ribo_pattern)
    return {"samplesheet": sheet, "registry": reg, "libraries": levels,
            "cells": int(A.n_obs), "design_rows_without_cells": left_out,
            "genes_mt": sum(1 for g in genes if g.startswith(mt_prefix)),
            "genes_ribo": sum(1 for g in genes if rx.match(g))}


def axes_reached(run_dir, *, assay: str, factors, mt_prefix: str, ribo_pattern: str,
                 genes_mt: int, genes_ribo: int) -> tuple:
    """(lines, misses): what the run's own records say each TOLD value reached. A miss names the
    value and where it did not arrive; an empty `misses` is the only pass."""
    import ast
    import json

    import pandas as pd

    run_dir = Path(run_dir)
    lines, misses = [], []
    try:
        described = json.loads((run_dir / "INPUTS.json").read_text())["described"]
        report = json.loads((run_dir / "reports" / "report.json").read_text())
    except (OSError, ValueError, KeyError) as e:
        return lines, [f"the run's records cannot be read ({e}), so no value can be shown to arrive"]
    rows = described.get("samplesheet") or []
    libs = [str(r.get("sample")) for r in rows]
    # WHAT EACH STEP RECORDED IT RECEIVED (engine/steps._received; step 5's `mito_assay`). The
    # samplesheet and the parameter record say what was TOLD; only these say what ARRIVED.
    rec = {}
    for st in report.get("steps") or []:
        for fd in st.get("found") or []:
            if fd.get("label") in ("assay_received", "mito_assay", "design_factors_received"):
                rec.setdefault(fd["label"], []).append((str(st.get("key")), str(fd.get("value"))))
    ran = rec.get("assay_received", []) + rec.get("mito_assay", [])
    gates = rec.get("design_factors_received", [])

    seen = sorted({str(r.get("assay") or "").strip().lower() for r in rows})
    said = str((report.get("run") or {}).get("assay") or "")
    other = [f"{k} ran under {v!r}" for k, v in ran if v != assay]
    if seen == [assay] and said == assay and ran and not other:
        lines.append(f"assay {assay!r}: every library of the run's samplesheet, its report, and the "
                     f"steps that record the assay they ran under ({', '.join(k for k, _ in ran)})")
    else:
        misses.append(f"assay {assay!r} was told; the run's samplesheet says {seen}, its report "
                      f"{said!r}, and " + ("; ".join(other) if other else
                                           "no step recorded the assay it ran under" if not ran
                                           else "every recording step agrees"))

    raw = str((described.get("parameters") or {}).get("design_factors") or "")
    try:
        declared = [str(x) for x in (ast.literal_eval(raw) if raw.startswith("[") else [])]
    except (ValueError, SyntaxError):
        declared = []
    found = [str(x) for x in ((report.get("confounding") or {}).get("factors") or [])]
    for f in factors:
        where = [w for w, ok in (("the run's parameters", f in declared),
                                 ("its samplesheet", bool(rows) and all(f in r for r in rows)),
                                 ("its design", f in found)) if not ok]
        if not gates:
            where.append("every gate's record (no gate recorded the factors it was given)")
        where += [f"the factors the {k} gate received ({v})" for k, v in gates
                  if f not in v.split(",")]
        if where:
            misses.append(f"design factor {f!r} was told and is absent from {', '.join(where)}")
        else:
            lines.append(f"design factor {f!r}: declared in the run's parameters, a column of its "
                         f"samplesheet, a factor of its design, and received by the gate(s) of "
                         f"{', '.join(k for k, _ in gates)}")

    for name, told, n_genes, col in (("mitochondrial prefix", mt_prefix, genes_mt, "pct_counts_mt"),
                                     ("ribosomal pattern", ribo_pattern, genes_ribo,
                                      "pct_counts_ribo")):
        if not n_genes:
            misses.append(f"{name} {told!r} was told and matches no gene of the object")
            continue
        zero, absent = [], []
        for s in libs:
            p = run_dir / "tables" / f"{_safe(s)}.percell.csv"
            if not p.is_file():
                p = run_dir / "tables" / f"{s}.percell.csv"
            try:
                v = pd.read_csv(p, usecols=[col])[col]
            except (OSError, ValueError):
                absent.append(s)
                continue
            if not (v.fillna(0) > 0).any():
                zero.append(s)
        if zero or absent:
            misses.append(f"{name} {told!r} matches {n_genes} gene(s) of the object, and the per-cell "
                          f"{col} is " + "; ".join(x for x in (
                              f"zero in every cell of {', '.join(zero)}" if zero else "",
                              f"absent for {', '.join(absent)}" if absent else "") if x))
        else:
            lines.append(f"{name} {told!r}: {n_genes} gene(s) of the object, a non-zero {col} in "
                         f"each of {len(libs)} libraries")
    return lines, misses
