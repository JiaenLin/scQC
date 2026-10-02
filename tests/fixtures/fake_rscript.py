#!/usr/bin/env python3
# A stand-in for `Rscript adapters/scdblfinder.R ...`, for tests that drive the CLI without R.
"""Honours the scDblFinder adapter's contract and nothing more: it reads the exported barcodes,
calls every 20th one a doublet (deterministically, so a run is reproducible), writes the calls CSV
the adapter reads, and prints the version and metric lines the adapter cross-checks - the rate,
the seed and the counts it was asked for. It scores nothing: a test using it checks the pipeline
around the detector, never the detector.

    fake_rscript.py SCRIPT MTX_DIR OUT_CSV DBR SD_TOKEN SEED [key=value ...]
"""
import gzip
import os
import sys

_, _script, mtx, out_csv, dbr, sd_token, seed, *_rest = sys.argv
path = os.path.join(mtx, "barcodes.tsv")
opener = open
if not os.path.exists(path):
    path, opener = path + ".gz", gzip.open
with opener(path, "rt") as fh:
    bcs = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
cls = ["doublet" if i % 20 == 19 else "singlet" for i in range(len(bcs))]
for name, ver in (("R", "fake 0"), ("Matrix", "0"), ("SingleCellExperiment", "0"),
                  ("scDblFinder", "0"), ("xgboost", "0")):
    print(f"##scqc-version\t{name}\t{ver}")
nd = cls.count("doublet")
sd_used = "package-default" if sd_token == "default" else sd_token
for k, v in (("n_cells", len(bcs)), ("n_genes", 0), ("n_doublets", nd),
             ("n_singlets", len(bcs) - nd), ("dbr_used", dbr), ("dbr_sd_used", sd_used),
             ("dbr_sd_formal_default", "none"), ("dbr_sd_formal_default_route", "fake"),
             ("seed_used", seed), ("threads_used", 1), ("out_csv", out_csv)):
    print(f"##scqc-metric\t{k}\t{v}")
os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
with open(out_csv + ".partial", "w") as fh:
    fh.write("barcode,doublet_score,doublet_class\n")
    for i, (b, c) in enumerate(zip(bcs, cls)):
        fh.write(f"{b},{0.9 if c == 'doublet' else 0.05 + (i % 7) / 100:.4f},{c}\n")
os.replace(out_csv + ".partial", out_csv)
