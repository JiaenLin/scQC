# Known issues

Defects that are measured, reproduced and not yet fixed. A pipeline that reports its own limits is
the point of this project, so they are written here rather than left in a run log.

## 1. Step 6 clustered the wrong population — FIXED in `02a422d`, verified

`_op_cluster` opened the denoised FULL DROPLET matrix while `cluster_flags.py` specified the
quality-filtered one. It now receives a `population` — cell call, count floors and the library's
mitochondrial ceiling — and applies it before clustering. Doublets stay ATTACHED and NOT applied,
because criterion D is only computable before a removal.

**Verified on the calibration cohort:** the `population clustered` REVIEW is gone; clusters per
library are 11–15 where they were 8–21 set by droplet contamination (one library went 8 → 14,
another 21 → 14); the deliverable is unchanged, since step 6 removes nothing.

**Two things this cost, both worth repeating to anyone touching this path.** The first fix
(`96f76d4`) *never applied*: it read the count floors from `thresholds_per_sample.csv`, which a
later task writes, so the read raised and an `except … return None` turned that into "cluster
everything" — while the run exited 0 and produced a byte-identical deliverable. The floors were
already in `results_by_key["05_quality"].metrics`. Both fallbacks are now refusals.

**The acceptance test for anything here is that the step-6 finding CHANGES, not that the run
succeeds.** The broken attempt succeeded.

Residual, and now a real question rather than an artefact: five libraries of ten report zero
flagged clusters. That is no longer explained by clustering empty droplets.

## 2. The run key does not cover the code — OPEN

`engine/runkey.py` hashes mode, the samplesheet's content and the DECLARED parameters. It does not
hash the pipeline version. Two runs of different scQC versions over the same inputs therefore
produce the SAME key unless a declared parameter also changed. On the calibration cohort two
genuinely different deliverables differed in key only because an interpreter path — itself a
declared parameter — had changed; that was luck, not design.

The header is honest that "a digest is not a guarantee of reproducibility" and that the key
answers *was this asked for in the same way*. That remains true, and it is still too easy to end
up with one identity over two results. Recording the commit inside `INPUTS.json`, or mixing it
into the digest, would close it.

## 3. Environments are not relocatable — OPEN, bit once

Conda R bakes absolute paths in at build time. Moving `env/rdoublet` left the old path in `bin/R`,
the `Rscript` binary and 32 other files, and **R could not start at all** — so the doublet step was
unrunnable, silently, until the next run was attempted. `R_HOME` does not override it.

`conf/env/install_rdoublet.sh` also pins `r-base` and `r-xgboost` but **not**
`bioconductor-scdblfinder`, so rebuilding to repair a move can change the doublet caller — which
is not something to do in the middle of measuring a different change. Pin it.

## 4. Smaller, all confirmed

- `bin/scqc` is `#!/usr/bin/env python3`; the orchestrator therefore runs under the system
  interpreter, which may lack `h5py`. `--python` governs subprocess tasks only. The failure
  arrives 40 lines into a run that looks healthy. Either document invoking it with the core
  interpreter, or have it re-exec.
- The CLI hardcodes `samplesheet.csv`. `docs/QUICKSTART.md` says `samplesheet.tsv`. The delimiter
  is auto-detected; the filename is not.
- `03_light_floor` has no `STEP_TEXT` entry, so the report's "what it does / what it cannot
  establish" is blank for step 3 — the one place a step's stated limit is absent rather than
  marked absent.
- Four `adapters/*.pyc.tmp` build artefacts are tracked in git.
