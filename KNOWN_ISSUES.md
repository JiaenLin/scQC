# Known issues

Defects that are measured, reproduced and not yet fixed. A pipeline that reports its own limits is
the point of this project, so they are written here rather than left in a run log.

## 1. Step 6 clusters the wrong population — STILL OPEN. A fix is committed and DOES NOT WORK

**Read this before assuming the mask is live.** Commit `96f76d4` added a `population` mask to
`_op_cluster` and wired `_population_for_cluster()` in `engine/steps.py` to supply it. On a full
re-run of the calibration cohort (2026-08-11) the mask **did not take effect**: step 6 reported
the identical "1413 of 1531 clusters have a median UMI of zero", and the deliverable came out
byte-identical. `_population_for_cluster()` returned None and the op fell back to clustering
everything.

The cause is in the fix, not the pipeline. `_population_for_cluster()` wraps its table reads in
`except (OSError, KeyError, ValueError): return None`, so a failed read is indistinguishable from
a deliberate "cluster everything" — the exact anti-pattern this project exists to prevent, sitting
inside the repair for it. It is why the run looked healthy while doing the old thing.

**Next step:** replace that `return None` with a refusal naming the file and the missing key, then
re-run. The op's own refusals — absent population, unknown cell-call column, over-strict mask —
are correct and covered by `tests/test_cluster_population.py`; only the supplier is broken. Check
`_tables(pipeline)` at step-6 execution time first: `mito_ceiling_per_sample.csv` has no scope row
while `thresholds_per_sample.csv` does, and the two readers treat them differently.

Everything below describes the defect itself and remains accurate.

## 1. Step 6 clusters the wrong population — OPEN, affects every cohort

**What happens.** `engine/steps.py::_cluster()` opens `results/objects/<sample>_ambient.h5`, the
denoised FULL DROPLET matrix: no cell call, no count floor. `modules/06_cluster_check/
cluster_flags.py` specifies the opposite in its own docstring — *"On the step-5 object:
quality-filtered, with the doublet flags ATTACHED and NOT applied."* The pipeline reports the
mismatch at runtime (`06_cluster_check / population clustered`) instead of refusing.

**Why it happens.** An ordering gap, not a typo. Step 6 needs a filtered population and nothing
builds one until step 7, and the doublet calls must still be attached-not-applied at that point
because criterion D is uncomputable after a removal.

**What it costs**, measured on the calibration cohort (10 libraries, 305,425 barcodes):

- **1,398 of 1,531 clusters (91.3%) contain no cell that reaches the deliverable.**
- Leiden spends a fixed resolution over what it is given, so the real cells are left
  UNDER-RESOLVED: one library received 245 clusters, of which **8** hold retained cells.
- Step 6's criteria are cluster MEDIANS and cannot fire on a blob that averages a library
  together. **The five libraries raising zero flags are the five with the fewest real clusters** —
  the dirtier the library, the fewer warnings it produces, which is backwards.
- On that cohort four of those five sat in one arm of the design, so the bias ran along a design
  factor and would have reached the study's primary readout.

Any figure of the form "N cells in M flagged clusters" inherits this and should not be quoted.

**The fix.** Step 5 already derives every floor and the per-library ceiling. Hand step 6 a
KEEP-MASK — cell call AND count floors AND mitochondrial ceiling — and apply it in memory before
`cluster()` runs. Doublets stay ATTACHED and are NOT applied. No object is written and nothing is
removed, so step 7 remains the only place a removal happens.

**Cost of fixing.** It regenerates the cluster flags, so any cohort already processed needs a
re-run before its flags are usable. Cell counts are unaffected — step 6 removes nothing.

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
