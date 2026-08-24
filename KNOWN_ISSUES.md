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

## 1b. Five of the then-fifteen figures were never drawn — FIXED, and three of the reasons were wrong

The report printed a reason where each missing figure should have been, and the reasons read as
statements about the pipeline's design: a step "computes something, uses it, and discards it".
For three of the five that was not what had happened.

| | what the report said | what was true |
|---|---|---|
| F1 | step 0 "records only the verdict, so the curve behind it is gone" | `adapters/matrix.barcode_rank()` existed and `--rank-points` was plumbed end to end. `run_summary_stats` COUNTED the pairs into `n_rank_points` and dropped the pairs themselves. Nothing ever asked step 0 for them |
| F5 | "it needs a sweep step this pipeline does not run" | `adapters/doublets.sweep()` existed, complete with its cross-setting version checks. No task called it |
| F6 | step 5 "fits a KDE, takes the minimum, records the valley position, and discards the curve" | `_op_valley` had always written `<sample>.valley_density.csv`. It went to the scratch directory, which nothing publishes and the report never looks in |
| F10, F11 | no embedding is computed | correct — that one really was absent |

All five are now produced. Two figures were added on the same reasoning: **F13**, because step 5
derives and applies TWO count floors and one density figure can carry only one of them, and
**F14/F15**, because the report's decision spine names three applied criteria and two of its
three rows read "no figure is produced for this axis" — a report showing a third of the filter
while looking complete.

**The lesson is the wrong reasons, not the missing figures.** A plausible explanation for an
absent result reads exactly like a correct one, it is more durable than the defect because it
tells the next reader not to look, and every one of these was written in good faith by someone
looking at the same code. Two things follow, and both are now in place:

- `report/collect.py` states this in its own header, so the next person to meet an absent figure
  checks the producing step before writing down why it cannot be produced.
- `tests/test_figure_collection.py` renders every figure from a fixture written as files. A
  builder can return a dict, be recorded as "assembled", and still be unrenderable, because
  `render_figures` calls `FIGURE_FUNCTIONS[id](**data)` and a key that is not a parameter of that
  function only fails at draw time. **The acceptance test for a figure is a figure.**

Residual, and stated because it is a limit rather than a defect: the sweep behind F5 is opt-in
(`--dbr-sd-sweep default,dbr,1`). It re-scores every library once per setting and applies
nothing, so it is not imposed on a run that did not ask for it; where the figure would be, the
report names the flag that produces it.

## 1c. The payload and the report beside it could disagree — FIXED

`scqc report <results>` rebuilds a document from `payload.json`, and `finish()` wrote the HTML and
the JSON while leaving the payload as whatever the report TASK had written. On a resumed run that
task is SKIPPED, so the payload on disk was the previous version's while the report beside it was
this one's — two files meant to be two views of one thing, quietly describing different runs.

Caught on this cohort with the worst possible timing: the stale payload predated the parameter
table and the freshness feed, so running `scqc report` would have regenerated a document with **no
section 2 and freshness NOT CHECKED**, over one that had both. It exits 0 and the result opens and
reads like a correct report.

`finish()` now writes the payload it built the document from.

## 1d. Two defects a resumed run walked straight into — both FIXED

**A run could be resumed exactly once.** `should_skip` requires the recorded status to be DONE,
and the orchestrator recorded the SKIPPED result over that DONE record. The next resume found a
non-done record and re-ran the task; the one after that re-ran everything. Silent: exit 0, same
deliverable, and the only visible difference is RUN where SKIP should be, in a log of seventy
near-identical lines. On the cohort that found it, a second resume re-copied 2.5 GB of matrices,
re-scored every library for doublets and re-clustered all ten — to rebuild one report. Nothing is
recorded for a skip now; the record that justified it is already there and is the better one.

**UMAP's numba cache is shared and was not serialised.** numba caches compiled functions into the
INSTALLED PACKAGE directory, which on a cluster is one NFS path every concurrent job writes to.
Ten libraries embedding at once raced, and one died with `OSError: [Errno 116] Stale file handle`
inside `numba/core/caching.py`. The traceback names numba, so it reads as an environment fault; it
is concurrency over a shared cache nothing was serialising, and it arrived with the embedding
(figures F10/F11). Each scanpy op now gets a `NUMBA_CACHE_DIR` unique to its prefix and op. The
cost is recompilation per job; the alternative fails one library in ten, and a cohort missing a
library is a different cohort.

**The second is why the first matters.** A full re-execution nobody asked for is not just slow —
it re-exercises every concurrency hazard in the pipeline. These two were found together because
one caused the other to be reached.

## 2. The run key did not cover the code — FIXED in `d573c20`

`engine/runkey.py` hashed mode, the samplesheet's content and the DECLARED parameters, and not the
pipeline version. Two runs of different scQC versions over the same inputs therefore produced the
SAME key unless a declared parameter also changed. On the calibration cohort two genuinely
different deliverables differed in key only because an interpreter path — itself a declared
parameter — had changed; that was luck, not design.

**Why it mattered more than it looked.** A DERIVED threshold is a function of the inputs *and* of
the code that derives it, so changing a derivation produces a different answer under an identical
key. Re-running would have resolved to the previous run's directory, found every task complete,
SKIPPED all of them, and republished the old numbers under the new code — exiting zero, with a
report that opens and reads like a correct one.

The commit is now part of the description and `results/INDEX.tsv` carries a `code` column, so a
code change writes BESIDE the old run exactly as a parameter change does. The cost is accepted and
real: any commit invalidates resume for every project, including one that changed only a docstring.

**Residual.** Two runs from the same commit with UNCOMMITTED edits still share a key. `dirty` is
recorded, but a modified tree has no identity to hash cheaply. Commit before a run whose output
anyone will rely on.

## 3. Environments are not relocatable — OPEN, bit once

Conda R bakes absolute paths in at build time. Moving `env/rdoublet` left the old path in `bin/R`,
the `Rscript` binary and 32 other files, and **R could not start at all** — so the doublet step was
unrunnable, silently, until the next run was attempted. `R_HOME` does not override it.

`conf/env/install_rdoublet.sh` also pins `r-base` and `r-xgboost` but **not**
`bioconductor-scdblfinder`, so rebuilding to repair a move can change the doublet caller — which
is not something to do in the middle of measuring a different change. Pin it.

## 4. Smaller, all confirmed

- ~~**`tests/test_figure_collection.py` FAILS where it should SKIP, on a bare interpreter.**~~ —
  **FIXED.** `from report import figures as figmod` sat inside the `try` whose `except
  ModuleNotFoundError` exists to skip the *rendering* block without matplotlib. The drawable
  check further down then called `figmod.FIGURE_FUNCTIONS` unguarded, so the suite raised
  `NameError: name 'figmod' is not defined` instead of skipping. Reproduced before the fix:
  `python3 bin/scqc selftest` under a system interpreter without matplotlib gave
  `24 passed, 1 failed, 2 skipped`; after, `25 passed, 0 failed, 2 skipped`.

  `report/figures.py` imports matplotlib **inside a function**, not at module level, so the
  module import never needed guarding — only the rendering did. `figmod` is now bound before the
  `try`, and the drawable check runs on every interpreter, because it does not need matplotlib
  at all.

  **The lesson is the one this file already records twice: a suite that crashes where it means to
  skip cannot tell a reader which of the two happened.** The suite enforcing "*the acceptance
  test for a figure is a figure*" was the one whose own outcome was unreadable, and a failure
  meaning "matplotlib is absent" looked in the summary line exactly like a figure that will not
  render.
- `bin/scqc` is `#!/usr/bin/env python3`; the orchestrator therefore runs under the system
  interpreter, which may lack `h5py`. `--python` governs subprocess tasks only. The failure
  arrives 40 lines into a run that looks healthy. Either document invoking it with the core
  interpreter, or have it re-exec.
- The CLI hardcodes `samplesheet.csv`. `docs/QUICKSTART.md` says `samplesheet.tsv`. The delimiter
  is auto-detected; the filename is not.
- ~~`03_light_floor` has no `STEP_TEXT` entry~~ — FIXED, and the first attempt was only half of
  it. Adding the entry was necessary and did nothing: `step_text()` is consulted only for steps
  that have TASKS, and no task carried the key `03_light_floor` — the floor is a declared
  parameter consumed inside step 4's export. The report went on saying "no record of this step
  at all", which was accurate. It now has a task that REPORTS what the floor did per library
  (`tables/light_floor.csv`) without applying or being able to change it. Step 5's entry was
  stale in the same table — it still described the mitochondrial ceiling as Tukey's fence, which
  0.2.0 replaced with a MAD fence at a derived k, keeping Tukey only as the calibrator.

  **The lesson is the same one this file already records about step 6:** a fix that leaves the
  symptom in place was not a fix, and "I edited the thing named in the defect" is not evidence.
  The acceptance test is that the REPORT changes.
- ~~Four `adapters/*.pyc.tmp` build artefacts are tracked in git.~~ — FIXED. Removed from
  the index and `*.pyc.tmp` added to `.gitignore`. They are half-written byte-code CPython
  leaves behind when a compile is interrupted, and they were being shipped to every user.
