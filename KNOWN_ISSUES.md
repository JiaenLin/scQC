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

The code is now part of the description and `results/INDEX.tsv` carries a `code` column, so a
code change writes BESIDE the old run exactly as a parameter change does.

**2026-11-09 (harness autodev W020, B1-106): keyed by content, not by commit.** The key held the
commit, so any commit - one that changed only a docstring - invalidated resume for every project.
It now holds a digest of the output-bearing code's bytes (`engine/runkey.py`, `OUTPUT_BEARING`): a
docs- or tests-only commit keeps the key, and an UNCOMMITTED edit to code that runs moves it, which
closes the residual this section used to carry (two runs from one commit with uncommitted edits
shared a key). The commit and dirty flag are recorded beside the key in INPUTS.json and in
STATUS.json.

## 3. Environments are not relocatable — OPEN, bit once

Conda R bakes absolute paths in at build time. Moving `env/rdoublet` left the old path in `bin/R`,
the `Rscript` binary and 32 other files, and **R could not start at all** — so the doublet step was
unrunnable, silently, until the next run was attempted. `R_HOME` does not override it.

`conf/env/install_rdoublet.sh` also pinned `r-base` and `r-xgboost` but **not**
`bioconductor-scdblfinder`, so rebuilding to repair a move could change the doublet caller — which
is not something to do in the middle of measuring a different change. **Pinned** by the harness's
autodev wave W014, at the versions the calibration run's own report records for the environment that
scored its doublets: R 4.3.3, scDblFinder 1.16.0 (Bioconductor 3.18), SingleCellExperiment 1.24.0,
xgboost and Matrix. A package has two numbers here: the specification names conda's (`r-xgboost=1.7.6`
as a cpu build, `r-matrix=1.6_5`), and the check after the build reads R's, which is what scDblFinder
runs against - xgboost **1.7.6.1** (conda-forge's 1.7.6 installs that CRAN release) and Matrix 1.6.5 - and
reads the cpu build from the environment's conda record. Any difference is refused; an environment that
existed is not modified, and the refusal says to remove it and re-run, or build into a new prefix.
W014's first version compared R's report with conda's 1.7.6 and so refused every correct environment;
`tests/test_the_doublet_environment_is_read_as_r_reports_it.py` fails on it. The rest of that conda
environment is not pinned (section 7).

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
- ~~The CLI hardcodes `samplesheet.csv`. `docs/QUICKSTART.md` says `samplesheet.tsv`.~~ — the
  quick start names `samplesheet.csv` throughout (2026-11-09). The delimiter is auto-detected; the
  filename is not, and `--samplesheet` names another.
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

## 5. Fixed by the harness's autodev wave W007 (2026-11-09) — each with a test that fails before it

- **Steps 6 and 7 applied different thresholds when a decisions file set them.** Step 6 clustered
  on the derived floors and ceiling while step 7 applied the adjudicated ones, so the cluster
  flags step 7 carries onto every delivered cell described a population the deliverable did not
  hold. One resolver now serves both (`steps._resolve_thresholds`), and step 6's tasks carry the
  decisions' quality block, so a changed decision re-runs them instead of being skipped as
  unchanged. `tests/test_one_floor_for_steps_6_and_7.py`.
- **The report called every observation a nucleus.** Captions, axis labels, F11's question and the
  step texts now name nuclei or cells by the declared assay (`report/units.py`); a cohort with no
  single assay is named in barcodes. The spine's floor class is read from the parameter table
  (DECLARED on cells, where it printed DERIVED), and step 5's text no longer says k is derived on
  nuclei, where it is declared. Single-cell suite, sections H and D4.
- **A first cohort's reference was every project's default.** The registry scQC shipped held one
  row, the first cohort's mouse reference, and ingest read it (after the project's parent
  directory) whenever a project gave no `--registry`. Both fallbacks are gone and the file with
  them; ingest and `validate` refuse without `--registry`, and `setup/init_project.sh` writes the
  project's registry header. Wiring section N; `tests/test_a_project_brings_its_own_registry.py`.
- **The decisions template pre-filled step 6's thresholds** with one calibration cohort's numbers
  (0.5, 15.0, 50.0), so a template copied unedited replaced what step 6 proposes from the
  project's own clusters. They are empty now, and the template's unread `doublets`, `resolution`
  and `d_doublet_pct` keys are gone. Same suite, section D.
- **A padded assay stopped the build.** Ingest strips `assay`; the graph only lowered it, so
  ` scrna` passed ingest and stopped at step 6 with a KeyError (a padded and a capitalised value
  were also read as two assays). Single-cell suite, A13 and D.
- **Declared floors named the aligner's call on a library supplied already denoised.**
  `count_floors.csv`'s `applied_to` and the parameter table's basis now follow each library's
  route. Single-cell suite, A10 and G.
- **The leak guard passed over 0 terms.** It read one path from `$SCQC_FORBIDDEN_TERMS`, so two
  lists joined as PATH joins them scanned nothing and printed PASS. It reads every
  `os.pathsep`-separated list and refuses one that is missing or holds no term.
  `tests/test_leak_guard_reads_every_list.py`.
- **Two local runs whose project folders shared a name shared their scratch** and failed each
  other (`$TMPDIR/scqc_<folder name>`). Local and `SCQC_SCRATCH` scratch is named for the run;
  the PBS scratch path is unchanged. `tests/test_runs_do_not_share_scratch.py`.
- **`refuse_if_stale` had no caller** while the README said a stale report is refused. It is
  removed, and the documents say what is true: staleness is flagged on the report's face and
  counted as a defect, never refused.

## 6. Open, found or confirmed while making W007

- **Single cells get no ambient estimate and no library-level metrics** (the third item of the
  single-cell audit). Nothing estimates the ambient fraction on a library this pipeline does not
  denoise, and the report says NOT MEASURED. A capability, not a defect fixed here.
- **`decisions.yml` is not in the run key.** An apply run with a changed decisions file runs in the
  same `results/<digest>/`, re-runs the tasks of steps 5 to 7 that carry it and replaces their
  outputs. Documented (README, user guide), not changed.
- **Step 6's declared thresholds need no attestation.** A `cluster_check` value in decisions.yml
  replaces step 6's proposal as written, with no `approved_by` or `verbatim`, unlike the quality
  thresholds step 7 applies.
- **F15's caption names the light floor as the ceiling's population.** On single cells the
  ceiling is derived over the called cells at or above the declared UMI floor; the figure draws
  every barcode above the light floor.
- **`scqc report --samplesheet` reads the delimiter from the file's first line, comment or not**
  (`report/collect.collect`), so a tab-separated samplesheet that opens with a `#` line is read as
  one column, and its design and assay are lost to the rebuilt report.

## 7. The environments, pinned by the harness's autodev wave W014 (2026-11-09)

Fixed, each with a test that fails before it:

- **The denoiser's torch and numpy were never named.** `setup/install_env.sh` installed
  `cellbender==0.3.2` and let pip choose torch, and `conf/env/install_cs_cb.sh` installed torch with no
  version at all: on one cohort's first build pip chose torch 2.14.1, so the recipe built a denoiser
  other than the one this pipeline was calibrated with. The environment is one lock now,
  `conf/env/cellbender.lock.txt`, built by one recipe, `conf/env/install_cellbender.sh`, which requires
  cellbender 0.3.2, torch 2.5.1 and numpy 1.26.4 of the lock and refuses any other before anything is
  built, and refuses an environment holding any other after. Of the three, cellbender 0.3.2 and torch
  2.5.1 are on record for the environment that made the calibration's denoised objects; numpy 1.26.4 is
  not - it is the pin drafted together with torch - and `--check` says which is which.
  `tests/test_the_denoiser_environment_is_pinned.py`.
- **The core lock built an environment its own suite failed in.** It had no pyarrow, which
  `tests/test_native_io.py` needs, and it was a workstation's freeze: a notebook stack, a Windows-only
  package, and not a closure on linux. It is the cluster environment's 55 pins now, interpreter 3.12.12,
  and `requirements.hpc.lock.txt`, which no script read, is gone. The interpreter was first pinned as
  3.12.13, read from an older run's records written by another environment; the fresh-prefix build
  measured the reference environment's 3.12.12 (2026-11-09), and on 3.12.13 the analysis call records
  differed from the reference run's in `versions.python` alone.
  `tests/test_the_core_lock_carries_what_the_suite_needs.py`.
- **Locks were resolved again at every build.** `pip install -r` re-resolved anything a lock left
  open. Every lock is applied with `--no-deps` now and every installed version compared with it; a lock
  that is not exact is refused; an existing environment is checked and never modified; `--check` reads
  the locks and builds nothing.
- **A recipe's refusal did not fail the installer.** `install_rdoublet.sh` refusing an xgboost >= 2 was
  printed as a WARNING, and the final check, which only loads packages, reported the environment ok.
  Every component whose recipe refuses now fails the installer.
- **`conf/env/build_reference.sh` defaulted to one cohort's reference** (organism, assembly, release).
  `SPECIES`, `ASSEMBLY` and `ENSEMBL_RELEASE` are required, and the refusal names all three. Its one
  bash-4 expansion is gone, so it also runs under macOS's bash 3.2.
  `tests/test_build_reference_requires_its_reference.py`.

Open:

- **The pin does not make CellBender 0.3.2 able to write its checkpoint.** pyro attaches a weakref
  (`.unconstrained`) to every parameter it manages; recent torch pickles a parameter's `__dict__`, so
  `torch.save` raises `cannot pickle 'weakref'` after the last epoch and the posterior step finds no
  checkpoint. The calibration runs met this in the environment whose torch the lock pins, and saved
  through a wrapper outside this repository that detaches the attribute around the save.
  `adapters/cellbender.py` refuses the empty result; it does not work around the cause. MEASURED
  (2026-11-09, torch 2.5.1, CPU build): `torch.save` of a `Parameter` carrying a weakref attribute raises
  `TypeError: cannot pickle 'weakref.ReferenceType' object` - the message the failed cohort run died with -
  while the bare `Parameter` pickles. So a denoiser task run inside scQC is expected to fail this way under
  the pinned torch too; that last step is an INFERENCE, not measured inside scQC.
- **The conda half is not locked.** The R environment pins R, scDblFinder (and with it its
  Bioconductor release), xgboost and Matrix; the other CRAN packages and the BLAS are solved when it is
  built. The interpreters are pinned to the patch, their conda builds are not, and the micromamba the
  installer fetches is the newest.
- **Linux x86_64 only, for the denoiser.** Its lock names torch's CUDA 12.1 build, which exists for no
  other platform; on another, pip stops at torch.

