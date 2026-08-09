# Acceptance test — reproduce a cohort's recorded QC result

**The pipeline is not finished until this passes.** It is a regression test, not a validation: it
proves the pipeline does what the reference cohort's recorded result says was done. It does not
prove that result was right.

## Why two tiers

Most of the chain is deterministic. **Ambient correction is not** — it is variational inference,
and its cell call sets the denominator for every number after it. Demanding an exact match
end-to-end would make the test fail for a reason that has nothing to do with a regression, and a
test that fails for the wrong reason gets switched off.

So the deterministic part is tested exactly, and the stochastic part is tested against a tolerance
that is **measured rather than guessed**.

## Tier 1 — deterministic, must be EXACT

**Input:** the stored denoised per-library objects and the recorded per-cell table.
**Covers:** steps 2–7. Every one is deterministic given its input; scDblFinder is seeded and
verified reproducible (its vignette's `RNGseed` caveat applies only to the `samples` argument,
which is not used).

Ten checkpoints are compared. The exact integers belong to the cohort you are checking against
and live only in your untracked `expected.local.tsv`; the percentages below are the calibration
cohort's, and are here to say what each checkpoint *means* and roughly where it should land.

| checkpoint | what it counts | calibration cohort |
|---|---|---|
| `droplets_analysed` | every droplet in the per-cell table, before any call | the denominator |
| `ambient_cells` | droplets the denoiser called a cell | 80.2% of droplets |
| `scored_for_doublets` | called cells above the light floor, so a detector could score them | 75.8% of called cells |
| `doublets_called` | cells the detector called a doublet at the chosen `dbr.sd` | 12.2% of those scored |
| `step5_after_quality` | cells surviving the count and mitochondrial filters — the doublet criterion excluded, because it is step 4's | 61.1% of called cells |
| `deliverable` | what survives everything | 51.9% of called cells |

Four more compare the **unique contribution** of each criterion — cells that criterion removes and
no other one would. This is what separates a criterion that is doing work from one that is
duplicating a neighbour. In the calibration cohort, the doublet criterion removed 9.2% of called
cells on its own, while the UMI floor removed 1.2%, the gene floor 1.0% and the mitochondrial
ceiling 0.4%.

**Any deviation is a regression.** Nothing in tier 1 has a legitimate reason to move.

**A checkpoint with a blank or missing expectation is a FAILURE, not a skip.** A check with
nothing to compare against reports success under every regression it exists to catch, so the
runner counts it as failed and says how many checkpoints were actually compared. The reverse also
fails: an expectation whose name no checkpoint reads is reported, because a row nobody reads is
how a file comes to describe coverage that does not exist.

## Tier 2 — from raw FASTQ, tolerance MEASURED first

**Input:** the FASTQ pairs. **Covers:** steps 0–7, including alignment and ambient correction.

Before this tier can have a pass condition, the tolerance has to be established:

1. Run the denoiser **twice on one library** with identical input, seed and environment.
2. Record how much the cell call moves.
3. The tolerance is that spread, stated — not a round number chosen afterwards.

Until step 3 is done, tier 2 has **no pass condition** and must not be claimed as passing. An
untested tolerance is worse than no tolerance: it reads like a criterion.

Two things already known to constrain it:

- **The aligner should be exactly reproducible** given the same reference and parameters, so any
  tier-2 drift is attributable to the denoiser rather than to alignment. That is worth verifying
  separately, because it makes the tolerance interpretable.
- **The reference must be the same index.** `references/_registry/registry.tsv` records the
  aligner and STAR version an index was built with, and whether it retains introns. A rebuilt
  environment could arrive with a different STAR; if it does, tier 2 is testing two things at once.

## What neither tier tests

- **That the thresholds are right.** The count floors, the mitochondrial ceiling and the doublet
  rule are adjudicated decisions recorded in the cohort's `decisions.yml`. Reproducing them proves
  the pipeline applies them faithfully, not that they were the right thresholds.
- **That the deliverable is biologically correct.** No cell type is identified anywhere in steps
  0–7, so the cost of the cut in lost populations is unknown — the calibration cohort records this
  itself.
- **Portability.** The calibration cohort is one tissue, one platform, one operator. Passing this
  says nothing about another platform's data. [`CALIBRATION.md`](../../CALIBRATION.md) classifies
  every measured parameter by whether it transfers, and
  [`docs/WORKFLOW.md`](../../docs/WORKFLOW.md) shows where each class enters the run; most of what
  the calibration cohort measured does not transfer as a number.
- **Threshold DERIVATION**, unless your cohort recorded its valleys. Tier 1 checks that the
  thresholds recorded were the thresholds applied. It re-derives them only if `valley_table`
  exists (see below), and otherwise says so rather than implying coverage it does not have.

## Files

```
tests/acceptance/
  run_tier1.py           the runner
  expected.tsv.template  the checkpoint names, values blank — copy and fill
  expected.local.tsv     your filled copy. GITIGNORED; never committed
  schema.local.tsv       optional table/column overrides. GITIGNORED
  README.md              this file
```

Inputs are not copied here — they are the reference index plus the cohort's own objects. Point the
runner at a cohort you hold:

```bash
COHORT_DIR=/path/to/cohort python tests/acceptance/run_tier1.py
```

With `COHORT_DIR` unset or pointing nowhere, and when pandas is absent, the runner prints a `SKIP:`
line and exits 0, so the suite stays runnable on a machine that hosts no cohort. Once `COHORT_DIR`
*is* set the check has been asked for, and a missing `expected.local.tsv` is a failure rather than
a skip — otherwise the absence of expectations would look exactly like a pass.

## Adapting it to your own tables

`run_tier1.py` reads a specific table layout. Its defaults are the calibration cohort's, and every
one of them is overridable from `schema.local.tsv` — one `key<TAB>value` per line, `#` comments
allowed, unknown keys rejected rather than ignored.

| key | default | what it names |
|---|---|---|
| `tables_dir` | `results/tables` | recorded tables, relative to `COHORT_DIR` |
| `objects_dir` | `results/objects/cellbender_h5ad` | per-library denoised objects |
| `denoised_glob` | `*_lr5e-5.h5ad` | which objects in that directory to read |
| `denoised_sample_sep` | `_cellbender` | filename separator; the library name is what precedes it |
| `per_cell_table` | `qc_filter_per_cell.csv.gz` | one row per droplet |
| `doublet_sweep_table` | `doublet_dbrsd_sweep_per_cell.csv.gz` | one row per scored cell |
| `cluster_profile_table` | `cluster_check_profile.csv.gz` | step 6's per-cluster profile |
| `valley_table` | `quality_valleys.csv` | per-library density valleys: `sample`, `metric`, `valley`, `bimodal` |
| `col_sample` | `sample` | library identifier column |
| `col_aligner_cell` | `celescope_cell` | the aligner's cell call, in the denoised object |
| `col_ambient_cell` | `cellbender_cell` | the denoiser's cell call |
| `col_doublet_score` | `doublet_score` | null where a cell was never scored |
| `col_total_counts` | `total_counts` | UMI per droplet |
| `col_keep` | `keep` | the recorded deliverable mask |
| `col_doublet_call` | `scdbl_sd0.06_call` | the doublet call in the sweep table |
| `fail_columns` | four `fail_*` columns | one boolean column per removal criterion |
| `doublet_fail_column` | `fail_scdblfinder_sd006` | which of those belongs to step 4, so step 5's checkpoint can exclude it |
| `design_factor` | `condition` | the name the design checks report under |
| `design_levels` | `treat=treat,ctrl=ctrl` | `level=substring`, first match wins; a library matching none is `unassigned` rather than silently joining an arm |
| `light_floor` | `200` | the technical floor for doublet scoring |
| `quality_floor` | `350` | the lowest count floor, used only to check the two do not collide |
| `cluster_resolution` | `1.0` | which resolution to read from the cluster profile |
| `cluster_algorithm` | `leiden` | which algorithm to read from it |

If you override `fail_columns`, replace the four `only_*` rows in your `expected.local.tsv` to
match: each is named `only_` followed by the column name with its `fail_` prefix removed. A name
that no longer matches produces a missing expectation, which fails — deliberately, and loudly.
