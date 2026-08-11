# Output reference

Every file a run writes, what is in it, and which files are meant to be read next.

Read [FILTERS.md](FILTERS.md) for how the numbers in these files are calculated.

---

## Where a run writes

Outputs are stored under a name derived from what produced them:

```
<project>/
├── results/
│   ├── INDEX.tsv                 one line per run: digest, mode, samples, parameters, first seen
│   ├── latest ->                 a pointer to the newest run
│   └── <digest>/
│       ├── INPUTS.json           what this directory was claimed for
│       ├── objects/              matrices
│       ├── tables/               CSV, small enough to open
│       ├── reports/              the rendered report and its JSON companion
│       └── figures/              (nothing writes one yet)
├── work/<digest>/                run manifest and intermediates
└── logs/<digest>/                one log per task
```

The digest is computed from the samplesheet's **content**, the declared parameters and the mode.
The same inputs go to the same directory, which is what lets a re-run reuse completed work. Change
a threshold and the digest changes with it, so the new run lands beside the old one rather than
over it.

`INDEX.tsv` is the file to read first — a directory of digests answers *where is it* and not
*which one do I want*.

> **`latest` is a convenience and nothing should depend on it.** It is rewritten by every run and
> is the only thing here that does not accumulate. Scripts should resolve a digest from
> `INDEX.tsv` and use it.

---

## objects/

| file | written by | contents |
|---|---|---|
| `<sample>_ambient.h5` | step 1 | the denoised matrix for one library. Every later step reads this. Unfiltered. |
| `cohort_per_sample/<sample>.filtered.h5ad` | step 7 | one library, filtered. **The primary deliverable.** |
| `cohort.deliverable.h5ad` | step 7 | all libraries, filtered, merged |

Step 7 runs in apply mode only. In evidence mode neither filtered object exists, and the task that
would write them is not in the task graph at all.

### The merged object is the merge

`cohort.deliverable.h5ad` is built by reading the per-library files back from disk and
concatenating them — not from copies held in memory while they were written. The two would
ordinarily be identical, which is the point: where they are not, the difference belongs in the
merged object rather than being hidden by a good copy in memory. What a reader can open is what
was merged.

The merge is checked against the sum of its parts. A concatenation that silently dropped an
observation would otherwise produce a smaller cohort that still looks complete.

**The per-library objects are primary and the merged one is derived.** Each library is filtered on
its own mitochondrial ceiling and cluster-checked on its own clustering; the cohort object is
those results put together. Use the per-library objects for annotation, and the merged one for
integration and differential testing.

### obs columns on a filtered object

| column | type | meaning |
|---|---|---|
| `sample` | categorical | the library this nucleus came from |
| `total_counts` | float | UMI in this nucleus — **the value the UMI floor was applied to** |
| `n_genes` | float | genes detected — the value the gene floor was applied to |
| `pct_counts_mt` | float | mitochondrial percentage — the value the ceiling was applied to |
| `doublet_score` | float / absent | the detector's score, where the call file carried one |
| `doublet_class` | label / absent | `singlet` / `doublet`, or absent where never scored |
| `cluster` | label | its cluster in that library's clustering. A **label**, not a number — cluster ids read as integers and are not |
| `cluster_FLAG` | `True` / `False` / absent | that cluster's flag verdict from step 6 |
| `cluster_WATCH` | `True` / `False` / absent | that cluster's watch verdict |
| `cluster_pct_doublet` | float / absent | the cluster's doublet percentage |
| `cluster_median_pct_mt` | float / absent | the cluster's median mitochondrial percentage |

Whatever the input object already carried is preserved alongside these.

**The first five are the values the filter read**, carried onto the object so a reader can see why
any nucleus survived without recomputing anything. They come from the per-cell table rather than
being recomputed at write time, so the object provably holds the numbers the criteria were
evaluated on — recomputing would give the same answers and could not prove it.

> **Absent is a third state and is preserved as one.** A barcode step 6 never labelled, or a
> cluster missing from the profile, leaves these `None` — never `False`, never `0`. *"This cluster
> was not flagged"* and *"this barcode was never examined"* are different facts, and a deliverable
> that conflates them claims a check that did not happen. Test with `.isna()`, not `== False`.

The flags are carried rather than recomputed because **re-clustering a filtered object does not
reproduce them**: criterion D is a doublet fraction, and once the doublets have been removed every
cluster is 0% doublet by construction.

### uns

`uns["scqc_apply"]` on the merged object records `library`, `library_n_in`, `library_n_kept`,
`n_delivered` and `built_from` — three parallel arrays rather than a list of records, because HDF5
has no record type and h5py fails on one.

---

## tables/

All CSV, all readable with the standard library alone.

### Per library

| file | one row per | key columns |
|---|---|---|
| `<sample>_doublets.csv` | scored barcode | `barcode`, `doublet_score`, `doublet_class` |
| `cell_calls.csv` | library | `aligner`, `denoiser`, `lost` |
| `valleys_umi.csv`, `valleys_genes.csv` | library | `valley`, `bimodal`, `note` |
| `<sample>.barcode_rank.csv` | **rank point** | the raw matrix's barcode-rank curve, downsampled log-uniformly: `rank`, `total_counts`, `n_barcodes`. Figure F1. Absent for a library rebuilt from FASTQ, which has no supplied matrix to measure |
| `<sample>.valley_density.csv` | grid point × metric | the KDE the valley was read off: `grid`, `grid_log10`, `density`, `is_valley`, `is_mode`. Figures F6 and F13 |
| `<sample>.embedding.csv` | **barcode** | the 2-D coordinates of every barcode the denoiser called, with `clustered` and `cluster` joined on: `x`, `y`. Figures F10 and F11 |
| `<sample>.light_floor.csv` | library | that library's export: how many barcodes cleared the floor, how many sat below it, how many were not selected |
| `<sample>.doublet_sweep.csv` | setting | the called rate at each swept `dbr.sd` for that library. `--dbr-sd-sweep` only |
| `mito_ceiling_per_sample.csv` | library | quartiles, `derived`, `ceiling`, `clamped`, and **the population it was derived over** |
| `<sample>.percell.csv` | **barcode** | every barcode the library held: the four measured values, the doublet score and class, and which of the five criteria fired. Apply mode only |
| `ambient_lr_diagnostics.csv` | library | `fraction_removed`, `convergence_indicator`, `measured` |
| `thresholds_per_sample.csv` | library | every threshold the run derived, each column marked *per library* or *cohort constant* |

`thresholds_per_sample.csv` carries a second header row giving each column's scope. It is the
quickest answer to *which numbers differ because the libraries differ*.

Beside `<sample>.percell.csv` you will find a small `<sample>.apply_measure.metrics.json`. It is
not a result: it records the size and modification time of each output **as the process that wrote
them finished writing**, which is how the caller tells a file this run produced from one that was
already there under the same name. Left in place because a proof that is deleted proves
nothing.

### Cohort

| file | one row per | notes |
|---|---|---|
| `cluster_profile.csv` | (sample, cluster) | the profile and the A/B/C/D/FLAG/WATCH verdicts. Not joinable to barcodes — use the `cluster` column on the object |
| `removal_ledger.csv` | **removed** barcode | every criterion that fired on it, not just the first |
| `ambient_summary.csv` | library | fraction removed, genes fully removed |
| `ambient_supplied.json` | — | provenance of any denoised object supplied rather than produced here |
| `doublet_sweep.csv` | (library, setting) | the called rate at each swept `dbr.sd`, assembled from the per-library tables. Figure F5. Written **only** when `--dbr-sd-sweep` was given; it applies nothing and changes no deliverable |
| `doublet_health.csv` | library | the rates the step-4 gate judged: `n_scored`, `n_called_doublet`, `rate_over_scored`. A gate whose evidence is not on disk cannot be re-checked |
| `light_floor.csv` | library | what the light floor left out: `n_exported`, `n_below_floor`, `n_not_selected`. The two reasons a barcode was never scored are kept apart — only the first is explained by a threshold |

> **The embedding is over the CELL-CALLED population, not the clustered one.** Step 6 clusters
> the cells that reach the deliverable — otherwise its flags describe empty droplets — and embeds
> every barcode the denoiser called, so the nuclei the count floors and the ceiling removed are
> still in the picture. That is the only way F11 can answer whether a removal took a coherent
> region. `clustered` is `False` for a barcode that was embedded and not clustered; it is never
> blank, because "embedded and not clustered" is a fact about that barcode and a blank would read
> as "not recorded".

> **The ledger lists what LEFT, not what stayed.** It is the record that makes a removal
> recoverable: re-read the input object with those barcodes and the removed population is back.
> It stores identifiers, not counts.
>
> `<sample>.percell.csv` is its complement and covers **every** barcode, kept and removed alike,
> with the value each criterion was evaluated on. It is the file to re-check a filter with, or to
> ask how close a retained nucleus came to a threshold — which the ledger cannot answer, because
> the ledger does not contain it.

---

## reports/

| file | |
|---|---|
| `qc_report.html` | the rendered report, self-contained — no request leaves the file |
| `report.json` | the same content as data. Every number in the HTML is in here beside it |

The report is organised by step. Each step block gives what the step does, **what it cannot
establish**, the numbers it produced and the file each came from.

**The report audits itself.** Anything the run should have recorded and did not is counted as a
defect on the front page and listed in the JSON, because a report that silently omits a section
reads exactly like a complete one. A defect count above zero is normal and is not a failure of the
run; it is the report saying what it was not given.

`report.json` also carries the run's wall-clock, the summed task time and the ratio between them.

---

## work/ and logs/

`work/<digest>/state.json` is the run manifest: one record per task, with its status, signature,
declared outputs, metrics and log path. It is what makes a re-run reuse completed work — delete it
to force everything to run again.

`work/<digest>/` also holds intermediates that no downstream step opens — the keep-lists and
annotation tables step 7 hands to the writer. The per-barcode criteria tables are **not** here:
they are a deliverable and live in `tables/`. `logs/<digest>/` holds one log per task, named for
the task.

---

## What is not produced

- **No figures.** `report/figures.py` exists and the report expects F1–F9; no step supplies one.
  Each absence is reported as a defect rather than omitted.
- **No freshness check.** The report can compare its own timestamp against its inputs, and no step
  supplies a newest-input time, so every report says `NOT CHECKED` rather than claiming to be
  current.
- **No modified inputs.** Nothing a run writes replaces anything it read. The denoised objects in
  `objects/` are the run's own copies; the matrices named in the samplesheet are never written to.
