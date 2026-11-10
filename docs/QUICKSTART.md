# Quick start

From nothing to a filtered object and a report. Ten minutes of your attention; the run itself
takes as long as the cohort takes.

If you only want to know what the pipeline decided and why, skip to
[reading the report](#5-read-the-report).

---

## 1 · Install

```bash
git clone https://github.com/JiaenLin/scQC.git && cd scQC
pip install -e .
./bin/scqc selftest
```

`selftest` runs the bundled suites. It needs no cohort and no third-party package: the decision
layer is standard library only, deliberately, because a gate you cannot install is a gate that
gets skipped.

To *run* a cohort you also need the tools the pipeline invokes:

```bash
setup/install_env.sh --prefix ~/scqc-env --all
```

Four environments, because the aligner, the denoiser and the analysis stack have mutually
incompatible pins. Each is built from its lock in `conf/env/` with nothing resolved on the day, and
checked against it afterwards; `setup/install_env.sh --check --all` reads the locks and builds nothing.

## 2 · Make a project

```bash
setup/init_project.sh --dir ~/projects/my-study --assay snrna --samples 10
```

scQC never writes inside its own directory. One installation serves any number of datasets, and
upgrading the pipeline cannot disturb a result you already have.

## 3 · Fill in the samplesheet

`<project>/samplesheet.csv` (the CLI looks for that name; `--samplesheet` names another). Five columns are required; everything else is either a
path the pipeline needs or a **design factor it discovers on its own**.

```tsv
sample	platform	species	reference	assay	matrix	mt_prefix	ribo_pattern	condition
S1	singleron	mus_musculus	mus_musculus/ensembl_112	snrna	/data/S1/outs/raw	mt-	^Rp[sl]	control
S2	singleron	mus_musculus	mus_musculus/ensembl_112	snrna	/data/S2/outs/raw	mt-	^Rp[sl]	treated
```

| column | required | what it is |
|---|---|---|
| `sample` | **yes** | unique library name; every output is keyed by it |
| `platform` | **yes** | declared, no default. To rebuild a matrix from FASTQ it must be `10x` or `singleron` (the aligners this pipeline runs), and anything else refuses rather than guessing; a verified supplied `matrix` from any platform is accepted, recorded `processor: none (supplied matrix)` |
| `species` | **yes** | for the record and the reference resolution |
| `reference` | **yes** | the genome/annotation the counts were made against, as a `<species>/<build>` key of the project's registry |
| `assay` | **yes** | `snrna` (nuclei) or `scrna` (whole cells). It sets the mitochondrial bound and the ambient policy, and it has **no default**: until 2026-10-02 it was optional, never reached the mitochondrial step, and every cohort was bounded as single-nucleus |
| `matrix` | to start from counts | a MatrixMarket directory, a CellRanger `.h5`, or an `.h5ad` |
| `fastq_r1` / `fastq_r2` | to start from reads | alignment runs first and produces the matrix |
| `mt_prefix`, `ribo_pattern` | for steps 5 and 6 | **species-specific; nothing guesses them** |
| `ambient_h5` | if already denoised | a CellBender object you produced elsewhere |
| `ambient_skip_reason`, `ambient_skip_approved_by` | to skip CellBender on a **nucleus** library | both, or the run is refused: the reason and the person who approved it. The library then takes the not-denoised route (raw counts, `aligner_cells` required); the skip is DECLARED in the report and STATUS.json. Never needed on `scrna`, where CellBender is never run |
| `aligner_cells` | for step 2; **required for `scrna` with a raw matrix** | the aligner's filtered matrix (CellRanger `filtered_feature_bc_matrix`, CeleScope `outs/filtered`). On nuclei it is compared with the denoiser's call; on single cells, which this pipeline does not denoise, it **is** the cell call, and a raw-droplet run without it is refused before anything is submitted. A single-cell `matrix` that holds called cells only - the authors' filtered matrix, as most public data is published - needs none: it is accepted as a **cell-called input**, its own barcodes are the cell call, and the run records that it has no empty droplets (no ambient estimate, no empty-droplet model, no F1) |
| anything else | no | a **design factor**, discovered automatically |

`mt_prefix` and `ribo_pattern` have no defaults on purpose. Mouse and human differ in case alone
(`mt-` against `MT-`), and a wrong prefix gives every cell `pct_counts_mt == 0` — which is
indistinguishable from a clean library and passes every mitochondrial gate silently.

**Design factors are discovered, and may be declared.** Any extra column becomes a factor if it
varies, has at most six levels, and leaves at least one level holding more than one sample. That
last rule excludes identifiers: a replicate id in a four-sample cohort has four levels and would
otherwise make every differential check a ratio between single libraries — arithmetic with no
evidence in it, reported in the same words as a real result. Every discovered factor is checked;
with `--design-factors genotype[,...]` a differential REFUSES only on a declared factor, and one on
a column discovered but not declared (a lane, a chemistry) is a REVIEW that says so. With nothing
declared, any discovered factor can refuse.

**The project's reference registry** is `<project>/registry.tsv`, tab-separated, one row per
reference the samplesheet names: `species`, `build`, `path`, then `built_by`, `star_version` and
`notes` ([what a complete row records](TOOLS_AND_REFERENCES.md)). `init_project.sh` writes its
header. scQC ships none, and both commands below refuse without `--registry`.

Check the samplesheet before computing anything:

```bash
scqc validate --samplesheet ~/projects/my-study/samplesheet.csv \
              --registry ~/projects/my-study/registry.tsv
```

## 4 · Run

```bash
scqc run --project ~/projects/my-study --registry ~/projects/my-study/registry.tsv --jobs 8
```

Independent tasks run concurrently; `--jobs 1` is serial, which is what you want when a failure
has to be read in one log. The default mode is `apply`, which writes filtered objects. A run asked
for differently writes beside the old one rather than over it; the one exception is a changed
`decisions.yml` — see [run keys](USER_GUIDE.md#run-keys-and-resuming).

Outputs land under `results/<digest>/`:

```
objects/   cohort.deliverable.h5ad and one filtered .h5ad per library
tables/    every number the report quotes, as CSV
reports/   qc_report.html — open this
figures/   (figures are embedded in the HTML; this holds any written separately)
```

## 5 · Read the report

Open `results/<digest>/reports/qc_report.html`. It is one self-contained file: no server, no
network, every figure embedded.

Read it in this order — it is arranged as an argument, not as a log:

1. **The decision strip.** Observations in, kept, removed, and *evenness across the design*. That
   last number is the one a conventional QC report does not carry: a filter that falls harder on
   one arm of the design puts a technical gradient exactly where the biology is measured.
2. **Where the design is confounded.** Every pair of design factors, classified `aliased` /
   `nested` / `crossed` by exact comparison of the partitions they induce over your libraries.
   Where two factors are aliased, every QC metric is drawn across the arms that leaves, with each
   library's own median inside its arm — an arm difference no larger than the spread between
   libraries of the same arm is a library effect. Then the removal rate per arm, per criterion. A
   design with no aliased pair says so in one line.
3. **Where each count floor came from, and what quality control did.** The density each valley
   was read off; then one row per criterion — the distribution before the cut, what survived it,
   and the threshold with how it was arrived at (DERIVED from this data, DECLARED before seeing
   it, or ADJUDICATED by a person after reading the evidence; per library, or one cohort
   constant).
4. **Findings.** Every REFUSE and REVIEW, with the step that raised it. REVIEW does not stop a
   run; it means a person has to look.
5. **What this run could not establish.** Each step states its own limit. An omitted limit reads
   as no limit.

A figure that could not be produced says so, in its place, with the reason. That is not the same
as a figure nobody wanted, and the report does not let the two look alike.

## 6 · Iterate on the report without re-running

```bash
scqc report ~/projects/my-study/results/<digest>
```

Rebuilds the document from the run's own files in seconds. It cannot reach the matrices, so no
number in it can change.

---

## Where to go next

- **[User guide](USER_GUIDE.md)** — every command, resuming, approvals, and what to do when a gate
  refuses.
- **[How each filter is calculated](FILTERS.md)** — the exact procedure behind each threshold.
- **[Output reference](OUTPUTS.md)** — every file and every column.
