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
incompatible pins.

## 2 · Make a project

```bash
setup/init_project.sh --dir ~/projects/my-study --assay snrna --samples 10
```

scQC never writes inside its own directory. One installation serves any number of datasets, and
upgrading the pipeline cannot disturb a result you already have.

## 3 · Fill in the samplesheet

`~/projects/my-study/samplesheet.tsv`. Four columns are required; everything else is either a
path the pipeline needs or a **design factor it discovers on its own**.

```tsv
sample	platform	species	reference	assay	matrix	mt_prefix	ribo_pattern	condition
S1	singleron	mus_musculus	refs/ensembl_112	snrna	/data/S1/outs/raw	mt-	^Rp[sl]	control
S2	singleron	mus_musculus	refs/ensembl_112	snrna	/data/S2/outs/raw	mt-	^Rp[sl]	treated
```

| column | required | what it is |
|---|---|---|
| `sample` | **yes** | unique library name; every output is keyed by it |
| `platform` | **yes** | `10x` or `singleron`. Anything else refuses rather than guessing |
| `species` | **yes** | for the record and the reference resolution |
| `reference` | **yes** | the genome/annotation the counts were made against |
| `matrix` | to start from counts | a MatrixMarket directory, a CellRanger `.h5`, or an `.h5ad` |
| `fastq_r1` / `fastq_r2` | to start from reads | alignment runs first and produces the matrix |
| `mt_prefix`, `ribo_pattern` | for steps 5 and 6 | **species-specific; nothing guesses them** |
| `ambient_h5` | if already denoised | a CellBender object you produced elsewhere |
| `aligner_cells` | for step 2 | the aligner's filtered matrix, to compare cell calls against |
| anything else | no | a **design factor**, discovered automatically |

`mt_prefix` and `ribo_pattern` have no defaults on purpose. Mouse and human differ in case alone
(`mt-` against `MT-`), and a wrong prefix gives every cell `pct_counts_mt == 0` — which is
indistinguishable from a clean library and passes every mitochondrial gate silently.

**Design factors are discovered, never declared.** Any extra column becomes a factor if it varies,
has at most six levels, and leaves at least one level holding more than one sample. That last rule
excludes identifiers: a replicate id in a four-sample cohort has four levels and would otherwise
make every differential check a ratio between single libraries — arithmetic with no evidence in
it, reported in the same words as a real result.

Check it before computing anything:

```bash
scqc validate ~/projects/my-study/samplesheet.tsv
```

## 4 · Run

```bash
scqc run --project ~/projects/my-study --jobs 8
```

Independent tasks run concurrently; `--jobs 1` is serial, which is what you want when a failure
has to be read in one log. The default mode is `apply`, which writes filtered objects. Nothing is
ever overwritten — see [run keys](USER_GUIDE.md#run-keys-and-resuming).

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
   and the threshold with how it was arrived at (DERIVED from this data, or DECLARED before
   seeing it; per library, or one cohort constant).
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
