# scQC

**Quality control for single-cell and single-nuclei RNA-seq that separates deriving a threshold
from applying it.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-0.4.0-blue.svg)](#status)

scQC measures what it can from the data, declines to guess what it cannot, and records who set
every threshold it applies — the data or a named person — so the two cannot be confused later.

📖 **[Quickstart](docs/QUICKSTART.md)** · **[User guide](docs/USER_GUIDE.md)** ·
**[Known issues](KNOWN_ISSUES.md)**

---

## Install

```bash
git clone https://github.com/JiaenLin/scQC.git && cd scQC
pip install -e .              # the decision layer: stdlib only, no dependencies
./bin/scqc selftest
```

The decision layer has no third-party dependencies, so a gate cannot be skipped for want of an
install. The tools a full run invokes are separate:

```bash
setup/install_env.sh --prefix ~/scqc-env --all
```

Four environments, because the aligner, the denoiser and the analysis stack have mutually
incompatible pins. Select with `--with-celescope`, `--with-cellbender`, `--with-doublet` or
`--all`. Requires conda, mamba or micromamba. Each Python environment is built from a lock in
`conf/env/` — every package at one exact version, the interpreter's patch included — with
`pip install --no-deps`, and every installed version is compared with the lock afterwards. The
denoiser's lock must pin torch 2.5.1 and numpy 1.26.4; any other is refused before anything is built.
`--check` reads the locks and builds nothing. An environment an interrupted build left behind -
`conda-meta/` and no interpreter - is refused; `--repair` removes it and builds it again.

## Run

```bash
setup/init_project.sh --dir ~/projects/my-study --assay snrna --samples 10
# fill in samplesheet.csv, and registry.tsv with the references it names

R=~/projects/my-study/registry.tsv
scqc run --project ~/projects/my-study --registry $R                  # apply, the default: derive, apply, deliver
scqc run --project ~/projects/my-study --registry $R --mode evidence  # optional: measure everything, remove nothing
# optional: override a threshold in decisions.yml (value, approved_by, verbatim) and run apply again
```

The reference registry is the project's: `--registry` is required, and scQC ships none.

scQC never writes inside its own directory. One installation serves any number of projects, and
upgrading it cannot disturb an existing result. Runs execute locally or on PBS with one job per
task.

## Parameter classes

Every parameter carries a class, and the class determines who may set it.

| class | meaning | set by |
|---|---|---|
| **FIXED** | true regardless of dataset | the pipeline; changing it is a code change |
| **DERIVED** | procedure fixed, value computed per dataset | the data — never a command-line argument |
| **DECLARED** | platform, species, design | you, in advance; no default exists |
| **ADJUDICATED** | not derivable; needs judgement | you, after reading the evidence, in your own words |

Adjudicated parameters have no defaults. A value recorded without an approver and their reasoning
is not treated as an approval, and the derived value is used instead.

## Two phases

**Apply mode**, the default, takes every measurement, applies the DERIVED and DECLARED
thresholds, records who decided each one, and performs the removal - at step 7, the only step that
removes. **Evidence mode** takes the same measurements, writes the report, and removes nothing.

`decisions.yml` is optional. An ADJUDICATED override of a threshold needs a `value`, an
`approved_by` and the approver's `verbatim` words; a value missing either attestation is not an
override, and the derived value is applied and reported as derived. Where a decisions file gives
approval words for the removal, they are matched against the thresholds they were given for: move a
threshold afterwards and they no longer apply.

## Steps

| # | step | does | removes |
|---|---|---|---|
| 0 | ingest | validate samplesheet, resolve reference, verify the input: raw values always, raw droplets where a denoiser needs them (nuclei) | — |
| 1 | ambient | denoise (CellBender) or accept a denoised object; audit the removal | — |
| 2 | cell call | compare aligner and denoiser calls; gate the loss | — |
| 3 | light floor | technical floor for doublet scoring — not a quality filter | — |
| 4 | doublets | score per sample, before quality filtering; flag only | — |
| 5 | quality | derive count floors and the per-library mitochondrial ceiling | — |
| 6 | cluster check | per-cluster flags: depth, mitochondrial, markers, doublet fraction | — |
| 7 | **apply** | measure, write the ledger, audit the removal (a failed check refuses), then write the filtered objects | **yes — only here** |

Single nuclei and single cells take different routes, decided by the declared `assay`.
Nuclei are denoised by CellBender (mandatory, unless a skip is declared in the samplesheet with
`ambient_skip_reason` and `ambient_skip_approved_by`, which is recorded and reported), and the
denoiser's call is the cell call. **CellBender is never run on single cells** (the PI, 2026-10-02);
no other ambient method is built for cells here yet, and one may be added per assay. On single
cells the raw counts are the object, the aligner's call (`aligner_cells`) is the cell call, and the
count bounds and the cluster check are the whole-cell ones. Every line calibrated on nuclei - the
clustering resolution, the doublet silence line, the cell-call gate lines, the design-differential
line - is a per-assay table with its provenance; the single-cell entries carry the same values and
say they were calibrated on nuclei. Apply mode runs on single cells too: step 7's cell criterion is
`fail_not_aligner_cell`, over the aligner's called cells. A single-cell library may also arrive
as called cells only - the authors' filtered matrix, as most public single-cell data is published:
it is accepted as a **cell-called input**, its own barcodes are the cell call, and the run records,
per library and in STATUS.json, that it holds no empty droplets - so no ambient estimate, no
empty-droplet model and no barcode-rank figure. On nuclei a cell-called matrix is still rejected:
the denoiser needs the empties. A library supplied already denoised (`ambient_h5`, with its
provenance) is used as supplied on either assay, and its cell call is the denoiser's.
Doublet scoring precedes quality filtering, as scDblFinder's documentation requires.

📄 **[Why each step exists](docs/RATIONALE.md)** · **[How each filter is
calculated](docs/FILTERS.md)** · **[Workflow diagrams](docs/WORKFLOW.md)**

## The confounding block

The report opens with what the design can and cannot support. Three figures, computed from the
samplesheet and the per-cell tables:

**F16** compares every pair of design factors and classifies the pair as `aliased` (each factor is
fixed once the other is known), `nested` (the implication runs one way) or `crossed` (separable).
The comparison is exact set arithmetic over the partitions each factor induces on the libraries —
there is no statistic and no p-value, because a confounded factor cannot be tested. The samplesheet
is redrawn beside the matrix so the classification can be checked by eye.

**F17** draws every QC metric across the arms the aliasing leaves: `total_counts`, `n_genes`,
complexity, `pct_counts_mt`, `pct_counts_ribo`, `nuclear_fraction` and `doublet_score`. A violin per
arm over the nuclei the filter kept, with **each library's own median drawn inside its arm**. Each
panel states the gap between the arm medians against the widest spread of library medians within a
single arm. A gap no larger than that spread is a library effect, not an arm effect.

**F18** gives the removal rate per arm, in total and per criterion, with the ratio between the
extreme arms printed for each. This is the part of the confounding the pipeline itself introduces.

Where no two factors are aliased, F17 and F18 say so instead of constructing a contrast the design
does not support.

📄 **[Report design](docs/REPORT_DESIGN.md)** — every figure, and what it is for.

## Output

Apply mode writes one filtered object per library plus a cohort object merged from those files
read back from disk. Every retained observation - nucleus or cell - carries `sample`, `cluster`,
`cluster_FLAG`, `cluster_WATCH` and the continuous values behind them, so nothing downstream needs
to re-cluster.

Every object carries `uns["scqc"]`: what the flag means, how many observations hold it, and a digest of
the exact set — so a consumer can act on the flag without inferring its meaning from a column
name, and can verify the column is still the one scQC wrote. `scqc stamp <objects>` adds it to
objects written before 0.4.0.

Per-library objects are primary; the cohort object is derived from them, since each library is
filtered on its own mitochondrial ceiling and cluster-checked on its own clustering.

📄 **[Output reference](docs/OUTPUTS.md)** — every file, every column.

### What a run leaves behind

Four states are distinguishable from the filesystem alone, by a reader who did not watch the run:

| file | written | says |
|---|---|---|
| `STATUS.json` | first as `partial`, last with the outcome | `ok` · `partial` (it died) · `refused` with a `fix` · `failed`; the commit, the job, the products actually on disk, `sees`, `cannot_show`, and every adjudicated threshold as an ask/decision pair |
| `RUNNING.txt` | at start | replaced at exit by one of the two below |
| `SEALED.txt` | exit 0 **and** every expected product present | a green light that checks the products, not just the exit status |
| `FAILED.txt` | anything else | names what is missing |

`scqc describe` prints the tool's declaration — needs, provides (one mask per criterion), sees,
gates and their verdicts, escapes, `cannot_show`, `state_version` — as JSON, for a host or a
reader that will not read the code.

## Reproducibility

All enforced today:

- **Nothing is overwritten by a run asked for differently - except through the decisions file.**
  Outputs land under `results/<digest>/`, the digest computed from the samplesheet content, the
  declared parameters, the mode and a digest of the pipeline's output-bearing code - a commit that
  changes only docs or tests keeps it, and the commit itself is recorded beside it. A directory
  records the inputs it was claimed for and refuses a run described differently. `decisions.yml` is NOT in the digest: an
  apply run with a changed decisions file runs in the same directory, re-runs the tasks of steps 5
  to 7 that carry it, and replaces their outputs.
- **Gates return a verdict, not a warning.** Each reduces to `REFUSE`, `REVIEW` or `PASS`, and the
  subcommands exit 2 on refusal. Reporting and enforcing are separate calls. A REFUSE can be
  accepted in `decisions.yml` (`accept: [{step, check, approved_by, verbatim}]`), by a named person
  in their own words; it is then recorded as an accepted REVIEW and in STATUS.json's `escapes`.
  A check that compares the run's own records with each other cannot fire on a correct run and is
  never accepted: step 7's ledger count, the count of the object written, and the removal audit's
  integrity checks (docs/WORKFLOW.md, "How a refusal is delivered").
- **Every removal is recoverable with the criterion that made it.** Removed observations are paired
  with the criteria that fired; the removal refuses if that record and the mask disagree on the
  count. Persisted as CSV, readable with the standard library alone.
- **Every threshold records who set it** — `DERIVED`, `DECLARED`, or `ADJUDICATED` with an
  approver and their own words — in the report's parameter table (`reports/report.json`); each
  adjudicated one is also in `STATUS.json`, and step 7 records who authorised the removal in the
  run's state. The removal ledger and the delivered objects do not carry it.
- **Reports carry the commit and tool versions**, obtained by asking each tool rather than reading
  a lockfile, and audit themselves: anything the run should have recorded and did not is counted
  as a defect on the front page.
- **A report older than its inputs is flagged, not refused.** Freshness is computed, printed on
  the report's face, counted as a defect when stale, and recorded as True, False or **None** —
  never False by default, because unchecked is not a pass. Nothing in scQC publishes a report, so
  nothing refuses one.

📄 **[Principles](docs/PRINCIPLES.md)** — the four enforced rules.

## Tests

```bash
scqc selftest     # every bundled suite including the adversarial one; non-zero exit if any failed
```

Each suite reports PASS, FAIL or SKIP. A SKIP is reported as its own outcome, not as a pass.
The decision-layer suites run on the standard library alone. The suites that drive the CLI, read
matrices or draw figures need the analysis stack (`conf/env/requirements.lock.txt`, or `core`
from `setup/install_env.sh`) and report SKIP without it; the cohort audit also needs a cohort
directory (`COHORT_DIR`).

Run `scqc` with `core`'s interpreter, not the system python — several steps read matrices in the
orchestrator's own process:

```bash
~/scqc-env/core/bin/python bin/scqc run --project my_project --mode evidence \
    --registry my_project/registry.tsv --python ~/scqc-env/core/bin/python
```

No data ships with this repository. The harness's development suite runs the pipeline end to end
on its synthetic two-shape fixture through `scqc fixture-run` (`DEVPOINTS.yaml`), which is told the
assay, the design's factor and the gene patterns and fails when the records the steps write show
one did not reach the step that reads it (each gate records the assay it ran under, the factors it
was given, and the factors of the design it computed on with their levels; a declared factor the
design leaves out is named with the reason, in the report's notices too; the report's records are
checked against the run's work state). A limit of the fixture, not of the tool: both of its shapes are cells, so it tells
only `scrna`, a step wired to `scrna` cannot be told apart from one that read it, and the nucleus
value never reaches the pipeline through this check; a real cohort is tested by reproducing its promoted run, not by a
runner shipped here.

## Part of a family

Four tools share one contract. Each works on its own, and they fit together:

| tool | does |
|---|---|
| [scQC](https://github.com/JiaenLin/scQC) | quality control, and who set each threshold |
| [scAnno](https://github.com/JiaenLin/scAnno) | cell-type labels, only as deep as the evidence goes |
| [scIntegrate](https://github.com/JiaenLin/scIntegrate) | whether you need batch integration, and which method |
| [scProfile](https://github.com/JiaenLin/scProfile) | communication, velocity, pseudotime, differential expression |

[single-cell-harness](https://github.com/JiaenLin/single-cell-harness) runs them as one stack, so a
decision made in one can be taken back later and every tool downstream is told what that costs.

## Changing this tool

Using scQC and changing it are different jobs, and the second one has its own suite.
`DEVPOINTS.yaml` in the root says what can be added here — criteria, steps and adapters — and what each one has to
declare. The suite lives in the harness:

```bash
git clone https://github.com/JiaenLin/single-cell-harness.git ../single-cell-harness
alias sch='PYTHONPATH=../single-cell-harness python3 -m sch'

sch dev map --root .                            # what plugs in here, and what it must declare
sch dev new criterion my_criterion       # a skeleton, a spec written before the code, and a test
sch dev check --point criterion --name my_criterion   # checks, cheapest first
```

Every check prints what it does **not** prove, so a green run is not read as more than it is. Two
of them run your addition against the same synthetic cohort written twice, with every column
renamed — code that asks for a role passes both, code that knows a column name passes one.

## Documentation

| document | answers |
|---|---|
| [QUICKSTART](docs/QUICKSTART.md) | install, samplesheet, run, reading the report |
| [USER_GUIDE](docs/USER_GUIDE.md) | every command, run keys, resuming, approvals, refusals |
| [FILTERS](docs/FILTERS.md) | how each filter is calculated and what it cannot establish |
| [OUTPUTS](docs/OUTPUTS.md) | every file a run writes, and the `obs` schema |
| [RATIONALE](docs/RATIONALE.md) | why each step exists and what goes wrong without it |
| [PRINCIPLES](docs/PRINCIPLES.md) | the four enforced rules |
| [WORKFLOW](docs/WORKFLOW.md) | diagrams of the pipeline, phases and parameter classes |
| [REPORT_DESIGN](docs/REPORT_DESIGN.md) | report layout and its figures |
| [TOOLS_AND_REFERENCES](docs/TOOLS_AND_REFERENCES.md) | tools, versions, reference resolution |
| [CALIBRATION](CALIBRATION.md) | what was measured, how much it varied, what one cohort cannot show |
| [KNOWN_ISSUES](KNOWN_ISSUES.md) | measured, reproduced, not yet fixed |
| [CONTRIBUTING](CONTRIBUTING.md) | what a pull request must answer before it can remove anything |

## Status

**0.4.0.** `scqc run` builds the task graph and executes it, locally or on PBS, invoking the
aligner, denoiser, doublet caller and analysis stack out of process. It writes the report and, in
apply mode, one filtered object per library plus a merged cohort object with a ledger naming every
barcode removed and why. All eighteen report figures are drawn. Per-step subcommands remain, for
judging tables produced elsewhere.

Read [KNOWN_ISSUES.md](KNOWN_ISSUES.md) before quoting a number from a run.

## Citing

If scQC contributes to published work, cite the repository and the version you ran. Every report
carries the commit.

## License

MIT — see [LICENSE](LICENSE).
