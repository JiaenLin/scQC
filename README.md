# scQC

**A quality-control pipeline for single-cell and single-nuclei RNA-seq that separates deriving a
threshold from applying it.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-0.1.0-blue.svg)](#status)

Most QC pipelines take thresholds as arguments. scQC treats them as findings: it measures what it
can from your data, refuses to guess what it cannot, and records who set every threshold it
applies — the data, or a person in their own words — so the two can never be confused later.

> **Read [Status](#status) before you plan a run.** At `0.1.0` scQC runs a cohort end to end —
> `scqc run` builds the task graph and executes it, locally or on a PBS scheduler with one job per
> task, invoking the aligner, the denoiser, the doublet caller and the analysis stack out of
> process. It writes a report and, in apply mode, a filtered object with its removal ledger.
> What it does **not** yet produce is any figure, and nothing feeds its freshness check. The
> per-step subcommands remain, for judging tables you produced elsewhere.

---

## Why

Three failures motivated this pipeline. Each is quiet, each survives review, and none is detected
by a run that completes successfully.

- **A delivered matrix that was already filtered.** Pre-filtered input cannot be un-filtered, and
  nothing downstream can tell. scQC verifies that its input is raw and refuses otherwise.
- **A doublet rate set by its prior rather than by the data.** The detector's default expected rate
  dominated the call. scQC sweeps that parameter and reports whether the rate is a measurement.
- **A threshold that transferred badly.** Across ten libraries of a single cohort — one tissue, one
  platform, one operator — the measured density valley ranged **274–473 UMI**. A hard-coded
  `--min_umi 500` sits above the valley in all ten. *The procedure transfers; the number does not.*

## The core idea

**Every parameter carries a class, and the class determines who is allowed to set it.**

| class | meaning | who sets it |
|---|---|---|
| **FIXED** | true regardless of dataset | the pipeline; changing it is a code change |
| **DERIVED** | procedure fixed, value computed per dataset | the data — never a command-line argument |
| **DECLARED** | platform, species, design | you, in advance; no default exists |
| **ADJUDICATED** | not derivable; needs judgement | you, *after* reading the evidence, in your own words |

A pipeline whose adjudicated parameters have defaults is reproducible only in the sense that it
does the same wrong thing every time.

## Two phases

The pipeline is built around two modes, and the separation between them is the design rather than
a workflow convenience:

- **evidence** — run every step, derive every DERIVED parameter, write the report, remove nothing,
  then stop. The apply task is not placed in the task graph at all, so there is no code path from
  this mode to a deletion.
- **apply** — the default. Run every step, then measure each removal criterion per barcode, write
  the removal ledger, and materialise the filtered object.

```bash
scqc run --project ./my-study --mode evidence   # measures everything, removes nothing
scqc run --project ./my-study                   # apply is the default
```

Thresholds are **DERIVED** by default: the pipeline applies what it measured and records it as
measured. To override one, put it in `decisions.yml` with an approver and your own words, and it
is recorded as **ADJUDICATED** instead. Both classes travel into the ledger, the written object
and the report, so a value the pipeline proposed is never later read as one a person chose.

```bash
cp decisions.template.yml decisions.yml         # optional; overrides what was derived
```

**Outputs are named after their inputs.** Each run writes under `results/<digest>/`, where the
digest is computed from the samplesheet's content, the declared parameters and the mode. The same
inputs go to the same directory, which is what lets a re-run reuse completed work; change a
threshold and the digest changes with it, so the new run lands **beside** the old one rather than
over it. `results/INDEX.tsv` says what each digest was and `results/latest` points at the newest.

Nothing is overwritten by a run that would have produced something different — that is a property
of the layout, not a rule anyone has to remember.

You can also judge one step at a time, from tables you produce yourself:

```bash
scqc validate --samplesheet samplesheet.csv       # every DECLARED field present, reference known
scqc verify --name mat --barcodes 12500 --genes 25000 \
            --min-counts 486 --max-counts 71204 --p98-counts 24900   # is this matrix raw?
scqc gate-cells --calls calls.csv                 # did the denoiser drop cells the aligner kept?
scqc doublet-health --rates rates.csv             # is the rate a measurement or the prior?
scqc quality --valleys valleys.csv --metric umi   # propose a count floor, or refuse to
scqc cluster-preflight --profile clusters.csv --kept 120000   # what step 6 found, before removal
scqc selftest                                     # run the bundled suites
```

**Exit code 2 is a refusal**, 0 is pass or review, 1 is an error, so a gate can stop a shell
script of your own. `--json` on the gate subcommands prints the findings as structured output.

Looking at the data and cutting it are separated in time and recorded separately.

## Install

```bash
git clone https://github.com/JiaenLin/scQC.git && cd scQC
pip install -e .              # the gates and the scqc command — stdlib only, no dependencies
./bin/scqc selftest           # or run it straight from the clone, without installing
```

The decision layer has no third-party dependencies at all, which is deliberate: a gate you cannot
install is a gate that gets skipped. `pip install -e '.[test]'` adds pandas, which only the
cohort-reading test needs.

The tools the pipeline invokes are a separate install, needed to run a cohort end to end:

```bash
setup/install_env.sh --prefix ~/scqc-env --all
```

Four separate environments, because the aligner, the denoiser and the analysis stack have
mutually incompatible pins. Options: `--with-celescope`, `--with-cellbender`, `--with-doublet`, or
`--all`. Requires conda, mamba or micromamba.

The installer finishes by suggesting an `SCQC_ENV_ROOT` export. One script reads it:
`conf/env/build_reference.sh` resolves `--celescope` from `$SCQC_ENV_ROOT/celescope/bin/celescope`
before falling back to `PATH`. Nothing else consults it; the environments are otherwise used
directly by path under the prefix you gave.

## Start a project

```bash
setup/init_project.sh --dir ~/projects/my-study --assay snrna --samples 10
```

scQC never writes into its own directory. One installation serves any number of datasets, and
upgrading the pipeline cannot disturb an existing result.

## The steps

| # | step | does | removes |
|---|---|---|---|
| 0 | ingest | validate samplesheet, resolve reference, verify input is raw | — |
| 1 | ambient | denoise (CellBender), or accept a denoised object; audit the removal; report any fit unlike its siblings | — |
| 2 | cell call | compare aligner and denoiser calls; gate the loss | — |
| 3 | light floor | technical floor for doublet scoring — *not* a quality filter | — |
| 4 | doublets | score per sample, before quality filtering; flag only | — |
| 5 | quality | derive count floors and the per-library mitochondrial ceiling | — |
| 6 | cluster check | per-cluster flags: depth, mitochondrial, markers, doublet | — |
| 7 | **apply** | measure every criterion, write the ledger, then write the filtered object | **yes — only here** |

🔬 **[How each filter is calculated](docs/FILTERS.md)** — the exact procedure for the UMI floor,
the gene floor, the mitochondrial ceiling, doublets and the cluster check: what each is derived
over, what it is applied to, and what it cannot establish.

Doublet scoring precedes quality filtering because scDblFinder's documentation requires it:
*"Further quality filtering should be performed downstream of doublet detection."* It is the only
one of the four common doublet tools that documents an ordering requirement at all.

Ambient correction is **mandatory for single-nuclei** and optional for single-cell.

📊 **[Workflow diagrams](docs/WORKFLOW.md)** — the full flow, the two phases, and where each gate
refuses.

## What makes it reproducible

Each entry says whether it is enforced today, because a claim about a property nobody has
implemented is the same class of defect this pipeline exists to catch.

- **Nothing is overwritten** — *built*. Outputs live under `results/<digest>/`, the digest computed
  from the samplesheet's content, the declared parameters and the mode. A run that would produce
  something different lands somewhere different. A directory records the inputs it was claimed for
  and refuses a run described differently.
- **Gates return a verdict, not a warning** — *built*. Each gate reduces its findings to `REFUSE`,
  `REVIEW` or `PASS`, and the `scqc` subcommands exit 2 on a refusal. A loss that falls ≥3× harder
  on one arm of your design is a `REFUSE` **where it is also material** — the worst arm at ≥1% —
  and a `REVIEW` below that floor, because a ratio between two near-zero rates is dominated by a
  single library and a gate that fires on correct behaviour gets switched off. Reporting and
  enforcing are separate calls: `verdict()` returns the string, the caller stops the run.
- **Every removal is recoverable, with the criterion that made it** — *built*.
  `build_removal_record()` pairs every removed observation with the criteria that fired on it, the
  removal refuses if that record and the mask disagree on the count, and `write_removal_record()`
  persists it as CSV — standard library only, so reading it back needs nothing this pipeline
  installed. The format is documented in `modules/07_apply/apply.py`.
- **Every threshold says who set it** — *built*. Each is recorded `DERIVED` (measured by the
  pipeline) or `ADJUDICATED` (declared in `decisions.yml` with an approver and their own words),
  in the ledger, the written object and the report. A value with a number but no approver and no
  words counts as neither, and the derived value is used instead: half an approval must not read
  as a whole one.
- **A declared removal is checked against the thresholds it was given for** — *built*. Where
  `decisions.yml` supplies an approval, it is matched against the action the current thresholds
  derive; move a threshold after approving and the approval no longer applies, by design. There is
  no force flag.
- **Reports carry the commit and the tool versions**, obtained by asking each tool rather than read
  from a lockfile — *built*. The report also audits itself: anything the run should have recorded
  and did not is counted as a defect on its own front page.
- **A report older than its inputs is refused rather than warned about** — *specified, not fed*.
  `freshness()` and `refuse_if_stale()` exist in `report/build.py`, and no step supplies a
  newest-input time, so every report says `NOT CHECKED` rather than claiming to be current.
  [docs/REPORT_DESIGN.md](docs/REPORT_DESIGN.md) is the specification.

See **[docs/PRINCIPLES.md](docs/PRINCIPLES.md)** for the four rules and why each exists.

## Documentation

| document | answers |
|---|---|
| [docs/FILTERS.md](docs/FILTERS.md) | **how each filter is calculated** — the exact procedure, the population it is derived over, and what it cannot establish |
| [docs/WORKFLOW.md](docs/WORKFLOW.md) | diagrams: the pipeline, the two phases, the parameter classes |
| [docs/PRINCIPLES.md](docs/PRINCIPLES.md) | the removal checklist and the three other enforced rules |
| [docs/REPORT_DESIGN.md](docs/REPORT_DESIGN.md) | the report layout and the nine figures. The layout is built; **the figures are not** |
| [docs/TOOLS_AND_REFERENCES.md](docs/TOOLS_AND_REFERENCES.md) | tools, versions, reference resolution |
| [CALIBRATION.md](CALIBRATION.md) | what was measured, how much it varied, what one cohort cannot establish |
| [CONTRIBUTING.md](CONTRIBUTING.md) | what a pull request has to answer before it can remove anything |
| [tests/acceptance/](tests/acceptance/) | how to check scQC against a dataset you already trust |

## Tests

```bash
scqc selftest            # every bundled suite, including the adversarial one; exit 1 if any failed
```

`selftest` reports each suite as PASS, FAIL or SKIP and returns a non-zero exit code if any
failed. **A SKIP is not a pass** and is reported as its own outcome — `tests/test_audit_ambient.py`
skips unless pandas is installed *and* `COHORT_DIR` points at a cohort whose tables it can read.

To run them by hand instead, accumulate the failures rather than looping bare — a plain
`for t in …; do python "$t"; done` exits with the status of the *last* suite, so a failure
anywhere before it is invisible:

```bash
fails=0
for t in tests/test_*.py; do python "$t" || fails=$((fails + 1)); done
python tests/adversarial.py || fails=$((fails + 1))
echo "$fails suite(s) failed"; [ "$fails" -eq 0 ]
```

Everything except the cohort audit runs on the standard library alone. **pandas** is needed by
`tests/test_audit_ambient.py` and by the acceptance harness, and by nothing else:

```bash
pip install -e '.[test]'
COHORT_DIR=/path/to/cohort python tests/acceptance/run_tier1.py
```

### Which interpreter runs `scqc`

**`core`'s.** Not the system python.

```bash
$SCQC_ENV_ROOT/core/bin/python bin/scqc run     --project my_project --mode evidence     --python $SCQC_ENV_ROOT/core/bin/python
```

`core` is installed unconditionally because the orchestrator needs the analysis stack in its own
process — several steps read matrices directly. Under a bare system python they fail with
`ModuleNotFoundError` partway through a run rather than at the start, which is the worst place
for an environment problem to appear.

The other environments exist for tools whose pins are incompatible with `core` and with each
other. They are reached through `--celescope`, `--cellbender` and `--rscript`, and are never
imported.


The acceptance harness is a regression test against a dataset **you** supply — no data ships with
this repository, so it cannot be run from a clone alone.

The adversarial suite exists because re-reading one's own code finds almost nothing: every defect
it has caught came from calling a function with a hostile input and looking at what happened.

## Status

**0.1.0.** Being precise about this, because a QC tool that overstates itself does damage
quietly — and one that understates itself is wrong in the same way, just harder to notice. Every
row below was checked against the tree rather than remembered.

| | |
|---|---|
| ✅ **Built and tested** | The decision layer: each step's policy, contract, threshold derivation and gate, importable as Python, plus the `scqc` CLI over it. `scqc selftest` on a clean clone reports **15 passed, 0 failed, 1 skipped** — the skip needs pandas or a cohort (`COHORT_DIR`) and is not evidence either way. |
| ✅ **Built** | The `scqc` command, as **per-step subcommands** over tables you supply: `validate`, `verify`, `gate-cells`, `doublet-health`, `quality`, `cluster-preflight`, `selftest`. Exit code 2 on a refusal, `--json` for structured output. |
| ✅ **Built** | Recoverable removal. Step 7 pairs each removed observation with the criteria that fired on it, refuses if that record and the mask disagree, and writes it as a CSV any reader can open. |
| ✅ **Built, and run end to end** | The two-mode driver. `scqc run --project … --mode evidence\|apply` builds the task graph and runs it, locally or on PBS with one job per task. A ten-library cohort completes as 47 tasks with nothing reused. In `evidence` mode the apply task is **not placed in the graph at all**, so there is no code path from it to a removal. |
| ✅ **Built** | The execution layer. Steps invoke the aligner, the denoiser, the doublet caller and the analysis stack out of process, each under its own interpreter, and read count matrices through the adapters. Versions are obtained by asking the tool, never read from a lockfile. |
| ✅ **Built** | The report. Every run writes `qc_report.html` and `report.json`, including a per-library table of every threshold the run derived with each column marked per-library or cohort constant. The report **audits itself**: anything the payload should have carried and did not is a defect counted on its own front page. |
| ✅ **Built** | Decisions are read. `decisions.yml` is parsed and validated, and `--mode apply` refuses on a missing or incomplete one, naming every problem at once rather than one per run. |
| ⚠️ **Not built — the figures** | `report/figures.py` exists and the report expects F1–F9, but no step supplies one. Each absence is reported as a defect rather than omitted. The five sections and the nine figures in [docs/REPORT_DESIGN.md](docs/REPORT_DESIGN.md) remain a specification. |
| ⚠️ **Not fed — freshness** | `freshness()` and `refuse_if_stale()` exist in `report/build.py`, and no step supplies a newest-input time, so every report says `NOT CHECKED` rather than claiming to be current. Of everything on this list it is the one whose absence is hardest to notice, because a stale artifact opens and reads exactly like a current one. |
| ✅ **Built** | Step 7, the only step that removes. It measures every criterion per barcode, writes the removal ledger, verifies the ledger against the mask, audits the result, and only then writes the filtered cohort object. Measure, record, write — in that order, so nothing is materialised before what left has been written down. Where `decisions.yml` supplies an approval, it is additionally matched against the action the current thresholds derive. |
| ✅ **Built** | Content-addressed outputs. Each run writes under `results/<digest>/`, named from the samplesheet's content, the declared parameters and the mode. A run that would produce something different lands somewhere different, so nothing is overwritten by one that disagrees with it. |
| ⚠️ **Not measured** | Run-to-run tolerance for the stochastic steps. Recorded as `UNMEASURED`; the pipeline will not assert a tolerance it has not measured. |
| ⚠️ **n = 1 cohort** | One tissue, one species, one platform. Every threshold is an existence proof, not a range. See [CALIBRATION.md](CALIBRATION.md). |

## Contributing

Issues and pull requests welcome. See **[CONTRIBUTING.md](CONTRIBUTING.md)**; the two expectations
that catch most reviews are that any change removing observations goes through step 7 and answers
the removal checklist in the pull request, and that a new gate carries a materiality bound.

## Citing

See [CITATION.cff](CITATION.cff).

## License

MIT — see [LICENSE](LICENSE).
