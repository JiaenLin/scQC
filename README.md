# scQC

**A quality-control pipeline for single-cell and single-nuclei RNA-seq that separates deriving a
threshold from applying it.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-0.0.1--dev-orange.svg)](#status)

Most QC pipelines take thresholds as arguments. scQC treats them as findings: it measures what it
can from your data, refuses to guess what it cannot, and requires a recorded human decision — in
that person's own words — before anything is removed.

> **Read [Status](#status) before you plan a run.** At `0.0.1-dev` scQC is a decision layer, not a
> workflow engine. The `scqc` command judges tables you have produced elsewhere, one step at a
> time; it does **not** read matrices, does **not** invoke CeleScope, Cell Ranger, CellBender,
> scDblFinder or scanpy, and does **not** write a report. The two-mode driver
> (`--mode evidence` / `--mode apply`) described below is the intended interface and is not built.

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

- **evidence** — run every step, derive every DERIVED parameter, write the report and the decision
  template, remove nothing, then stop.
- **apply** — re-run with your `decisions.yml`, verify that every ADJUDICATED value carries your
  own words, and only then remove.

```bash
scqc --project ./my-study --mode evidence   # NOT BUILT — measures everything, removes nothing
cp decisions.template.yml decisions.yml     # NOT BUILT — record your decisions
scqc --project ./my-study --mode apply      # NOT BUILT — applies them; refuses if any is missing
```

> **The two-mode driver does not exist.** `scqc` takes subcommands, not `--mode`; nothing runs the
> steps in sequence and nothing writes a report or a decision template. The block above is the
> interface the modules are written against, and it is also what `setup/init_project.sh` writes
> into a new project's README, where it will not run either. See [Status](#status).

What `scqc` does today is judge one step at a time, from tables you produce yourself:

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

**Exit code 2 is a refusal**, 0 is pass or review, 1 is an error — so a gate can stop a shell
script today even though nothing drives the steps end to end. `--json` on the gate subcommands
prints the findings as structured output.

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

The tools the pipeline reasons *about* are a separate install, and only needed once the execution
layer exists:

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
| 1 | ambient | denoise (CellBender); audit the result; halve the learning rate if degenerate | — |
| 2 | cell call | compare aligner and denoiser calls; gate the loss | — |
| 3 | light floor | technical floor for doublet scoring — *not* a quality filter | — |
| 4 | doublets | score per sample, before quality filtering; flag only | — |
| 5 | quality | derive count floors and mitochondrial ceiling | — |
| 6 | cluster check | per-cluster flags: depth, mitochondrial, markers, doublet | — |
| 7 | **apply** | pre-flight, verify approval, remove | **yes — only here** |

Doublet scoring precedes quality filtering because scDblFinder's documentation requires it:
*"Further quality filtering should be performed downstream of doublet detection."* It is the only
one of the four common doublet tools that documents an ordering requirement at all.

Ambient correction is **mandatory for single-nuclei** and optional for single-cell.

📊 **[Workflow diagrams](docs/WORKFLOW.md)** — the full flow, the two phases, and where each gate
refuses.

## What makes it reproducible

Three of these are enforced in code today and two are design commitments the execution layer will
have to honour. Each says which, because a claim about a property nobody has implemented is the
same class of defect this pipeline exists to catch.

- **Gates return a verdict, not a warning** — *built*. Each gate reduces its findings to `REFUSE`,
  `REVIEW` or `PASS`, and the `scqc` subcommands exit 2 on a refusal. A loss that falls ≥3× harder
  on one arm of your design is a `REFUSE` **where it is also material** — the worst arm at ≥1% —
  and a `REVIEW` below that floor, because a ratio between two near-zero rates is dominated by a
  single library and a gate that fires on correct behaviour gets switched off. Reporting and
  enforcing are separate calls: `verdict()` returns the string, the caller stops the run.
- **A removal needs the operator's own words for that exact action** — *built*.
  `modules/07_apply/apply.py` raises unless a recorded approval matches the action text verbatim,
  and there is no force flag.
- **Every removal is recoverable, with the criterion that made it** — *built*.
  `build_removal_record()` pairs every removed observation with the criteria that fired on it, the
  gate refuses if that record and the mask disagree on the count, and `write_removal_record()`
  persists it as CSV — standard library only, so reading it back needs nothing this pipeline
  installed. The format is documented in `modules/07_apply/apply.py`.
- **Decisions are a versioned input**, not a command line, so that the same data and the same
  `decisions.yml` are the same run — *design*. Nothing reads or hashes a decisions file yet.
- **Reports carry the commit**, the tool versions, the reference registry entry and the hash of the
  decisions file, and a report older than its inputs is refused rather than warned about —
  *design*. No code writes a report and nothing compares timestamps.
  [docs/REPORT_DESIGN.md](docs/REPORT_DESIGN.md) is the specification.

See **[docs/PRINCIPLES.md](docs/PRINCIPLES.md)** for the four rules and why each exists.

## Documentation

| document | answers |
|---|---|
| [docs/WORKFLOW.md](docs/WORKFLOW.md) | diagrams: the pipeline, the two phases, the parameter classes |
| [docs/PRINCIPLES.md](docs/PRINCIPLES.md) | the removal checklist and the three other enforced rules |
| [docs/REPORT_DESIGN.md](docs/REPORT_DESIGN.md) | **specification, not shipped behaviour** — the report layout, the nine figures and the reproducibility block the report writer must produce |
| [docs/TOOLS_AND_REFERENCES.md](docs/TOOLS_AND_REFERENCES.md) | tools, versions, reference resolution |
| [CALIBRATION.md](CALIBRATION.md) | what was measured, how much it varied, what one cohort cannot establish |
| [CONTRIBUTING.md](CONTRIBUTING.md) | what a pull request has to answer before it can remove anything |
| [tests/acceptance/](tests/acceptance/) | how to check scQC against a dataset you already trust |

## Tests

```bash
scqc selftest            # ten unit suites + the adversarial suite; exit 1 if any failed
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
quietly — and one that *understates* itself is also wrong, which half this table was until the
driver it says does not exist ran a cohort end to end. Every row below was checked against the
tree rather than remembered.

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
| ⚠️ **Not built — the deliverable object** | Step 7 validates the decisions, runs its pre-flight, and then refuses: there is no combined object for it to filter. **Nothing in this pipeline has ever removed an observation, in either mode.** |
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
