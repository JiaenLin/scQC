# scQC

**A quality-control pipeline for single-cell and single-nuclei RNA-seq that separates deriving a
threshold from applying it.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-0.3.1-blue.svg)](#status)

Most QC pipelines take thresholds as arguments. scQC treats them as findings: it measures what it
can from your data, refuses to guess what it cannot, and records who set every threshold it
applies — the data, or a person in their own words — so the two can never be confused later.

> **Read [Status](#status) before you plan a run.** At `0.3.1` scQC runs a cohort end to end —
> `scqc run` builds the task graph and executes it, locally or on a PBS scheduler with one job per
> task, invoking the aligner, the denoiser, the doublet caller and the analysis stack out of
> process. It writes a report and, in apply mode, one filtered object per library plus a merged
> cohort object, with a ledger naming every barcode removed and why.
> It draws all fifteen of the report's figures, and it now supplies its own newest-input time, so
> the report states whether it is fresh instead of declining to say. The per-step subcommands
> remain, for judging tables you produced elsewhere.
>
> **[KNOWN_ISSUES.md](KNOWN_ISSUES.md) lists what is measured, reproduced and not yet fixed.**
> Read it before quoting a number: it is where a defect lives between being found and being
> repaired, and one entry there is the difference between a figure that means something and one
> that does not.

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
  the removal ledger, and materialise the filtered objects — **one per library and one combined**,
  each retained nucleus carrying its cluster and that cluster's flags from step 6.

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
scqc stamp results/*/objects/*.h5ad               # declare what the flag means, on objects
                                                  # written before 0.3.1
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
| 6 | cluster check | per-cluster flags: depth, mitochondrial, markers, doublet — at resolution **2.0**, with each `--extra-resolutions` value profiled into a sibling table beside it | — |
| 7 | **apply** | measure every criterion, write the ledger, then write the filtered objects | **yes — only here** |

### Why each step is there

None of these is a stage in a conventional pipeline that happens to be numbered. Each exists
because something specific goes wrong without it, and each was placed where it is by a constraint
rather than by preference.

- **0 · ingest** — a delivered matrix is very often the aligner's `outs/filtered`: genuine integer
  counts that have *already* been through cell calling, and nothing about the file says so — not
  its name, not its shape, not a header field. That matters immediately, because an ambient model
  learns its background profile **from the empty droplets**. Hand it a cell-called matrix and
  "ambient correction is not possible with these files" becomes available as a conclusion, and it
  reads on the page exactly like a real finding.
- **1 · ambient** — a nucleus holds roughly an order of magnitude less RNA than the cell it came
  from, and the ambient pool in a nuclear prep *is* the cytoplasm of the cells lysed to make it.
  A pipeline that leaves that in place for snRNA-seq is not applying a lighter touch; it is
  analysing a mixture and calling it a cell. The audit exists because the headline fraction cannot
  see failure: ten good libraries spanned 15.5–23.7% while the one genuinely degenerate fit sat at
  8.3% — inside a plausible range, and *low* rather than high.
- **2 · cell call** — the error here is not symmetric. A nucleus wrongly called empty never reaches
  the steps that would have judged it; a droplet wrongly called a cell still has to pass the count
  floors, the ceiling, the doublet call and the cluster check. Permissive is the safe direction,
  and the gate enforces exactly that asymmetry. Its third check exists because a loss can be small
  and still be structured: all 617 cells lost in the calibration cohort fell on **one** level of
  the design factor, which no per-sample percentage would have shown.
- **3 · light floor** — doublet detectors build their null by summing pairs of observed
  transcriptomes. If the pool contains near-empty droplets, the artificial doublets are
  debris+cell or debris+debris, which are not doublets, and the null is calibrated on the wrong
  thing. It is emphatically not a quality filter: a barcode below it is recorded **UNSCORED**,
  which is unknown, not a singlet.
- **4 · doublets** — a consensus of detectors was measured and rejected: scds returned 6.00% in
  all ten libraries (it is a fixed quantile), DoubletDetection ranged 1.82–20.90%. *A vote does
  not protect against a bad detector; it averages one in.* The protection is a diagnostic instead.
  `dbr` has no default because scDblFinder's own default is the 10x loading formula, which applied
  to a Singleron run tracked library size at r = 0.872 and would have imported a 19.65% vs 14.98%
  condition differential onto the study's primary readout.
- **5 · quality** — the clearest case in the pipeline of a procedure transferring where a number
  cannot. The count floors are *measured* as the density valley between the debris and cell modes,
  and bimodality is tested rather than assumed, because a density minimum exists in any smooth
  curve. On snRNA the mitochondrial `k` is **declared**, not derived, and the reason is biology: a
  nucleus contains no mitochondria, so `pct_counts_mt` there measures cytoplasmic carry-over —
  deriving `k` from the cohort's own tail would let a badly-prepared cohort earn itself a wider
  licence.
- **6 · cluster check** — this is the **only** point at which the doublet-fraction criterion is
  computable. After step 7 every cluster is 0% doublet by construction and that flag becomes a
  tautology which passes for the wrong reason. The flags are a conjunction rather than a
  disjunction because a mitochondrial-marker cluster at normal coverage is the signature of a
  mitochondria-rich cell type, not of damage; that change took 22 flagged clusters to 8.
- **7 · apply** — confining every removal to one step is what lets a reviewer establish, quickly,
  that no other step can silently cost them cells. The internal ordering is load-bearing too: an
  object filtered before its ledger was checked would have performed the removal before anything
  verified it, so the ledger is written **before** the object — a crash between the two leaves a
  record of a removal that did not happen, rather than a removal with no record.

Apply mode writes one filtered object per library, and a cohort object that is the **merge of
those files read back from disk**. Every retained nucleus carries `sample`, `cluster`,
`cluster_FLAG`, `cluster_WATCH` and the continuous values behind them, so nothing downstream has
to re-cluster to recover what step 6 found — and a re-clustering of the filtered object would not
give the same answer anyway, since criterion D is a tautology once the doublets are gone.

**Every object also carries `uns["scqc"]`: a declaration of what the flag MEANS**, how many
nuclei carry it, and a digest of the exact set. Without one, a downstream tool has two bad
options — ignore the flag and annotate nuclei this pipeline flagged as technical, or guess from
the column name, which is that tool deciding what is technical on scQC's behalf. The declaration
is the third option: scQC states its decision, the consumer decides what to do about it, and the
digest lets the consumer prove the column is still the one scQC wrote. `scqc stamp <objects>`
adds it to anything written before `0.3.1` without a re-run.

The per-library objects are primary and the cohort object is derived from them: each library is
filtered on its own mitochondrial ceiling and cluster-checked on its own clustering.
📄 **[Output reference](docs/OUTPUTS.md)** — every file, every column.

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
- **A report older than its inputs is refused rather than warned about** — *built*.
  `engine/pipeline.py` supplies the newest-input time over the files the run actually read, and
  `freshness()` compares it against the generated time, recording `stale` as True, False **or
  None** — never False by default, because None is NOT CHECKED and NOT CHECKED is not a pass.
  `refuse_if_stale()` acts on it, and a run that could not check says so with its reason.
  [docs/REPORT_DESIGN.md](docs/REPORT_DESIGN.md) is the specification.

See **[docs/PRINCIPLES.md](docs/PRINCIPLES.md)** for the four rules and why each exists.

## Documentation

| document | answers |
|---|---|
| [docs/QUICKSTART.md](docs/QUICKSTART.md) | **start here** — install, samplesheet, run, and how to read the report |
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | **operating it** — every command, run keys and resuming, approvals, and what a refusal means |
| [docs/FILTERS.md](docs/FILTERS.md) | **how each filter is calculated** — the exact procedure, the population it is derived over, and what it cannot establish |
| [docs/OUTPUTS.md](docs/OUTPUTS.md) | **every file a run writes** — the object `obs` schema, every table's columns, and which files are meant to be read next |
| [docs/WORKFLOW.md](docs/WORKFLOW.md) | diagrams: the pipeline, the two phases, the parameter classes |
| [docs/PRINCIPLES.md](docs/PRINCIPLES.md) | the removal checklist and the three other enforced rules |
| [docs/REPORT_DESIGN.md](docs/REPORT_DESIGN.md) | the report layout and its thirteen figures. All are built and all are drawn from a finished run's own tables; the one that is not drawn by default (F5) needs a sweep that re-scores every library, and the report names the flag that requests it |
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

**0.3.1.** Being precise about this, because a QC tool that overstates itself does damage
quietly — and one that understates itself is wrong in the same way, just harder to notice. Every
row below was checked against the tree rather than remembered.

| | |
|---|---|
| ✅ **Built and tested** | The decision layer: each step's policy, contract, threshold derivation and gate, importable as Python, plus the `scqc` CLI over it. **27 suites** — 26 `tests/test_*.py` plus the adversarial one. **Run `scqc selftest` yourself and read the counts off your own run**; this table does not quote a number, because a test count copied into prose is a number that goes stale between releases while still looking authoritative. **Use an interpreter that has the scientific stack**: `bin/scqc` is `#!/usr/bin/env python3`, so under a bare system python the suites needing numpy/anndata SKIP rather than run, and a skip counted as a pass is the defect this project exists to catch. |
| ✅ **Built** | The `scqc` command, as **per-step subcommands** over tables you supply: `validate`, `verify`, `gate-cells`, `doublet-health`, `quality`, `cluster-preflight`, `selftest`. Exit code 2 on a refusal, `--json` for structured output. |
| ✅ **Built** | Recoverable removal. Step 7 pairs each removed observation with the criteria that fired on it, refuses if that record and the mask disagree, and writes it as a CSV any reader can open. |
| ✅ **Built, and run end to end** | The two-mode driver. `scqc run --project … --mode evidence\|apply` builds the task graph and runs it, locally or on PBS with one job per task. A ten-library cohort completes as 47 tasks with nothing reused. In `evidence` mode the apply task is **not placed in the graph at all**, so there is no code path from it to a removal. |
| ✅ **Built** | The execution layer. Steps invoke the aligner, the denoiser, the doublet caller and the analysis stack out of process, each under its own interpreter, and read count matrices through the adapters. Versions are obtained by asking the tool, never read from a lockfile. |
| ✅ **Built** | The report. Every run writes `qc_report.html` and `report.json`, including a per-library table of every threshold the run derived with each column marked per-library or cohort constant. The report **audits itself**: anything the payload should have carried and did not is a defect counted on its own front page. |
| ✅ **Built** | Decisions are read. `decisions.yml` is parsed and validated, and `--mode apply` refuses on a missing or incomplete one, naming every problem at once rather than one per run. |
| ✅ **Built** | The figures — all fifteen ids, `F1`–`F15`, drawn by ten functions. Four ids deliberately share a function with another: applying the same treatment twice reads as a comparison, so the gene axis and the two other applied axes get their own ids rather than being folded into one. `report/collect.py` assembles each from a NAMED table a reader can open, and `tests/test_figure_collection.py` renders every one of them from a fixture, because a builder that returns a dict can still be unrenderable. A figure that genuinely cannot be drawn is shown as a named absence saying what would produce it — not as a gap. |
| ⚠️ **Known defect** | Step 6's cluster flags were computed on the wrong population until `02a422d`; cohorts processed before it must be re-run before their flags are read. See [KNOWN_ISSUES.md](KNOWN_ISSUES.md). |
| ✅ **Built** | Freshness. `engine/pipeline.py::_newest_input_block()` supplies the newest-input time over the files the run read — listing an absent file rather than skipping it, so a time computed over a missing input cannot silently read as fresh — and `report/build.py` records `generated`, `newest_input`, `checked`, `stale`, a reason and a margin in `provenance.freshness`. `stale` is three-valued and never False by default. This was the last row on this table to move from *specified* to *built*, and it was the one whose absence would have been hardest to notice: a stale artifact opens and reads exactly like a current one. |
| ✅ **Built** | Step 7, the only step that removes. It measures every criterion per barcode, writes the removal ledger, verifies the ledger against the mask, audits the result, and only then writes the filtered cohort object. Measure, record, write — in that order, so nothing is materialised before what left has been written down. Where `decisions.yml` supplies an approval, it is additionally matched against the action the current thresholds derive. |
| ✅ **Built** | The delivered object declares its own flag. `uns["scqc"]` records the flag column, what it MEANS, how many carry it, how many were examined, and a digest of the exact mask - so a downstream tool can act on it knowing whose decision it is rather than guessing from a column name. The digest is a cross-tool contract with scAnno, held by a known-answer vector asserted in both suites. `scqc stamp` backfills objects written earlier. |
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
