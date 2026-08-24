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
`--all`. Requires conda, mamba or micromamba.

## Run

```bash
setup/init_project.sh --dir ~/projects/my-study --assay snrna --samples 10

scqc run --project ~/projects/my-study --mode evidence   # measure, decide nothing
# read the report, record approvals in decisions.yml
scqc run --project ~/projects/my-study --mode apply      # filter, once
```

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

**Evidence mode** measures every criterion, writes the report, and removes nothing.
**Apply mode** takes the same measurements plus the recorded approvals and performs the removal.
An approval is matched against the thresholds it was given for; move a threshold afterwards and
the approval no longer applies. There is no force flag.

## Steps

| # | step | does | removes |
|---|---|---|---|
| 0 | ingest | validate samplesheet, resolve reference, verify input is raw | — |
| 1 | ambient | denoise (CellBender) or accept a denoised object; audit the removal | — |
| 2 | cell call | compare aligner and denoiser calls; gate the loss | — |
| 3 | light floor | technical floor for doublet scoring — not a quality filter | — |
| 4 | doublets | score per sample, before quality filtering; flag only | — |
| 5 | quality | derive count floors and the per-library mitochondrial ceiling | — |
| 6 | cluster check | per-cluster flags: depth, mitochondrial, markers, doublet fraction | — |
| 7 | **apply** | measure, write the ledger, then write the filtered objects | **yes — only here** |

Ambient correction is mandatory for single-nuclei and optional for single-cell. Doublet scoring
precedes quality filtering, as scDblFinder's documentation requires.

📄 **[Why each step exists](docs/RATIONALE.md)** · **[How each filter is
calculated](docs/FILTERS.md)** · **[Workflow diagrams](docs/WORKFLOW.md)**

## Output

Apply mode writes one filtered object per library plus a cohort object merged from those files
read back from disk. Every retained nucleus carries `sample`, `cluster`, `cluster_FLAG`,
`cluster_WATCH` and the continuous values behind them, so nothing downstream needs to re-cluster.

Every object carries `uns["scqc"]`: what the flag means, how many nuclei hold it, and a digest of
the exact set — so a consumer can act on the flag without inferring its meaning from a column
name, and can verify the column is still the one scQC wrote. `scqc stamp <objects>` adds it to
objects written before 0.4.0.

Per-library objects are primary; the cohort object is derived from them, since each library is
filtered on its own mitochondrial ceiling and cluster-checked on its own clustering.

📄 **[Output reference](docs/OUTPUTS.md)** — every file, every column.

## Reproducibility

All enforced today:

- **Nothing is overwritten.** Outputs land under `results/<digest>/`, the digest computed from the
  samplesheet content, the declared parameters and the mode. A directory records the inputs it was
  claimed for and refuses a run described differently.
- **Gates return a verdict, not a warning.** Each reduces to `REFUSE`, `REVIEW` or `PASS`, and the
  subcommands exit 2 on refusal. Reporting and enforcing are separate calls.
- **Every removal is recoverable with the criterion that made it.** Removed observations are paired
  with the criteria that fired; the removal refuses if that record and the mask disagree on the
  count. Persisted as CSV, readable with the standard library alone.
- **Every threshold records who set it** — `DERIVED` or `ADJUDICATED` with an approver and their
  own words — in the ledger, the object and the report.
- **Reports carry the commit and tool versions**, obtained by asking each tool rather than reading
  a lockfile, and audit themselves: anything the run should have recorded and did not is counted
  as a defect on the front page.
- **A report older than its inputs is refused.** Freshness is recorded as True, False or **None** —
  never False by default, because unchecked is not a pass.

📄 **[Principles](docs/PRINCIPLES.md)** — the four enforced rules.

## Tests

```bash
scqc selftest     # every bundled suite including the adversarial one; non-zero exit if any failed
```

Each suite reports PASS, FAIL or SKIP. A SKIP is reported as its own outcome, not as a pass.
Everything except the cohort audit runs on the standard library alone; `pip install -e '.[test]'`
adds pandas for the rest.

Run `scqc` with `core`'s interpreter, not the system python — several steps read matrices in the
orchestrator's own process:

```bash
$SCQC_ENV_ROOT/core/bin/python bin/scqc run --project my_project --mode evidence \
    --python $SCQC_ENV_ROOT/core/bin/python
```

The acceptance harness regression-tests against a dataset you supply; no data ships with this
repository.

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
barcode removed and why. All eighteen report figures are drawn, three of them a confounding
block: which design factors these libraries cannot tell apart, how far apart the arms that
leaves sit on every QC metric, and whether the filter widened the gap. Per-step subcommands
remain, for judging tables produced elsewhere.

Read [KNOWN_ISSUES.md](KNOWN_ISSUES.md) before quoting a number from a run.

## Citing

If scQC contributes to published work, cite the repository and the version you ran. Every report
carries the commit.

## License

MIT — see [LICENSE](LICENSE).
