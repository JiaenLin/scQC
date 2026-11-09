# User guide

How to operate scQC: the commands, how a run resumes, what a refusal means, and how to answer one.

New here? Start with the **[quick start](QUICKSTART.md)**.

---

## Commands

| command | what it does |
|---|---|
| `scqc validate --samplesheet <file> --registry <file>` | check a samplesheet before anything is computed |
| `scqc run --project <dir>` | run the pipeline over a project |
| `scqc report <results-dir>` | rebuild a finished run's report from its own files |
| `scqc verify` | is a matrix raw? — from summary statistics you supply |
| `scqc gate-cells` | did the denoiser drop cells the aligner kept? |
| `scqc doublet-health` | is a doublet rate a measurement, or the prior it was given? |
| `scqc quality` | propose a count floor from per-library valleys |
| `scqc cluster-preflight` | contradictions to settle before a removal is applied |
| `scqc selftest` | run the bundled suites |

The middle five take numbers you already have and answer one question each. They exist so a
decision can be checked without running a pipeline — and so the same check gives the same answer
whether it is reached through the pipeline or by hand.

### `scqc run`

```bash
scqc run --project ~/projects/my-study [--jobs N] [--mode evidence|apply] [--decisions FILE] \
          [--registry FILE] [--executor pbs --queue Q [--gpu-queue G]]
```

- `--registry` — the project's reference registry (TSV: `species`, `build`, `path`, ...). The
  samplesheet's `reference` must be a `species/build` key in it. **Required:** ingest refuses
  every library without it, and scQC ships no registry to fall back to.
- `--gpu-queue` — under `--executor pbs`, the queue a task that declares a GPU (the denoiser)
  is sent to; every other task goes to `--queue`. Without it, every task goes to `--queue`.

- `--jobs` — independent tasks run at once. `0` (default) uses the machine, capped at 16. Use `1`
  when a failure has to be read in a single log.
- `--mode` — `apply` (default) applies the DERIVED and DECLARED thresholds, records who decided
  each one, and writes filtered objects; `evidence` computes and reports every threshold but writes
  no filtered object.
- `--decisions` — optional; a YAML file recording overrides a person has approved, and their words.
  See [approvals](#approvals-and-adjudicated-values).

> This guide documented the second mode as `measure` until 2026-08-11. The CLI has only ever
> accepted `evidence`, so `--mode measure` failed at argument parsing — a reader following this
> page could not run the command it gave them. `tests/test_docs.py` did not catch it because it
> validates the numeric constants quoted in `docs/FILTERS.md`, not prose or option names.

## Run keys and resuming

**Outputs are stored under a name derived from what produced them.** The run key is a digest of
the samplesheet, the declared parameters, the mode and the pipeline's OUTPUT-BEARING CODE - the
bytes of `engine/`, `modules/`, `adapters/`, `report/`, `lib/`, `scqc_cli.py`, `conf/env/` and
`VERSION` (`engine/runkey.py`, `OUTPUT_BEARING`). A commit that changes only docs or tests keeps
the key, so completed work is reused; an edit to code that runs moves it, committed or not. The
commit and dirty flag are recorded beside the key in INPUTS.json (`code_records`) and in
STATUS.json:

```
results/8f2a91c4d7e0/              ← this samplesheet, these parameters, this mode, this code
results/8f2a91c4d7e0/INPUTS.json   ← what that digest was computed from
```

Same inputs → same directory, and completed work is reused. Change a declared parameter, or the
code, and the digest changes with it, so the new run writes **beside** the old one rather than
over it. **The decisions file is the exception:** it is not in the digest, so an apply run with a
changed `decisions.yml` runs in the same directory, re-runs the tasks that carry it (steps 5 to 7)
and replaces their outputs.

A task is skipped when all three hold: it completed before, its signature is unchanged, and every
output it recorded is still on disk. **Delete an output to force its task to re-run** — deleting
`reports/qc_report.html` re-runs the report step and nothing else, which takes seconds.

Gate findings are stored in `work/<digest>/state.json` and survive a skip, so a resumed run's
report carries the findings the original raised. A step that re-runs replaces its own.

## Modes

**`evidence`** derives every threshold, writes every table and figure, and removes nothing. Use it
to see what the data proposes before committing.

**`apply`** does all of that and then writes the filtered objects. It is the default, and it is
the only mode that removes anything. Before it writes, it:

1. re-derives every criterion per cell and writes the **removal ledger** — one row per removed
   barcode naming the criteria that removed it;
2. runs an audit over the completed decision (does each criterion remove anything uniquely? does
   `removed` decompose exactly into the recorded criteria? is anything unexamined reaching the
   deliverable?);
3. writes one object per library, then builds the cohort object by **reading those files back**.

The per-library objects are primary; the cohort object is derived from them. Each library is
filtered on its own mitochondrial ceiling and cluster-checked on its own clustering, so a pooled
object is a merge of per-library results, not a re-analysis.

## Approvals and ADJUDICATED values

Every parameter in the report carries a class:

| class | meaning |
|---|---|
| `FIXED` | true regardless of dataset; changing it is a code change |
| `DERIVED` | computed from this dataset; reproducible from the data alone |
| `DECLARED` | supplied by the operator up front, before seeing the result |
| `ADJUDICATED` | a person decided, after seeing the result, in their own words |

`decisions.yml` is optional: without it apply mode applies what the data proposed (DERIVED) and
what was declared (DECLARED), and says so. An ADJUDICATED override needs a `value`, an
`approved_by` and the approver's `verbatim` words. A threshold nobody can attribute is a threshold
nobody can review, and "someone approved this" is not a record of who or why: a value missing
either attestation is not an override, and the derived value is applied instead.

## When a gate refuses

A REFUSE stops the run. A REVIEW does not — it means a person must look, not that the run is
wrong, and the distinction survives all the way into the report rather than being collapsed into
pass/fail. Collapsing it trains a reader to ignore both.

**A run that refuses still writes its report.** Refusing to write the document would destroy the
only record of why the run stopped.

Read what the refusal says before deciding it is wrong. Each one names the numbers behind it and
what would have to be true for it to pass. The common ones:

- **"the design differential exceeds 3×"** — a removal is falling much harder on one arm of the
  design than another. This is the check that exists to catch a technical gradient being laid
  down where the biology is measured.
- **"`mt_prefix` is not declared"** — species-specific and never guessed. See the quick start.
- **"the declared `cellbender_barcodes` disagrees with the object"** — the CSV and the denoised
  object are from different runs, so comparing the aligner against either says nothing.
- **"a criterion removes nothing another does not"** — not necessarily wrong, but that criterion
  is not carrying its own weight, and the report says so rather than leaving it to be noticed.

## Reading the outputs

```
results/<digest>/
├── objects/
│   ├── cohort.deliverable.h5ad          the merge of the per-library objects
│   └── cohort_per_sample/<S>.filtered.h5ad
├── tables/
│   ├── removal_ledger.csv               one row per removed barcode, and why
│   ├── thresholds_per_sample.csv        every threshold, with its scope
│   ├── <S>.percell.csv                  every value each criterion was evaluated on
│   └── ...
└── reports/
    ├── qc_report.html                   the document
    ├── report.json                      every number in it, machine-readable
    └── payload.json                     what the document was built from
```

Every retained observation - nucleus or cell - carries `sample`, `cluster`, `cluster_FLAG`,
`cluster_WATCH` and the continuous values behind them, so nothing downstream has to re-cluster to recover what step 6
found. Unknown is kept distinct from absent throughout: a nullable boolean is `pd.NA`, not
`False`; a missing number is `NaN`, not `0.0`.

**If gene symbols are not unique in your reference**, `var` is re-indexed by the identifier and
the symbols are kept in `var["gene_symbol"]`. Symbol matching on a delivered object — `mt-`,
`^Rp[sl]` — should read that column. A reference whose symbols are unique is left as it arrived.

📄 Full column-by-column reference: **[OUTPUTS.md](OUTPUTS.md)**.

## Running on a cluster

scQC's executor submits each task as its own job where a scheduler is configured, and each job
writes its own log. On PBS Pro, `#PBS -o` resolves against the *submitting* host, so a job that
lands on another node can exit 0 having delivered nothing; scQC redirects inside the job script
instead and treats the exit marker in the log as authoritative over `qstat`.

Point the orchestrator at a project on shared storage and let it submit. Two practical notes:

- **Give the orchestrator enough memory.** It reads matrices for some audits.
- **Outputs are content-addressed**, so two runs with different parameters can share a project
  directory without either overwriting the other - but not two runs that differ only in their
  decisions file, which share one directory.

## Troubleshooting

| symptom | cause |
|---|---|
| report says NOT DETERMINED | no gate was evaluated in this run and none was stored — not a pass |
| a figure says NOT PRODUCED | the run did not record the data it needs; the reason names the file or the flag that would produce it. F5 is the expected one: its sweep is opt-in, and `--dbr-sd-sweep default,dbr,1` draws it |
| a task will not re-run | it completed, its signature is unchanged and its outputs exist — delete an output |
| every cell has `pct_counts_mt == 0` | wrong `mt_prefix` for the species; matching is case-sensitive |
| `scqc report` says there is nothing to rebuild | that run finished before `payload.json` was written; re-run the report step |

---

- **[How each filter is calculated](FILTERS.md)**
- **[Output reference](OUTPUTS.md)**
- **[Report design](REPORT_DESIGN.md)**
- **[Principles](PRINCIPLES.md)** — why the gates are shaped the way they are
