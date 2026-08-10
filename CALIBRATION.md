# The calibration cohort

Every threshold in this pipeline that is not a published default was measured against one real
dataset. This file says what that dataset was, what it can support, and — more importantly — what
it cannot.

## What it was

Ten single-nuclei RNA-seq libraries of a solid tissue, in a two-factor design with five libraries
per arm of the factor used in the differential checks. Droplet-based 3' chemistry, two chemistry
versions within the cohort, intron-inclusive quantification. 305,425 droplets analysed, 244,968
called cells before QC.

It is **unpublished**, so it is not named here: its identifying details — the study, the sample
names, the species, the tissue, and what its design factor actually was — are absent from this
repository. The libraries are referred to as `ctrl_01…ctrl_05` and `treat_01…treat_05`, and the
design factor is a generic two-level `condition`. That renaming is cosmetic: the measured values
below are the real ones.

One thing that does name a species is `references/_registry/registry.tsv`, and it is worth saying
what it is so it is not read as a leak. That row is a **worked example of a registry entry** — it
exists to show what a complete one has to record: species, build, the command that built it, the
aligner version, the gene count and whether introns are retained. It is not a statement about which
reference this cohort was quantified against, and the index itself is not distributed.

## Why the numbers are here and the identity is not

Removing the measurements along with the identity would have left the test suite asserting against
invented values — and a test with fabricated fixtures passes exactly like a real one, which is the
worst of both worlds. What identifies a study is its name, its samples, its species and tissue and
its design semantics. A cell count attached to `treat_02` identifies nothing.

The design *structure* is kept, as a generic two-level factor, because two of the pipeline's gates
exist specifically to catch a loss that falls on one arm of a design and cannot be exercised by a
fixture that has no arms.

## What was measured, and how stable it was

| parameter | value here | varied across libraries | class |
|---|---|---|---|
| ambient learning rate | 5e-5 | one cohort-wide decision | DERIVED per cohort |
| UMI floor | 350 | **274–473** per library | DERIVED per dataset |
| gene floor | 250 | bimodal in every library | DERIVED per dataset |
| mito ceiling | **per library, 6.12–25.00%** | **4× spread in Q3** | **DERIVED per library**; bound DECLARED |
| doublet `dbr.sd` | 0.06 | default was prior-driven | DERIVED per dataset |
| cluster A / B / C | 0.5× / 15% / 50% | proposed against distributions | **ADJUDICATED** |
| cell-call REVIEW / REFUSE | 5% / 10% | one library at 5.0% | DECLARED |
| design differential refusal | 3× | observed 1.04–1.08× | DECLARED |

The column that matters is the last one. **Only the DECLARED rows are safe to carry to another
dataset unchanged.** Everything marked DERIVED must be re-derived, and the per-library spread in
column three is why: a UMI floor that ranged 274–473 *within a single cohort* has no business
being a constant across cohorts.

### The mitochondrial ceiling changed class, and it is worth saying why

It was **ADJUDICATED — "not derivable, no valley"** until the calibration cohort was measured
against its own delivered object. The fixed 40% ceiling in use there removed **zero** nuclei: the
count floors had already taken every nucleus it would have taken. A ceiling that removes nothing
is not permissive, it is absent, and it had been reported as protection.

"No valley" was true and the conclusion drawn from it was wrong. The valley route does not apply
to a unimodal distribution — but that rules out *one* derivation, not all of them. Each library's
own upper Tukey fence, `Q3 + 1.5 × IQR`, needs no valley and no free parameter.

What remains adjudicated is narrower and more honest: whether a mitochondria-high **population**
is damage or a mitochondria-rich cell type. That needs an identity, so the pipeline emits
cluster-level medians and stops.

Two things this row should teach a reader of another cohort:

- **The 6.12–25.00% range is not transferable.** It is what ten libraries of one cohort produced.
  The **procedure** ships; the numbers do not.
- **The bound is the only DECLARED part**, and it must be re-declared in the analyst's own words.
  On this cohort it bound in 3 of 10 libraries — a guard rail. `derive_mito_ceiling` refuses if it
  binds in most, because a declared number that overrides most libraries has become the threshold
  while still being reported as derived.

## What one cohort cannot establish

Stated plainly, because this is the section that is usually omitted:

- **Whether any of it generalises.** n = 1 cohort, one tissue, one species, one platform. Every
  number above is an existence proof that the method works somewhere, not evidence of a range.
- **The tier-2 tolerance.** The stochastic steps — ambient correction, doublet scoring — have a
  run-to-run spread that has **never been measured**. Until a repeat run exists, the pipeline
  records that tolerance as `UNMEASURED` and will not assert one.
- **Whether the gate thresholds are correctly placed.** The 5% / 10% cell-call lines and the 3×
  design-differential line were set against distributions from this cohort and one accepted
  outcome. A threshold that has never fired in anger is a hypothesis.
- **Anything about single-cell (as opposed to single-nuclei) data.** The ambient step is mandatory
  for nuclei and optional for cells; that branch has not been exercised on real cell data.
- **Droplet-free or plate-based protocols.** Not supported and not tested.

## How to calibrate on your own data

Run `--mode evidence`. It derives every DERIVED parameter, refuses to apply anything, and writes
the report described in `docs/REPORT_DESIGN.md`. Read it, decide the ADJUDICATED rows, write them
into `decisions.yml` with your own words, then run `--mode apply`.

The two modes exist so that the act of looking at the data and the act of cutting it are separated
in time and recorded separately. That separation is the whole design.

**That is the design and not yet the software.** At `0.0.1-dev` there is no `--mode` driver and no
report writer, and nothing runs the steps in sequence. What exists is one `scqc` subcommand per
gate, judging tables you produce yourself — `scqc quality --valleys … --metric umi` will propose
the count floor above, or refuse to. Until the driver exists, calibrating on your own data means
running those steps in order by hand. The Status table in [README.md](README.md) says exactly what
is missing.
