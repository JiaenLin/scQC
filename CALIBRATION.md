# The calibration cohort

The nucleus thresholds in this pipeline that are not published defaults were first measured against
one real dataset; the single-cell count bounds and floors are declared, not measured
(`docs/FILTERS.md`). This file says what kind of dataset that was, what was calibrated on it and how, what
it can support, and - more importantly - what it cannot.

## What it was

Single-nuclei RNA-seq libraries of a solid tissue, in a two-factor design with replicate libraries
in each arm of the factor used in the differential checks. Droplet-based 3' chemistry, two
chemistry versions within the cohort, intron-inclusive quantification.

It is **unpublished**, so it is not named here: its study, sample names, species, tissue and design
levels are absent. Its measurements are that cohort's own and are recorded in the cohort's project
records (its stage decisions and run logs), outside this repository. The tool's rules - the
procedures, bounds and lines - do not depend on them, and the tool reads none of them when it runs.
Some sentences here and in the code still say, in words or as a rounded illustration, what those
measurements showed, to explain why a rule has its shape; none changes a computed value. A
generated word list catches the names and the numbers it
holds; a measured number no list holds is caught only by reading, so every number this repository
keeps should say what it is: a default, a bound with its reason, a published value, or synthetic.

One thing that does name a species is the worked registry row in `docs/TOOLS_AND_REFERENCES.md`,
and it is worth saying what it is so it is not read as a leak. It is an **example of a registry
entry** - it shows what a complete one has to record: species, build, the command that built it,
the aligner version and, in its notes, the gene count and whether introns are retained. No registry ships with scQC:
each project passes its own with `--registry`, and the tool falls back to none.

## Why the numbers are not here

Until 2026-11-11 this file and the test suite carried the cohort's real measurements - its cell
and droplet counts, per-library ceilings, rates and spreads - on the argument that a cell count
attached to an anonymous library identifies nothing. It does: a generated word list made from the
cohort's own records finds those numbers, and a value that belongs to one dataset is exactly what
a portable tool must not carry (harness autodev W048). The measurements are the calibration
cohort's and live in its records; this file keeps what was calibrated, how, and what the
measurement showed, in words. The tests assert against synthetic values whose expected results
are computed, not against one cohort's numbers; the last test still running on the cohort's own
valleys (`tests/test_quality.py`) moved to synthetic ones in W048's fix round. One suite still
reads a cohort, and only one you hold: `tests/test_audit_ambient.py` audits the directory
`COHORT_DIR` names, and skips without it.

The design *structure* is kept, as a generic two-level factor, because two of the pipeline's gates
exist specifically to catch a loss that falls on one arm of a design and cannot be exercised by a
fixture that has no arms.

## What was calibrated, and how stable it was

| parameter | how it is set | how it varied across the cohort's libraries | class |
|---|---|---|---|
| ambient learning rate | one cohort-wide decision, a person's, from the per-library fits the tool compares | - | decided per cohort; the tool assesses the fits, it derives no rate |
| UMI floor | the density valley (nuclei; declared on cells) | substantially, per library | DERIVED per dataset |
| gene floor | the density valley (nuclei; declared on cells) | bimodal in every library | DERIVED per dataset |
| mito ceiling | **per library**, inside a declared bound | **severalfold in Q3** | **DERIVED per library**; bound and, on snRNA, `k` both DECLARED |
| doublet `dbr.sd` | declared (`--dbr-sd`, or `dbr_sd` per library); the sweep (`--dbr-sd-sweep`) shows what each setting would call and applies nothing | the default was prior-driven | DECLARED per run: it describes this cohort's loading and is re-declared, not carried |
| cluster A / B / C | B and C proposed from the cohort's own cluster distributions; A a convention (0.5) | the values are the cohort's | proposed, not approved; a value in `decisions.yml` replaces it |
| cell-call REVIEW / REFUSE | 5% / 10% | a library on the REVIEW line | DECLARED |
| design differential refusal | 3x | observed far below it | DECLARED |

The column that matters is the last one. **Only the DECLARED rows are safe to carry to another
dataset unchanged.** Everything marked DERIVED must be re-derived, and the per-library spread in
column three is why: a UMI floor that moved that much *within a single cohort* has no business
being a constant across cohorts.

### The mitochondrial ceiling changed class, and it is worth saying why

It was **ADJUDICATED - "not derivable, no valley"** until the calibration cohort was measured
against its own delivered object. The fixed ceiling in use there removed **no** nuclei: the
count floors had already taken every nucleus it would have taken. A ceiling that removes nothing
is not permissive, it is absent, and it had been reported as protection.

"No valley" was true and the conclusion drawn from it was wrong. The valley route does not apply
to a unimodal distribution - but that rules out *one* derivation, not all of them. Each library's
own upper Tukey fence, `Q3 + 1.5 x IQR`, needs no valley and no free parameter.

What remains adjudicated is narrower and more honest: whether a mitochondria-high **population**
is damage or a mitochondria-rich cell type. That needs an identity, so the pipeline emits
cluster-level medians and stops.

Two things this row should teach a reader of another cohort:

- **No ceiling range is transferable.** Each is what the libraries of one cohort produced under
  one rule. The **procedure** ships; the numbers do not.
- **The bound is DECLARED and must be re-declared** in the analyst's own words - and on snRNA so is
  `k`. See the 2026-08-13 note below: what stays adjudicated differs by assay, because on nuclei a
  mitochondria-high population is not a cell type question at all.

### Revised 2026-08-11 - Tukey calibrates, MAD applies, and `k` is derived

The account above describes the fence as `Q3 + 1.5 x IQR` applied directly. That is no longer what
is applied, and the reasoning for the change is worth more than the change itself.

**What was wrong with applying Tukey.** Nothing statistically - it is the better estimator for a
right-skewed tail, and on this cohort it gave the most even *removal* of anything tested
(unbounded). The objection is biological. Tukey adapts without limit to a wide tail, and on one
library it adapted far past any bound. A ceiling that high is a statement about that prep, not
about what a nucleus can be. Where a tissue's dominant cell type is also an unusually
mitochondria-rich and fragile one - common in solid tissue - letting the ceiling track a library's
own tail that far means the filter's severity is set by prep quality, and set loosest exactly where
the population at risk is largest.

**Why the spread is prep, not biology.** The decisive measurement is available because the design
has replicates: the fence varied **more between two animals of the same group than between the
design groups**. A biological replicate cannot differ by more than the design and still have that
variation be the design. The UMI valley is the control and behaves oppositely - it varied less
within a group than between groups - which is why *it* collapses safely to one cohort constant and
the mitochondrial fence does not.

**Tukey is retained as the *calibration instrument* and never applied:**

```
k_i = (Q3_i + 1.5*IQR_i - median_i) / (1.4826 * MAD_i)      per library
k   = round(median(k_i))                                     scRNA: ONE cohort constant, DERIVED,
                                                               then held within 3 to 10; outside,
                                                               the edge applies, said as a REVIEW
                                                             snRNA: computed and REPORTED, not used
raw = median + k * 1.4826 * MAD                              applied
```

On this cohort the per-library `k_i` spread widely. The per-library k is a perfect monotone
function of skew (`IQR / 1.4826*MAD`, which is 1.349 for a normal distribution and ran well above
it here) - which is also why no single k reproduces Tukey, and why a derived cohort value is
rounded to an integer rather than carried to a decimal it cannot support.

> **What this cohort calibrated is no longer what snRNA applies.** On **2026-08-13** the two assays
> were separated: `snrna` moved to a **5-10% bound at a DECLARED k = 3**, `scrna` kept 10-30% with k
> derived. A nucleus contains no mitochondria, so `pct_counts_mt` on that assay measures
> cytoplasmic carry-over rather than the cell - deriving k from the tail would widen the fence for
> exactly the cohorts whose preparations were dirtiest. See `docs/FILTERS.md`. **Everything below
> describes the retired rule and is kept because it is the evidence the change was made against.**

**Under the retired rule** - k = 4, bound 10-25% - the bound bound in a minority of this cohort's
libraries, all lower, so the result classified as `derived` rather than `bound_dominated`. Tukey
flagged **the most skewed library** as the one place the two independent routes disagreed by more
than 1.5x.

**Under the applied rule** - k = 3, bound 5-10%, with the `mt < 50%` derivation cut also in place -
the bound binds in about half of the same libraries, at both ends, and the result is still
classified `derived`, narrowly. The Tukey-implied k stays inside the 2x review
line of the applied 3. Delivered nuclei fall by a few per cent against the unbounded fence.

### The trade this makes, measured, because it is not visible from the number

Threshold evenness and removal evenness pull against each other, monotonically, across every rule
tested on this cohort. Each column ranks the six rules, 1 = the most even (ranks, not values: the
measured values are the cohort's):

| rule | ceiling spread | removal differential (interaction arm) |
|---|---|---|
| one constant for the cohort | 1 (none) | 6 |
| **MAD k=3, bound 5-10 - applied (snRNA)** | **2** | **3** |
| MAD k=4, bound 10-25 - *retired 2026-08-13* | 3 | 5 |
| MAD k=5, bound 10-25 | 4 | 4 |
| Tukey, bound 10-25 | 5 | 2 |
| Tukey, unbounded | 6 | 1 |

The applied rule is more even than the retired one on **both** columns; down the other rows one
column is bought with the other.

The applied row carries the `mt < 50%` derivation cut and the rows below it do not, so it is not
exactly comparable with them; it is shown in the same table because the trade is the point.

Every step toward a more uniform *threshold* costs a less uniform *effect*. The reason is that the
quantity being fenced genuinely varies - Q3 and the IQR themselves spread severalfold across these
libraries - so an estimator whose output spreads less is not being more consistent, it is
under-adapting. The choice of where to sit on that curve is the analyst's, and it is a judgement
about which error is worse for the study at hand. Here it was made deliberately in favour of the
biological ceiling.

**Two costs recorded rather than smoothed over.** The nuclei a ceiling just above the old 5% bound
was removing were indistinguishable from those it kept in most of the affected libraries (nearly
the median depth of retained nuclei) - which is what motivated the floor. But in one the band
*did* look damaged, so the floor admits some questionable nuclei there. And because nuclei above
10% also carry normal depth in this cohort, nothing here establishes that 10% is where to stop -
only that the old floor was too low.

## What one cohort cannot establish

Stated plainly, because this is the section that is usually omitted:

- **Whether any of it generalises.** n = 1 cohort, one tissue, one species, one platform. Every
  result above is an existence proof that the method works somewhere, not evidence of a range.
- **The tier-2 tolerance.** The stochastic steps - ambient correction, doublet scoring - have a
  run-to-run spread that has **never been measured**. Until a repeat run exists, no tolerance
  is asserted: the code holds none and the report prints none.
- **Whether the gate thresholds are correctly placed.** The 5% / 10% cell-call lines and the 3x
  design-differential line were set against distributions from this cohort and one accepted
  outcome. A threshold that has never fired in anger is a hypothesis.
- **Anything about single-cell (as opposed to single-nuclei) data.** CellBender runs on nuclei
  (unless a library declares its correction skipped) and never on cells (`CELLBENDER_RUNS`); the
  single-cell route has run on one real cell cohort only, and the nucleus-calibrated lines it reads
  from per-assay tables (the resolution, the doublet silence line, the cell-call lines, the
  design-differential line) carry a single-cell entry that says it was calibrated on nuclei and is
  unmeasured on cells. Two more were set on nuclei and are single constants on both assays: the
  doublet step's spread lines (`SPREAD_IMPOSED`, `SPREAD_UNSTABLE`), which only REVIEW.
- **Droplet-free or plate-based protocols.** Not supported and not tested.

## How to calibrate on your own data

Run `--mode evidence` to see what the data proposes: it derives every DERIVED parameter, applies
nothing, and writes the report described in `docs/REPORT_DESIGN.md`. Apply mode, the default,
applies the DERIVED and DECLARED values and records them as such. Where you disagree with one, write
an ADJUDICATED override into `decisions.yml` - the value, who approved it and their own words - and
run apply again; the report's parameter table then shows the override as ADJUDICATED (who approved
it, and their words, are in `report.json` and in STATUS.json's `escapes`).

Both modes are built and run end to end, locally or on a scheduler. The per-step subcommands remain
for judging tables you produced elsewhere — `scqc quality --valleys … --metric umi --assay snrna` will
propose the count floor above, or refuse to. The Status section of [README.md](README.md) says what
is implemented.

*(This paragraph said "at `0.0.1-dev` there is no `--mode` driver and no report writer" until
0.3.0. Both had existed since 0.2.0. A document describing software that has moved on is the same
defect this pipeline exists to catch, arriving in prose.)*

## The clustering resolution - measured, 2026-08-13

The applied resolution is **2.0**, and the reason was criterion D: a cluster's doublet frequency
above a threshold. **D was removed on 2026-08-15** (`tests/test_cluster_flags.py` keeps it out), so
step 6 flags no doublet pocket today: the cluster's doublet percentage is reported, and no criterion
reads it. The sweep is kept as the evidence the value was chosen on, and is still its recorded
source (`cluster_flags.RESOLUTION_SOURCE`). A doublet pocket is small; at a
coarse resolution it dissolves into the lineage cluster around it and its doublet frequency is
diluted to the library average. Swept over five resolutions (1.0, 1.5, 2.0, 2.5, 3.0) on the
calibration cohort, criterion D flagged nothing below 2.0 and the same pockets from 2.0 up, while
the number of clusters grew with every step. The pockets D found at 2.0 are mostly doublets, deep,
and mitochondria-poor - the heterotypic signature, invisible one step coarser. Nothing else in the
check changed with resolution in a way that argued for a different default.

**What the sweep did NOT support**, recorded because it was proposed and tested:

- **A stopping rule on marker strength.** "Split until a cluster has no gene above AUC 0.6", used
  by two published workflows for this tissue, almost never fires anywhere below resolution 3.0 on
  this cohort. It does not constrain, and it is not adopted.
- **Higher resolution as a way to flag less.** It flags *more*: the share of the deliverable in a
  flagged cluster rose from the coarsest resolution to the finest. Resolution redistributes
  flagging into smaller units; it does not shrink it.

`cluster_profile()` will not rank a cluster under `min_cells_for_markers` (3 cells, the tool's own
line). At 2.0 the smallest cluster in this cohort passes that line narrowly, and a different seed
could put a library under it.
