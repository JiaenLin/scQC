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

### Revised 2026-08-11 — Tukey calibrates, MAD applies, and `k` is derived

The account above describes the fence as `Q3 + 1.5 × IQR` applied directly. That is no longer what
is applied, and the reasoning for the change is worth more than the change itself.

**What was wrong with applying Tukey.** Nothing statistically — it is the better estimator for a
right-skewed tail, and on this cohort it gave the most even *removal* of anything tested (1.09×
across the interaction arm, unbounded). The objection is biological. Tukey adapts without limit to
a wide tail, and on one library it adapted to **25.88%**. A ceiling that high is a statement about
that prep, not about what a nucleus can be. In heart in particular, where cardiomyocytes are the
most mitochondria-rich cell type in the body, letting the ceiling track a library's own tail that
far means the filter's severity is set by prep quality.

**Why the spread is prep, not biology.** The decisive measurement is available because the design
has replicates: the fence varied **more between two mice of the same group (3.87×) than between
the design groups (2.59×)**. A biological replicate cannot differ by more than the design and
still have that variation be the design. The UMI valley is the control and behaves oppositely —
within-group 1.33× against between-group 1.52× — which is why *it* collapses safely to one cohort
constant and the mitochondrial fence does not.

**The design now applied.** Tukey is retained as the *calibration instrument* and never applied:

```
k_i = (Q3_i + 1.5*IQR_i - median_i) / (1.4826 * MAD_i)      per library
k   = round(median(k_i))                                     ONE cohort constant, DERIVED
raw = median + k * 1.4826 * MAD                              applied
```

On this cohort `k_i` ran **3.44 to 6.56** (1.90×) for a median of **4.26**, giving **k = 4**. The
per-library k is a perfect monotone function of skew (`IQR / 1.4826·MAD`, which is 1.349 for a
normal distribution and ran 1.59–2.82 here) — which is also why no single k reproduces Tukey, and
why the cohort value is rounded to an integer rather than carried to a decimal it cannot support.

Applied at bound 10–25%, this cohort gives ceilings **10.00–20.51% (2.05×)**, with the bound
binding in 4 of 10 libraries — all lower, none upper — so the result classifies as `derived`
rather than `bound_dominated`. Tukey then flags **Young1**, the most skewed library, as the one
place the two independent routes disagree by more than 1.5×.

### The trade this makes, measured, because it is not visible from the number

Threshold evenness and removal evenness pull against each other, monotonically, across every rule
tested on this cohort:

| rule | ceiling spread | removal differential (interaction arm) |
|---|---|---|
| cohort constant 12.65% | 1.00× | 2.77× |
| **MAD k=4, bound 10–25 — applied** | **2.05×** | **1.75×** |
| MAD k=5, bound 10–25 | 2.46× | 1.60× |
| Tukey, bound 10–25 | 2.50× | 1.37× |
| Tukey, unbounded | 4.23× | 1.09× |

Every step toward a more uniform *threshold* costs a less uniform *effect*. The reason is that the
quantity being fenced genuinely varies — Q3 spreads 3.91× and IQR 4.75× across these libraries —
so an estimator whose output spreads less is not being more consistent, it is under-adapting. The
choice of where to sit on that curve is the analyst's, and it is a judgement about which error is
worse for the study at hand. Here it was made deliberately in favour of the biological ceiling.

**Two costs recorded rather than smoothed over.** The nuclei a ~6–7% ceiling was removing were
indistinguishable from those it kept in three of four affected libraries (median depth 0.88–0.96×
that of retained nuclei) — which is what motivated the floor. But in the fourth, `Aging_HFD1`, the
band *did* look damaged (0.52×), so the floor admits ~173 questionable nuclei there. And because
nuclei above 10% also carry normal depth in this cohort, nothing here establishes that 10% is
where to stop — only that 6% was too low.

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
