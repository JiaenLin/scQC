# How each filter is calculated

One section per filter. Each gives the exact procedure, the population the threshold is
**derived over**, the population it is **applied to**, what it can establish, and what it cannot.

Two conventions run through all of them and are worth stating once.

**Derived over ≠ applied to.** A threshold tuned on observations it can never act on is set by
data that will be discarded anyway. Where the two populations differ, the section says so and says
why — the difference is deliberate in every case below, and the reason is never the same twice.

**Unknown is never a value.** A quantity that could not be computed is emitted as `None`, never as
`0`, `False` or `NaN`. Every gate is three-valued, and a blank must not read as a pass. This
matters most where it is least visible: `NaN >= 50` is `False` and nothing objects, so a missing
input would otherwise turn into a silent clearance.

---

## Table of contents

| Filter | Step | Threshold is | Scope |
|---|---|---|---|
| [Light floor](#light-floor) | 3 | declared, defaults to 200 UMI | cohort |
| [UMI floor](#umi-floor-and-gene-floor) | 5 | derived from a density valley | one cohort constant |
| [Gene floor](#umi-floor-and-gene-floor) | 5 | derived from a density valley | one cohort constant |
| [Mitochondrial ceiling](#mitochondrial-ceiling) | 5 | derived per library, bounded | per library |
| [Doublets](#doublets) | 4 | the detector's call | per library |
| [Cluster check](#cluster-check) | 6 | proposed from the cohort | per cluster, flags only |

---

## Light floor

**What it is.** The minimum UMI count a barcode must have before a doublet detector is asked to
score it. It is not a quality filter and it removes nothing from the deliverable.

**How it is set.** Declared, via `--light-floor`, defaulting to **200**. That default is
scDblFinder's own documented value and is the only one of the supported detectors to publish one;
for the others the pipeline records it as *assumed*, not documented, so the distinction survives
into the report.

**Why a floor exists at all.** Doublet detectors simulate artificial doublets by summing pairs of
observed barcodes. At very low counts that simulation is miscalibrated — a 30-UMI droplet is not a
cell whose profile can be meaningfully combined with another — and the resulting scores are noise
presented in the same units as signal.

**What happens to barcodes below it.** They are **never scored**, and are reported as *unknown*,
never as singlets. This is load-bearing: folding an unexamined barcode into the denominator as a
negative deflates every doublet rate the pipeline reports.

**A refusal you may meet.** If the light floor is not strictly below the applied UMI floor, the
run stops. Otherwise barcodes between the two reach the deliverable having never been examined for
doublets, and in the ledger that is indistinguishable from having been examined and found clean.

---

## UMI floor and gene floor

Both are derived the same way, on different columns: `total_counts` for the UMI floor,
`n_genes_by_counts` for the gene floor.

### The procedure, exactly

1. Take the metric for every barcode in the denoised object, keeping the **positive** values.
2. Fit a Gaussian kernel density estimate (`scipy.stats.gaussian_kde`, Scott's bandwidth) to
   **log10** of those values. Log10 because the two populations being separated — near-empty
   droplets and cells — differ by orders of magnitude, and on a linear axis the lower one occupies
   a handful of pixels.
3. Evaluate on a linear grid of **512** points spanning the observed range.
4. Find the interior local maxima. The valley is the **minimum between the two tallest modes**.
5. The floor is `10 ** grid[valley]`, in the metric's own units. Its precision is the grid
   spacing, reported as `grid_step_umi_at_valley` — a floor is not more precise than the grid it
   was found on.

The reported density is a density with respect to log10(metric). A figure drawing it must label
the axis that way, or the area under the curve means nothing.

### The bimodality test

Six conditions, **all** of which must hold, each reported individually so a refusal can be read
rather than guessed at. These are shape thresholds, not hypothesis tests; none of the numbers is a
p-value.

| # | Condition | Default |
|---|---|---|
| 1 | two interior modes exist at the reference bandwidth | — |
| 2 | a grid point lies strictly between them | — |
| 3 | both modes carry mass — the smaller side holds ≥ 0.1% of observations *and* ≥ 8 of them | `min_mode_mass`, `min_mode_observations` |
| 4 | the valley is deep enough relative to the lower peak | `min_valley_depth = 0.10` |
| 5 | the modes are separated on the log axis | `min_mode_separation_log10 = 0.30` |
| 6 | the answer is stable under a bandwidth change | `bw_stability_factor = 1.5` |

### What a failed test does — and does not — do

A shallow valley **changes what the floor is called, not whether you get one.**

If a library's minimum exists but is a shoulder rather than a dip, the floor is classified
`declared_informed` — a judgement the density did not contradict — and that library is **named in
the proposal**. It is not silently promoted to a measurement, and it is not thrown away either:
the depth test answers *"is this number a measurement or a judgement?"*, which is a question about
provenance, not about whether the number is usable.

**If no library in the cohort is bimodal, the derivation refuses.** There is then nothing to take
a median *of*, and the result would be a judgement wearing a derivation's clothes. Choose the cut
explicitly and record it as one.

### From per-library valleys to one cohort constant

The per-library valleys are combined into a single constant. A per-library count floor would make
the filter a technical property that varies across the design — which is the failure the
design-differential checks exist to catch — so one constant is preferred even where it fits some
libraries better than others.

The proposal is refused if it falls outside the bounds:

| Metric | Bounds |
|---|---|
| UMI | **200 – 1,000** |
| genes | **100 – 600** |

The bounds are the real guard. On one cohort, two libraries' minima wandered to ~1,040 UMI under a
narrower kernel; the upper bound rejects that on sight, and no depth test is needed to catch it.
Where per-library valleys differ by more than **2.0×**, the spread is reported for review: a
single constant fits some libraries better than others and the poorest-fit library should be
looked at.

**Derived over:** every barcode in the denoised object, with no prior floor. A valley is a
boundary *between two modes* and needs both present. Pre-cutting at any count floor deletes most
of the debris mode — the mode the boundary is defined against — and the minimum then lands inside
the nucleus mode or stops existing altogether.

**Applied to:** every barcode, in step 7.

**Cannot establish:** whether the cut is *right*. It establishes that the distribution has two
modes and that the cut sits between them and inside plausible bounds. A tight real population can
fail the dispersion test, and a large enough artifact can pass it.

---

## Mitochondrial ceiling

**What it is.** An upper bound on each barcode's mitochondrial fraction, `pct_counts_mt` —
the percentage of a barcode's counts falling in genes matched by the declared `mt_prefix`.

### The procedure, exactly

Per library, over the population defined below:

```
# 1. calibrate: what multiple of the scaled MAD does Tukey's fence sit at, per library?
IQR   = Q3 - Q1
tukey = Q3 + 1.5 x IQR                       # the calibration instrument, never applied
k_i   = (tukey - median) / (1.4826 x MAD)
k     = round(median(k_i))                   # ONE cohort constant, DERIVED

# 2. apply
raw     = median + k x 1.4826 x MAD          # the applied fence
ceiling = min(max(raw, lo), hi)              # clamped to the assay's bounds
```

**Two estimators, two jobs.** Tukey decides *how far out* a fence should sit — it is anchored at
Q3 and adapts to a right-skewed tail, which is what makes it a good calibrator. The MAD fence is
what gets applied, because it is anchored at the median and compresses the extremes, so no single
library's ceiling runs away from the cohort. Tukey is then retained per library as an
**independent cross-check** and is never applied.

`k` is **derived, not declared** — a property of the cohort's tail shape rather than a number
someone chose. It is rounded to an integer because the per-library `k_i` typically span far more
than a decimal's worth: on the calibration cohort they ran **3.44 to 6.56** (1.90×) for a median
of 4.26, applied as **k = 4**. `MAD_K_BOUNDS = (2, 10)` are sanity limits; a cohort landing outside
them is reported, not silently clamped to the edge. A library with `MAD = 0` is excluded from the
calibration, because its implied k is undefined and a library that cannot contribute must not be
able to change the answer by being counted as a default.

Quartiles and the MAD are linearly interpolated. The bounds are a **declared statement about what
a nucleus can be**, not a measurement:

| Assay | Bounds |
|---|---|
| `snrna` | 10 – 25% |
| `scrna` | 10 – 30% |

An unknown assay with no explicit bounds is refused — guessing a bound for an unknown assay is
guessing what a cell can be. Non-default bounds require `declared_by`: the analyst's own words for
why they are what they are, because the bound is the one part of this the data cannot supply.

> **The snRNA lower bound was 5% until 2026-08-11.** It was raised because 5% permitted an applied
> ceiling to vary **4.08×** across ten libraries of one cohort (6.12% to 25.00%) — one filter doing
> very different things to samples of one experiment. A library whose fence lands at 6% is being
> cut close to its own third quartile, and nothing establishes that a nuclei prep separating at 6%
> differs biologically from one separating at 25% by that factor. At 10–25% the same cohort spreads
> 2.50×.
>
> **The trade this makes**, measured, because it is not obvious: raising the floor acts only on
> libraries with a low fence, so where those correlate with a design arm it makes the *applied
> threshold* more even while making the *removal rate* less even. On that cohort, mitochondrial
> removal by diet moved 1.10× → 1.30×, and within the arm carrying the interaction 0.91× → 0.73×
> — both far below the 3× refusal line — while 1,228 more nuclei of 117,021 were retained.
> Threshold-evenness and removal-evenness are different properties that pull against each other.

**When the bound binds, the per-library derivation was not used.** Each library's `clamped` field
records `upper`, `lower` or blank, and the two directions are counted and named **separately**,
because they are opposite events: an **upper** clamp removes nuclei the library's own fence would
have kept — the only direction that can delete signal — while a **lower** clamp retains nuclei the
fence would have cut.

If the bound decides the ceiling in more than half the cohort, the result is **reclassified rather
than refused**: `provenance` becomes `bound_dominated` instead of `derived`, and a REVIEW note says
so. A ceiling that is really the declared constant for most libraries must not be *described* as
derived — but refusing to run does not fix the classification, it only stops the analysis, and it
fires hardest when an analyst has deliberately narrowed the bound, which is legitimate.

**The spread of the applied ceiling is reported, and reviewed above `CEILING_SPREAD_REVIEW`
(3.0×).** This is not implied by the design differential: the differential asks whether removal
falls evenly across the arms of the design, the spread asks whether one filter is treating samples
of one experiment alike at all. A cohort can pass the first while its ceiling varies fourfold —
which is exactly how a 6.12% library and a 25.00% library came to sit in one deliverable with
neither check objecting.

### Which population the quartiles are taken over

**Derived over:** called cells **at or above the light floor**. Both restrictions matter:

- *Called cells*, because a barcode with no counts has no meaningful percentage — `0/0`.
- *At or above the light floor*, because a percentage needs a denominator large enough to mean
  something. A 30-UMI droplet with 10 mitochondrial counts reads 33%, and a fence built on Q3 is
  set by exactly those barcodes. At ~100 counts the binomial standard error on a fraction near 18%
  is close to four percentage points — noise, sitting precisely in the tail that determines Q3.

Deriving it over every called barcode instead produces ceilings materially higher, and looser on
exactly the cells the ceiling is supposed to police.

**Note the asymmetry with the count floors**, which take *no* floor. The two quantities need
opposite treatment: a mode boundary needs both modes present, a percentage needs a denominator.
Both are computed in the same pass over the same object, and the population each was taken over is
recorded beside it — `pop_floor_umi` and `pop_n_all_called` in the ceiling table. A threshold whose
derivation population is unrecorded cannot be compared with anyone else's.

**Applied to:** every barcode, in step 7, per library.

**Cannot establish:** whether a mitochondria-high *population* is damage or a mitochondria-rich
cell type. The ceiling establishes only that each library was cut at its own outlier fence and that
the result is even across the design. Telling damage from biology needs cell identity, which this
pipeline does not establish.

---

## Doublets

**What it is.** A per-barcode call from an external detector — scDblFinder is the supported one —
run per library, never pooled.

### The procedure

Barcodes at or above the light floor are exported and scored. The detector is given:

| Parameter | Meaning |
|---|---|
| `dbr` | expected doublet rate |
| `dbr_sd` | uncertainty on that rate |
| `seed` | fixed, and recorded |

**`dbr` and `dbr_sd` are declared and have no default.** The expected doublet rate is a property
of how the libraries were loaded and cannot be read off the data; the adapter refuses a missing
one rather than substituting a value nobody chose. They may be given per library — loading
concentration can differ — or once for the cohort.

The rate is reported as **doublets / scored**, and the denominator is stated every time. The same
call over all called cells gives a materially different number, and a rate whose denominator is
unstated is not comparable with a published one.

### The health checks

The point is to establish that the rate is a *measurement* and not the prior handed to the
detector:

| Check | Fires when |
|---|---|
| effectively silent | a library calls < 0.5% |
| rate imposed, not measured | spread across libraries < 1.05× — the fraction is being set, not found |
| unstable | spread > 5× |
| design-differential | worst arm / best arm ≥ 3×, **and** the worst arm calls ≥ 1% |

The materiality condition on the last one is deliberate: a ratio between two near-zero rates is
dominated by single libraries, and is reported rather than refused.

**Derived over / applied to:** barcodes at or above the light floor, per library. Barcodes below
it are unknown, not singlets.

**Cannot establish:** whether a called doublet *is* one. It checks that the rate is a measurement
rather than the prior, and that it does not fall unevenly across the design.

---

## Cluster check

**What it is.** Step 6 clusters each library, profiles every cluster, and **flags**. It removes
nothing. A cluster-level removal deletes a whole population and is the most destructive class of
removal available; it requires an explicit human decision taken elsewhere.

### The four criteria

| | Criterion | Fires when |
|---|---|---|
| **A** | low RNA | cluster median UMI < *f* × that **sample's own** median |
| **B** | high mitochondrial | cluster median `pct_counts_mt` above a threshold |
| **C** | uninformative | share of the cluster's top-20 markers falling in the locked mt+ribo set |
| **D** | doublet | cluster doublet frequency above a threshold |

```
FLAG  = (A and C) or (B and C) or D
WATCH = C alone
```

A alone, B alone, and every continuous value are **reported, not flagged**.

**Why the conjunction.** C alone is not a QC failure — a cluster whose markers are mitochondrial
at normal coverage is the signature of a mitochondria-rich cell type rather than damage. A alone
is low depth with informative markers: a genuine low-RNA population, since some cell types carry
less RNA. The conjunction removes both false-positive classes. Its cost is that the result now
depends on two cut-points at once, so **the continuous values ship beside every boolean** — a
`FLAG` column summarises two thresholds and a reader cannot see from it which one a cluster failed.

**A is sample-relative** because median depth differs between libraries; an absolute floor would
flag whole libraries rather than clusters. Note that criterion A compares against a *fraction*, so
a percentage supplied there would silently never fire.

### Where the thresholds come from

**Proposed from the cohort in hand**, not inherited. B is the p95 of cluster mitochondrial
content. C is the midpoint of the gap between the bulk and the tail *if the distribution is
bimodal*, and the p95 otherwise — with the proposal stating which, because a percentile is a
weaker basis than a valley. Anything declared in `decisions.yml` overrides by name, and the
recorded source says which values came from where.

Thresholds measured on one cohort **do not transfer**. They are chosen against that cohort's own
distributions, and the module proposes rather than remembers.

### Three-valued, deliberately

Markers are computed at the default resolution only, so elsewhere C is **unknown** — and an
unknown must never be written as `False`. `NaN >= 50` is `False` and nothing objects, which would
quietly turn `FLAG` into the full rule at one resolution and "D alone" at every other, under one
name. The flag counts would then differ across a sweep because of what was *calculated*, not
because of anything in the data.

The same holds inside the conjunction: `A and C` with A unknown is **unknown**, not False. Where
markers are absent, C, FLAG and WATCH are all withheld rather than reduced to their D term — D
itself is still reported, so nothing is hidden.

**Cannot establish:** whether a flagged cluster is technical. A cluster flagged for mitochondrial
content cannot be told from a mitochondria-rich cell type without an identity.

---

## What step 7 does with all of them

Step 7 is the only step that removes anything, and it is deliberately three stages:

1. **Measure.** Every criterion is evaluated for every barcode and written to a table — kept and
   removed alike. Nothing is filtered and no object is written.
2. **Record and check.** The ledger names each removed barcode and **all** the criteria that fired
   on it, never just the first. The arithmetic is verified against the mask, and a record that
   disagrees with the mask stops the run rather than reporting whichever number was asked for
   first.
3. **Write.** The kept barcodes are materialised into one object.

An object filtered before its ledger was checked would have performed the removal before anything
verified it. The write stage takes a **keep-list**, not a threshold, so the object cannot be
filtered by anything other than what was recorded.

Every criterion is `True` or `False` for every barcode — never unknown. An observation removed on
a criterion nobody evaluated has not been judged, and the ledger refuses one.

The applied criteria are, in ledger order:

```
fail_not_cellbender_cell   the denoiser left it no counts
fail_umi_floor             below the applied UMI floor
fail_gene_floor            below the applied gene floor
fail_mito_ceiling          above that library's mitochondrial ceiling
fail_doublet               called a doublet by the detector
```

### What apply mode writes

Both shapes, because they answer different questions:

| path | what it is |
|---|---|
| `objects/cohort.deliverable.h5ad` | every library, filtered, in one object |
| `objects/cohort_per_sample/<sample>.filtered.h5ad` | one filtered object per library |
| `tables/removal_ledger.csv` | one row per **removed** barcode, with every criterion that fired on it |

The combined object is what integration and differential testing read. The per-library objects are
what annotation reads — and this pipeline's cluster check is per library precisely because
identity is decided there.

Every retained nucleus carries what step 6 found about its cluster:

| `obs` column | |
|---|---|
| `sample` | the library it came from |
| `cluster` | its cluster in that library's clustering, as a **label** — a cluster id reads as a number and is not one |
| `cluster_FLAG`, `cluster_WATCH` | that cluster's verdicts, as `True` / `False` / **absent** |
| `cluster_pct_doublet`, `cluster_median_pct_mt` | the continuous values behind them |

Absent is a third state and is preserved as such. A barcode step 6 never labelled, or a cluster
missing from the profile, leaves these `None` — not `False` and not `0`. *"This cluster was not
flagged"* and *"this barcode was never examined"* are different facts, and a deliverable that
conflates them claims a cluster check that did not happen.

Carrying the flags forward is not a convenience. Without them the next stage has to re-cluster to
ask a question step 6 already answered — and re-clustering a *filtered* object does not give the
same answer, because criterion D is a tautology once the doublets have been removed.

Each threshold is recorded as **ADJUDICATED** (declared in `decisions.yml`, with an approver and
their own words) or **DERIVED** (what the pipeline measured). A run without a decisions file uses
the derived values and says so; a proposal is never later read as a decision.
