# The QC report — SPECIFICATION

> **This document specified a report that was not built. It is built now, and this note is kept
> because the header outlived the fact by several releases.** `report/build.py` writes the HTML
> and the JSON, `report/figures.py` draws all eighteen figures and `report/collect.py` assembles
> their data from a finished run's tables. The freshness check is fed: the run supplies its newest
> input's time (`engine/pipeline.py`, `_newest_input_block`), and a report says `NOT CHECKED` only
> where a time cannot be read or compared. The Status section of
> [README.md](../README.md) remains the inventory.

What the pipeline is to hand a human: one self-contained HTML file, no external requests, openable
from a filesystem five years from now.

## What this borrows from nf-core, and where it deliberately differs

nf-core and MultiQC solved the parts everybody needs: one document per run, every sample in it,
software versions collected automatically, consistent section shapes so a reader who has seen one
report can read any of them. All of that is kept.

The difference is what the document is *for*. A MultiQC report answers **what happened**. A QC
report has to answer **what was decided, on what evidence, and what is still unresolved** —
because QC is not a measurement, it is a sequence of choices about what to discard, and the
choices are what a reviewer, a collaborator, or you-in-a-year actually need.

| | MultiQC / nf-core | here |
|---|---|---|
| organising unit | the tool that ran | the decision that was made |
| a threshold | appears as a value | appears with **who chose it and against what** |
| a warning | a coloured cell | a **verdict** that can block the run |
| provenance | versions | versions + reference + the decision file + input verification |
| what it can't do | not addressed | **a required section in every step** |
| staleness | not addressed | **flagged** on the document's face and in `report.json` — an artifact older than its inputs is wrong, not old. Refusing to publish one is not built (section 4) |

Three things this report does that a metrics report structurally cannot:

- **It distinguishes a number the data produced from a number a person chose.** Both look
  identical on a page and they carry completely different weight.
- **It states what each step cannot establish.** The omission is what makes an honest report read
  like an over-confident one.
- **It shows never-examined as its own category**, never folded into a negative. A barcode that
  was never scored for doublets is not a barcode that passed.

---

## Layout

The order is the argument: a reader who stops after ten seconds should stop having read the thing
that matters most. What that is turned out not to be the verdict.

The document opens with a **masthead** (cohort, library count, mode, verdict chip) and a
**decision strip**: observations in, kept, removed, and *evenness across the design* — the widest
removal ratio any gate measured between arms of the design. That last figure leads because it is
the one a conventional QC report does not carry, and the one that decides whether a technical
gradient has been laid down where the biology is measured.

Then **the design itself** (the confounding block: F16, F17 and F18), and **where each count floor
came from** (F6 and F13).

Then **what quality control did**: five fixed rows - the UMI floor (F7), the gene floor (F14), the
mitochondrial ceiling (F15), the doublet call (F10) and the survivors on the same embedding (F11) -
each showing the distribution before the cut and what survived it, with the threshold and how it
was arrived at (DERIVED / DECLARED / ADJUDICATED, per library / cohort constant) on a rule down the
right. The page's own sub-heading reads "one row per criterion"; the rows are four of the applied
criteria, not all of them (see "The spine holds five fixed rows" under Figures). A row whose figure
was not produced states the cut anyway, beside a NOT PRODUCED note carrying the reason.

After the spine: the same distributions on a linear axis; what each criterion removed uniquely;
what went in and what was examined; the cluster check; every threshold per library; what was
decided and by whom; what this run could not establish; the findings; provenance; open items; and
the report's own defects.

**What each criterion removed carries F9 and a table, and the table is not a caption for the
figure.** F9 draws the cohort, and the cohort is the wrong unit for the question the section is
asked: a criterion taking 1% of one library and 30% of another (an illustration) sums to an unremarkable cohort bar,
and that is a technical gradient sitting exactly where the biology is measured. So the table is one
row per library — `fired` and `alone` per criterion, the cohort as a marked last row — with the
widest per-criterion ratio across libraries stated beneath it, because that ratio is the finding
and nobody should have to scan the whole table for it. `fired` and `alone` are both shown: printing
only the first makes every criterion look load-bearing, only the second makes a criterion that
agrees with its neighbours look inert.

**A figure is never omitted because it is empty.** An absent block reads exactly like a block
nobody needed, so a figure that was not produced is a bordered NOT PRODUCED note carrying the
reason. Sections are not all held to that: the per-library thresholds and the parameter table are
left out of the page when the run supplied nothing for them, and so is the confounding block when
the design was not stated; the absence is then a row of the defects section at the foot.

The sections below specify the content each block carries.

### 1 · Verdict

The deliverable in one line, then every REFUSE and REVIEW raised by any gate, with the step that
raised it. This block is `verdict.lines` in `report.json`. On the page the masthead's chip carries
the run's overall verdict and the decision strip the deliverable's counts, so a reader meets a
REFUSE or a REVIEW before the figures; the findings themselves, each with its step, are the
Findings section, which sits below the figures and the limits.

```
(illustrative numbers)
DELIVERABLE   N of M called cells   (x% removed)
REFUSE        none
REVIEW        step 2  one library lost y% of aligner calls to the denoiser
              step 7  n nuclei survive inside a cluster that is mostly doublets
OPEN          tier-2 tolerance UNMEASURED - run-to-run spread never quantified
```

A run with a REFUSE still produces a report. Refusing to write the document would destroy the
only record of why the run stopped.

### 2 · What was decided, and by whom

One row per parameter the run records. **The class column is the point of this report.**
(Illustrative values; the rows are the ones `engine/pipeline.py`, `_parameters_block`, writes.)

| parameter | value | class | basis |
|---|---|---|---|
| cell caller | the denoiser's output on nuclei; the aligner's call on single cells | FIXED | pipeline contract |
| variable-gene selection | every gene, no class excluded | FIXED | no exclusion term exists |
| removal point | step 7 only | FIXED | no other step drops an observation |
| clustering resolution | 2.0 | FIXED | the assay's entry, where `--resolution` was not given |
| each declared flag that was given: `--dbr`, `--dbr-sd`, `--dbr-sd-sweep`, `--light-floor`, `--resolution`, `--extra-resolutions`, `--seed`, `--cpu` | as given | DECLARED | the flag that supplied it |
| UMI floor | *u* | DERIVED | step 5: the density valley per library, one cohort constant, bounded (DECLARED on single cells) |
| gene floor | *g* | DERIVED | the same, on the gene axis: on nuclei bounded to [100, 600] |
| mitochondrial ceiling | *m*% per library | DERIVED | step 5: median + k × 1.4826 × MAD, bounded (DECLARED where the bound decides most libraries) |
| any of those three, overridden in `decisions.yml` | the override | **ADJUDICATED** | declared in `decisions.yml`; the row carries the approver and their words |

The ambient learning rate and the cluster check's A/B/C thresholds have no row in this table
today, and `dbr.sd` is a DECLARED row, not a derived one.

- **FIXED** — the pipeline's own contract; changing it is changing the pipeline.
- **DERIVED** — computed from this dataset. Reproducible from the data alone.
- **DECLARED** — supplied by the operator up front, before seeing the result.
- **ADJUDICATED** — a human decided, after seeing the result, and the row carries their own words.

An ADJUDICATED row without verbatim text is a bug, and the row is withheld: its name stays, its cells are empty, and the defects section says why.

### 3 · The steps

One record per step in `report.json` (`steps`), each with the same four parts:

1. **What this step does** (`purpose`).
2. **The figures** (`figures`).
3. **What it found** (`found`) — numbers, with the files they came from (`sources`).
4. **What this step cannot establish** (`cannot_establish`).

The page is not laid out by step (see Layout): it places the figures by question and prints the
fourth part of every step in one section, "What this run could not establish". That fourth part
is required. A step that ran and supplied none carries
`cannot establish: NOT STATED — this is a defect in the report, not an absence of limits.`

### 4 · Provenance

```
pipeline        scQC <tag>  commit <sha>   (dirty: no)
invocation      the exact command line
decisions       decisions.yml  sha256 <...>
reference       registry name, genome build, annotation, intron handling
tools           every executable, its version, and how the version was obtained
input check     P1 raw VALUES: PASS    P2 raw DROPLETS: PASS
generated       <ISO timestamp>   newest input: <ISO timestamp>
```

These lines are `provenance.lines` in `report.json`, with the project, the runtime and a freshness
line beside them. The page's Provenance section prints less today: the invocation, the generation
time with its source, and the wall clock. The commit, the dirty flag, the decisions hash, the
reference, the tools and the input check are in the JSON only.

The last line is the freshness contract, and it is machine-checked: `generated` earlier than
`newest input` is a row of the page's defects section and counted as a defect. The **refusal to publish**
this specification asked for is not built — nothing in scQC publishes a report — and it stays the
most load-bearing open item, because staleness has no symptom and nothing else will catch it.

### 5 · Open items

Named, each with what would close it and who it is blocked on (`open_items.items` in
`report.json`; `NOT STATED` where the run did not say). A run names the steps that were blocked,
and a run with none says `none recorded` explicitly, because a missing section and an empty one
are not the same claim. The page's Open items section does not print the items today: it prints
the block's two field names (`stated`, `items`), so read them in the JSON, or in `verdict.lines`.

---

## Figures

Eighteen, each answering a single question a reader would otherwise have to ask. Fifteen follow
the steps; F10/F11, F7/F12 and F6/F13 are pairs — the same treatment applied twice, because a
comparison is what carries the meaning. F16, F17 and F18 belong to no step: they describe the
design of the experiment rather than the behaviour of the pipeline.

**F10 and F11 must share one embedding.** Re-embedding the retained nuclei changes the layout, and
a reader comparing the two panels would be looking at a difference that may be the projection
rather than the data. `fig_f10_umap_per_library` therefore takes coordinates and never computes
them.

**The embedding is built over every called barcode - the denoiser's call on nuclei, the aligner's
on single cells - and not over the population step 6 clusters.** The clustering must run on the cells that reach the deliverable, or its flags
describe empty droplets. The embedding must contain the cells a criterion REMOVED, or F11 is drawn
over a population every one of them has already left, and it can only ever answer "there were
none". Those are different populations, so they are two passes, and step 6 does both.

**A flag with three values gets three colours.** A barcode below the light floor was never handed
to the doublet detector: it is UNKNOWN, not a singlet. F10 draws it in its own colour and counts
it out of the percentage, because the alternative is a "not a doublet" cloud that is a quarter
unexamined with nothing on the page saying so.

**A capped axis states what it hides.** F12 clips a long tail so the bulk is legible, and prints
the share of nuclei that fall outside. A truncated axis that does not say it is truncated reports
a distribution nobody drew.

**The report names what was observed by the declared assay.** Every caption, axis label, step text
and figure question below that says "nucleus" says "cell" on a single-cell cohort: `{one}` and
`{many}` in the table are filled from the samplesheet's `assay` (`report/units.py`), and a cohort
with no single assay is named in barcodes, never in nuclei by default.

| # | step | question it answers | form |
|---|---|---|---|
| F1 | 0 ingest | is this really raw, unfiltered input? | barcode-rank per library, log–log, the aligner's cut marked |
| F2 | 1 ambient | how much was removed, from what, and evenly? | removal per library + per gene; **removal rate by design level** |
| F3 | 2 cell call | did the denoiser drop cells the aligner kept? | aligner vs denoiser per library, gate lines drawn |
| F4 | 3 floor | what was never examined? | stacked bar: scored / below floor / **never scored** |
| F5 | 4 doublets | is the rate a measurement or the prior? | sweep curve per library (a published band is shaded only if one is passed; a run passes none) |
| F6 | 5 quality | where is the UMI cut and why there? | per-library density, valley marked, cut drawn |
| F7 | 5 quality | what did the cut change? | before/after violins on a log axis (F12 is the linear one) |
| F8 | 6 clusters | are any clusters technical? | depth × mito scatter, flags coloured, thresholds drawn |
| F9 | 7 apply | what did each criterion remove *uniquely*? | unique vs shared contribution per criterion |
| F10 | 4 doublets | where in the manifold did the doublets sit? | one embedding per library, doublet calls coloured |
| F11 | 7 apply | did the removed {many} leave as a population, or scattered? | **the same embeddings**, removed {many} coloured |
| F12 | 5 quality | the same count distributions, on the scale people work in | F7's data on a linear axis, with the share above the cap stated |
| F13 | 5 quality | where is the GENE cut and why there? - step 5 derives two floors and applies both | F6's form on the gene axis, with the gene floor drawn |
| F14 | 5 quality | what did the GENE floor change? | F7's form on the gene axis |
| F15 | 5 quality | what did the MITOCHONDRIAL ceiling change, library by library? | F7's form on the mitochondrial axis, over the barcodes above the light floor, with each library's OWN ceiling drawn as a segment |
| F16 | the design | what is confounded with what, before any number is read? | the samplesheet as a grid, beside every pair of factors classified `aliased` / `nested` / `crossed` |
| F17 | the design | how far apart do the confounded arms sit on each QC metric — further than the libraries inside an arm? | one panel per QC metric: a violin per arm over the kept nuclei, each library's own median drawn on it, and the between-arm gap against the widest within-arm spread |
| F18 | the design | did the filter remove the same share from each confounded arm? | removal rate per arm, in total and per criterion |

**F2 and F6 are the two that matter most.** F2 is the only figure that can show a technical
removal has become an apparent biological difference. F6 is the only one that shows a derived
threshold sitting where the data actually separates — or not.

### The confounding block — F16, F17, F18

These three figures describe the design of the experiment, not the behaviour of the pipeline. They
are computed from the samplesheet and the per-cell tables, and no step of the pipeline can change
what they report.

**F16 — which factors cannot be separated.** For every pair of design factors, F16 compares the
partition each induces over the libraries and classifies the pair:

| | |
|---|---|
| `aliased` | each factor is fixed once the other is known — the partitions are identical |
| `nested` | one factor is fixed by the other, but not the reverse |
| `crossed` | both levels of each occur with both levels of the other |

The comparison is exact set arithmetic. There is no statistic and no p-value: a confounded factor
cannot be tested, which is why it is reported here instead. The panel beside the matrix redraws
the samplesheet as a grid, so the classification can be checked by eye.

**Arms.** Where factors are aliased, the arm is the partition they share, and its label carries
every aliased factor's level — `arm-1 / chem-v2 / batch-A`. F17 and F18 compare arms.

**F17 — every QC metric across the arms.** One panel per metric, drawn over the nuclei the filter
kept. Each panel shows:

- a violin per arm, over every nucleus in that arm (capped at 20,000 per arm, at an even stride);
- one marker per library at that library's own median, in the same colour, shape and position in
  every panel;
- the span those library medians cover, shaded;
- the arm's median of them.

The heading states the gap between the arm medians, the widest spread of library medians within a
single arm, and their ratio. A gap no larger than the within-arm spread indicates a library effect
rather than an arm effect. The ratio is an estimate, not a test.

Metrics drawn: `total_counts`, `n_genes`, complexity (`log10 genes / log10 UMI`), `pct_counts_mt`,
`pct_counts_ribo`, `nuclear_fraction`, `doublet_score`. A metric the per-cell tables do not carry
is drawn as a stated absence naming the missing column.

**Gene-class percentages name their pattern.** `pct_counts_mt` and `pct_counts_ribo` are decided by
a samplesheet pattern, and the axis prints it: `ribosomal % per nucleus (pattern ^Rp[sl])`. A
percentage whose gene class is not on the page cannot be checked, and a ribosomal pattern has
matched a ribosomal protein kinase before (PRINCIPLES.md section 1).

**F18 — what the filter removed from each arm.** The removal rate per arm, in total and per
criterion, over every barcode each arm held. The ratio between the extreme arms is printed for the
filter as a whole and for each criterion separately. This is the part of the confounding the
pipeline itself introduces.

**A design with no aliased pair produces a stated absence.** F16 is still drawn; F17 and F18 report
that there is no arm to compare, rather than constructing a contrast the design does not support.

**A per-library threshold is never drawn as one line.** F15's ceiling is one rule per library,
and a single line across the panel would assert a cohort constant that was never applied -
with nothing on the page to say it is not one. Each library carries its own segment, and the
panel states that the value is per library and gives its range.

**The spine holds five fixed rows, and they are not every applied criterion.** The report's
decision spine names the count floor (F7), the gene floor (F14), the mitochondrial ceiling (F15)
and the doublet call (F10), then shows the survivors on the same embedding (F11): four of the six
criteria step 7 applies. The cell-call criterion (`fail_not_cellbender_cell`, or
`fail_not_aligner_cell` on single cells) and the joint mitochondrial / nuclear-fraction criterion
(`fail_mito_nf`) have no row; what they removed is in the section "What each criterion removed"
(F9 and its per-library table). The gene-floor and ceiling rows read "no figure is produced for
this axis" until F14 and F15 existed: a report showing a third of the filter while looking complete.

**F13 exists because step 5 applies TWO floors.** `cut` is one cohort constant per figure and the
two floors are in different units, so a single density figure can carry only one of them — and the
UMI floor drawn over a gene density would be a line in the wrong units on every panel. A report
that showed the derivation of one applied threshold and not the other would be showing half the
filter while looking complete.

### Rules the figures are built to

Each of these exists because its absence produced a wrong reading at least once. They are rules for
whoever draws a figure, kept by reading: no test enforces them, and where the eighteen figures
depart from one it is said here.

- **No legend states a number the code did not compute.** The drawing functions take their cuts,
  counts and denominators as arguments.
- **Never-examined is its own colour.** Never merged into a negative.
- **A rate carries its denominator**, in the axis label or the annotation; F10 and F11 state
  theirs in the figure's subtitle.
- **Same scale across libraries** in any panel meant for comparison. A per-panel scale makes a
  4× difference look identical to a 4% one.
- **n on the panel**, because a figure with no n reads as universal whether or not that was
  meant: F1 to F6, F8, F9, F13, F16's map, F17 and F18's first panel print one. The before-and-after
  violins (F7, F12, F14, F15) print each library's median and the cut, and the embeddings (F10,
  F11) the share highlighted and the count never determined; neither prints an n.
- **A threshold is drawn, never described.** If a cut is at a value the figure was given, there is
  a line at that value.
- **Colour is not the only encoding** where a flag marks few objects: F8's four states have four
  markers and F4's never-scored bar a hatch. F10 and F11 separate their three states by colour.
- **Log axes labelled in original units**, with real ticks.

---

## Reproducibility

Beyond version capture, four properties nf-core does not currently give you. The first, second and
fourth are implemented; the third is computed and flagged, not enforced:

1. **Decisions are an input, not an argument.** `decisions.yml` is hashed, and its
   hash is in `report.json` (`provenance.decisions`). Two runs with the same data and the same
   decisions file are the same run; a command line does not have that property.
2. **The commit is in `report.json`.** A report from a dirty tree says `dirty: yes` there
   (`provenance.pipeline`); the page's Provenance section prints neither today.
3. **Freshness is enforced, not advised.** Stale artifacts have no symptom — they open, render
   and read exactly like correct ones. So this is machine-checked; blocking is not built.
4. **Removal is recoverable** — *implemented*. Every removed observation keeps the criteria that
   removed it, written beside the deliverable as CSV by `modules/07_apply/apply.py`, so any
   question of the form "what if we had not dropped those" is answerable without a re-run. The
   report has only to quote the file; the record does not depend on the report existing.

## Format

Single-file HTML: inline CSS, no script, data URIs for figures, no external requests. It must
open with no network and no server. A machine-readable `report.json` carrying every number in the
document is written beside it — the HTML is for people, the JSON is for the next pipeline.
