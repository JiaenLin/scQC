# The QC report — SPECIFICATION

> **This document specified a report that was not built. It is built now, and this note is kept
> because the header outlived the fact by several releases.** `report/build.py` writes the HTML
> and the JSON, `report/figures.py` draws all fifteen figures and `report/collect.py` assembles
> their data from a finished run's tables. What is still NOT fed is the freshness check: nothing
> supplies a newest-input time, so every report says `NOT CHECKED` rather than claiming to be
> current — see [KNOWN_ISSUES.md](../KNOWN_ISSUES.md). The Status table in
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
| staleness | not addressed | **refused** — an artifact older than its inputs is wrong, not old |

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

Then **what quality control did**: one row per criterion, each showing the distribution before the
cut and what survived it, with the threshold and how it was arrived at (DERIVED / DECLARED, per
library / cohort constant) on a rule down the right. A criterion with no figure states the cut
anyway and says no figure was produced for that axis.

After the spine: the same distributions on a linear axis; what each criterion removed uniquely;
what went in and what was examined; the cluster check; every threshold per library; what was
decided and by whom; what this run could not establish; the findings; provenance; open items; and
the report's own defects.

**A section is never omitted because it is empty.** An absent block reads exactly like a block
nobody needed, so each states what it could not establish instead.

The sections below specify the content each block carries.

### 1 · Verdict

The deliverable in one line, then every REFUSE and REVIEW raised by any gate, with the step that
raised it. Not an appendix. If a run has an unresolved contradiction, the reader meets it before
the figures.

```
DELIVERABLE   127,050 of 244,968 called cells   (48.1% removed)
REFUSE        none
REVIEW        step 2  one library lost 5.0% of aligner calls to the denoiser
              step 7  41 nuclei survive inside a cluster that is 74.4% doublet
OPEN          tier-2 tolerance UNMEASURED - run-to-run spread never quantified
```

A run with a REFUSE still produces a report. Refusing to write the document would destroy the
only record of why the run stopped.

### 2 · What was decided, and by whom

One row per parameter. **The class column is the point of this report.**

| parameter | value | class | basis |
|---|---|---|---|
| ambient learning rate | 5e-5 | DERIVED | default rejected: cohort assessment halved it |
| cell caller | denoiser output | FIXED | pipeline contract |
| UMI floor | 350 | DERIVED | KDE valley; per-library 274–473, cohort median 346 |
| gene floor | 250 | DERIVED | KDE valley; bounded to [100, 600] |
| mito ceiling | 40% | **ADJUDICATED** | *"apply 40% as the pre QC steps"* — operator, 2026-08-09 |
| doublet `dbr.sd` | 0.06 | DERIVED | sweep; the default was prior-driven |
| cluster A/B/C | 0.5× / 15% / 50% | **ADJUDICATED** | proposed against the distributions, approved |

- **FIXED** — the pipeline's own contract; changing it is changing the pipeline.
- **DERIVED** — computed from this dataset. Reproducible from the data alone.
- **DECLARED** — supplied by the operator up front, before seeing the result.
- **ADJUDICATED** — a human decided, after seeing the result, and the row carries their own words.

An ADJUDICATED row without verbatim text is a bug, and the report says so in place of the row.

### 3 · The steps

One section per step, each with the same four parts, always in this order:

1. **What this step does** — two sentences.
2. **The figure.**
3. **What it found** — numbers, each traceable to a named file.
4. **What this step cannot establish.**

That fourth part is required. A section that omits it renders as
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

The last line is the freshness contract, and it is to be machine-checked: `generated` earlier than
`newest input` must be a **refusal to publish**, not a warning. No such check exists in the code
yet — this is the single most load-bearing unimplemented item in this specification, because
staleness has no symptom and so nothing else will catch it.

### 5 · Open items

Named, each with what would close it and who it is blocked on. A run with no open items says
`none` explicitly, because a missing section and an empty one are not the same claim.

---

## Figures

Fifteen, each answering a single question a reader would otherwise have to ask. Nine follow the
steps; F10/F11, F7/F12 and F6/F13 are pairs — the same treatment applied twice, because a
comparison is what carries the meaning.

**F10 and F11 must share one embedding.** Re-embedding the retained nuclei changes the layout, and
a reader comparing the two panels would be looking at a difference that may be the projection
rather than the data. `fig_f10_umap_per_library` therefore takes coordinates and never computes
them.

**The embedding is built over every barcode the denoiser called, and not over the population step
6 clusters.** The clustering must run on the cells that reach the deliverable, or its flags
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

| # | step | question it answers | form |
|---|---|---|---|
| F1 | 0 ingest | is this really raw, unfiltered input? | barcode-rank per library, log–log, the aligner's cut marked |
| F2 | 1 ambient | how much was removed, from what, and evenly? | removal per library + per gene; **removal rate by design level** |
| F3 | 2 cell call | did the denoiser drop cells the aligner kept? | aligner vs denoiser per library, gate lines drawn |
| F4 | 3 floor | what was never examined? | stacked bar: scored / below floor / **never scored** |
| F5 | 4 doublets | is the rate a measurement or the prior? | sweep curve per library, published band shaded |
| F6 | 5 quality | where is the UMI cut and why there? | per-library density, valley marked, cut drawn |
| F7 | 5 quality | what did the cut change? | before/after violins, log and linear |
| F8 | 6 clusters | are any clusters technical? | depth × mito scatter, flags coloured, thresholds drawn |
| F9 | 7 apply | what did each criterion remove *uniquely*? | unique vs shared contribution per criterion |
| F10 | 4 doublets | where in the manifold did the doublets sit? | one embedding per library, doublet calls coloured |
| F11 | 7 apply | did the removed nuclei leave as a population, or scattered? | **the same embeddings**, removed nuclei coloured |
| F12 | 5 quality | the same count distributions, on the scale people work in | F7's data on a linear axis, with the share above the cap stated |
| F13 | 5 quality | where is the GENE cut and why there? - step 5 derives two floors and applies both | F6's form on the gene axis, with the gene floor drawn |
| F14 | 5 quality | what did the GENE floor change? | F7's form on the gene axis |
| F15 | 5 quality | what did the MITOCHONDRIAL ceiling change, library by library? | F7's form on the mitochondrial axis, over the barcodes above the light floor, with each library's OWN ceiling drawn as a segment |

**F2 and F6 are the two that matter most.** F2 is the only figure that can show a technical
removal has become an apparent biological difference. F6 is the only one that shows a derived
threshold sitting where the data actually separates — or not.

**A per-library threshold is never drawn as one line.** F15's ceiling is ten different rules,
and a single line across the panel would assert a cohort constant that was never applied -
with nothing on the page to say it is not one. Each library carries its own segment, and the
panel states that the value is per library and gives its range.

**The spine lists every applied criterion, so every applied criterion needs a figure.** The
report's decision spine names the count floor, the gene floor and the mitochondrial ceiling.
Two of those three rows read "no figure is produced for this axis" until F14 and F15 existed:
a report showing a third of the filter while looking complete.

**F13 exists because step 5 applies TWO floors.** `cut` is one cohort constant per figure and the
two floors are in different units, so a single density figure can carry only one of them — and the
UMI floor drawn over a gene density would be a line in the wrong units on every panel. A report
that showed the derivation of one applied threshold and not the other would be showing half the
filter while looking complete.

### Rules every figure obeys

Each of these exists because its absence produced a wrong reading at least once.

- **No legend states a number the code did not compute.** Every annotation is passed in, never
  typed.
- **Never-examined is its own colour.** Never merged into a negative.
- **Every rate carries its denominator**, in the axis label or the annotation.
- **Same scale across libraries** in any panel meant for comparison. A per-panel scale makes a
  4× difference look identical to a 4% one.
- **n on every panel.** A figure with no n reads as universal whether or not that was meant.
- **A threshold is drawn, never described.** If a cut is at 350, there is a line at 350.
- **Colour is not the only encoding.** Every flag also has a shape or a label.
- **Log axes labelled in original units**, with real ticks.

---

## Reproducibility

Beyond version capture, four properties nf-core does not currently give you. The first three are
requirements on the report writer and are not implemented; the fourth is:

1. **Decisions are an input, not an argument.** `decisions.yml` is versioned, hashed, and its
   hash is in the report. Two runs with the same data and the same decisions file are the same
   run; a command line does not have that property.
2. **The commit is in the document.** A report from a dirty tree says `dirty: yes` on its face.
3. **Freshness is enforced, not advised.** Stale artifacts have no symptom — they open, render
   and read exactly like correct ones. So this is machine-checked and blocking.
4. **Removal is recoverable** — *implemented*. Every removed observation keeps the criteria that
   removed it, written beside the deliverable as CSV by `modules/07_apply/apply.py`, so any
   question of the form "what if we had not dropped those" is answerable without a re-run. The
   report has only to quote the file; the record does not depend on the report existing.

## Format

Single-file HTML: inline CSS, inline JS, data URIs for figures, no external requests. It must
open with no network and no server. A machine-readable `report.json` carrying every number in the
document is written beside it — the HTML is for people, the JSON is for the next pipeline.
