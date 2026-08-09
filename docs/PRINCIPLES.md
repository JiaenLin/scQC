# Principles

Four rules the pipeline enforces in code. They are here because each was learned by getting it
wrong, and because a principle that lives only in a document is a principle that gets skipped.

---

## 1 · The removal checklist

**Before anything is removed, five questions must be answered in writing.** This applies to any
action that removes, excludes, filters, corrects, collapses or transforms cells, genes, samples or
clusters.

The reason is simple and uncomfortable: *a technical decision that deletes signal is
indistinguishable, downstream, from a biological finding that the signal was absent.*

1. **What exactly does this remove?** The actual list — gene symbols, barcode counts, sample names
   — not a description of the category. Print it and read it. A regex intended to match ribosomal
   genes matched a ribosomal protein *kinase*, an mTOR signalling gene, in a metabolic study. It
   was called "ribosomal genes" until somebody printed the list.
2. **Does the removed signal exist anywhere else in the data?** If yes, removal is cheap. If no,
   removal is a decision not to look.
3. **Is the removal differential across the design?** Compute it per group. A filter that removes
   53% of one group and 6% of another has converted a technical property into an apparent
   biological difference, and no downstream analysis can undo that.
4. **Is there a reversible version?** Excluding a gene from variable-gene selection changes only
   the embedding and leaves it testable; excluding it from differential testing does not. Prefer
   the reversible form.
5. **If it must be irreversible, is what was removed still on disk?** Removed barcodes, genes and
   clusters must remain recoverable and be recorded, not merely counted.

**Question 3 is checked in code at three steps.** Each computes a rate per level of every design
factor and records a `REFUSE` finding on a differential, rather than a warning about one:

| module | the rate it groups by design level | bound |
|---|---|---|
| `modules/01_ambient/audit_ambient.py` | fraction of counts removed by the denoiser | ≥3× **and** ≥1% on the worst arm |
| `modules/02_cells/cellcall_gate.py` | aligner-called cells lost at the cell call | ≥3× **and** ≥1% on the worst arm |
| `modules/04_doublets/doublet_health.py` | the rate at which doublets are **called** | ≥3× **and** ≥1% on the worst arm |

The third is not a removal rate. Step 4 removes nothing; it scores and flags. What is checked is
whether the detector *calls* at a different rate on one arm of the design than another, because a
detector that fails on one arm imports that failure into every count downstream of it — and it can
do so while every per-sample number looks unremarkable, which is why the check is on the ratio
across arms rather than on any library's own rate.

All three bound the ratio by materiality, for the reason set out under rule 3, and all three
report rather than refuse when one arm is at zero: a ratio against zero is **undefined**, not
infinite, and manufacturing a huge number out of it refuses runs that removed almost nothing.
A one-sided loss is still raised — as a REVIEW, with the arms printed — because it is exactly
the shape the check exists to catch even when the ratio cannot express it.

**Reporting a refusal and enforcing one are separate calls.** These three modules build findings
and reduce them with `verdict()` to `REFUSE`, `REVIEW` or `PASS`. Returning the string `"REFUSE"`
does not stop anything — the caller does, and in the shipped caller that is the `scqc` subcommand
exiting 2. The modules that stop a run themselves are the ones that would otherwise hand back a
usable value alongside the verdict: a derived threshold (`modules/05_quality/quality.py`) or a
removal (`modules/07_apply/apply.py`) raise, because returning a number and a complaint together
is an invitation to take the number and drop the complaint.

Two corollaries:

- **A published filter's citation establishes that it was reasonable for someone else's question.**
  On your data it is a hypothesis, not a result. Run the five questions on it anyway.
- **Re-audit a classification whenever a new consumer starts using it.** Errors a gene list
  tolerates while it feeds a summary percentage become destructive the moment it starts deciding
  what to discard.

---

## 2 · Only one step removes anything

Steps 0–6 measure, score, flag and refuse. **Step 7 is the only code path permitted to drop an
observation**, it requires a recorded approval containing the operator's own words for that
specific action, and **there is no force flag**. The approval is matched against the action text,
not against a description of it: words spoken about a different action are not consent to this one.

This is why the modules read as if they are over-explaining what they do not do. A reviewer needs
to be able to establish, quickly, that a step cannot silently cost them cells.

Question 5 of the checklist is answered in the same place. Step 7 builds a record pairing every
removed observation with the criteria that fired on it, refuses if that record and the mask
disagree about how many went, and writes it as a CSV readable with nothing installed. A removal
you cannot enumerate afterwards is a removal nobody can question.

---

## 3 · A gate that fires on correct behaviour gets switched off

Every refusal in this pipeline is bounded so it cannot fire on a legitimate run. Three concrete
forms that bound takes:

- **Bound a ratio by materiality.** A 5.42× differential computed from 541 observations against 76
  is arithmetic, not evidence. Gates carry a materiality floor as well as a ratio, so a run is not
  stopped over half a percentage point of compositional distortion.
- **Never block the remedy.** A freshness check that refuses the command that would rebuild the
  stale artifact is worse than no check.
- **Compare approvals on stripped text.** A trailing newline is not a difference of intent, and a
  refusal over one trains the operator to route around the gate. The words must still be the
  operator's own — but whitespace is not consent.

A gate with no escape hatch gets disabled wholesale. Every deliberate bypass in this pipeline is
therefore *possible and logged*, never silent and never impossible.

---

## 4 · Unknown is not a value

A missing measurement must never be silently read as zero, false, or passing. The bug has three
recognisable shapes, and all three have been found in this pipeline:

- a blank marker percentage compared with `>=` evaluated to `False` and quietly meant "clean";
- an unknown doublet fraction read through `or 0` counted every cell in the cluster as surviving;
- a robust scale of zero, guarded with `or 1e-12`, turned a 0.07% difference into ~10¹² standard
  deviations.

Each looked exactly like a correct result. The pipeline now reports never-examined as its own
category everywhere it appears — in the tables, in the gates, and as its own colour in the
figures.

---

## What this costs

These rules make the pipeline slower to run end-to-end and more likely to stop and ask. That is
the intended trade. The failure mode they prevent — a defensible-looking result that quietly
deleted the effect you were measuring — does not announce itself, and is usually discovered, if at
all, long after it has been built upon.
