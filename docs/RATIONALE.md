# Rationale

Why each part of the pipeline exists, and what goes wrong without it. Moved here from the
README so that document can describe the tool; the reasoning is unchanged.

See also [PRINCIPLES.md](PRINCIPLES.md) for the four enforced rules and
[FILTERS.md](FILTERS.md) for how each threshold is calculated.

---

## What motivated it

Three failures motivated this pipeline. Each is quiet, each survives review, and none is detected
by a run that completes successfully.

- **A delivered matrix that was already filtered.** Pre-filtered input cannot be un-filtered, and
  nothing downstream can tell. scQC verifies that its input is raw and refuses otherwise.
- **A doublet rate set by its prior rather than by the data.** The detector's default expected rate
  dominated the call. scQC sweeps that parameter and reports whether the rate is a measurement.
- **A threshold that transferred badly.** Across ten libraries of a single cohort — one tissue, one
  platform, one operator — the measured density valley ranged **274–473 UMI**. A hard-coded
  `--min_umi 500` sits above the valley in all ten. *The procedure transfers; the number does not.*

---

## Why each step is there

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
