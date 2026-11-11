# Decides whether ambient correction runs and refuses an unsafe skip.
# It removes no observation itself; CellBender adjusts counts and deletes no barcode.
"""Step 1 — ambient RNA removal. MANDATORY for single-nuclei, NOT RUN on single cells.

WHY THE ASYMMETRY IS STRUCTURAL, NOT A PREFERENCE

A nucleus holds roughly an order of magnitude less RNA than the cell it came from, and the
ambient pool in a nuclear prep IS the cytoplasm of the cells that were lysed to make it. The same
absolute background is therefore a far larger fraction of a nuclear transcriptome than of a whole
cell's, and where the dominant cell type is large, RNA-rich and mitochondria-dense, lysing it
floods the pool.

Measured across the calibration cohort's ten solid-tissue libraries, CellBender 0.3.2:

    fraction of all counts removed as ambient    15.5% - 23.7%
    genes removed entirely                       12 - 240 per library

Between a sixth and a quarter of the data. A pipeline that leaves that in place for snRNA-seq is
not applying a lighter touch; it is analysing a mixture and calling it a cell.

For scRNA-seq the cell retains its cytoplasm and the ambient fraction is ordinarily much lower.
This read "default-on-but-skippable for `scrna`" until 2026-10-02, and the run route ignored even
that: it scheduled CellBender for every library whatever the assay, so the first single-cell
cohort spent its GPU hours denoising whole cells and failed there (single-cell-harness ADR-0027).
The PI's ruling that day: "do not run cellbender for single cell" - "single cell and single
nuclei is fundamentally different!". So the rule is a POLICY per assay, `CELLBENDER_RUNS` below:
CellBender is mandatory for `snrna` and is never run on `scrna`. The words name CellBender, and
the policy is read no wider than them (harness autodev W020, B1-075): no other ambient method is
built for cells here yet, and one may be added, per assay, without touching this rule. A
single-cell library arrives already corrected (SUPPLIED, with its provenance) or takes its raw
counts as the object, and the report says which. Its cell call is then the aligner's, because
there is no denoiser to take one from.

A NUCLEUS LIBRARY'S CORRECTION MAY BE SKIPPED BY DECLARATION (harness autodev W020, B1-105).
"There is no flag for this" made the skip impossible, against the project's rule that a bypass is
possible and logged, never impossible (docs/PRINCIPLES.md section 3); a site without a GPU queue
could not run nuclei at all. A skip needs a reason AND a named approver (samplesheet columns
`ambient_skip_reason`, `ambient_skip_approved_by`); it is recorded DECLARED and stated in the
report, and the library then takes the not-denoised route - raw counts, the aligner's call - with
every other nucleus rule unchanged. Half a declaration skips nothing.

CELLBENDER IS A DENOISER, NOT THIS PIPELINE'S CELL CALLER

Its cell/empty split is a by-product. Cell selection belongs to the count thresholds, doublet
detection and annotation downstream. Judge a run on its denoising, not on its cell count - and
prefer permissive calling, because a nucleus wrongly called empty is gone before any later step
can look at it. The pipeline therefore carries the cell call forward as a COLUMN and uses it only
as the empty-droplet boundary the doublet tools require.

THE DECLARATION IS CHECKED, NOT TRUSTED

`assay` is DECLARED because the pipeline cannot safely infer it - but it is also MEASURABLE, and
the two must agree. snRNA-seq carries a high intronic fraction and scRNA-seq a low one, so if
step 0 emits both intron-inclusive and exon-only counts the fraction is a direct read-out. A
declaration of `scrna` on data with a nuclear intronic fraction means either the wrong assay was
declared - and mandatory correction is about to be skipped - or the wrong reference was used.
Either way it must stop the run rather than be resolved by whichever value happened to be typed.
"""

from __future__ import annotations

from dataclasses import dataclass

ASSAYS = ("snrna", "scrna")

# WHETHER CELLBENDER RUNS, PER ASSAY. A policy, not a default a flag can flip: a nucleus is
# denoised because the ambient pool is a sixth to a quarter of its counts, and CellBender is never
# run on a whole cell, by the PI's ruling of 2026-10-02 (see above). Every assay in ASSAYS has an
# entry, and a new assay must state its own rather than inherit either. CellBender-specific: an
# ambient method built for cells would have its own table beside this one.
CELLBENDER_RUNS = {"snrna": True, "scrna": False}
#: Why the single-cell policy is what it is, in the words that set it.
CELLBENDER_POLICY_SOURCE = ("the PI, 2026-10-02: \"do not run cellbender for single cell\" - "
                            "\"single cell and single nuclei is fundamentally different!\"")
# CellBender package defaults. A departure from either - a halved learning rate, for instance -
# is a stated choice, not a setting, and must be recorded as one. `lr_policy` decides when the
# learning rate has to move and requires the halved run to be RE-MEASURED before it is adopted.
DEFAULTS = {"fpr": 0.0, "learning_rate": 1e-4}
# Intronic fraction is bimodal between the assays. These bounds are deliberately WIDE - they
# exist to catch a declaration that is plainly wrong, not to adjudicate a borderline case.
INTRONIC_EXPECT = {"snrna": (0.35, 1.00), "scrna": (0.00, 0.35)}

@dataclass
class AmbientPlan:
    sample: str
    assay: str
    run: bool
    mandatory: bool
    reason: str = ""
    params: dict = None
    # SUPPLIED is a third state and NOT a kind of skip. See plan_ambient.
    supplied: dict = None
    # A skip DECLARED for a nucleus library: {reason, approved_by}. None otherwise (B1-105).
    declared_skip: dict = None

    @property
    def state(self) -> str:
        return "RUN" if self.run else ("SUPPLIED" if self.supplied else "SKIP")

    def __str__(self) -> str:
        tag = "mandatory" if self.mandatory else "not denoised here"
        s = f"[{self.state:8s}] {self.sample:14s} assay={self.assay:6s} ({tag})"
        if self.reason:
            s += f"\n {self.reason}"
        return s

class AmbientRefusal(RuntimeError):
    """Raised when a skip is requested where correction is mandatory."""

# Provenance a supplied correction must carry. Without these the object is an anonymous matrix
# asserted to be denoised, which is not a claim anything downstream can check or a reader trust.
SUPPLIED_REQUIRED = ("tool", "version", "params", "produced_by")


def plan_ambient(sample, assay, skip=False, skip_reason="", intronic_fraction=None,
                 params=None, supplied=None, approved_by="") -> AmbientPlan:
    """Decide whether CellBender runs for this sample, and refuse an unsafe skip.

    THREE STATES, AND THE THIRD IS NOT A KIND OF SKIP

      RUN       this pipeline corrects the matrix.
      SKIP      no correction happens here. For scRNA it is the policy (CELLBENDER_RUNS), not a
                request; for snRNA it is refused unless DECLARED with a reason and a named
                approver (`skip_reason`, `approved_by`), and is then recorded as declared.
      SUPPLIED  the correction ALREADY HAPPENED, elsewhere, and its output is the input here.

    Conflating SUPPLIED with SKIP is the failure this distinction exists to prevent, and it can
    go wrong in both directions. Treating a supplied object as a skip refuses a perfectly valid
    re-entry into the pipeline at step 2 - the common case of "I already denoised these, do the
    rest". Treating it as a run lets the report claim this pipeline performed a correction it
    never performed, and attaches this run's parameters to someone else's output.

    So SUPPLIED demands provenance: which tool, which version, which parameters, and who ran it.
    A matrix merely asserted to be denoised is an anonymous matrix. Downstream cannot verify the
    claim - the ambient AUDIT needs the raw counts to measure what was removed, and by
    construction a supplied run does not have them - so what cannot be measured must at least be
    attributed.
    """
    assay = str(assay).strip().lower()
    if assay not in ASSAYS:
        raise AmbientRefusal(
            f"{sample}: assay '{assay}' is not one of {ASSAYS}. It is DECLARED and has no "
            f"default - the pipeline cannot infer it, and guessing wrong either skips a "
            f"mandatory correction or applies an unnecessary one.")

    mandatory = assay == "snrna"
    p = {**DEFAULTS, **(params or {})}

    # Check 2: the declaration against the data. Wide bounds, so this fires only on a plain
    # contradiction rather than on a borderline library.
    if intronic_fraction is not None:
        lo, hi = INTRONIC_EXPECT[assay]
        if not (lo <= intronic_fraction <= hi):
            other = "snrna" if assay == "scrna" else "scrna"
            raise AmbientRefusal(
                f"{sample}: declared assay '{assay}' but the intronic fraction is "
                f"{intronic_fraction:.2f}, outside the {lo:.2f}-{hi:.2f} expected for it and "
                f"consistent with '{other}'. Either the assay is mis-declared - in which case "
                f"{'a mandatory correction is about to be skipped' if assay == 'scrna' else 'an unnecessary one is about to run'} "
                f"- or the reference is wrong. Resolve it; do not pick one.")

    if supplied is not None:
        if skip:
            raise AmbientRefusal(
                f"{sample}: `skip` and `supplied` are contradictory. A supplied correction is "
                f"one that already happened; a skip is one that never does. Pass one.")
        if not isinstance(supplied, dict):
            raise AmbientRefusal(
                f"{sample}: `supplied` must be a provenance record, not {type(supplied).__name__}. "
                f"Required keys: {', '.join(SUPPLIED_REQUIRED)}.")
        missing = [k for k in SUPPLIED_REQUIRED
                   if not str(supplied.get(k) or "").strip()]
        if missing:
            raise AmbientRefusal(
                f"{sample}: a SUPPLIED ambient correction is missing {', '.join(missing)}. "
                f"This pipeline did not perform the correction, so it cannot describe it, and an "
                f"object merely asserted to be denoised is an anonymous matrix. The ambient audit "
                f"cannot check the claim either - measuring what was removed needs the raw counts, "
                f"which a supplied run does not have. What cannot be measured must be attributed.")
        return AmbientPlan(
            sample, assay, False, mandatory,
            f"SUPPLIED: corrected by {supplied['tool']} {supplied['version']} "
            f"({supplied['params']}), run by {supplied['produced_by']}. NOT performed by this "
            f"run - the fraction removed is NOT MEASURABLE here and must not be reported as if "
            f"it were.",
            p, supplied=dict(supplied))

    if not CELLBENDER_RUNS[assay]:
        # NOT a skip somebody asked for, so it needs no reason from the caller: it is this
        # pipeline's rule for the assay, and the rule is the reason. A reason given anyway is
        # kept beside it.
        return AmbientPlan(
            sample, assay, False, mandatory,
            f"NOT RUN: CellBender is not run on {assay} ({CELLBENDER_POLICY_SOURCE}). No other "
            f"ambient method is built for cells here yet; one may be added per assay. The "
            f"ambient fraction is NOT MEASURED and must not be reported as low; the cell call is "
            f"the aligner's."
            + (f" Also recorded: {skip_reason}" if str(skip_reason).strip() else ""),
            p)

    # Only an assay CellBender runs on reaches here, and on every such assay it is mandatory. A
    # skip is DECLARED or it is refused: a reason without a named approver, or an approver
    # without a reason, is half a declaration and skips nothing.
    if skip:
        reason, who = str(skip_reason or "").strip(), str(approved_by or "").strip()
        if reason and who:
            return AmbientPlan(
                sample, assay, False, mandatory,
                f"SKIPPED BY DECLARATION: ambient correction is mandatory for {assay} and was "
                f"skipped for this library - reason: {reason}; approved by {who}. DECLARED, not "
                f"the policy: the ambient fraction is NOT MEASURED and must not be reported as "
                f"low, and the cell call is the aligner's.",
                p, declared_skip={"reason": reason, "approved_by": who})
        missing = ", ".join(x for x, v in (("a reason (ambient_skip_reason)", reason),
                                           ("a named approver (ambient_skip_approved_by)", who))
                            if not v)
        raise AmbientRefusal(
            f"{sample}: ambient correction is mandatory for snRNA-seq, and a skip must be "
            f"DECLARED with a reason and a named approver - this one lacks {missing}. A nucleus "
            f"holds roughly an order of magnitude less RNA than its cell and the ambient pool IS "
            f"the lysed cytoplasm, so the background is a far larger fraction of the signal - "
            f"measured as a large share of all counts in solid-tissue libraries. Skipping "
            f"it analyses a mixture, and a person has to say so, in writing, before it is done.")

    return AmbientPlan(sample, assay, True, mandatory, "ambient correction will run", p)
