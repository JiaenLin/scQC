# Decides whether ambient correction runs and refuses an unsafe skip.
# It removes no observation itself; CellBender adjusts counts and deletes no barcode.
"""Step 1 — ambient RNA removal. MANDATORY for single-nuclei, optional for single-cell.

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

For scRNA-seq the cell retains its cytoplasm, the ambient fraction is ordinarily much lower, and
whether to correct is a judgement about a particular experiment. So: mandatory for `snrna`,
default-on-but-skippable for `scrna`, and a skip must be RECORDED rather than silent.

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

    def __str__(self) -> str:
        head = "RUN" if self.run else "SKIP"
        tag = "mandatory" if self.mandatory else "optional"
        s = f"[{head:4s}] {self.sample:14s} assay={self.assay:6s} ({tag})"
        if self.reason:
            s += f"\n {self.reason}"
        return s

class AmbientRefusal(RuntimeError):
    """Raised when a skip is requested where correction is mandatory."""

def plan_ambient(sample, assay, skip=False, skip_reason="", intronic_fraction=None,
                 params=None) -> AmbientPlan:
    """Decide whether CellBender runs for this sample, and refuse an unsafe skip."""
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

    if not skip:
        return AmbientPlan(sample, assay, True, mandatory,
                           "ambient correction will run", p)

    if mandatory:
        raise AmbientRefusal(
            f"{sample}: ambient correction cannot be skipped for snRNA-seq. There is no flag "
            f"for this. A nucleus holds roughly an order of magnitude less RNA than its cell "
            f"and the ambient pool IS the lysed cytoplasm, so the background is a far larger "
            f"fraction of the signal - measured at 15.5-23.7% of all counts across ten "
            f"solid-tissue libraries. Skipping it does not analyse nuclei; it analyses a "
            f"mixture.")

    if not str(skip_reason).strip():
        raise AmbientRefusal(
            f"{sample}: skipping ambient correction for scRNA-seq is permitted but must be "
            f"RECORDED. Give a reason - it is written to the run manifest, the same way every "
            f"other deliberate bypass in this pipeline is. A skip that is possible and logged "
            f"stays reviewable; a silent one is indistinguishable from an oversight.")

    return AmbientPlan(sample, assay, False, mandatory,
                       f"SKIPPED, recorded reason: {skip_reason}", p)
