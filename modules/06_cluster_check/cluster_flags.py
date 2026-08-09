# Flags clusters against thresholds proposed from the cohort in hand, and returns a table.
# It removes nothing. A cluster-level removal deletes a whole population and is the most
# destructive class of removal in this pipeline; it requires an explicit assessment and a human
# decision, taken elsewhere.
"""Step 6 - the cluster check: flags, thresholds proposed from the cohort, nothing removed.

WHERE IT RUNS, AND WHY THERE IS NO CHOICE ABOUT IT

On the step-5 object: quality-filtered, with the doublet flags ATTACHED and NOT applied. After
step 7 every cluster is 0% doublet by construction and flag D becomes a tautology that passes
for the wrong reason. This is the only point at which D is computable.

THE FLAGS

  A  low RNA        cluster median UMI < f x that sample's OWN median. Sample-relative because
                    median depth differs across libraries; an absolute floor flags whole
                    libraries rather than clusters.
  B  high mito      cluster median %mt above a threshold. This is NOT the per-cell
                    mitochondrial ceiling and cannot be: in the calibration cohort no cluster
                    reached 27.5% against a 40% cut, so the per-cell value never fires at
                    cluster level.
  C  uninformative  share of the top-N markers in the locked mt+ribo set. Reported SPLIT into
                    C_mt and C_ribo - in the calibration cohort the ribosomal half was empty
                    (max 10%), which is expected for nuclei, where cytoplasmic ribosomal
                    transcripts are depleted, and folding the two together would have hidden
                    which one fired.
  D  doublet        cluster doublet frequency above a threshold. The default of 70% is the
                    cited method's number (doi:10.1038/s41467-025-64774-4), carried over rather
                    than measured here.

  FLAG  = (A and C) or (B and C) or D
  WATCH = C alone
  A alone, B alone and every continuous value are REPORTED, not flagged.

WHY A CONJUNCTION, AND WHAT IT COSTS

C alone is not a QC failure. In the calibration cohort every C-alone cluster was
mitochondrial-marker-driven at normal-to-deep coverage - the signature of a mitochondria-rich
cell type rather than of damage. A alone is low depth with INFORMATIVE markers, i.e. a genuine
low-RNA population: some cell types carry less RNA by biology. The conjunction removes both
false-positive classes - 22 clusters became 8.

The cost is that the result now depends on two cut-points at once. In the calibration cohort one
cluster sat at 15.67% mito (over B) and 45% markers (under C) and dropped out by five points on
a single axis. So the continuous values ship beside the boolean, always: a FLAG column is a
summary of two thresholds and a reader cannot see from it which one a cluster failed.

C, FLAG AND WATCH ARE MISSING WHERE MARKERS WERE NOT COMPUTED - NEVER FALSE

Markers are computed at the default resolution only; they are the expensive part. Everywhere
else C is UNKNOWN, and an unknown must never be written as False. `NaN >= 50` is False and
nothing objects, so a missing marker table quietly turns FLAG into the full rule at the default
resolution and into "D alone" at every other one; the flag counts then differ across a sweep
because of what was CALCULATED, not because of anything in the DATA. A gap must not read as a
cliff.

The same holds inside the conjunction, which is why the logic below is three-valued: `A and C`
with A unknown is UNKNOWN, not False, and `not A` with A unknown is not "A did not fire". A
blank must never be readable as a pass.

THE THRESHOLDS DO NOT TRANSFER

The calibration cohort's A < 0.5x, B > 15%, C >= 50% were chosen against its own distributions.
B = 15 is the p95 of a single cohort. C >= 50 works because that distribution was bimodal - the
median cluster had ZERO mt/ribo markers in its top 20 and the tail jumped to 43-55% - which is a
property of that data, not a guarantee. So this module PROPOSES thresholds from the cohort it is
given, and applies only what is passed to it.
"""

from __future__ import annotations

from dataclasses import dataclass, field


def _unknown(v) -> bool:
    """True when a value carries no information: None, or a NaN from a blank table cell.

    `nan is not None` is True and `nan >= x` is False, so a NaN slips through an `is not None`
    guard and then reads as a value that failed the test. Both spellings of unknown must be
    caught at the same place or the three-valued logic has a hole in it.
    """
    return v is None or (isinstance(v, float) and v != v)


D_DEFAULT = 70.0 # the cited method's number, not one measured here
TOPN = 20

class ClusterRefusal(RuntimeError):
    """Raised when the flags cannot be computed honestly."""

def _tri(v):
    """Normalise one operand to True, False or None, whatever type carries it.

    Identity is the wrong test and it inverted this module. `numpy.bool_(False) is False` is
    False, so `_and(np.False_, np.False_)` fell through the False branch and the unknown branch
    and returned True: a cluster that failed both criteria was reported as FLAGGED. Every profile
    this module reads comes from a numpy-backed table, so the wrong branch was the usual one.

    A string is read from its VOCABULARY, never by truthiness. `bool("False")` is True, so a
    profile round-tripped through CSV would have flagged every cluster; a word this function does
    not recognise is unknown, because guessing at it is how the previous defect happened.
    """
    if v is None:
        return None
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        s = v.strip().lower()
        if s in ("true", "t", "yes", "y", "1"):
            return True
        if s in ("false", "f", "no", "n", "0"):
            return False
        return None
    if isinstance(v, float) and v != v:                 # NaN, including numpy's
        return None
    if type(v).__name__ in ("NAType", "NaTType", "MaskedConstant"):
        return None
    try:
        if v != v:                                      # numpy/pandas NaN of any width
            return None
    except Exception:                                   # noqa: BLE001
        return None
    try:
        return bool(v)
    except Exception:                                   # noqa: BLE001
        return None


def _and(*vals):
    """Three-valued AND. False beats unknown; unknown beats True.

    False if any operand is False, None if any operand is unknown, True otherwise. `True and
    None` is None because the conjunction has not been evaluated, which is not the same as
    having been evaluated and passed.
    """
    vals = [_tri(v) for v in vals]
    if any(v is False for v in vals):
        return False
    if any(v is None for v in vals):
        return None
    return True

def _or(*vals):
    """Three-valued OR. True beats unknown; unknown beats False."""
    vals = [_tri(v) for v in vals]
    if any(v is True for v in vals):
        return True
    if any(v is None for v in vals):
        return None
    return False

def _not(v):
    """Three-valued NOT. The negation of an unknown is unknown, not True."""
    v = _tri(v)
    return None if v is None else (not v)

@dataclass
class Thresholds:
    a_umi_frac: float
    b_pct_mt: float
    c_uninformative: float
    d_doublet: float = D_DEFAULT
    source: str = "PROPOSED from the cohort - not approved"

    def __str__(self) -> str:
        return (f"A < {self.a_umi_frac}x sample median · B > {self.b_pct_mt:g}% mito · "
                f"C >= {self.c_uninformative:g}% markers · D > {self.d_doublet:g}% doublet"
                f"\n {self.source}")

def propose(profile, d_doublet=D_DEFAULT) -> Thresholds:
    """Derive candidate thresholds from THIS cohort's distributions, the way step 5 does.

    `profile` is a list of dicts with umi_frac_of_sample, median_pct_mt, pct_uninformative.
    """
    mt = sorted(r["median_pct_mt"] for r in profile if not _unknown(r.get("median_pct_mt")))
    un = sorted(r["pct_uninformative"] for r in profile
                if not _unknown(r.get("pct_uninformative")))
    if not mt or not un:
        raise ClusterRefusal("cannot propose thresholds without mitochondrial and marker "
                             "profiles for the default resolution")

    def q(v, p):
        return v[min(len(v) - 1, int(p * len(v)))]

    b = round(q(mt, 0.95), 1)
    # C is proposed at the midpoint of the gap between the bulk and the tail, IF there is a gap.
    # If the distribution is not bimodal the proposal is refused rather than fudged - the same
    # rule the count floors follow.
    nonzero = [x for x in un if x > 0]
    if not nonzero:
        raise ClusterRefusal(
            "no cluster has any uninformative markers - criterion C has nothing to threshold. "
            "Report it as absent rather than picking a cut")
    if q(un, 0.75) == 0 and max(un) > 0:
        # The MIDPOINT of the gap between the bulk and the tail. Check the algebra, not the
        # shape: an expression such as `(0 + m)/2 + m/2` reads as a midpoint and simplifies to
        # `m`, returning the tail's MINIMUM, so a cluster sitting exactly at the start of the
        # tail is flagged. Arithmetic that looks like an average and is not is worse than no
        # comment at all.
        bulk_top = max((x for x in un if x < min(nonzero)), default=0.0)
        c = float(round((bulk_top + min(nonzero)) / 2, 0))
        note = (f"marker share is bimodal - the bulk tops out at {bulk_top:.0f}% and the tail "
                f"starts at {min(nonzero):.0f}%; C proposed at the midpoint, {c:.0f}%")
    else:
        c = round(q(un, 0.95), 0)
        note = (f"marker share is NOT clearly bimodal (p75 = {q(un, 0.75):.0f}%); C proposed at "
                f"the p95, {c:.0f}% - a percentile, not a valley, so it is a weaker basis than "
                f"the count floors have")
    return Thresholds(0.5, b, c, d_doublet,
                      f"PROPOSED from this cohort - B at the p95 of cluster mito ({b:g}%); "
                      f"{note}. NOT approved; A is a convention, not a measurement")

@dataclass
class Flagged:
    rows: list = field(default_factory=list)

    def counts(self) -> dict:
        out = {}
        for k in ("A", "B", "C", "D", "FLAG", "WATCH"):
            out[k] = sum(1 for r in self.rows if r.get(k) is True)
        return out

    def unknown_counts(self) -> dict:
        """How many clusters each criterion was NOT evaluated on.

        `counts()` answers "how many fired". Without this a reader cannot tell a criterion that
        cleared every cluster from one that was never computed on any of them.
        """
        out = {}
        for k in ("A", "B", "C", "D", "FLAG", "WATCH"):
            out[k] = sum(1 for r in self.rows if _unknown(r.get(k)))
        return out

def apply_flags(profile, thr: Thresholds, markers_computed=True) -> Flagged:
    """Compute A/B/C/D and the conjunction. Each verdict is True, False or None - never coerced.

    None means NOT EVALUATED - the input it needs was absent - and it propagates: a conjunction
    with an unknown operand that is not already decided by a False is itself unknown, and so is
    its negation. C, FLAG and WATCH are None wherever markers were not computed.
    """
    out = []
    for r in profile:
        row = dict(r)
        row["A"] = (r["umi_frac_of_sample"] < thr.a_umi_frac
                    if not _unknown(r.get("umi_frac_of_sample")) else None)
        row["B"] = (r["median_pct_mt"] > thr.b_pct_mt
                    if not _unknown(r.get("median_pct_mt")) else None)
        row["D"] = (r["pct_doublet"] > thr.d_doublet
                    if not _unknown(r.get("pct_doublet")) else None)

        if not markers_computed or _unknown(r.get("pct_uninformative")):
            # MISSING, never False. A blank must not read as "evaluated and passed".
            #
            # FLAG is withheld here rather than reduced to its D term. Three-valued OR would
            # return True for a >70%-doublet cluster whatever C is, which is defensible in
            # isolation and wrong in a sweep: the column would then mean the full rule at the
            # resolution where markers exist and "D alone" everywhere else, under one name. D
            # itself is reported, so nothing is hidden by withholding the conjunction.
            row["C"] = row["C_mt"] = row["C_ribo"] = None
            row["FLAG"] = row["WATCH"] = None
        else:
            row["C"] = r["pct_uninformative"] >= thr.c_uninformative
            # A missing split is missing too: defaulting it to 0 would report "evaluated, below
            # threshold" for a quantity no file contains.
            row["C_mt"] = (r["pct_mt_markers"] >= thr.c_uninformative
                           if not _unknown(r.get("pct_mt_markers")) else None)
            row["C_ribo"] = (r["pct_ribo_markers"] >= thr.c_uninformative
                             if not _unknown(r.get("pct_ribo_markers")) else None)
            row["FLAG"] = _or(_and(row["A"], row["C"]),
                              _and(row["B"], row["C"]),
                              row["D"])
            row["WATCH"] = _and(row["C"], _not(row["A"]), _not(row["B"]))
        out.append(row)
    return Flagged(out)

def sweep_summary(by_resolution) -> list:
    """Only quantities computed at EVERY resolution. C, FLAG and WATCH are never included.

    Reporting FLAG across a sweep makes a gap in what was CALCULATED look like a cliff in the
    DATA. Each count is of clusters where the criterion FIRED, and each is shipped beside the
    number of clusters it could be evaluated on, so a small count cannot be misread as a
    quiet resolution when it is really a missing input.
    """
    rows = []
    for res, prof in sorted(by_resolution.items()):
        v = list(prof)
        row = {"resolution": res, "clusters": len(v)}
        for k in ("A", "B", "D"):
            row[k] = sum(1 for r in v if r.get(k) is True)
            row[f"{k}_evaluated"] = sum(1 for r in v if not _unknown(r.get(k)))
        for label, key in (("max_pct_mt", "median_pct_mt"),
                           ("max_pct_doublet", "pct_doublet")):
            # A max over an iterable containing None raises; a max over an absent key raises a
            # different way. Both are legitimate here - these fields are genuinely missing at
            # some resolutions - so the unknowns are dropped and counted, never coerced.
            known = [r[key] for r in v if r.get(key) is not None]
            row[label] = max(known) if known else None
            row[f"{label}_n"] = len(known)
        rows.append(row)
    return rows
