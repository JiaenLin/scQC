# Applies the agreed removals and writes the record of what left - the only step in this
# pipeline that removes anything. Every removal goes through `preflight()` and requires a
# recorded approval with the operator's verbatim words; there is no force flag. The masks below
# select what LEAVES, which is why they are gated.
"""Step 7 - apply: remove the flagged nuclei, and refuse to report a clean result over a
contradiction step 6 already found.

WHY STEP 6 FEEDS THIS STEP AT ALL

Built naively, step 6 is a dead end: it flags clusters, removes nothing, and step 7 applies
per-cell criteria that never consult it. Measured on the calibration cohort, the gap that leaves
is not small - 7,188 of the 7,734 nuclei in flagged clusters survive into the deliverable, 5.7%
of it.

The sharpest case is 41 nuclei in a cluster that is 74.4% doublet and were NOT called doublets
individually. That is precisely the residual the cited method's cluster rule exists to catch -
doublets that cluster TOGETHER and so escape a per-cell test - and the per-cell criteria cannot
see it by construction.

So step 6's output is not a filter of its own. It is EVIDENCE ABOUT STEPS 4 AND 5, and step 7
uses it three ways:

  1 PRE-FLIGHT. A cluster-level contradiction blocks a clean report. If a >70%-doublet cluster
     has survivors, the per-cell doublet call is under-removing where doublets cluster together,
     and saying "127,050 kept" without that is saying less than is known.
  2 CARRIED FORWARD. Every retained nucleus keeps its cluster's flags. It costs nothing, and
     without it the next stage recomputes a clustering to ask a question already answered.
  3 OFFERED, NOT TAKEN. A cluster-level removal is prepared in full - what it removes, listed;
     the design differential; what is lost if it is wrong - and then stops for a decision. The
     cited method removes these clusters; this pipeline requires a human to.

WHAT THE PRE-FLIGHT CANNOT DO

It cannot say the flagged nuclei are bad. Eight of the calibration cohort's nine flags are B&C -
high mitochondrial content with mitochondrial markers - and that is as consistent with a
mitochondria-rich cell type as with damage. Step 6 becomes useful by making the contradiction
visible and actionable, not by resolving it. Resolving it needs an identity, which is
annotation.

THE REMOVAL RECORD

A removal is recoverable only if what left can be named afterwards, so the gate does not accept
a count on its own. `build_removal_record()` takes the per-observation identifiers and one
boolean mask per criterion and returns one row per REMOVED observation, listing every criterion
that removed it; `write_removal_record()` persists it as CSV with the standard library alone.

  It guarantees   one row per removed observation, keyed by the identifier the caller
                  supplied; ALL criteria that fired for it, not just the first; and that the
                  record's own total and the count handed to the gate agree, or the removal is
                  refused rather than reported.
  It does not     recover the data. It stores identifiers, not counts - re-reading the input
                  with those identifiers is what recovers an observation. Nor can it know about
                  a criterion that was applied upstream and not passed in: the record then
                  reads as though the observation left for fewer reasons than it did, which is
                  the one failure mode to watch for.
  It refuses      duplicate identifiers, because a record whose key is ambiguous cannot answer
                  "what if that one had not been dropped", and an unknown mask entry, because a
                  removal decided on an unevaluated criterion is not a decision.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field


def _unknown(v) -> bool:
    """True when a value carries no information: None, or a NaN from a blank table cell.

    `nan is not None` is True and `nan >= x` is False, so a NaN slips through an `is not None`
    guard and then reads as a value that failed the test. Both spellings of unknown must be
    caught at the same place or the three-valued logic has a hole in it.
    """
    return v is None or (isinstance(v, float) and v != v)


D_CLUSTER = 70.0

class ApplyRefusal(RuntimeError):
    """Raised when a removal may not proceed."""

@dataclass
class PreflightFinding:
    check: str
    severity: str
    message: str
    detail: list = field(default_factory=list)

    def __str__(self) -> str:
        s = f"[{self.severity:6s}] {self.check}\n {self.message}"
        for d in self.detail:
            s += f"\n - {d}"
        return s

def preflight(cluster_profile, kept_total, d_cluster=D_CLUSTER) -> list:
    """Read step 6's table and report what the per-cell criteria did not catch.

    Removes nothing: it counts what SURVIVES.
    """
    out = []
    flagged = [c for c in cluster_profile if c.get("FLAG") is True]

    # 1. the doublet residual - the contradiction the per-cell call cannot see
    resid, rows = 0, []
    for c in cluster_profile:
        if not _unknown(c.get("pct_doublet")) and c["pct_doublet"] > d_cluster:
            surv = round(c["n"] * (1 - c["pct_doublet"] / 100))
            resid += surv
            # A median that was never computed is printed as unknown. Substituting 0 puts a
            # number in the report that no file on disk contains.
            mu = c.get("median_umi")
            mu_txt = f"{mu:,.0f} UMI" if mu is not None else "unknown"
            rows.append(f"{c.get('sample','?')} c{c.get('cluster','?')}: {surv:,} of {c['n']:,} "
                        f"survive a cluster that is {c['pct_doublet']:.1f}% doublet "
                        f"(median {mu_txt})")
    out.append(PreflightFinding(
        "doublet residual in cluster-level doublets",
        "REVIEW" if resid else "ok",
        (f"{resid:,} nuclei survive inside clusters above {d_cluster:g}% doublet. The per-cell "
         f"call cannot see doublets that cluster together; the cited method "
         f"(doi:10.1038/s41467-025-64774-4) removes such clusters entirely, and this pipeline "
         f"requires a decision rather than doing so silently"
         if resid else
         f"no cluster exceeds {d_cluster:g}% doublet"),
        rows))

    # 2. how much of the deliverable sits in a flagged cluster
    if flagged:
        # An UNKNOWN doublet fraction must not be folded in at 0% with `or 0`: that counts every
        # nucleus in the cluster as surviving, the same class of error as a blank C reading as
        # False. Clusters with no doublet fraction are reported separately instead.
        known = [c for c in flagged if not _unknown(c.get("pct_doublet"))]
        unknown = [c for c in flagged if _unknown(c.get("pct_doublet"))]
        surv = sum(round(c["n"] * (1 - c["pct_doublet"] / 100)) for c in known)
        if unknown:
            out.append(PreflightFinding(
                "flagged clusters with no doublet fraction", "REVIEW",
                f"{len(unknown)} flagged cluster(s) covering "
                f"{sum(c['n'] for c in unknown):,} nuclei have no doublet fraction. They are "
                f"NOT counted as 0% - an unknown is not a zero - and how many of them reach "
                f"the deliverable is unknown too"))
        out.append(PreflightFinding(
            "retained nuclei in flagged clusters", "REVIEW",
            f"{surv:,} of {kept_total:,} retained nuclei ({100*surv/max(kept_total,1):.1f}%) sit "
            f"in one of {len(flagged)} clusters flagged by step 6. This is NOT a claim that they "
            f"are bad - a cluster flagged for mitochondrial content cannot be told from a "
            f"mitochondria-rich cell type without an identity - but a deliverable reported "
            f"without it says less than is known"))

    # 3. criteria that never engaged at population level
    mt = [c["median_pct_mt"] for c in cluster_profile if c.get("median_pct_mt") is not None]
    if mt:
        out.append(PreflightFinding(
            "mitochondrial ceiling at cluster level", "ok",
            f"highest cluster median is {max(mt):.2f}% - the per-cell ceiling trims individual "
            f"nuclei and removes no population. Whether mitochondria-high POPULATIONS are cells "
            f"or damage is untouched by it"))
    return out

def annotate_kept(obs_rows, cluster_profile, cluster_key="cluster", sample_key="sample"):
    """Carry each nucleus's cluster flags into the deliverable.

    Adds columns; nothing is dropped. Without this the next stage recomputes a clustering to ask
    a question step 6 answered.
    """
    idx = {(c.get(sample_key), c.get(cluster_key)): c for c in cluster_profile}
    for r in obs_rows:
        c = idx.get((r.get(sample_key), r.get(cluster_key)))
        r["cluster_FLAG"] = None if c is None else c.get("FLAG")
        r["cluster_WATCH"] = None if c is None else c.get("WATCH")
        r["cluster_pct_doublet"] = None if c is None else c.get("pct_doublet")
        r["cluster_median_pct_mt"] = None if c is None else c.get("median_pct_mt")
    return obs_rows

def propose_cluster_removal(cluster_profile, design=None, d_cluster=D_CLUSTER) -> dict:
    """Prepare a cluster-level removal in full, and stop. Never applies it.

    Everything a decision needs: what goes, how much, whether it is even across the design.
    """
    victims = [c for c in cluster_profile
               if not _unknown(c.get("pct_doublet")) and c["pct_doublet"] > d_cluster]
    n = sum(c["n"] for c in victims)
    prop = {
        "rule": f"remove clusters with doublet frequency > {d_cluster:g}%",
        "source": "doi:10.1038/s41467-025-64774-4, the half of the cited method not applied "
                  "per cell",
        "clusters": [f"{c.get('sample')} c{c.get('cluster')} "
                     f"(n={c['n']:,}, {c['pct_doublet']:.1f}% doublet)" for c in victims],
        "n_removed": n,
        "status": "PROPOSED - not applied. A cluster-level removal deletes a whole population "
                  "and is the most destructive class of removal in this pipeline",
    }
    if design and victims:
        by, no_level = {}, 0
        for c in victims:
            lvl = design.get(c.get("sample"))
            # A level of "", 0 or False is a level. Testing truthiness drops it, the tally then
            # looks even or empty, and the one warning this function exists to raise never
            # fires. Only an ABSENT sample - `None` - is unknown, and it is counted as unknown.
            if not _unknown(lvl):
                by[lvl] = by.get(lvl, 0) + c["n"]
            else:
                no_level += c["n"]
        prop["by_design_level"] = by
        if no_level:
            prop["n_no_design_level"] = no_level
            prop["design_note"] = (
                f"{no_level:,} of the {n:,} nuclei belong to samples absent from the design "
                f"map. They are not assigned to a level and the tally below does not cover them")
        if len(by) == 1:
            prop["warning"] = (
                f"every removed nucleus with a known level falls on a single design level "
                f"({list(by)[0]!r}). A one-sided removal converts a technical property into an "
                f"apparent biological difference, which is what the design-differential check "
                f"exists to catch, and this one cannot be checked by ratio because the other "
                f"level is zero")
    return prop

@dataclass
class RemovalRecord:
    """One row per removed observation, with every criterion that removed it.

    `rows` is a list of (identifier, [criterion, ...]) in input order. `criteria` is the full
    set of criteria the record was built over, including any that removed nothing - a criterion
    with zero rows is evidence it was evaluated, which a record listing only what fired cannot
    give.
    """
    criteria: list = field(default_factory=list)
    rows: list = field(default_factory=list)
    n_in: int = 0

    @property
    def n_removed(self) -> int:
        return len(self.rows)

    @property
    def n_kept(self) -> int:
        return self.n_in - self.n_removed

    def by_criterion(self) -> dict:
        """Removed observations per criterion. The sum exceeds n_removed when criteria overlap."""
        out = {k: 0 for k in self.criteria}
        for _, why in self.rows:
            for k in why:
                out[k] += 1
        return out

    def identifiers(self) -> list:
        return [i for i, _ in self.rows]

def build_removal_record(identifiers, criterion_masks) -> RemovalRecord:
    """Name what leaves. Takes the per-observation identifiers and one boolean mask per
    criterion, and returns the removed observations with the criteria that removed them.

    `identifiers` is a sequence of length N whose entries are unique. `criterion_masks` maps a
    criterion name to a sequence of N booleans, True where that criterion removes the
    observation. An observation is removed if ANY mask is True, and it carries ALL of them.
    """
    ids = list(identifiers)
    n = len(ids)
    if not criterion_masks:
        raise ApplyRefusal(
            "a removal record needs at least one named criterion. A record that cannot say WHY "
            "an observation left is a list of casualties, not a record")

    seen, dupes = set(), []
    for i in ids:
        if i in seen:
            dupes.append(i)
        seen.add(i)
    if dupes:
        raise ApplyRefusal(
            f"{len(dupes)} duplicate identifier(s) among {n} observations, e.g. "
            f"{dupes[:5]!r}. A record keyed on an ambiguous identifier cannot answer 'what if "
            f"that one had not been dropped'. Make the key unique - barcodes repeat across "
            f"samples, so combine the sample with the barcode - rather than merging them here")

    names = list(criterion_masks)
    masks = {}
    for name in names:
        m = list(criterion_masks[name])
        if len(m) != n:
            raise ApplyRefusal(
                f"criterion {name!r} has {len(m)} mask entries for {n} observations. A mask "
                f"that does not line up with the identifiers records the wrong observations, "
                f"and nothing downstream can detect that it did")
        unknown = sum(1 for v in m if v is None)
        if unknown:
            raise ApplyRefusal(
                f"criterion {name!r} is unevaluated for {unknown} of {n} observations. An "
                f"unknown is not a False: removing on a criterion that was never computed, or "
                f"recording it as not having fired, both state more than was measured")
        masks[name] = [bool(v) for v in m]

    rows = []
    for pos, ident in enumerate(ids):
        why = [name for name in names if masks[name][pos]]
        if why:
            rows.append((ident, why))
    return RemovalRecord(criteria=names, rows=rows, n_in=n)

def write_removal_record(record: RemovalRecord, path) -> str:
    """Persist the record as CSV: one row per removed observation, one column per criterion.

    Standard library only, so reading it back needs nothing this pipeline installed. Columns are
    `identifier`, `n_criteria`, `criteria` (the names, `|`-separated) and one 0/1 column per
    criterion in the order the criteria were supplied.
    """
    cols = ["identifier", "n_criteria", "criteria"] + list(record.criteria)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(cols)
        for ident, why in record.rows:
            fired = set(why)
            w.writerow([ident, len(why), "|".join(why)]
                       + [1 if k in fired else 0 for k in record.criteria])
    return str(path)

class Kept(int):
    """The number of observations that remain - an int that carries what left.

    It compares, formats and arithmetics exactly as the plain count did, so an existing caller
    reading it as a number is unaffected; `.record` and `.record_path` are there for one that
    needs to say which observations went and why.
    """

    def __new__(cls, value, record=None, record_path=None):
        obj = super().__new__(cls, value)
        obj.record = record
        obj.record_path = record_path
        return obj

def _as_count(value, name, action):
    """An observation count is a non-negative whole number or the removal does not proceed."""
    try:
        n = int(value)
    except (TypeError, ValueError) as exc:
        raise ApplyRefusal(
            f"{name} is not a count - refused: {action!r}\n got {value!r}, which is not a "
            f"number of observations") from exc
    if n != value or isinstance(value, bool):
        raise ApplyRefusal(
            f"{name} is not a whole number - refused: {action!r}\n got {value!r}; an "
            f"observation count must be whole, and silently truncating one hides an upstream "
            f"mistake")
    return n

def apply_removal(n_in, removed_mask_sum, action, user_verbatim, approvals,
                  record: RemovalRecord = None, record_path=None) -> Kept:
    """The gate. Removes nothing without the operator's own words recorded against THIS action.

    Returns the number kept, carrying `record` and the path it was written to when one is
    supplied. Pass the record from `build_removal_record()` to make the removal recoverable:
    the gate then checks it against the count it was given and refuses if the two disagree,
    because a number derived one way is a number nothing has checked.
    """
    # Arithmetic first: an approval for an impossible removal is not an approval of anything.
    n_in = _as_count(n_in, "n_in", action)
    removed = _as_count(removed_mask_sum, "removed_mask_sum", action)
    if removed < 0 or removed > n_in:
        raise ApplyRefusal(
            f"impossible arithmetic - refused: {action!r}\n {removed:,} observations cannot be "
            f"removed from {n_in:,}. A mask sum outside 0..n_in means the mask and the object "
            f"have come apart - subtracting anyway yields a plausible-looking count of a "
            f"population that does not exist.")
    if not str(user_verbatim).strip():
        raise ApplyRefusal(
            f"no verbatim words - refused: {action!r}\n a removal requires the operator's own "
            f"words. Do not summarise them, do not paraphrase them, and do not reuse words "
            f"spoken about a different action - an approval recorded against a description "
            f"that has drifted from the action is not an approval.")
    if action not in approvals:
        raise ApplyRefusal(
            f"no approval recorded - refused: {action!r}\n no approval exists for this exact "
            f"action text. The action string is what an approval is matched against; if the "
            f"text changed, the approval must be given again. There is no force flag.")
    # Compared on stripped text. An approval stored with a trailing newline would otherwise
    # refuse a correct CONFIRM, and a gate that fires on correct behaviour gets switched off.
    # Stripping whitespace does not weaken it: the words must still be the operator's own, for
    # THIS action.
    if approvals[action].strip() != str(user_verbatim).strip():
        raise ApplyRefusal(
            f"approval text does not match - refused: {action!r}\n the recorded approval is "
            f"{approvals[action]!r}, not {user_verbatim!r}. Reusing words spoken about "
            f"something else is how an instruction becomes a fabricated consent.")

    if record is not None:
        if record.n_in != n_in:
            raise ApplyRefusal(
                f"record describes another object - refused: {action!r}\n the removal record "
                f"covers {record.n_in:,} observations and the gate was given {n_in:,}.")
        if record.n_removed != removed:
            raise ApplyRefusal(
                f"record disagrees with the mask - refused: {action!r}\n the removal record "
                f"names {record.n_removed:,} removed observations and the mask sums to "
                f"{removed:,}. Two independent routes to the same number disagree; the removal "
                f"stops rather than reporting whichever was asked for first.")
    elif record_path is not None:
        raise ApplyRefusal(
            f"nothing to record - refused: {action!r}\n a record path was given with no record "
            f"to write. The file would be created empty and read as 'nothing was removed'.")

    written = write_removal_record(record, record_path) if record_path is not None else None
    return Kept(n_in - removed, record=record, record_path=written)
