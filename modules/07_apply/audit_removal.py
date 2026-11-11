# This module AUDITS a completed removal. It removes nothing and it writes no object.
# rule-one: no-removal - every mask here selects rows to COUNT.
"""Step 7 audit: was the removal exactly what the pipeline says it was?

WHY THIS EXISTS

Steps 2 and 3 restrict a population without removing an observation, and step 4 flags without
removing. That is the design, and it is also the design's weak point: a restriction that drops
nothing leaves no gap for anyone to notice. Three failures are invisible without an audit -

  * a nucleus never EXAMINED for doublets counted as a nucleus examined and found clean;
  * a rate computed over all droplets when its denominator should be called cells;
  * a removal that happened somewhere other than step 7, so the ledger under-reports it.

None of those makes an object smaller than expected. They make the NUMBERS DESCRIBING it wrong
while every count still reconciles, which is why they need a check that reconciles the counts
deliberately rather than a reader who notices something missing.

WHAT IT CHECKS, AND WHAT IT CANNOT

It checks that the deliverable is EXACTLY the recorded criteria applied jointly, that nothing was
removed without a recorded reason, that no unexamined nucleus survived carrying a fabricated
"clean" verdict, and that the only object with a different cell set is the cluster-check
intermediate. It cannot check whether the criteria were the RIGHT ones - that is the operator's
judgement and steps 5 and 6 exist to inform it.

THE ONE SANCTIONED EXCEPTION

The cluster check needs an object where doublets are still present, because once the doublet
criterion is applied every cluster is 0% doublet by construction and the cluster-level doublet
test becomes a tautology that passes for the wrong reason. That object is a STRICT SUPERSET of
the deliverable and differs from it by doublets ALONE. Both halves are checked: a superset that
differs by anything else is not the sanctioned exception, it is a second, unaudited filter.

A FAILED CHECK STOPS THE RUN (harness autodev W043, GAP-119)

Every check that finds the removal is not what the pipeline recorded is a REFUSE, the word every
other gate writes, and step 7 hands `verdict()` to `Pipeline.gate`, which stops the run there -
after the ledger, before any matrix is written. This module wrote FAIL, and the gate stops on
REFUSE only, so a failed audit was recorded, read by the report as an unrecognised severity, and
the deliverable was written anyway.

ONE OF THESE REFUSALS CAN BE ACCEPTED, AND THE REST NEVER (W043, fix round 1)

Only the floor ordering (ACCEPTABLE) is a gate in docs/PRINCIPLES.md's sense: it compares two
parameters, so it can fire on a correct run, and an `accept` entry naming this step, the check,
who accepted it and their words lifts it. Every other REFUSE (INTEGRITY) compares the per-cell
table with itself. A table 07_measure writes cannot fail one, and accepting one would deliver an
object that contradicts its own ledger, so `integrity()` raises them before the gate and no
accept entry lifts them. They were lifted by name: a check accepted while it held lifted it
later, when an edited table made the deliverable hold a cell its ledger lists as removed.
"""

from __future__ import annotations

from dataclasses import dataclass, field


class AuditFailure(RuntimeError):
    """Raised when the removal is not what the pipeline recorded."""


STEP = "07_apply"
# THE ONE REFUSAL HERE AN ACCEPT ENTRY MAY LIFT. It compares two parameters, so it fires on a run
# whose floor was adjudicated under the light floor with every called cell above both, and what it
# protects is checked row by row: by "step 3: no unexamined nucleus reaches the deliverable" here,
# and by 07_measure's own stop on a kept cell nothing scored.
ACCEPTABLE = ("step 3: the light floor is strictly below the quality floor",)
# THE INTEGRITY CHECKS: each compares the per-cell table with itself, and none can fire on a table
# 07_measure wrote (keep is the complement of removed, removed the OR of every criterion column).
# `integrity()` fails closed: a REFUSE check added later and not named in ACCEPTABLE is one of these.
INTEGRITY = (
    "keep and removed are complementary",
    "removal decomposes into the recorded criteria",
    "step 2: no uncalled droplet reaches the deliverable",
    "step 3: no unexamined nucleus reaches the deliverable",
    "step 4: no doublet call outside the examined set",
    "step 6 object is a superset of the deliverable",
    "step 6 object differs from the deliverable by doublets ALONE",
)


def _norm(x) -> str:
    # The comparison Pipeline._apply_accepts makes, so an entry is refused exactly where the gate
    # would have matched it.
    return " ".join(str(x or "").split()).lower()


@dataclass
class Finding:
    check: str
    severity: str          # "ok" | "REVIEW" | "REFUSE"
    message: str
    detail: list = field(default_factory=list)

    def __str__(self) -> str:
        s = f"[{self.severity:6s}] {self.check}\n         {self.message}"
        for d in self.detail[:12]:
            s += f"\n           - {d}"
        if len(self.detail) > 12:
            s += f"\n           - ... and {len(self.detail) - 12} more"
        return s


def _col(rows, name):
    if name not in rows:
        raise AuditFailure(
            f"the per-cell table has no column {name!r}. The audit cannot pass by skipping a "
            f"check it could not run - an absent column is a failed audit, not a silent one. "
            f"Columns present: {', '.join(sorted(rows))}")
    return rows[name]


def audit(rows, *, cell_col="cellbender_cell", keep_col="keep", removed_col="removed",
          criteria=(), scored_col=None, light_floor=None, quality_floor=None,
          umi_col="total_counts", predoublet_keep=None, doublet_criterion=None) -> list:
    """Audit a completed removal from the per-cell table.

    `rows` maps column name -> a sequence, one entry per DROPLET (not per kept cell) - the table
    must cover everything the pipeline saw, or the checks below can only confirm the survivors
    agree with themselves.

    `criteria` names the columns holding each applied criterion, without the `fail_` prefix
    convention imposed: pass the actual column names. They are what "removed" must decompose into.
    """
    out: list = []
    cell = [bool(x) for x in _col(rows, cell_col)]
    keep = [bool(x) for x in _col(rows, keep_col)]
    removed = [bool(x) for x in _col(rows, removed_col)]
    n = len(cell)
    crit = {c: [bool(x) for x in _col(rows, c)] for c in criteria}

    # --- A. keep and removed are complements, not two independently maintained columns ----------
    both = [i for i in range(n) if keep[i] and removed[i]]
    neither = [i for i in range(n) if not keep[i] and not removed[i]]
    out.append(Finding(
        "keep and removed are complementary", "REFUSE" if (both or neither) else "ok",
        (f"{len(both)} row(s) are BOTH kept and removed and {len(neither)} are neither - the two "
         f"columns disagree, so at least one of them is not what the object contains"
         if (both or neither) else
         f"every one of {n:,} rows is exactly one of kept or removed")))

    # --- B. the deliverable is EXACTLY the criteria applied jointly ------------------------------
    # Recomputed from the criterion columns rather than trusted. This is the check that catches a
    # removal performed somewhere other than step 7: such a nucleus is `removed` with no criterion
    # explaining it, and nothing else in the pipeline would notice.
    if crit:
        rebuilt = [(not cell[i]) or any(crit[c][i] for c in crit) for i in range(n)]
        mismatch = [i for i in range(n) if rebuilt[i] != removed[i]]
        unexplained = [i for i in range(n)
                       if removed[i] and cell[i] and not any(crit[c][i] for c in crit)]
        out.append(Finding(
            "removal decomposes into the recorded criteria", "REFUSE" if mismatch else "ok",
            (f"{len(mismatch):,} row(s) differ between the recorded `removed` and the join of "
             f"{len(crit)} criteria. Of those, {len(unexplained):,} are called cells removed with "
             f"NO criterion set - a removal that happened outside step 7 looks exactly like this"
             if mismatch else
             f"`removed` reproduces exactly from {len(crit)} criteria over {n:,} rows"),
            [f"row {i}" for i in mismatch[:12]]))

    # --- C. step 2 did not contaminate: no non-cell survives -------------------------------------
    leak = [i for i in range(n) if keep[i] and not cell[i]]
    out.append(Finding(
        "step 2: no uncalled droplet reaches the deliverable", "REFUSE" if leak else "ok",
        (f"{len(leak):,} droplet(s) not called a cell are in the deliverable"
         if leak else
         f"all {sum(keep):,} kept observations were called cells; the {sum(1 for c in cell if not c):,} "
         f"uncalled droplets are accounted for as a recorded criterion, not dropped silently")))

    # --- D. step 3 did not contaminate ------------------------------------------------------------
    if scored_col is not None:
        scored = [bool(x) for x in _col(rows, scored_col)]
        kept_unscored = [i for i in range(n) if keep[i] and cell[i] and not scored[i]]
        # THE INVARIANT THAT MAKES STEP 3 SAFE. An unscored nucleus carries "not a doublet" because
        # nothing examined it, not because it was examined and found clean. That is harmless ONLY
        # while every unscored nucleus is removed by the count floor anyway - i.e. while the light
        # floor sits strictly below the quality floor. Lower the quality floor under the light
        # floor and unexamined nuclei enter the deliverable wearing a verdict nobody computed.
        out.append(Finding(
            "step 3: no unexamined nucleus reaches the deliverable",
            "REFUSE" if kept_unscored else "ok",
            (f"{len(kept_unscored):,} nuclei in the deliverable were never scored for doublets. "
             f"They carry 'not a doublet' because nothing looked, which is not the same claim"
             if kept_unscored else
             f"every kept nucleus was examined; the {sum(1 for i in range(n) if cell[i] and not scored[i]):,} "
             f"unexamined cells are all removed by other criteria")))

        if light_floor is not None and quality_floor is not None:
            bad = light_floor >= quality_floor
            out.append(Finding(
                "step 3: the light floor is strictly below the quality floor",
                "REFUSE" if bad else "ok",
                (f"light floor {light_floor} is NOT below quality floor {quality_floor}. The "
                 f"previous check may pass today and stop passing on the next cohort: this "
                 f"inequality is the only reason unexamined nuclei cannot survive"
                 if bad else
                 f"light floor {light_floor} < quality floor {quality_floor}, so an unexamined "
                 f"nucleus is removed by the count floor by construction")))

        # The doublet RATE must be reported over the examined set. Both denominators are computed
        # so the difference is visible rather than a footnote.
        if doublet_criterion and doublet_criterion in crit:
            d = crit[doublet_criterion]
            n_sc = sum(1 for i in range(n) if scored[i])
            n_cell = sum(1 for i in range(n) if cell[i])
            n_d = sum(1 for i in range(n) if d[i])
            r_sc = n_d / n_sc if n_sc else 0.0
            r_cell = n_d / n_cell if n_cell else 0.0
            out.append(Finding(
                "step 4: the doublet rate names its denominator", "REVIEW",
                f"{n_d:,} called: {100*r_sc:.2f}% of the {n_sc:,} EXAMINED, "
                f"{100*r_cell:.2f}% of all {n_cell:,} called cells. Published rates use the "
                f"first; dividing by all cells puts never-examined nuclei in the denominator and "
                f"deflates it. Quote the denominator with the number."))
            outside = [i for i in range(n) if d[i] and not scored[i]]
            out.append(Finding(
                "step 4: no doublet call outside the examined set",
                "REFUSE" if outside else "ok",
                (f"{len(outside):,} nuclei are called doublets but were never scored - the join "
                 f"assigned a call to a barcode the detector never saw"
                 if outside else "every doublet call belongs to a nucleus that was examined")))

    # --- E. the cluster-check object: superset, and differing by doublets ALONE --------------------
    if predoublet_keep is not None:
        pre = [bool(x) for x in predoublet_keep]
        not_superset = [i for i in range(n) if keep[i] and not pre[i]]
        extra = [i for i in range(n) if pre[i] and not keep[i]]
        out.append(Finding(
            "step 6 object is a superset of the deliverable",
            "REFUSE" if not_superset else "ok",
            (f"{len(not_superset):,} nuclei are in the deliverable but NOT in the cluster-check "
             f"object - it is not a superset, so it is a second filter nobody audited"
             if not_superset else
             f"contains all {sum(keep):,} delivered nuclei plus {len(extra):,} more")))
        if doublet_criterion and doublet_criterion in crit:
            d = crit[doublet_criterion]
            wrong = [i for i in extra if not d[i]]
            out.append(Finding(
                "step 6 object differs from the deliverable by doublets ALONE",
                "REFUSE" if wrong else "ok",
                (f"{len(wrong):,} of its {len(extra):,} extra nuclei are not doublet-flagged. The "
                 f"sanctioned exception is 'doublets retained'; anything else retained is an "
                 f"unaudited difference between the two objects"
                 if wrong else
                 f"all {len(extra):,} extra nuclei are doublet-flagged, as intended")))

    # --- F. every criterion is load-bearing, and overlaps are visible ------------------------------
    if crit:
        detail = []
        for c, m in crit.items():
            flagged = sum(1 for x in m if x)
            others = [k for k in crit if k != c]
            alone = sum(1 for i in range(n)
                        if m[i] and cell[i] and not any(crit[k][i] for k in others))
            detail.append(f"{c}: flagged {flagged:,}, removed ALONE {alone:,}")
        dead = [c for c in crit
                if sum(1 for i in range(n)
                       if crit[c][i] and cell[i]
                       and not any(crit[k][i] for k in crit if k != c)) == 0]
        out.append(Finding(
            "every applied criterion removes something no other does",
            "REVIEW" if dead else "ok",
            (f"{len(dead)} criterion(a) remove nothing another does not: {', '.join(dead)}. That "
             f"is not necessarily wrong, but a criterion with no unique contribution is not "
             f"protecting the deliverable and should not be described as if it were"
             if dead else "each criterion is solely responsible for some removal"),
            detail))
    return out


def verdict(findings) -> str:
    """"REFUSE", "REVIEW" or "PASS" - the closed vocabulary `Pipeline.gate` acts on.

    It returned "FAIL" and "ok", words the gate does not stop on, so a failed audit stopped
    nothing (harness autodev W043, GAP-119).
    """
    if any(f.severity == "REFUSE" for f in findings):
        return "REFUSE"
    return "REVIEW" if any(f.severity == "REVIEW" for f in findings) else "PASS"


def enforce(findings) -> None:
    """Raise on any REFUSE, outside the pipeline. REVIEW findings do not stop a run.

    Step 7 does not call this: it calls `integrity()`, then passes `verdict()` to `Pipeline.gate`,
    which stops the run and can lift only the floor ordering. This is for a caller with no pipeline.
    """
    bad = [f for f in findings if f.severity == "REFUSE"]
    if bad:
        raise AuditFailure(
            "the removal is not what the pipeline recorded:\n"
            + "\n".join(str(f) for f in bad))


def integrity(findings, accept=()) -> list:
    """Raise on a failed integrity check, whatever decisions.yml accepts; report an entry naming one.

    `accept` is decisions.yml's `accept` list; entries of other steps are not read. Any REFUSE
    finding whose check is not in ACCEPTABLE raises `AuditFailure`, naming the check, refusing by
    name each accept entry of this step that names an integrity check, and saying how to fix the
    cause. It prints no entry to paste, because no entry would count. When every integrity check
    held, returns one REVIEW finding per entry of this step that names one: the entry lifted nothing
    and would lift nothing when the check fails, and the report says so.
    """
    liftable = {_norm(c) for c in ACCEPTABLE}
    mine = [e for e in (accept or []) if isinstance(e, dict) and _norm(e.get("step")) == _norm(STEP)]
    failed = [f for f in findings if f.severity == "REFUSE" and _norm(f.check) not in liftable]
    named = {_norm(c) for c in INTEGRITY} | {_norm(f.check) for f in failed}
    refused = [e for e in mine if _norm(e.get("check")) in named]

    def entry(e) -> str:
        who = str(e.get("approved_by") or "").strip() or "no one named"
        return f"the accept entry for {str(e.get('check'))!r} (approved by {who}) is refused"

    if failed:
        also = [f.check for f in findings if f.severity == "REFUSE" and _norm(f.check) in liftable]
        raise AuditFailure(
            "the per-cell tables contradict themselves, so a deliverable written from them would "
            "not be the population its own ledger records:\n"
            + "\n".join(f"  {f}" for f in failed)
            + "".join(f"\n  ({entry(e)}: an integrity check is never accepted)" for e in refused)
            + "\n\nNo object was written. These are integrity checks, not gates: a table 07_measure "
              "writes cannot fail one, so a failure means a per-cell table was changed or stitched "
              "after 07_measure wrote it, and an accept entry does not lift one - accepting it would "
              "deliver an object that contradicts its own ledger. Fix the cause: remove each changed "
              "tables/<library>.percell.csv so that 07_measure measures the library again (or run "
              "with --force), and run again."
            + "".join(f"\n  Also refused, and acceptable by name once the tables hold: {c!r}."
                      for c in also))
    return [Finding(
        "an accept entry names an integrity check", "REVIEW",
        f"{entry(e)}: that check is an integrity check, which is never accepted - the entry "
        f"lifted nothing on this run, and would lift nothing if the check failed. Remove it "
        f"from the decisions file.") for e in refused]
