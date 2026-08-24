"""Assemble `payload["figures"]` from the files a finished run left on disk.

Every figure in `report/figures.py` had a data contract and a drawing routine, and nothing ever
supplied the data: `payload["figures"]` was `{}` on every run, so a completed pipeline produced a
text-only report while the report itself dutifully printed "NOT PRODUCED" for every one. This
module is the missing half.

TWO RULES IT FOLLOWS, BOTH FROM THE PROJECT'S EVIDENCE RULE

1.  A figure is built from a NAMED FILE in `results/tables/`, and the file is recorded in the
    figure's `source`, so every mark on a plot can be traced to something a reader can open.
2.  A figure whose data the run did not produce is NOT drawn from a substitute. It is left out
    with a stated reason, and the report prints the reason in its place. A plausible figure drawn
    from the nearest available number is worse than no figure, because nothing on the page says
    it is not the thing it appears to be.

WHAT CANNOT BE BUILT FROM A FINISHED RUN, AND WHY

`UNAVAILABLE` holds the text the report prints where a figure would have been. It names the step
that would have to record something rather than merely reporting that a file is absent.

THE ENTRIES IT USED TO HOLD WERE MOSTLY WRONG, AND WRONGLY IN THE SAME DIRECTION

Four of the five said a step "computes something, uses it and discards it". Three of them did no
such thing:

  F1   `adapters/matrix.barcode_rank()` existed, `--rank-points` was plumbed end to end, and
       `run_summary_stats` COUNTED the pairs into `n_rank_points` and dropped the pairs. Step 0
       was never asked for them.
  F6   `_op_valley` has always written `<sample>.valley_density.csv`. It went to the scratch
       directory, which nothing publishes and this module never looks in.
  F5   `adapters/doublets.sweep()` existed, complete with its cross-setting version checks. No
       task called it.

Only F10/F11 were what the note claimed: no embedding was computed anywhere. All five are now
built, and the lesson is recorded because it is the expensive one - A PLAUSIBLE EXPLANATION FOR
A MISSING FIGURE READS EXACTLY LIKE A CORRECT ONE, and it is more durable than the defect,
because it tells the next reader not to look.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

# What each absent figure would need. The text reaches the reader, so it names what to do about
# it. Every entry here is REMOVED as soon as its builder succeeds - see `_try`.
UNAVAILABLE = {
    "F1": "the barcode-rank curve of each raw matrix. Step 0 measures it while verifying that "
          "the input is unfiltered and writes `tables/<sample>.barcode_rank.csv`; this run has "
          "none, which means either no library was ingested from a supplied matrix (a cohort "
          "rebuilt from FASTQ has no such curve) or the run predates the table.",
    "F5": "a sweep of the doublet prior. Step 4 scores each library at ONE dbr.sd, and one point "
          "cannot say whether the called rate tracks the prior or the data. The sweep is not run "
          "by default because it re-scores every library once per setting; request it with "
          "`--dbr-sd-sweep default,dbr,1` and this figure is drawn from "
          "`tables/doublet_sweep.csv`.",
    "F6": "the density curve the UMI valley was found on, `tables/<sample>.valley_density.csv`.",
    "F10": "a per-library embedding, `tables/<sample>.embedding.csv`, written by step 6.",
    "F11": "a per-library embedding, as F10.",
    "F13": "the density curve the gene valley was found on, "
           "`tables/<sample>.valley_density.csv`.",
}


# What each figure shows. A caption states what is drawn and what would be read off it; it does
# not interpret the result, which belongs to whoever reads the run.
CAPTIONS = {
    "F1": "Every library's barcode-rank curve, from the RAW matrix step 0 verified, on shared "
          "log-log axes. A curve that runs down into the single digits still holds its empty "
          "droplets; one that stops at a floor has already been cell-called, and the floor is "
          "usually a round number somebody chose. Neither is visible in the file's name or shape.",
    "F5": "The called doublet rate at each swept setting of the detector's prior, all libraries "
          "on one axis. A curve that follows the prior is a rate the prior set; one that stays "
          "put was measured from the library.",
    "F6": "Per library, the density the UMI valley was measured on, with that library's own "
          "valley and the single cohort floor drawn at the same place on every panel. A library "
          "the constant fits poorly shows as a cut sitting away from its own valley.",
    "F13": "The same as F6 for the GENE axis. Step 5 derives and applies two count floors, and "
           "each is only checkable against the density it came from.",
    "F10": "Where the doublet calls sit in each library's manifold. Doublets on the boundaries "
           "between clusters are behaving as doublets should; a coherent region called doublet "
           "is a population, and that is a different finding.",
    "F11": "The SAME coordinates as F10, coloured by whether any criterion removed the nucleus. "
           "The embedding is built over every barcode the denoiser called, so the nuclei the "
           "count floors and the ceiling removed are still in the picture - which is the only "
           "way to see whether a removal took a coherent region or scattered points.",
    "F2": "Ambient counts removed per library, the distribution across genes, and the removal "
          "rate split by each design factor. A rate that differs across an arm of the design is "
          "a technical removal that will read as biology downstream.",
    "F3": "Cells called by the aligner against cells called by the denoiser, per library. Points "
          "below the equality line are libraries where the denoiser called fewer cells than the "
          "aligner - a cell-selection decision made by a tool chosen for denoising.",
    "F4": "Per library, how many nuclei the doublet detector examined, how many sat below the "
          "light floor and were skipped for that stated reason, and how many were not examined "
          "for any recorded reason.",
    "F7": "Counts per nucleus before and after the filter, one pair of violins per library, on a "
          "log axis. The log axis is what shows two populations and a cut falling between them.",
    "F8": "Every cluster from step 6, positioned by its share of its library's UMI against its "
          "median mitochondrial percentage, marked by whether the cluster check flagged it. "
          "Clusters whose markers were never computed are drawn in the never-examined colour.",
    "F9": "What each criterion removed on its own, separated from what it removed jointly with "
          "another. A criterion whose removals are all shared changed nothing by itself.",
    "F12": "The same count distributions as F7 on a linear axis - the scale a reader works in - "
           "capped so the bulk is legible.",
    "F14": "F7's form on the GENE axis: genes per nucleus before and after the filter, with the "
           "cohort gene floor across every library. Step 7 applies this floor as well as the "
           "count floor, and a report showing only one of them shows part of the filter.",
    "F15": "Mitochondrial % per nucleus before and after the filter, over the barcodes above the "
           "light floor - the population the ceiling was derived over, because a percentage of a "
           "30-UMI droplet is not a measurement. The ceiling is PER LIBRARY, so each library "
           "carries its own segment and there is no cohort line to read across.",
}


def _rows(path: Path) -> list:
    with Path(path).open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


from engine.task import TaskFailure  # noqa: E402


def _num(value):
    """A number, or None for a blank. Never 0.0 for a blank - that is a measurement."""
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _flag(value):
    """A CSV cell as True/False/None. `"False"` is truthy as a string; this is why."""
    text = str(value).strip().lower()
    if text in ("true", "1"):
        return True
    if text in ("false", "0"):
        return False
    return None


def _mapping_or_empty(value) -> dict:
    return value if isinstance(value, dict) else {}


def _per_sample_tables(tables: Path, suffix: str) -> dict:
    """`{sample: rows}` for every `<sample><suffix>` in tables/, keyed by the file's OWN column.

    Keyed by the column and not by the filename, because the filename is a rendering of the
    sample and the column is the sample. A library whose name contains a dot, or one renamed
    between the run and the read, breaks the stem and does not break the column. The stem is used
    only when the file carries no `sample` column at all.
    """
    out: dict = {}
    for path in sorted(tables.glob(f"*{suffix}")):
        rows = _rows(path)
        if not rows:
            continue
        name = str(rows[0].get("sample") or "").strip() or path.name[: -len(suffix)]
        out.setdefault(name, []).extend(rows)
    return out


def _f1(tables: Path) -> dict:
    """The raw barcode-rank curve per library, and the aligner's cut where step 2 recorded one."""
    curves = _per_sample_tables(tables, ".barcode_rank.csv")
    if not curves:
        raise FileNotFoundError("no <sample>.barcode_rank.csv in tables/")
    libraries = {}
    for s, rows in curves.items():
        pairs = [(_num(r.get("rank")), _num(r.get("total_counts"))) for r in rows]
        pairs = [(r, t) for r, t in pairs if r is not None and t is not None]
        pairs.sort(key=lambda p: p[0])
        # `n_barcodes` is the matrix's own barcode count, not the number of points plotted: the
        # curve is downsampled log-uniformly, so the two differ by orders of magnitude and the
        # figure labels the panel with the one that describes the FILE.
        n_bc = next((_num(r.get("n_barcodes")) for r in rows
                     if _num(r.get("n_barcodes")) is not None), None)
        libraries[s] = {"ranks": [p[0] for p in pairs], "counts": [p[1] for p in pairs],
                        "n_barcodes": int(n_bc) if n_bc is not None else None}
    data = {"libraries": libraries}
    calls = tables / "cell_calls.csv"
    if calls.exists():
        called = {r["sample"]: _num(r.get("aligner")) for r in _rows(calls)}
        if any(v is not None for v in called.values()):
            data["called_cells"] = called
    return data


def _f5(tables: Path) -> dict:
    """The doublet-prior sweep: called rate against the swept dbr.sd, one line per library."""
    rows = _rows(tables / "doublet_sweep.csv")
    if not rows:
        raise ValueError("tables/doublet_sweep.csv is empty, so no setting was swept")
    sweep: dict = {}
    for r in rows:
        x, y = _num(r.get("dbr_sd_value")), _num(r.get("rate_over_scored"))
        if x is None or y is None:
            continue
        e = sweep.setdefault(str(r["sample"]), {"x": [], "y": []})
        e["x"].append(x)
        e["y"].append(y)
    for e in sweep.values():
        order = sorted(range(len(e["x"])), key=lambda i: e["x"][i])
        e["x"] = [e["x"][i] for i in order]
        e["y"] = [e["y"][i] for i in order]
    if not sweep:
        raise ValueError("no swept setting carried both a dbr.sd value and a rate. The token "
                         "'default' has no numeric value by design - it means the argument was "
                         "not passed - so a sweep of that alone cannot be drawn against an axis")
    applied = {_num(r.get("dbr_sd_applied")) for r in rows}
    applied.discard(None)
    data = {"sweep": sweep, "param_label": "dbr.sd (the detector's prior on the expected rate)",
            "rate_denominator": "of the nuclei actually scored; those below the light floor are "
                                "in no denominator"}
    if len(applied) == 1:
        data["chosen"] = applied.pop()
    return data


#: Which metric each density figure draws, and what the axis is called. The two ids exist because
#: `cut` is one cohort constant per figure and the two floors are in different units.
DENSITY_FIGURES = {"F6": ("umi", "UMI per nucleus", "umi_floor_proposed"),
                   "F13": ("genes", "genes per nucleus", "gene_floor_proposed")}


def _f6(tables: Path, thresholds: list, fig_id: str, metric: str, metric_label: str,
        cut_key: str, n_by_sample=None) -> dict:
    """One metric's density per library, with that library's valley and the cohort cut.

    The bounds come from `modules/05_quality`, imported rather than restated: they are the range
    outside which a derived floor is REFUSED, and a figure shading a different range from the one
    the code enforces would be reassuring about a guard that is not there.
    """
    curves = _per_sample_tables(tables, ".valley_density.csv")
    if not curves:
        raise FileNotFoundError("no <sample>.valley_density.csv in tables/")

    # `n` is the population the density was ESTIMATED OVER, which is every barcode in the object
    # - the valleys are measured before any floor, because a valley is a boundary between two
    # modes and pre-cutting deletes the one made of debris. `<sample>.percell.csv` covers exactly
    # that population, so its row count is the n; None when it was not read, never 0.
    n_by_sample = _mapping_or_empty(n_by_sample)
    densities = {}
    for s, rows in curves.items():
        pts = [(_num(r.get("grid")), _num(r.get("density"))) for r in rows
               if str(r.get("metric", "")).strip() == metric]
        pts = [(g, d) for g, d in pts if g is not None and d is not None]
        if pts:
            densities[s] = {"x": [p[0] for p in pts], "y": [p[1] for p in pts],
                            "n": n_by_sample.get(s)}
    if not densities:
        raise ValueError(f"no library recorded a density for metric {metric!r}")

    # The valley from the file written for it, at full precision - not the rounded copy in the
    # per-library threshold table, which exists to be read by a human.
    valleys = {}
    vfile = tables / f"valleys_{metric}.csv"
    if vfile.exists():
        valleys = {r["sample"]: _num(r.get("valley")) for r in _rows(vfile)}

    data = {"densities": densities, "metric_label": metric_label, "fig_id": fig_id}
    if valleys:
        data["valleys"] = valleys
    cuts = {_num(r.get(cut_key)) for r in thresholds}
    cuts.discard(None)
    if len(cuts) == 1:
        data["cut"] = cuts.pop()
    try:
        mod_dir = Path(__file__).resolve().parent.parent / "modules" / "05_quality"
        if str(mod_dir) not in sys.path:
            sys.path.insert(0, str(mod_dir))
        import quality  # noqa: PLC0415

        data["bounds"] = list(quality.UMI_BOUNDS if metric == "umi" else quality.GENE_BOUNDS)
    except Exception:                                                     # noqa: BLE001
        # A shaded margin nobody can attribute is worse than no margin. Absent is absent.
        pass
    return data


def _f10(tables: Path, percell: dict, *, colour_by: str, fig_id: str) -> dict:
    """The per-library embedding, with each barcode's doublet call and removal joined onto it.

    JOINED BY BARCODE, NEVER BY POSITION. The embedding covers the barcodes the denoiser called;
    the per-cell table covers every barcode the library held. The two files are written by
    different steps in different orders, and zipping them would colour the wrong points with no
    symptom on the page.

    A barcode with no per-cell row, and one whose flag the run never determined, both arrive as
    None. F10 draws those in their own colour and counts them out of the percentage, which is why
    this must not turn either of them into False.
    """
    embeddings = _per_sample_tables(tables, ".embedding.csv")
    if not embeddings:
        raise FileNotFoundError("no <sample>.embedding.csv in tables/")

    out = {}
    for s, rows in embeddings.items():
        by_barcode = {str(r.get("barcode")): r for r in (percell.get(s) or [])}
        xs, ys, marks = [], [], []
        for r in rows:
            x, y = _num(r.get("x")), _num(r.get("y"))
            if x is None or y is None:
                continue
            xs.append(x)
            ys.append(y)
            cell = by_barcode.get(str(r.get("barcode")))
            if cell is None:
                marks.append(None)
            elif colour_by == "doublet":
                # `doublet_class` is the call; `doublet_scored` says whether there was one to
                # make. A barcode below the light floor was never examined and is UNKNOWN - it is
                # not a singlet, and reading its blank class as one is the whole reason the
                # scored flag is written.
                scored = _flag(cell.get("doublet_scored"))
                cls = str(cell.get("doublet_class") or "").strip().lower()
                marks.append(None if scored is not True or not cls else cls == "doublet")
            else:
                marks.append(_flag(cell.get("removed")))
        if xs:
            out[s] = {"x": xs, "y": ys, colour_by: marks}
    if not out:
        raise ValueError("every embedding table was present and held no usable coordinate")
    return {"embeddings": out, "colour_by": colour_by, "fig_id": fig_id}


def _f2(tables: Path, sheet: list) -> dict:
    """How much ambient was removed, and whether evenly across the design."""
    rows = _rows(tables / "ambient_summary.csv")
    per_library = {r["sample"]: _num(r.get("fraction_removed_overall")) for r in rows}

    # The design panel is the one that can show a technical removal masquerading as biology, so
    # it is built from whatever factors the samplesheet actually carries rather than a fixed list.
    #
    # It said that before, and then used a fixed list: ("age", "diet", "chemistry", "batch") -
    # the four factor names of the cohort this pipeline was calibrated on. On that cohort the
    # panel was therefore correct, and on any other it rendered EMPTY while the comment above it
    # claimed otherwise. The gates in steps 1, 2 and 4 had always used the discovered factors, so
    # a run could refuse on a design differential the report then declined to draw.
    #
    # `_design()` is the one discovery routine, in `engine/steps.py`, and it is imported rather
    # than reimplemented here: it excludes constant columns, columns with too many levels, and -
    # the one that is easy to get wrong - identifier columns with a single sample per level.
    by_design = {}
    if sheet:
        from engine.steps import _design  # noqa: PLC0415  (avoids a circular import at module load)

        of_sample = {r.get("sample"): r for r in sheet}
        factors = list(_design(sheet))
        for factor in factors:
            levels = {}
            for s, frac in per_library.items():
                level = (of_sample.get(s) or {}).get(factor)
                if not level or frac is None:
                    continue
                levels.setdefault(str(level), []).append(float(frac))
            if levels:
                by_design[factor] = {
                    lv: {"rate": sum(v) / len(v), "n_samples": len(v),
                         "denominator": "mean of the per-library removed fractions"}
                    for lv, v in levels.items()}
    data = {"per_library": per_library}
    if by_design:
        data["by_design"] = by_design
    return data


def _f3(tables: Path) -> dict:
    """Aligner against denoiser cell calls. The CSV column is `denoiser`; F3 names it
    `cellbender`, and the rename happens here rather than in either of them."""
    calls = {}
    for r in _rows(tables / "cell_calls.csv"):
        calls[r["sample"]] = {"aligner": _num(r.get("aligner")),
                              "cellbender": _num(r.get("denoiser")),
                              "lost": _num(r.get("lost"))}
    return {"calls": calls}


def _percell(tables: Path, samples, objects=None, want=()) -> dict:
    """The per-cell rows, optionally completed from the objects the run itself wrote.

    A COLUMN THE TABLES LACK IS NOT A COLUMN THE RUN LACKS. `pct_counts_ribo` is computed into
    `obs` by the metrics step on every run and was written to the per-cell table only after that
    column was added, so a finished run holds the values while its CSV does not. Re-running the
    pipeline to recover a number it already computed is the expensive answer to the wrong
    question.

    OPT-IN and READ-ONLY: nothing is written back, the tables are not modified, and without
    `--objects` this behaves exactly as it did. `anndata` is imported only when the flag is used,
    so a report of a run whose columns are all present needs nothing extra installed.
    """
    per = {s: _rows(tables / f"{s}.percell.csv") for s in samples}
    missing = [c for c in want if not any(c in (r or {}) for rows in per.values() for r in rows[:1])]
    if not objects or not missing:
        return per

    import anndata as ad  # noqa: PLC0415  (only when the caller asked for this)

    found = sorted(Path(objects).glob("*.h5ad"))
    if not found:
        return per
    for path in found:
        try:
            obs = ad.read_h5ad(path, backed="r").obs
        except Exception:                                                   # noqa: BLE001
            continue
        have = [c for c in missing if c in obs.columns]
        if not have:
            continue
        # BY BARCODE, never by position. The per-cell table and the object are both ordered by
        # the run, and "both were written by the same run" is not a guarantee that a filter did
        # not reorder one of them.
        for col in have:
            lookup = dict(zip(obs.index.astype(str), obs[col].tolist()))
            for rows in per.values():
                for r in rows:
                    if col not in r:
                        v = lookup.get(str(r.get("barcode")))
                        r[col] = "" if v is None else v
    return per


def _f4(percell: dict, floors: dict) -> dict:
    """What was never examined - and the three categories are never merged.

    `scored` is a nucleus the detector examined. `below_floor` is one it did not examine for a
    STATED reason: it sat under the light floor. `never_scored` is one it did not examine for
    some other reason, and that is the category the figure exists to make visible.
    """
    coverage = {}
    for s, rows in percell.items():
        floor = floors.get(s)
        scored = below = never = 0
        for r in rows:
            if _flag(r.get("doublet_scored")):
                scored += 1
                continue
            umi = _num(r.get("total_counts"))
            if floor is not None and umi is not None and umi < floor:
                below += 1
            else:
                never += 1
        coverage[s] = {"scored": scored, "below_floor": below, "never_scored": never}
    return {"coverage": coverage}


def _f7(percell: dict, cut, *, scale: str, fig_id: str, column: str = "total_counts",
        metric_label: str = "UMI per nucleus", cut_word: str = "floor",
        only_at_or_above=None) -> dict:
    """What one applied criterion changed: each library's values before it, and what survived.

    ONE BUILDER FOR EVERY APPLIED AXIS, and the reason is the report's own spine: step 7 applies a
    count floor, a gene floor and a mitochondrial ceiling, and the spine lists all three. Drawing
    only the first left two rows of that table reading "no figure is produced for this axis" - the
    report showing a third of the filter while looking complete.

    `only_at_or_above` restricts the BEFORE population to barcodes at or above that `total_counts`.
    It exists for the mitochondrial axis and must not be used on the count axes. A percentage needs
    a denominator big enough to mean something - a 30-UMI droplet with 10 mitochondrial counts
    reads 33% - so the ceiling is DERIVED over the barcodes above the light floor, and a figure
    that drew it against every droplet would show the threshold sitting in a cloud it was never
    computed from. The count axes need the opposite: their whole point is that the debris mode is
    visible and the cut falls between it and the nuclei.
    """
    distributions = {}
    n_excluded = 0
    for s, rows in percell.items():
        before, after = [], []
        for r in rows:
            v = _num(r.get(column))
            if v is None:
                continue
            if only_at_or_above is not None:
                depth = _num(r.get("total_counts"))
                if depth is None or depth < float(only_at_or_above):
                    n_excluded += 1
                    continue
            before.append(v)
            if _flag(r.get("keep")):
                after.append(v)
        distributions[s] = {"before": before, "after": after}
    data = {"distributions": distributions, "metric_label": metric_label,
            "scale": scale, "fig_id": fig_id, "cut_word": cut_word}
    if cut is not None:
        data["cut"] = cut
    if scale == "linear" and only_at_or_above is None:
        # The linear panel is unreadable at full range - a handful of very deep nuclei set the
        # axis and the bulk collapses into the first pixel. The figure states that it is capped;
        # an uncapped-looking axis that had been capped would be the problem.
        #
        # NOT for a percentage: it is already bounded at 100, so there is no long tail to clip and
        # a cap would hide the high-mitochondrial population the ceiling exists to remove - the
        # one part of that figure a reader has come to look at.
        pooled = sorted(v for d in distributions.values() for v in d["before"])
        if pooled:
            data["cap"] = pooled[min(len(pooled) - 1, int(0.99 * len(pooled)))]
    return data


def _f8(tables: Path) -> dict:
    """The step-6 cluster profile, read through the module that knows its column types.

    Not `csv.DictReader` alone: every value would arrive as text, `"False"` is truthy, and the
    figure would draw every cluster as flagged. `read_profile_csv` is the reader that already
    exists for this and it is the reason the three-valued flags survive the round trip.
    """
    import sys

    root = Path(__file__).resolve().parent.parent
    mod_dir = root / "modules" / "06_cluster_check"
    if str(mod_dir) not in sys.path:
        sys.path.insert(0, str(mod_dir))
    import cluster_flags  # noqa: PLC0415

    clusters = cluster_flags.read_profile_csv(tables / "cluster_profile.csv")
    return {"clusters": clusters,
            "x_key": "umi_frac_of_sample",
            "x_label": "share of the library's UMI (%)"}


def _f9(tables: Path, n_in, n_removed) -> dict:
    """What each criterion removed UNIQUELY, from the ledger of removed observations."""
    from . import figures as figmod  # noqa: PLC0415

    path = tables / "removal_ledger.csv"
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        # The criteria are the ledger's own `fail_*` columns, so a criterion that fired on
        # nothing still appears - evidence it was evaluated, which a list of what fired cannot
        # give - and a criterion added later needs no change here.
        criteria = [c for c in (reader.fieldnames or []) if c.startswith("fail_")]
        rows = [(r.get("identifier"),
                 [c for c in criteria if _flag(r.get(c)) or str(r.get(c)).strip() == "1"])
                for r in reader]
    # criterion_contributions() returns {"per_criterion": {...}, "n_removed": n,
    # "n_multi_criterion": n} - the mapping is one key INSIDE it, not the return value. Passing
    # the whole thing through made F9 iterate "n_removed" as though it were a criterion and die
    # on `int.get`, which is what it looks like when a docstring describes the nested shape.
    contrib = figmod.criterion_contributions(rows, criteria=criteria)
    data = {"per_criterion": contrib["per_criterion"],
            "n_removed": n_removed if n_removed is not None else contrib["n_removed"]}
    if n_in is not None:
        data["n_in"] = n_in
    return data


def _f16(sheet):
    """{factor: {sample: level}} for every factor the samplesheet actually carries."""
    from engine.steps import _design  # noqa: PLC0415  (circular import at module load)
    factors = list(_design(sheet))
    if not factors:
        raise TaskFailure(
            "the samplesheet carries no design factor with two or more levels, so there is no "
            "pair of factors that could be confounded. That is a property of this cohort, not "
            "a failure.")
    of_sample = {r.get("sample"): r for r in sheet}
    out = {}
    for f in factors:
        levels = {s: (r or {}).get(f) for s, r in of_sample.items() if (r or {}).get(f)}
        if levels:
            out[f] = levels
    return {"design": out}


#: The per-library columns F17 draws, and why each one is readable across a confounded split:
#: every one is set by dissociation, chemistry or the cell call rather than by the condition
#: under test. NOT a fixed list of factor names - these are scQC's own metric columns, which
#: exist on every cohort.
_F17_METRICS = (("ribo_median", "ribosomal %, median"),
                ("mito_q3", "mitochondrial %, Q3"),
                ("doublet_rate_pct", "doublet rate %"),
                ("cells_lost", "cells lost at the call"))


def _median(values):
    """The middle value, or None. Local rather than imported from `figures`, which is the module
    that needs matplotlib; this one must stay importable without it."""
    v = sorted(x for x in values if x is not None)
    if not v:
        return None
    mid = len(v) // 2
    return v[mid] if len(v) % 2 else (v[mid - 1] + v[mid]) / 2.0


#: The per-cell columns F17 draws. Every one is set by dissociation, chemistry or sequencing
#: depth rather than by the condition under test, which is what makes them readable across a
#: split whose factors cannot be told apart. NOT a list of factor names - these are scQC's own
#: measurement columns and exist on every cohort.
_F17_COLUMNS = (("pct_counts_ribo", "ribosomal % per nucleus"),
                ("n_genes", "genes detected per nucleus"),
                ("pct_counts_mt", "mitochondrial % per nucleus"),
                ("total_counts", "UMI per nucleus"),
                ("nuclear_fraction", "nuclear fraction"))


def _confounded_pair(sheet):
    """The two factors whose levels determine each other for the most libraries, and by how much.

    FOUND, NOT NAMED. `_design()` discovers the factors - the same routine the gates use - and
    the pair is chosen by measurement. A fixed list of factor names is how the ambient design
    panel came to render empty on every cohort but the one it was written beside.
    """
    from engine.steps import _design  # noqa: PLC0415

    factors = list(_design(sheet))
    of_sample = {r.get("sample"): r for r in sheet}
    best, best_score = None, -1.0
    for i_, a in enumerate(factors):
        for b in factors[i_ + 1:]:
            both = [s for s in of_sample if of_sample[s].get(a) and of_sample[s].get(b)]
            if len(both) < 2:
                continue
            fwd = {}
            for s in both:
                fwd.setdefault(of_sample[s][a], set()).add(of_sample[s][b])
            score = sum(1 for s in both if len(fwd[of_sample[s][a]]) == 1) / len(both)
            if score > best_score:
                best, best_score = (a, b), score
    return best, best_score, factors, of_sample


def _arm_label(level, factor, partners, of_sample, samples):
    """`aged - V2 - batch A`: every factor that shares this partition, under one violin.

    An axis reading `aged` / `young` invites the reading that AGE is what differs. When three
    factors are the same split of the libraries, no honest axis can name one of them, and this
    is the whole reason the panel exists.
    """
    parts = [str(level)]
    for other in partners:
        vals = {str(of_sample[s].get(other)) for s in samples if of_sample[s].get(other)}
        if len(vals) == 1:
            parts.append(vals.pop())
    return "\n".join(parts)


def _f17(sheet, percell):
    """Per-cell distributions by arm, plus each library's median over them."""
    best, score, _factors, of_sample = _confounded_pair(sheet)
    if not best:
        raise TaskFailure("fewer than two design factors carry a level on two or more libraries")
    a, b = best
    arms = {}
    for s in percell:
        lv = (of_sample.get(s) or {}).get(a)
        if lv:
            arms.setdefault(str(lv), []).append(s)
    if len(arms) < 2:
        raise TaskFailure(f"{a!r} has fewer than two levels among the libraries with per-cell "
                          f"tables, so there is nothing to compare across")

    distributions, per_library, present = {}, {}, []
    for col, label in _F17_COLUMNS:
        by_sample, meds = {}, {}
        for s, rows in percell.items():
            vals = [_num(r.get(col)) for r in rows]
            vals = [v for v in vals if v is not None]
            if vals:
                by_sample[s] = vals
                meds[s] = _median(vals)
        if sum(len(v) for v in by_sample.values()) >= 10:
            distributions[label] = by_sample
            per_library[label] = meds
            present.append(col)
    if not distributions:
        raise TaskFailure(
            "no per-cell column in tables/<sample>.percell.csv carries values. Ribosomal content "
            "is written per cell only on runs made after that column was added.")

    partners = [f for f in (b,) if f]
    labels = {lv: _arm_label(lv, a, partners, of_sample, ss) for lv, ss in arms.items()}
    return {"distributions": distributions, "arms": arms, "arm_labels": labels,
            "per_library": per_library, "pair": f"{a} / {b}",
            "agreement": round(score, 3), "columns": present}


def _f18(sheet, percell):
    """What fraction of each arm the filter removed, and each library inside it."""
    best, _score, _factors, of_sample = _confounded_pair(sheet)
    if not best:
        raise TaskFailure("fewer than two design factors carry a level on two or more libraries")
    a, b = best
    arms = {}
    for s, rows in percell.items():
        lv = (of_sample.get(s) or {}).get(a)
        if not lv:
            continue
        removed = sum(1 for r in rows if _flag(r.get("removed")))
        entry = arms.setdefault(str(lv), {"removed": 0, "n": 0, "libraries": {}})
        entry["removed"] += removed
        entry["n"] += len(rows)
        entry["libraries"][s] = (removed, len(rows))
    if len(arms) < 2:
        raise TaskFailure(f"{a!r} has fewer than two levels among the libraries with per-cell "
                          f"tables, so removal cannot be compared across arms")
    return {"arms": arms, "pair": f"{a} / {b}"}


def collect(tables, *, samplesheet_rows=None, samplesheet=None, objects=None,
            n_in=None, n_removed=None) -> tuple:
    """Every figure this run's files can support, plus a note for every one they cannot.

    Returns `(figures, notes)`. `figures` is the `payload["figures"]` dict - `{id: {"data":
    {...}, "source": "..."}}`. `notes` maps a figure id to why it is absent, for the report to
    print where the figure would have been.

    A builder that raises is caught and becomes a note. The report is the only record of the run
    and a figure that could not be assembled must not be able to prevent it being written - the
    same reasoning `render_figures` applies one layer down.
    """
    tables = Path(tables)
    notes = dict(UNAVAILABLE)
    figures = {}

    # The rows the run was built from, when the caller has them in hand; the file only when it
    # does not. A rebuild from a finished directory has no Pipeline to ask, which is the case
    # the path argument exists for.
    sheet = [dict(r) for r in (samplesheet_rows or [])]
    if not sheet and samplesheet and Path(samplesheet).exists():
        text = Path(samplesheet).read_text(encoding="utf-8")
        delim = "\t" if "\t" in text.splitlines()[0] else ","
        sheet = [r for r in csv.DictReader(
            [ln for ln in text.splitlines() if not ln.startswith("#")], delimiter=delim)]

    thresholds, floors = [], {}
    umi_cut = gene_cut = light_floor = None
    mito_cuts: dict = {}
    tp = tables / "thresholds_per_sample.csv"
    if tp.exists():
        # Row 1 of this table is a SCOPE row ("per library" / "cohort constant"), not a sample.
        thresholds = [r for r in _rows(tp) if str(r.get("sample", "")).strip() != "scope"]
        floors = {r["sample"]: _num(r.get("light_floor_umi")) for r in thresholds}

        def _cohort(key):
            """The one value a cohort-constant column holds, or None if it holds several.

            None rather than the first or the mean: a column that was supposed to be one
            constant and is not is a finding, and a figure drawing one of its values as though
            it were the constant would state a rule the run did not apply."""
            seen = {_num(r.get(key)) for r in thresholds}
            seen.discard(None)
            return seen.pop() if len(seen) == 1 else None

        umi_cut = _cohort("umi_floor_proposed")
        gene_cut = _cohort("gene_floor_proposed")
        light_floor = _cohort("light_floor_umi")
        # PER LIBRARY by design - kept as a mapping and never collapsed. See `_violin_row`.
        for r in thresholds:
            v = _num(r.get("mito_ceiling_pct"))
            if v is not None:
                mito_cuts[r["sample"]] = v

    samples = [r["sample"] for r in thresholds] or \
              sorted(p.name.split(".percell.csv")[0] for p in tables.glob("*.percell.csv"))

    percell = {}

    def _try(fid, source, build):
        try:
            data = build()
        except FileNotFoundError as exc:
            # THE STANDING NOTE SURVIVES A MISSING FILE, and this is not a detail.
            #
            # Every builder here opens a file, so the ordinary way for a figure to be absent is
            # FileNotFoundError - and overwriting the note with the exception replaced the one
            # sentence that tells a reader what to do ("request it with --dbr-sd-sweep") with a
            # path and an errno, which tells them the pipeline is broken. A file that was never
            # written is the case `UNAVAILABLE` exists to describe.
            #
            # A figure with no standing note has nothing better to say, so it keeps the error.
            notes.setdefault(fid, f"could not be assembled - {type(exc).__name__}: {exc}")
            return
        except Exception as exc:                                            # noqa: BLE001
            # Anything else IS a defect - a column that changed name, a value that will not
            # parse - and it replaces the standing note, because the standing note would now be
            # a wrong explanation of a real fault.
            notes[fid] = f"could not be assembled - {type(exc).__name__}: {exc}"
            return
        figures[fid] = {"data": data, "source": source, "caption": CAPTIONS.get(fid, "")}
        notes.pop(fid, None)

    _try("F1", "tables/<sample>.barcode_rank.csv + tables/cell_calls.csv", lambda: _f1(tables))
    _try("F2", "tables/ambient_summary.csv", lambda: _f2(tables, sheet))
    _try("F3", "tables/cell_calls.csv", lambda: _f3(tables))
    _try("F5", "tables/doublet_sweep.csv", lambda: _f5(tables))
    _try("F8", "tables/cluster_profile.csv", lambda: _f8(tables))

    if samples:
        try:
            percell = _percell(tables, samples, objects=objects,
                               want=[c for c, _l in _F17_COLUMNS])
        except Exception as exc:                                            # noqa: BLE001
            for fid in ("F4", "F7", "F12"):
                notes[fid] = (f"needs the per-cell tables, which could not be read - "
                              f"{type(exc).__name__}: {exc}")
            for fid in ("F10", "F11"):
                notes[fid] = (f"the embedding has no flags to colour without the per-cell "
                              f"tables, which could not be read - {type(exc).__name__}: {exc}")
    if percell:
        src = "tables/<sample>.percell.csv"
        _try("F4", src + " + tables/thresholds_per_sample.csv",
             lambda: _f4(percell, floors))
        _try("F7", src, lambda: _f7(percell, umi_cut, scale="log", fig_id="F7"))
        _try("F12", src, lambda: _f7(percell, umi_cut, scale="linear", fig_id="F12"))
        # THE OTHER TWO APPLIED AXES. Same builder, same form, different column - so a reader
        # comparing the three is comparing the filter and not three chart styles.
        _try("F14", src + " + tables/thresholds_per_sample.csv",
             lambda: _f7(percell, gene_cut, scale="log", fig_id="F14", column="n_genes",
                         metric_label="genes detected per nucleus"))
        # PER LIBRARY, and passed as a mapping so the figure draws ten segments rather than
        # asserting a cohort constant that was never applied. Restricted to the barcodes above
        # the light floor, which is the population the ceiling was DERIVED over.
        _try("F15", src + " + tables/thresholds_per_sample.csv",
             lambda: _f7(percell, mito_cuts or None, scale="linear", fig_id="F15",
                         column="pct_counts_mt", metric_label="mitochondrial % per nucleus",
                         cut_word="ceiling", only_at_or_above=light_floor))
        # F10 and F11 are the SAME coordinates under two colourings, and that is the point of
        # storing an embedding rather than recomputing one: a reader comparing them is looking at
        # the flags, never at the projection.
        emb_src = "tables/<sample>.embedding.csv + " + src
        _try("F10", emb_src, lambda: _f10(tables, percell, colour_by="doublet", fig_id="F10"))
        _try("F11", emb_src, lambda: _f10(tables, percell, colour_by="removed", fig_id="F11"))

    # THE CONFOUNDING PAIR, from the samplesheet and the per-library table this report already
    # reads. F16 is the design itself; F17 is how far apart its arms sit on quantities the
    # condition under test did not set. Both come from files a finished run has, so a reader can
    # get them from `scqc report` without the pipeline running again.
    #
    # `_design()` DISCOVERS the factors - the same routine the gates use. A fixed list of factor
    # names is how the ambient design panel came to render empty on every cohort but the one it
    # was written beside, and this is the same trap one figure over.
    if sheet:
        _try("F16", "the samplesheet", lambda: _f16(sheet))
    if sheet and percell:
        src16 = "the samplesheet + tables/<sample>.percell.csv"
        _try("F17", src16, lambda: _f17(sheet, percell))
        _try("F18", src16, lambda: _f18(sheet, percell))

    # After the per-cell tables, because that is where the denominator comes from: F9 draws each
    # criterion's unique removals against the population they were removed FROM, and a bar chart
    # with no denominator cannot say whether a criterion removed much or little.
    n_considered = sum(len(r) for r in percell.values()) if percell else n_in
    _try("F9", "tables/removal_ledger.csv", lambda: _f9(tables, n_considered, n_removed))

    # Also after them, for the same reason in a different place: the density figures label each
    # panel with the population the estimate was made over, and the per-cell table is what counts
    # it. The figures are still drawn without it - the n reads NOT SUPPLIED and the curve is
    # unaffected - so this is an ordering preference, not a dependency.
    n_by_sample = {s: len(rows) for s, rows in percell.items()} if percell else {}
    for fid, (metric, label, cut_key) in DENSITY_FIGURES.items():
        _try(fid, f"tables/<sample>.valley_density.csv + tables/valleys_{metric}.csv",
             lambda f=fid, m=metric, lb=label, ck=cut_key:
                 _f6(tables, thresholds, f, m, lb, ck, n_by_sample))

    return figures, notes
