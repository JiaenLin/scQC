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
import math
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
    # THE CONFOUNDING BLOCK. Its inputs are the samplesheet and the per-cell tables, both of
    # which an apply-mode run has - so an absence here is a property of the DESIGN, and the note
    # says which property rather than naming a file.
    "F16": "the design. `collect()` was given no samplesheet rows, or the samplesheet carries "
           "no factor: `_design()` excludes a column that is constant, one with more levels "
           "than the libraries can support, and one with a single library per level.",
    "F17": "a confounded pair. No two design factors partition these libraries identically, so "
           "there is no arm to compare - which is a clean design, not a missing figure. F16 "
           "shows every pair and how it was classified.",
    "F18": "a confounded pair, as F17; or the per-cell tables, which measure mode does not "
           "write.",
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
    "F16": "Every library against every design factor the samplesheet carries, and then every "
           "pair of those factors compared exactly. Two factors that partition the libraries "
           "identically are ALIASED: no analysis of these data can say which of the two a "
           "difference belongs to. This is a property of the experiment and no step of this "
           "pipeline can change it.",
    "F17": "Each QC metric across the arms the aliasing leaves, over the nuclei the filter "
           "kept. The violin is the arm; the dots are the individual libraries' own medians "
           "and the amber band is the span they cover. An arm difference no larger than the "
           "spread between libraries of the same arm is a library effect, whatever it is "
           "labelled. The ratio printed on each panel is that comparison and is an estimate, "
           "not a test - a confounded factor cannot be tested.",
    "F18": "The share of each arm's barcodes the filter removed, in total and per criterion. "
           "This is the part of the confounding the pipeline itself creates: a criterion that "
           "takes a larger share out of one arm has made the arms differ for a reason that has "
           "nothing to do with the factor they are named after.",
    "F15": "Mitochondrial % per nucleus before and after the filter, over the barcodes above the "
           "light floor - the population the ceiling was derived over, because a percentage of a "
           "30-UMI droplet is not a measurement. The ceiling is PER LIBRARY, so each library "
           "carries its own segment and there is no cohort line to read across.",
}


def _rows(path: Path) -> list:
    with Path(path).open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


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


def _percell(tables: Path, samples) -> dict:
    return {s: _rows(tables / f"{s}.percell.csv") for s in samples}


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


# ---------------------------------------------------- F16-F18 · the confounding estimation
#
# WHAT THIS BLOCK IS FOR. Every other builder here assembles evidence about what the PIPELINE
# did. These three assemble evidence about what the EXPERIMENT is: which design factors cannot
# be told apart from each other in these libraries, how far apart the resulting arms sit on
# every QC metric this run measured, and whether the filter widened that gap.
#
# NOTHING HERE IS A TEST. A confounded factor cannot be tested - that is what confounded means -
# so no p-value is computed and none would mean anything. What is computed is exact: the
# partition each factor induces over the libraries, compared with every other factor's, and the
# distributions of the measured values within each arm.

#: At most this many values per arm reach a violin. A violin drawn from 20,000 nuclei and one
#: drawn from 200,000 are the same picture, and the payload is written to disk. The library
#: MEDIANS are always taken over every value, before any cap - a median of a subsample is a
#: different number and this figure's whole argument rests on those medians.
CONFOUNDING_MAX_POINTS = 20_000

#: The per-cell QC metrics, in the order the panels are drawn. `(column, label, log)`.
#: `complexity` is not a column: it is derived below from two that are, and is labelled as such.
#:
#: EVERY METRIC THE PER-CELL TABLE CARRIES IS HERE. A metric left out of this list is one whose
#: arm difference nobody will look at, and the arm difference is the entire question.
CONFOUNDING_METRICS = (
    ("total_counts", "UMI per nucleus", True),
    ("n_genes", "genes detected per nucleus", True),
    ("complexity", "complexity (log10 genes / log10 UMI)", False),
    ("pct_counts_mt", "mitochondrial % per nucleus", False),
    ("pct_counts_ribo", "ribosomal % per nucleus", False),
    ("nuclear_fraction", "nuclear fraction (intronic / intronic+exonic)", False),
    ("doublet_score", "doublet score", False),
)


def _partition(mapping: dict) -> dict:
    """`{level: frozenset(samples)}` for one factor's `{sample: level}`."""
    out: dict = {}
    for s, lv in mapping.items():
        if lv is None or not str(lv).strip():
            continue
        out.setdefault(str(lv), set()).add(str(s))
    return {k: frozenset(v) for k, v in out.items()}


def _relation(a_map: dict, b_map: dict) -> tuple:
    """How two factors sit against each other, by EXACT comparison of their partitions.

    Returns `(kind, detail)` with kind in `aliased` / `nested` / `crossed`.

    This is arithmetic on sets, not a statistic. Two factors are aliased when knowing one tells
    you the other and knowing the other tells you the one - their partitions of the libraries are
    identical, and no quantity of data separates them. One is nested in the other when the
    implication runs one way only. Anything else is crossed and can be estimated separately.

    Restricted to the libraries that carry BOTH values, and the count is returned in the detail:
    a pair compared over three of ten libraries is a different claim from one compared over all
    ten, and the figure prints which.
    """
    shared = sorted(set(a_map) & set(b_map))
    pairs = [(str(a_map[s]), str(b_map[s])) for s in shared
             if str(a_map[s]).strip() and str(b_map[s]).strip()]
    if not pairs:
        return "crossed", "no library carries a value for both"
    b_of_a: dict = {}
    a_of_b: dict = {}
    for a, b in pairs:
        b_of_a.setdefault(a, set()).add(b)
        a_of_b.setdefault(b, set()).add(a)
    b_from_a = all(len(v) == 1 for v in b_of_a.values())
    a_from_b = all(len(v) == 1 for v in a_of_b.values())
    n = len(pairs)
    if b_from_a and a_from_b:
        return "aliased", f"one-to-one over {n} librar{'y' if n == 1 else 'ies'}"
    if b_from_a:
        return "nested", f"the second is fixed by the first, over {n} libraries"
    if a_from_b:
        return "nested", f"the first is fixed by the second, over {n} libraries"
    return "crossed", f"both levels of each occur with both of the other, over {n} libraries"


def design_relations(sheet: list) -> dict:
    """The design, every pair of its factors compared, and the arms that survive the aliasing.

    Returns `{"factors", "levels", "relations", "arms", "aliased_group", "samples"}`.
    `arms` is None when no two factors are aliased - in which case there is nothing for F17 and
    F18 to compare, and saying so is the correct answer rather than inventing a contrast.

    THE ARM IS THE ALIASED GROUP'S OWN PARTITION. Where age, chemistry and batch all induce the
    same split, the arm is that split and its label carries all three names - because a reader
    looking at a difference between those arms must see, in the label itself, that it could be
    any of the three.
    """
    from engine.steps import _design  # noqa: PLC0415  (avoids a circular import at module load)

    levels = _design(sheet or [])
    factors = sorted(levels)
    samples = [str(r.get("sample")) for r in (sheet or []) if str(r.get("sample") or "").strip()]
    out = {"factors": factors, "levels": levels, "relations": [], "arms": None,
           "aliased_group": [], "samples": samples}
    if len(factors) < 2:
        return out

    edges = []
    for i, a in enumerate(factors):
        for b in factors[i + 1:]:
            kind, detail = _relation(levels[a], levels[b])
            out["relations"].append({"a": a, "b": b, "kind": kind, "detail": detail})
            if kind == "aliased":
                edges.append((a, b))
    if not edges:
        return out

    # The connected components of the aliased relation. Aliasing is transitive over exact
    # partitions, so a component is a set of factors that are all one factor as far as these
    # libraries can tell; the largest one is what the arms are built from.
    parent = {f: f for f in factors}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    groups: dict = {}
    for f in factors:
        groups.setdefault(find(f), []).append(f)
    group = max((g for g in groups.values() if len(g) > 1), key=lambda g: (len(g), g))
    group = sorted(group)
    out["aliased_group"] = group

    arms: dict = {}
    for s in (samples or sorted({x for f in group for x in levels[f]})):
        vals = [str(levels[f].get(s) or "") for f in group]
        if not all(v.strip() for v in vals):
            continue
        arms.setdefault("\n".join(vals), []).append(str(s))
    out["arms"] = {k: arms[k] for k in sorted(arms)} if len(arms) > 1 else None
    return out


def _arm_of_sample(arms: dict) -> dict:
    return {s: label for label, members in (arms or {}).items() for s in members}


def _ordered_samples(relations: dict) -> list:
    """Libraries ordered so an arm's rows are contiguous - F16 draws a bracket only if they are."""
    arms = relations.get("arms")
    known = relations.get("samples") or []
    if not arms:
        return sorted(known) if known else []
    out = []
    for label in sorted(arms):
        out.extend(sorted(arms[label]))
    out.extend(sorted(s for s in known if s not in set(out)))
    return out


def _f16(relations: dict) -> dict:
    if not relations.get("factors"):
        raise ValueError(
            "the samplesheet carries no design factor. `_design()` excludes a column that is "
            "constant, one with more levels than libraries can support, and one with a single "
            "library per level - so a cohort of one or two libraries has no factor by "
            "construction and there is nothing here to be confounded")
    return {"levels": relations["levels"], "samples": _ordered_samples(relations),
            "factors": relations["factors"], "relations": relations["relations"],
            "arms": relations.get("arms")}


def _gene_pattern_note(sheet: list, column: str, what: str) -> str:
    """`(prefix mt-)` for the axis label, or a statement that the libraries disagree.

    The gene class behind a percentage is decided by a pattern in the samplesheet, and the
    percentage on the page is only as good as that pattern. scQC's own docs/PRINCIPLES.md
    section 1 is about a ribosomal pattern that matched a ribosomal protein KINASE, so the
    pattern is printed beside the number rather than left for a reader to go and look up.
    """
    seen = sorted({str(r.get(column) or "").strip() for r in (sheet or [])
                   if str(r.get(column) or "").strip()})
    if not seen:
        return f" ({what} NOT STATED in the samplesheet)"
    if len(seen) == 1:
        return f" ({what} {seen[0]})"
    return f" ({what} DIFFERS across libraries: {', '.join(seen)})"


def _f17(percell: dict, relations: dict, sheet: list) -> dict:
    """Every QC metric the per-cell tables carry, per confounded arm, over the KEPT nuclei.

    THE KEPT POPULATION, deliberately. This figure is about what a confound will contaminate
    downstream, and downstream sees what the filter kept. What the filter itself did differently
    between the arms is F18's question and is not mixed into this one.
    """
    arms = relations.get("arms")
    if not arms:
        raise ValueError(
            "no two design factors are aliased in this samplesheet, so there is no confounded "
            "arm to compare. F16 shows every pair and how it was classified")
    of_sample = _arm_of_sample(arms)
    order = sorted(arms)

    values: dict = {key: {a: {} for a in order} for key, _lab, _log in CONFOUNDING_METRICS}
    for s, rows in percell.items():
        arm = of_sample.get(str(s))
        if arm is None:
            continue
        for key, _lab, _log in CONFOUNDING_METRICS:
            values[key][arm][s] = []
        for r in rows:
            if _flag(r.get("keep")) is not True:
                continue
            counts, genes = _num(r.get("total_counts")), _num(r.get("n_genes"))
            for key, _lab, _log in CONFOUNDING_METRICS:
                if key == "complexity":
                    # log10 genes / log10 UMI - the standard novelty score. Undefined at or
                    # below 10 UMI, where the denominator approaches 1 and the ratio explodes;
                    # such a nucleus contributes nothing rather than a large invented number.
                    v = (None if counts is None or genes is None or counts < 10 or genes < 1
                         else math.log10(genes) / math.log10(counts))
                else:
                    v = _num(r.get(key))
                if v is not None:
                    values[key][arm][s].append(float(v))

    metrics, absent = [], {}
    for key, label, log in CONFOUNDING_METRICS:
        if key == "pct_counts_mt":
            label += _gene_pattern_note(sheet, "mt_prefix", "prefix")
        elif key == "pct_counts_ribo":
            label += _gene_pattern_note(sheet, "ribo_pattern", "pattern")
        by_arm, medians, n_total = {}, {}, 0
        for arm in order:
            per_lib = values[key][arm]
            pooled = [v for s in sorted(per_lib) for v in per_lib[s]]
            n_total += len(pooled)
            if not pooled:
                continue
            step = max(1, -(-len(pooled) // CONFOUNDING_MAX_POINTS))
            by_arm[arm] = pooled[::step][:CONFOUNDING_MAX_POINTS]
            medians[arm] = {s: sorted(v)[len(v) // 2] for s, v in per_lib.items() if v}
        if not n_total:
            absent[label] = (f"the per-cell tables carry no value for `{key}`, so this metric "
                             f"was not measured by this run")
            continue
        metrics.append({"key": key, "label": label, "log": bool(log),
                        "by_arm": by_arm, "library_medians": medians})
    if not metrics:
        raise ValueError("no QC metric in the per-cell tables carried a value for any arm")
    return {"metrics": metrics, "arms": order, "absent": absent or None,
            "cap": CONFOUNDING_MAX_POINTS,
            "population": "the nuclei the filter KEPT (`keep` is true)"}


def _f18(percell: dict, relations: dict) -> dict:
    """What each criterion removed from each arm, over every barcode the arm held."""
    arms = relations.get("arms")
    if not arms:
        raise ValueError(
            "no two design factors are aliased in this samplesheet, so there is no confounded "
            "arm whose removal rates could be compared. F16 shows every pair")
    of_sample = _arm_of_sample(arms)
    order = sorted(arms)

    criteria: list = []
    for rows in percell.values():
        for r in rows[:1]:
            criteria = [c for c in r if str(c).startswith("fail_")]
        if criteria:
            break
    if not criteria:
        raise ValueError("no `fail_*` column in the per-cell tables, so no criterion can be "
                         "attributed to an arm")

    counts = {a: {"n_in": 0, "n_removed": 0, **{c: 0 for c in criteria}} for a in order}
    for s, rows in percell.items():
        arm = of_sample.get(str(s))
        if arm is None:
            continue
        bucket = counts[arm]
        for r in rows:
            bucket["n_in"] += 1
            if _flag(r.get("removed")) is True:
                bucket["n_removed"] += 1
            for c in criteria:
                if _flag(r.get(c)) is True:
                    bucket[c] += 1

    rates, overall = {}, {}
    for a in order:
        b = counts[a]
        n = b["n_in"]
        overall[a] = {"n_in": n, "n_removed": b["n_removed"],
                      "pct": (100.0 * b["n_removed"] / n) if n else None}
        rates[a] = {c: (100.0 * b[c] / n) if n else None for c in criteria}
    return {"rates": rates, "arms": order, "overall": overall, "criteria": criteria}


def collect(tables, *, samplesheet_rows=None, samplesheet=None,
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

    # THE DESIGN, computed once and shared by all three confounding figures so they cannot
    # disagree about which factors are aliased. A samplesheet that will not parse is caught here
    # rather than three times: `relations` degrades to an empty design and each `_try` below
    # records its own reason.
    try:
        relations = design_relations(sheet)
    except Exception as exc:                                            # noqa: BLE001
        relations = {"factors": [], "levels": {}, "relations": [], "arms": None,
                     "aliased_group": [], "samples": []}
        notes["F16"] = (f"the design could not be read from the samplesheet - "
                        f"{type(exc).__name__}: {exc}")
    _try("F16", "the run's samplesheet", lambda: _f16(relations))

    if samples:
        try:
            percell = _percell(tables, samples)
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
        # THE CONFOUNDING ESTIMATE. Both read the same per-cell tables every other applied-axis
        # figure reads, so a number on them is the number the criteria were evaluated on.
        _try("F17", src + " + the run's samplesheet",
             lambda: _f17(percell, relations, sheet))
        _try("F18", src + " + the run's samplesheet", lambda: _f18(percell, relations))

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
