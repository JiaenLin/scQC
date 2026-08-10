"""Assemble `payload["figures"]` from the files a finished run left on disk.

Every figure in `report/figures.py` had a data contract and a drawing routine, and nothing ever
supplied the data: `payload["figures"]` was `{}` on every run, so a completed pipeline produced a
text-only report while the report itself dutifully printed "NOT PRODUCED" twelve times. This
module is the missing half.

TWO RULES IT FOLLOWS, BOTH FROM THE PROJECT'S EVIDENCE RULE

1.  A figure is built from a NAMED FILE in `results/tables/`, and the file is recorded in the
    figure's `source`, so every mark on a plot can be traced to something a reader can open.
2.  A figure whose data the run did not produce is NOT drawn from a substitute. It is left out
    with a stated reason, and the report prints the reason in its place. A plausible figure drawn
    from the nearest available number is worse than no figure, because nothing on the page says
    it is not the thing it appears to be.

WHAT CANNOT BE BUILT FROM A FINISHED RUN, AND WHY

Four of the twelve need data no step currently writes down. They are listed in `UNAVAILABLE`
with what each would need, and that text is what the report shows. They are not defects in this
module; they are steps that compute something, use it, and discard it.
"""
from __future__ import annotations

import csv
from pathlib import Path

# What each absent figure would need. The text reaches the reader, so it names the step that
# would have to record something, not merely the fact that a file is missing.
UNAVAILABLE = {
    "F1": "the barcode-rank curve of each RAW matrix. Step 0 reads those matrices to verify "
          "they are unfiltered but records only the verdict, so the curve behind it is gone by "
          "the time the report is built. Recording rank/count pairs during ingest would make "
          "this figure free.",
    "F5": "a sweep of the doublet prior. Step 4 scores each library at ONE prior (dbr), which "
          "is a single point; the figure asks whether the called rate tracks the prior across "
          "its range, and one point cannot answer it. It needs a sweep step this pipeline does "
          "not run.",
    "F6": "the density curve each valley was found on. Step 5 fits a KDE, takes the minimum, "
          "records the valley position, and discards the curve.",
    "F10": "a per-library embedding. Step 6 clusters each library on a neighbour graph but "
           "computes no 2-D embedding and stores no coordinates.",
    "F11": "a per-library embedding, as F10.",
}


# What each figure shows. A caption states what is drawn and what would be read off it; it does
# not interpret the result, which belongs to whoever reads the run.
CAPTIONS = {
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


def _f2(tables: Path, sheet: list) -> dict:
    """How much ambient was removed, and whether evenly across the design."""
    rows = _rows(tables / "ambient_summary.csv")
    per_library = {r["sample"]: _num(r.get("fraction_removed_overall")) for r in rows}

    # The design panel is the one that can show a technical removal masquerading as biology, so
    # it is built from whatever factors the samplesheet actually carries rather than a fixed list.
    by_design = {}
    if sheet:
        of_sample = {r.get("sample"): r for r in sheet}
        factors = [c for c in ("age", "diet", "chemistry", "batch") if c in (sheet[0] or {})]
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


def _f7(percell: dict, cut, *, scale: str, fig_id: str) -> dict:
    """What the cut changed: every library's counts before it and what survived it."""
    distributions = {}
    for s, rows in percell.items():
        before, after = [], []
        for r in rows:
            umi = _num(r.get("total_counts"))
            if umi is None:
                continue
            before.append(umi)
            if _flag(r.get("keep")):
                after.append(umi)
        distributions[s] = {"before": before, "after": after}
    data = {"distributions": distributions, "metric_label": "UMI per nucleus",
            "scale": scale, "fig_id": fig_id}
    if cut is not None:
        data["cut"] = cut
    if scale == "linear":
        # The linear panel is unreadable at full range - a handful of very deep nuclei set the
        # axis and the bulk collapses into the first pixel. The figure states that it is capped;
        # an uncapped-looking axis that had been capped would be the problem.
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

    thresholds, floors, umi_cut = [], {}, None
    tp = tables / "thresholds_per_sample.csv"
    if tp.exists():
        # Row 1 of this table is a SCOPE row ("per library" / "cohort constant"), not a sample.
        thresholds = [r for r in _rows(tp) if str(r.get("sample", "")).strip() != "scope"]
        floors = {r["sample"]: _num(r.get("light_floor_umi")) for r in thresholds}
        cuts = {_num(r.get("umi_floor_proposed")) for r in thresholds}
        cuts.discard(None)
        umi_cut = cuts.pop() if len(cuts) == 1 else None

    samples = [r["sample"] for r in thresholds] or \
              sorted(p.name.split(".percell.csv")[0] for p in tables.glob("*.percell.csv"))

    percell = {}

    def _try(fid, source, build):
        try:
            data = build()
        except Exception as exc:                                            # noqa: BLE001
            notes[fid] = f"could not be assembled - {type(exc).__name__}: {exc}"
            return
        figures[fid] = {"data": data, "source": source, "caption": CAPTIONS.get(fid, "")}
        notes.pop(fid, None)

    _try("F2", "tables/ambient_summary.csv", lambda: _f2(tables, sheet))
    _try("F3", "tables/cell_calls.csv", lambda: _f3(tables))
    _try("F8", "tables/cluster_profile.csv", lambda: _f8(tables))

    if samples:
        try:
            percell = _percell(tables, samples)
        except Exception as exc:                                            # noqa: BLE001
            for fid in ("F4", "F7", "F12"):
                notes[fid] = (f"needs the per-cell tables, which could not be read - "
                              f"{type(exc).__name__}: {exc}")
    if percell:
        src = "tables/<sample>.percell.csv"
        _try("F4", src + " + tables/thresholds_per_sample.csv",
             lambda: _f4(percell, floors))
        _try("F7", src, lambda: _f7(percell, umi_cut, scale="log", fig_id="F7"))
        _try("F12", src, lambda: _f7(percell, umi_cut, scale="linear", fig_id="F12"))

    # After the per-cell tables, because that is where the denominator comes from: F9 draws each
    # criterion's unique removals against the population they were removed FROM, and a bar chart
    # with no denominator cannot say whether a criterion removed much or little.
    n_considered = sum(len(r) for r in percell.values()) if percell else n_in
    _try("F9", "tables/removal_ledger.csv", lambda: _f9(tables, n_considered, n_removed))

    return figures, notes
