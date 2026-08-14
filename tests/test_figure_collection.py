"""Every figure the report declares must be assembled from a finished run's tables, and DRAWN.

WHAT THIS GUARDS, AND WHY IT IS TWO CHECKS AND NOT ONE

`report/collect.py` builds `payload["figures"]` as `{id: {"data": {...}}}` and `render_figures`
calls `FIGURE_FUNCTIONS[id](**data)`. So a builder can succeed, produce a dict, be recorded as
"assembled", and still be unrenderable - because one key it emits is not a parameter of the
function that receives it. `collect()` returning thirteen entries is therefore NOT evidence that
thirteen figures exist. This test renders every one of them.

It is the same lesson KNOWN_ISSUES records about step 6: the broken attempt SUCCEEDED. The
acceptance test for a figure is a figure.

WHY THE FIXTURE IS WRITTEN OUT AS FILES

Because that is the contract. `collect()` takes a directory and nothing else - a rebuild by
`scqc report <results>` has no Pipeline to ask - so a test that handed it objects in memory would
be testing a route no run takes. Every file below is written with the columns the producing step
writes, and the test fails if a builder needs a column no step produces.

Run: python tests/test_figure_collection.py
"""
from __future__ import annotations

import csv
import math
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from report.collect import DENSITY_FIGURES, UNAVAILABLE, collect  # noqa: E402

PASS, FAIL = [], []
SAMPLES = ("libA", "libB")


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}"
          + (f"  --  {detail}" if detail and not cond else ""))


def write(path: Path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def build_tables(tables: Path) -> None:
    """One finished run's tables/, as the steps write them."""
    tables.mkdir(parents=True, exist_ok=True)

    # The mitochondrial ceiling DIFFERS between the libraries, deliberately: it is the one
    # threshold this pipeline derives per library, and a fixture giving both the same value would
    # pass whether or not F15 collapsed ten rules into one cohort line.
    write(tables / "thresholds_per_sample.csv",
          ["sample", "light_floor_umi", "umi_valley", "umi_floor_proposed",
           "gene_valley", "gene_floor_proposed", "mito_ceiling_pct"],
          [["scope", "cohort constant", "per library", "cohort constant",
            "per library", "cohort constant", "per library"]]
          + [[s, 200, 340, 350, 250, 260, ceil]
             for s, ceil in zip(SAMPLES, (12.5, 19.75))])

    write(tables / "ambient_summary.csv", ["sample", "fraction_removed_overall"],
          [[s, 0.11] for s in SAMPLES])
    write(tables / "cell_calls.csv", ["sample", "aligner", "denoiser", "lost"],
          [[s, 900, 1000, 0] for s in SAMPLES])

    # --- F1: the raw barcode-rank curve, one file per library (step 0)
    for s in SAMPLES:
        write(tables / f"{s}.barcode_rank.csv",
              ["sample", "rank", "total_counts", "n_barcodes"],
              [[s, r, max(1.0, 20000.0 / r), 5000] for r in (1, 5, 25, 120, 600, 3000, 5000)])

    # --- F6/F13: the density each valley was measured on, both metrics (step 5)
    for s in SAMPLES:
        rows = []
        for metric, centre in (("umi", 2.6), ("genes", 2.4)):
            for i in range(60):
                g10 = 1.3 + i * 0.035
                rows.append([s, metric, 10 ** g10, g10,
                             math.exp(-((g10 - centre) ** 2) / 0.02), False, False])
        write(tables / f"{s}.valley_density.csv",
              ["sample", "metric", "grid", "grid_log10", "density", "is_valley", "is_mode"], rows)
    write(tables / "valleys_umi.csv", ["sample", "metric", "valley", "bimodal"],
          [[s, "umi", 340.5, True] for s in SAMPLES])
    write(tables / "valleys_genes.csv", ["sample", "metric", "valley", "bimodal"],
          [[s, "genes", 250.5, True] for s in SAMPLES])

    # --- F5: the dbr.sd sweep (step 4b, only when requested)
    write(tables / "doublet_sweep.csv",
          ["sample", "setting", "dbr_sd_value", "n_scored", "n_called", "rate_over_scored",
           "dbr_sd_applied"],
          [[s, lab, val, 900, int(900 * rate), rate, 0.06]
           for s in SAMPLES
           for lab, val, rate in (("dbr", 0.06, 0.11), ("1", 1.0, 0.24))])

    # --- F8: the cluster profile (step 6)
    write(tables / "cluster_profile.csv",
          ["sample", "cluster", "n", "median_umi", "umi_frac_of_sample", "median_pct_mt",
           "pct_uninformative", "pct_doublet", "A", "B", "C", "D", "FLAG", "WATCH"],
          [[s, c, 100, 900, 12.5, 4.0, 1.0, 3.0, "False", "False", "False", "False",
            "False", "False"] for s in SAMPLES for c in (0, 1)])

    # --- F4/F7/F12/F10/F11: the per-cell record, and the coordinates joined onto it
    # EVERY COLUMN AN APPLIED CRITERION IS MEASURED ON. Step 7 applies a count floor, a gene
    # floor and a mitochondrial ceiling, so a fixture carrying only `total_counts` lets F14 and
    # F15 "render" as blank panels and report PASS - a probe that cannot fail.
    for s in SAMPLES:
        cells = []
        emb = []
        for i in range(300):
            bc = f"{s}_{i:04d}"
            scored = i % 4 != 0                      # a quarter below the light floor: UNKNOWN
            doublet = scored and i % 9 == 0
            removed = doublet or i % 11 == 0
            cells.append([bc, s, 120 + 7 * i, 40 + 3 * (i % 210), 1.5 + (i % 37),
                          scored, "doublet" if doublet else ("singlet" if scored else ""),
                          removed, not removed])
            emb.append([bc, s, math.cos(i / 7.0) * 4, math.sin(i / 5.0) * 4, True, i % 2])
        write(tables / f"{s}.percell.csv",
              ["barcode", "sample", "total_counts", "n_genes", "pct_counts_mt",
               "doublet_scored", "doublet_class", "removed", "keep"], cells)
        write(tables / f"{s}.embedding.csv",
              ["barcode", "sample", "x", "y", "clustered", "cluster"], emb)

    # --- F9: the ledger of what left, and which criterion took it. EVERY row carries at least
    # one criterion, because the builder refuses a removal with no recorded reason - correctly,
    # and a fixture that trips that refusal is testing the fixture.
    ledger = []
    for s in SAMPLES:
        for n, i in enumerate(range(0, 300, 11)):
            fails = [n % 3 == 0, n % 5 == 0, n % 4 == 0]
            if not any(fails):
                fails[0] = True
            ledger.append([f"{s}_{i:04d}", *fails])
    write(tables / "removal_ledger.csv",
          ["identifier", "fail_umi_floor", "fail_mito_ceiling", "fail_doublet"], ledger)


print("\nassembling every figure from a finished run's tables/")
with tempfile.TemporaryDirectory() as tmp:
    tables = Path(tmp) / "tables"
    build_tables(tables)
    sheet = [{"sample": s, "age": a, "diet": d}
             for s, a, d in (("libA", "young", "chow"), ("libB", "aged", "hfd"))]
    figures, notes = collect(tables, samplesheet_rows=sheet)

    # The five that had a "not produced" note are exactly the ones this change is about, so they
    # are named individually: a count would pass while the wrong five were present.
    for fid in ("F1", "F5", "F6", "F10", "F11", "F13", "F14", "F15"):
        check(f"{fid} is assembled", fid in figures,
              notes.get(fid, "absent with no note, which is worse"))
    for fid in ("F2", "F3", "F4", "F7", "F8", "F9", "F12"):
        check(f"{fid} still assembles", fid in figures, notes.get(fid, ""))
    check("no figure is left with a note", not notes, str(sorted(notes)))

    print("\nthe applied axes carry values, not empty panels")
    # A figure whose distributions are all empty still RENDERS - as a blank panel saying so - and
    # would report PASS below. Each applied axis is checked for values before that happens.
    for fid, axis in (("F7", "counts"), ("F14", "genes"), ("F15", "mitochondrial %")):
        dists = (figures.get(fid, {}).get("data") or {}).get("distributions") or {}
        n_before = sum(len(v.get("before") or []) for v in dists.values())
        n_after = sum(len(v.get("after") or []) for v in dists.values())
        check(f"{fid} ({axis}) has values before and after the cut",
              n_before > 0 and n_after > 0, f"before={n_before} after={n_after}")

    # F15's cut is the one threshold this pipeline derives PER LIBRARY. Handed over as a scalar it
    # would draw one line across every panel, asserting a cohort constant that was never applied.
    f15_cut = (figures.get("F15", {}).get("data") or {}).get("cut")
    check("F15's ceiling reaches the figure as a per-library mapping",
          isinstance(f15_cut, dict) and len(f15_cut) == len(SAMPLES), str(f15_cut))
    check("...and the libraries genuinely differ, so a collapse would be visible",
          isinstance(f15_cut, dict) and len(set(f15_cut.values())) > 1, str(f15_cut))
    check("F15 says it is a ceiling, not a floor",
          (figures.get("F15", {}).get("data") or {}).get("cut_word") == "ceiling")
    # The mitochondrial axis is measured over the barcodes above the light floor - the population
    # the ceiling was derived over - so it must hold FEWER values than the count axis.
    f7_n = sum(len(v.get("before") or [])
               for v in ((figures.get("F7", {}).get("data") or {}).get("distributions")
                         or {}).values())
    f15_n = sum(len(v.get("before") or [])
                for v in ((figures.get("F15", {}).get("data") or {}).get("distributions")
                          or {}).values())
    check("F15 is restricted to the population its ceiling was derived over",
          0 < f15_n < f7_n, f"F15={f15_n} F7={f7_n}")

    print("\nthe data each builder emits is what its function accepts")
    # Bound OUTSIDE the guard below: report.figures imports matplotlib inside a function, not at
    # module level, so this import cannot raise ModuleNotFoundError - while the RENDERING does
    # need matplotlib. Binding it inside the try left `figmod` undefined on the skip path, and
    # the drawable check further down then raised NameError instead of skipping. A suite that
    # crashes where it means to skip cannot tell a reader which of the two happened.
    from report import figures as figmod

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for fid in sorted(figures):
            fn = figmod.FIGURE_FUNCTIONS.get(fid)
            if fn is None:
                check(f"{fid} has a drawing function", False, "not in FIGURE_FUNCTIONS")
                continue
            try:
                fig = fn(**figures[fid]["data"])
                check(f"{fid} renders", fig is not None)
                plt.close(fig)
            except Exception as e:                                        # noqa: BLE001
                check(f"{fid} renders", False, f"{type(e).__name__}: {e}")
    except ModuleNotFoundError as e:                                      # noqa: BLE001
        print(f"  SKIP rendering: needs {e.name}")

    print("\nevery id the report positions has a caption and an availability note")
    from report.build import FIGURE_QUESTIONS

    for fid in sorted(FIGURE_QUESTIONS):
        check(f"{fid} is drawable", fid in figmod.FIGURE_FUNCTIONS)
    check("F13 is declared by a step", "F13" in FIGURE_QUESTIONS)
    check("both density ids are collected", set(DENSITY_FIGURES) == {"F6", "F13"},
          str(sorted(DENSITY_FIGURES)))

print("\nabsence is still reported as absence")
with tempfile.TemporaryDirectory() as tmp:
    empty = Path(tmp) / "tables"
    empty.mkdir(parents=True)
    figures, notes = collect(empty)
    check("an empty run produces no figures", not figures, str(sorted(figures)))
    # The point of the note is that it tells a reader what to do. A run with no sweep must say
    # how to ask for one, not that the pipeline cannot do it.
    check("F5's note names the flag that requests it", "--dbr-sd-sweep" in notes.get("F5", ""),
          notes.get("F5", "")[:60])
    check("every unavailable figure carries a note",
          all(fid in notes for fid in UNAVAILABLE), str(sorted(set(UNAVAILABLE) - set(notes))))

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
sys.exit(1 if FAIL else 0)
