# Builds graphs, converts tiny synthetic matrices and runs the CLI on a synthetic single-cell
# cohort in a temporary directory; runs no denoiser and needs no GPU.
"""Single cells are not nuclei: the route a whole-cell cohort takes through this pipeline.

Until 2026-10-02 the assay changed two numbers - the mitochondrial bounds and how k is chosen - and
nothing else. The first single-cell cohort (single-cell-harness ADR-0027) went down the nuclear
route whole: CellBender was scheduled for every library, its run failed on all four, and had it
succeeded the steps after it would have taken "has counts after denoising" as the cell call, bounded
the count floors by nuclear valleys, counted ribosomal markers against whole cells, and reached step
6 asking for a `cellbender_cell` column nothing in this tool writes. The PI's ruling: "do not run
cellbender for single cell" - "single cell and single nuclei is fundamentally different!".

Each check below fails on the code that preceded it:

  A  the graph: a single-cell library is not denoised, takes the aligner's call, and every step that
     selects cells reads that call; `aligner_cells` is required before anything is submitted; apply
     mode builds, with step 7's cell criterion named for the aligner (A8).
  B  the conversion writes the cell call as a column - the aligner's for cells, the denoiser's for
     nuclei (Q10: no conversion in this tool ever wrote `cellbender_cell`).
  C  step 5's measurement takes the call from that column, not from "has counts", and places the
     mitochondrial quartiles on called cells only.
  D  the CLI, on a synthetic single-cell cohort: no denoiser is invoked, the object keeps every
     droplet with the aligner's call marked, and step 2 says it had nothing to compare.
"""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
fails: list[str] = []

from engine import graph, steps  # noqa: E402
from engine.task import Refusal  # noqa: E402


class _P:
    results = Path("/tmp/r"); work = Path("/tmp/w"); project = Path("/tmp/p")
    decisions: dict = {}
    mode = "evidence"
    samples: list = []


def _row(s, assay, **extra):
    return {"sample": s, "platform": "10x", "species": "homo_sapiens",
            "reference": "homo_sapiens/GRCh38", "assay": assay, "matrix": f"/d/{s}/raw",
            "genotype": "WT" if s.endswith("1") else "Mut", **extra}


def _build(rows, mode="evidence"):
    _P.samples, _P.mode, _P.decisions = rows, mode, {}
    return {t.key: t for t in graph.main_stage(_P, "python", {}, {r["sample"]: {"mode": "accept"}
                                                                 for r in rows})}


# ---- A. the graph
CELLS = [_row("A1", "scrna", aligner_cells="/d/A1/filtered"),
         _row("B2", "scrna", aligner_cells="/d/B2/filtered")]
g = _build(CELLS)
gpu = sorted(k for k, t in g.items() if getattr(t, "gpu", False))
if gpu:
    fails.append(f"A1: a single-cell cohort scheduled GPU (denoiser) tasks: {gpu}")
for s in ("A1", "B2"):
    t = g.get(f"01_ambient/{s}")
    if t is None or t.fn is not getattr(steps, "_ambient_none", object()):
        fails.append(f"A2: 01_ambient/{s} is {getattr(t, 'fn', None)!r}, not the no-denoiser route")
    for k in (f"05_quality/{s}", f"04_doublets/{s}", f"06_cluster/{s}"):
        got = g[k].params.get("cell_call_key") if k in g else None
        if got != "aligner_cell":
            fails.append(f"A3: {k} reads the cell call from {got!r}; it must read the aligner's "
                         f"column, or it takes 'has counts' on a raw matrix where every droplet has")
    cls = g[f"06_cluster/{s}"].params.get("uninformative_classes")
    if cls != ["mt"]:
        fails.append(f"A4: criterion C on whole cells counts {cls!r}; ribosomal transcript is a "
                     f"cell's cytoplasm, not carry-over")
audit = g["01_ambient_audit"].params
if set(audit.get("not_denoised") or []) != {"A1", "B2"} or audit.get("per_sample"):
    fails.append(f"A5: the ambient audit would measure a fraction for libraries nothing denoised "
                 f"(per_sample {sorted(audit.get('per_sample') or {})}, not_denoised "
                 f"{audit.get('not_denoised')})")
if (g["02_cells"].params.get("call_by") or {}) != {"A1": "aligner", "B2": "aligner"}:
    fails.append(f"A6: step 2 is not told the call is the aligner's: {g['02_cells'].params.get('call_by')}")

try:
    _build([_row("A1", "scrna"), _row("B2", "scrna")])
    fails.append("A7: a single-cell cohort with no aligner_cells built; it has no cell call at all")
except Refusal as e:
    if "aligner_cells" not in str(e):
        fails.append(f"A7: refused, but without naming aligner_cells: {str(e)[:120]}")
# Apply mode builds for single cells, and step 7 reads the aligner's call: its cell criterion is
# named for the aligner, not the denoiser (ALIGNER_CELL_CRITERION in the adapter).
try:
    ga = _build(CELLS, mode="apply")
    for s in ("A1", "B2"):
        if ga[f"07_measure/{s}"].params.get("cell_call_key") != "aligner_cell":
            fails.append(f"A8: 07_measure/{s} would take 'has counts' as the cell call on a raw "
                         f"matrix, calling every droplet a cell")
except Refusal as e:
    fails.append(f"A8: apply mode refused for single cells: {str(e)[:120]}")
# Step 5 on single cells: no valley over raw droplets, the ceiling derived at the declared floor.
for s in ("A1", "B2"):
    q5 = g[f"05_quality/{s}"].params
    if not q5.get("valleys_not_used") or q5.get("mito_floor_umi") != 500:
        fails.append(f"A11: 05_quality/{s} would measure a valley over raw droplets, or derive the "
                     f"ceiling below the floor it is applied above: {q5.get('valleys_not_used')!r}, "
                     f"{q5.get('mito_floor_umi')!r}")
# The nuclear-fraction criterion points the other way in a whole cell: on scrna the column is
# IGNORED and the barrier is told which libraries carried it (harness autodev W020, B1-076) -
# the run is not refused over a column a shared samplesheet carries for every row.
try:
    gi = _build([_row("A1", "scrna", aligner_cells="/d/f", cellreads_stats="/d/CellReads.stats"),
                 _row("B2", "scrna", aligner_cells="/d/f")])
    if gi["05_quality/A1"].params.get("cellreads_stats"):
        fails.append("A12: the nucleus-calibrated nuclear-fraction criterion was armed on single cells")
    if (gi["05_quality"].params.get("cellreads_stats_ignored") or {}).get("A1") != "/d/CellReads.stats":
        fails.append(f"A12: the ignored column is not recorded for the report: "
                     f"{gi['05_quality'].params.get('cellreads_stats_ignored')}")
except Refusal as e:
    fails.append(f"A12: a single-cell cohort was refused over cellreads_stats: {str(e)[:120]}")

# Nuclei are untouched: denoised on the GPU, the call left to the denoiser, both marker classes.
NUC = [_row("A1", "snrna"), _row("B2", "snrna")]
gn = _build(NUC)
if not all(getattr(gn.get(f"01_ambient/{s}"), "gpu", False) for s in ("A1", "B2")):
    fails.append("A9: nuclei are no longer denoised")
if any("cell_call_key" in t.params for t in gn.values()):
    fails.append("A9: a nuclear graph names a cell-call column; its call is the denoiser's")
if any("valleys_not_used" in gn[f"05_quality/{s}"].params for s in ("A1", "B2")):
    fails.append("A11: nuclei lost their valley")
if gn["06_cluster/A1"].params.get("uninformative_classes") != ["mt", "ribo"]:
    fails.append(f"A9: nuclei lost a marker class: {gn['06_cluster/A1'].params.get('uninformative_classes')}")
# A single-cell library that arrives DENOISED ELSEWHERE is still a whole cell.
SUP = {"ambient_h5": "/cb/x.h5", "ambient_tool": "CellBender", "ambient_version": "0.3.2",
       "ambient_params": "--fpr 0", "ambient_produced_by": "a collaborator"}
gs = _build([_row("A1", "scrna", **SUP), _row("B2", "scrna", **SUP)])
if gs["06_cluster/A1"].params.get("uninformative_classes") != ["mt"]:
    fails.append("A10: a supplied single-cell library was profiled with the nuclear marker classes")
# Step 5 is told whose call each library's cells are: the denoiser's for a library supplied
# already denoised, the aligner's for one this pipeline did not denoise (C04).
if (gs["05_quality"].params.get("call_by") or {}) != {"A1": "denoiser", "B2": "denoiser"} \
        or (g["05_quality"].params.get("call_by") or {}) != {"A1": "aligner", "B2": "aligner"}:
    fails.append(f"A10: step 5 is not told whose called cells its floors judge: supplied "
                 f"{gs['05_quality'].params.get('call_by')}, raw {g['05_quality'].params.get('call_by')}")
# A13. The assay as a person types it: padded, or in capitals. Ingest strips it before checking,
# so it passes there; the graph lowered it without stripping, and ` scrna` then stopped the build
# with a KeyError at step 6's marker classes - before anything was submitted, but after ingest
# had passed it (C03).
try:
    gp = _build([_row("A1", " scrna", aligner_cells="/d/A1/filtered"),
                 _row("B2", "SCRNA ", aligner_cells="/d/B2/filtered")])
    got = {gp["05_quality"].params["assay"].get(s) for s in ("A1", "B2")}
    if got != {"scrna"} or gp["06_cluster/A1"].params.get("uninformative_classes") != ["mt"] \
            or gp["01_ambient/A1"].fn is not steps._ambient_none:
        fails.append(f"A13: a padded or capitalised 'scrna' reached the steps as {got}, or took "
                     f"the nuclear route")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"A13: a padded or capitalised 'scrna' stopped the build: "
                 f"{type(e).__name__}: {e}")
_P.mode = "evidence"
print("A. the graph: single cells take the aligner's call, nuclei the denoiser's")

# ---- G. step 5 says whose called cells its declared floors were applied to (C04)
import types as _types  # noqa: E402

from engine.pipeline import Pipeline as _Pipeline  # noqa: E402


class _Done:
    ok = True

    def __init__(self, metrics):
        self.metrics = metrics


def _floors_table(call_by):
    """Run step 5's barrier over two measured single-cell libraries; (applied_to, metrics)."""
    with tempfile.TemporaryDirectory() as td:
        pp = _types.SimpleNamespace(results=Path(td) / "results", results_by_key={
            f"05_quality/{s}": _Done({"valleys": {}, "bimodal": {},
                                      "mito_quartiles": {"n": 160, "median": 2.0, "q1": 1.5,
                                                         "q3": 2.6, "mad": 0.4},
                                      "mito_population": {"floor_umi": 500}})
            for s in ("A1", "B2")})
        (pp.results / "tables").mkdir(parents=True)
        t = _types.SimpleNamespace(key="05_quality", params={
            "samples": ["A1", "B2"], "assay": {"A1": "scrna", "B2": "scrna"},
            "light_floor": 200, "call_by": call_by})
        out = steps._quality_stage(t, pp, None)
        with open(pp.results / "tables" / "count_floors.csv", encoding="utf-8") as fh:
            applied = {r["applied_to"] for r in csv.DictReader(fh)}
        return applied, out["metrics"]


try:
    sup, sup_m = _floors_table({"A1": "denoiser", "B2": "denoiser"})
    raw_, raw_m = _floors_table({"A1": "aligner", "B2": "aligner"})
    mix, _ = _floors_table({"A1": "aligner", "B2": "denoiser"})
    if sup != {"the denoiser's called cells"}:
        fails.append(f"G1: floors on libraries supplied already denoised say they judged {sup}")
    if raw_ != {"the aligner's called cells"}:
        fails.append(f"G2: floors on libraries no denoiser touched say they judged {raw_}")
    if not any("aligner" in a and "denoiser" in a and "A1" in a and "B2" in a for a in mix):
        fails.append(f"G3: a cohort of both routes is described as {mix}")
    basis = {r["name"]: r.get("basis", "") for r in _Pipeline._parameters_block(
        _types.SimpleNamespace(FIXED_PARAMETERS=_Pipeline.FIXED_PARAMETERS,
                               DECLARED_PARAMETERS=_Pipeline.DECLARED_PARAMETERS, decisions={},
                               tools={}, samples=[],
                               results_by_key={"05_quality": _Done(sup_m),
                                               "07_apply": _Done({})}))}
    if "the denoiser's called cells" not in basis.get("UMI floor", "") \
            or "aligner" in basis.get("UMI floor", ""):
        fails.append(f"G4: the parameter table's basis for the floor: {basis.get('UMI floor')!r}")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"G: step 5's floors could not be measured here: {type(e).__name__}: {e}")
print("G. step 5 names the caller whose cells its declared floors judged, by each library's route")

# ---- H. the report names cells as cells (GAP-024)
#
# Written so that it also runs against the code before the words were per assay - every caption
# fixed text, every question one string, the step texts one table - and fails THERE on what it
# asserts rather than on a name that did not exist yet.
from report import build as _rb  # noqa: E402
from report import collect as _rc  # noqa: E402

try:
    from report.units import fill as _fill  # noqa: E402
except ImportError:
    def _fill(text, assay):
        return text
_question = getattr(_rb, "_question", None) or (
    lambda fid, payload: _rb.FIGURE_QUESTIONS.get(fid, ""))


def _step_text(step, assay):
    try:
        return steps.step_text(step, assay)
    except TypeError:
        return steps.step_text(step)


_NUC = ("nucleus", "nuclei")
for fid, cap in sorted(_rc.CAPTIONS.items()):
    if any(w in _fill(cap, "scrna").lower() for w in _NUC):
        fails.append(f"H1: {fid}'s caption on a single-cell cohort still says nucleus: "
                     f"{_fill(cap, 'scrna')[:120]}")
if "per nucleus" not in _fill(_rc.CAPTIONS["F7"], "snrna"):
    fails.append("H1: F7's caption no longer names nuclei on a nucleus cohort")
q11 = _question("F11", {"run": {"assay": "scrna"}})
if "cells" not in q11 or "nucle" in q11:
    fails.append(f"H2: F11 asks of a single-cell cohort: {q11!r}")
if "nuclei" not in _question("F11", {"run": {"assay": "snrna"}}):
    fails.append("H2: F11 no longer asks about nuclei on a nucleus cohort")
for step in steps.STEP_TEXT:
    what, cannot = _step_text(step, "scrna")
    if any(w in (what + cannot).lower() for w in _NUC):
        fails.append(f"H3: step {step}'s text on a single-cell cohort says nucleus: {what[:100]}")
if "denoise" not in _step_text("01_ambient", "scrna")[0] \
        or _step_text("01_ambient", "scrna") == _step_text("01_ambient", "snrna"):
    fails.append("H3: step 1's text is the nuclear one on single cells")
print("H. captions, F11's question and the step texts name cells on a single-cell cohort")

# ---- E. step 5's floors and step 6's thresholds, as cells and as nuclei
from engine.pipeline import step_module  # noqa: E402
qm = step_module("quality")
try:
    pu, pg = qm.declared_floor("umi", "scrna", light_floor=200), qm.declared_floor("genes", "scrna")
    if (pu.constant, pg.constant, pu.provenance) != (500, 250, "declared"):
        fails.append(f"E1: declared single-cell floors {pu.constant}/{pg.constant} {pu.provenance}")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"E1: no declared floor for single cells: {type(e).__name__}: {e}")
try:
    qm.declared_floor("umi", "snrna")
    fails.append("E2: nuclei were handed a declared floor; theirs is measured")
except Exception:                                                     # noqa: BLE001
    pass
cf = step_module("cluster_flags")
flat = [{"umi_frac_of_sample": 1.0, "median_pct_mt": m, "pct_uninformative": 0.0,
         "pct_mt_markers": 0.0, "pct_ribo_markers": None} for m in (2.0, 2.5, 3.0, 3.5, 4.0, 5.0)]
try:
    thr = cf.propose(flat, b_floor=cf.B_FLOOR_PCT["scrna"])
    if thr.c_uninformative is not None or thr.b_pct_mt < 10.0:
        fails.append(f"E3: thresholds on healthy cells: {thr}")
    fl = cf.apply_flags(flat, thr)
    if any(r.get("FLAG") or r.get("B") for r in fl.rows):
        fails.append("E3: a healthy whole-cell cluster (median mito 2-5%) was flagged")
except Exception as e:                                                # noqa: BLE001
    fails.append(f"E3: a cohort with no uninformative marker stopped step 6: {type(e).__name__}: {e}")
print("E. declared floors on cells only; step 6 stands on healthy cells")

# ---- B, C need the analysis stack.
try:
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    HAVE = True
except ImportError:
    HAVE = False
    print("SKIP B-D: anndata / numpy / scipy not importable here")

GENES = (["MT-CO1", "MT-ND1", "MT-ND2", "RPL3", "RPS6", "RPL13"]
         + [f"GENE{i}" for i in range(54)])


def _cohort_counts(seed):
    """Debris, called cells, and droplets deep enough to clear the light floor the caller rejected."""
    rng = np.random.default_rng(seed)
    n_g = len(GENES)
    rows, bcs, called = [], [], []

    def droplet(total, mito_frac):
        p = np.full(n_g, (1 - mito_frac) / (n_g - 3))
        p[:3] = mito_frac / 3
        return rng.multinomial(int(total), p)

    for i in range(400):                                   # debris: 1-15 UMI
        rows.append(droplet(rng.integers(1, 16), 0.05)); bcs.append(f"D{i:05d}-1")
    for i in range(160):                                   # cells the aligner called
        rows.append(droplet(rng.integers(1500, 6000), rng.uniform(0.02, 0.12)))
        bcs.append(f"C{i:05d}-1"); called.append(bcs[-1])
    for i in range(40):                                    # 400-900 UMI, NOT called, mostly mito
        rows.append(droplet(rng.integers(400, 900), 0.45)); bcs.append(f"U{i:05d}-1")
    return sp.csr_matrix(np.vstack(rows)), bcs, called


def _write_10x(dirpath: Path, X, bcs):
    dirpath.mkdir(parents=True, exist_ok=True)
    G = sp.coo_matrix(X.T)                                 # features x barcodes
    with gzip.open(dirpath / "matrix.mtx.gz", "wt") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{G.shape[0]} {G.shape[1]} {G.nnz}\n")
        for r, c, v in zip(G.row, G.col, G.data):
            fh.write(f"{r + 1} {c + 1} {int(v)}\n")
    with gzip.open(dirpath / "barcodes.tsv.gz", "wt") as fh:
        fh.write("".join(f"{b}\n" for b in bcs))
    with gzip.open(dirpath / "features.tsv.gz", "wt") as fh:
        fh.write("".join(f"ENSG{i:011d}\t{g}\tGene Expression\n" for i, g in enumerate(GENES)))


if HAVE:
    from adapters import scanpy_ops as so

if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        X, bcs, called = _cohort_counts(0)
        _write_10x(td / "raw", X, bcs)
        _write_10x(td / "filtered", X[[bcs.index(b) for b in called]], called)

        # B1. the aligner's call, as a column, every droplet kept
        dest = td / "S_ambient.h5"
        try:
            _, met = so._op_convert(None, {"source": str(td / "raw"), "dest": str(dest),
                                           "cell_key": "aligner_cell",
                                           "cell_barcodes": str(td / "filtered")}, str(td / "c"))
            obj = ad.read_h5ad(dest)
            if obj.n_obs != len(bcs):
                fails.append(f"B1: the object kept {obj.n_obs} of {len(bcs)} droplets; the valleys "
                             f"need the debris")
            got = set(obj.obs_names[obj.obs["aligner_cell"].to_numpy().astype(bool)])
            if got != set(called):
                fails.append(f"B1: aligner_cell marks {len(got)} droplets, the aligner called "
                             f"{len(called)}")
            if met.get("cell_call") != "aligner" or met.get("cells_called") != len(called):
                fails.append(f"B1: metrics do not say whose call it is: {met}")
        except Exception as e:                                        # noqa: BLE001
            fails.append(f"B1: {type(e).__name__}: {e}")

        # B2. a call naming droplets the matrix does not hold is not this library's call
        _write_10x(td / "foreign", X[:3], ["ZZZ0001-1", "ZZZ0002-1", called[0]])
        try:
            so._op_convert(None, {"source": str(td / "raw"), "dest": str(td / "x.h5"),
                                  "cell_key": "aligner_cell", "cell_barcodes": str(td / "foreign")},
                           str(td / "c2"))
            fails.append("B2: a call from another run was accepted")
        except Exception as e:                                        # noqa: BLE001
            if "not droplets of" not in str(e):
                fails.append(f"B2: refused for the wrong reason: {e}")

        # B3. Q10 - the denoiser's route writes the column step 6 reads
        try:
            Xd = X.tolil()
            for i in range(400, len(bcs)):
                if bcs[i].startswith("U"):
                    Xd[i, :] = 0                                  # what a denoiser calls empty
            _write_10x(td / "denoised", Xd.tocsr(), bcs)
            _, met = so._op_convert(None, {"source": str(td / "denoised"),
                                           "dest": str(td / "N_ambient.h5"),
                                           "cell_key": "cellbender_cell"}, str(td / "c3"))
            objn = ad.read_h5ad(td / "N_ambient.h5")
            if "cellbender_cell" not in objn.obs.columns:
                fails.append("B3: the denoiser's route wrote no cellbender_cell; step 6 refuses")
            elif int(objn.obs["cellbender_cell"].sum()) != int((np.asarray(Xd.sum(axis=1)).ravel() > 0).sum()):
                fails.append("B3: cellbender_cell is not 'has denoised counts'")
        except Exception as e:                                        # noqa: BLE001
            fails.append(f"B3: {type(e).__name__}: {e}")
        print("B. the conversion writes the cell call as a column, for either caller")

        # C. step 5's measurement takes the call from the column
        try:
            obj = ad.read_h5ad(dest)
            outs, m = so._op_valley(obj, {"sample": "S", "metrics": ["umi", "genes"],
                                          "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]",
                                          "mito_floor_umi": 200,
                                          "cell_call_key": "aligner_cell"}, str(td / "v"))
            csvp = next(o for o in outs if str(o).endswith(".called_barcodes.csv"))
            with open(csvp, encoding="utf-8") as fh:
                got = {r["barcode"] for r in csv.DictReader(fh)}
            if got != set(called):
                fails.append(f"C1: step 5 recorded {len(got)} called barcodes; the aligner called "
                             f"{len(called)}. 'Has counts' on a raw matrix calls every droplet")
            if m.get("n_called_by_denoiser") is not None or m.get("n_called") != len(called):
                fails.append(f"C2: the count is filed under the denoiser: {m.get('n_called_by_denoiser')}"
                             f" / n_called {m.get('n_called')}")
            pop = m.get("mito_population") or {}
            # The 40 rejected droplets clear the 200-UMI floor and are ~45% mitochondrial; on cells
            # only, the population above the floor is the 160 called cells.
            if pop.get("n_above_floor") != len(called) or not pop.get("called_only"):
                fails.append(f"C3: the mitochondrial quartiles were placed on {pop.get('n_above_floor')} "
                             f"droplets, not on the {len(called)} called cells")
        except Exception as e:                                        # noqa: BLE001
            fails.append(f"C: {type(e).__name__}: {e}")
        print("C. step 5 measures the aligner's cells, not every droplet with counts")

        # C4. no valley on cells, but only with the reason in hand
        try:
            _, m4 = so._op_valley(ad.read_h5ad(dest), {
                "sample": "S", "metrics": [], "valleys_not_used": "declared floors",
                "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]", "mito_floor_umi": 500,
                "cell_call_key": "aligner_cell"}, str(td / "v4"))
            if m4.get("valleys"):
                fails.append(f"C4: valleys were measured anyway: {m4.get('valleys')}")
            if (m4.get("mito_population") or {}).get("floor_umi") != 500:
                fails.append("C4: the ceiling's population is not taken at the declared floor")
        except Exception as e:                                        # noqa: BLE001
            fails.append(f"C4: {type(e).__name__}: {e}")
        try:
            so._op_valley(ad.read_h5ad(dest), {"sample": "S", "metrics": [], "mt_prefix": "MT-",
                                               "ribo_pattern": "^RP[SL]", "mito_floor_umi": 200},
                          str(td / "v5"))
            fails.append("C5: an empty metric list with no reason was accepted")
        except Exception:                                             # noqa: BLE001
            pass

        # F. step 7's measurement: the aligner's criterion, over called cells only
        try:
            outs, m7 = so._op_apply_measure(ad.read_h5ad(dest), {
                "sample": "S", "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]",
                "umi_floor": 500, "gene_floor": 250, "mito_ceiling_pct": 20.0, "light_floor": 200,
                "doublet_csv": None, "nf_csv": None, "nf_floor": None, "nf_trigger_pct": None,
                "cell_call_key": "aligner_cell"}, str(td / "S"))
            with open(outs[0], encoding="utf-8") as fh:
                rd = csv.DictReader(fh)
                hdr, body = rd.fieldnames, list(rd)
            if "fail_not_aligner_cell" not in hdr or "fail_not_cellbender_cell" in hdr \
                    or "aligner_cell" not in hdr:
                fails.append(f"F1: step 7 names the cell criterion {[h for h in hdr if 'cell' in h]}")
            if len(body) != len(called) or m7.get("n_not_tabled") != len(bcs) - len(called):
                fails.append(f"F2: the per-cell table holds {len(body)} rows for {len(called)} "
                             f"called cells (left out: {m7.get('n_not_tabled')})")
        except Exception as e:                                        # noqa: BLE001
            fails.append(f"F: {type(e).__name__}: {e}")
        print("F. step 7 measures called cells under the aligner's criterion")

# ---- D. the CLI, on a synthetic single-cell cohort
if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sheet = td / "samplesheet.tsv"
        cols = ["sample", "platform", "species", "reference", "assay", "matrix", "aligner_cells",
                "genotype", "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]
        rows = []
        for i, s in enumerate(("WT_1", "MUT_1")):
            X, bcs, called = _cohort_counts(10 + i)
            _write_10x(td / s / "raw", X, bcs)
            _write_10x(td / s / "filtered", X[[bcs.index(b) for b in called]], called)
            # The assay as typed, padded on one row and in capitals on the other (C03): ingest
            # strips it, and the graph must read it the same way or the build stops at step 6.
            rows.append([s, "10x", "homo_sapiens", "homo_sapiens/GRCh38",
                         " scrna" if i == 0 else "SCRNA",
                         str(td / s / "raw"), str(td / s / "filtered"),
                         "WT" if s.startswith("WT") else "Mut", "MT-", "^RP[SL]", "0.05", "0.05"])
        sheet.write_text("\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
        reg = td / "registry.tsv"
        reg.write_text("species\tbuild\tpath\nhomo_sapiens\tGRCh38\t-\n")
        # A denoiser that cannot run, and that leaves a mark if anything calls it. Its output
        # goes to a task log, not to the CLI's, so the mark is a file rather than a message.
        mark = td / "DENOISER_INVOKED"
        fake = td / "no_denoiser.sh"
        fake.write_text(f"#!/bin/sh\ntouch '{mark}'\nexit 3\n")
        fake.chmod(0o755)
        rfake = td / "no_rscript.sh"
        rfake.write_text("#!/bin/sh\nexit 3\n")
        rfake.chmod(0o755)
        p = subprocess.run(
            [sys.executable, str(ROOT / "scqc_cli.py"), "run", "--project", str(td / "proj"),
             "--samplesheet", str(sheet), "--registry", str(reg), "--mode", "evidence",
             "--executor", "local", "--python", sys.executable, "--cellbender", str(fake),
             "--rscript", str(rfake), "--jobs", "2", "--seed", "0"],
            capture_output=True, text=True, timeout=900)
        log = p.stdout + p.stderr
        if mark.exists():
            fails.append("D1: the CLI invoked a denoiser on a single-cell cohort")
        res = [d for d in (td / "proj" / "results").iterdir()
               if d.is_dir() and not d.is_symlink() and d.name != "latest"]
        if not res:
            fails.append(f"D: no results directory; the run stopped before step 1:\n{log[-1500:]}")
        else:
            R = res[0]
            for s in ("WT_1", "MUT_1"):
                h5 = R / "objects" / f"{s}_ambient.h5"
                if not h5.exists():
                    fails.append(f"D2: {h5.name} was not written:\n{log[-1500:]}")
                    continue
                o = ad.read_h5ad(h5, backed="r")
                if "aligner_cell" not in o.obs.columns or o.n_obs != 600:
                    fails.append(f"D2: {s}: {o.n_obs} droplets, columns {list(o.obs.columns)}")
                o.file.close()
            cc = R / "tables" / "cell_calls.csv"
            if not cc.exists():
                fails.append(f"D3: step 2 wrote no cell_calls.csv:\n{log[-1500:]}")
            else:
                with open(cc, encoding="utf-8") as fh:
                    got = {r["sample"]: r for r in csv.DictReader(fh)}
                for s in ("WT_1", "MUT_1"):
                    r = got.get(s) or {}
                    if r.get("aligner") != "160" or r.get("denoiser") != "" or r.get("lost") != "":
                        fails.append(f"D3: {s} cell_calls row {r}; the aligner's 160 with the "
                                     f"denoiser and the loss BLANK, since neither exists")
            # D4. The report names what it observed: the cohort's assay, and step texts written
            # for cells (GAP-024). Read from the payload the report was built from.
            import json as _json
            pl = R / "reports" / "payload.json"
            if not pl.exists():
                fails.append(f"D4: no payload.json beside the report:\n{log[-800:]}")
            else:
                payload = _json.loads(pl.read_text(encoding="utf-8"))
                if (payload.get("run") or {}).get("assay") != "scrna":
                    fails.append(f"D4: the report's run block names the assay "
                                 f"{(payload.get('run') or {}).get('assay')!r}")
                for st in payload.get("steps") or []:
                    text = f"{st.get('what_it_does', '')} {st.get('cannot_establish', '')}"
                    if any(w in text.lower() for w in ("nucleus", "nuclei")):
                        fails.append(f"D4: the report's {st.get('key')} text says nucleus on "
                                     f"single cells: {text[:100]}")
        print("D. the CLI: no denoiser invoked, every droplet kept, step 2 compares nothing")

# ---- I. every nucleus-calibrated line is declared per assay, and the policy names CellBender
# (harness autodev W020, B1-074, B1-075, B1-105)
import importlib.util  # noqa: E402


def _load(rel, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


_am = _load("modules/01_ambient/ambient.py", "w020_ambient")
_cg = _load("modules/02_cells/cellcall_gate.py", "w020_cg")
_dh = _load("modules/04_doublets/doublet_health.py", "w020_dh")
_aa = _load("modules/01_ambient/audit_ambient.py", "w020_aa")
_qm = _load("modules/05_quality/quality.py", "w020_qm")
_cf = _load("modules/06_cluster_check/cluster_flags.py", "w020_cf")
ASSAYS = set(_am.ASSAYS)
TABLES = {"cluster_flags.RESOLUTION": (_cf, "RESOLUTION"),
          "doublet_health.ZERO_RATE": (_dh, "ZERO_RATE"),
          "doublet_health.DESIGN_REFUSE": (_dh, "DESIGN_REFUSE"),
          "doublet_health.MATERIAL": (_dh, "MATERIAL"),
          "cellcall_gate.REVIEW_LOST": (_cg, "REVIEW_LOST"),
          "cellcall_gate.REFUSE_LOST": (_cg, "REFUSE_LOST"),
          "cellcall_gate.REFUSE_RATIO": (_cg, "REFUSE_RATIO"),
          "cellcall_gate.MATERIAL_LOSS": (_cg, "MATERIAL_LOSS"),
          "cellcall_gate.DESIGN_REFUSE": (_cg, "DESIGN_REFUSE"),
          "audit_ambient.DESIGN_REFUSE": (_aa, "DESIGN_REFUSE"),
          "audit_ambient.MATERIAL_REMOVAL": (_aa, "MATERIAL_REMOVAL"),
          "quality.DESIGN_REFUSE": (_qm, "DESIGN_REFUSE"),
          "quality.MATERIAL": (_qm, "MATERIAL"),
          "ambient.CELLBENDER_RUNS": (_am, "CELLBENDER_RUNS")}
for label, (mod, name) in TABLES.items():
    table = getattr(mod, name, None)
    if not isinstance(table, dict) or set(table) != ASSAYS:
        fails.append(f"I1: {label} is not a table keyed by exactly {sorted(ASSAYS)}: {table!r}")
for label, mod in (("cluster_flags", _cf), ("doublet_health", _dh), ("cellcall_gate", _cg),
                   ("audit_ambient", _aa), ("quality", _qm)):
    src = getattr(mod, "RESOLUTION_SOURCE" if label == "cluster_flags" else "LINES_SOURCE", None)
    if not isinstance(src, dict) or set(src) != ASSAYS or "nuclei" not in str(src.get("scrna")):
        fails.append(f"I1: {label} carries no per-assay provenance saying the single-cell line was "
                     f"calibrated on nuclei: {src!r}")
# No default assay: an unknown one is refused by every gate, never read as either.
for label, call in (
        ("cellcall_gate.gate", lambda: _cg.gate({"A": {"aligner": 10, "cellbender": 12, "lost": 0}},
                                                assay="spatial")),
        ("doublet_health.health", lambda: _dh.health({"A": 0.05}, assay="spatial")),
        ("audit_ambient.audit", lambda: _aa.audit([{"sample": "A", "fraction_removed_overall": 0.1,
                                                    "genes_fully_removed": 0}], None, None,
                                                  assay="spatial")),
        ("quality.assess_mito_removal", lambda: _qm.assess_mito_removal(
            {"A": [1.0, 2.0, 3.0, 4.0]}, {"A": 5.0}, assay="spatial"))):
    try:
        call()
        fails.append(f"I2: {label} accepted assay 'spatial' instead of refusing it")
    except (ValueError, _qm.ThresholdRefusal):
        pass
    except TypeError as e:
        fails.append(f"I2: {label} takes no assay: {e}")
# The policy names CellBender, and nothing wider.
r = _am.plan_ambient("x", "scrna").reason
if "CellBender" not in r or "never run by this pipeline" in r or "not corrected at all" in r:
    fails.append(f"I3: plan_ambient's single-cell reason does not name CellBender, or bans more: {r[:160]}")
# A nucleus library's correction can be skipped by declaration - with a reason AND an approver.
sk = _am.plan_ambient("n", "snrna", skip=True, skip_reason="no GPU queue at this site; pilot ambient 4%",
                      approved_by="the PI")
if sk.state != "SKIP" or not sk.declared_skip or sk.declared_skip.get("approved_by") != "the PI":
    fails.append(f"I4: a declared skip on snrna is {sk.state}, declared_skip {sk.declared_skip}")
for kw in ({"skip_reason": "no GPU queue"}, {"approved_by": "the PI"}, {}):
    try:
        _am.plan_ambient("n", "snrna", skip=True, **kw)
        fails.append(f"I4: a snrna skip with {kw} was allowed - half a declaration skips nothing")
    except _am.AmbientRefusal:
        pass
# ...and through the graph: the skipped library takes the not-denoised route, the step-1 task and
# the audit carry the declaration, and the resolution comes from the assay's table.
gsk = _build([_row("A1", "snrna", aligner_cells="/d/A1/filtered",
                   ambient_skip_reason="no GPU queue at this site; pilot ambient 4%",
                   ambient_skip_approved_by="the PI"),
              _row("B2", "snrna")])
t1 = gsk.get("01_ambient/A1")
if t1 is None or t1.fn is not steps._ambient_none or t1.params.get("approved_by") != "the PI":
    fails.append(f"I5: a declared skip did not send the library down the not-denoised route with "
                 f"its declaration: {getattr(t1, 'params', None)}")
if getattr(gsk.get("01_ambient/B2"), "gpu", False) is not True:
    fails.append("I5: the other nucleus library is no longer denoised")
if (gsk["01_ambient_audit"].params.get("declared_skips") or {}).get("A1", {}).get("approved_by") != "the PI":
    fails.append(f"I5: the audit does not carry the declared skip: "
                 f"{gsk['01_ambient_audit'].params.get('declared_skips')}")
try:
    _build([_row("A1", "snrna", ambient_skip_reason="no GPU queue"), _row("B2", "snrna")])
    fails.append("I5: a skip with a reason and no approver built")
except Refusal as e:
    if "approver" not in str(e):
        fails.append(f"I5: refused without naming the missing approver: {str(e)[:120]}")
for s, want in (("A1", _cf.RESOLUTION["snrna"]),):
    got = gsk[f"06_cluster/{s}"].params.get("resolution")
    if got != want:
        fails.append(f"I6: with no --resolution, 06_cluster/{s} takes {got!r}, not the assay's {want}")
if g["06_cluster/A1"].params.get("resolution") != _cf.RESOLUTION["scrna"]:
    fails.append(f"I6: single cells take {g['06_cluster/A1'].params.get('resolution')!r} rather than "
                 f"their assay's entry")
if {t.params.get("assay", {}).get("A1") for k, t in gsk.items()
        if k in ("02_cells", "04_doublet_health", "01_ambient_audit")} != {"snrna"}:
    fails.append("I7: the gating tasks are not told the assay whose lines they judge by")
print("I. every nucleus-calibrated line is a per-assay table with its provenance; the policy names")
print("   CellBender; a nucleus library's correction is skipped only by a declared, approved reason")

print("\n" + "=" * 74)
if fails:
    print("FAILED:")
    for f in fails:
        print(" -", f)
    raise SystemExit(1)
print("proved: a single-cell cohort is not denoised, takes the aligner's call into every step that")
print("selects cells, is bounded and profiled as cells, and cannot reach a removal under a denoiser's")
print("name; nuclei take the route they took before.")
