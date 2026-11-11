# Runs `scqc fixture-run` - the real CLI, then the real `scqc run` it calls - on a small synthetic
# object written twice under different names, in temporary directories; needs the analysis stack, not R.
"""The development suite's fixture reaches the pipeline, in both of its shapes (harness autodev W020,
B1-041).

The harness's development suite hands a tool ONE object of cells, written twice: the same numbers
under different column names, layer names and gene spellings. scQC's fixture tier ran `describe`,
because scQC reads a project, not an object - so no step ever ran on it. `scqc fixture-run` splits
the object into a project from the names it is TOLD and runs `scqc run --mode apply`. Checked here:

  1  each shape completes: a report (`reports/report.json` in --out, copied flat from the run), a
     deliverable in the project beside it (`<out>.project`), and the run's STATUS.json in --out,
     saying ok;
  2  the two shapes - the same numbers under other names - keep the same cells: a step that read a
     column, a layer or a gene spelling by name would make them differ;
  3  a name the object does not carry is refused by name, not guessed;
  4  every value it is TOLD reaches the step that reads it, by the run's own records (harness
     autodev W048, GAP-023): the assay, the design's factor, the mitochondrial prefix and the
     ribosomal pattern each print an AXIS line, and none is NOT REACHED;
  5  a value that did not arrive is named: the same run's records read against another assay,
     another factor, a prefix matching no gene, and a per-cell fraction that is zero in a library;
  6  the assay has no default: a fixture run that does not declare one is refused by the parser.

The object is this file's own synthetic one, not the harness's (which the suite builds where it
runs); the doublet detector is the stand-in, so nothing here says a call is right.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
fails: list[str] = []

try:
    import anndata as ad
    import numpy as np
    import pandas as pd
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")

#: shape -> (sample column, counts layer, mitochondrial prefix, ribosomal pattern)
SHAPES = {"a": ("sample", "counts", "MT-", "^RP[SL]"),
          "b": ("library_id", "raw_counts", "mt-", "^Rp[sl]")}


def spell(gene: str, shape: str) -> str:
    if shape == "a":
        return gene
    if gene.startswith("MT-"):
        return "mt-" + gene[3:].capitalize()
    return gene.capitalize() if gene.startswith(("RPL", "RPS")) else gene


def write_object(td: Path, shape: str) -> tuple:
    """Three libraries of cells, the same counts in both shapes."""
    rng = np.random.default_rng(7)
    genes = ["MT-CO1", "MT-ND1", "MT-ND2", "RPL3", "RPS6"] + [f"GENE{i:04d}" for i in range(295)]
    libs = {"L1": 80, "L2": 70, "L10": 60}             # one name a prefix of another
    rows, sample = [], []
    for s, n in libs.items():
        mu = rng.gamma(1.5, 1.0, size=len(genes)) + 0.5
        rows.append(rng.poisson(np.outer(rng.gamma(8.0, 1 / 8.0, size=n), mu * 4.0)))
        sample += [s] * n
    X = np.vstack(rows).astype("float32")
    key, layer, _mt, _ribo = SHAPES[shape]
    arm = {"L1": "ctrl", "L2": "treated", "L10": "treated"}
    obs = pd.DataFrame({key: pd.Categorical(sample), "arm_of": [arm[s] for s in sample]},
                       index=[f"{s}_C{i:05d}-1" for i, s in enumerate(sample)])
    A = ad.AnnData(X=X.copy(), obs=obs, var=pd.DataFrame(index=[spell(g, shape) for g in genes]))
    A.layers[layer] = X.copy()
    h5 = td / f"fixture_{shape}.h5ad"
    A.write_h5ad(h5)
    design = td / f"design_{shape}.csv"
    pd.DataFrame({key: list(libs) + ["L99"], "condition_of": [arm[s] for s in libs] + ["ctrl"]}
                 ).to_csv(design, index=False)
    return h5, design


def fixture_run(td: Path, shape: str, h5: Path, design: Path, **override):
    key, layer, mt, ribo = SHAPES[shape]
    args = {"--sample-key": key, "--counts-layer": layer, "--mt-prefix": mt, "--ribo-pattern": ribo,
            "--assay": "scrna", "--design-factors": "condition_of"}
    args.update(override)
    out = td / f"run_{shape}"
    cmd = [sys.executable, str(ROOT / "scqc_cli.py"), "fixture-run", "--observations", str(h5),
           "--design", str(design), "--dbr", "0.05", "--jobs", "4", "--out", str(out)]
    for k, v in args.items():
        if v is not None:
            cmd += [k, v]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    return p.returncode, p.stdout + p.stderr, out


if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        kept = {}
        run_a = None
        for shape in SHAPES:
            h5, design = write_object(td, shape)
            rc, log, out = fixture_run(td, shape, h5, design)
            report = out / "reports" / "report.json"
            status = out / "STATUS.json"
            if rc != 0 or not report.is_file():
                tail = "\n".join(ln for ln in log.splitlines() if ln.strip())[-900:]
                fails.append(f"1 shape {shape}: fixture-run did not complete (exit {rc}):\n{tail}")
                continue
            st = json.loads(status.read_text()) if status.is_file() else {}
            if st.get("status") != "ok":
                fails.append(f"1 shape {shape}: STATUS.json beside the project says "
                             f"{st.get('status')!r}: {st.get('headline')}")
            project = out.parent / (out.name + ".project")
            d = ad.read_h5ad(project / "results" / "latest" / "objects" / "cohort.deliverable.h5ad",
                             backed="r")
            kept[shape] = int(d.n_obs)
            d.file.close()
            print(f"  shape {shape}: completed, {kept[shape]} cells kept")
            axes = [ln for ln in log.splitlines() if ln.startswith("AXIS ")]
            want = ("AXIS assay 'scrna'", "AXIS design factor 'condition_of'",
                    f"AXIS mitochondrial prefix '{SHAPES[shape][2]}'",
                    f"AXIS ribosomal pattern '{SHAPES[shape][3]}'")
            if any(ln.startswith("AXIS NOT REACHED") for ln in axes) or \
                    not all(any(ln.startswith(w) for ln in axes) for w in want):
                fails.append(f"4 shape {shape}: the told values did not all reach their steps: {axes}")
            if shape == "a":
                run_a = (project / "results" / "latest").resolve()
        if len(kept) == 2 and kept["a"] != kept["b"]:
            fails.append(f"2: the two shapes kept different numbers of cells {kept} - a step read "
                         f"a name rather than the role it was told")
        print("1-2. both shapes reach a report and a deliverable, and keep the same cells")

        h5, design = write_object(td, "b")
        rc, log, _ = fixture_run(td, "b", h5, design, **{"--counts-layer": "counts"})
        if rc == 0 or "'counts' is absent" not in log:
            fails.append(f"3: a counts layer the object does not carry was not refused by name "
                         f"(exit {rc}): {log[-300:]}")
        print("3. a name the object does not carry is refused by name")
        print("4. every told value reaches the step that reads it, in both shapes")

        sys.path.insert(0, str(ROOT))
        try:
            from adapters.fixture_project import axes_reached
        except ImportError as e:
            fails.append(f"5: nothing reads the told values back from the run's records ({e})")
            run_a = None
        if run_a is None:
            fails.append("5: no completed shape-a run to read the told values back from")
        told = dict(assay="scrna", factors=["condition_of"], mt_prefix="MT-",
                    ribo_pattern="^RP[SL]", genes_mt=3, genes_ribo=2)
        if run_a is not None:
            lines, misses = axes_reached(run_a, **told)
            if misses or len(lines) != 4:
                fails.append(f"5: the run as told must show every value arriving: {lines} {misses}")
            zeroed = td / "zeroed"
            for sub in ("reports", "tables"):
                (zeroed / sub).mkdir(parents=True, exist_ok=True)
            (zeroed / "INPUTS.json").write_text((run_a / "INPUTS.json").read_text())
            (zeroed / "reports" / "report.json").write_text(
                (run_a / "reports" / "report.json").read_text())
            for p in (run_a / "tables").glob("*.percell.csv"):
                t = pd.read_csv(p)
                if p.name.startswith("L2."):
                    t["pct_counts_mt"] = 0.0
                t.to_csv(zeroed / "tables" / p.name, index=False)
            for case, run, change, says in (
                    ("another assay", run_a, {"assay": "snrna"}, "assay 'snrna' was told"),
                    ("another factor", run_a, {"factors": ["arm_of"]},
                     "design factor 'arm_of' was told"),
                    ("a prefix matching no gene", run_a, {"mt_prefix": "mt-", "genes_mt": 0},
                     "matches no gene of the object"),
                    ("a fraction that is zero in one library", zeroed, {}, "zero in every cell of L2")):
                _, m = axes_reached(run, **dict(told, **change))
                if not any(says in x for x in m):
                    fails.append(f"5 {case}: not named as a value that did not arrive: {m}")
        print("5. a told value that did not arrive is named, axis by axis")

        rc, log, _ = fixture_run(td, "a", *write_object(td, "a"), **{"--assay": None})
        if rc == 0 or "--assay" not in log:
            fails.append(f"6: a fixture run that declares no assay must be refused (exit {rc}): "
                         f"{log[-300:]}")
        print("6. the assay is declared, never defaulted")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the fixture reaches the pipeline in both shapes")
