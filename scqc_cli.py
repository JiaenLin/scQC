"""scQC command-line interface.

WHAT THIS CAN AND CANNOT DO — read this before using it.

scQC does not read matrices, run aligners, denoise, cluster, or write filtered objects. Every
module in this repository consumes numbers that YOU computed, and returns a decision about them.
That is the design: the judgement is the part worth reviewing, and it is separable from the
computation.

So this CLI audits a quality-control run you performed elsewhere. Point it at your own tables and
it will tell you whether your cell calls lost a population, whether an ambient correction fell
unevenly across your design, whether a doublet rate is a measurement or a prior, whether a
threshold sits somewhere defensible, and whether a removal you are about to make has a recorded
approval. It will refuse where refusing is warranted.

It cannot produce those tables for you. Commands that would require it to are absent rather than
stubbed, because a stub that prints "not implemented" still appears in `--help` as a capability.

Import mechanics: the step directories are named `00_ingest`, `01_ambient` and so on, which are
not valid Python identifiers, and there are no `__init__.py` files. Modules are therefore loaded
by file path through importlib rather than imported normally.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VERSION = (ROOT / "VERSION").read_text(encoding="utf-8").strip() if (ROOT / "VERSION").exists() else "unknown"

# Findings contain typographic characters. A Windows console defaults to cp1252 and renders them
# as mojibake, which makes a correct refusal look like a corrupted one - so force UTF-8 on the
# streams we own. Guarded because `reconfigure` needs Python 3.7+ and a real text stream, and
# stdout may be a pipe or a captured buffer.
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):                              # pragma: no cover
        pass

_MODULES = {
    "verify_raw": ROOT / "lib" / "verify_raw.py",
    "ingest": ROOT / "modules" / "00_ingest" / "ingest.py",
    "ambient": ROOT / "modules" / "01_ambient" / "ambient.py",
    "lr_policy": ROOT / "modules" / "01_ambient" / "lr_policy.py",
    "audit_ambient": ROOT / "modules" / "01_ambient" / "audit_ambient.py",
    "cellcall_gate": ROOT / "modules" / "02_cells" / "cellcall_gate.py",
    "light_floor": ROOT / "modules" / "03_light_floor" / "light_floor.py",
    "doublet": ROOT / "modules" / "04_doublets" / "doublet.py",
    "doublet_health": ROOT / "modules" / "04_doublets" / "doublet_health.py",
    "quality": ROOT / "modules" / "05_quality" / "quality.py",
    "cluster_flags": ROOT / "modules" / "06_cluster_check" / "cluster_flags.py",
    "apply": ROOT / "modules" / "07_apply" / "apply.py",
}

_cache: dict = {}


def load(name: str):
    """Load a step module by path. Cached, so repeated calls do not re-execute it."""
    if name in _cache:
        return _cache[name]
    path = _MODULES.get(name)
    if path is None or not path.exists():
        raise SystemExit(f"scqc: internal error - module '{name}' not found at {path}")
    # `lib/` must be importable before ingest, which imports verify_raw at module level.
    lib = str(ROOT / "lib")
    if lib not in sys.path:
        sys.path.insert(0, lib)
    spec = importlib.util.spec_from_file_location(f"scqc_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    _cache[name] = mod
    return mod


# --------------------------------------------------------------------------------------------
# helpers


def read_csv(path: Path) -> list[dict]:
    """Read a CSV, skipping comment lines. Returns a list of dicts."""
    if not path.exists():
        raise SystemExit(f"scqc: no such file: {path}")
    with path.open(encoding="utf-8", newline="") as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    if not lines:
        raise SystemExit(f"scqc: {path} contains no data")
    rows = list(csv.DictReader(lines))
    if not rows:
        raise SystemExit(f"scqc: {path} has a header but no data rows")
    return rows


def design_from(rows: list[dict], skip: set[str], max_levels: int = 6) -> dict:
    """Every column with 2..max_levels levels is a design factor, unless it is an identifier.

    Discovered, not declared. Three things are excluded: a constant column, which separates
    nothing; a column with more levels than the pipeline can test; and a column with ONE SAMPLE
    PER LEVEL. The last is the one that slips through a bare `<= max_levels` test - a replicate
    id in a four-sample cohort has four levels, and every differential check then computes a
    ratio between single libraries, which is arithmetic rather than evidence and is reported in
    the same words as a real design differential. A factor must leave at least one level holding
    more than one sample.
    """
    if not rows:
        return {}
    out = {}
    ceiling = min(max_levels, max(len(rows) - 1, 1))
    for col in rows[0]:
        if col in skip:
            continue
        vals = {r.get(col, "") for r in rows if (r.get(col) or "").strip()}
        if 2 <= len(vals) <= ceiling:
            out[col] = {r["sample"]: r[col] for r in rows if r.get(col)}
    return out


def note(msg: str, as_json: bool) -> None:
    """Print a human preamble. Under --json it goes to stderr.

    A header on stdout makes the document unparseable, so `scqc ... --json | jq` failed at
    character 0. Anything that is not the document itself belongs on stderr.
    """
    print(msg, file=sys.stderr if as_json else sys.stdout)


def emit(findings, as_json: bool) -> str:
    """Print findings and return the overall verdict."""
    sev_rank = {"REFUSE": 2, "REVIEW": 1, "ok": 0}
    worst = 0
    payload = []
    for f in findings:
        sev = getattr(f, "severity", "ok")
        worst = max(worst, sev_rank.get(sev, 0))
        payload.append({"check": getattr(f, "check", "?"), "severity": sev,
                        "message": getattr(f, "message", str(f))})
        if not as_json:
            print(f)
    verdict = {2: "REFUSE", 1: "REVIEW", 0: "PASS"}[worst]
    if as_json:
        print(json.dumps({"verdict": verdict, "findings": payload}, indent=2))
    else:
        print(f"\nverdict: {verdict}")
    return verdict


def code_for(verdict: str) -> int:
    """REFUSE is a failure exit. REVIEW is not - it needs a human, but it is not an error."""
    return 2 if verdict == "REFUSE" else 0


# --------------------------------------------------------------------------------------------
# commands


def cmd_validate(a) -> int:
    ing = load("ingest")
    sheet = a.samplesheet or (Path(a.project) / "samplesheet.csv" if a.project else None)
    if sheet is None:
        raise SystemExit("scqc validate: give --samplesheet FILE or --project DIR")
    rows = read_csv(Path(sheet))
    reg_path = Path(a.registry) if a.registry else ROOT / "references" / "_registry" / "registry.tsv"
    registry = ing.read_registry(reg_path)
    print(f"samplesheet : {sheet}  ({len(rows)} row(s))")
    print(f"registry    : {reg_path}  ({len(registry)} entr{'y' if len(registry) == 1 else 'ies'})")
    if not registry:
        print("  NOTE: the registry is empty or missing, so every declared reference will fail.")
    bad = 0
    results = []
    for r in rows:
        errs = ing.validate_row(r, registry)
        results.append({"sample": r.get("sample", "?"), "errors": errs})
        if errs:
            bad += 1
            if not a.json:
                print(f"\n[INVALID] {r.get('sample', '(no sample name)')}")
                for e in errs:
                    print(f"   - {e}")
        elif not a.json:
            print(f"[ok     ] {r.get('sample')}")
    factors = design_from(rows, skip={"sample", "platform", "species", "reference",
                                      "assay", "fastq_r1", "fastq_r2", "matrix"})
    if a.json:
        print(json.dumps({"rows": len(rows), "invalid": bad, "results": results,
                          "design_factors": {k: sorted(set(v.values())) for k, v in factors.items()}},
                         indent=2))
    else:
        print(f"\ndesign factors discovered: "
              f"{ {k: sorted(set(v.values())) for k, v in factors.items()} or 'none'}")
        if not factors:
            print("  Without a design factor the differential checks CANNOT RUN. They are the")
            print("  checks that catch a technical loss masquerading as biology.")
        print(f"\n{len(rows) - bad} of {len(rows)} row(s) valid")
    return 2 if bad else 0


def cmd_verify(a) -> int:
    v = load("verify_raw")
    verdict = v.verify(
        name=a.name, n_barcodes=a.barcodes, n_genes=a.genes, min_counts=a.min_counts,
        max_counts=a.max_counts, p98_counts=a.p98_counts,
        expected_genes=a.expected_genes, integer_counts=not a.non_integer)
    print(verdict)
    return 0 if verdict.usable else 2


def cmd_gate_cells(a) -> int:
    g = load("cellcall_gate")
    rows = read_csv(Path(a.calls))
    need = {"sample", "aligner", "cellbender"}
    missing = need - set(rows[0])
    if missing:
        raise SystemExit(f"scqc gate cells: {a.calls} is missing column(s): {sorted(missing)}\n"
                         f"  required: sample, aligner, cellbender; optional: lost")
    calls = {}
    for r in rows:
        aligner, cb = int(r["aligner"]), int(r["cellbender"])
        calls[r["sample"]] = {
            "aligner": aligner, "cellbender": cb,
            "lost": int(r["lost"]) if r.get("lost") not in (None, "") else max(aligner - cb, 0),
        }
    design = design_from(rows, skip=need | {"lost"})
    note(f"libraries: {len(calls)}   design factors: {list(design) or 'none'}\n", a.json)
    return code_for(emit(g.gate(calls, design), a.json))


def cmd_doublet_health(a) -> int:
    dh = load("doublet_health")
    rows = read_csv(Path(a.rates))
    if not {"sample", "rate"} <= set(rows[0]):
        raise SystemExit(f"scqc doublet health: {a.rates} needs columns: sample, rate")
    rates = {r["sample"]: float(r["rate"]) for r in rows}
    if any(v > 1 for v in rates.values()):
        rates = {k: v / 100.0 for k, v in rates.items()}
        note("rates looked like percentages; interpreted as such\n", a.json)
    design = design_from(rows, skip={"sample", "rate"})
    note(f"libraries: {len(rates)}   design factors: {list(design) or 'none'}\n", a.json)
    findings = dh.health(rates, design, n_kept_unscored=a.unscored,
                         detector_name=a.detector, reproducible=not a.non_reproducible)
    return code_for(emit(findings, a.json))


def cmd_quality(a) -> int:
    q = load("quality")
    rows = read_csv(Path(a.valleys))
    if not {"sample", "valley"} <= set(rows[0]):
        raise SystemExit(f"scqc quality: {a.valleys} needs columns: sample, valley"
                         f"  (optional: bimodal, true/false)")
    valleys = [q.Valley(r["sample"], a.metric, float(r["valley"]),
                        str(r.get("bimodal", "true")).strip().lower() not in ("false", "0", "no"))
               for r in rows]
    try:
        prop = q.derive(valleys, a.metric, light_floor=a.light_floor)
    except Exception as e:                                            # noqa: BLE001
        print(f"REFUSED: {e}")
        return 2
    print(prop)
    return 0


def cmd_cluster_preflight(a) -> int:
    ap = load("apply")
    rows = read_csv(Path(a.profile))

    def num(v):
        if v is None or str(v).strip() == "":
            return None                     # an unknown, never a zero
        try:
            return float(v)
        except ValueError:
            return None

    prof = []
    for r in rows:
        prof.append({
            "sample": r.get("sample"), "cluster": r.get("cluster"),
            "n": int(float(r["n"])) if r.get("n") else 0,
            "median_umi": num(r.get("median_umi")),
            "median_pct_mt": num(r.get("median_pct_mt")),
            "pct_doublet": num(r.get("pct_doublet")),
            "FLAG": str(r.get("FLAG", "")).strip().lower() in ("true", "1", "yes"),
            "WATCH": str(r.get("WATCH", "")).strip().lower() in ("true", "1", "yes"),
        })
    note(f"clusters: {len(prof)}   kept total: {a.kept:,}\n", a.json)
    findings = ap.preflight(prof, kept_total=a.kept)
    return code_for(emit(findings, a.json))


def cmd_selftest(a) -> int:
    suites = sorted((ROOT / "tests").glob("test_*.py")) + [ROOT / "tests" / "adversarial.py"]
    passed, failed, skipped = [], [], []
    for t in suites:
        if not t.exists():
            continue
        r = subprocess.run([sys.executable, str(t)], capture_output=True, text=True,
                           cwd=str(ROOT), encoding="utf-8", errors="replace")
        out = (r.stdout or "") + (r.stderr or "")
        # A suite that skipped exits 0. Counting that as a pass is precisely the defect this
        # project exists to catch - a check that never ran, reported as a check that succeeded.
        #
        # The marker has to be a WHOLE-SUITE one. Matching a bare "SKIP" anywhere caught
        # `[SKIP] <case>` lines, which several suites print as ordinary per-case output for a
        # case that exercises a skip path - so a fully-executed suite was reported as not run.
        # A suite announces its own skip with a line beginning `SKIP:`.
        did_skip = any(ln.strip().startswith("SKIP:") for ln in out.splitlines()) \
            or "ModuleNotFoundError" in out
        if r.returncode != 0:
            status, bucket = "FAIL", failed
        elif did_skip:
            status, bucket = "SKIP", skipped
        else:
            status, bucket = "PASS", passed
        bucket.append(t.name)
        print(f"  {status}  {t.name}")
        if status == "SKIP":
            why = next((ln.strip() for ln in out.splitlines() if ln.strip().startswith("SKIP:")), "needs pandas")
            print(f"        {why[:100]}")
        if status == "FAIL" and a.verbose:
            print("\n".join("        " + ln for ln in out.strip().splitlines()[-12:]))
    print(f"\n{len(passed)} passed, {len(failed)} failed, {len(skipped)} skipped")
    if skipped:
        print("A skipped suite has NOT been checked. It needs pandas (`pip install scqc[test]`)")
        print("or a cohort directory (COHORT_DIR); it is not evidence of anything either way.")
    return 1 if failed else 0


# --------------------------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="scqc", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Audit a single-cell/nuclei RNA-seq quality-control run.",
        epilog=(
            "scQC judges tables you produced elsewhere; it does not read matrices or run tools.\n"
            "Exit codes: 0 pass or review, 2 refuse, 1 error.\n"
            "See docs/PRINCIPLES.md for what each gate enforces and why."))
    p.add_argument("--version", action="version", version=f"scQC {VERSION}")
    sub = p.add_subparsers(dest="cmd", metavar="COMMAND")

    v = sub.add_parser("validate", help="check a samplesheet before anything is computed")
    v.add_argument("--project"); v.add_argument("--samplesheet"); v.add_argument("--registry")
    v.add_argument("--json", action="store_true"); v.set_defaults(fn=cmd_validate)

    r = sub.add_parser("verify", help="is a matrix raw? (from summary statistics you supply)")
    r.add_argument("--name", required=True)
    r.add_argument("--barcodes", type=int, required=True)
    r.add_argument("--genes", type=int, required=True)
    r.add_argument("--min-counts", dest="min_counts", type=float, required=True)
    r.add_argument("--max-counts", dest="max_counts", type=float, required=True)
    r.add_argument("--p98-counts", dest="p98_counts", type=float, required=True)
    r.add_argument("--expected-genes", dest="expected_genes", type=int)
    r.add_argument("--non-integer", action="store_true")
    r.set_defaults(fn=cmd_verify)

    g = sub.add_parser("gate-cells", help="did the denoiser drop cells the aligner kept?")
    g.add_argument("--calls", required=True,
                   help="CSV: sample,aligner,cellbender[,lost][,<design columns>]")
    g.add_argument("--json", action="store_true"); g.set_defaults(fn=cmd_gate_cells)

    d = sub.add_parser("doublet-health", help="is a doublet rate a measurement or the prior?")
    d.add_argument("--rates", required=True, help="CSV: sample,rate[,<design columns>]")
    d.add_argument("--detector", default="scDblFinder")
    d.add_argument("--unscored", type=int, default=0,
                   help="cells kept that were never scored (an unknown, not a zero)")
    d.add_argument("--non-reproducible", action="store_true")
    d.add_argument("--json", action="store_true"); d.set_defaults(fn=cmd_doublet_health)

    q = sub.add_parser("quality", help="propose a count floor from per-library valleys")
    q.add_argument("--valleys", required=True, help="CSV: sample,valley[,bimodal]")
    q.add_argument("--metric", choices=["umi", "genes"], required=True)
    q.add_argument("--light-floor", dest="light_floor", type=int, default=200)
    q.set_defaults(fn=cmd_quality)

    c = sub.add_parser("cluster-preflight", help="contradictions before a removal is applied")
    c.add_argument("--profile", required=True,
                   help="CSV: sample,cluster,n[,median_umi,median_pct_mt,pct_doublet,FLAG,WATCH]")
    c.add_argument("--kept", type=int, required=True, help="cells that would survive the filter")
    c.add_argument("--json", action="store_true"); c.set_defaults(fn=cmd_cluster_preflight)

    s = sub.add_parser("selftest", help="run the bundled test suites")
    s.add_argument("-v", "--verbose", action="store_true"); s.set_defaults(fn=cmd_selftest)
    return p


def main(argv=None) -> int:
    p = build_parser()
    a = p.parse_args(argv)
    if not getattr(a, "fn", None):
        p.print_help()
        return 1
    try:
        return a.fn(a)
    except SystemExit:
        raise
    except KeyboardInterrupt:
        return 130
    except Exception as e:                                            # noqa: BLE001
        print(f"scqc: {type(e).__name__}: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
