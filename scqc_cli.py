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
import os
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
    """Read a delimited table, skipping comment lines. Returns a list of dicts.

    THE DELIMITER IS DETECTED, AND AN UNDETECTABLE ONE IS REFUSED.

    This read comma-separated only. Handed a tab-separated samplesheet - the extension this
    project's own docs use - csv.DictReader did not fail: it produced one column per row whose
    single key was the entire header line. The run then reported "samples: 10", started stage 1,
    and died on `KeyError: 'sample'` inside the task graph, three layers from the cause.

    That is the failure worth engineering against. A wrong delimiter does not produce a parse
    error, it produces a table of the right LENGTH and the wrong SHAPE, and every check that
    counts rows agrees with it.
    """
    if not path.exists():
        raise SystemExit(f"scqc: no such file: {path}")
    with path.open(encoding="utf-8", newline="") as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    if not lines:
        raise SystemExit(f"scqc: {path} contains no data")

    header = lines[0].rstrip("\r\n")
    counts = {d: header.count(d) for d in ("\t", ",", ";")}
    best = max(counts, key=lambda d: counts[d])
    if counts[best] == 0:
        raise SystemExit(
            f"scqc: {path} has a single-column header ({header[:60]!r}...). No tab, comma or "
            f"semicolon separates it, so there is nothing to read columns from.")
    # Ambiguity is refused rather than guessed. A header carrying both tabs and commas in
    # comparable numbers parses two different ways into two different tables, and picking one
    # silently is how the wrong shape gets used.
    rivals = [d for d, n in counts.items() if d != best and n >= counts[best]]
    if rivals:
        names = {"\t": "tab", ",": "comma", ";": "semicolon"}
        raise SystemExit(
            f"scqc: {path} is ambiguous - its header contains "
            + ", ".join(f"{n} {names[d]}(s)" for d, n in counts.items() if n)
            + ". Quote the fields or use one delimiter; guessing produces a table of the right "
              "length and the wrong shape.")

    rows = list(csv.DictReader(lines, delimiter=best))
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


def cmd_stamp(a) -> int:
    """Add `uns["scqc"]` to objects written before this pipeline declared what its flag means.

    The declaration exists so a downstream tool can act on `cluster_FLAG` knowing whose decision
    it is. Objects written earlier carry the flag and not the declaration, and re-running a
    cohort to add one would be hours of compute to write a dictionary. This adds it in place.

    What it will NOT do is invent the run identity. `--run-key` and `--commit` are yours to
    supply and default to empty; a stamped object that names a run it did not come from is worse
    than one that names none, because the wrong provenance is the kind a reader trusts.

    The digest is computed HERE, from the column as it stands in the file. So a stamp added
    afterwards describes the object in hand rather than the object scQC once wrote - and if the
    flag has been altered since, the stamp records the altered column and says nothing false.
    """
    from pathlib import Path as _P

    try:
        import anndata as ad
    except ImportError:
        print("scqc stamp: needs anndata.  Use the core environment's interpreter.",
              file=sys.stderr)
        return 1
    sys.path.insert(0, str(_P(__file__).resolve().parent))
    from adapters import declaration as decl

    version = ""
    try:
        version = (_P(__file__).resolve().parent / "VERSION").read_text(encoding="utf-8").strip()
    except OSError:
        pass

    rc = 0
    for raw in a.objects:
        p = _P(raw)
        if not p.exists():
            print(f"  MISSING  {p}", file=sys.stderr)
            rc = 1
            continue
        A = ad.read_h5ad(p)
        existing = A.uns.get(decl.KEY)
        if existing and not a.force:
            ok, why = decl.verify(A)
            print(f"  already declared  {p.name}  ({why})")
            print("                    --force re-stamps it")
            continue
        d = decl.build(A, sample=str(existing.get("sample", "") if existing else ""),
                       run_key=a.run_key, commit=a.commit, version=version)
        flag = d["flag_column"] or "no flag column"
        n = d["n_flagged"]
        print(f"  {'would stamp' if a.dry_run else 'stamped'}  {p.name}  "
              f"{A.n_obs:,} obs  {flag}"
              + (f"  {n:,} flagged  digest {d['flag_digest']}" if d["flag_column"] else ""))
        if not a.dry_run:
            A.uns[decl.KEY] = d
            A.write_h5ad(p)
    if a.dry_run:
        print("\n  --dry-run: nothing was written.")
    return rc


def cmd_report(a) -> int:
    """Rebuild a finished run's report from the files in its own results directory.

    The pipeline writes the report at the end of a run, which means a change to the report or to
    the figures used to need the whole run again - hours, to redraw a plot. This reads the
    tables the run already wrote and rebuilds from them, so the report is iterable on its own.

    It rebuilds ONLY what the directory can support: the payload's steps and per-library table
    come from the existing report.json, and the figures from the tables. Nothing is recomputed
    from the matrices, so a number in the rebuilt report is the same number as in the original.
    """
    from report.build import build_report
    from report.collect import collect as collect_figures

    results = Path(a.results).expanduser().resolve()
    reports, tables = results / "reports", results / "tables"
    source = reports / "payload.json"
    if not tables.is_dir():
        print(f"scqc: {tables} does not exist - not a results directory from `scqc run`")
        return 2
    if not source.exists():
        # NOT report.json. That file is the rendered document, and building a document from a
        # document silently produces a report with empty steps rather than an error.
        print(f"scqc: {source} does not exist, so there is nothing to rebuild from. It is "
              f"written by the report step; a run finished before that step wrote it has only "
              f"report.json, which is the rendered document and not a payload.")
        return 2

    # Read back, not rebuilt: this is the payload report_payload() produced during the run, so
    # the document below is the same one that run would have written, with the figures added.
    payload = json.loads((reports / "payload.json").read_text(encoding="utf-8"))
    figures, notes = collect_figures(tables, samplesheet=a.samplesheet)
    payload["figures"] = figures
    payload["figure_notes"] = notes
    print(f"figures assembled: {len(figures)}  ({', '.join(sorted(figures)) or 'none'})")
    for fid in sorted(notes):
        print(f"  {fid} NOT PRODUCED - {notes[fid]}")
    source.write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")
    doc = build_report(payload, reports / "qc_report.html", reports / "report.json")
    # doc["figures"] IS the render index - {id: {"rendered": bool, "source": str}} - not a block
    # containing one. Reading it as though it were nested reports "rendered 0" for a build that
    # drew every figure it was given, which is a lie in the direction nobody checks.
    index = doc.get("figures") or {}
    rendered = sum(1 for v in index.values() if v.get("rendered"))
    print(f"rendered {rendered} of {len(index)} figure(s) into {reports / 'qc_report.html'}")
    for fid, entry in sorted(index.items()):
        if not entry.get("rendered"):
            print(f"  {fid} did not draw: {entry.get('source')}")
    return 0


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


def cmd_run(a) -> int:
    """Run the pipeline over a project, in two stages.

    Stage one is ingest alone, because what it decides - accept this matrix, or rebuild it from
    FASTQ - determines which tasks stage two contains. Building the whole graph first would put
    an alignment in it that may not be needed, or omit one that is, and the manifest would record
    both outcomes identically.

    A refusal from any gate stops the run and is reported as a refusal, not as an error: it is a
    correct outcome that needs reading, and collapsing it into a failure teaches a user to treat
    it as flaky infrastructure. The report is written either way, because a run that stopped is
    exactly the run whose reasons somebody needs.
    """
    sys.path.insert(0, str(ROOT))
    from engine import graph
    from engine.decisions import load as load_decisions
    from engine.executor import make_executor
    from engine.pipeline import Pipeline
    from engine.task import Refusal, Status, first_line

    # Argument contradictions are settled BEFORE any file is opened. Otherwise a misuse of the
    # flags surfaces as a missing-file error about something else entirely, and the reader fixes
    # the wrong thing.
    if a.mode == "evidence" and a.decisions:
        raise SystemExit(
            "scqc run: --decisions is only read in apply mode.\n"
            "       Evidence mode exists to produce the evidence those decisions are made\n"
            "       against; supplying them here would mean they were chosen beforehand.")

    project = Path(a.project).expanduser().resolve()
    sheet = Path(a.samplesheet) if a.samplesheet else project / "samplesheet.csv"
    rows = read_csv(sheet)

    decisions = {}
    if a.mode == "apply":
        # Optional. Apply mode runs on the thresholds the pipeline derived and records them as
        # DERIVED; a decisions file, when present, overrides them and is recorded as ADJUDICATED.
        # The distinction survives into the deliverable, so nothing later reads a proposal as a
        # decision.
        dpath = Path(a.decisions) if a.decisions else project / "decisions.yml"
        if a.decisions or dpath.exists():
            decisions = load_decisions(dpath)

    executor = make_executor(a.executor, **({"queue": a.queue, "project": a.pbs_project}
                                            if a.executor == "pbs" else {}))
    jobs = a.jobs if a.jobs and a.jobs > 0 else min(16, (os.cpu_count() or 4))
    tools = {k: v for k, v in {
        "celescope": a.celescope, "cellranger": a.cellranger, "cellbender": a.cellbender,
        "rscript": a.rscript, "device": ("cpu" if a.cpu else "cuda"),
        "light_floor": a.light_floor, "seed": a.seed, "resolution": a.resolution,
        # DECLARED, with no default anywhere. The doublet adapter refuses a missing dbr rather
        # than substituting one, because an expected doublet rate is a property of how the
        # libraries were loaded and cannot be read off the data. It reaches the graph either
        # per-sample from the samplesheet - loading concentration can differ between libraries -
        # or cohort-wide from these flags.
        "dbr": a.dbr, "dbr_sd": a.dbr_sd,
        # Split here rather than in the graph, so a malformed list is a usage error before any
        # file is opened. An empty token is dropped: `default,,1` is a typo, not a request to
        # sweep a setting with no name, and `resolve_dbr_sd` would refuse it three steps later.
        "dbr_sd_sweep": ([t.strip() for t in str(a.dbr_sd_sweep).split(",") if t.strip()]
                         if a.dbr_sd_sweep else None),
        # Same treatment, same reason: a malformed list is a usage error before any file is
        # opened. The default resolution is dropped if it appears here, so `--resolution 2
        # --extra-resolutions 1,2,3` profiles 2 once rather than twice under two names.
        "extra_resolutions": (sorted({r for r in
                                      (float(t.strip()) for t in
                                       str(a.extra_resolutions).split(",") if t.strip())
                                      if r != float(a.resolution)})
                              if a.extra_resolutions else None),
    }.items() if v is not None}
    # AFTER `tools`, because the run's output directory is named from the parameters it was given.
    pipe = Pipeline(project=project, mode=a.mode, executor=executor, samples=rows,
                    decisions=decisions, force=a.force, jobs=jobs, tools=tools)
    python_exe = a.python or sys.executable

    # Which code is running, printed by scQC rather than by whatever script launched it. The job
    # script echoed `git log --oneline -1`, and a compute node has no git: the banner read
    # `code :` with nothing after it, on every run that mattered.
    from engine.provenance import git_provenance
    _g = git_provenance(ROOT)
    print(f"code     : {_g.get('commit')}"
          + (f" ({_g['branch']})" if _g.get("branch") else "")
          + ("  tree: MODIFIED" if _g.get("dirty") else
             ("  tree: clean" if _g.get("dirty") is False else "  tree: not determined")))
    print(f"project  : {project}")
    print(f"mode     : {a.mode}      executor: {getattr(executor, 'name', '?')}")
    print(f"samples  : {len(rows)}")
    print(f"jobs     : {jobs} concurrent")
    print()

    def finish(code: int, stopped=None) -> int:
        # UNDER THE RUN'S OWN DIRECTORY, not beside it. This wrote to
        # `<project>/results/reports/`, outside the content-addressed directory the rest of the
        # run writes into - so the report escaped the layout that guarantees nothing is
        # overwritten, and every subsequent run replaced it whatever its parameters had been. The
        # one artifact a reader opens first was the one artifact with no run identity.
        html = pipe.results / "reports" / "qc_report.html"
        html.parent.mkdir(parents=True, exist_ok=True)
        try:
            from report.build import build_report
            # The reason the run stopped is the single most useful line in the document, so it
            # is passed rather than left None. A report that records a stop without its cause
            # sends the reader to the logs for the one thing the report exists to tell them.
            if stopped is None:
                bad = [f"{k}: {r.status.value}" for k, r in sorted(pipe.results_by_key.items())
                       if r.status.value in ("refused", "failed")]
                stopped = "; ".join(bad) or None
            payload = pipe.report_payload(stopped=stopped)
            # THE PAYLOAD IS REWRITTEN WITH THE DOCUMENT BUILT FROM IT.
            #
            # This wrote the HTML and the JSON and left `payload.json` as whatever the report TASK
            # had written - and on a resumed run that task is SKIPPED, so the payload on disk was
            # the previous version's while the report beside it was this one's. `scqc report`
            # rebuilds from `payload.json`, so running it would have regenerated an OLDER document
            # over a newer one: on this cohort, a report with no parameter table and freshness
            # NOT CHECKED, replacing one that had both. Staleness with no symptom, between two
            # files that are meant to be two views of the same thing.
            (pipe.results / "reports" / "payload.json").write_text(
                json.dumps(payload, indent=1, default=str), encoding="utf-8")
            build_report(payload, html, pipe.results / "reports" / "report.json")
            print(f"\nreport   : {html}")
        except Exception as e:                                        # noqa: BLE001
            print(f"\nreport   : COULD NOT BE WRITTEN - {type(e).__name__}: {e}",
                  file=sys.stderr)
            print("           The run's own record is missing; treat its outcome as unrecorded.",
                  file=sys.stderr)
        return code

    try:
        print("stage 1 - ingest")
        pipe.run(graph.ingest_stage(pipe, python_exe))
        bad = {k: r for k, r in pipe.results_by_key.items() if not r.ok}
        if bad:
            for k, r in sorted(bad.items()):
                print(f"  {r.status.value.upper():8s} {k}: {first_line(r.message, 110)}")
            return finish(2)

        ingest = {r.key.split("/", 1)[1]: r.metrics
                  for r in pipe.results_by_key.values() if r.key.startswith("00_ingest/")}
        print("\nstage 2 - everything decided by stage 1")
        pipe.run(graph.main_stage(pipe, python_exe, tools, ingest))
    except Refusal as e:
        print(f"\nREFUSED\n{e}")
        return finish(2)

    refused = [k for k, r in pipe.results_by_key.items() if r.status is Status.REFUSED]
    failed = [k for k, r in pipe.results_by_key.items() if r.status is Status.FAILED]
    print()
    for label, keys in (("refused", refused), ("failed", failed)):
        if keys:
            print(f"  {label}: {', '.join(sorted(keys))}")
    return finish(2 if refused else (1 if failed else 0))


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

    r2 = sub.add_parser("run", help="run the pipeline over a project")
    r2.add_argument("--project", required=True)
    r2.add_argument("--mode", choices=["evidence", "apply"], default="apply",
                    help="apply (default) derives the thresholds, applies them and writes the "
                         "deliverable; evidence measures and applies nothing")
    r2.add_argument("--samplesheet")
    r2.add_argument("--decisions", help="apply mode only; evidence mode refuses it")
    r2.add_argument("--executor", choices=["local", "pbs"], default="local")
    r2.add_argument("--queue"); r2.add_argument("--pbs-project", dest="pbs_project")
    r2.add_argument("--python", help="interpreter that has scanpy/anndata")
    r2.add_argument("--celescope"); r2.add_argument("--cellranger")
    r2.add_argument("--cellbender"); r2.add_argument("--rscript")
    r2.add_argument("--cpu", action="store_true", help="no GPU available")
    r2.add_argument("--light-floor", dest="light_floor", type=int, default=200)
    # No defaults. A doublet rate nobody chose is the threshold this pipeline most objects to,
    # and scDblFinder's own default is derived from a 10x loading model that does not transfer.
    r2.add_argument("--dbr", type=float,
                    help="expected doublet rate (DECLARED; per-sample `dbr` in the samplesheet "
                         "overrides it). No default: it describes how the libraries were loaded.")
    r2.add_argument("--dbr-sd", dest="dbr_sd", type=float,
                    help="uncertainty on --dbr (DECLARED; per-sample `dbr_sd` overrides it)")
    # OFF by default because it re-scores every library once per setting. It changes no
    # deliverable; it is the evidence behind figure F5, which asks whether the rate the run
    # applied was measured or was the prior's, and one setting cannot answer that.
    r2.add_argument("--dbr-sd-sweep", dest="dbr_sd_sweep",
                    help="comma-separated dbr.sd settings to sweep for figure F5, e.g. "
                         "'default,dbr,1'. 'default' omits the argument, 'dbr' sets dbr.sd = dbr, "
                         "anything else is a number. Costs one extra scDblFinder run per library "
                         "per setting and applies nothing.")
    # 2.0, not 1.0. At 1.0 a droplet-based library of this size resolves to major lineages only,
    # and criterion D - a cluster's doublet frequency - cannot fire, because a doublet pocket is
    # small and dissolves into the lineage cluster around it. Measured on the calibration cohort:
    # D fired on nothing at 1.0 or 1.5 and on two clusters at 2.0, both above 72% doublet.
    r2.add_argument("--resolution", type=float, default=2.0,
                    help="the clustering resolution step 6 flags on, and the one every downstream "
                         "artifact carries (default 2.0)")
    # Profiled and flagged ALONGSIDE the default, into sibling files. They change no deliverable
    # and nothing downstream reads them; they exist so the flags at one resolution can be read
    # against the flags at another, which is the only way to tell a contamination POCKET (isolates
    # as resolution rises) from a contamination GRADIENT (never isolates, at any resolution).
    r2.add_argument("--extra-resolutions", dest="extra_resolutions", default="1.0,3.0",
                    help="comma-separated additional resolutions to profile and flag beside the "
                         "default, into sibling files nothing downstream reads "
                         "(default '1.0,3.0'; pass '' for none)")
    r2.add_argument("--seed", type=int, default=0)
    r2.add_argument("--force", action="store_true", help="re-run every task, ignoring the manifest")
    # Independent tasks run concurrently. 1 is the old serial behaviour, which is what you want
    # when a failure has to be read in one log; the default uses the machine.
    r2.add_argument("--jobs", "-j", type=int, default=0,
                    help="independent tasks to run at once (0 = auto: cores, capped at 16)")
    r2.set_defaults(fn=cmd_run)

    rp = sub.add_parser("report",
                        help="rebuild the report of a finished run, from that run's own files")
    rp.add_argument("results", help="a results/<digest> directory written by `scqc run`")
    rp.add_argument("--samplesheet", default=None,
                    help="the run's samplesheet, for the design panel of F2")
    rp.set_defaults(fn=cmd_report)

    st = sub.add_parser("stamp",
                        help="add the uns['scqc'] declaration to objects already written")
    st.add_argument("objects", nargs="+", help=".h5ad file(s) this pipeline produced")
    st.add_argument("--run-key", default="", help="the run that produced them, if known")
    st.add_argument("--commit", default="", help="the commit that produced them, if known")
    st.add_argument("--force", action="store_true",
                    help="re-stamp an object that already carries a declaration")
    st.add_argument("--dry-run", action="store_true", help="report, write nothing")
    st.set_defaults(fn=cmd_stamp)

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
        # The TRACEBACK is printed, not just the type and message.
        #
        # Refusals and SystemExit are the designed, user-facing errors and they say what to do.
        # Anything reaching here is a BUG, and for a bug the message alone is often useless:
        # `KeyError: 'sample'` names a key that appears in a dozen places and points at none of
        # them. A report that omits what a reader needs in order to act sends them somewhere
        # else - a whole diagnosis was spent on a healthy environment for exactly this reason.
        import traceback
        print(f"scqc: {type(e).__name__}: {e}", file=sys.stderr)
        print("scqc: this is a bug in the pipeline, not a refusal. Traceback:", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
