# Runs the real CLI on small synthetic single-cell cohorts in temporary directories, with the stand-in
# doublet detector (tests/fixtures/fake_rscript.py), and calls two gate modules directly; needs the
# analysis stack, not R.
"""A refusal can be accepted, and a silence is judged against the library's own declaration
(harness autodev W020, B1-103).

Every refusal said "record a decision that accepts it", and no code read one: a real 3x biological
difference between two arms stopped the run for good, against the project's own rule that every
deliberate bypass is possible and logged (docs/PRINCIPLES.md section 3). Checked here:

  1  through the real CLI, apply mode: two genotypes loaded at different densities give a 3.5x
     doublet differential and the step-4 gate REFUSES; with an `accept` entry in decisions.yml
     naming that check, a person and their words, the run completes, the finding reads ACCEPTED,
     and STATUS.json's `escapes` carries the acceptance; half an approval (a bare "ok") accepts
     nothing and the refusal says why
  2  the cell-call gate: a denoiser calling fewer cells than the aligner while losing under the
     REVIEW line of the aligner's cells is REVIEW; losing more is still REFUSED
  3  a library under the silent-library line is REFUSED only where its declared rate predicted
     doublets: too few scored cells to expect one, or a rate declared under the line, is REVIEW
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
fails: list[str] = []


def module(rel: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m                   # a dataclass module must be registered to be built
    spec.loader.exec_module(m)
    return m


# ---- 2. the cell-call gate's ratio check
cg = module("modules/02_cells/cellcall_gate.py", "w020_cellcall_gate")
for lost, want in ((200, "REVIEW"), (600, "REFUSE")):
    calls = {"A": {"aligner": 10_000, "cellbender": 9_900, "lost": lost},
             "B": {"aligner": 10_000, "cellbender": 14_000, "lost": 0}}
    f = next(x for x in cg.gate(calls, assay="snrna")
             if x.check.startswith("CellBender must not be stricter"))
    if f.severity != want:
        fails.append(f"2: stricter than the aligner, losing {lost / 100:.0f}%: {f.severity}, "
                     f"expected {want}")
print("2. stricter than the aligner is REVIEW under the loss line, REFUSE above it")

# ---- 3. a silence judged against the library's declaration
dh = module("modules/04_doublets/doublet_health.py", "w020_doublet_health")
CASES = [
    # (label, rate, declared, n scored, expected severity of 'silent libraries')
    ("seven cells at a declared 5%", 0.0, 0.05, 7, "REVIEW"),
    ("a rate declared under the line", 0.003, 0.004, 10_000, "REVIEW"),
    ("ten thousand cells at a declared 6%", 0.003, 0.06, 10_000, "REFUSE"),
]
for label, rate, declared, n, want in CASES:
    rates = {"S": rate, "T": 0.05}
    try:
        got = next(x for x in dh.health(rates, declared_rate={"S": declared, "T": 0.05},
                                        n_scored={"S": n, "T": 5_000}, assay="snrna")
                   if x.check == "silent libraries").severity
    except TypeError as e:
        got = f"not judged ({e})"
    if got != want:
        fails.append(f"3: {label}: 'silent libraries' is {got}, expected {want}")
old = next(x for x in dh.health({"S": 0.0, "T": 0.05}, assay="snrna")
           if x.check == "silent libraries").severity
if old != "REFUSE":
    fails.append(f"3: with no declaration to judge against, a silent library is {old}, not REFUSE")
print("3. a silence is refused only where the library's own declaration expected doublets")

# ---- 1. through the real CLI
try:
    import anndata  # noqa: F401
    from fixtures import cohorts as C
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")

COLS = ["sample", "platform", "species", "reference", "assay", "matrix", "aligner_cells",
        "genotype", "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]
CHECK = "design differential: genotype"


def cohort(td: Path) -> Path:
    rows = []
    for i, (s, g, dbr) in enumerate((("CT_1", "ctl", 0.02), ("TR_1", "trt", 0.07),
                                     ("CT_2", "ctl", 0.02), ("TR_2", "trt", 0.07))):
        raw, filt, _ = C.cellranger_library(td / "data", s, 30 + i, C.HUMAN)
        rows.append({"sample": s, "platform": "10x", "species": "homo_sapiens",
                     "reference": "homo_sapiens/ref", "assay": "scrna", "matrix": raw,
                     "aligner_cells": filt, "genotype": g, "mt_prefix": "MT-",
                     "ribo_pattern": "^RP[SL]", "dbr": dbr, "dbr_sd": dbr})
    sheet = td / "samplesheet.tsv"
    sheet.write_text("\t".join(COLS) + "\n"
                     + "".join("\t".join(str(r[c]) for c in COLS) + "\n" for r in rows))
    (td / "registry.tsv").write_text("species\tbuild\tpath\nhomo_sapiens\tref\t-\n")
    return sheet


def run(td: Path, sheet: Path, label: str, accept: list | None):
    proj = td / f"proj_{label}"
    cmd = [sys.executable, str(ROOT / "scqc_cli.py"), "run", "--project", str(proj),
           "--samplesheet", str(sheet), "--registry", str(td / "registry.tsv"), "--mode", "apply",
           "--executor", "local", "--python", sys.executable,
           "--rscript", str(ROOT / "tests" / "fixtures" / "fake_rscript.py"),
           "--jobs", "4", "--seed", "0"]
    if accept is not None:
        import yaml
        d = td / f"decisions_{label}.yml"
        d.write_text(yaml.safe_dump({"accept": accept}))
        cmd += ["--decisions", str(d)]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    res = [x for x in (proj / "results").iterdir() if x.is_dir() and not x.is_symlink()]
    return p.returncode, p.stdout + p.stderr, (res[0] if res else None)


def refusal_words(R) -> str:
    """What STATUS.json says the refusal was - the gate's own words, where a host reads them."""
    if R is None or not (R / "STATUS.json").is_file():
        return ""
    return str((json.loads((R / "STATUS.json").read_text()).get("refusal") or {}).get("reason"))


if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sheet = cohort(td)
        rc, log, R = run(td, sheet, "none", None)
        if rc != 2 or CHECK not in refusal_words(R):
            fails.append(f"1: without an accept entry the {CHECK!r} REFUSE did not stop the run "
                         f"(exit {rc})")
        print("1a. a 3.5x doublet differential refuses without an accept entry")

        words = ("the two genotypes were loaded at different densities on purpose, so their "
                 "doublet rates differ by design")
        rc, log, R = run(td, sheet, "half", [{"step": "04_doublets", "check": CHECK,
                                              "approved_by": "the PI", "verbatim": "ok"}])
        if rc != 2 or "does not count" not in refusal_words(R):
            fails.append(f"1: half an approval (verbatim 'ok') was taken as one, or the refusal "
                         f"did not say why it does not count (exit {rc})")
        print("1b. half an approval accepts nothing, and says why")

        rc, log, R = run(td, sheet, "accepted", [{"step": "04_doublets", "check": CHECK,
                                                  "approved_by": "the PI", "verbatim": words}])
        if rc != 0 or R is None:
            fails.append(f"1: with a complete accept entry the run did not complete (exit {rc}):\n"
                         + "\n".join(log.splitlines()[-12:]))
        else:
            st = json.loads((R / "STATUS.json").read_text())
            esc = [e for e in st.get("escapes") or []
                   if (e.get("ask") or {}).get("check") == CHECK]
            if not esc or (esc[0].get("decision") or {}).get("by") != "the PI" \
                    or (esc[0].get("decision") or {}).get("why") != words:
                fails.append(f"1: STATUS.json's escapes do not carry the acceptance: "
                             f"{st.get('escapes')}")
            pay = json.loads((R / "reports" / "payload.json").read_text())
            text = json.dumps(pay)
            if "ACCEPTED by the PI" not in text:
                fails.append("1: the report does not say the refusal was accepted, or by whom")
            if not (R / "objects" / "cohort.deliverable.h5ad").exists():
                fails.append("1: the accepted run wrote no deliverable")
        print("1c. with an accept entry the run completes, and STATUS.json and the report say so")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("a refusal can be accepted, by a named person in their own words, and nothing else lifts it")
