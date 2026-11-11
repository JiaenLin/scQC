# Runs the real CLI on a small synthetic single-cell cohort in temporary directories, with the stand-in
# doublet detector (tests/fixtures/fake_rscript.py), and calls the step-7 audit and Pipeline.gate
# directly; needs the analysis stack, not R.
"""A failed removal audit stops the run before any matrix is written (harness autodev W043, GAP-119).

Step 7's audit (modules/07_apply/audit_removal.py) checks that the deliverable is exactly the
recorded criteria applied jointly. It wrote FAIL, and Pipeline.gate stops on REFUSE only, so a
failed audit was recorded - the report counted it as an unrecognised severity under a REVIEW
verdict - and the objects were written anyway, while docs/WORKFLOW.md said the caller stops.
Checked here:

  1  the audit's verdict on a removal table whose keep and removed columns disagree is REFUSE, and
     Pipeline.gate raises on it; on a correct table it does not
  2  through the real CLI, apply mode: a UMI floor recorded in decisions.yml below the light floor
     (every called cell sits above both, so no unexamined cell is kept today - the inequality is
     what guarantees it on the next cohort) fails 'the light floor is strictly below the quality
     floor'. The run stops at 07_apply with exit 2, STATUS.json's refusal names the step and the
     check, no deliverable object, per-library object or write record exists, and the ledger does
  3  the same refusal is lifted the way every refusal is: an accept entry naming 07_apply, the
     check, a person and their words completes the run, recorded in STATUS.json's escapes and in
     the report; half an approval (a bare "ok") accepts nothing and the refusal says why
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
fails: list[str] = []
CHECK = "step 3: the light floor is strictly below the quality floor"


def module(rel: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m                   # a dataclass module must be registered to be built
    spec.loader.exec_module(m)
    return m


# ---- 1. the audit's verdict, and the gate acting on it
ar = module("modules/07_apply/audit_removal.py", "w043_audit_removal")
ROWS = {"aligner_cell": [1, 1, 1, 0], "keep": [1, 1, 0, 0], "removed": [0, 0, 1, 1],
        "fail_not_aligner_cell": [0, 0, 0, 1], "fail_umi_floor": [0, 0, 1, 0]}
BROKEN = dict(ROWS, removed=[1, 0, 1, 1])                # row 0 is both kept and removed
CRIT = ["fail_not_aligner_cell", "fail_umi_floor"]
try:
    from engine.executor import make_executor
    from engine.pipeline import Pipeline
    from engine.task import Refusal
    HAVE_ENGINE = True
except ImportError as e:
    HAVE_ENGINE = False
    print(f"SKIP: the engine is not importable here ({e})")

if HAVE_ENGINE:
    for label, rows, stops in (("a correct removal", ROWS, False),
                               ("keep and removed disagreeing", BROKEN, True)):
        f = ar.audit(rows, criteria=CRIT, cell_col="aligner_cell")
        v = ar.verdict(f)
        with tempfile.TemporaryDirectory() as td:
            p = Pipeline(project=Path(td) / "proj", mode="apply", executor=make_executor("local"),
                         samples=[{"sample": "S1", "assay": "scrna"}], decisions={}, jobs=1)
            try:
                p.gate("07_apply", f, v)
                raised = None
            except Refusal as e:
                raised = str(e)
            if stops and (v != "REFUSE" or raised is None or "complementary" not in raised):
                fails.append(f"1: {label}: verdict {v!r}, and the gate "
                             + ("did not stop the run" if raised is None else f"said {raised[:80]!r}"))
            if not stops and raised is not None:
                fails.append(f"1: {label} was refused: {raised[:120]}")
            if stops:
                got = [d["severity"] for d in p.findings if "complementary" in d["check"]]
                if got != ["REFUSE"]:
                    fails.append(f"1: {label} is recorded as {got}, not REFUSE")
print("1. a failed audit is REFUSE, and the gate stops on it; a correct removal passes")

# ---- 2 and 3. through the real CLI
try:
    import anndata  # noqa: F401
    import yaml
    from fixtures import cohorts as C
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")

COLS = ["sample", "platform", "species", "reference", "assay", "matrix", "aligner_cells",
        "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]
FLOOR = {"umi_floor": {"value": 150, "approved_by": "the analyst",
                       "verbatim": "a floor of 150 UMI keeps every cell the aligner called"}}
WORDS = ("every called cell is above the light floor here, so no unexamined cell can be kept; "
         "the floor is recorded on purpose")


def cohort(td: Path) -> Path:
    rows = []
    for i, s in enumerate(("LIB_A", "LIB_B")):
        raw, filt, _ = C.cellranger_library(td / "data", s, 40 + i, C.HUMAN)
        rows.append({"sample": s, "platform": "10x", "species": "homo_sapiens",
                     "reference": "homo_sapiens/ref", "assay": "scrna", "matrix": raw,
                     "aligner_cells": filt, "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]",
                     "dbr": 0.05, "dbr_sd": 0.05})
    sheet = td / "samplesheet.tsv"
    sheet.write_text("\t".join(COLS) + "\n"
                     + "".join("\t".join(str(r[c]) for c in COLS) + "\n" for r in rows))
    (td / "registry.tsv").write_text("species\tbuild\tpath\nhomo_sapiens\tref\t-\n")
    return sheet


def run(td: Path, sheet: Path, label: str, accept: list | None):
    proj = td / f"proj_{label}"
    d = td / f"decisions_{label}.yml"
    d.write_text(yaml.safe_dump({"quality": FLOOR, **({"accept": accept} if accept else {})}))
    cmd = [sys.executable, str(ROOT / "scqc_cli.py"), "run", "--project", str(proj),
           "--samplesheet", str(sheet), "--registry", str(td / "registry.tsv"), "--mode", "apply",
           "--executor", "local", "--python", sys.executable,
           "--rscript", str(ROOT / "tests" / "fixtures" / "fake_rscript.py"),
           "--jobs", "2", "--seed", "0", "--decisions", str(d)]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    res = [x for x in (proj / "results").iterdir() if x.is_dir() and not x.is_symlink()]
    return p.returncode, p.stdout + p.stderr, (res[0] if res else None)


def status(R) -> dict:
    return json.loads((R / "STATUS.json").read_text()) if R and (R / "STATUS.json").is_file() else {}


WRITTEN = ("cohort.deliverable.h5ad", "cohort_per_sample", "cohort.apply_write.metrics.json")

if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sheet = cohort(td)
        rc, log, R = run(td, sheet, "refused", None)
        why = str((status(R).get("refusal") or {}).get("reason"))
        if rc != 2 or "07_apply" not in why or CHECK not in why:
            fails.append(f"2: a failed removal audit did not stop the run at 07_apply (exit {rc}; "
                         f"refusal {why[:120]!r}):\n" + "\n".join(log.splitlines()[-8:]))
        made = [w for w in WRITTEN if R is not None and (R / "objects" / w).exists()]
        if made:
            fails.append(f"2: the refused run wrote {', '.join(made)}")
        if R is None or not (R / "tables" / "removal_ledger.csv").is_file():
            fails.append("2: the refused run left no removal ledger, which the docs say it writes")
        print("2. a failed removal audit stops the run at 07_apply, before any object is written")

        rc, log, R = run(td, sheet, "half", [{"step": "07_apply", "check": CHECK,
                                              "approved_by": "the PI", "verbatim": "ok"}])
        why = str((status(R).get("refusal") or {}).get("reason"))
        if rc != 2 or "does not count" not in why:
            fails.append(f"3: half an approval (verbatim 'ok') lifted the audit's refusal, or the "
                         f"refusal did not say why it does not count (exit {rc})")
        rc, log, R = run(td, sheet, "accepted", [{"step": "07_apply", "check": CHECK,
                                                  "approved_by": "the PI", "verbatim": WORDS}])
        if rc != 0 or R is None or not (R / "objects" / "cohort.deliverable.h5ad").exists():
            fails.append(f"3: with a complete accept entry the run did not deliver (exit {rc}):\n"
                         + "\n".join(log.splitlines()[-8:]))
        else:
            esc = [e for e in status(R).get("escapes") or []
                   if (e.get("ask") or {}).get("check") == CHECK]
            if not esc or (esc[0].get("decision") or {}).get("why") != WORDS:
                fails.append(f"3: STATUS.json's escapes do not carry the acceptance: "
                             f"{status(R).get('escapes')}")
            if "ACCEPTED by the PI" not in (R / "reports" / "payload.json").read_text():
                fails.append("3: the report does not say the audit's refusal was accepted, or by whom")
        print("3. an accept entry lifts it by name and is recorded; half an approval does not")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("a failed removal audit stops the run, and only a named acceptance lifts it")
