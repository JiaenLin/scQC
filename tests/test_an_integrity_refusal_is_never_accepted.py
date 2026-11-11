# Runs the real CLI on a small synthetic single-cell cohort in temporary directories, with the stand-in
# doublet detector (tests/fixtures/fake_rscript.py), and calls the step-7 audit directly; needs the
# analysis stack, not R.
"""An integrity check of step 7's removal audit is never accepted by name (harness autodev W043, fix round 1).

W043 made a failed removal audit stop the run, and an accept entry matching step and check name
lifted any of its refusals. Seven of the eight compare the per-cell tables with themselves, and a
table 07_measure wrote cannot fail one. An entry accepting "keep and removed are complementary",
written while it held, lifted it on the next run after one library's table was edited (a removed
doublet kept, a kept cell dropped, the count unchanged): exit 0, and the deliverable held a cell its
own ledger lists as removed. Checked here:

  1  every REFUSE the audit can emit is named in INTEGRITY or in ACCEPTABLE, which is the floor
     ordering alone; integrity() raises on a failed integrity check whatever is accepted, refusing
     the entry by name and printing no entry to paste, and reports an entry naming one that held
  2  the real CLI, apply mode, two synthetic single-cell libraries, the UMI floor adjudicated under
     the light floor, decisions accepting the floor ordering and "keep and removed are
     complementary": the run delivers, the floor is in STATUS.json's escapes and the other entry is
     not, and a REVIEW of 07_apply in the report refuses that entry by name
  3  one library's per-cell table edited as above, the deliverable deleted, the same decisions: exit
     2, no deliverable; the refusal names the check, refuses the entry by name, and prints no entry
     to paste for it
  4  the remedy is not blocked: the edited table removed, the run measures that library again and
     delivers, and the barcode the ledger lists as removed is not in the deliverable
"""
from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
fails: list[str] = []
FLOOR = "step 3: the light floor is strictly below the quality floor"
A = "keep and removed are complementary"
STANZA = f'check: "{A}"'


def module(rel: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m                   # a dataclass module must be registered to be built
    spec.loader.exec_module(m)
    return m


# ---- 1. the names, and integrity() on its own
ar = module("modules/07_apply/audit_removal.py", "w043f_audit_removal")
ACC, INT = tuple(getattr(ar, "ACCEPTABLE", ()) or ()), tuple(getattr(ar, "INTEGRITY", ()) or ())
# One table that fails every check at once (row 0 kept and removed with no criterion; row 1 kept and
# never scored, outside the cluster object; row 2 a doublet call nothing scored; row 3 an extra
# cluster-object cell that is not a doublet; row 4 an uncalled droplet kept), floors in the wrong order.
BAD = {"aligner_cell": [1, 1, 1, 1, 0, 0], "keep": [1, 1, 0, 0, 1, 0],
       "removed": [1, 0, 1, 1, 0, 1], "fail_doublet": [0, 0, 1, 0, 0, 0],
       "fail_umi_floor": [0, 0, 0, 1, 0, 0], "doublet_scored": [1, 0, 0, 1, 0, 0]}
PRE = [1, 0, 0, 1, 1, 0]
GOOD = {"aligner_cell": [1, 1, 1, 0], "keep": [1, 0, 0, 0], "removed": [0, 1, 1, 1],
        "fail_doublet": [0, 1, 0, 0], "fail_umi_floor": [0, 0, 1, 0], "doublet_scored": [1, 1, 1, 0]}
KW = dict(criteria=["fail_doublet", "fail_umi_floor"], cell_col="aligner_cell",
          scored_col="doublet_scored", doublet_criterion="fail_doublet")
ENTRY = {"step": "07_apply", "check": A, "approved_by": "the PI",
         "verbatim": "accepted for one library on an earlier run after reading its two rows by hand"}
bad = ar.audit(BAD, light_floor=200, quality_floor=150.0, predoublet_keep=PRE, **KW)
refusing = {f.check for f in bad if f.severity == "REFUSE"}
if ACC != (FLOOR,):
    fails.append(f"1: ACCEPTABLE is {ACC!r}, not the floor ordering alone")
if refusing != set(INT) | {FLOOR}:
    fails.append(f"1: the audit's refusals {sorted(refusing)} are not INTEGRITY plus the floor "
                 f"ordering ({sorted(INT)})")
if not hasattr(ar, "integrity"):
    fails.append("1: the audit has no integrity(): every refusal it gives goes to the gate, and an "
                 "accept entry lifts any of them by name")
else:
    try:
        ar.integrity(bad, [ENTRY])
        fails.append("1: integrity() did not raise on a table failing seven integrity checks")
    except ar.AuditFailure as e:
        why = str(e)
        if A not in why or "is refused" not in why or STANZA in why:
            fails.append(f"1: the integrity refusal does not refuse the entry by name, or prints an "
                         f"entry to paste: {why[:300]!r}")
    good = ar.audit(GOOD, light_floor=200, quality_floor=350.0, **KW)
    got = ar.integrity(good, [ENTRY, dict(ENTRY, check=FLOOR)])
    if [(f.check, f.severity) for f in got] != [("an accept entry names an integrity check", "REVIEW")] \
            or A not in got[0].message or "refused" not in got[0].message:
        fails.append(f"1: on a correct table an entry naming an integrity check is not reported once, "
                     f"as a REVIEW refusing it by name: {[(f.check, f.severity) for f in got]}")
print("1. the integrity checks are named, never lifted, and an entry naming one is reported")

# ---- 2 to 4. through the real CLI
try:
    import anndata
    import yaml
    from fixtures import cohorts as C
    HAVE = True
except ImportError as e:
    HAVE = False
    print(f"SKIP: the analysis stack is not importable here ({e})")

COLS = ["sample", "platform", "species", "reference", "assay", "matrix", "aligner_cells",
        "mt_prefix", "ribo_pattern", "dbr", "dbr_sd"]
QUALITY = {"umi_floor": {"value": 150, "approved_by": "the analyst",
                         "verbatim": "a floor of 150 UMI keeps every cell the aligner called"}}
ACCEPT = [{"step": "07_apply", "check": FLOOR, "approved_by": "the PI",
           "verbatim": "every called cell is above the light floor here; the floor is recorded on purpose"},
          ENTRY]


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


def run(td: Path, sheet: Path):
    proj, d = td / "proj", td / "decisions.yml"
    d.write_text(yaml.safe_dump({"quality": QUALITY, "accept": ACCEPT}))
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


def delivered(R) -> list:
    a = anndata.read_h5ad(R / "objects" / "cohort.deliverable.h5ad", backed="r")
    try:
        return list(a.obs_names)
    finally:
        a.file.close()


def true(r, c) -> bool:
    return str(r.get(c, "")).strip().lower() == "true"


if HAVE:
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        sheet = cohort(td)
        rc, log, R = run(td, sheet)
        obj = R / "objects" / "cohort.deliverable.h5ad" if R else None
        if rc != 0 or obj is None or not obj.exists():
            fails.append(f"2: the consistent run did not deliver (exit {rc}):\n"
                         + "\n".join(log.splitlines()[-8:]))
            R = None
        else:
            esc = [(e.get("ask") or {}).get("check") for e in status(R).get("escapes") or []]
            if FLOOR not in esc or A in esc:
                fails.append(f"2: STATUS.json's escapes are {esc}: the floor ordering must be accepted "
                             f"and the integrity check must not")
            gates = json.loads((R / "reports" / "payload.json").read_text()).get("gates") or []
            said = [g for g in gates if isinstance(g, dict) and g.get("step") == "07_apply"
                    and g.get("severity") == "REVIEW" and A in json.dumps(g) and "refused" in json.dumps(g)]
            if not said:
                fails.append("2: the report does not say the accept entry naming an integrity check "
                             "was refused")
        print("2. a consistent run delivers; only the floor ordering is accepted, the other entry is reported")

        if R is not None:
            t = R / "tables" / "LIB_B.percell.csv"
            with open(t, newline="") as fh:
                rows = list(csv.DictReader(fh))
            hdr = list(rows[0])
            r1 = next(r for r in rows if true(r, "removed") and true(r, "aligner_cell")
                      and true(r, "doublet_scored"))
            r2 = next(r for r in rows if true(r, "keep"))
            r1["keep"], r2["keep"] = "True", "False"
            with open(t, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=hdr)
                w.writeheader()
                w.writerows(rows)
            (R / "objects" / "cohort.deliverable.h5ad").unlink()
            rc, log, R = run(td, sheet)
            why = str((status(R).get("refusal") or {}).get("reason"))
            if rc != 2 or "07_apply" not in why or A not in why:
                fails.append(f"3: an edited table with the check accepted was not refused at 07_apply "
                             f"(exit {rc}; refusal {why[:160]!r})")
            if "is refused" not in why or STANZA in why:
                fails.append(f"3: the refusal does not refuse the accept entry by name, or prints an "
                             f"entry to paste for it: {why[:400]!r}")
            if R is not None and (R / "objects" / "cohort.deliverable.h5ad").exists():
                fails.append(f"3: the refused run wrote a deliverable that holds LIB_B_{r1['barcode']}, "
                             f"which its ledger lists as removed")
            print("3. a table edited to contradict itself is refused whatever the decisions file accepts; "
                  "nothing is delivered")

            t.unlink()
            rc, log, R = run(td, sheet)
            obj = R / "objects" / "cohort.deliverable.h5ad" if R else None
            if rc != 0 or obj is None or not obj.exists():
                fails.append(f"4: with the edited table removed the run did not deliver (exit {rc}):\n"
                             + "\n".join(log.splitlines()[-8:]))
            elif f"LIB_B_{r1['barcode']}" in delivered(R):
                fails.append(f"4: the deliverable holds LIB_B_{r1['barcode']}, which the ledger removes")
            print("4. the remedy is not blocked: the library is measured again and delivered as recorded")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("an integrity check of the removal audit is never accepted, and its remedy is not blocked")
