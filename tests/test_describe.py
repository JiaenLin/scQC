#!/usr/bin/env python3
"""`scqc describe`: the tool's declaration, and its agreement with the code it describes.

WHAT IS CHECKED

  1. `scqc describe` prints JSON with the contract's fields: needs, provides, sees, cannot_show,
     state_version, gates, verdicts, escapes, criteria
  2. the criteria it names are exactly APPLY_CRITERIA in adapters/scanpy_ops.py, and every one
     has a mask in `provides`
  3. every gate module it names exists and defines the three-verdict vocabulary
  4. cannot_show has at least three sentences; cannot_establish carries every STEP_TEXT step
  5. state_version is an integer and sees is a list (an empty list would be a positive claim)

Stdlib only.
"""
from __future__ import annotations

import io
import json
import re
import sys
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f"  {detail}" if detail and not ok else ""))
    if not ok:
        FAILS.append(name)


def main():
    import scqc_cli
    out = io.StringIO()
    with redirect_stdout(out):
        rc = scqc_cli.main(["describe"])
    d = json.loads(out.getvalue())
    check("1 describe exits 0 and prints JSON", rc == 0 and isinstance(d, dict))
    for f in ("needs", "provides", "sees", "cannot_show", "state_version", "gates", "verdicts",
              "escapes", "criteria", "name", "version", "commit"):
        check(f"1 field {f}", f in d)
    src = (ROOT / "adapters" / "scanpy_ops.py").read_text(encoding="utf-8", errors="replace")
    m = re.search(r"APPLY_CRITERIA\s*=\s*\((.*?)\)", src, re.S)
    crit = re.findall(r"\"([a-z_]+)\"", m.group(1))
    check("2 criteria are APPLY_CRITERIA", d["criteria"] == crit, f"{d['criteria']} vs {crit}")
    check("2 one mask per criterion in provides", all(f"mask/{c}" in d["provides"] for c in crit))
    for name, g in d["gates"].items():
        mod = ROOT / g["module"]
        text = mod.read_text(encoding="utf-8", errors="replace") if mod.exists() else ""
        check(f"3 gate {name}: module exists, verdicts within the vocabulary, words as written are in it",
              mod.exists() and set(g["verdicts"]) <= {"PASS", "REVIEW", "REFUSE"}
              and all(w in text for w in g["as_written"]), f"{g}")
    check("3 verdict vocabulary is closed", d["verdicts"] == ["PASS", "REVIEW", "REFUSE"])
    check("4 cannot_show has at least three sentences", len(d["cannot_show"]) >= 3)
    from engine import steps
    check("4 cannot_establish covers every step", set(d["cannot_establish"]) == set(steps.STEP_TEXT))
    check("5 state_version is an integer", isinstance(d["state_version"], int) and not isinstance(d["state_version"], bool))
    check("5 sees is a list", isinstance(d["sees"], list))
    check("5 version agrees with VERSION", d["version"] == (ROOT / "VERSION").read_text().strip())
    print(f"\n{'PASS' if not FAILS else 'FAIL'}: {len(FAILS)} failing")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
