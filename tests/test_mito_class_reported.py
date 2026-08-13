# The report must state the class the run APPLIED, not the class the code was written expecting.
"""Does the report say how the mitochondrial ceiling was arrived at, truthfully?

WHY THIS TEST EXISTS

Two places named the mitochondrial ceiling's parameter class in fixed text - the parameter table
("DERIVED", "k derived from each library's Tukey fence") and the spine rule ("DERIVED · per
library"). Both were correct when written and both became false on 2026-08-13, when snRNA moved to
a DECLARED k and a bound narrow enough to dominate the result.

That is the worst failure this report can have and the hardest to see: the section whose entire
job is to say who set each number would have said "the data did", about a number an analyst
declared, on a page that otherwise renders perfectly. Nothing downstream could detect it, because
a wrong class reads exactly like a right one.

So the class is now READ FROM THE RUN in both places, and this checks that it travels: quality.py
decides it, step 5 records it, the parameter table reports it, and the spine takes it from there
rather than restating it.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "modules" / "05_quality"))

import quality  # noqa: E402
from report import build as B  # noqa: E402

fails: list[str] = []
print("Mitochondrial ceiling: the class the run applied is the class the report states")
print("=" * 78)


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'ok    ' if ok else 'FAILED'} {name}{('   ' + detail) if detail else ''}")
    if not ok:
        fails.append(f"{name}: {detail}")


def cohort(spread: float, n: int = 6, base: float = 2.0) -> dict:
    """n libraries whose mitochondrial values differ by `spread`, as a quartile summary."""
    out = {}
    for i in range(n):
        v = sorted([base * (1 + i * spread) * (1 + (j % 11) / 6) for j in range(300)])
        q = quality.fence(v)
        med = v[len(v) // 2]
        dev = sorted(abs(x - med) for x in v)
        out[f"lib{i}"] = {"n": len(v), "median": med, "q1": q[0], "q3": q[1],
                          "mad": dev[len(dev) // 2]}
    return out


# ---------------------------------------------------------------- quality.py decides it
tight = cohort(0.02)
snr = quality.derive_mito_ceiling_from_quartiles(tight, assay="snrna")
sc = quality.derive_mito_ceiling_from_quartiles(tight, assay="scrna")
check("snrna reports its k as declared by the assay", snr["k_source"] == "declared_assay",
      f"got {snr['k_source']!r}")
check("scrna reports its k as derived", sc["k_source"] == "derived", f"got {sc['k_source']!r}")
check("a caller's k outranks both", quality.derive_mito_ceiling_from_quartiles(
    tight, assay="snrna", k=6)["k_source"] == "declared_caller")

# A declared k must never silence the Tukey comparison, and the gap must be reported.
check("a declared k still carries the Tukey comparison",
      any("Tukey-implied" in n for n in snr["notes"]),
      "the declaration would otherwise be marking its own homework")

# A cohort where every library is flat: the derived route must REFUSE, and the declared route must
# proceed while saying it has no cross-check. Silence in the second case is the trap - it looks
# exactly like agreement.
flat = {f"lib{i}": {"n": 100, "median": 3.0, "q1": 3.0, "q3": 3.0, "mad": 0.0} for i in range(4)}
refused = False
try:
    quality.derive_mito_ceiling_from_quartiles(flat, assay="scrna")
except quality.ThresholdRefusal:
    refused = True
check("a flat cohort refuses where k must be derived", refused,
      "there is nothing to calibrate against")
flat_dec = quality.derive_mito_ceiling_from_quartiles(flat, assay="snrna")
check("...but proceeds where k is declared", flat_dec["k"] == quality.MITO_MAD_K["snrna"])
check("...and says it has no cross-check rather than staying silent",
      any("could NOT be computed" in n for n in flat_dec["notes"]),
      "no cross-check is not agreement")

# ---------------------------------------------------------------- and the report carries it
# PARAM_CLASSES is a CLOSED vocabulary of four. A row outside it is dropped by the validator, so
# a fifth class invented to carry "k was declared" would have removed the parameter from the table
# rather than mislabelling it in it - which is why the class stays one of the four and the detail
# lives in `basis`. Checked here so the next person to reach for a fifth finds this first.
check("the class vocabulary is exactly the four", set(B.PARAM_CLASSES) ==
      {"FIXED", "DERIVED", "DECLARED", "ADJUDICATED"}, f"{sorted(B.PARAM_CLASSES)}")

for klass, want_spine in (("DERIVED", "DERIVED · per library"),
                          ("DECLARED", "DECLARED · per library")):
    payload = {
        "run": {"project": "t", "mode": "apply"},
        "deliverable": {"n_in": 100, "n_kept": 90, "unit": "observations"},
        "parameters": [{"name": "mitochondrial ceiling", "value": "per library",
                        "class": klass, "basis": "step 5"}],
        "per_sample": {"source": "t.csv",
                       "columns": [{"key": "mito_ceiling_pct", "label": "mito",
                                    "scope": "per library"}],
                       "rows": [{"sample": "A", "mito_ceiling_pct": 7.5}]},
    }
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "r.html"
        B.build_report(payload, out, Path(td) / "r.json")
        page = out.read_text(encoding="utf-8")
    check(f"spine shows {klass!r}, taken from the parameter table", want_spine in page,
          "" if want_spine in page else f"{want_spine!r} absent; the spine restated a class")

# The one that would have shipped: a bound-dominated ceiling must NOT read as derived anywhere.
payload = {
    "run": {"project": "t", "mode": "apply"},
    "deliverable": {"n_in": 100, "n_kept": 90, "unit": "observations"},
    "parameters": [{"name": "mitochondrial ceiling", "value": "per library", "class": "DECLARED",
                    "basis": "step 5: ...; the bound decided the ceiling in most libraries, so "
                             "this is BOUND-DOMINATED and is reported as declared, not derived"}],
}
with tempfile.TemporaryDirectory() as td:
    out = Path(td) / "r.html"
    B.build_report(payload, out, Path(td) / "r.json")
    page = out.read_text(encoding="utf-8")
check("a bound-dominated ceiling never reads as DERIVED on the spine",
      "DERIVED · per library" not in page and "BOUND-DOMINATED" in page,
      "the spine must not label a declaration a measurement")

# A declared k keeps the class DERIVED - the applied number really is computed from each library's
# own distribution - but the BASIS must say the k was declared, or the report claims the whole
# derivation came from the data.
payload = {
    "run": {"project": "t", "mode": "apply"},
    "deliverable": {"n_in": 100, "n_kept": 90, "unit": "observations"},
    "parameters": [{"name": "mitochondrial ceiling", "value": "per library", "class": "DERIVED",
                    "basis": "step 5: median + k*1.4826*MAD over the barcodes above the light "
                             "floor and below 50% mitochondrial, with k = 3 DECLARED for assay "
                             "snrna, and the whole bounded to 5.0-10.0%"}],
}
with tempfile.TemporaryDirectory() as td:
    out = Path(td) / "r.html"
    B.build_report(payload, out, Path(td) / "r.json")
    page = out.read_text(encoding="utf-8")
check("a declared k stays DERIVED but says so in the basis",
      "DERIVED · per library" in page and "DECLARED for assay" in page
      and "below 50% mitochondrial" in page,
      "the class is one of four; the detail belongs in the basis")

print("=" * 78)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("class reporting OK - quality decides it, the report states it, nothing restates it")
