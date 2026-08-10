"""The mitochondrial ceiling: does it derive, and does it REFUSE when it should?

A derivation is only as good as what it declines to return. Every refusal below corresponds to a
way the calibration cohort could have been filtered wrongly, and three of them are things that
actually happened there before being caught:

  * a fence derived on an already-filtered set (measures the previous fence, not the data)
  * a bound that binds in most libraries (declared number masquerading as a derived one)
  * a removal that is uneven across the design (technical property read as biology)

Run: python tests/test_mito_ceiling.py
"""
from __future__ import annotations

import importlib.util
import pathlib
import random
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("quality", ROOT / "modules/05_quality/quality.py")
quality = importlib.util.module_from_spec(spec)
# Registered BEFORE exec: @dataclass resolves annotations via sys.modules[cls.__module__], which
# is None for a spec-loaded module that was never registered, and fails with a bare AttributeError.
sys.modules[spec.name] = quality
spec.loader.exec_module(quality)

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}" + (f"  --  {detail}" if detail and not cond
                                                       else ""))


def refuses(name, fn, must_mention=()):
    try:
        fn()
    except quality.ThresholdRefusal as e:
        missing = [m for m in must_mention if m.lower() not in str(e).lower()]
        check(name, not missing, f"message omits {missing}")
        return
    except Exception as e:  # noqa: BLE001
        check(name, False, f"raised {type(e).__name__}, not ThresholdRefusal: {e}")
        return
    check(name, False, "did not refuse")


def lognormal(rng, n, mu, sigma):
    """Mitochondrial percentages are right-skewed and bounded at 100."""
    return [min(100.0, rng.lognormvariate(mu, sigma)) for _ in range(n)]


rng = random.Random(0)

# --- the fence itself ---------------------------------------------------------------------------
print("\nfence()")
q1, q3, f = quality.fence([1, 2, 3, 4, 5, 6, 7, 8])
check("quartiles by linear interpolation", abs(q1 - 2.75) < 1e-9 and abs(q3 - 6.25) < 1e-9,
      f"got q1={q1} q3={q3}")
check("fence is Q3 + 1.5*IQR", abs(f - (6.25 + 1.5 * 3.5)) < 1e-9, f"got {f}")
refuses("fewer than 4 values refuses", lambda: quality.fence([1, 2, 3]), ("at least 4",))
check("non-finite values are dropped, not propagated",
      quality.fence([1, 2, 3, 4, 5, float("nan")])[2] == quality.fence([1, 2, 3, 4, 5])[2])

# --- derivation ----------------------------------------------------------------------------------
print("\nderive_mito_ceiling()")
tight = {f"lib{i}": lognormal(rng, 3000, 0.7, 0.8) for i in range(6)}
d = quality.derive_mito_ceiling(tight, assay="snrna")
cs = d["ceilings"]
check("one ceiling per library", len(cs) == len(tight))
check("every ceiling lies inside the bound",
      all(d["bounds"][0] <= m.ceiling <= d["bounds"][1] for m in cs.values()))
check("ceilings differ between libraries (not a cohort constant)",
      len({round(m.ceiling, 6) for m in cs.values()}) > 1)
check("clamped is labelled, not silent",
      all(m.clamped in ("", "lower", "upper") for m in cs.values()))
check("scrna gets a wider bound than snrna",
      quality.MITO_BOUNDS["scrna"][1] > quality.MITO_BOUNDS["snrna"][1])

refuses("unknown assay with no explicit bounds refuses",
        lambda: quality.derive_mito_ceiling(tight, assay="spatial"), ("unknown assay",))
refuses("a custom bound with no declared_by refuses",
        lambda: quality.derive_mito_ceiling(tight, bounds=(2.0, 60.0)), ("declared_by",))
check("a custom bound WITH declared_by is accepted",
      quality.derive_mito_ceiling(tight, bounds=(2.0, 60.0),
                                  declared_by="pilot tissue, cytoplasm retained")["bounds"]
      == (2.0, 60.0))
refuses("inverted bounds refuse",
        lambda: quality.derive_mito_ceiling(tight, bounds=(30.0, 5.0), declared_by="x"),
        ("lo < hi",))

# The bound must be a guard rail. Libraries with a huge spread push every fence past the bound;
# when that happens for most of them the "derived" ceiling is really the declared number.
wild = {f"lib{i}": lognormal(rng, 3000, 2.4, 1.6) for i in range(6)}
refuses("a bound that binds in most libraries refuses",
        lambda: quality.derive_mito_ceiling(wild, assay="snrna"),
        ("binds in", "guard rail"))

# --- the differential check -----------------------------------------------------------------------
print("\nassess_mito_removal()")
design = {"diet": {f"lib{i}": ("chow" if i % 2 else "hfd") for i in range(6)}}
a = quality.assess_mito_removal(tight, d, design=design)
check("reports an overall removal rate", 0.0 <= a["overall_rate"] <= 1.0)
check("reports a per-arm ratio", "diet" in a["arms"] and a["arms"]["diet"]["ratio"] is not None)
check("passes on an even cohort", a.get("worst_ratio", 0) < quality.DESIGN_REFUSE)

no_design = quality.assess_mito_removal(tight, d)
check("missing design is reported, not passed over",
      any("NO DESIGN" in n for n in no_design["notes"]))

# One arm mitochondria-heavy, the other clean: the ceilings follow each library, so the REMOVAL
# stays even. That is the property being claimed for a per-library rule, so it is tested.
split = {}
for i in range(6):
    split[f"lib{i}"] = lognormal(rng, 3000, 1.6 if i % 2 else 0.4, 0.8)
ds = quality.derive_mito_ceiling(split, assay="snrna")
asx = quality.assess_mito_removal(split, ds, design=design)
check("a per-library fence keeps removal even when arms genuinely differ",
      asx["worst_ratio"] < quality.DESIGN_REFUSE,
      f"ratio {asx.get('worst_ratio')}")

# A FIXED ceiling on the same data is the counter-example: it removes far more from the
# mitochondria-heavy arm. This is why the rule is per library.
fixed = {s: 12.0 for s in split}
try:
    quality.assess_mito_removal(split, fixed, design=design)
    check("a fixed ceiling on the same data is refused", False, "it passed")
except quality.ThresholdRefusal:
    check("a fixed ceiling on the same data is refused", True)

refuses("a library with no ceiling refuses rather than borrowing one",
        lambda: quality.assess_mito_removal(split, {k: v for k, v in list(ds["ceilings"].items())[:2]},
                                            design=design),
        ("no derived ceiling",))

# Exactly-on-the-line: 3.00x must refuse, not slip under by one ULP.
check("_at_least catches exactly-on-the-line", quality._at_least(0.30 / 0.10, 3.0))

# --- the note ------------------------------------------------------------------------------------
print("\nmito_ceiling_note()")
note = quality.mito_ceiling_note()
check("no longer claims the ceiling is underivable", "NOT DERIVABLE" not in note)
check("still records what remains adjudicated", "adjudicat" in note.lower())
check("names the bound as declared", "declared" in note.lower())

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
if FAIL:
    for f_ in FAIL:
        print(f"  FAILED: {f_}")
    sys.exit(1)
