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
check("clamped is labelled, not silent",
      all(m.clamped in ("", "lower", "upper") for m in cs.values()))

# A CLEAN cohort gets the declared floor, and that is the correct answer, not a failure. `tight`
# is a low-mitochondrial cohort - what a good nuclei prep looks like - so every fence lands below
# the 10% snRNA floor and the ceiling becomes one number. The pipeline must SAY that rather than
# present a flat constant as if each library had derived it.
check("a clean cohort lands on the floor, and is not called derived",
      d["provenance"] == "bound_dominated", f"got {d['provenance']!r}")
check("...and the flat ceiling is the floor itself",
      {round(m.ceiling, 6) for m in cs.values()} == {float(d["bounds"][0])})

# Per-library variation is a property of the ESTIMATOR, so it is tested where the bound is not
# deciding: with a wide declared bound the same cohort must produce ten different fences.
free = quality.derive_mito_ceiling(tight, bounds=(0.5, 60.0),
                                   declared_by="probe: isolate the estimator from the bound")
check("unbounded, ceilings differ between libraries (not a cohort constant)",
      len({round(m.ceiling, 6) for m in free["ceilings"].values()}) > 1)
check("unbounded, the derivation is reported as derived",
      free["provenance"] == "derived", f"got {free['provenance']!r}")

# The applied fence is the MAD route; Tukey is carried as an independent second derivation and is
# never applied. Both must be present on every library, and they must actually differ in general -
# a cross-check that always equals the applied number is not a check.
m0 = next(iter(free["ceilings"].values()))
check("both fences are recorded per library",
      m0.derived > 0 and m0.tukey > 0 and m0.mad >= 0)
check("the APPLIED fence is the MAD route, not Tukey",
      abs(m0.derived - (m0.median + free["k"] * quality.MAD_SCALE * m0.mad)) < 1e-9,
      f"derived={m0.derived} median={m0.median} mad={m0.mad} k={free['k']}")
check("the cross-check is Tukey and is a different number",
      abs(m0.tukey - (m0.q3 + quality.IQR_MULT * m0.iqr)) < 1e-9
      and m0.tukey != m0.derived)
check("skew ratio is reported (1.349 would be normal)", m0.skew_ratio > 0)
# k is DERIVED from the Tukey fence, not a module constant. That is the whole point: Tukey
# calibrates how far out the fence sits, MAD applies it, and neither is a number someone typed.
check("k is derived from the cohort, not declared", free.get("k_source") == "derived")
check("...and the selection is reported, not just its result",
      free["k_selection"]["per_library"] and free["k_selection"]["median_k"] > 0)
check("...an explicitly passed k is recorded as declared instead",
      quality.derive_mito_ceiling(tight, bounds=(0.5, 60.0), declared_by="probe",
                                  k=7)["k_source"] == "declared")
check("a cohort whose libraries share a shape derives a tight implied k",
      free["k_selection"]["spread"] < 2.0, f"spread {free['k_selection']['spread']}")
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
#
# Until 2026-08-11 that RAISED. It now RECLASSIFIES: refusing did not stop the number being
# wrongly classified, it stopped the run, and it fired hardest when an analyst had deliberately
# narrowed the bound. The class has to be honest; the run does not have to die for it to be.
wild = {f"lib{i}": lognormal(rng, 3000, 2.4, 1.6) for i in range(6)}
w = quality.derive_mito_ceiling(wild, assay="snrna")
check("a bound that binds in most libraries does NOT refuse", w["ceilings"] != {})
check("...it is reclassified as bound_dominated, not derived",
      w["provenance"] == "bound_dominated", f"got {w['provenance']!r}")
check("...and says so where a reader will meet it",
      any("BOUND-DOMINATED" in n for n in w["notes"]), f"notes={w['notes']}")
check("a cohort the bound only rails stays derived", free["provenance"] == "derived")

# The two clamp directions are opposite events and must not share a counter: an upper clamp
# removes nuclei the library's own fence would have kept, a lower clamp retains nuclei it would
# have cut. Only the first can delete signal.
check("clamp directions are reported separately",
      set(w["clamped_lower"]) | set(w["clamped_upper"])
      == {s for s, m in w["ceilings"].items() if m.clamped},
      f"lower={w['clamped_lower']} upper={w['clamped_upper']}")
check("a library is never counted in both directions",
      not (set(w["clamped_lower"]) & set(w["clamped_upper"])))

# Spread of the APPLIED ceiling - not implied by the design differential, and the property that
# went unmeasured while a 6.12% library and a 25.00% library sat in one deliverable.
check("ceiling spread is reported", isinstance(w["ceiling_spread"], float))
check("ceiling spread is max/min of the applied ceilings",
      abs(w["ceiling_spread"]
          - max(m.ceiling for m in w["ceilings"].values())
          / min(m.ceiling for m in w["ceilings"].values())) < 1e-9)
check("a spread above the review line is called out",
      (w["ceiling_spread"] <= quality.CEILING_SPREAD_REVIEW)
      or any("differs" in n and "least and most permissive" in n for n in w["notes"]),
      f"spread={w['ceiling_spread']} notes={w['notes']}")

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

# One arm mitochondria-heavy, the other clean. A per-library fence follows each library, so the
# REMOVAL stays even - that is the property being claimed for a per-library rule, so it is tested
# where the estimator is actually free to express it.
split = {}
for i in range(6):
    split[f"lib{i}"] = lognormal(rng, 3000, 1.6 if i % 2 else 0.4, 0.8)
ds_free = quality.derive_mito_ceiling(split, bounds=(0.5, 60.0),
                                      declared_by="probe: the per-library property, unbounded")
asx = quality.assess_mito_removal(split, ds_free, design=design)
check("an unbounded per-library fence keeps removal even when arms genuinely differ",
      asx["worst_ratio"] < quality.DESIGN_REFUSE,
      f"ratio {asx.get('worst_ratio')}")

# THE COST OF THE FLOOR, asserted rather than discovered later.
#
# With the snRNA floor at 10%, the CLEAN arm can no longer track its own distribution - its fences
# are pushed up to 10 while the heavy arm keeps a fence of its own. Removal then diverges across
# the arms, and on a cohort split this hard it diverges past the refusal line. That is a real
# consequence of choosing a floor for biological plausibility, and the pipeline refuses rather
# than applying it quietly. The test exists so that the day someone widens or narrows the floor,
# this is what tells them what they changed.
ds = quality.derive_mito_ceiling(split, assay="snrna")
try:
    quality.assess_mito_removal(split, ds, design=design)
    check("the 10% floor on a clean-vs-heavy cohort is refused, not applied quietly",
          False, "it passed - the floor no longer costs anything, which needs explaining")
except quality.ThresholdRefusal as e:
    check("the 10% floor on a clean-vs-heavy cohort is refused, not applied quietly",
          "removes unevenly" in str(e), str(e)[:90])

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
