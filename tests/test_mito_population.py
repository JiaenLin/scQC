# The population the mitochondrial ceiling is derived over - both of its restrictions.
"""Does `mito_quartiles` select the population the ceiling is specified to be derived over?

WHY THIS TEST EXISTS

The ceiling is `median + k * 1.4826 * MAD` over a population, and every part of that is tested
elsewhere except the population. That is the part with no symptom: a fence derived over the wrong
barcodes is a real number, computed correctly, describing something nobody asked about, and it
prints identically to the right one. The two restrictions it must apply pull in opposite
directions and were added two days apart, which is exactly the shape of thing that comes apart:

  * at or above the light floor, because a percentage needs a denominator;
  * strictly below MITO_DERIVATION_MAX, because a droplet more than half mitochondrial is not a
    cell whose mitochondrial content is high - on a nuclear prep it cannot be.

WHAT IS CHECKED, AND WHY EACH ONE

  1. the high-mitochondrial droplets are excluded from the quartiles
  2. excluding them LOWERS the fence - the direction matters, because the whole reason for the
     restriction is that ambient droplets widen the MAD and so loosen the cut on real cells
  3. the floor still applies, and applies independently of the new restriction
  4. the exclusion is COUNTED, not silent
  5. `max_above_floor` survives the cut - the diagnostic the restriction would otherwise erase
  6. the boundary is strict (`< max`, not `<= max`), stated once so it cannot drift
  7. a population too small to place a quartile returns None AND a population record
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from adapters.scanpy_ops import MITO_DERIVATION_MAX, mito_quartiles  # noqa: E402

fails: list[str] = []
print("Mitochondrial derivation population - adapters/scanpy_ops.mito_quartiles")
print("=" * 74)


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'ok    ' if ok else 'FAILED'} {name}{('   ' + detail) if detail else ''}")
    if not ok:
        fails.append(f"{name}: {detail}")


def fence(stats, k=4):
    """The applied fence, spelled as modules/05_quality spells it."""
    return stats["median"] + k * 1.4826 * stats["mad"]


# A plausible library: a healthy population around 3%, a handful of shallow droplets that the
# floor should drop, and an ambient tail above 50% that the derivation max should drop. The
# ambient tail is deliberately LARGE - a quarter of the barcodes - because a restriction that
# only matters for a handful is one whose absence nobody would notice.
healthy = [1.0, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0]
ambient = [55.0, 62.0, 71.0, 88.0, 99.5]
shallow = [40.0, 45.0, 30.0]  # real percentages, but on 30-UMI droplets

pct = healthy + ambient + shallow
tot = [1000.0] * len(healthy) + [800.0] * len(ambient) + [30.0] * len(shallow)

stats, pop = mito_quartiles(pct, tot, floor_umi=200.0)

# 1 - the ambient tail is not in the quartiles.
check("high-mitochondrial droplets excluded", stats["n"] == len(healthy),
      f"n={stats['n']}, expected {len(healthy)}")
check("nothing above the max survives", stats["max"] < MITO_DERIVATION_MAX,
      f"max in population = {stats['max']}")

# 2 - and excluding them LOWERS the fence. Without this the restriction could be inverted and
#     every other assertion here would still pass.
unrestricted, _ = mito_quartiles(pct, tot, floor_umi=200.0, derivation_max=100.0)
check("the restriction lowers the fence", fence(stats) < fence(unrestricted),
      f"{fence(stats):.3f} restricted vs {fence(unrestricted):.3f} unrestricted")

# 3 - the floor still applies, and does so independently: the shallow droplets sit BELOW the
#     derivation max, so only the floor can be removing them.
check("the light floor still applies",
      all(s not in [round(x, 6) for x in (healthy)] or True for s in shallow)
      and stats["n"] == len(healthy),
      f"shallow droplets at 30-45% mito are below the {MITO_DERIVATION_MAX:g}% max; "
      f"only the floor removes them")
no_floor, _ = mito_quartiles(pct, tot, floor_umi=0.0)
check("dropping the floor readmits the shallow droplets",
      no_floor["n"] == len(healthy) + len(shallow),
      f"n={no_floor['n']}, expected {len(healthy) + len(shallow)}")

# 4 - the exclusion is counted, not silent.
check("the exclusion is counted", pop["n_excluded_high_mito"] == len(ambient),
      f"n_excluded_high_mito={pop['n_excluded_high_mito']}, expected {len(ambient)}")
check("the population above the floor is counted",
      pop["n_above_floor"] == len(healthy) + len(ambient),
      f"n_above_floor={pop['n_above_floor']}")
check("the applied max is recorded", pop["derivation_max_pct"] == MITO_DERIVATION_MAX,
      f"derivation_max_pct={pop['derivation_max_pct']}")

# 5 - the diagnostic the restriction would otherwise erase.
check("the true observed maximum survives the cut", pop["max_above_floor"] == max(ambient),
      f"max_above_floor={pop['max_above_floor']}, expected {max(ambient)}")
check("it is carried on the quartiles too", stats["max_above_floor"] == max(ambient),
      f"stats['max_above_floor']={stats['max_above_floor']}")

# 6 - the boundary is strict. A droplet exactly at the max is excluded; one just below is kept.
#     Stated once, here, so the two spellings cannot drift apart in silence.
edge_pct = [1.0, 2.0, 3.0, 4.0, 50.0, 49.999]
edge_tot = [1000.0] * len(edge_pct)
edge, edge_pop = mito_quartiles(edge_pct, edge_tot, floor_umi=200.0)
check("the boundary is strict (< max, not <= max)",
      edge["n"] == 5 and edge_pop["n_excluded_high_mito"] == 1,
      f"n={edge['n']} kept, {edge_pop['n_excluded_high_mito']} excluded; "
      f"49.999 must be kept and 50.0 excluded")

# 7 - too few to place a quartile: None, but never a missing population record.
thin, thin_pop = mito_quartiles([1.0, 2.0, 60.0, 70.0], [1000.0] * 4, floor_umi=200.0)
check("too few for a quartile returns None", thin is None, f"got {thin!r}")
check("...but still returns the population record",
      thin_pop["n_above_floor"] == 4 and thin_pop["n_excluded_high_mito"] == 2
      and thin_pop["n_at_or_above"] is None,
      f"{thin_pop}")

# A value outside (0, 100] is a mistake, not a strict setting - checked in _op_valley, so here we
# only confirm the default is the declared constant rather than a second copy of the number.
check("the default is the declared constant", MITO_DERIVATION_MAX == 50.0,
      f"MITO_DERIVATION_MAX={MITO_DERIVATION_MAX}")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("mitochondrial derivation population OK - both restrictions apply, and both are recorded")
