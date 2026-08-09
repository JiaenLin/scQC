# Plans a doublet sweep, recommends a setting and checks call health, then prints; flags only.
"""Step 4 test: a real dbr.sd sweep and a real detector failure, against one tuned detector."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "04_doublets"))
from doublet import DoubletRefusal, SweepResult, check_detector, recommend # noqa: E402
from doublet_health import health, verdict # noqa: E402

SAMPLES = ["treat_03", "treat_05", "ctrl_02", "treat_04", "ctrl_05",
           "treat_02", "treat_01", "ctrl_04", "ctrl_01", "ctrl_03"]
# dbr = 0.06 throughout; only dbr.sd moves between the three settings
SWEEP_RATES = {
    "default": [.1004, .1031, .1032, .1022, .1062, .1030, .1038, .1074, .1079, .1068],
    "dbr": [.1185, .1207, .1186, .1167, .1234, .1174, .1244, .1245, .1281, .1223],
    "1": [.1281, .1562, .1929, .2965, .1567, .1807, .1825, .1766, .1934, .3169],
}
DEEP = {"default": 0.328, "dbr": 0.369, "1": 0.637}
BAND = {"default": 2, "dbr": 9, "1": 0}
CONDITION = {s: ("treat" if "treat" in s else "ctrl") for s in SAMPLES}

class ScDblFinder:
    name = "scDblFinder"
    reproducible = True
    needs_empty_drops_removed = True
    min_umi_floor = 200
    imports_rate_prior = True

    def score(self, matrix, sample, seed, **kw):
        raise NotImplementedError

class Sloppy:
    name = "MyDetector"
    reproducible = False
    needs_empty_drops_removed = True
    min_umi_floor = None
    imports_rate_prior = True

fails = []
print("Step 4 - one detector, tuned\n" + "=" * 74)

print("A. contract: dbr not declared, and the detector imports a prior")
try:
    check_detector(ScDblFinder(), dbr=None, light_floor=200)
    fails.append("A: must refuse without a declared dbr")
except DoubletRefusal as e:
    print(f"[REFUSED] {str(e)[:165]}...")

print("\nB. contract satisfied")
for n in check_detector(ScDblFinder(), dbr=0.06, light_floor=200):
    print(f" note: {n}")
print(" ok")

print("\nC. the dbr.sd sweep")
res = []
for k, v in SWEEP_RATES.items():
    r = SweepResult(k, dict(zip(SAMPLES, v)), deep_decile_rate=DEEP[k],
                    in_published_band=BAND[k], n_samples_in_band=10)
    m = {}
    for s, rate in r.per_sample_rate.items():
        m.setdefault(CONDITION[s], []).append(rate)
    mm = {a: sum(x) / len(x) for a, x in m.items()}
    r.worst_arm_ratio = max(mm.values()) / min(mm.values())
    r.worst_arm_factor = "condition"
    res.append(r)
    print(" " + str(r))
rec = recommend(res)
print("\n" + str(rec))
if rec.setting != "dbr":
    fails.append(f"C: expected dbr.sd=dbr, got {rec.setting}")

print("\n" + "-" * 74)
print("D. health of the chosen calls (dbr.sd = dbr)")
f = health(res[1].per_sample_rate, {"condition": CONDITION}, deep_decile_rate=DEEP["dbr"],
           n_kept_unscored=0, detector_name="scDblFinder", reproducible=True)
for x in f:
    print(x)
print(f"\n verdict: {verdict(f)}")
if verdict(f) == "REFUSE":
    fails.append("D: the accepted calls must not be refused")

print("\n" + "-" * 74)
print("E. Scrublet's real failure, replayed through the health check")
SCRUB = {"ctrl_03": .2064, "ctrl_04": .0729, "ctrl_05": .0000, "treat_03": .0000,
         "treat_04": .0001, "treat_05": .0000, "ctrl_01": .1139, "ctrl_02": .0001,
         "treat_01": .0000, "treat_02": .0000}
f2 = health(SCRUB, {"condition": CONDITION}, detector_name="Scrublet", reproducible=True)
for x in f2:
    if x.severity != "ok":
        print(x)
if verdict(f2) != "REFUSE":
    fails.append("E: Scrublet's failure must REFUSE")

print("\nF. a supplied detector that declares itself non-reproducible")
for n in check_detector(Sloppy(), dbr=0.06, light_floor=200):
    print(f" note: {n}")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: C reaches the chosen dbr.sd from the measurements rather than from the")
print("conclusion, and E shows one detector plus a health check catching the failure that")
print("a four-tool consensus is usually assembled to catch")
