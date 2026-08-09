# Plans ingestion for a set of inputs and prints each plan; it processes no data.
"""Step 0 test: every decision the ingest planner can make.

Five cases, chosen so each exercises one branch and three of them are real files with a known
right answer:

  1 vendor counts_matrix.tar passes P1, FAILS P2 (EmptyDrops-called) -> must REJECT, rebuild
  2 vendor .h5ad fails both -> must REJECT, rebuild
  3 a matrix reprocessed from FASTQ passes both -> must ACCEPT
  4 declared platform 'bgi' not implemented -> must BLOCK, not guess
  5 missing 'species' DECLARED field with no default -> must BLOCK
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "00_ingest"))
from ingest import plan_one, read_registry # noqa: E402

REG = read_registry(HERE.parent / "references" / "_registry" / "registry.tsv")
# An EXAMPLE species and reference, in the registry's `<species>/<build>` key form. Species is
# a DECLARED parameter with no default, so a value has to appear somewhere in this test; this
# one is illustrative and nothing in the planner is specific to it.
SPECIES = "mus_musculus"
REF = f"{SPECIES}/ensembl_112_filtered"

# Summary statistics for three real artifacts of one library, as recorded by the prior-filtering
# audit: the vendor's delivered matrix, the vendor's delivered analysis object, and the
# unfiltered droplet matrix reprocessed from FASTQ.
STATS = {
    "vendor_tar": dict(n_barcodes=17685, n_genes=34290, min_counts=500.0,
                       max_counts=21238.0, p98_counts=5421.0, expected_genes=34290),
    "vendor_h5ad": dict(n_barcodes=16564, n_genes=20000, min_counts=498.0,
                        max_counts=5420.0, p98_counts=5421.0, expected_genes=34290),
    "ours_raw": dict(n_barcodes=35915, n_genes=34290, min_counts=1,
                     max_counts=21238.0, p98_counts=None, expected_genes=34290),
}
CASES = [
    ("1 vendor counts_matrix.tar (CeleScope outs/filtered)", "vendor_tar", "run"),
    ("2 vendor .h5ad (delivered analysis object)", "vendor_h5ad", "run"),
    ("3 unfiltered droplet matrix reprocessed from FASTQ", "ours_raw", "accept"),
]

fails = []
print("Step 0 - ingest planning, library ctrl_01\n" + "=" * 74)

for label, key, expect in CASES:
    row = {"sample": "ctrl_01", "platform": "singleron", "species": SPECIES,
           "reference": REF, "matrix": __file__, "fastq_r1": "ctrl_01_R1.fq.gz"}
    p = plan_one(row, REG, stats_fn=lambda _p, k=key: STATS[k])
    p.sample = label
    print(p)
    if p.mode != expect:
        fails.append(f"{label}: expected {expect}, got {p.mode}")
    print()

print("-" * 74)
row = {"sample": "4 platform bgi (declared, not implemented)", "platform": "bgi",
       "species": SPECIES, "reference": REF, "fastq_r1": "x.fq.gz"}
p = plan_one(row, REG, stats_fn=lambda _p: STATS["ours_raw"])
print(p)
if p.mode != "blocked":
    fails.append("bgi: expected blocked")

print()
row = {"sample": "5 species omitted", "platform": "singleron", "species": "",
       "reference": REF, "fastq_r1": "x.fq.gz"}
p = plan_one(row, REG, stats_fn=lambda _p: STATS["ours_raw"])
print(p)
if p.mode != "blocked":
    fails.append("missing species: expected blocked")

print("\n" + "=" * 74)
if fails:
    print("FAILED:")
    for f in fails:
        print(" -", f)
    raise SystemExit(1)
print("proved: all 5 cases behave as specified - the two pre-filtered vendor artifacts are")
print("REJECTED and rebuilt, and an undeclared species blocks rather than defaulting")
