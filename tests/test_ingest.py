# Plans ingestion for a set of inputs and prints each plan; it processes no data.
"""Step 0 test: every decision the ingest planner can make.

Five cases, chosen so each exercises one branch and three of them are real files with a known
right answer:

  1 vendor counts_matrix.tar passes P1, FAILS P2 (EmptyDrops-called) -> must REJECT, rebuild
  2 vendor .h5ad fails both -> must REJECT, rebuild
  3 a matrix reprocessed from FASTQ passes both -> must ACCEPT
  4 declared platform 'bgi' not implemented -> must BLOCK, not guess
  5 missing 'species' DECLARED field with no default -> must BLOCK
  6 missing 'assay' -> must BLOCK (it had a silent default of single-nucleus)
  7 an assay with no policy -> must BLOCK, not be read as either
  8 a CELL-CALLED matrix (min UMI 500) on SINGLE CELLS -> ACCEPTED and named cell-called; the same
    row on single NUCLEI -> BLOCKED without FASTQ, rebuilt with it (harness autodev W020, B1-046:
    P2 is required where an empty-droplet model consumes it, and none runs on cells); a P1 failure
    is rejected on both; verify() takes no default assay
  9 a platform with no aligner here ('bd') and a VERIFIED supplied matrix -> ACCEPTED, recorded
    'processor: none (supplied matrix)'; the same platform with FASTQ only, or with a rejected
    matrix to rebuild -> BLOCKED, naming the missing aligner (harness autodev W020, B1-048)
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "00_ingest"))
from ingest import plan_one, read_registry # noqa: E402

# An EXAMPLE species and reference, in the registry's `<species>/<build>` key form. Species is
# a DECLARED parameter with no default, so a value has to appear somewhere in this test; this
# one is illustrative and nothing in the planner is specific to it.
SPECIES = "mus_musculus"
REF = f"{SPECIES}/ensembl_112_filtered"
# The test's OWN registry, written as a project writes one. scQC ships none: the registry it
# used to ship held one cohort's reference, read for any project that brought none.
import tempfile  # noqa: E402

with tempfile.TemporaryDirectory() as _td:
    _reg_file = Path(_td) / "registry.tsv"
    _reg_file.write_text(f"species\tbuild\tpath\n{SPECIES}\tensembl_112_filtered\t-\n",
                         encoding="utf-8")
    REG = read_registry(_reg_file)

# Summary statistics in the shape of three artifacts of one library, as a prior-filtering audit
# records them (synthetic values): the vendor's delivered matrix, the vendor's delivered analysis
# object, and the unfiltered droplet matrix reprocessed from FASTQ.
STATS = {
    "vendor_tar": dict(n_barcodes=15000, n_genes=30000, min_counts=500.0,
                       max_counts=20000.0, p98_counts=5000.0, expected_genes=30000),
    "vendor_h5ad": dict(n_barcodes=14000, n_genes=20000, min_counts=470.0,
                        max_counts=4999.0, p98_counts=5000.0, expected_genes=30000),
    "ours_raw": dict(n_barcodes=36000, n_genes=30000, min_counts=1,
                     max_counts=20000.0, p98_counts=None, expected_genes=30000),
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
           "reference": REF, "assay": "snrna", "matrix": __file__, "fastq_r1": "ctrl_01_R1.fq.gz"}
    p = plan_one(row, REG, stats_fn=lambda _p, k=key: STATS[k])
    p.sample = label
    print(p)
    if p.mode != expect:
        fails.append(f"{label}: expected {expect}, got {p.mode}")
    print()

print("-" * 74)
row = {"sample": "4 platform bgi (declared, not implemented)", "platform": "bgi",
       "species": SPECIES, "reference": REF, "assay": "snrna", "fastq_r1": "x.fq.gz"}
p = plan_one(row, REG, stats_fn=lambda _p: STATS["ours_raw"])
print(p)
if p.mode != "blocked":
    fails.append("bgi: expected blocked")

print()
row = {"sample": "5 species omitted", "platform": "singleron", "species": "",
       "reference": REF, "assay": "snrna", "fastq_r1": "x.fq.gz"}
p = plan_one(row, REG, stats_fn=lambda _p: STATS["ours_raw"])
print(p)
if p.mode != "blocked":
    fails.append("missing species: expected blocked")

# 6 and 7. THE ASSAY HAD A DEFAULT, AND IT WAS THE FIRST COHORT'S. It was not required here and
# the mitochondrial step read it as `params.get("assay", "snrna")`, so a single-cell cohort was
# bounded as nuclei and nothing said so. Absent or unknown, it now blocks at step 0.
print()
for label, assay, says in (("6 assay omitted", "", "'assay' is required"),
                           ("7 assay with no policy", "spatial", "assay 'spatial' unknown")):
    row = {"sample": label, "platform": "10x", "species": SPECIES, "reference": REF,
           "assay": assay, "fastq_r1": "x.fq.gz"}
    p = plan_one(row, REG, stats_fn=lambda _p: STATS["ours_raw"])
    print(p)
    if p.mode != "blocked" or says not in p.reason:
        fails.append(f"{label}: expected blocked saying {says!r}, got {p.mode}: {p.reason}")

# 8. A CELL-CALLED MATRIX ON SINGLE CELLS (harness autodev W020, B1-046). The register's CHECK:
# a scrna row with min_counts 500 is accepted; the same snrna row is blocked.
print()
for assay, fastq, expect in (("scrna", "", "accept"), ("snrna", "", "blocked"),
                             ("scrna", "x_R1.fq.gz", "accept"), ("snrna", "x_R1.fq.gz", "run")):
    row = {"sample": f"8 {assay} cell-called{' +FASTQ' if fastq else ''}", "platform": "10x",
           "species": SPECIES, "reference": REF, "assay": assay, "matrix": __file__,
           "fastq_r1": fastq}
    p = plan_one(row, REG, stats_fn=lambda _p: STATS["vendor_tar"])
    print(p)
    if p.mode != expect:
        fails.append(f"8 {assay}{' +FASTQ' if fastq else ''}: expected {expect}, got {p.mode}: "
                     f"{p.reason}")
    if assay == "scrna":
        kind = getattr(p, "input_kind", None)
        if kind != "called_cells" or "CELL-CALLED" not in p.reason:
            fails.append(f"8 {assay}: accepted without being named a cell-called input "
                         f"(input_kind {kind!r}): {p.reason}")
        if p.verdict is None or p.verdict.p2_raw_droplets:
            fails.append(f"8 {assay}: the P2 failure was not kept in the verdict")
    print()
# P1 is required on every assay: the vendor's capped, gene-subset object is rejected for cells too.
row = {"sample": "8 scrna capped object", "platform": "10x", "species": SPECIES, "reference": REF,
       "assay": "scrna", "matrix": __file__}
p = plan_one(row, REG, stats_fn=lambda _p: STATS["vendor_h5ad"])
print(p)
if p.mode != "blocked":
    fails.append(f"8 scrna with a P1 failure: expected blocked, got {p.mode}")
sys.path.insert(0, str(HERE.parent / "lib"))
import verify_raw as _vr  # noqa: E402
try:
    _vr.verify(name="x", **STATS["vendor_tar"])
    fails.append("8: verify() ran with no assay - a default assay would lend one assay's rule to the other")
except TypeError:
    pass
for bad in ("", "spatial"):
    try:
        _vr.verify(name="x", assay=bad, **STATS["vendor_tar"])
        fails.append(f"8: verify() accepted assay {bad!r}")
    except ValueError:
        pass

# 9. THE PLATFORM SELECTS AN ALIGNER, WHICH ONLY THE FASTQ ROUTE RUNS (harness autodev W020, B1-048).
print()
for label, matrix, fastq, key, expect in (
        ("9 bd + verified matrix", __file__, "", "ours_raw", "accept"),
        ("9 bd + FASTQ only", "", "x_R1.fq.gz", "ours_raw", "blocked"),
        ("9 bd + rejected matrix + FASTQ", __file__, "x_R1.fq.gz", "vendor_tar", "blocked")):
    row = {"sample": label, "platform": "bd", "species": SPECIES, "reference": REF,
           "assay": "snrna", "matrix": matrix, "fastq_r1": fastq}
    p = plan_one(row, REG, stats_fn=lambda _p, k=key: STATS[k])
    print(p)
    if p.mode != expect:
        fails.append(f"{label}: expected {expect}, got {p.mode}: {p.reason}")
    if expect == "accept" and p.processor != "none (supplied matrix)":
        fails.append(f"{label}: the record names processor {p.processor!r} for a matrix nothing "
                     f"aligned")
    if expect == "blocked" and "no aligner" not in p.reason:
        fails.append(f"{label}: blocked without naming the missing aligner: {p.reason}")
    print()

# THE THREE LISTS OF ASSAYS ARE ONE LIST. Ingest admits an assay, ambient chooses a policy for it
# and quality bounds it; an assay admitted by one and unknown to another refuses three steps late.
import ast  # noqa: E402
import ingest as _ing  # noqa: E402

def _literal(rel, name):
    for node in ast.walk(ast.parse((HERE.parent / rel).read_text(encoding="utf-8"))):
        if isinstance(node, ast.Assign) and any(getattr(x, "id", "") == name for x in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found in {rel}")

_lists = {"ingest.ASSAYS": set(_ing.ASSAYS),
          "verify_raw.P2_REQUIRED": set(_vr.P2_REQUIRED),
          "ambient.ASSAYS": set(_literal("modules/01_ambient/ambient.py", "ASSAYS")),
          "quality.MITO_BOUNDS": set(_literal("modules/05_quality/quality.py", "MITO_BOUNDS"))}
if len({frozenset(v) for v in _lists.values()}) != 1:
    fails.append(f"the assays disagree between modules: {_lists}")

print("\n" + "=" * 74)
if fails:
    print("FAILED:")
    for f in fails:
        print(" -", f)
    raise SystemExit(1)
print("proved: all 9 cases behave as specified - the two pre-filtered vendor artifacts are")
print("REJECTED and rebuilt for nuclei, a cell-called matrix is ACCEPTED and named for cells, and an")
print("undeclared species or assay blocks rather than defaulting")
