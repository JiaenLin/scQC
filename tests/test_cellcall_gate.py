# Gates recorded aligner and denoiser cell counts and prints the findings; it calls no cells.
"""Step 2 test: a constructed set of cell calls, plus the faults the gate exists to stop.

Every count here is CONSTRUCTED from a unit and the lines the gate documents, and every expected
result is COMPUTED below from those lines, so no number stands for any cohort's: REVIEW when at
least REVIEW_AT of the aligner's cells are lost (the line inclusive), REFUSE when more than
REFUSE_ABOVE are, and "stricter" when the denoiser calls fewer cells than the aligner, refused
when it also loses at least REVIEW_AT. Each line is probed on both sides of it.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "modules" / "02_cells"))
from cellcall_gate import gate, verdict # noqa: E402

# THE LINES THIS TEST HOLDS THE GATE TO, stated here rather than read from the module, so a line
# that moves fails the test instead of moving the test with it.
REVIEW_AT, REFUSE_ABOVE = 0.05, 0.10

# CONSTRUCTED calls. Each aligner count is a multiple of UNIT squared, so a loss ON a line is a
# whole number of cells; the denoiser calls half again as many as the aligner everywhere.
UNIT = 40
LIBS = ["ctrl_01", "ctrl_02", "ctrl_03", "ctrl_04", "treat_01", "treat_02", "treat_03", "treat_04"]
COHORT = {s: dict(aligner=UNIT ** 2 * (i + 3), lost=0) for i, s in enumerate(LIBS)}
for v in COHORT.values():
    v["cellbender"] = v["aligner"] * 3 // 2


def on_line(lib, line, calls):
    """The number of the aligner's cells that is exactly `line` of them."""
    return round(line * calls[lib]["aligner"])


# Every loss on ONE arm: treat_01 exactly ON the review line, treat_02 one cell under it, treat_03
# a handful. The second factor splits the two largest losses, so a loss balanced on it is seen.
COHORT["treat_01"]["lost"] = on_line("treat_01", REVIEW_AT, COHORT)
COHORT["treat_02"]["lost"] = on_line("treat_02", REVIEW_AT, COHORT) - 1
COHORT["treat_03"]["lost"] = UNIT // 8
CONDITION = {s: s.split("_")[0] for s in LIBS}
BATCH = {s: ("b1" if int(s.split("_")[1]) % 2 else "b2") for s in LIBS}
DESIGN = {"condition": CONDITION, "batch": BATCH}


def expected_lost(calls):
    """What the documented lines say about each library's loss: {library: "REFUSE" | "review"}."""
    out = {}
    for s, v in calls.items():
        f = v["lost"] / v["aligner"]
        if f > REFUSE_ABOVE:
            out[s] = "REFUSE"
        elif f >= REVIEW_AT:
            out[s] = "review"
    return out


def listed(findings):
    """{library: tag} as the gate's loss finding lists them."""
    hit = [x for x in findings if x.check.startswith("aligner cells lost")]
    return {d.split(":")[0]: d.rsplit(" ", 1)[1] for d in hit[0].detail}, hit[0].severity


def stricter(findings):
    return [x for x in findings if "stricter" in x.check][0]


fails = []
print("Step 2 - cell-call gate\n" + "=" * 74)
print("A. a constructed nuclei cohort, as accepted")
f = gate(COHORT, DESIGN, assay="snrna")
for x in f:
    print(x)
v = verdict(f)
print(f"\n verdict: {v}")
got, sev = listed(f)
want = expected_lost(COHORT)
if want != {"treat_01": "review"}:
    fails.append(f"A: the construction is wrong - the lines put {want} at or over them")
if got != want:
    fails.append(f"A: the gate listed {got}; the documented lines list {want} (the 5% line is inclusive)")
if sev != "REVIEW" or v != "REVIEW":
    fails.append(f"A: a loss on the review line and nothing over the refuse line is REVIEW, got {sev}/{v}")
if not any(x.severity == "REVIEW" and x.check == "loss differential: condition" for x in f):
    fails.append("A: every lost cell on one arm of the condition must be surfaced")
if any(x.severity == "REFUSE" and x.check == "loss differential: batch" for x in f):
    fails.append("A: a loss balanced across the batch must not refuse on it")

print("\n" + "-" * 74)
print("B. injected: CellBender stricter than the aligner")
# STRICTER IS A FILTER WHEN IT COSTS CELLS (harness autodev W020, B1-103): a ratio below 1.0 with
# under REVIEW_LOST of the aligner's cells lost is REVIEW; with REVIEW_LOST or more it is REFUSED.
b = {k: dict(x) for k, x in COHORT.items()}
b["ctrl_04"]["cellbender"] = b["ctrl_04"]["aligner"] - UNIT
for lost, want_sev, why in ((0, "REVIEW", "a ratio below 1.0 that loses nothing"),
                            (on_line("ctrl_04", REVIEW_AT, b) - 1, "REVIEW",
                             "a ratio below 1.0 that loses one cell under the review line"),
                            (on_line("ctrl_04", REVIEW_AT, b), "REFUSE",
                             "a ratio below 1.0 that loses exactly the review line (inclusive)")):
    b["ctrl_04"]["lost"] = lost
    hit = stricter(gate(b, DESIGN, assay="snrna"))
    print(f" lost {lost}: {hit.severity}")
    if hit.severity != want_sev:
        fails.append(f"B: {why} must be {want_sev}, got {hit.severity}")

print("\n" + "-" * 74)
print("C. injected: losses at, just over, and well over the refuse line in one sample")
for lost_of, label in ((lambda c: on_line("ctrl_02", REFUSE_ABOVE, c), "exactly the refuse line"),
                       (lambda c: on_line("ctrl_02", REFUSE_ABOVE, c) + 1, "one cell over it"),
                       (lambda c: int(0.14 * c["ctrl_02"]["aligner"]), "14% of the aligner's cells")):
    c = {k: dict(x) for k, x in COHORT.items()}
    c["ctrl_02"]["lost"] = lost_of(c)
    got, sev = listed(gate(c, DESIGN, assay="snrna"))
    want = expected_lost(c)
    print(f" {label}: listed {got.get('ctrl_02')}, finding {sev}")
    if got != want:
        fails.append(f"C: {label}: the gate listed {got}; the documented lines list {want}")
    if sev != ("REFUSE" if "REFUSE" in want.values() else "REVIEW"):
        fails.append(f"C: {label}: the finding is {sev}")

print("\n" + "=" * 74)
if fails:
    print("FAILED:"); [print(" -", x) for x in fails]; raise SystemExit(1)
print("proved: each line holds on both sides - a loss ON the review line is REVIEWed, one cell under")
print(" it is not; ON the refuse line is REVIEWed, one cell over it REFUSED; a stricter denoiser is")
print(" REFUSED from the review line up - and a loss on one arm is seen where no per-sample % can see it")
