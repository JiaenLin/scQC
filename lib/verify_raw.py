# Input verification for step 0: reads matrix summary statistics and returns a verdict.
# It removes nothing, writes nothing, and subsets nothing; its whole purpose is to refuse an
# input rather than to alter one.
"""Step 0 of the pipeline: is this matrix actually usable as the pipeline's input?

WHY THIS EXISTS. "Raw" is ambiguous, and the ambiguity is expensive because it is discovered, if
at all, only after a matrix has been analysed. Two INDEPENDENT properties hide under the word:

    P1  RAW VALUES      unnormalised integer counts, no ceiling, no gene subsetting
    P2  RAW DROPLETS    every barcode, including the empties

A matrix can pass P1 and fail P2, and that combination is the common one: a delivered count
matrix is frequently an aligner's `outs/filtered` - genuine integer counts (P1 pass) that have
already been through cell calling (P2 fail). Nothing in the file says so. There is no header
field, no flag, and no difference in shape; separating the two takes an audit of the value
distributions, which is what this module is.

P2 is the one that matters for what comes next: an ambient-RNA model learns the background
profile FROM the empty droplets. Hand it a cell-called matrix and it has nothing to learn from -
and the conclusion "ambient correction is not possible with these files" is then available,
wrong, and indistinguishable on the page from a real finding.

THE VERDICT DECIDES THE MODE, NOT THE USER. A supplied matrix that fails is not a warning to be
clicked through: step 0 rejects the matrix and falls back to running the processor from FASTQ.
"""

from __future__ import annotations

from dataclasses import dataclass, field

@dataclass
class Verdict:
    name: str
    p1_raw_values: bool = True
    p2_raw_droplets: bool = True
    reasons: list = field(default_factory=list)

    @property
    def usable(self) -> bool:
        return self.p1_raw_values and self.p2_raw_droplets

    def __str__(self) -> str:
        head = "USABLE" if self.usable else "REJECTED"
        bits = [f"[{head:8s}] {self.name}",
                f" P1 raw values : {'pass' if self.p1_raw_values else 'FAIL'}",
                f" P2 raw droplets : {'pass' if self.p2_raw_droplets else 'FAIL'}"]
        bits += [f" - {r}" for r in self.reasons]
        return "\n".join(bits)

def verify(name, n_barcodes, n_genes, min_counts, max_counts, p98_counts,
           expected_genes=None, integer_counts=True) -> Verdict:
    """Verdict on one matrix from summary statistics a caller can compute cheaply."""
    v = Verdict(name)

    # ---- P1: are the VALUES raw?
    if not integer_counts:
        v.p1_raw_values = False
        v.reasons.append("values are not integers - normalised or scaled")
    # A ceiling shows up as the maximum sitting on a percentile of itself. A genuine maximum is
    # an outlier and does not coincide with p98 to four significant figures.
    if p98_counts and max_counts and abs(max_counts - p98_counts) / max(p98_counts, 1) < 0.01:
        v.p1_raw_values = False
        v.reasons.append(
            f"max UMI {max_counts:,.0f} sits on the 98th percentile {p98_counts:,.0f} - "
            f"an upper cap was applied. A UMI ceiling is what doi:10.1101/140848 used AS a "
            f"doublet filter, so this silently pre-removes the deepest nuclei")
    if expected_genes and n_genes < expected_genes:
        v.p1_raw_values = False
        v.reasons.append(
            f"{n_genes:,} genes against {expected_genes:,} in the declared reference - "
            f"{expected_genes - n_genes:,} genes are absent and nothing downstream recovers them")

    # ---- P2: are the DROPLETS raw?
    # The decisive tell. A raw droplet matrix contains barcodes with a handful of UMI; a
    # cell-called one has a floor, and that floor is usually a round number somebody chose.
    if min_counts is not None and min_counts >= 100:
        v.p2_raw_droplets = False
        v.reasons.append(
            f"minimum UMI is {min_counts:,.0f} - a raw droplet matrix contains near-empty "
            f"barcodes with single-digit counts. This has been cell-called, so the empty "
            f"droplets CellBender needs to learn the ambient profile are already gone")
    return v

def demo() -> None:
    """Three worked verdicts on INVENTED statistics: no files, no data, no dependencies.

    The numbers below are illustrative and are not measurements of anything. What they
    demonstrate is the shape of each verdict, and in particular the middle case - the one that
    passes P1 and fails P2, which is the failure this module exists for.
    """
    ref_genes = 25_000 # an illustrative declared reference size

    print("verify_raw demo - ALL NUMBERS BELOW ARE INVENTED\n" + "=" * 74)

    print(verify("A delivered analysis object: capped and gene-subset",
                 n_barcodes=12_000, n_genes=20_000,
                 min_counts=500, max_counts=18_400, p98_counts=18_400,
                 expected_genes=ref_genes))
    print()
    print(verify("A count matrix from the aligner's filtered output",
                 n_barcodes=12_500, n_genes=ref_genes,
                 min_counts=486, max_counts=71_204, p98_counts=24_900,
                 expected_genes=ref_genes))
    print()
    print(verify("An unfiltered droplet matrix",
                 n_barcodes=310_000, n_genes=ref_genes,
                 min_counts=1, max_counts=71_204, p98_counts=None,
                 expected_genes=ref_genes))

    print("\n" + "=" * 74)
    print("The middle verdict is the one to read. Its values are genuine integer counts, so it "
          "passes P1\nand looks like a usable input - but its floor of 486 UMI says the empty "
          "droplets are already\ngone, and those droplets are what an ambient model learns from.")

if __name__ == "__main__":
    import sys

    args = sys.argv[1:]
    if args == ["--demo"]:
        demo()
    elif not args:
        print("verify_raw is a library: import `verify` and pass it summary statistics you have "
              "computed.\nFor a self-contained illustration of the three verdicts, run:\n"
              "    python lib/verify_raw.py --demo")
    else:
        sys.stderr.write(f"unrecognised argument(s): {' '.join(args)}\n"
                         f"usage: python lib/verify_raw.py [--demo]\n")
        raise SystemExit(2)
