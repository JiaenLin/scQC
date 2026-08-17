"""The nuclear fraction, read from an aligner's per-barcode read summary.

    nuclear fraction = intronic / (intronic + exonic)

WHAT IT MEASURES, AND WHY IT IS PAIRED WITH THE MITOCHONDRIAL PERCENTAGE
------------------------------------------------------------------------
On single-NUCLEUS data a high mitochondrial percentage cannot mean what it means on whole cells: a
nucleus contains no mitochondria, so the reads are cytoplasmic carry-over and ambient RNA in the
droplet, not damage to the thing being measured. `pct_counts_mt` therefore describes the DROPLET,
and on its own it cannot separate

    an intact nucleus in a dirty suspension        from        debris and ambient RNA

The nuclear fraction is the axis that can, because nuclear pre-mRNA is intronic: a droplet whose
reads are largely intronic contains a nucleus whatever else is floating around it. The two are
useful together and misleading apart, which is why nothing here is a filter by itself.

WHAT THIS MODULE DOES NOT DO
----------------------------
It defines no threshold, proposes none, and removes nothing. It reads a file and returns numbers.

An UNDEFINED nuclear fraction — zero intronic and zero exonic — is returned as `None` and must
stay that way. It is not a low fraction and not a high one; it is a barcode the measurement could
not be made on, and a caller that coerces it to 0.0 has invented evidence of a droplet with no
nucleus in it.

THE FILE
--------
STARsolo's `CellReads.stats`: one row per barcode, tab-separated, with a header naming the columns.
Column POSITIONS are not assumed — they are read from the header — because a column order that
happens to hold today is a hypothesis about the aligner's next version.

Two things in it that are not barcodes and must not be summed as if they were: STARsolo writes
aggregate rows such as `CBnotInPasslist` for reads it could not assign. They are named and skipped.

Antisense counts (`exonicAS`, `intronicAS`) are read only when asked for. Including them silently
would change what the number means without changing its name.
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

#: Rows STARsolo writes that are aggregates rather than barcodes. Summed as barcodes they would
#: contribute a single enormous "cell" to any pooled statistic.
NOT_A_BARCODE = frozenset({
    "CBnotInPasslist", "noCBmatch", "noCB", "CBmatch", "unmapped", "TOTAL", "Total",
})

#: The columns summed for each half of the ratio, sense-only and with antisense.
SENSE = {"intronic": ("intronic",), "exonic": ("exonic",)}
WITH_ANTISENSE = {"intronic": ("intronic", "intronicAS"), "exonic": ("exonic", "exonicAS")}


class NuclearFractionError(RuntimeError):
    """The nuclear fraction could not be read. Never raised for a merely undefined value."""


def nuclear_fraction(intronic, exonic):
    """`intronic / (intronic + exonic)`, or None where the denominator is zero.

    None is the answer for a barcode with no assignable reads, and it is a THIRD state. Returning
    0.0 would say the droplet held no nuclear signal, which is a measurement nobody made.
    """
    try:
        i, e = float(intronic), float(exonic)
    except (TypeError, ValueError):
        return None
    total = i + e
    if total <= 0:
        return None
    return i / total


def _open_text(path: Path):
    """Open plain or gzipped text. The aligner writes either, depending on the wrapper."""
    p = Path(path)
    if p.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")
    return p.open("r", encoding="utf-8", errors="replace")


def _bare(barcode: str, sample) -> str:
    """The barcode as the ALIGNER spells it: without the `<sample>_` prefix the object carries.

    The object's barcodes are prefixed so that libraries can be concatenated without collision;
    `CellReads.stats` predates that and is not. Joining the two spellings directly gives an
    intersection of exactly ZERO barcodes while both files describe the same library — an error
    that looks like missing data rather than like a mismatch, and one this pipeline has already
    made once between a matrix and a cell-barcode CSV.
    """
    if sample:
        pre = f"{sample}_"
        if barcode.startswith(pre):
            return barcode[len(pre):]
    return barcode


def read_cellreads(path, target_barcodes, *, sample=None, antisense=False) -> tuple:
    """`(values_by_target_barcode, stats)` for the barcodes asked for, and nothing else.

    `target_barcodes` is the barcodes as the OBJECT spells them; the returned dict is keyed the
    same way, so the caller never handles two spellings.

    STREAMED AND FILTERED ON READ, and that is the whole performance argument. The file holds one
    row per DROPLET — on a typical run 1.7 million of them, 131 MB — while an object holds tens of
    thousands of called cells. Building a dict of every row and then selecting from it costs ~45x
    the memory for the same answer. Each line is split ONCE to look at the barcode, and split the
    rest of the way only for the ~2% that are wanted.
    """
    p = Path(path)
    if not p.exists():
        raise NuclearFractionError(f"nuclear-fraction source does not exist: {p}")

    cols = WITH_ANTISENSE if antisense else SENSE
    wanted = {_bare(str(b), sample): str(b) for b in target_barcodes}
    if not wanted:
        raise NuclearFractionError(
            f"{p}: no target barcodes were given, so there is nothing to join against. An empty "
            f"request is a caller error, not a library with no cells.")

    out: dict = {}
    n_rows = n_skipped = 0
    with _open_text(p) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        need = ["CB"] + [c for pair in cols.values() for c in pair]
        missing = [c for c in need if c not in idx]
        if missing:
            raise NuclearFractionError(
                f"{p}: header has no column(s) {missing}. Present: {header[:12]}... "
                f"Column positions are read from the header and never assumed, so this is a "
                f"different file format rather than a different column order.")
        cb_i = idx["CB"]
        i_idx = [idx[c] for c in cols["intronic"]]
        e_idx = [idx[c] for c in cols["exonic"]]
        # The barcode is field 0 in every STARsolo build seen, but the header is authoritative:
        # the cheap pre-split is only valid when it really is first.
        cheap = cb_i == 0

        for line in fh:
            n_rows += 1
            if cheap:
                cb = line.split("\t", 1)[0]
                if cb not in wanted:
                    continue
                parts = line.rstrip("\n").split("\t")
            else:
                parts = line.rstrip("\n").split("\t")
                cb = parts[cb_i] if cb_i < len(parts) else ""
                if cb not in wanted:
                    continue
            if cb in NOT_A_BARCODE:
                n_skipped += 1
                continue
            try:
                intronic = sum(float(parts[i]) for i in i_idx)
                exonic = sum(float(parts[i]) for i in e_idx)
            except (IndexError, ValueError):
                # A malformed row is not a zero. It is recorded as unjoined and the barcode keeps
                # no value, so it fails nothing downstream.
                continue
            out[wanted[cb]] = nuclear_fraction(intronic, exonic)

    n_defined = sum(1 for v in out.values() if v is not None)
    stats = {
        "source": str(p),
        "columns": "+".join(cols["intronic"]) + " / (" + "+".join(cols["intronic"]) + " + "
                   + "+".join(cols["exonic"]) + ")",
        "antisense": bool(antisense),
        "rows_read": n_rows,
        "aggregate_rows_skipped": n_skipped,
        "n_target": len(wanted),
        "n_joined": len(out),
        "n_defined": n_defined,
        "join_pct": (100.0 * len(out) / len(wanted)) if wanted else 0.0,
    }
    return out, stats


def median(values) -> float | None:
    """The median of the DEFINED values, or None if there are none.

    Written out rather than imported from `statistics` so that a None in the input is a visible
    filter here rather than a TypeError three frames away.
    """
    xs = sorted(float(v) for v in values if v is not None)
    if not xs:
        return None
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])
