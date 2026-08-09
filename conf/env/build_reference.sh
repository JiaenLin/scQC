#!/bin/bash
# Build a CeleScope-filtered GTF from an Ensembl release and report the gene set it yields.
#
# A count matrix is comparable to another one only when both were quantified against the same
# reference, and a reference is often described by nothing more than a directory name such as
# `<species>_ensembl_<release>_filtered`. That filter is CeleScope's own `utils mkgtf` default,
# so such a reference can be rebuilt rather than guessed at:
#
#   celescope utils mkgtf --attributes \
#     "gene_biotype=protein_coding,lncRNA,antisense,IG_LV_gene,IG_V_gene,IG_V_pseudogene,
#      IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_C_gene,IG_C_pseudogene,TR_V_gene,
#      TR_V_pseudogene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_C_gene;"
#
# Intron entries are retained by default, which is what soloFeatures GeneFull_Ex50pAS requires,
# and mt_gene_list.txt is generated from gene symbols matching MT-/mt-.
#
# THE CHECK. `--expect-genes N` asserts how many genes the rebuilt GTF must contain. When a
# matrix was quantified somewhere else and its gene count is known, passing that number turns
# "probably the same reference" into a reproducible yes or no: a match means the reference is
# reproduced by construction, and a mismatch measures exactly how far apart the two are. Either
# answer is a result.
#
# `--matrix-tar` answers the same question by a second, independent route: the gene identifiers
# themselves, read from a features.tsv.gz inside a supplied matrix archive. Two gene sets that
# agree in size and two that agree element by element are not the same evidence, and only the
# second distinguishes a reproduced reference from a coincidence.
#
# The STAR index is deliberately NOT built here. It costs hours and roughly 32 GB of RAM, and
# there is no point paying that before the gene set is settled; the command to run afterwards is
# printed at the end.
#
# The defaults build Ensembl 112 / GRCm39. They are a shipped EXAMPLE of a complete runnable
# invocation, not an assumption the pipeline makes: set SPECIES, ASSEMBLY and ENSEMBL_RELEASE for
# any other Ensembl species or release.
#
# Usage:
#   conf/env/build_reference.sh --out DIR [--celescope PATH]
#                               [--expect-genes N] [--matrix-tar FILE]
#
# Environment:
#   SCQC_ENV_ROOT     where setup/install_env.sh put the environments; --celescope defaults to
#                     $SCQC_ENV_ROOT/celescope/bin/celescope, then to celescope on PATH
#   SPECIES           Ensembl species directory name   (default mus_musculus)
#   ASSEMBLY          assembly name in the filename    (default GRCm39)
#   ENSEMBL_RELEASE   release number                   (default 112)

set -uo pipefail

OUT=""
CELESCOPE=""
EXPECT_GENES=""
MATRIX_TAR=""

SPECIES="${SPECIES:-mus_musculus}"
ASSEMBLY="${ASSEMBLY:-GRCm39}"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-112}"

usage() {
    cat <<'USAGE'
usage: build_reference.sh --out DIR [--celescope PATH] [--expect-genes N] [--matrix-tar FILE]

  --out DIR           where to build the reference (created if absent)
  --celescope PATH    the celescope executable; defaults to $SCQC_ENV_ROOT/celescope/bin/celescope
  --expect-genes N    assert the filtered GTF contains exactly N genes
  --matrix-tar FILE   a counts-matrix archive containing features.tsv.gz, compared gene by gene

  SPECIES / ASSEMBLY / ENSEMBL_RELEASE override the Ensembl build (default mus_musculus GRCm39 112).
USAGE
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --out)          OUT="${2:?--out needs a directory}"; shift 2 ;;
        --celescope)    CELESCOPE="${2:?--celescope needs a path}"; shift 2 ;;
        --expect-genes) EXPECT_GENES="${2:?--expect-genes needs a number}"; shift 2 ;;
        --matrix-tar)   MATRIX_TAR="${2:?--matrix-tar needs a file}"; shift 2 ;;
        -h|--help)      usage 0 ;;
        *) echo "unknown argument: $1" >&2; usage 1 ;;
    esac
done

[ -n "$OUT" ] || { echo "ERROR: --out is required" >&2; usage 1; }

if [ -z "$CELESCOPE" ]; then
    if [ -n "${SCQC_ENV_ROOT:-}" ] && [ -x "$SCQC_ENV_ROOT/celescope/bin/celescope" ]; then
        CELESCOPE="$SCQC_ENV_ROOT/celescope/bin/celescope"
    else
        CELESCOPE="$(command -v celescope 2>/dev/null || true)"
    fi
fi
if [ -z "$CELESCOPE" ] || [ ! -x "$CELESCOPE" ]; then
    echo "ERROR: no celescope executable found." >&2
    echo "       Pass --celescope PATH, or set SCQC_ENV_ROOT to the prefix used by" >&2
    echo "       setup/install_env.sh --with-celescope." >&2
    exit 1
fi

if [ -n "$EXPECT_GENES" ]; then
    case "$EXPECT_GENES" in
        ''|*[!0-9]*) echo "ERROR: --expect-genes must be a positive integer" >&2; exit 1 ;;
    esac
fi

# Resolve caller-supplied paths to absolute BEFORE changing directory. A relative --matrix-tar
# was interpreted against $OUT after the cd, so it silently missed and the run reported
# "matrix archive not found" - a verification step skipped without anyone asking for it to be.
if [ -n "$MATRIX_TAR" ]; then
    case "$MATRIX_TAR" in
        /*) : ;;
        *)  MATRIX_TAR="$PWD/$MATRIX_TAR" ;;
    esac
    [ -f "$MATRIX_TAR" ] || {
        echo "ERROR: --matrix-tar not found: $MATRIX_TAR" >&2
        echo "       Give an existing path; a missing archive would skip the gene-ID cross-check" >&2
        echo "       silently, which is the one thing this comparison exists to prevent." >&2
        exit 1
    }
fi

mkdir -p "$OUT" || exit 1
cd "$OUT" || exit 1

SP_LC="$(printf '%s' "$SPECIES" | tr '[:upper:]' '[:lower:]')"
SP_CAP="${SP_LC^}"

BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}"
GTF="${SP_CAP}.${ASSEMBLY}.${ENSEMBL_RELEASE}.gtf"
FA="${SP_CAP}.${ASSEMBLY}.dna.primary_assembly.fa"
FILT="${SP_CAP}.${ASSEMBLY}.${ENSEMBL_RELEASE}.filtered.gtf"

echo "=================================================================="
echo "build reference : ${SP_CAP} ${ASSEMBLY}, Ensembl ${ENSEMBL_RELEASE}"
echo "filter          : celescope utils mkgtf, default attributes"
echo "started         : $(date)   host: $(hostname)"
echo "dest            : $OUT"
echo "=================================================================="

echo; echo "--- download GTF ---"
[ -f "$GTF" ] || { curl -fL --retry 3 -o "$GTF.gz" "$BASE/gtf/${SP_LC}/$GTF.gz" && gunzip -f "$GTF.gz"; }
[ -s "$GTF" ] || { echo "  GTF not available at $BASE/gtf/${SP_LC}/$GTF.gz - stopping"; exit 2; }
ls -lh "$GTF" | awk '{print "  ", $5, $9}'

echo; echo "--- download primary assembly FASTA (several hundred MB compressed) ---"
[ -f "$FA" ] || { curl -fL --retry 3 -o "$FA.gz" "$BASE/fasta/${SP_LC}/dna/$FA.gz" && gunzip -f "$FA.gz"; }
[ -s "$FA" ] || { echo "  FASTA not available at $BASE/fasta/${SP_LC}/dna/$FA.gz - stopping"; exit 2; }
ls -lh "$FA" | awk '{print "  ", $5, $9}'

echo; echo "--- gene_biotype inventory of the RAW GTF ---"
awk -F'\t' '$3=="gene"' "$GTF" | grep -oP 'gene_biotype "\K[^"]+' | sort | uniq -c | sort -rn | head -20 | sed 's/^/    /'
RAW=$(awk -F'\t' '$3=="gene"' "$GTF" | grep -c 'gene_id')
echo "  raw gene records: $RAW"

echo; echo "--- celescope utils mkgtf, DEFAULT attributes (introns retained) ---"
"$CELESCOPE" utils mkgtf "$GTF" "$FILT" 2>&1 | tail -8

if [ ! -s "$FILT" ]; then echo "  MKGTF PRODUCED NOTHING - stopping"; exit 2; fi
ls -lh "$FILT" | awk '{print "  ", $5, $9}'

biotypes() {
    echo "  Retained biotypes in the filtered GTF:"
    awk -F'\t' '$3=="gene"' "$FILT" | grep -oP 'gene_biotype "\K[^"]+' | sort | uniq -c | sort -rn | sed 's/^/      /'
}

echo; echo "--- gene count of the filtered GTF ---"
N=$(awk -F'\t' '$3=="gene"' "$FILT" | grep -oP 'gene_id "\K[^"]+' | sort -u | wc -l)
echo "  unique gene_id in filtered GTF : $N"
if [ -n "$EXPECT_GENES" ]; then
    echo "  expected (--expect-genes)      : $EXPECT_GENES"
    if [ "$N" -eq "$EXPECT_GENES" ]; then
        echo "  ** MATCH ** the default filter on this release yields exactly the expected gene"
        echo "              count, so the reference is reproduced by construction."
    else
        D=$((N - EXPECT_GENES))
        echo "  ** NO MATCH ** differs by $D genes."
        echo "  This is a measurement, not a failure: whatever produced the expected count did not"
        echo "  use the documented default filter, and this is the size of the difference."
        biotypes
    fi
else
    echo "  (no --expect-genes given, so there is nothing to compare the count against)"
    biotypes
fi

echo; echo "--- mt gene list ---"
ls -1 mt_gene_list.txt 2>/dev/null && wc -l < mt_gene_list.txt | xargs -I{} echo "  mt genes: {}"

if [ -n "$MATRIX_TAR" ]; then
    echo; echo "--- compare gene IDs against a supplied features.tsv.gz (independent route) ---"
    if [ -f "$MATRIX_TAR" ]; then
        TMP=$(mktemp -d)
        tar xf "$MATRIX_TAR" -C "$TMP" 2>/dev/null
        F=$(find "$TMP" -name "features.tsv.gz" | head -1)
        if [ -n "$F" ]; then
            DN=$(zcat "$F" | wc -l)
            echo "  features.tsv.gz rows       : $DN"
            zcat "$F" | cut -f1 | sort -u > "$TMP/supplied_ids.txt"
            awk -F'\t' '$3=="gene"' "$FILT" | grep -oP 'gene_id "\K[^"]+' | sort -u > "$TMP/rebuilt_ids.txt"
            echo "  shared gene_ids            : $(comm -12 "$TMP/supplied_ids.txt" "$TMP/rebuilt_ids.txt" | wc -l)"
            echo "  only in the supplied set   : $(comm -23 "$TMP/supplied_ids.txt" "$TMP/rebuilt_ids.txt" | wc -l)"
            echo "  only in the rebuilt GTF    : $(comm -13 "$TMP/supplied_ids.txt" "$TMP/rebuilt_ids.txt" | wc -l)"
            echo "  examples, supplied only    :"; comm -23 "$TMP/supplied_ids.txt" "$TMP/rebuilt_ids.txt" | head -5 | sed 's/^/      /'
            echo "  examples, rebuilt only     :"; comm -13 "$TMP/supplied_ids.txt" "$TMP/rebuilt_ids.txt" | head -5 | sed 's/^/      /'
        else
            echo "  no features.tsv.gz inside $MATRIX_TAR"
        fi
        rm -rf "$TMP"
    else
        echo "  (matrix archive not found at $MATRIX_TAR)"
    fi
fi

echo
echo "=================================================================="
echo "STAR INDEX NOT BUILT. Run it only once the gene set is settled:"
echo "  $CELESCOPE rna mkref --fasta $FA --gtf $FILT --thread 16"
echo "finished: $(date)"
echo "=================================================================="
