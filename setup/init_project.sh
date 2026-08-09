#!/usr/bin/env bash
# Create a project directory for one dataset.
#
# scQC never writes into its own directory. A project is a separate tree that scQC reads from and
# writes into, so one installation serves any number of datasets and upgrading the pipeline can
# never disturb a result.
#
# Usage:
#   setup/init_project.sh --dir ~/projects/my-study --assay snrna
#   setup/init_project.sh --dir ~/projects/my-study --assay scrna --samples 8
#
# Creates the layout, a samplesheet template, a decisions template, and a README that records
# what the directory is - so it is still interpretable by someone who did not create it.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
VERSION="$(cat "$ROOT/VERSION" 2>/dev/null || echo unknown)"

DIR=""
ASSAY=""
NSAMPLES=4
PLATFORM=""
SPECIES=""
REFERENCE=""

usage() { sed -n '2,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir)       DIR="${2:?--dir needs a path}"; shift 2 ;;
        --assay)     ASSAY="${2:?--assay needs snrna or scrna}"; shift 2 ;;
        --samples)   NSAMPLES="${2:?--samples needs a count}"; shift 2 ;;
        --platform)  PLATFORM="${2:?--platform needs 10x or singleron}"; shift 2 ;;
        --species)   SPECIES="${2:?--species needs a species name}"; shift 2 ;;
        --reference) REFERENCE="${2:?--reference needs a registry key}"; shift 2 ;;
        -h|--help)   usage 0 ;;
        *) echo "unknown argument: $1" >&2; usage 1 ;;
    esac
done

[[ -n "$DIR"   ]] || { echo "ERROR: --dir is required" >&2; usage 1; }
[[ -n "$ASSAY" ]] || { echo "ERROR: --assay is required (snrna or scrna)" >&2; usage 1; }
[[ "$ASSAY" == "snrna" || "$ASSAY" == "scrna" ]] || {
    echo "ERROR: --assay must be 'snrna' or 'scrna', got '$ASSAY'" >&2; exit 1; }
[[ "$NSAMPLES" =~ ^[0-9]+$ ]] && [[ "$NSAMPLES" -ge 1 ]] || {
    echo "ERROR: --samples must be a positive integer, got '$NSAMPLES'" >&2; exit 1; }

# platform is validated against the same whitelist the ingest step uses. Pre-filling a template
# with a value that step 0 will reject helps nobody, so an unsupported platform is refused here.
if [[ -n "$PLATFORM" ]]; then
    case "$PLATFORM" in
        10x|singleron) : ;;
        bgi) echo "ERROR: platform 'bgi' is recognised but not implemented; step 0 will refuse" >&2
             exit 1 ;;
        *) echo "ERROR: --platform must be 10x or singleron, got '$PLATFORM'" >&2; exit 1 ;;
    esac
fi

DIR="${DIR/#\~/$HOME}"

# Refuse to write into a directory that already holds a project. Overwriting a samplesheet or a
# decisions file silently would destroy the record of what a previous run was told to do.
if [[ -e "$DIR/samplesheet.csv" || -e "$DIR/decisions.yml" ]]; then
    echo "ERROR: $DIR already contains a project (samplesheet.csv or decisions.yml)." >&2
    echo "       Refusing to overwrite. Choose another --dir, or move the existing files aside." >&2
    exit 1
fi

mkdir -p "$DIR"/{data,results/{tables,figures,reports,objects},work,logs}

# --- samplesheet -------------------------------------------------------------------------------
# The four REQUIRED columns are sample, platform, species and reference; they have no defaults
# because every one of them changes how a later step behaves, and a wrong guess is silent. An
# earlier version of this script omitted platform and species, so `scqc validate` rejected every
# row of the file this script had just written - a generated template must pass the validator that
# ships beside it.
{
    echo "# scQC samplesheet. One row per library."
    echo "#"
    echo "# REQUIRED - no defaults exist for these:"
    echo "#   sample     unique library name; becomes the key in every output table"
    echo "#   platform   10x | singleron   (bgi is recognised but not implemented and will refuse)"
    echo "#   species    e.g. mus_musculus; gene-class patterns are species-specific"
    echo "#   reference  a '<species>/<build>' key from references/_registry/registry.tsv"
    echo "#"
    echo "# One of these two is required:"
    echo "#   matrix     absolute path to a RAW, UNFILTERED matrix"
    echo "#   fastq_r1   absolute path; supply when no raw matrix exists"
    echo "#"
    echo "# Optional:"
    echo "#   fastq_r2   absolute path"
    echo "#   assay      snrna | scrna. Ambient correction is mandatory for snrna."
    echo "#"
    echo "# Design factors: add one column per factor (condition, batch, timepoint, ...)."
    echo "# Any column with 2-6 distinct values is treated as a design factor and every"
    echo "# differential check is computed across it. Columns with one value are ignored."
    echo "#"
    echo "# Chemistry is deliberately absent: it is detected per sample and recorded, never"
    echo "# declared, because a declared chemistry that disagrees with the data is a silent error."
    echo "sample,platform,species,reference,assay,fastq_r1,fastq_r2,matrix,condition"
    for i in $(seq 1 "$NSAMPLES"); do
        printf 'sample_%02d,%s,%s,%s,%s,,,,\n' \
            "$i" "$PLATFORM" "$SPECIES" "$REFERENCE" "$ASSAY"
    done
} > "$DIR/samplesheet.csv"

# --- decisions ---------------------------------------------------------------------------------
cat > "$DIR/decisions.template.yml" <<'YAML'
# Decisions for this dataset.
#
# Copy to decisions.yml and fill in after reading the evidence-mode report. Do NOT fill it in
# beforehand: the point of the two-phase design is that these values are chosen against the
# report, and a value chosen before the evidence is a guess wearing a number's clothes.
#
# Every ADJUDICATED entry needs `value` AND `approved_by` AND `verbatim`. The verbatim field must
# be your own words about THIS decision. The pipeline compares them and refuses without them;
# there is no override.

quality:
  umi_floor:
    value:            # integer, bounded to [200, 1000]
    class: ADJUDICATED
    approved_by:
    verbatim:
  gene_floor:
    value:            # integer, bounded to [100, 600]
    class: ADJUDICATED
    approved_by:
    verbatim:
  mito_ceiling_pct:
    value:            # percentage; rarely derivable - usually a judgement
    class: ADJUDICATED
    approved_by:
    verbatim:

doublets:
  detector: scDblFinder
  dbr:                # expected doublet rate, or null to use the detector's default
  dbr_sd:             # uncertainty on that rate

cluster_check:
  resolution: 1.0
  a_umi_fraction: 0.5     # cluster median UMI below this fraction of the sample median
  b_mito_pct: 15.0        # cluster median mitochondrial percentage above this
  c_uninformative_pct: 50.0
  d_doublet_pct: 70.0
  approved_by:
  verbatim:

apply:
  # The exact action string this approval authorises. Changing any threshold above changes this
  # string, which invalidates the approval - by design.
  action:
  approved_by:
  verbatim:
YAML

# --- README ------------------------------------------------------------------------------------
cat > "$DIR/README.md" <<EOF
# $(basename "$DIR")

Created by scQC $VERSION on $(date -u +%Y-%m-%d).

## What this directory is

A project directory for one dataset. scQC reads \`samplesheet.csv\` and \`decisions.yml\` from
here and writes everything it produces into \`results/\`. The pipeline's own installation is
elsewhere and is never modified by a run.

## Layout

| path | contents |
|---|---|
| \`data/\` | inputs. Treat as read-only: never modify a delivered file |
| \`samplesheet.csv\` | one row per library, plus design factor columns |
| \`decisions.yml\` | operator decisions. Copy from \`decisions.template.yml\` |
| \`results/tables/\` | CSV/TSV, small enough to open and read |
| \`results/figures/\` | figures |
| \`results/reports/\` | rendered HTML reports and \`report.json\` |
| \`results/objects/\` | large binaries (\`.h5ad\`). Exclude from any sync client |
| \`work/\` | scratch. Safe to delete |
| \`logs/\` | run logs |

## Assay

\`$ASSAY\`. $( [[ "$ASSAY" == "snrna" ]] \
  && echo "Ambient RNA correction is MANDATORY: nuclear preparations carry substantial ambient
signal from lysed cells, and skipping it produces spurious cross-population expression." \
  || echo "Ambient RNA correction is optional but recommended." )

## How to run

\`\`\`bash
# 1. Fill in samplesheet.csv, then look at the data without touching it:
scqc --project . --mode evidence

# 2. Read results/reports/, then record your decisions:
cp decisions.template.yml decisions.yml && \$EDITOR decisions.yml

# 3. Apply them:
scqc --project . --mode apply
\`\`\`

## Status

- [ ] samplesheet.csv filled in
- [ ] evidence mode run
- [ ] report read
- [ ] decisions.yml written
- [ ] apply mode run
EOF

cat <<EOF
Created $DIR

  samplesheet.csv         $NSAMPLES row(s), assay=$ASSAY  <- fill this in first
  decisions.template.yml  copy to decisions.yml AFTER reading the evidence report
  README.md               what this directory is
  data/ results/ work/ logs/

Next:
  \$EDITOR $DIR/samplesheet.csv
  scqc --project $DIR --mode evidence
EOF
