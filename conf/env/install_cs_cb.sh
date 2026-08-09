#!/bin/bash
# Install the CeleScope aligner and the CellBender denoiser as two separate conda environments.
#
# Usage:
#   conf/env/install_cs_cb.sh <env-root> [celescope] [cellbender]
#
# <env-root> is the directory given to setup/install_env.sh --prefix. Environment names are fixed
# by the pipeline: <env-root>/celescope and <env-root>/cellbender. Naming no component builds
# both. An environment that already exists is left alone.
#
# WHY TWO ENVIRONMENTS. The aligner's solve pins an older numpy than the denoiser's torch build
# tolerates. One environment holding both resolves only by relaxing a pin, and a relaxed pin is
# how a pipeline stops reproducing.
#
# WHY THE ALIGNER VERSION IS PINNED. Quantification is only comparable across runs when the
# aligner version is fixed; with the pin in place, a difference between two matrices is a
# property of the run rather than of the software that produced it. 2.7.3 is the default here
# because it is the version the shipped reference recipe and argument checks were written
# against. Changing it is a deliberate act, and it belongs in the record of the analysis.
#
# The conda spec comes from the v2.7.3 TAG, not from master. They differ: 2.7.3 pins trust4=1.0.7
# and leaves snpeff unpinned, while master pins snpeff=5.2 and adds leidenalg, igraph and pyarrow.
# A spec fetched from master therefore installs a different environment under the same version
# number. The tag's spec is shipped alongside this script as conda_pkgs_v2.7.3.txt and is used
# in preference to any download.
#
# WHY micromamba IS BOOTSTRAPPED WHEN NOTHING ELSE IS FOUND. conda releases predating the libmamba
# solver take hours on a bioconda environment this size, or fail to solve it at all. micromamba is
# a single static binary, needs no privileges, and installs into <env-root> like everything else.
#
# Environment:
#   SCQC_CONDA        package manager to use; setup/install_env.sh passes the one it discovered
#   TORCH_INDEX_URL   wheel index for torch (default the CUDA 12.1 build)

set -uo pipefail

ROOT="${1:-}"
if [ -z "$ROOT" ] || [ "$ROOT" = "-h" ] || [ "$ROOT" = "--help" ]; then
    echo "usage: install_cs_cb.sh <env-root> [celescope] [cellbender]" >&2
    exit 2
fi
case "$ROOT" in
    /*) : ;;
    *)  echo "ERROR: <env-root> must be an absolute path (got '$ROOT')" >&2; exit 2 ;;
esac
shift

DO_CELESCOPE=0
DO_CELLBENDER=0
if [ "$#" -eq 0 ]; then
    DO_CELESCOPE=1
    DO_CELLBENDER=1
else
    for c in "$@"; do
        case "$c" in
            celescope)  DO_CELESCOPE=1 ;;
            cellbender) DO_CELLBENDER=1 ;;
            *) echo "ERROR: unknown component '$c' (expected celescope or cellbender)" >&2; exit 2 ;;
        esac
    done
fi

BIN="$ROOT/bin"
LOGDIR="$ROOT/logs"
LOG="$LOGDIR/install_cs_cb_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$BIN" "$LOGDIR" || exit 1
exec > >(tee "$LOG") 2>&1

echo "=================================================================="
echo "CeleScope / CellBender install"
echo "host      : $(hostname)"
echo "started   : $(date)"
echo "env root  : $ROOT"
echo "components: $([ $DO_CELESCOPE -eq 1 ] && printf 'celescope ')$([ $DO_CELLBENDER -eq 1 ] && printf 'cellbender')"
echo "log       : $LOG"
echo "=================================================================="

# ------------------------------------------------------------------ package manager
MM="${SCQC_CONDA:-}"
if [ -z "$MM" ]; then
    for c in micromamba mamba conda; do
        if command -v "$c" >/dev/null 2>&1; then MM="$(command -v "$c")"; break; fi
    done
fi
if [ -z "$MM" ]; then
    echo; echo "--- no conda, mamba or micromamba found; bootstrapping micromamba into $BIN ---"
    ( cd "$ROOT" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba )
    chmod +x "$BIN/micromamba" 2>/dev/null
    MM="$BIN/micromamba"
fi
[ -x "$MM" ] || { echo "ERROR: package manager '$MM' is not executable" >&2; exit 1; }
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$ROOT/mamba}"
echo "package manager: $MM  ($("$MM" --version 2>&1 | head -1))"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPEC_SHIPPED="$HERE/conda_pkgs_v2.7.3.txt"
SPEC="$ROOT/conda_pkgs_v2.7.3.txt"

# ------------------------------------------------------------------ CeleScope
CEL="$ROOT/celescope"
if [ $DO_CELESCOPE -eq 1 ]; then
    if [ ! -d "$CEL" ]; then
        if [ -s "$SPEC_SHIPPED" ]; then
            echo; echo "--- using the shipped v2.7.3 conda spec ---"
            cp "$SPEC_SHIPPED" "$SPEC"
        else
            echo; echo "--- fetching the v2.7.3 conda spec (NOT master) ---"
            curl -fsSL -o "$SPEC" \
                https://raw.githubusercontent.com/singleron-RD/CeleScope/v2.7.3/conda_pkgs.txt
        fi
        echo "spec:"; sed 's/^/    /' "$SPEC"

        echo; echo "--- creating celescope env ---"
        "$MM" create -y -p "$CEL" -c conda-forge -c bioconda --file "$SPEC"
        echo "conda layer exit: $?"

        echo; echo "--- pip install celescope==2.7.3 ---"
        "$CEL/bin/pip" install --no-cache-dir "celescope==2.7.3"
        echo "pip exit: $?"
    else
        echo; echo "--- celescope env already present at $CEL, skipping create ---"
    fi
fi

# ------------------------------------------------------------------ CellBender
CB="$ROOT/cellbender"
if [ $DO_CELLBENDER -eq 1 ]; then
    if [ ! -d "$CB" ]; then
        echo; echo "--- creating cellbender env (python 3.9, per the upstream dev install docs) ---"
        "$MM" create -y -p "$CB" -c conda-forge python=3.9 pytables

        # The torch wheel has to match the CUDA runtime on the machine that will RUN the denoiser,
        # which is not necessarily the machine installing it. The CUDA 12.1 build works against any
        # 12.x driver stack; override TORCH_INDEX_URL for a different one, or point it at plain
        # PyPI for a CPU-only install.
        TORCH_INDEX_URL="${TORCH_INDEX_URL:-https://download.pytorch.org/whl/cu121}"
        echo; echo "--- torch from $TORCH_INDEX_URL ---"
        "$CB/bin/pip" install --no-cache-dir torch --index-url "$TORCH_INDEX_URL"

        echo; echo "--- cellbender 0.3.2 ---"
        "$CB/bin/pip" install --no-cache-dir "cellbender==0.3.2"
    else
        echo; echo "--- cellbender env already present at $CB, skipping create ---"
    fi
fi

# ------------------------------------------------------------------ VERIFY
echo
echo "=================================================================="
echo "VERIFICATION"
echo "=================================================================="
FAIL=0
chk() { if [ "$2" = "0" ]; then echo "  PASS  $1"; else echo "  FAIL  $1"; FAIL=$((FAIL+1)); fi; }

if [ $DO_CELESCOPE -eq 1 ]; then
    echo "--- CeleScope ---"
    V=$("$CEL/bin/celescope" --version 2>&1 | tr -d '\r'); echo "  celescope --version : $V"
    echo "$V" | grep -q "2\.7\.3"; chk "celescope version is 2.7.3" $?
    S=$("$CEL/bin/STAR" --version 2>&1 | head -1); echo "  STAR --version      : $S"
    echo "$S" | grep -q "2\.7\.11a"; chk "STAR is 2.7.11a (the v2.7.3 pin)" $?
    "$CEL/bin/python" -c "import celescope; print('  celescope import OK', getattr(celescope, '__version__', ''))"; chk "celescope imports" $?
    "$CEL/bin/multi_rna" -h >/dev/null 2>&1; chk "multi_rna -h runs" $?
    "$CEL/bin/celescope" rna mkref -h >/dev/null 2>&1; chk "celescope rna mkref -h runs" $?
    for t in samtools subread featureCounts picard gatk; do
        command -v "$CEL/bin/$t" >/dev/null 2>&1 && echo "    present: $t"
    done

    # The alignment step sets these by name. Printing them here means a spelling or a removal is
    # caught at install time rather than inside a scheduled job.
    echo
    echo "--- the multi_rna arguments the pipeline sets ---"
    "$CEL/bin/multi_rna" -h 2>&1 | grep -iE "chemistry|soloFeatures|soloCellFilter|genomeDir|mapfile|thread|outdir" | sed 's/^/    /' | head -20
fi

if [ $DO_CELLBENDER -eq 1 ]; then
    echo
    echo "--- CellBender ---"
    CV=$("$CB/bin/cellbender" --version 2>&1 | head -1); echo "  cellbender --version: $CV"
    "$CB/bin/cellbender" remove-background --help >/dev/null 2>&1; chk "cellbender remove-background --help runs" $?
    "$CB/bin/python" -c "
import torch
print('  torch               :', torch.__version__)
print('  torch CUDA build    :', torch.version.cuda)
print('  cuda available HERE :', torch.cuda.is_available())
"; chk "torch imports" $?

    echo
    echo "--- the remove-background flags the pipeline sets ---"
    "$CB/bin/cellbender" remove-background --help 2>&1 | grep -E "expected-cells|total-droplets|fpr|epochs|cuda|input|output" | sed 's/^/    /' | head -14
fi

echo
echo "=================================================================="
echo "checks failed: $FAIL"
if [ $DO_CELLBENDER -eq 1 ]; then
    echo "NOTE: a machine with no visible GPU reports cuda available = False. That proves nothing"
    echo "      about the machine the denoiser will run on; check it again inside a GPU job."
fi
echo "finished: $(date)"
echo "=================================================================="
[ "$FAIL" -eq 0 ] || exit 1
