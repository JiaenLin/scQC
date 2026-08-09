#!/bin/bash
# Install R with scDblFinder and scds as a separate conda environment, from precompiled packages.
#
# Usage:
#   conf/env/install_rdoublet.sh <env-root>
#
# <env-root> is the directory given to setup/install_env.sh --prefix. The environment name is
# fixed by the pipeline: <env-root>/rdoublet. If it already exists, only the xgboost pin below is
# enforced; nothing else is touched.
#
# WHY A SEPARATE ENVIRONMENT. Installing these packages from source into an existing conda R
# fails on every Bioconductor package containing C or C++ code, all for one reason: a conda R
# build's Makeconf points at a conda toolchain compiler that is not installed unless it was asked
# for, so each package stops with `x86_64-conda-linux-gnu-cc: command not found`. Dozens of
# distinct-looking failures, one cause. Adding a compiler to an environment that has already
# produced results changes that environment; taking the packages from bioconda instead, where
# they arrive already compiled, changes nothing that exists.
#
# WHY r-xgboost IS PINNED TO 1.7.6. Solved without the pin, conda-forge supplies xgboost 3.x,
# and both packages here are Bioconductor 3.18 releases written against the 1.7.x R API. The
# consequences differ, and only one of them is loud:
#
#   scds 1.18.0         Fails outright. bcds calls xgboost(mm, label = ..., nrounds = ...);
#                       `label` was removed and `y` became required, so the run stops with
#                       `argument "y" is missing, with no default`.
#
#   scDblFinder 1.16.0  Survives quietly. It passes max_depth, eta, subsample, nthread and
#                       eval_metric as top-level arguments, which xgboost 3.x accepts only
#                       through a deprecation shim that folds them into `params` and warns.
#                       It completes and returns scores.
#
# The second case is what the pin is for. Scores produced through an API-conversion shim are not
# the scores the method was characterised on, and nothing downstream can tell the difference:
# they look entirely normal. The version mismatch is visible at all only because the other package
# crashes, so the pin protects precisely the tool that would otherwise fail in silence.
#
# 1.7.6 has a cpu_r43 build, matching r-base 4.3 exactly. The verification at the end refuses the
# install if a 2.x or newer xgboost ends up in the environment anyway.
#
# Environment:
#   SCQC_CONDA   package manager to use; setup/install_env.sh passes the one it discovered

set -uo pipefail

ROOT="${1:-}"
if [ -z "$ROOT" ] || [ "$ROOT" = "-h" ] || [ "$ROOT" = "--help" ]; then
    echo "usage: install_rdoublet.sh <env-root>" >&2
    exit 2
fi
case "$ROOT" in
    /*) : ;;
    *)  echo "ERROR: <env-root> must be an absolute path (got '$ROOT')" >&2; exit 2 ;;
esac

ENVDIR="$ROOT/rdoublet"
BIN="$ROOT/bin"
mkdir -p "$BIN" || exit 1

MM="${SCQC_CONDA:-}"
if [ -z "$MM" ]; then
    for c in micromamba mamba conda; do
        if command -v "$c" >/dev/null 2>&1; then MM="$(command -v "$c")"; break; fi
    done
fi
if [ -z "$MM" ]; then
    echo "--- no conda, mamba or micromamba found; bootstrapping micromamba into $BIN ---"
    ( cd "$ROOT" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba )
    chmod +x "$BIN/micromamba" 2>/dev/null
    MM="$BIN/micromamba"
fi
[ -x "$MM" ] || { echo "ERROR: package manager '$MM' is not executable" >&2; exit 1; }
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$ROOT/mamba}"
echo "package manager: $MM  ($("$MM" --version 2>&1 | head -1))"

if [ -d "$ENVDIR" ]; then
    echo "$ENVDIR exists - enforcing the xgboost pin only"
    "$MM" install -y -p "$ENVDIR" -c conda-forge -c bioconda \
        r-xgboost=1.7.6 2>&1 | tail -20
else
    "$MM" create -y -p "$ENVDIR" \
        -c conda-forge -c bioconda \
        r-base=4.3 r-matrix r-xgboost=1.7.6 \
        bioconductor-scdblfinder bioconductor-scds \
        bioconductor-singlecellexperiment 2>&1 | tail -25
fi

echo
echo "=== verify ==="
"$ENVDIR/bin/Rscript" -e '
for (p in c("Matrix","SingleCellExperiment","scDblFinder","scds"))
  cat(sprintf("  %-24s %s  %s\n", p, requireNamespace(p, quietly=TRUE),
              tryCatch(as.character(packageVersion(p)), error=function(e) "?")))
cat(sprintf("  %-24s %s  %s\n", "xgboost", requireNamespace("xgboost", quietly=TRUE),
            as.character(packageVersion("xgboost"))))
cat("  R:", R.version.string, "\n")
v <- packageVersion("xgboost")
if (v >= "2.0.0") { cat("  FATAL: xgboost", as.character(v),
                        "breaks scds and silently shims scDblFinder\n"); quit(status=4) }
' 2>&1 | tail -12
rc=${PIPESTATUS[0]}
exit "$rc"
