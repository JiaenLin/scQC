#!/bin/bash
# Install R with scDblFinder as a separate conda environment, from precompiled packages.
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
# while scDblFinder 1.16.0 is a Bioconductor 3.18 release written against the 1.7.x R API. It
# does not fail. It passes max_depth, eta, subsample, nthread and eval_metric as top-level
# arguments, which xgboost 3.x accepts only through a deprecation shim that folds them into
# `params` and warns - so it completes and returns scores.
#
# That is the whole reason for the pin. Scores produced through an API-conversion shim are not
# the scores the method was characterised on, and nothing downstream can tell the difference:
# they look entirely normal.
#
# HOW THIS WAS FOUND, and why the finder is no longer installed. The mismatch surfaced because a
# sibling package, scds 1.18.0, crashed outright on the same xgboost - `bcds` calls
# xgboost(label = ...), `label` was removed and `y` became required. scds was then kept in this
# environment as though it were a canary. It was not one: nothing in the pipeline calls it, and
# the verification only did `requireNamespace("scds")`, which loads the package without ever
# touching the xgboost API. A canary that is never asked to sing is just weight - and a required
# package the pipeline never invokes turns an unrelated breakage into a BROKEN environment.
#
# The guard that actually works is the explicit version check at the end of this script, which
# reads packageVersion("xgboost") directly and refuses >= 2.0.0. scds was removed 2026-08-10;
# this note stays because the reasoning behind the pin is worth keeping even once the package
# that revealed it is gone.
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
        bioconductor-scdblfinder \
        bioconductor-singlecellexperiment 2>&1 | tail -25
fi

echo
echo "=== verify ==="
# The package list here is exactly what adapters/scdblfinder.R declares in NEEDED. Verifying
# anything else makes the environment report BROKEN for a package the pipeline never calls, and
# a guard that fails for reasons unrelated to whether the run can proceed is a guard that gets
# switched off.
"$ENVDIR/bin/Rscript" -e '
for (p in c("Matrix","SingleCellExperiment","scDblFinder"))
  cat(sprintf("  %-24s %s  %s\n", p, requireNamespace(p, quietly=TRUE),
              tryCatch(as.character(packageVersion(p)), error=function(e) "?")))
cat(sprintf("  %-24s %s  %s\n", "xgboost", requireNamespace("xgboost", quietly=TRUE),
            as.character(packageVersion("xgboost"))))
cat("  R:", R.version.string, "\n")
v <- packageVersion("xgboost")
if (v >= "2.0.0") { cat("  FATAL: xgboost", as.character(v),
                        "silently shims scDblFinder onto a converted API\n"); quit(status=4) }
' 2>&1 | tail -12
rc=${PIPESTATUS[0]}
exit "$rc"
