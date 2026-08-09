#!/bin/bash
# SUPERSEDED: install scDblFinder and scds from source into a shared R library tree.
#
# Use conf/env/install_rdoublet.sh instead. That is the route setup/install_env.sh takes and the
# one the pipeline is tested against. This script is kept only for the case where an existing R
# must be reused and a second conda environment is not an option.
#
# Two reasons it is not the supported route:
#
#   1. It builds from source. A conda-provided R has no compiler unless one was installed with
#      it, and every Bioconductor package containing C or C++ code then stops with
#      `x86_64-conda-linux-gnu-cc: command not found`.
#
#   2. It cannot pin xgboost. BiocManager resolves xgboost from CRAN, which supplies 3.x. Against
#      that version scds fails outright, while scDblFinder keeps running through a deprecation
#      shim and returns scores that look entirely normal but were not computed the way the method
#      was characterised. install_rdoublet.sh pins r-xgboost=1.7.6 for exactly this reason. Here
#      the version can only be checked after the fact, which the final block does: it exits
#      non-zero rather than leave a silently shimmed install in place.
#
# Usage:
#   conf/env/install_scdblfinder.sh <env-root> [Rscript]
#
# Packages are installed into <env-root>/Rlib. Rscript defaults to
# <env-root>/celescope/bin/Rscript, then to Rscript on PATH. No Python environment is modified.
#
# Environment:
#   BIOC_VERSION   Bioconductor release to install (default 3.19, which matches R 4.4.x)

set -uo pipefail

ROOT="${1:-}"
if [ -z "$ROOT" ] || [ "$ROOT" = "-h" ] || [ "$ROOT" = "--help" ]; then
    echo "usage: install_scdblfinder.sh <env-root> [Rscript]" >&2
    echo "       (superseded; prefer conf/env/install_rdoublet.sh)" >&2
    exit 2
fi
case "$ROOT" in
    /*) : ;;
    *)  echo "ERROR: <env-root> must be an absolute path (got '$ROOT')" >&2; exit 2 ;;
esac

RS="${2:-}"
if [ -z "$RS" ]; then
    if [ -x "$ROOT/celescope/bin/Rscript" ]; then
        RS="$ROOT/celescope/bin/Rscript"
    else
        RS="$(command -v Rscript 2>/dev/null || true)"
    fi
fi
[ -n "$RS" ] && [ -x "$RS" ] || {
    echo "ERROR: no Rscript found. Pass one as the second argument." >&2
    exit 1
}

export BIOC_VERSION="${BIOC_VERSION:-3.19}"
export R_LIBS_USER="$ROOT/Rlib"
mkdir -p "$R_LIBS_USER" || exit 1

echo "SUPERSEDED script - conf/env/install_rdoublet.sh is the supported route."
echo "Rscript      : $RS"
echo "library      : $R_LIBS_USER"
echo "Bioconductor : $BIOC_VERSION"

"$RS" -e '
.libPaths(Sys.getenv("R_LIBS_USER"))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
bioc <- Sys.getenv("BIOC_VERSION", "3.19")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = bioc, ask = FALSE, update = FALSE)
cat("BiocManager version:", as.character(BiocManager::version()), "\n")
BiocManager::install(c("scDblFinder", "scds"), ask = FALSE, update = FALSE, Ncpus = 8)
for (p in c("scDblFinder", "scds", "SingleCellExperiment", "Matrix")) {
  cat(sprintf("  %-22s %s\n", p, requireNamespace(p, quietly = TRUE)))
}
if (!requireNamespace("xgboost", quietly = TRUE)) {
  cat("  FATAL: xgboost is not installed; scds cannot run\n"); quit(status = 4)
}
v <- packageVersion("xgboost")
cat(sprintf("  %-22s %s\n", "xgboost", as.character(v)))
if (v >= "2.0.0") {
  cat("  FATAL: xgboost", as.character(v), "breaks scds and silently shims scDblFinder.\n")
  cat("  Use conf/env/install_rdoublet.sh, which pins r-xgboost=1.7.6.\n")
  quit(status = 4)
}
' 2>&1 | tail -40
rc=${PIPESTATUS[0]}
exit "$rc"
