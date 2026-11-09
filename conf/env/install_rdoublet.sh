#!/bin/bash
# Install R with scDblFinder as a separate conda environment, from precompiled packages.
#
# Usage:
#   conf/env/install_rdoublet.sh <env-root>
#   conf/env/install_rdoublet.sh --check      print the pinned specification; build nothing
#
# <env-root> is the directory given to setup/install_env.sh --prefix. The environment name is
# fixed by the pipeline: <env-root>/rdoublet. If it already exists it is checked against the pins
# below and never modified: an environment something may be using is not changed underneath it,
# and one holding another version is refused rather than reused.
#
# THE PINS are the versions the environment that scored the calibration cohort's doublets reported in
# that run's own report (its provenance block, "observed"): R 4.3.3, scDblFinder 1.16.0,
# SingleCellExperiment 1.24.0, xgboost 1.7.6.1 and Matrix 1.6.5. scDblFinder 1.16.0 is a Bioconductor
# 3.18 release, and bioconda's build of it pins every Bioconductor dependency to that release, which
# admits SingleCellExperiment 1.24.0 and no other.
#
# A PACKAGE HAS TWO VERSION NUMBERS HERE, conda's and R's, and they are not always the same.
# conda-forge's r-xgboost 1.7.6 installs the CRAN release of the R package, which R reports as 1.7.6.1;
# conda's r-matrix 1.6_5 is R's Matrix 1.6.5. So the conda specification names conda's numbers, and the
# check after the build reads R's: what R reports is what scDblFinder runs against. The first version of
# this check compared R's report with conda's number and refused every correct environment (harness
# autodev W014's verifier). The CPU build of r-xgboost is read from the environment's conda record.
#
# Everything else the solver picks - the other CRAN packages, the BLAS - is NOT pinned: the conda half
# of this environment is not locked, and a build months later can differ there. The versions R reports
# are printed at the end; read them before comparing two builds.
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
# reads packageVersion("xgboost") directly - 1.7.6.1, the R package conda-forge's 1.7.6 installs - and
# refuses any other, naming the shim when it is >= 2. scds was removed 2026-08-10; this note stays
# because the reasoning behind the pin is worth keeping even once the package that revealed it is gone.
#
# 1.7.6 has a cpu_r43 build, matching r-base 4.3 exactly, and the CPU build is what is asked for: a
# GPU build chosen on a GPU node would be another scorer. The verification at the end refuses the
# install if any pinned version is not the one R reports, or the conda record names another build.
#
# Environment:
#   SCQC_CONDA   package manager to use; setup/install_env.sh passes the one it discovered

set -uo pipefail

# conda's numbers, for the specification
R_VERSION=4.3.3
SCDBLFINDER_VERSION=1.16.0          # Bioconductor 3.18
SCE_VERSION=1.24.0                  # the only one Bioconductor 3.18 admits
XGBOOST_VERSION=1.7.6               # conda-forge's r-xgboost, a cpu_ build ...
MATRIX_VERSION=1.6_5                # conda-forge's r-matrix
# R's numbers, for the check: what `packageVersion()` reports once those are installed
XGBOOST_R_VERSION=1.7.6.1           # ... whose R package is the CRAN release 1.7.6.1
MATRIX_R_VERSION=1.6.5
SPEC=(r-base=$R_VERSION "r-matrix=$MATRIX_VERSION" "r-xgboost=$XGBOOST_VERSION=cpu_*"
      bioconductor-scdblfinder=$SCDBLFINDER_VERSION
      bioconductor-singlecellexperiment=$SCE_VERSION)

ROOT="${1:-}"
if [ "$ROOT" = "--check" ]; then
    echo "  rdoublet: conda-forge then bioconda: ${SPEC[*]}"
    echo "    read back after the build, as R reports them: R $R_VERSION, scDblFinder $SCDBLFINDER_VERSION," \
         "SingleCellExperiment $SCE_VERSION, xgboost $XGBOOST_R_VERSION (conda's $XGBOOST_VERSION)," \
         "Matrix $MATRIX_R_VERSION (conda's $MATRIX_VERSION); and r-xgboost's conda build a cpu_ one." \
         "The other packages are not pinned."
    exit 0
fi
if [ -z "$ROOT" ] || [ "$ROOT" = "-h" ] || [ "$ROOT" = "--help" ]; then
    echo "usage: install_rdoublet.sh <env-root> | --check" >&2
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

if [ -d "$ENVDIR" ] && [ ! -x "$ENVDIR/bin/Rscript" ]; then
    # EXISTING IS NOT COMPLETE (single-cell-harness ADR-0027): a build killed mid-transaction
    # leaves the directory without its interpreter, and "exists" then skipped the build for ever.
    echo "ERROR: $ENVDIR exists but has no bin/Rscript - the remains of an interrupted build." >&2
    echo "       Remove it and run again:  rm -rf '$ENVDIR'" >&2
    exit 1
fi
made=0
if [ -d "$ENVDIR" ]; then
    echo "$ENVDIR exists - checked against the pins below, not modified"
else
    # conda-forge and bioconda alone: a channel list from an rc file must not decide what is solved.
    "$MM" create -y -p "$ENVDIR" --override-channels \
        -c conda-forge -c bioconda "${SPEC[@]}" 2>&1 | tail -25
    [ -x "$ENVDIR/bin/Rscript" ] \
        || { echo "ERROR: the package manager did not create $ENVDIR with an Rscript (above)" >&2; exit 1; }
    made=1
fi

echo
echo "=== verify ==="
# The package list here is exactly what adapters/scdblfinder.R declares in NEEDED. Verifying
# anything else makes the environment report BROKEN for a package the pipeline never calls, and
# a guard that fails for reasons unrelated to whether the run can proceed is a guard that gets
# switched off.
# R REPORTS what it holds, one `name<TAB>version` line each, and the comparison is made here, in the
# shell: it then reads the same whoever reports, and R's numbers are compared with R's.
REPORT="$("$ENVDIR/bin/Rscript" -e '
ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
cat(sprintf("R\t%s.%s\n", R.version$major, R.version$minor))
for (p in c("scDblFinder", "SingleCellExperiment", "xgboost", "Matrix", "scran", "scuttle",
            "BiocNeighbors", "BiocSingular", "BiocParallel"))
  cat(sprintf("%s\t%s\n", p, ver(p)))
' 2>&1)"
rc=$?
printf '%s\n' "$REPORT" | awk -F'\t' 'NF == 2 { printf "  %-24s %s\n", $1, $2; next } { print "  " $0 }'
if [ "$rc" -ne 0 ]; then
    echo "ERROR: $ENVDIR/bin/Rscript could not report its packages (exit $rc, above)" >&2
    exit 1
fi
reported() { printf '%s\n' "$REPORT" | awk -F'\t' -v k="$1" '$1 == k { print $2; exit }'; }
bad=""
pinned() {
    local got
    got="$(reported "$1")"
    [ "$got" = "$2" ] || bad="$bad; $1 ${got:-not reported} (pinned $2)"
}
pinned R "$R_VERSION"
pinned scDblFinder "$SCDBLFINDER_VERSION"
pinned SingleCellExperiment "$SCE_VERSION"
pinned xgboost "$XGBOOST_R_VERSION"
pinned Matrix "$MATRIX_R_VERSION"
x="$(reported xgboost)"
case "${x%%.*}" in
    ''|*[!0-9]*) : ;;
    *) [ "${x%%.*}" -lt 2 ] || echo "  FATAL: xgboost $x silently shims scDblFinder onto a converted API" ;;
esac
# The build, from the environment's own conda record: <name>-<version>-<build>.json
rec="$(ls "$ENVDIR"/conda-meta/r-xgboost-*.json 2>/dev/null | head -1)"
if [ -z "$rec" ]; then
    bad="$bad; r-xgboost has no conda record in $ENVDIR/conda-meta"
else
    rec="${rec##*/}"; rec="${rec%.json}"; rec="${rec#r-xgboost-}"
    [ "${rec%%-*}" = "$XGBOOST_VERSION" ] || bad="$bad; r-xgboost ${rec%%-*} in conda (pinned $XGBOOST_VERSION)"
    case "${rec#*-}" in
        cpu_*) echo "  r-xgboost build           ${rec#*-}" ;;
        *) bad="$bad; r-xgboost build ${rec#*-} (pinned a cpu_ build)" ;;
    esac
fi
if [ -n "$bad" ]; then
    echo "  REFUSED - not the pinned versions:${bad#;}"
    if [ "$made" -eq 0 ]; then
        echo "  It existed before this run and was not modified. Remove it and re-run, or build into a" \
             "new prefix." >&2
    else
        echo "  The solver did not give the pinned versions. Remove $ENVDIR and re-run; if it gives them" \
             "again, read the specification above against conda-forge and bioconda." >&2
    fi
    exit 5
fi
echo "  every pinned version is the one installed"
