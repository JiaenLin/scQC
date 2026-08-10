#!/usr/bin/env bash
# Create the environments scQC needs.
#
# Four separate environments, deliberately. They have hard, mutually incompatible pins - the
# aligner needs an older numpy than the analysis stack, and the denoiser needs its own torch
# build. One environment holding all of them resolves only by relaxing pins, and a relaxed pin
# is how a pipeline stops reproducing.
#
#   core        analysis stack: scanpy, anndata, pandas       (always)
#   celescope   CeleScope aligner                             (--with-celescope)
#   cellbender  CellBender ambient-RNA denoiser               (--with-cellbender)
#   rdoublet    R: scDblFinder doublet scoring                (--with-doublet)
#
# Usage:
#   setup/install_env.sh --prefix ~/scqc-env --all
#   setup/install_env.sh --prefix ~/scqc-env --with-cellbender --with-doublet
#
# Requires conda, mamba or micromamba; the search order is below. Nothing is installed outside
# --prefix.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
CONF="$ROOT/conf/env"

PREFIX=""
WITH_CELESCOPE=0
WITH_CELLBENDER=0
WITH_DOUBLET=0

usage() {
    # Print the header block itself rather than a copy of it: a fixed line range goes stale the
    # first time the header is edited, and a usage message that has drifted from the code is
    # indistinguishable from one that has not.
    awk 'NR == 1 { next } /^#/ { sub(/^# ?/, ""); print; next } { exit }' "${BASH_SOURCE[0]}"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)          PREFIX="${2:?--prefix needs a path}"; shift 2 ;;
        --with-celescope)  WITH_CELESCOPE=1; shift ;;
        --with-cellbender) WITH_CELLBENDER=1; shift ;;
        --with-doublet)    WITH_DOUBLET=1; shift ;;
        --all)             WITH_CELESCOPE=1; WITH_CELLBENDER=1; WITH_DOUBLET=1; shift ;;
        -h|--help)         usage 0 ;;
        *) echo "unknown argument: $1" >&2; usage 1 ;;
    esac
done

[[ -n "$PREFIX" ]] || { echo "ERROR: --prefix is required" >&2; usage 1; }

# --- find a package manager --------------------------------------------------------------------
# Three places, in order of preference. The third matters: on a module-based HPC conda is
# installed but deliberately absent from the default PATH, so a PATH-only check reports "not
# installed" on a machine where it is one `module load` away. That is a false negative on exactly
# the class of system this pipeline is meant to run on.
CONDA=""
CONDA_VIA="PATH"

find_on_path() {
    # Resolve to an ABSOLUTE PATH, never a bare command name. The helper scripts receive this
    # through SCQC_CONDA and guard it with `[ -x "$MM" ]`; `-x conda` tests the relative path
    # ./conda, which does not exist, so a bare name made every optional component abort with
    # "package manager 'conda' is not executable" on a machine where conda WAS on PATH.
    local c p
    for c in mamba micromamba conda; do
        p="$(command -v "$c" 2>/dev/null)" || continue
        [[ -n "$p" && -x "$p" ]] && { CONDA="$p"; return 0; }
    done
    return 1
}

find_installed() {
    local p
    for p in "$HOME"/{miniforge3,mambaforge,miniconda3,anaconda3}/bin \
             /opt/{miniforge3,miniconda3,conda,anaconda3}/bin \
             /apps/{anaconda3,miniconda3}/bin \
             /usr/local/{miniconda3,anaconda3}/bin; do
        for c in mamba conda; do
            if [[ -x "$p/$c" ]]; then
                CONDA="$p/$c"; CONDA_VIA="$p"; return 0
            fi
        done
    done
    return 1
}

find_via_module() {
    # `module` is a shell function, so it must be sourced into THIS shell - and the load must not
    # be piped or run in a subshell, or the PATH change is discarded with the subshell.
    local init m
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash; do
        [[ -r "$init" ]] && { . "$init"; break; }
    done
    command -v module >/dev/null 2>&1 || type module >/dev/null 2>&1 || return 1
    for m in $(module avail 2>&1 | tr ' ' '\n' \
               | grep -iE '^(anaconda|miniconda|miniforge|conda|mamba)' | sort -u); do
        if module load "$m" >/dev/null 2>&1 && find_on_path; then
            CONDA_VIA="module load $m"
            echo "note: loaded environment module '$m' to find conda" >&2
            return 0
        fi
    done
    return 1
}

find_on_path || find_installed || find_via_module || {
    cat >&2 <<'ERR'
ERROR: no conda, mamba or micromamba could be found.

Searched:
  1. $PATH
  2. common install prefixes (~/miniforge3, /opt/miniconda3, /apps/anaconda3, ...)
  3. environment modules matching anaconda/miniconda/miniforge/conda/mamba

On an HPC, conda is often present but not on the default PATH. Try:
  module avail 2>&1 | grep -i conda
  module load <the module it lists>

Otherwise install Miniforge: https://github.com/conda-forge/miniforge
ERR
    exit 1
}

case "$PREFIX" in
    /*|~*) : ;;
    *) echo "ERROR: --prefix must be an absolute path (got '$PREFIX')" >&2; exit 1 ;;
esac
PREFIX="${PREFIX/#\~/$HOME}"

# Environments are NOT relocatable: conda writes the absolute prefix into script shebangs and
# into compiled shared-library paths. Moving one afterwards produces an environment whose
# executables silently resolve to the wrong interpreter, so the path is fixed at creation.
if [[ "$PREFIX" == *" "* ]]; then
    echo "ERROR: --prefix must not contain spaces; conda shebangs break." >&2
    exit 1
fi

echo "installer : $CONDA  (found via: $CONDA_VIA)"
echo "version   : $("$CONDA" --version 2>&1 | head -1)"
echo "prefix    : $PREFIX"
mkdir -p "$PREFIX"

created=()

make_env() {
    local name="$1" py="$2" path="$PREFIX/$1"
    if [[ -d "$path" ]]; then
        echo "  [skip] $name already exists at $path"
        return 0
    fi
    echo "  [make] $name (python $py)"
    "$CONDA" create -y -p "$path" -c conda-forge "python=$py" >/dev/null
    created+=("$name")
}

# The interpreter version a lock file was resolved against, read FROM the lock file.
#
# It was hardcoded below as 3.11 while the lock pinned a package requiring >= 3.12, and nothing
# compared the two. The install died on the first package with a wall of resolver output listing
# every anndata ever released; `set -e` then aborted before any other component was built, and
# what remained was an environment that existed as a directory and contained nothing. A lock file
# is only valid for the interpreter it was resolved against, so that interpreter belongs in it.
lock_python() {
    local f="$1" v
    v="$(sed -n 's/^#[[:space:]]*scqc-python:[[:space:]]*\([0-9][0-9.]*\).*/\1/p' "$f" | head -1)"
    if [[ -z "$v" ]]; then
        echo "ERROR: $f carries no '# scqc-python: X.Y' line." >&2
        echo "       A lock file without the interpreter it was resolved against cannot be" >&2
        echo "       installed safely - that version is part of the lock." >&2
        exit 1
    fi
    printf '%s' "$v"
}

# --- core --------------------------------------------------------------------------------------
echo
echo "== core =="
CORE_PY="$(lock_python "$CONF/requirements.lock.txt")"
make_env core "$CORE_PY"
# Verified BEFORE installing. Conda resolves to the nearest satisfiable build, so the environment
# that exists is not necessarily the one requested, and the mismatch would otherwise surface as
# an unreadable dependency error rather than as the version conflict it is.
GOT_PY="$("$PREFIX/core/bin/python" -c 'import sys; print("%d.%d" % sys.version_info[:2])')"
if [[ "$GOT_PY" != "$CORE_PY" ]]; then
    echo "ERROR: the core environment has python $GOT_PY, but" >&2
    echo "       conf/env/requirements.lock.txt was resolved against $CORE_PY." >&2
    echo "       Installing it here would fail outright, or silently install versions other" >&2
    echo "       than the ones locked. Remove $PREFIX/core and re-run, or re-lock for $GOT_PY." >&2
    exit 1
fi
echo "  python $GOT_PY, matching the lock"
"$PREFIX/core/bin/pip" install -q --no-input -r "$CONF/requirements.lock.txt"
echo "  installed from conf/env/requirements.lock.txt"

# --- optional components -------------------------------------------------------------------------
if [[ $WITH_CELESCOPE -eq 1 ]]; then
    echo
    echo "== celescope =="
    echo "  see conf/env/install_cs_cb.sh - it pins the aligner and its reference tooling"
    SCQC_CONDA="$CONDA" bash "$CONF/install_cs_cb.sh" "$PREFIX" celescope || {
        echo "  WARNING: celescope install failed; core is still usable" >&2
    }
fi

if [[ $WITH_CELLBENDER -eq 1 ]]; then
    echo
    echo "== cellbender =="
    make_env cellbender 3.10
    "$PREFIX/cellbender/bin/pip" install -q --no-input cellbender==0.3.2
    echo "  cellbender 0.3.2"
fi

if [[ $WITH_DOUBLET -eq 1 ]]; then
    echo
    echo "== rdoublet =="
    # install_rdoublet.sh, not install_scdblfinder.sh: the latter builds from source, which needs
    # a compiler the conda R does not ship, and it cannot pin xgboost - and an unpinned xgboost
    # returns scDblFinder scores that look normal but were computed through a deprecation shim.
    SCQC_CONDA="$CONDA" bash "$CONF/install_rdoublet.sh" "$PREFIX" || {
        echo "  WARNING: scDblFinder install failed; supply your own detector instead" >&2
    }
fi

# --- verify: report what is actually importable, not what was requested --------------------------
echo
echo "== verification =="
ok=1
check() {
    local label="$1" py="$2"; shift 2
    if [[ ! -x "$py" ]]; then
        printf '  %-12s MISSING\n' "$label"; ok=0; return
    fi
    # Joined with a subshell IFS, NOT with "${*// /,}".
    #
    # That form substitutes into EACH positional parameter and then joins the results with a
    # space - it never joins first and substitutes after. So it produced
    # `import scanpy anndata pandas numpy`, a SyntaxError, and every environment with more than
    # one package to check reported BROKEN however healthy it was. `core` is the only such
    # environment here, so the verification step has never once passed. The single-package envs
    # verified correctly, which is why it looked selective rather than wrong.
    local mods err
    mods="$(IFS=,; printf '%s' "$*")"
    # stderr CAPTURED, not discarded. `2>/dev/null` turned a real diagnosis into the words
    # "imports failed", and a failure report that omits the failure sends the reader to the wrong
    # place - here, to a perfectly good environment.
    if err="$("$py" -c "import $mods" 2>&1)"; then
        printf '  %-12s ok   %s\n' "$label" "$("$py" --version 2>&1)"
    else
        printf '  %-12s BROKEN: %s\n' "$label" "$(printf '%s' "$err" | tail -1)"
        printf '  %-12s tried: import %s\n' "" "$mods"
        ok=0
    fi
}
check core       "$PREFIX/core/bin/python"       scanpy anndata pandas numpy
[[ $WITH_CELLBENDER -eq 1 ]] && check cellbender "$PREFIX/cellbender/bin/python" cellbender
[[ $WITH_CELESCOPE  -eq 1 ]] && check celescope  "$PREFIX/celescope/bin/python"  celescope

# rdoublet is R, so the python import check above does not apply to it. Load the packages that
# are ACTUALLY USED - the same list adapters/scdblfinder.R declares in NEEDED. A present Rscript
# with an absent scDblFinder is the failure mode this catches.
#
# scds was in this list until 2026-08-10 and nothing in the pipeline has ever called it. It made
# a missing package that no step needs report the environment BROKEN, which is the same shape as
# the multi-package import bug fixed the same day: a guard failing for a reason unrelated to
# whether the run can proceed.
if [[ $WITH_DOUBLET -eq 1 ]]; then
    rs="$PREFIX/rdoublet/bin/Rscript"
    if [[ ! -x "$rs" ]]; then
        printf '  %-12s MISSING\n' rdoublet; ok=0
    elif "$rs" -e 'q(status = if (all(sapply(c("scDblFinder"), requireNamespace, quietly=TRUE))) 0 else 1)' >/dev/null 2>&1; then
        printf '  %-12s ok   %s\n' rdoublet "$("$rs" --version 2>&1 | head -1)"
    else
        printf '  %-12s BROKEN (scDblFinder did not load)\n' rdoublet; ok=0
    fi
fi

cat <<EOF

Environments created: ${#created[@]}
Point scQC at them with:

  export SCQC_ENV_ROOT="$PREFIX"

Add that line to your shell profile, or pass --env-root on every run.
EOF

[[ $ok -eq 1 ]] || {
    echo
    echo "One or more environments did not verify. Read the lines marked BROKEN or MISSING above" >&2
    echo "before running the pipeline - a half-installed environment fails deep in a run, not at" >&2
    echo "the start." >&2
    exit 1
}
