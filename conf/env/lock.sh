# Reading, checking and applying a LOCK. Sourced by setup/install_env.sh and conf/env/install_cellbender.sh;
# never run on its own.
#
# A lock is a pip requirements file holding every package of one environment's closure at one exact
# version. It names the interpreter it was resolved against on a line of its own, `# scqc-python: X.Y.Z`,
# patch included. A line may name an index (`--index-url URL`, `--extra-index-url URL`): where a wheel
# comes from is part of what it is, as torch's CUDA build shows. Every other line is `name==version`,
# optionally followed by an environment marker.
#
# A range, a bare name, a URL or an editable is not a pin. pip would resolve it on the day of the build,
# and the environment would be whatever the index held that day - which is how a denoiser environment
# acquired a torch nobody had chosen. So a lock is applied with `pip install --no-deps`, and pip resolves
# nothing; every installed version is then compared with the lock; and a lock is checked before a package
# manager is looked for, so a refusal costs nothing.

# The interpreter a lock was resolved against, read FROM the lock file, as major.minor.patch.
#
# It was once hardcoded in the installer as 3.11 while the lock pinned a package requiring >= 3.12, and
# nothing compared the two. The patch is required as well: an analysis call record writes the
# interpreter's full version down (`versions.python`), so two environments built from one lock on
# different days must not differ in it.
lock_python() {
    local f="$1" v
    v="$(sed -n 's/^#[[:space:]]*scqc-python:[[:space:]]*\([0-9][0-9.]*\).*/\1/p' "$f" | head -1)"
    if [[ -z "$v" ]]; then
        echo "ERROR: $f carries no '# scqc-python: X.Y.Z' line." >&2
        echo "       A lock file without the interpreter it was resolved against cannot be" >&2
        echo "       installed safely - that version is part of the lock." >&2
        exit 1
    fi
    case "$v" in
        [0-9]*.[0-9]*.[0-9]*) : ;;
        *)  echo "ERROR: $f names python $v. A lock names the interpreter's patch as well (X.Y.Z):" >&2
            echo "       two builds from one lock must hold the same interpreter." >&2
            exit 1 ;;
    esac
    printf '%s' "$v"
}

# The full version of an environment's interpreter, X.Y.Z.
env_python() {
    "$1" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])'
}

# `name version` for every exact pin of a lock, the name normalised as pip compares names: lower case,
# every run of - _ . as one -.
lock_pins() {
    awk '
        { sub(/\r$/, "") }
        /^[[:space:]]*(#|$)/ { next }
        /^[[:space:]]*-/ { next }
        {
            line = $0
            sub(/[[:space:]]*;.*/, "", line)
            gsub(/[[:space:]]/, "", line)
            if (split(line, part, "==") == 2) {
                name = tolower(part[1]); gsub(/[-_.]+/, "-", name)
                print name, part[2]
            }
        }
    ' "$1"
}

# The requirement line(s) of a lock that name a package, whatever form they take.
lock_line() {
    awk -v want="$2" '
        { sub(/\r$/, "") }
        /^[[:space:]]*(#|$|-)/ { next }
        {
            if (match($0, /^[[:space:]]*[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                name = tolower(substr($0, RSTART, RLENGTH))
                gsub(/[[:space:]]/, "", name); gsub(/[-_.]+/, "-", name)
                if (name == want) print $0
            }
        }
    ' "$1"
}

# Refuse a lock that is not exact; otherwise print one line: what it is, how many pins, which interpreter.
check_lock() {
    local f="$1" label="$2" bad dup n
    if [[ ! -s "$f" ]]; then
        echo "ERROR: $label: there is no lock at $f." >&2
        return 1
    fi
    lock_python "$f" >/dev/null
    bad="$(awk '
        { sub(/\r$/, "") }
        /^[[:space:]]*(#|$)/ { next }
        /^[[:space:]]*--(extra-)?index-url([[:space:]]+|=)[^[:space:]]+[[:space:]]*$/ { next }
        {
            line = $0
            sub(/[[:space:]]*;.*/, "", line)
            if (line !~ /^[[:space:]]*[A-Za-z0-9][A-Za-z0-9._-]*[[:space:]]*==[[:space:]]*[A-Za-z0-9][A-Za-z0-9.+!_-]*[[:space:]]*$/)
                printf "         line %d: %s\n", NR, $0
        }
    ' "$f")"
    if [[ -n "$bad" ]]; then
        echo "ERROR: $label: $f holds requirements that are not one package at one exact version:" >&2
        printf '%s\n' "$bad" >&2
        echo "       pip would resolve them on the day of the build, so the environment would be" >&2
        echo "       whatever the index held that day. Pin each one (name==version)." >&2
        return 1
    fi
    dup="$(lock_pins "$f" | awk '{ c[$1]++ } END { for (n in c) if (c[n] > 1) printf " %s", n }')"
    if [[ -n "$dup" ]]; then
        echo "ERROR: $label: $f pins these more than once:$dup" >&2
        return 1
    fi
    n="$(lock_pins "$f" | wc -l | tr -d ' ')"
    echo "  $label: $n exact pins, python $(lock_python "$f"), from $f"
}

# The lock pins each `name==version` given - whatever else it says. A version is compared without its
# local label (`+cu121` names the build of a wheel, not its version), and nothing else is accepted: a
# lock is regenerated by a resolver, these are decisions, and a lock that disagrees with one is refused.
require_pins() {
    local f="$1" label="$2" req name want line got bad=0
    shift 2
    for req in "$@"; do
        name="${req%%==*}"; want="${req#*==}"
        line="$(lock_line "$f" "$name" | head -1)"
        got="$(lock_pins "$f" | awk -v n="$name" '$1 == n { print $2; exit }')"
        if [[ -z "$line" ]]; then
            echo "ERROR: $label: the lock does not name $name. It must pin $name $want." >&2
            bad=1
        elif [[ -z "$got" ]]; then
            echo "ERROR: $label: the lock names $name as '$line', which is not one version." >&2
            echo "       It must pin $name $want." >&2
            bad=1
        elif [[ "${got%%+*}" != "$want" ]]; then
            echo "ERROR: $label: the lock pins $name $got. It must pin $name $want; no other is accepted." >&2
            bad=1
        fi
    done
    return "$bad"
}

# Apply a lock in an environment's interpreter, resolving nothing.
apply_lock() {
    "$1" -m pip install --no-deps --no-input --disable-pip-version-check -q -r "$2"
}

# Compare every version a lock pins with the version installed, in the environment's own interpreter.
# One line per difference and a summary; a difference fails.
verify_lock() {
    "$1" - "$2" "$3" <<'PY'
import sys
from importlib import metadata

lock, label = sys.argv[1], sys.argv[2]
try:
    from packaging.markers import Marker
    from packaging.version import Version
except Exception:                       # an interpreter without packaging compares the text
    Marker = Version = None


def same(got, want):
    if got == want:
        return True
    try:
        return Version is not None and Version(got) == Version(want)
    except Exception:
        return False


pins = []
for raw in open(lock, encoding="utf-8"):
    line = raw.strip()
    if not line or line.startswith(("#", "-")):
        continue
    spec, _, marker = line.partition(";")
    if marker.strip() and Marker is not None and not Marker(marker.strip()).evaluate():
        continue
    name, _, version = spec.strip().partition("==")
    pins.append((name.strip(), version.strip()))

wrong = []
for name, want in pins:
    try:
        got = metadata.version(name)
    except metadata.PackageNotFoundError:
        got = None
    if got is None or not same(got, want):
        wrong.append(name)
        print(f"  {label}: {name} {got or 'is not installed'}{' is installed' if got else ''}; "
              f"the lock pins {name}=={want}")
if wrong:
    print(f"  {label}: {len(wrong)} of the {len(pins)} locked versions are not the ones installed")
    sys.exit(1)
print(f"  {label}: every one of the {len(pins)} locked versions is the one installed")
PY
}
