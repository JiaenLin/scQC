# Execution adapter: FASTQ -> raw unfiltered count matrix, via CeleScope multi_rna.
# It removes no observation. STARsolo's cell call is a by-product of this step and is never this
# pipeline's cell call; what the step delivers downstream is the RAW droplet matrix.
"""CeleScope adapter - run the aligner, then hand back the unfiltered matrix or refuse.

WHY THIS ADAPTER EXISTS AT ALL

The pipeline needs a matrix that still contains its empty droplets. `lib/verify_raw.py` exists
because a delivered matrix usually does not: it is the aligner's `outs/filtered`, cell-called,
and nothing in the file says so. When step 0 rejects a supplied matrix the only remedy is to
regenerate it from FASTQ, and that is what this module does.

TWO PHASES, BECAUSE `multi_rna` DOES NOT RUN THE PIPELINE

`multi_rna --mod shell` is a PLANNER. It reads the mapfile and writes `shell/<sample>.sh` into
the working directory; it aligns nothing. A caller that runs `multi_rna` and then looks for a
matrix finds an empty tree and a zero exit status - the exact shape of failure `engine/task.py`
is written against. So `run_celescope` runs the planner, locates the script it wrote, refuses if
it is absent, and then executes it as a second command with its own log.

THREE HAZARDS THAT PRODUCE WRONG OUTPUT RATHER THAN AN ERROR

  1. `--chemistry` is per SAMPLE and must be passed explicitly. Chemistry has been observed to
     differ inside a single cohort - one version for some libraries and an earlier one for the
     rest - so a run-level value is wrong for part of the cohort without ever saying so. A
     library processed under the wrong chemistry yields near-empty output because the barcode
     and UMI positions do not match the read structure, and the run still exits zero. `auto` is
     refused for the same reason a default is: whatever it detects is not the value the report
     will be able to state was requested. See `require_chemistry`.

  2. CeleScope invokes STAR by BARE NAME through a shell, with no check that it was found. If
     the environment's `bin` is not on `PATH` the job dies inside the generated script with a
     missing-file message that names neither STAR nor the environment. The environment bin is
     therefore prepended to `PATH` in the env dict handed to the executor, and STAR is resolved
     on that constructed `PATH` BEFORE the aligner is launched, so the failure arrives in
     seconds instead of after a scheduler allocation.

  3. A PREVIOUS RUN'S OUTPUT SATISFIES EVERY CHECK THIS ADAPTER MAKES. Locating a matrix after
     the command returns proves that a matrix is there, not that this invocation wrote it: a run
     that exits zero having written nothing passes every one of those checks whenever leftovers
     are lying at the paths they look at, and the previous run's numbers are then recorded under
     this run's parameters. `run_celescope` therefore REFUSES to start when `<work_dir>/<sample>`
     or `shell/<sample>.sh` already exists - the same refusal the Cell Ranger adapter makes for an
     existing pipestance, and for the stronger reason that phase two EXECUTES `shell/<sample>.sh`:
     a leftover script runs the earlier run's chemistry and reference under this run's recorded
     parameters. Because the tree provably did not exist beforehand, anything found afterwards
     was written by this invocation, and that - not the mtime - is the proof. The mtimes are
     recorded as corroboration only: they come from the compute node's clock while the refusal is
     made on the orchestrator's, and a gate that fires on clock skew is a gate that gets removed.

PATHS ARE NOT ROUND-TRIPPED

Anything that reaches the command line, the mapfile or the executor's `cwd` is the TEXT the
caller supplied. `str(Path("/refs/mm39"))` is `\\refs\\mm39` on a Windows orchestrator, so every
`Path()` round trip is a silent rewrite into a path no POSIX cluster can open, arriving as an
error about a missing reference. A RELATIVE path cannot be left alone - it would resolve against
the work directory rather than the caller's - so it is resolved against this host's working
directory and the fact is recorded in `metrics["path_resolution"]`, per path, rather than done
quietly. `Path` is still used freely for local existence checks; what it must not do is decide
the spelling of an argument.

WHAT IS REPRODUCED, EXACTLY

The reference cohort was processed with:

    multi_rna --mapfile MAP --genomeDir REF --chemistry CHEM \\
              --soloFeatures "Gene GeneFull_Ex50pAS Velocyto" \\
              --report_soloFeature GeneFull_Ex50pAS \\
              --soloCellFilter EmptyDrops_CR --thread N --mod shell

Every flag is preserved and each is a module constant. `--soloFeatures` is ONE argv element
containing spaces, not three elements: the quoting in the line above is the shell's, and
splitting it would hand `Gene` to `--soloFeatures` and leave the other two as stray positional
arguments. `GeneFull_Ex50pAS` retains intronic reads, which for single-nuclei data is where most
of the signal is; `Gene` is kept alongside it because the exon-only counts are what a per-cell
nuclear fraction is computed from, and `Velocyto` because spliced/unspliced layers cannot be
recovered later without re-running the aligner. There is no `--outdir`: CeleScope writes
relative to the working directory, which is why `work_dir` is a parameter and is passed to the
executor as `cwd`.

WHERE THE RAW MATRIX IS

`<work_dir>/<sample>/outs/raw` is the stable public interface and is searched first. The
STARsolo path underneath it - `01.starsolo/<sample>_Solo.out/<feature>/raw` - is searched as a
fallback, because a pipeline that hard-codes the internal path breaks on a version that renames
a numbered step directory, and reports the successful runs as failures when it does.

A directory that exists is not a matrix. Every candidate must hold a MatrixMarket triplet, the
header dimensions must be non-zero, and they are cross-checked against the line counts of
`features` and `barcodes` - two independent routes to the same two numbers, which is what
distinguishes a complete matrix from a truncated one.

WHAT COULD NOT BE VERIFIED HERE

CeleScope is not installed on the machine this module was written on. Argument construction,
mapfile composition, path arithmetic and output parsing are pure functions with no tool
dependency and are exercisable anywhere; everything that depends on the tool's own behaviour -
whether `multi_rna --version` is accepted, the precise name of the generated script, the exact
output tree of a given release - is asserted at runtime and refuses with the searched paths
named rather than guessing.
"""

from __future__ import annotations

import gzip
import os
import shutil
import sys
import time
from pathlib import Path, PurePosixPath, PureWindowsPath

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from engine.provenance import NOT_INVOKED, tool_version  # noqa: E402
from engine.task import TaskFailure  # noqa: E402

# ---------------------------------------------------------------------------------------------
# The reference invocation. Changing any of these changes what the matrix IS, so they are
# constants with names rather than defaults buried in a signature.
# ---------------------------------------------------------------------------------------------
DEFAULT_EXE = "multi_rna"
SOLO_FEATURES = "Gene GeneFull_Ex50pAS Velocyto"   # one argv element - see the module docstring
REPORT_SOLO_FEATURE = "GeneFull_Ex50pAS"
SOLO_CELL_FILTER = "EmptyDrops_CR"
MOD = "shell"
STAR_EXE = "STAR"

#: STAR index members whose absence means `--genomeDir` is not an index at all. Checked before
#: launch because STAR discovers it hours in, after a scheduler has granted the allocation.
STAR_INDEX_MEMBERS = ("Genome", "SA", "SAindex")
#: A leftover STAR scratch directory means the index build was interrupted. The files that
#: remain load, up to a point, and then fail in a way that reads as a data problem.
STAR_INDEX_INCOMPLETE = "_STARtmp"

#: Spellings that mean "nobody supplied a chemistry". `nan` is in the list because a blank cell
#: read from a samplesheet by pandas arrives as float nan and str()s to exactly that.
CHEMISTRY_NOT_A_VALUE = frozenset(
    {"", "auto", "none", "null", "na", "n/a", "nan", "unknown", "unspecified", "?", "-"})

FASTQ_SUFFIXES = (".fq.gz", ".fastq.gz", ".fq", ".fastq")
#: Trailing read markers, longest first so `_R1` is stripped before `_1` can match anything.
READ_MARKERS = (("_R1", "1"), ("_R2", "2"), ("_r1", "1"), ("_r2", "2"),
                ("_1", "1"), ("_2", "2"))

MATRIX_NAMES = ("matrix.mtx.gz", "matrix.mtx")
BARCODE_NAMES = ("barcodes.tsv.gz", "barcodes.tsv")
FEATURE_NAMES = ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")


# ---------------------------------------------------------------------------------------------
# Unknown is not a value
# ---------------------------------------------------------------------------------------------
#: Type NAMES of the missing-value scalars that are neither None nor a float. Matched by name so
#: that the check costs nothing, and works, when pandas and numpy are not installed - the CLI has
#: to stay importable on a bare clone, so neither may be imported at module scope.
_MISSING_TYPE_NAMES = frozenset({"NAType", "NaTType", "MaskedConstant", "NullScalar"})


def _is_null_arrow_scalar(value) -> bool:
    """True for a pyarrow scalar that holds no value. Duck-typed; pyarrow is never imported.

    A null Arrow scalar is not one type but one per Arrow type - `StringScalar`, `Int64Scalar`,
    `TimestampScalar` - so the type NAME carries the TYPE and says nothing about the nullness,
    and only the untyped `pyarrow.scalar(None)` is a `NullScalar`. Every other route in
    `is_missing` misses it: it is not `None`, it is not a `str`, `x != x` is False for it, and
    `pandas.isna` does not recognise it. What it DOES do is `str()` to the four characters
    `None`, which is how it reached a command line as a directory named after itself.

    It arrives whenever a samplesheet was read through Arrow directly rather than through a
    pandas nullable dtype - pandas hands back `pandas.NA`, pyarrow hands back this - and it is
    the shape the previous review of this predicate missed.

    Asked about itself it answers plainly: `is_valid` is a `bool` and `as_py` is a method. BOTH
    are required before the answer is trusted, so an unrelated object that happens to carry an
    `is_valid` attribute is not read as an Arrow scalar. Duck-typed rather than imported for the
    same reason as everything else here - a bare clone has no pyarrow, and importing one inside
    a predicate that runs on every argument would cost the CLI its start-up time.
    """
    if getattr(type(value), "as_py", None) is None:
        return False
    valid = getattr(value, "is_valid", None)
    return isinstance(valid, bool) and not valid


def is_missing(value) -> bool:
    """True when a value carries no information, in every shape one actually arrives in.

    `None` is the shape everyone remembers and the rarest one in practice. A blank samplesheet
    cell is `float('nan')` when pandas read it through numpy, `pandas.NA` when it read it through
    the nullable or pyarrow-backed dtypes, `pandas.NaT` when the column was parsed as a date, and
    `numpy.ma.masked` when it came through a masked array. None of those is `None`, only the
    first is a `float`, and every one of them survives an `is not None` guard and then compares
    False against every threshold - which downstream is indistinguishable from a value that was
    measured and did not exceed the cut. That reads as a PASS, which is why this predicate exists
    once per module and every "is this unknown?" question in the file goes through it.

    Five routes, cheapest first:

      * identity against `None`;
      * blank and whitespace-only text, in `str` and in `bytes`;
      * the type NAME, which catches `pandas.NA`, `pandas.NaT` and `numpy.ma.masked` without
        importing anything;
      * an Arrow scalar's own `is_valid`, which catches the pyarrow-backed nulls that have no
        shared type name to match - see `_is_null_arrow_scalar`;
      * `value != value`, which catches `float('nan')`, `numpy.float64('nan')`,
        `numpy.float32('nan')` - not a `float` subclass - and `numpy.datetime64('NaT')`.

    pandas is then consulted, but only if it is ALREADY in `sys.modules`: a value cannot be a
    pandas scalar in a process that never imported pandas, so looking there is both sufficient
    and free, while importing pandas inside a predicate that runs on every argument would cost a
    second of start-up to a CLI that does not otherwise need it.

    For a value that may be a numpy boolean the rule in this module is `bool(x)` AFTER
    `is_missing(x)` has been checked, never `x is True` - `numpy.bool_(True) is True` is False,
    so identity reads a genuinely flagged row as unflagged.
    """
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip()
    if isinstance(value, bytes):
        return not value.decode("utf-8", "replace").strip()
    if type(value).__name__ in _MISSING_TYPE_NAMES:
        return True
    if _is_null_arrow_scalar(value):
        return True
    try:
        if bool(value != value):
            return True
    except (TypeError, ValueError):
        pass                      # an object whose __ne__ refuses is not thereby missing
    pandas = sys.modules.get("pandas")
    if pandas is not None and not hasattr(value, "__len__"):
        try:
            verdict = pandas.isna(value)
        except (TypeError, ValueError):
            return False
        if not hasattr(verdict, "__len__"):
            try:
                return bool(verdict)
            except (TypeError, ValueError):
                return False
    return False


def is_true(value, name: str = "flag") -> bool:
    """Truthiness of a flag that may be a numpy boolean, with unknown refused rather than read.

    `numpy.bool_(True) is True` is False, so identity misreads every flag arriving from a
    numpy-backed table; `bool(pandas.NA)` raises. The rule is therefore: refuse unknown, then
    `bool()`.
    """
    if is_missing(value):
        raise TaskFailure(
            f"{name} is {value!r}, which is not True or False. A flag that was never set is not "
            f"the same as one that was set to False, and reading it as False here would record a "
            f"decision nobody made.")
    return bool(value)


def _require_present(value, name: str, sample: str = "?"):
    """Refuse an unknown BEFORE it reaches `Path()`. Returns the value unchanged.

    `command_path` is where a path is spelled and where an unknown one is refused, but several
    checks legitimately run before it - a directory has to be proven to exist before the command
    is built - and each of those constructs a `Path` first. Two things happen there that are not
    a refusal: `Path(None)` raises a `TypeError` naming neither the argument nor the sample, and
    `Path("")` is the CURRENT DIRECTORY, so an unsupplied `work_dir` silently becomes "run here".
    This makes the refusal arrive with the argument's name on it, at the top, either way.
    """
    if is_missing(value):
        raise TaskFailure(
            f"{sample}: {name} is {value!r}, which is not a path, and there is no default.\n"
            f"  Unknown is not a value here. Spelled out it becomes a directory named after the "
            f"sentinel\n  - `None`, `<NA>`, `nan`, `NaT` - and the failure then arrives as a "
            f"missing file at a path\n  nobody wrote, several messages away from the argument "
            f"that was never supplied.")
    return value


def require_chemistry(chemistry, sample: str = "?") -> str:
    """Return the chemistry as a string, or refuse. There is no default and `auto` is not one.

    A cohort whose libraries were built with different chemistry versions, processed under a
    single `--chemistry`, produces near-empty output for the mismatched half and exits zero. The
    resulting matrix is not corrupt in any way a downstream check would notice; it is simply
    small, which reads as a poor library rather than as a wrong argument.

    `auto` is refused as well. Detection may well pick correctly, but the value that governed
    the run is then decided inside the tool and is not the value the samplesheet declares, so
    two runs of the same command are not guaranteed to be the same run.
    """
    if is_missing(chemistry):
        raise TaskFailure(
            f"{sample}: --chemistry is required and has no default.\n"
            f"  It is a per-SAMPLE property. Chemistry has been observed to differ within one\n"
            f"  cohort, and a library processed under the wrong chemistry yields near-empty\n"
            f"  output WITHOUT an error, because the barcode and UMI offsets no longer match\n"
            f"  the read structure. Declare it per sample in the samplesheet.")
    text = str(chemistry).strip()
    if text.lower() in CHEMISTRY_NOT_A_VALUE:
        raise TaskFailure(
            f"{sample}: --chemistry {text!r} is not a chemistry.\n"
            f"  'auto' asks the tool to decide, which means the value that governed the run is\n"
            f"  not the value the samplesheet declares and cannot be reproduced from it. Pass\n"
            f"  the version this library was built with (for example GEXSCOPE-V2 or\n"
            f"  GEXSCOPE-V3 - the accepted names are the ones your CeleScope release lists, and\n"
            f"  this adapter deliberately does not police them).")
    return text


def _require_positive_int(value, name: str, sample: str = "?") -> int:
    """An integer >= 1, or refuse. Never a fallback: a silently halved thread count is a bill."""
    if is_missing(value):
        raise TaskFailure(f"{sample}: {name} is required and has no default")
    try:
        n = int(str(value).strip())
    except (TypeError, ValueError):
        raise TaskFailure(f"{sample}: {name}={value!r} is not an integer") from None
    if n < 1:
        raise TaskFailure(f"{sample}: {name}={n} must be at least 1")
    return n


def is_absolute_path(path) -> bool:
    """Anchored under POSIX **or** Windows rules. Pure.

    `Path.is_absolute()` answers for the host running this code, and that is the wrong host: an
    orchestrator on Windows submitting to a POSIX cluster has `/refs/mm39` reported as relative,
    and a correct run is then refused - a gate firing on correct behaviour, which is how gates
    get switched off.

    A rooted-but-driveless path such as `\\scratch\\run` counts as anchored too. The hazard being
    guarded is a RELATIVE path silently resolving against the work directory; a rooted one does
    not do that, so refusing it would buy nothing and cost a legitimate run.
    """
    text = str(path)
    return (PurePosixPath(text).is_absolute() or PureWindowsPath(text).is_absolute()
            or bool(PureWindowsPath(text).root))


def _require_absolute(path, name: str) -> str:
    """Absolute paths only, because the command runs with `cwd` set to the working directory.

    A relative `--genomeDir` given against the caller's own directory silently resolves against
    the work directory instead. It does not error; it points somewhere else.

    The value is returned as the TEXT it arrived as. Round-tripping it through `Path` would
    rewrite `/refs/mm39` as `\\refs\\mm39` when the orchestrator runs on Windows, which is a
    path no POSIX cluster can open and one no error message would explain.

    An unknown is refused before `str()` is taken of it. Left to the absoluteness test it would
    be reported as "must be an absolute path, got 'None'", which reads as a caller who passed a
    relative path rather than as one who passed nothing.
    """
    _require_present(path, name)
    text = str(path)
    if not is_absolute_path(text):
        raise TaskFailure(
            f"{name} must be an absolute path, got {text!r}.\n"
            f"  CeleScope is launched with the working directory set to the run directory, so a\n"
            f"  relative path is resolved against THAT, not against the caller's directory. The\n"
            f"  failure is silent: the path exists somewhere, just not where it was meant.")
    return text


def command_path(value, name: str, sample: str = "?") -> tuple:
    """`(text, how)` - the spelling a path must have on the command line, and how it was obtained.

    This is the one place a path may be rewritten, and it rewrites only what it must. An ALREADY
    ABSOLUTE path is returned as the caller's own text, untouched: `str(Path('/refs/mm39'))` is
    `\\refs\\mm39` on a Windows orchestrator, so a `Path` round trip is a silent rewrite into a
    path the POSIX cluster cannot open, and it arrives as an error about a missing reference
    rather than as an error about a mangled argument.

    A RELATIVE path cannot be left alone: CeleScope runs with `cwd` set to the work directory, so
    a relative argument resolves against THAT and quietly points somewhere else. It is resolved
    against this host's working directory - the only host that can - and `how` says so, so the
    rewrite appears in the run's recorded parameters instead of happening silently.

    AN UNKNOWN IS REFUSED BEFORE ANY OF THAT HAPPENS, and the order is the whole point. This
    function used to take `str(value)` on its first line, which is where a sentinel stops being
    recognisable: `None` becomes the four characters `None`, `pandas.NA` becomes `<NA>`,
    `pandas.NaT` becomes `NaT`, and `float('nan')` becomes `nan`. Each is then a perfectly
    ordinary RELATIVE path, so the function resolved it against this host's working directory and
    handed back something like `<cwd>/None` with `how` reporting, accurately, that it had been
    resolved. An unsupplied `--genomeDir` reached the aligner as a directory named after the
    sentinel, and the run failed hours later as a missing file at a path nobody had written -
    a message pointing at the filesystem instead of at the argument that was never given.

    Every argv element, mapfile column, `cwd` and recorded metric in this module comes through
    here, so this is the one check that has to be made in the right order.
    """
    _require_present(value, name, sample)
    text = str(value)
    if is_absolute_path(text):
        return text, "verbatim: absolute as supplied"
    if is_missing(text):
        # Not dead: `value` can be a non-string object - a `Path(' ')`, a wrapper - that is not
        # itself unknown but whose text is blank, and `Path('').resolve()` is the cwd.
        raise TaskFailure(f"{sample}: {name} is empty; a path has no default here")
    resolved = str(Path(text).resolve())
    return resolved, (f"resolved against the orchestrator's working directory "
                      f"{Path.cwd()} because it was supplied relative")


# ---------------------------------------------------------------------------------------------
# Mapfile
# ---------------------------------------------------------------------------------------------
def derive_fastq_prefixes(names) -> dict:
    """Group FASTQ file names by the prefix that remains once the read marker is removed.

    Pure: it takes names, not a directory, so the grouping rules can be tested against any
    naming convention without a filesystem. Returns
    `{prefix: {"1": [...], "2": [...], "unpaired": [...]}}`.

    `_001` is stripped before the read marker because the bcl2fastq convention puts the chunk
    number last. A name carrying no recognised marker is recorded under `unpaired` rather than
    dropped, so a caller can print what was actually there when the grouping is ambiguous.
    """
    groups: dict = {}
    for raw in names:
        name = str(raw)
        stem = name
        for suffix in sorted(FASTQ_SUFFIXES, key=len, reverse=True):
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        else:
            continue                       # not a FASTQ - checksums and manifests live here too
        if stem.endswith("_001"):
            stem = stem[:-4]
        read = "unpaired"
        for marker, which in READ_MARKERS:
            if stem.endswith(marker):
                stem, read = stem[: -len(marker)], which
                break
        slot = groups.setdefault(stem, {"1": [], "2": [], "unpaired": []})
        slot[read].append(name)
    return groups


def discover_fastq_prefix(fastq_dir, sample: str = "?") -> str:
    """The single FASTQ prefix in a directory, or a refusal naming what was found.

    Ambiguity is not resolved by picking one. Two prefixes in a directory means either two
    libraries were staged together - in which case a mapfile naming one of them quietly analyses
    half the data - or the naming does not follow a convention this function recognises. Both
    are for the caller to settle by passing `fastq_prefix` explicitly.
    """
    _require_present(fastq_dir, "fastq_dir", sample)
    d = Path(fastq_dir)
    if not d.is_dir():
        raise TaskFailure(f"{sample}: FASTQ directory does not exist: {d}")
    entries = sorted(p.name for p in d.iterdir() if p.is_file())
    groups = derive_fastq_prefixes(entries)
    if not groups:
        listing = "\n".join(f"    {n}" for n in entries[:20]) or "    (directory is empty)"
        raise TaskFailure(
            f"{sample}: no FASTQ files in {d}\n"
            f"  looked for names ending in {', '.join(FASTQ_SUFFIXES)}; the directory holds:\n"
            f"{listing}")
    if len(groups) > 1:
        found = ", ".join(sorted(groups))
        raise TaskFailure(
            f"{sample}: {len(groups)} FASTQ prefixes in {d}: {found}\n"
            f"  A mapfile names ONE prefix. Naming one of several here would analyse part of\n"
            f"  the input and report a complete run. Pass fastq_prefix= explicitly, or stage\n"
            f"  one library per directory.")
    prefix, reads = next(iter(groups.items()))
    if not reads["1"] or not reads["2"]:
        have = {k: v for k, v in reads.items() if v}
        raise TaskFailure(
            f"{sample}: prefix {prefix!r} in {d} is not a read pair - found {have}\n"
            f"  Both reads are required: one carries the barcode and UMI, the other the cDNA.")
    return prefix


def mapfile_line(fastq_prefix: str, fastq_dir, sample: str) -> str:
    """One CeleScope mapfile row: fastq_prefix, fastq_dir, sample - tab separated.

    Pure. The three fields are checked for tabs and newlines because either would shift every
    later field by one column, and a mapfile whose sample column holds a directory path is read
    without complaint.

    Each field is tested for unknown BEFORE `str()` is taken of it, for the reason spelled out in
    `command_path`: the blank test below sees `str(None)`, which is four non-blank characters and
    passes it. A mapfile is a data file, not a command line - the sentinel would be written into
    the column, CeleScope would resolve it as a directory name, and the run would report an empty
    sample rather than a bad argument.
    """
    for name, value in (("fastq_prefix", fastq_prefix), ("fastq_dir", fastq_dir),
                        ("sample", sample)):
        if is_missing(value):
            raise TaskFailure(
                f"mapfile {name} is {value!r}; a mapfile column has no default. Written out it "
                f"would become a column reading `{value!s}`, which CeleScope reads as an ordinary "
                f"value.")
    fields = {"fastq_prefix": str(fastq_prefix), "fastq_dir": str(fastq_dir),
              "sample": str(sample)}
    for name, value in fields.items():
        if not value.strip():
            raise TaskFailure(f"mapfile {name} is empty; a mapfile column has no default")
        if "\t" in value or "\n" in value or "\r" in value:
            raise TaskFailure(
                f"mapfile {name}={value!r} contains a tab or newline. The mapfile is tab "
                f"separated, so this would silently shift every following column.")
    return "{fastq_prefix}\t{fastq_dir}\t{sample}\n".format(**fields)


def write_mapfile(sample: str, fastq_dir, out, fastq_prefix=None) -> Path:
    """Write the one-row mapfile CeleScope reads, and verify it round-trips.

    `out` may be the file to write or a directory to write `<sample>.mapfile` into. The file is
    read back and compared after writing: a mapfile that was not written where it was believed
    to be produces a `multi_rna` run over an empty sample set, which exits zero.

    When `fastq_prefix` is omitted it is DISCOVERED from the directory, never assumed - see
    `discover_fastq_prefix`, which refuses rather than choosing between candidates.

    The directory column holds the caller's own text when that text is absolute, and a locally
    resolved form only when it was relative - see `command_path`. CeleScope resolves this column
    against ITS working directory, so a relative value here points somewhere else; a `Path` round
    trip of an absolute one rewrites a POSIX path with backslashes on a Windows orchestrator, and
    the mapfile then names a directory that exists on neither machine.
    """
    if is_missing(sample):
        raise TaskFailure("sample is required and has no default")
    sample = str(sample).strip()
    _require_present(fastq_dir, "fastq_dir", sample)
    _require_present(out, "out", sample)
    d = Path(fastq_dir)
    if not d.is_dir():
        raise TaskFailure(
            f"{sample}: FASTQ directory does not exist: {d}\n"
            f"  CeleScope resolves the mapfile's directory column at run time and reports an "
            f"empty sample rather than a missing path.")
    fastq_text, _ = command_path(fastq_dir, "fastq_dir", sample)

    prefix = str(fastq_prefix).strip() if not is_missing(fastq_prefix) \
        else discover_fastq_prefix(d, sample)

    target = Path(out)
    if target.is_dir():
        target = target / f"{sample}.mapfile"
    target.parent.mkdir(parents=True, exist_ok=True)

    line = mapfile_line(prefix, fastq_text, sample)
    target.write_text(line, encoding="utf-8")

    if not target.is_file():
        raise TaskFailure(f"{sample}: mapfile was not created at {target}")
    written = target.read_text(encoding="utf-8")
    if written != line:
        raise TaskFailure(
            f"{sample}: mapfile at {target} does not read back as it was written.\n"
            f"  wrote: {line!r}\n  read : {written!r}")
    return target


# ---------------------------------------------------------------------------------------------
# Argument construction - pure, testable without the tool
# ---------------------------------------------------------------------------------------------
def build_command(mapfile, genome_dir, chemistry, thread, sample: str = "?",
                  exe: str = DEFAULT_EXE, solo_features: str = SOLO_FEATURES,
                  report_solo_feature: str = REPORT_SOLO_FEATURE,
                  solo_cell_filter: str = SOLO_CELL_FILTER, mod: str = MOD,
                  extra=()) -> list:
    """The `multi_rna` argv, in the reference cohort's order. Pure - no filesystem, no tool.

    `solo_features` is emitted as a single argv element. It contains spaces and that is correct:
    the reference invocation quotes it, and splitting it into three elements would pass only
    `Gene` to the flag.

    `extra` is appended verbatim for flags this adapter does not model. Anything placed there is
    a departure from the reference invocation and belongs in the run's recorded parameters, not
    only in the command line. Verbatim does not extend to unknown: an element that is `None` or
    `pandas.NA` would be spelled onto the command line as `None` or `<NA>` and read by the tool
    as a value, so it is refused with its position named.
    """
    chem = require_chemistry(chemistry, sample)
    n = _require_positive_int(thread, "thread", sample)
    mapfile = _require_absolute(mapfile, "mapfile")
    genome_dir = _require_absolute(genome_dir, "genomeDir")
    for name, value in (("solo_features", solo_features),
                        ("report_solo_feature", report_solo_feature),
                        ("solo_cell_filter", solo_cell_filter), ("mod", mod), ("exe", exe)):
        if is_missing(value):
            raise TaskFailure(f"{sample}: {name} is required; it defines what the matrix IS")
    extra = list(extra)
    for i, item in enumerate(extra):
        if is_missing(item):
            raise TaskFailure(
                f"{sample}: extra[{i}] is {item!r}, which is not a command-line argument. "
                f"str() would place `{item!s}` on the argv and the tool would read it as a value.")
    return [str(exe),
            "--mapfile", mapfile,
            "--genomeDir", genome_dir,
            "--chemistry", chem,
            "--soloFeatures", str(solo_features),
            "--report_soloFeature", str(report_solo_feature),
            "--soloCellFilter", str(solo_cell_filter),
            "--thread", str(n),
            "--mod", str(mod),
            *[str(x) for x in extra]]


def build_env(env_bin=None, base_path=None, pathsep=None, extra=None) -> dict:
    """The environment overlay for the executor, with the tool environment's bin FIRST on PATH.

    This is the second hazard from the module docstring made concrete. CeleScope shells out to
    STAR by bare name and does not check that it was found, so an environment that is installed
    but not on `PATH` fails inside the generated script with a message that names neither.

    `pathsep` is a parameter rather than always `os.pathsep` because an orchestrator running on
    one platform may be submitting to another: `os.pathsep` is `;` on Windows and would build a
    PATH no POSIX cluster can read. Pass `':'` explicitly when submitting across platforms.

    A `PATH` inside `extra` is MERGED, never overwritten and never discarded: it becomes the base
    that `env_bin` is prepended to. Either half winning outright loses the other silently - a
    PATH is still present and still plausible afterwards - and the symptom is the bare-name STAR
    failure this function exists to prevent.
    """
    env = {str(k): str(v) for k, v in (extra or {}).items()}
    if is_missing(env_bin):
        return env
    sep = os.pathsep if pathsep is None else str(pathsep)
    if "PATH" in env:
        base = env["PATH"]
    elif base_path is None:
        base = os.environ.get("PATH", "")
    else:
        base = str(base_path)
    # Not routed through Path: on a Windows orchestrator that rewrites a POSIX bin directory
    # with backslashes, and the resulting PATH silently matches nothing on the cluster.
    head = str(env_bin).strip()
    env["PATH"] = f"{head}{sep}{base}" if base else head
    return env


def resolve_on_path(names, env_bin=None, base_path=None, pathsep=None) -> dict:
    """Which of `names` is findable on the PATH `build_env` constructs. Values are paths or None.

    None means NOT FOUND and is reported as its own outcome by the caller; it is never allowed
    to stand in for a path. The lookup happens on the host that builds the command, which for a
    scheduler is the submit host - stated so a mismatch between submit and compute hosts is
    read as a limitation of this check rather than as evidence about the compute node.
    """
    env = build_env(env_bin, base_path=base_path, pathsep=pathsep)
    search = env.get("PATH") if "PATH" in env else base_path
    out = {}
    for name in names:
        out[str(name)] = shutil.which(str(name), path=search)
    return out


def stage_log(log, stage: str) -> Path:
    """`<log>` for the first command and `<log stem>.<stage><suffix>` for the others.

    Two commands run here and each gets its own captured log. Sharing one path would leave the
    planner's output overwritten by the aligner's, and the planner's output is where a mapfile
    that matched no sample is visible.

    An unknown `log` is refused rather than turned into `None.plan` beside the working directory:
    a log written somewhere nobody looks is the same as no log, and it is discovered only when
    the run has already failed for a different reason.
    """
    _require_present(log, "log")
    p = Path(log)
    return p.with_name(f"{p.stem}.{stage}{p.suffix}")


# ---------------------------------------------------------------------------------------------
# Output location and parsing - pure where it can be, filesystem-only where it cannot
# ---------------------------------------------------------------------------------------------
def raw_matrix_candidates(work_dir, sample: str,
                          feature: str = REPORT_SOLO_FEATURE) -> list:
    """Where the raw matrix could be, most stable first. Pure path arithmetic - nothing is read.

    `outs/raw` is CeleScope's published output location and does not carry a numbered step
    directory in its name. Hard-coding the STARsolo path underneath it has been observed to
    report successful runs as failures after a release renamed the step directory, so it is the
    fallback rather than the first choice.
    """
    root = Path(work_dir) / str(sample)
    return [root / "outs" / "raw",
            root / "01.starsolo" / f"{sample}_Solo.out" / feature / "raw",
            root / "01.starsolo" / "Solo.out" / feature / "raw"]


def matrix_triplet(directory) -> dict:
    """The three MatrixMarket files in a directory, or `{}` if the set is not complete.

    An empty result is deliberately falsy and deliberately not a partial dict: two files out of
    three is not most of a matrix, it is a directory that will fail to load.
    """
    d = Path(directory)
    if not d.is_dir():
        return {}
    found = {}
    for key, names in (("matrix", MATRIX_NAMES), ("barcodes", BARCODE_NAMES),
                       ("features", FEATURE_NAMES)):
        for name in names:
            candidate = d / name
            if candidate.is_file():
                found[key] = candidate
                break
    return found if len(found) == 3 else {}


def find_raw_matrix(work_dir, sample: str, feature: str = REPORT_SOLO_FEATURE) -> Path:
    """The directory holding the raw, unfiltered droplet matrix, or a refusal naming the search.

    A candidate that EXISTS but holds no complete triplet is reported separately from one that
    does not exist, because they mean different things: the first is a run that produced a
    directory and no matrix - the failure mode where a step "succeeds" and writes nothing - and
    the second is usually a version whose tree this adapter has not been taught.
    """
    searched, empty = [], []
    for candidate in raw_matrix_candidates(work_dir, sample, feature):
        searched.append(str(candidate))
        if not candidate.is_dir():
            continue
        if matrix_triplet(candidate):
            return candidate
        empty.append(str(candidate))

    root = Path(work_dir) / str(sample)
    if root.is_dir():
        for candidate in sorted(root.glob(f"*/*Solo.out/{feature}/raw")):
            searched.append(str(candidate))
            if matrix_triplet(candidate):
                return candidate
            empty.append(str(candidate))

    detail = "\n".join(f"    {s}" for s in searched)
    note = ""
    if empty:
        note = ("\n  These EXIST but hold no complete matrix triplet, which means the run "
                "produced a\n  directory and no counts:\n"
                + "\n".join(f"    {e}" for e in empty))
    raise TaskFailure(
        f"{sample}: no raw (unfiltered) matrix found under {Path(work_dir) / str(sample)}\n"
        f"  searched, in order:\n{detail}{note}\n"
        f"  A triplet is one of {MATRIX_NAMES} plus one of {BARCODE_NAMES} plus one of "
        f"{FEATURE_NAMES}.\n"
        f"  The raw matrix is not optional here: the ambient model learns its background from "
        f"the\n  empty droplets, and `outs/filtered` no longer contains them.")


def parse_mtx_header(lines, name: str = "<matrix>") -> dict:
    """Dimensions from a MatrixMarket header. Pure - it takes lines, so it is testable on text.

    Returns rows as `n_genes` and columns as `n_barcodes`, which is the features-by-barcodes
    orientation STARsolo and Cell Ranger both write. A zero in either dimension, or zero stored
    entries, is refused rather than returned: an empty matrix loads, plots and summarises
    without complaint, and is the output a wrong `--chemistry` produces.
    """
    banner = None
    for raw in lines:
        line = str(raw).strip()
        if not line:
            continue
        if banner is None:
            if not line.startswith("%%MatrixMarket"):
                raise TaskFailure(
                    f"{name}: not a MatrixMarket file - first line is {line[:60]!r}")
            banner = line
            continue
        if line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) != 3:
            raise TaskFailure(
                f"{name}: expected three integers on the dimension line, got {line[:60]!r}")
        try:
            rows, cols, entries = (int(p) for p in parts)
        except ValueError:
            raise TaskFailure(
                f"{name}: dimension line is not three integers: {line[:60]!r}") from None
        if rows <= 0 or cols <= 0 or entries <= 0:
            raise TaskFailure(
                f"{name}: matrix declares {rows} genes x {cols} barcodes with {entries} stored "
                f"entries.\n  An empty matrix is what a mismatched --chemistry produces, and it "
                f"exits zero.")
        return {"banner": banner, "n_genes": rows, "n_barcodes": cols, "n_entries": entries}
    raise TaskFailure(f"{name}: no dimension line found in the MatrixMarket header")


def _open_text(path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(str(p), "rt", encoding="utf-8", errors="replace")
    return open(str(p), "rt", encoding="utf-8", errors="replace")


def read_mtx_header(path, max_lines: int = 64) -> dict:
    """`parse_mtx_header` over the first lines of a file, gzipped or not."""
    p = Path(path)
    if not p.is_file():
        raise TaskFailure(f"matrix file does not exist: {p}")
    lines = []
    with _open_text(p) as fh:
        for i, line in enumerate(fh):
            lines.append(line)
            if i + 1 >= max_lines:
                break
    return parse_mtx_header(lines, name=str(p))


def count_lines(path) -> int:
    """Line count of a text file, gzipped or not, counted in binary and without loading it.

    A final line with no trailing newline is counted. Under-counting by one here would look like
    a one-row disagreement between the matrix header and its features file, which is exactly the
    kind of near-miss that gets explained away.
    """
    p = Path(path)
    if not p.is_file():
        raise TaskFailure(f"file does not exist: {p}")
    total, last = 0, b""
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(str(p), "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            total += chunk.count(b"\n")
            last = chunk[-1:]
    if last and last != b"\n":
        total += 1
    return total


def verify_matrix_dims(triplet: dict, header: dict, sample: str = "?",
                       check_barcodes: bool = True) -> dict:
    """Cross-check the header's dimensions against the two label files, and refuse a mismatch.

    The same two numbers by two independent routes. A truncated transfer, a half-written file or
    a transposed matrix all produce a header that parses cleanly on its own; only the comparison
    against the labels shows it. `check_barcodes` exists because a raw droplet matrix names
    millions of barcodes and a caller may have a reason to skip that read - it is never skipped
    silently, and the outcome is recorded as its own value.
    """
    n_features = count_lines(triplet["features"])
    if n_features != header["n_genes"]:
        raise TaskFailure(
            f"{sample}: matrix header declares {header['n_genes']:,} genes but "
            f"{triplet['features'].name} has {n_features:,} rows.\n"
            f"  These must agree. They disagree when a file was truncated in transfer, when the "
            f"matrix\n  is transposed, or when the two files came from different runs.")
    out = {"features_rows": n_features}
    if not check_barcodes:
        out["barcodes_rows"] = "not counted"
        return out
    n_barcodes = count_lines(triplet["barcodes"])
    if n_barcodes != header["n_barcodes"]:
        raise TaskFailure(
            f"{sample}: matrix header declares {header['n_barcodes']:,} barcodes but "
            f"{triplet['barcodes'].name} has {n_barcodes:,} rows.")
    out["barcodes_rows"] = n_barcodes
    return out


def has_crlf(data) -> bool:
    """True if this content carries Windows line endings. Pure. It has to be given BYTES.

    A generated shell script delivered with CRLF fails at the shebang with `bad interpreter`,
    the job leaves the queue in under a second, and the driver reports it identically to a run
    that executed and produced nothing - which sends the diagnosis into the chemistry or the
    reference instead of into the file.

    This function was correct and its only call site made it useless: it was fed
    `Path.read_text()`, which opens in UNIVERSAL-NEWLINE mode and has already translated every
    CRLF to LF before the caller sees a single character. Fed that, it cannot return True on any
    file, so the guard could never fire and a CRLF script would have gone to the queue with the
    check reporting it clean. The call site now reads bytes. A `str` is still accepted so the
    rule can be tested on literal text, but a `str` that came from `read_text()` is not evidence
    about a file - only `read_bytes()` is.
    """
    if isinstance(data, (bytes, bytearray, memoryview)):
        return b"\r\n" in bytes(data)
    return "\r\n" in str(data)


def find_shell_script(work_dir, sample: str) -> Path:
    """The script `multi_rna --mod shell` wrote for this sample, or a refusal listing the tree.

    `multi_rna` plans and does not run. If this file is absent the planner accepted the mapfile
    and produced nothing from it - usually a mapfile whose sample column does not match, which
    is not an error to CeleScope because an empty sample set is a legal one.
    """
    shell_dir = Path(work_dir) / "shell"
    script = shell_dir / f"{sample}.sh"
    if script.is_file():
        return script
    if not shell_dir.is_dir():
        raise TaskFailure(
            f"{sample}: {shell_dir} was not created.\n"
            f"  `multi_rna --mod {MOD}` writes the per-sample scripts there and runs nothing "
            f"itself.\n  A zero exit with no shell/ directory means the mapfile named no sample "
            f"the planner\n  recognised - check the mapfile's third column against --sample.")
    present = sorted(p.name for p in shell_dir.iterdir())
    listing = "\n".join(f"    {n}" for n in present[:20]) or "    (empty)"
    raise TaskFailure(
        f"{sample}: {script} was not written.\n"
        f"  {shell_dir} contains:\n{listing}\n"
        f"  The script name is taken from the mapfile's sample column; if the names above look "
        f"like\n  the sample under a different spelling, the mapfile and the task disagree.")


def _mtime(path):
    """The modification time of a path, or None. Never raises: this is corroboration, not proof."""
    try:
        return Path(path).stat().st_mtime
    except OSError:
        return None


def refuse_leftovers(work_dir, sample: str) -> dict:
    """Refuse to start when a previous run's tree or planner script is in the way.

    This is the check that makes every later one mean something. `find_raw_matrix` proves a
    matrix is present; it cannot prove that THIS invocation produced it, and a CeleScope run that
    exits zero having written nothing - a mapfile matching no sample, a chemistry that yields
    near-empty output, a generated script that died at its shebang - passes every check in this
    module the moment leftovers are lying at the paths those checks look at. The previous run's
    numbers are then recorded under this run's parameters, and nothing in the output says so.

    Nothing is deleted. Two reasons, and the second is the load-bearing one:

      * a previous tree is the record of what happened to that run, exactly as in the Cell Ranger
        adapter, and deleting evidence to make a re-run convenient is not this adapter's call;
      * `shell/<sample>.sh` is EXECUTED by phase two. A leftover script runs the earlier run's
        chemistry and reference while this run records the parameters it was asked for.

    Returns what it verified, so the run's metrics can state that the paths were empty rather
    than leave a reader to assume it.
    """
    work = Path(work_dir)
    sample_dir = work / str(sample)
    script = work / "shell" / f"{sample}.sh"
    if sample_dir.exists():
        raise TaskFailure(
            f"{sample}: {sample_dir} already exists.\n"
            f"  CeleScope writes its per-sample output tree there. This adapter neither writes "
            f"into an\n  existing one nor deletes it: a previous tree is the record of what "
            f"happened to that run,\n  and - the reason this refusal is load-bearing - a run that "
            f"exits 0 having written nothing\n  is indistinguishable from a successful one while "
            f"an earlier run's matrix sits at the path\n  the check looks at. Move it aside, or "
            f"pass a different work_dir. Whether this step needs\n  to run at all is a resume "
            f"question and belongs to the run manifest, not here.")
    if script.exists():
        raise TaskFailure(
            f"{sample}: {script} already exists.\n"
            f"  `multi_rna --mod {MOD}` writes that script and phase two EXECUTES it. A leftover "
            f"from an\n  earlier run would be executed instead - the earlier run's chemistry, the "
            f"earlier run's\n  reference - while this run records the parameters it was asked "
            f"for, and the two would\n  never be compared. Move it aside, or pass a different "
            f"work_dir.")
    return {"sample_dir_absent_before_run": True, "shell_script_absent_before_run": True,
            "checked": [str(sample_dir), str(script)]}


# ---------------------------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------------------------
def version(exe: str = DEFAULT_EXE, args=("--version",)) -> str:
    """The observed version string, or `not invoked` - never a value read from a lockfile.

    Whether `multi_rna` itself answers `--version` has not been confirmed against the tool. If
    it does not, `engine.provenance.tool_version` returns `not invoked`, and that is reported as
    such rather than filled in from the environment name, which would be a fabricated
    provenance record and indistinguishable from a real one.
    """
    return tool_version(exe, tuple(args))


# ---------------------------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------------------------
def run_celescope(sample: str, fastq_dir, genome_dir, chemistry, work_dir, log,
                  *, threads, env_bin, fastq_prefix=None, mapfile=None,
                  exe: str = DEFAULT_EXE, solo_features: str = SOLO_FEATURES,
                  report_solo_feature: str = REPORT_SOLO_FEATURE,
                  solo_cell_filter: str = SOLO_CELL_FILTER, mod: str = MOD,
                  require_on_path=(STAR_EXE,), require_index_files=STAR_INDEX_MEMBERS,
                  pathsep=None, check_barcode_count: bool = True, timeout_s=None,
                  extra=(), shell_exe: str = "bash", executor=None) -> dict:
    """Run CeleScope for one sample and return the raw matrix, or fail with what was searched.

    Two commands, two logs: the planner (`multi_rna --mod shell`), then the script it wrote.
    Nothing about the sample, the chemistry or the reference is decided here - the orchestrator
    decides and this executes.

    `threads` and `env_bin` are keyword arguments WITHOUT defaults on purpose. A thread count
    that silently defaults wastes a scheduler allocation; an environment bin that silently
    defaults reintroduces the PATH hazard this adapter exists partly to close. Pass `env_bin=None`
    to assert deliberately that the environment is already on `PATH` - the assertion is then
    checked, because `require_on_path` is resolved either way.

    EVERY DECLARED OUTPUT IS PROVEN ABSENT BEFORE THE FIRST COMMAND RUNS. `refuse_leftovers`
    refuses when `<work_dir>/<sample>` or `shell/<sample>.sh` exists, so a file found after the
    commands return was written BY THIS INVOCATION - which existence alone can never show. That
    refusal, not a timestamp, is the proof: mtimes are stamped by the compute node's clock while
    the refusal is made on the orchestrator's, and a freshness gate that fires on clock skew
    between the two is a gate someone switches off. The mtimes are recorded beside the outputs as
    corroboration.

    Paths reach the command line, the mapfile and `cwd` as the caller's own text when that text is
    absolute; only a relative path is resolved, and `metrics["path_resolution"]` says which and
    against what. See `command_path`.

    Returns `{"outputs": [...], "metrics": {...}, "versions": {...}}`. `outputs` lists only files
    that were checked to exist after the commands returned; the raw matrix DIRECTORY is reported
    as `metrics["raw_matrix"]`.
    """
    if executor is None:
        raise TaskFailure("run_celescope requires an executor; there is no in-process fallback")
    if is_missing(sample):
        raise TaskFailure("sample is required and has no default")
    sample = str(sample).strip()
    if str(mod) != MOD:
        raise TaskFailure(
            f"{sample}: this adapter implements --mod {MOD} only, got {mod!r}.\n"
            f"  The other modes submit the work themselves, so the second phase here - locating "
            f"and\n  running shell/{sample}.sh - would either duplicate the run or wait on "
            f"nothing.")

    n_threads = _require_positive_int(threads, "threads", sample)
    chem = require_chemistry(chemistry, sample)
    count_barcodes = is_true(check_barcode_count, "check_barcode_count")

    # Every path argument is tested for unknown before it reaches `Path()`. `command_path` refuses
    # these too, but it runs AFTER the existence checks below, and `Path(None)` raises a TypeError
    # there that names neither the argument nor the sample while `Path("")` is the current
    # directory. The refusal has to carry the argument's name to be worth anything.
    for _value, _name in ((fastq_dir, "fastq_dir"), (genome_dir, "--genomeDir"),
                          (work_dir, "work_dir"), (log, "log")):
        _require_present(_value, _name, sample)

    fq = Path(fastq_dir)
    if not fq.is_dir():
        raise TaskFailure(f"{sample}: FASTQ directory does not exist: {fq}")
    ref = Path(genome_dir)
    if not ref.is_dir():
        raise TaskFailure(f"{sample}: --genomeDir does not exist: {ref}")
    missing_index = [n for n in require_index_files if not (ref / n).exists()]
    if missing_index:
        raise TaskFailure(
            f"{sample}: --genomeDir {ref} is missing {', '.join(missing_index)}.\n"
            f"  A STAR index that is absent is discovered inside STAR, after the scheduler has "
            f"granted\n  the allocation, and reported as a file error rather than as a wrong "
            f"reference.")
    if (ref / STAR_INDEX_INCOMPLETE).exists():
        raise TaskFailure(
            f"{sample}: --genomeDir {ref} still contains {STAR_INDEX_INCOMPLETE}, so the index "
            f"build did not finish.\n  The files that are there load up to a point and then fail "
            f"in a way that reads as a data problem.")

    # The command's own spelling of each path. Absolute text is passed through untouched; only a
    # relative path is resolved, and how each was obtained is recorded rather than assumed.
    fq_text, fq_route = command_path(fastq_dir, "fastq_dir", sample)
    ref_text, ref_route = command_path(genome_dir, "--genomeDir", sample)

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    work_text, work_route = command_path(work_dir, "work_dir", sample)

    # Pattern: a previous run's output must never be accepted as this run's. Checked BEFORE
    # anything is written, so that what is found afterwards can only have come from this run.
    leftovers = refuse_leftovers(work, sample)
    started = time.time()

    # `is_missing`, not `is None`: an unsupplied mapfile arrives as `None` from a caller, as
    # `pandas.NA` or `nan` from a samplesheet column, and as `""` from an argument parser. All
    # three mean the same thing - write it into the work directory - and `Path(pandas.NA)` raises
    # while `Path("")` would put the mapfile in the current directory instead.
    map_target = work if is_missing(mapfile) else Path(mapfile)
    map_path = write_mapfile(sample, fastq_dir, map_target, fastq_prefix)
    map_text, map_route = command_path(map_path, "mapfile", sample)

    cmd = build_command(map_text, ref_text, chem, n_threads, sample=sample, exe=exe,
                        solo_features=solo_features, report_solo_feature=report_solo_feature,
                        solo_cell_filter=solo_cell_filter, mod=mod, extra=extra)
    env = build_env(env_bin, pathsep=pathsep)

    # The PATH hazard, checked before anything expensive starts.
    wanted = [str(x) for x in require_on_path]
    resolved = resolve_on_path([str(exe), *wanted], env_bin, pathsep=pathsep)
    unresolved = [n for n in wanted if resolved[n] is None]
    if unresolved:
        where = env["PATH"] if "PATH" in env else os.environ.get("PATH", "")
        raise TaskFailure(
            f"{sample}: {', '.join(unresolved)} not found on the PATH this run would use.\n"
            f"  CeleScope invokes STAR by BARE NAME through a shell and does not check that it "
            f"was\n  found, so this failure would otherwise surface hours later, inside the "
            f"generated\n  script, as a missing-file message naming neither STAR nor the "
            f"environment.\n"
            f"  Pass env_bin=<environment>/bin. PATH searched:\n    {where}\n"
            f"  (resolved on the host building this command; on a scheduler that is the submit "
            f"host.)")

    plan_log = stage_log(log, "plan")
    executor.shell(cmd, log=plan_log, env=env, cwd=work_text, timeout_s=timeout_s)

    # `shell/<sample>.sh` was proven absent above, so finding it here means the planner wrote it
    # on this invocation. That is the whole of the freshness argument; the mtime is recorded, not
    # relied on, because it comes from a different clock than the one `started` was read from.
    script = find_shell_script(work, sample)
    # READ AS BYTES. `read_text()` opens in universal-newline mode and translates CRLF to LF
    # before the caller sees it, so a CRLF check over that text can never fire - which is exactly
    # what this guard used to do.
    body = script.read_bytes()
    if has_crlf(body):
        raise TaskFailure(
            f"{sample}: {script} has Windows line endings.\n"
            f"  The shebang then reads as `/bin/bash^M` and the job dies before it starts, in "
            f"under a\n  second, leaving output that is indistinguishable from a run that "
            f"executed and\n  produced nothing. Fix the transfer, not the script.")

    run_log = stage_log(log, "run")
    executor.shell([str(shell_exe), str(script)], log=run_log, env=env, cwd=work_text,
                   timeout_s=timeout_s)

    raw_dir = find_raw_matrix(work, sample, report_solo_feature)
    triplet = matrix_triplet(raw_dir)
    if not triplet:
        raise TaskFailure(f"{sample}: {raw_dir} lost its matrix triplet between checks")
    header = read_mtx_header(triplet["matrix"])
    dims = verify_matrix_dims(triplet, header, sample, check_barcodes=count_barcodes)

    outputs = [map_path, script, triplet["matrix"], triplet["barcodes"], triplet["features"]]
    for p in outputs:
        if not Path(p).is_file():
            raise TaskFailure(f"{sample}: promised output is absent after the run: {p}")
    # Every one of these sits under a path that did not exist when this function started, so its
    # presence is proof of authorship rather than of luck. Recorded as ages so a reader can see
    # for themselves; a negative or wildly large age means the two hosts disagree about the time,
    # which is worth knowing and is not on its own a reason to fail a run that produced output.
    output_ages = {}
    for p in outputs:
        m = _mtime(p)
        output_ages[Path(p).name] = None if m is None else round(started - m, 3)

    star_path = resolved[STAR_EXE] if STAR_EXE in resolved else None
    exe_path = resolved[str(exe)] if resolved[str(exe)] else str(exe)
    versions = {
        "celescope": version(exe_path),
        "STAR": tool_version(star_path) if star_path else NOT_INVOKED,
    }

    metrics = {
        "sample": sample,
        "chemistry": chem,
        "threads": n_threads,
        "soloFeatures": str(solo_features),
        "report_soloFeature": str(report_solo_feature),
        "soloCellFilter": str(solo_cell_filter),
        "mod": str(mod),
        "fastq_dir": fq_text,
        "genomeDir": ref_text,
        "work_dir": work_text,
        "path_resolution": {"fastq_dir": fq_route, "genomeDir": ref_route,
                            "work_dir": work_route, "mapfile": map_route},
        "outputs_absent_before_run": leftovers,
        "output_age_s_at_start": output_ages,
        "started_epoch_s": started,
        "freshness_proof": ("the sample tree and shell/<sample>.sh were proven absent before the "
                            "planner ran, so every output found afterwards was written by this "
                            "invocation; the ages above are corroboration from the compute "
                            "node's clock, not the proof"),
        "mapfile": str(map_path),
        "shell_script": str(script),
        "raw_matrix": str(raw_dir),
        "n_genes": header["n_genes"],
        "n_barcodes": header["n_barcodes"],
        "n_entries": header["n_entries"],
        "features_rows": dims["features_rows"],
        "barcodes_rows": dims["barcodes_rows"],
        "env_bin": str(env_bin) if not is_missing(env_bin) else "not set",
        "star_resolved": str(star_path) if star_path else "not on PATH",
        "logs": [str(plan_log), str(run_log)],
    }
    return {"outputs": outputs, "metrics": metrics, "versions": versions}
