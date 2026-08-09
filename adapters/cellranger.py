# Execution adapter: FASTQ -> raw unfiltered count matrix, via `cellranger count`.
# It removes no observation. Cell Ranger's own cell call is a by-product of this step and is
# never this pipeline's cell call; what the step delivers downstream is outs/raw_feature_bc_matrix.
"""Cell Ranger adapter - run the aligner, then hand back the unfiltered matrix or refuse.

NOTHING IN THIS MODULE HAS BEEN RUN AGAINST CELL RANGER

Cell Ranger is installed nowhere this module could reach, so no line below has been executed
against the real tool. That is stated first because a reader who assumes otherwise will trust
the wrong things. What CAN be relied on is the split the module is built around: argument
construction, FASTQ-name matching, MatrixMarket header parsing and metrics-summary parsing are
pure functions over strings and take no tool, no filesystem and no dependencies, so they are
testable anywhere. Everything that depends on the tool's behaviour is asserted at run time and
refuses with the searched paths named. The two places where a version difference will bite are
called out below; both are argument construction, and both are parameters rather than
assumptions.

THE COMMAND

    cellranger count --id ID --transcriptome REF --fastqs DIR[,DIR...] --sample NAME \\
                     --create-bam=false --localcores N --localmem GB

`--create-bam` is Cell Ranger 8.0 syntax and is REQUIRED there - the flag has no default and
the run refuses without it. Releases up to 7.x do not accept it at all and spell the same
intent `--no-bam`. Emitting the wrong one is a hard failure rather than a silent one, which is
the good case, so `bam_style` selects between them instead of the module picking.

Cell Ranger writes its pipestance into the CURRENT WORKING DIRECTORY under `--id`, which is why
`work_dir` is a parameter and is handed to the executor as `cwd`. An `--id` directory that
already exists is refused, never removed: only step 7 of this pipeline removes anything, and a
half-finished pipestance is evidence about what happened, not litter.

WHY THE FASTQ NAMES ARE CHECKED BEFORE THE RUN

`--sample` selects input by filename. The convention is
`{sample}_S{n}_L{lane}_{R1,R2,I1,I2}_001.fastq.gz`, and a name that does not match it is not an
error to the caller - it is simply not selected. A `--sample` that matches nothing produces a
fast, clear failure from Cell Ranger itself; the case worth catching earlier is a directory
whose files are named by some other convention, because the run then either refuses after the
job has been scheduled or proceeds on a subset. `matching_fastqs` is pure and separates a strict
convention match from a loose prefix match, and the loose case is RECORDED in the returned
metrics rather than silently accepted.

WHAT THE RAW MATRIX IS

`<work_dir>/<id>/outs/raw_feature_bc_matrix` - every barcode, including the empties. It is the
one this pipeline consumes: an ambient model learns its background profile FROM the empty
droplets, and `filtered_feature_bc_matrix` no longer has them. The filtered matrix is located
and reported when present, so a caller can compare the two cell calls, but it is never returned
as the raw one.

DELIBERATE DUPLICATION

The MatrixMarket helpers here are near-copies of the ones in the CeleScope adapter. That is
intended: an adapter describes one tool's output tree, and sharing the description would mean a
change made for one aligner silently altering how the other's output is validated.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import re
import shutil
import sys
import time
from pathlib import Path, PurePosixPath, PureWindowsPath

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from engine.provenance import tool_version  # noqa: E402
from engine.task import TaskFailure  # noqa: E402

DEFAULT_EXE = "cellranger"
SUBCOMMAND = "count"

#: Relative to `<work_dir>/<id>`. The raw one is what this pipeline consumes.
RAW_MATRIX_SUBPATH = ("outs", "raw_feature_bc_matrix")
FILTERED_MATRIX_SUBPATH = ("outs", "filtered_feature_bc_matrix")
METRICS_SUMMARY_SUBPATH = ("outs", "metrics_summary.csv")
WEB_SUMMARY_SUBPATH = ("outs", "web_summary.html")

#: 8.0+ requires `--create-bam`; <= 7.x accepts only `--no-bam`. Neither is guessed.
BAM_STYLES = ("create-bam", "no-bam")

MATRIX_NAMES = ("matrix.mtx.gz", "matrix.mtx")
BARCODE_NAMES = ("barcodes.tsv.gz", "barcodes.tsv")
FEATURE_NAMES = ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")

#: The bcl2fastq / mkfastq naming convention `--sample` selects on.
FASTQ_CONVENTION = re.compile(
    r"^(?P<sample>.+)_S\d+(_L\d{3})?_(?P<read>[RI][12])_001\.fastq(\.gz)?$")

#: An `--id` may not carry a path separator, whitespace or a comma: the first two are rejected
#: by the tool and the third would be read as a list elsewhere on the command line.
ID_FORBIDDEN = ("/", "\\", ",", " ", "\t", "\n")


# ---------------------------------------------------------------------------------------------
# Unknown is not a value
# ---------------------------------------------------------------------------------------------
#: Type NAMES of the missing-value scalars that are neither None nor a float. Matched by name so
#: that the check costs nothing, and works, when pandas and numpy are not installed - this module
#: must stay importable on a bare clone, so neither may be imported at module scope.
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
    False against every threshold - indistinguishable, downstream, from a value that was measured
    and did not exceed the cut. That reads as a PASS, which is why this predicate exists once per
    module and every "is this unknown?" question in the file goes through it.

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
    second of start-up to a CLI that does not otherwise need it. `metrics_summary.csv` is read by
    pandas in many callers, and its blanks arrive here as whichever of those five shapes that
    caller's dtype backend produces.

    For a value that may be a numpy boolean the rule in this module is `bool(x)` AFTER
    `is_missing(x)` has been checked, never `x is True` - `numpy.bool_(True) is True` is False,
    so identity reads a genuinely flagged row as unflagged. See `is_true`.
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
    numpy-backed table; `bool(pandas.NA)` raises rather than returning False. The rule is
    therefore: refuse unknown, then `bool()`.
    """
    if is_missing(value):
        raise TaskFailure(
            f"{name} is {value!r}, which is not True or False. A flag that was never set is not "
            f"the same as one set to False, and reading it as False here would record a decision "
            f"nobody made.")
    return bool(value)


def _require_present(value, name: str, sample: str = "?"):
    """Refuse an unknown BEFORE it reaches `Path()`. Returns the value unchanged.

    `command_path` is where a path is spelled and where an unknown one is refused, but several
    checks legitimately run before it - a reference directory has to be proven to exist before
    the command is built - and each of those constructs a `Path` first. Two things happen there
    that are not a refusal: `Path(None)` raises a `TypeError` naming neither the argument nor the
    sample, and `Path("")` is the CURRENT DIRECTORY, so an unsupplied `work_dir` silently becomes
    "run here" - which for Cell Ranger means writing a pipestance into it. This makes the refusal
    arrive with the argument's name on it, at the top, either way.
    """
    if is_missing(value):
        raise TaskFailure(
            f"{sample}: {name} is {value!r}, which is not a path, and there is no default.\n"
            f"  Unknown is not a value here. Spelled out it becomes a directory named after the "
            f"sentinel\n  - `None`, `<NA>`, `nan`, `NaT` - and the failure then arrives as a "
            f"missing file at a path\n  nobody wrote, several messages away from the argument "
            f"that was never supplied.")
    return value


def _require_positive_int(value, name: str, sample: str = "?") -> int:
    """An integer >= 1, or refuse. `--localcores` and `--localmem` govern what the job can do."""
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

    `Path.is_absolute()` answers for the host running this code, which is not necessarily the
    host that will run the command: an orchestrator on Windows submitting to a POSIX cluster has
    `/refs/GRCh38` reported as relative, and refusing a correct run is how a gate gets switched
    off. A rooted-but-driveless path counts as anchored for the same reason - the hazard is a
    RELATIVE path resolving against the pipestance directory, which a rooted one does not do.
    """
    text = str(path)
    return (PurePosixPath(text).is_absolute() or PureWindowsPath(text).is_absolute()
            or bool(PureWindowsPath(text).root))


def _require_absolute(path, name: str) -> str:
    """Absolute paths only: the command runs with `cwd` set to the pipestance directory.

    A relative `--fastqs` resolves against the working directory rather than against the
    caller's, which does not error - it points somewhere else, and Cell Ranger reports it as
    finding no input for the sample.

    The value is returned as the TEXT it arrived as. Round-tripping through `Path` would rewrite
    `/refs/GRCh38` as `\\refs\\GRCh38` on a Windows orchestrator - a path the cluster cannot
    open, arriving as an error about missing input.

    An unknown is refused before `str()` is taken of it. Left to the absoluteness test it would
    be reported as "must be an absolute path, got 'None'", which reads as a caller who passed a
    relative path rather than as one who passed nothing.
    """
    _require_present(path, name)
    text = str(path)
    if not is_absolute_path(text):
        raise TaskFailure(
            f"{name} must be an absolute path, got {text!r}.\n"
            f"  Cell Ranger is launched with the working directory set to the run directory, so "
            f"a\n  relative path resolves against THAT and is reported as missing input rather "
            f"than as a\n  wrong path.")
    return text


def command_path(value, name: str, sample: str = "?") -> tuple:
    """`(text, how)` - the spelling a path must have on the command line, and how it was obtained.

    The one place a path may be rewritten, and it rewrites only what it must. An ALREADY ABSOLUTE
    path is returned as the caller's own text: `str(Path('/refs/GRCh38'))` is `\\refs\\GRCh38` on
    a Windows orchestrator, so a `Path` round trip is a silent rewrite into a path the POSIX
    cluster cannot open, and Cell Ranger reports it as finding no input rather than as a mangled
    argument. `_require_absolute` promises the value is used as the text it arrived as, and this
    is what makes that promise true of `run_cellranger` and not only of `build_command`.

    A RELATIVE path cannot be left alone - it would resolve against the pipestance directory - so
    it is resolved against this host's working directory, the only host that can do it, and `how`
    says so. The rewrite then appears in the run's recorded parameters instead of happening
    quietly.

    AN UNKNOWN IS REFUSED BEFORE ANY OF THAT HAPPENS, and the order is the whole point. This
    function used to take `str(value)` on its first line, which is where a sentinel stops being
    recognisable: `None` becomes the four characters `None`, `pandas.NA` becomes `<NA>`,
    `pandas.NaT` becomes `NaT`, and `float('nan')` becomes `nan`. Each is then a perfectly
    ordinary RELATIVE path, so the function resolved it against this host's working directory and
    handed back something like `<cwd>/None` with `how` reporting, accurately, that it had been
    resolved. An unsupplied `--transcriptome` reached the aligner as a directory named after the
    sentinel, and the run failed as a missing file at a path nobody had written - a message
    pointing at the filesystem instead of at the argument that was never given.

    Every argv element, `cwd` and recorded metric in this module comes through here, so this is
    the one check that has to be made in the right order.
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


def require_id(run_id, sample: str = "?") -> str:
    """A usable `--id`, or refuse.

    `--id` names the pipestance directory. A separator or a comma in it either fails inside the
    tool or, on the comma, is read as a list by a neighbouring argument - so it is checked here,
    where the message can say which character is the problem.
    """
    if is_missing(run_id):
        raise TaskFailure(f"{sample}: --id is required and has no default")
    text = str(run_id).strip()
    bad = [c for c in ID_FORBIDDEN if c in text]
    if bad:
        shown = ", ".join(repr(c) for c in bad)
        raise TaskFailure(
            f"{sample}: --id {text!r} contains {shown}.\n"
            f"  --id names a directory created in the working directory; a separator, a space "
            f"or a\n  comma there is either rejected by the tool or parsed as something else.")
    return text


# ---------------------------------------------------------------------------------------------
# FASTQ selection - pure, testable on names alone
# ---------------------------------------------------------------------------------------------
def matching_fastqs(names, sample: str) -> dict:
    """Split file names into those `--sample` selects strictly, loosely, or not at all. Pure.

    `strict` follows the mkfastq convention `{sample}_S1_L001_R1_001.fastq.gz`, which is what
    `--sample` matches on. `loose` starts with `{sample}_` but does not follow it - the shape a
    renamed public download usually has, and one Cell Ranger will not select. `other` is
    everything else in the directory, kept so a failure message can print what WAS there instead
    of only what was missing.

    An unknown `sample` is refused rather than stringified: `str(None)` would select files named
    `None_S1_L001_R1_001.fastq.gz`, find none, and report "no FASTQ matches --sample 'None'" -
    a message about the directory's contents when the fault is that nothing was passed.
    """
    if is_missing(sample):
        raise TaskFailure(f"--sample is {sample!r}; it selects the input files and has no default")
    sample = str(sample)
    out = {"strict": [], "loose": [], "other": []}
    for raw in names:
        name = str(raw)
        m = FASTQ_CONVENTION.match(name)
        if m and m.group("sample") == sample:
            out["strict"].append(name)
        elif name.startswith(sample + "_"):
            out["loose"].append(name)
        else:
            out["other"].append(name)
    for key in out:
        out[key].sort()
    return out


def scan_fastq_dirs(fastq_dirs, sample: str) -> dict:
    """`matching_fastqs` over every `--fastqs` directory, merged. Reads names only, not files."""
    merged = {"strict": [], "loose": [], "other": []}
    for i, d in enumerate(fastq_dirs):
        _require_present(d, f"--fastqs[{i}]", str(sample))
        path = Path(d)
        if not path.is_dir():
            raise TaskFailure(f"{sample}: --fastqs directory does not exist: {path}")
        found = matching_fastqs(sorted(p.name for p in path.iterdir() if p.is_file()), sample)
        for key in merged:
            merged[key].extend(f"{path.name}/{n}" for n in found[key])
    return merged


# ---------------------------------------------------------------------------------------------
# Argument construction - pure, testable without the tool
# ---------------------------------------------------------------------------------------------
def build_command(run_id, transcriptome, fastqs, sample, localcores, localmem_gb,
                  exe: str = DEFAULT_EXE, create_bam: bool = False,
                  bam_style: str = "create-bam", extra=()) -> list:
    """The `cellranger count` argv. Pure - no filesystem, no tool, no environment.

    `fastqs` is a sequence of directories and is joined with commas, which is how Cell Ranger
    takes more than one. A directory whose own path contains a comma is refused: it would split
    into two paths that do not exist, and the error names neither of them.

    `bam_style` picks between the two spellings of "do not write a BAM" - `--create-bam=false`
    on 8.0+, where the flag is mandatory, and `--no-bam` on 7.x and earlier, which does not know
    the newer one. The module does not detect the release, because detecting it would mean
    running the tool to decide how to run the tool.

    `extra` is appended verbatim for flags this adapter does not model. Verbatim does not extend
    to unknown: an element that is `None` or `pandas.NA` would be spelled onto the command line
    as `None` or `<NA>` and read by the tool as a value, so it is refused with its position named.
    """
    rid = require_id(run_id, str(sample))
    if is_missing(sample):
        raise TaskFailure("--sample is required and has no default")
    if is_missing(exe):
        raise TaskFailure("the cellranger executable is required and has no default")
    ref = _require_absolute(transcriptome, "--transcriptome")
    cores = _require_positive_int(localcores, "--localcores", str(sample))
    mem = _require_positive_int(localmem_gb, "--localmem", str(sample))

    if isinstance(fastqs, (str, bytes, Path)):
        fastqs = [fastqs]
    dirs = [_require_absolute(f, "--fastqs") for f in fastqs]
    if not dirs:
        raise TaskFailure("--fastqs is required: at least one directory, and no default exists")
    for d in dirs:
        if "," in d:
            raise TaskFailure(
                f"--fastqs path contains a comma: {d!r}\n"
                f"  Cell Ranger separates multiple directories with commas, so this would be "
                f"split into\n  two paths that do not exist.")

    if str(bam_style) not in BAM_STYLES:
        raise TaskFailure(
            f"bam_style={bam_style!r} is not one of {BAM_STYLES}.\n"
            f"  Cell Ranger 8.0+ requires --create-bam and refuses without it; 7.x and earlier "
            f"accept\n  only --no-bam. There is no spelling that works on both, and this adapter "
            f"does not guess\n  the release.")

    cmd = [str(exe), SUBCOMMAND,
           "--id", rid,
           "--transcriptome", ref,
           "--fastqs", ",".join(dirs),
           "--sample", str(sample).strip()]
    if str(bam_style) == "create-bam":
        cmd += ["--create-bam=" + ("true" if create_bam else "false")]
    elif not create_bam:
        cmd += ["--no-bam"]
    cmd += ["--localcores", str(cores), "--localmem", str(mem)]
    extra = list(extra)
    for i, item in enumerate(extra):
        if is_missing(item):
            raise TaskFailure(
                f"extra[{i}] is {item!r}, which is not a command-line argument. str() would place "
                f"`{item!s}` on the argv and Cell Ranger would read it as a value.")
    return cmd + [str(x) for x in extra]


def build_env(env_bin=None, base_path=None, pathsep=None, extra=None) -> dict:
    """The environment overlay for the executor, with a tool directory FIRST on PATH.

    Cell Ranger ships its own dependencies and does not have CeleScope's bare-name STAR problem,
    so this is here for the ordinary case of an installation that is not on the login shell's
    PATH rather than to close a specific hazard.

    `pathsep` is explicit for the same reason as in the CeleScope adapter: `os.pathsep` is `;`
    on Windows and would build a PATH no POSIX cluster can read when an orchestrator submits
    across platforms.

    A `PATH` inside `extra` is MERGED, never overwritten and never discarded: it becomes the base
    that `env_bin` is prepended to. Either half winning outright loses the other silently, since
    a PATH is still present and still plausible afterwards.
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

    None means NOT FOUND and is reported as its own outcome; it never stands in for a path. The
    lookup happens on the host building the command, which for a scheduler is the submit host.
    """
    env = build_env(env_bin, base_path=base_path, pathsep=pathsep)
    search = env.get("PATH") if "PATH" in env else base_path
    return {str(n): shutil.which(str(n), path=search) for n in names}


# ---------------------------------------------------------------------------------------------
# Output location and parsing
# ---------------------------------------------------------------------------------------------
def matrix_triplet(directory) -> dict:
    """The three MatrixMarket files in a directory, or `{}` if the set is not complete.

    Deliberately all-or-nothing: two files out of three is not most of a matrix, it is a
    directory that fails to load, and returning a partial dict invites a caller to use it.
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


def find_raw_matrix(work_dir, run_id: str) -> Path:
    """`<work_dir>/<id>/outs/raw_feature_bc_matrix`, checked to hold a matrix, or a refusal.

    A directory that exists and holds no triplet is reported differently from one that does not
    exist, because the first is the failure this repository keeps meeting: a step that returned
    zero and wrote nothing, discovered several steps later as a confusing error about a
    different file.
    """
    root = Path(work_dir) / str(run_id)
    raw = root.joinpath(*RAW_MATRIX_SUBPATH)
    if matrix_triplet(raw):
        return raw
    outs = root / "outs"
    if raw.is_dir():
        present = sorted(p.name for p in raw.iterdir())
        raise TaskFailure(
            f"{run_id}: {raw} exists but holds no complete matrix triplet.\n"
            f"  contains: {', '.join(present) or '(empty)'}\n"
            f"  A triplet is one of {MATRIX_NAMES} plus one of {BARCODE_NAMES} plus one of "
            f"{FEATURE_NAMES}.")
    if outs.is_dir():
        present = sorted(p.name for p in outs.iterdir())
        raise TaskFailure(
            f"{run_id}: {raw} does not exist.\n"
            f"  {outs} contains: {', '.join(present) or '(empty)'}\n"
            f"  The RAW matrix is what this pipeline consumes; filtered_feature_bc_matrix has "
            f"already\n  been cell-called and no longer holds the empty droplets an ambient "
            f"model learns from.")
    raise TaskFailure(
        f"{run_id}: no outs/ directory under {root}.\n"
        f"  Cell Ranger writes its pipestance into the working directory under --id. Either the "
        f"run\n  did not reach the output stage, or it ran with a different working directory "
        f"than the one\n  searched here.")


def parse_mtx_header(lines, name: str = "<matrix>") -> dict:
    """Dimensions from a MatrixMarket header. Pure - it takes lines, so it is testable on text.

    Rows are returned as `n_genes` and columns as `n_barcodes`, the features-by-barcodes
    orientation Cell Ranger writes. A zero dimension or zero stored entries is refused rather
    than returned: an empty matrix loads, plots and summarises without complaint.
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
                f"entries.")
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

    A final line with no trailing newline is counted, because under-counting by one would look
    like a one-row disagreement between a matrix header and its features file - the size of
    near-miss that gets explained away rather than investigated.
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


def verify_matrix_dims(triplet: dict, header: dict, run_id: str = "?",
                       check_barcodes: bool = True) -> dict:
    """Cross-check the header's dimensions against the two label files, and refuse a mismatch.

    The same two numbers by two independent routes. A truncated transfer, a half-written file
    and a transposed matrix all produce a header that parses cleanly on its own. `check_barcodes`
    is a parameter because a raw droplet matrix names millions of barcodes; when it is skipped
    the outcome is recorded as `not counted` rather than as a number.
    """
    n_features = count_lines(triplet["features"])
    if n_features != header["n_genes"]:
        raise TaskFailure(
            f"{run_id}: matrix header declares {header['n_genes']:,} genes but "
            f"{triplet['features'].name} has {n_features:,} rows.\n"
            f"  These must agree; they do not when a file was truncated, when the matrix is "
            f"transposed,\n  or when the two files came from different runs.")
    out = {"features_rows": n_features}
    if not check_barcodes:
        out["barcodes_rows"] = "not counted"
        return out
    n_barcodes = count_lines(triplet["barcodes"])
    if n_barcodes != header["n_barcodes"]:
        raise TaskFailure(
            f"{run_id}: matrix header declares {header['n_barcodes']:,} barcodes but "
            f"{triplet['barcodes'].name} has {n_barcodes:,} rows.")
    out["barcodes_rows"] = n_barcodes
    return out


def parse_metrics_summary(text: str) -> dict:
    """`outs/metrics_summary.csv` as `{"fields": {...}, "blank": [...]}`. Pure - takes the text.

    Values are returned as the STRINGS Cell Ranger wrote. They arrive thousands-separated and
    percent-suffixed (`"1,234"`, `"91.2%"`), and converting them here would mean deciding what a
    blank cell means - which is the decision this pipeline refuses to make quietly. Blank fields
    are listed by name instead of being coerced to zero, so a consumer must handle never-measured
    as its own category.
    """
    rows = list(csv.DictReader(io.StringIO(text)))
    if len(rows) != 1:
        raise TaskFailure(
            f"metrics_summary.csv has {len(rows)} data rows; exactly one was expected. "
            f"A multi-row summary is a different pipeline's output and its columns do not mean "
            f"the same thing.")
    row = rows[0]
    fields, blank = {}, []
    for key, value in row.items():
        if key is None:
            continue
        name = str(key).strip()
        if is_missing(value):
            blank.append(name)
            continue
        fields[name] = str(value).strip()
    return {"fields": fields, "blank": sorted(blank)}


# ---------------------------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------------------------
def version(exe: str = DEFAULT_EXE, args=("--version",)) -> str:
    """The observed version string, or `not invoked` - never inferred from a path or a module.

    The release matters here more than for most tools, because `--create-bam` and `--no-bam`
    divide at 8.0. This records what answered; it does not choose the flag from it, since the
    version is observed on the host that builds the command and the job may run on another.
    """
    return tool_version(exe, tuple(args))


# ---------------------------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------------------------
def run_cellranger(sample: str, fastq_dirs, transcriptome, work_dir, log,
                   *, run_id=None, localcores, localmem_gb, env_bin=None,
                   exe: str = DEFAULT_EXE, create_bam: bool = False,
                   bam_style: str = "create-bam", pathsep=None,
                   check_barcode_count: bool = True, require_strict_fastq_names: bool = False,
                   timeout_s=None, extra=(), executor=None) -> dict:
    """Run `cellranger count` for one sample and return the raw matrix, or fail with the search.

    Nothing about the sample, the reference or the resources is decided here - the orchestrator
    decides and this executes. `run_id` defaults to the sample name because the pipestance
    directory is an execution detail rather than a scientific choice; every parameter that
    changes what the matrix IS is required.

    `localcores` and `localmem_gb` are keyword arguments without defaults. A silently defaulted
    core count wastes a scheduler allocation, and a silently defaulted memory cap is worse: Cell
    Ranger treats `--localmem` as a ceiling and a job that exceeds it is killed in a way that
    reads as a crash rather than as a resource setting.

    An existing `--id` directory is REFUSED, not removed. Resume belongs to `engine/state.py`,
    which decides whether a task needs to run at all; deleting a previous pipestance here would
    destroy the evidence of what went wrong with it.

    THAT REFUSAL IS ALSO WHAT MAKES THE OUTPUT CHECKS MEAN ANYTHING. Finding a matrix after the
    command returns proves a matrix is there, not that this invocation wrote it; a run that exits
    zero having written nothing passes every such check whenever a previous pipestance is lying
    at the path. Because the pipestance is proven ABSENT before the command starts, and every
    output this function reports lives inside it, anything found afterwards was written by this
    invocation. The mtimes are recorded as corroboration only - they come from the compute node's
    clock while the refusal is made on the orchestrator's, and a freshness gate that fires on
    clock skew is a gate someone switches off.

    Paths reach the command line as the caller's own text when that text is absolute; only a
    relative path is resolved, and `metrics["path_resolution"]` says which and against what. See
    `command_path`.

    Returns `{"outputs": [...], "metrics": {...}, "versions": {...}}`. `outputs` lists only files
    checked to exist after the command returned; the raw matrix DIRECTORY is reported as
    `metrics["raw_matrix"]`.
    """
    if executor is None:
        raise TaskFailure("run_cellranger requires an executor; there is no in-process fallback")
    if is_missing(sample):
        raise TaskFailure("sample is required and has no default")
    sample = str(sample).strip()
    rid = require_id(sample if is_missing(run_id) else run_id, sample)
    want_bam = is_true(create_bam, "create_bam")
    count_barcodes = is_true(check_barcode_count, "check_barcode_count")
    strict_names = is_true(require_strict_fastq_names, "require_strict_fastq_names")

    # Every path argument is tested for unknown before it reaches `Path()`. `command_path` refuses
    # these too, but it runs AFTER the existence checks below, and `Path(None)` raises a TypeError
    # there that names neither the argument nor the sample while `Path("")` is the current
    # directory - which for `work_dir` means writing a pipestance into wherever the orchestrator
    # happens to be. The refusal has to carry the argument's name to be worth anything.
    for _value, _name in ((transcriptome, "--transcriptome"), (work_dir, "work_dir"),
                          (log, "log")):
        _require_present(_value, _name, sample)

    ref = Path(transcriptome)
    if not ref.is_dir():
        raise TaskFailure(
            f"{sample}: --transcriptome does not exist: {ref}\n"
            f"  It is a Cell Ranger reference directory, not a FASTA or a GTF.")
    ref_text, ref_route = command_path(transcriptome, "--transcriptome", sample)

    if isinstance(fastq_dirs, (str, bytes, Path)):
        fastq_dirs = [fastq_dirs]
    dirs, dir_texts, dir_routes = [], [], {}
    for i, d in enumerate(fastq_dirs):
        _require_present(d, f"--fastqs[{i}]", sample)
        p = Path(d)
        if not p.is_dir():
            raise TaskFailure(f"{sample}: --fastqs directory does not exist: {p}")
        dirs.append(p)
        text, route = command_path(d, "--fastqs", sample)
        dir_texts.append(text)
        dir_routes[text] = route
    if not dirs:
        raise TaskFailure(f"{sample}: --fastqs is required and no directory was given")

    found = scan_fastq_dirs(dirs, sample)
    if not found["strict"] and not found["loose"]:
        listing = "\n".join(f"    {n}" for n in found["other"][:20]) or "    (no files)"
        raise TaskFailure(
            f"{sample}: no FASTQ in {', '.join(str(d) for d in dirs)} matches --sample "
            f"{sample!r}.\n"
            f"  --sample selects by filename: {{sample}}_S1_L001_R1_001.fastq.gz. The "
            f"directories hold:\n{listing}\n"
            f"  A --sample that matches nothing is caught here rather than after the job has "
            f"been\n  scheduled.")
    naming = "convention" if found["strict"] else "non-standard"
    if naming == "non-standard" and strict_names:
        listing = "\n".join(f"    {n}" for n in found["loose"][:20])
        raise TaskFailure(
            f"{sample}: files begin with {sample!r} but none follows the "
            f"{{sample}}_S1_L001_R1_001.fastq.gz convention that --sample selects on:\n{listing}")

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    work_text, work_route = command_path(work_dir, "work_dir", sample)
    pipestance = work / rid
    if pipestance.exists():
        raise TaskFailure(
            f"{sample}: {pipestance} already exists.\n"
            f"  Cell Ranger refuses to write into an existing --id directory, and this adapter "
            f"does not\n  delete it: a previous pipestance is the record of what happened to "
            f"that run. Move it\n  aside, or pass a different run_id. Whether the step needs to "
            f"run at all is a resume\n  question and belongs to the run manifest, not here.")

    cmd = build_command(rid, ref_text, dir_texts, sample, localcores, localmem_gb, exe=exe,
                        create_bam=want_bam, bam_style=bam_style, extra=extra)
    env = build_env(env_bin, pathsep=pathsep)
    resolved = resolve_on_path([str(exe)], env_bin, pathsep=pathsep)
    exe_path = resolved[str(exe)]

    started = time.time()
    executor.shell(cmd, log=Path(log), env=env, cwd=work_text, timeout_s=timeout_s)

    raw_dir = find_raw_matrix(work, rid)
    triplet = matrix_triplet(raw_dir)
    if not triplet:
        raise TaskFailure(f"{sample}: {raw_dir} lost its matrix triplet between checks")
    header = read_mtx_header(triplet["matrix"])
    dims = verify_matrix_dims(triplet, header, rid, check_barcodes=count_barcodes)

    outputs = [triplet["matrix"], triplet["barcodes"], triplet["features"]]
    metrics_path = pipestance.joinpath(*METRICS_SUMMARY_SUBPATH)
    if metrics_path.is_file():
        summary = parse_metrics_summary(metrics_path.read_text(encoding="utf-8", errors="replace"))
        summary_status = "parsed"
        outputs.append(metrics_path)
    else:
        # Recorded as its own outcome. Absent is not empty, and neither is zero.
        summary = {"fields": {}, "blank": []}
        summary_status = "absent"

    web = pipestance.joinpath(*WEB_SUMMARY_SUBPATH)
    if web.is_file():
        outputs.append(web)
    filtered = pipestance.joinpath(*FILTERED_MATRIX_SUBPATH)

    for p in outputs:
        if not Path(p).is_file():
            raise TaskFailure(f"{sample}: promised output is absent after the run: {p}")
    # Each of these sits inside the pipestance, which was proven absent before the command ran, so
    # its presence is proof of authorship rather than of luck. The ages are corroboration and are
    # recorded rather than enforced; see the docstring on why a clock is not the proof here.
    output_ages = {}
    for p in outputs:
        try:
            output_ages[Path(p).name] = round(started - Path(p).stat().st_mtime, 3)
        except OSError:
            output_ages[Path(p).name] = None

    metrics = {
        "sample": sample,
        "run_id": rid,
        "transcriptome": ref_text,
        "fastqs": list(dir_texts),
        "path_resolution": {"--transcriptome": ref_route, "work_dir": work_route, **dir_routes},
        "pipestance_absent_before_run": True,
        "output_age_s_at_start": output_ages,
        "started_epoch_s": started,
        "freshness_proof": ("the pipestance was proven absent before the command ran, so every "
                            "output found afterwards was written by this invocation; the ages "
                            "above are corroboration from the compute node's clock, not the "
                            "proof"),
        "fastq_naming": naming,
        "fastq_files_matched": len(found["strict"]) + len(found["loose"]),
        "localcores": _require_positive_int(localcores, "--localcores", sample),
        "localmem_gb": _require_positive_int(localmem_gb, "--localmem", sample),
        "bam_style": str(bam_style),
        "create_bam": want_bam,
        "work_dir": work_text,
        "pipestance": str(pipestance),
        "raw_matrix": str(raw_dir),
        "filtered_matrix": str(filtered) if matrix_triplet(filtered) else "absent",
        "n_genes": header["n_genes"],
        "n_barcodes": header["n_barcodes"],
        "n_entries": header["n_entries"],
        "features_rows": dims["features_rows"],
        "barcodes_rows": dims["barcodes_rows"],
        "metrics_summary": summary_status,
        "metrics_summary_fields": summary["fields"],
        "metrics_summary_blank": summary["blank"],
        "env_bin": str(env_bin) if not is_missing(env_bin) else "not set",
        "log": str(Path(log)),
    }
    versions = {"cellranger": version(exe_path if exe_path else str(exe))}
    return {"outputs": outputs, "metrics": metrics, "versions": versions}
