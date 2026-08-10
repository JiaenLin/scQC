# Execution adapter for CellBender remove-background: builds the command, runs it through an
# executor, verifies that it actually produced something, and measures what it removed.
# It removes no observation of its own and writes nothing into its inputs. Every number it
# returns is read back off a file that was checked to exist.
"""Step 1a - run the denoiser, then prove it ran.

WHAT THIS FILE IS FOR, AND WHAT IT DELIBERATELY IS NOT

`modules/01_ambient/` decides whether CellBender runs (`ambient.py`), whether its learning rate
has to move (`lr_policy.py`) and whether the result is fit to use (`audit_ambient.py`). None of
those touch a file. This adapter is the other half: it turns a decision into a process and a
process back into numbers, and it takes no decisions itself. There is no policy here, no
samplesheet, no default cohort, no threshold. Paths and parameters arrive as arguments and the
caller owns every one of them.

A ZERO EXIT IS NOT SUCCESS. THE HAZARD THIS ADAPTER EXISTS FOR

CellBender 0.3.2 paired with some pyro/torch builds cannot write its checkpoint, and its
posterior step declines to run without one. The observable result is a job that trains for its
full epoch budget, logs nothing that reads as fatal, exits 0, and leaves no output matrix. Every
scheduler, every wrapper and every resume manifest records that run as a success; the next step
then fails on a missing file, three hours and one queue wait later, with an error naming the
wrong thing.

So `run_remove_background()` treats the exit code as necessary and not sufficient. After the
command returns it checks that the output `.h5` exists, that it begins with the HDF5 signature,
that it is not trivially small, and - where h5py is importable - that it contains a count matrix
with at least one barcode and at least one stored value. A failure of any of those raises
`TaskFailure` naming the checkpoint/posterior failure mode, because the message a stranger reads
at 2am is part of the code.

THE REPORT STEP SHELLS OUT BY BARE NAME

The end of `remove-background` renders its HTML report by invoking `jupyter nbconvert` from
`PATH`, with no check on the result. Under a scheduler, or from a wrapper that trimmed the
environment, `jupyter` is frequently not on `PATH` even though CellBender itself is - so the
report silently does not appear while the run reports success. `build_env()` puts the
environment's own `bin` directory at the FRONT of `PATH` for exactly this, and the returned
metrics record whether the report was written rather than assuming either answer.

WHAT "FRACTION REMOVED" MEANS HERE, AND WHY IT NEEDS THE RAW MATRIX

The denoised object does not carry the counts it started from, so no function in this file can
report a removal fraction from the CellBender output alone. `parse_metrics()` therefore requires
the raw input as well and compares the two over ONE explicit barcode scope - by default the
called cells. Computed over all droplets the same quantity is dominated by the empties, which are
reduced to near nothing by design, and it is then not the quantity the cohort ranges in
`modules/01_ambient/ambient.py` describe. The scope is an argument, it is recorded in the
returned metrics, and it must not vary within a cohort.

COLUMN NAMES ARE NOT NEGOTIABLE

`summary_row()` and `per_gene_rows()` emit exactly the columns `audit_ambient.audit()` reads -
`sample`, `fraction_removed_overall`, `genes_fully_removed`; and `sample`, `symbol`,
`fraction_removed`, `raw_detection_frac`, `denoised_detection_frac`. Extra columns travel
alongside them because a table a human opens should answer more than the gate asks, but the
named five are produced by this file and consumed by that one, and renaming either end silently
turns a check into a no-op.

UNKNOWN IS NOT A VALUE (docs/PRINCIPLES.md section 4)

Nothing here substitutes zero for a measurement that was not taken. A gene with no counts in the
raw matrix has no removal fraction - not 0.0, not 1.0 - so it is excluded from the per-gene table
and counted separately, and no NaN is ever written into a column a gate reads. A diagnostic that
could not be derived from the artifact it was asked for is ABSENT from the returned dict rather
than present as None, so that a caller merging two dicts cannot have a real measurement
overwritten by a blank one.

`is_missing()` is the single predicate that decides what "unknown" means here, and it covers
every shape one arrives in - None, float nan, numpy nan, pandas.NA, pandas.NaT, numpy.ma.masked -
because a sentinel that slips through is compared with `>=` or passed to `int()` and reads as a
value that failed the test, which is a PASS. Guarding on the PRESENCE of a key rather than on its
VALUE is the same defect wearing a different hat: `{"fraction_removed_overall": None}` satisfied
`if required not in metrics` and put the None straight into the audit table.

A PREVIOUS RUN'S OUTPUT IS NOT THIS RUN'S RESULT

`verify_products()` proves a usable matrix is at the path. It cannot prove this invocation put it
there, and the checkpoint/posterior failure above produces exactly a run that exits 0 and writes
nothing - so with a previous run's `.h5` in place, every check passes and the old numbers are
recorded under the new parameters. `run_remove_background()` therefore DELETES every declared
product before it launches the command, so that what exists afterwards can only have been written
by this invocation. See its docstring for what is deleted and what is deliberately not.
"""

from __future__ import annotations

import csv
import os
import re
import sys
import time
from pathlib import Path
from typing import Sequence

try:
    from engine.task import TaskFailure
except ImportError:  # loaded by path from outside the repository root
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from engine.task import TaskFailure
from engine.provenance import NOT_INVOKED, tool_version

# ---------------------------------------------------------------------------- constants

SUBCOMMAND = "remove-background"

#: CellBender's package default for the target false positive rate, restated here only so that
#: the emitted command line is explicit. `modules/01_ambient/ambient.py` owns the policy.
DEFAULT_FPR = 0.0

#: CellBender 0.3.2 selects the CPU by the ABSENCE of `--cuda`; it has no `--cpu-only` flag, and
#: passing one aborts the run in argparse before a single epoch. `device="cpu"` therefore emits
#: no device flag at all. A site build that does accept an explicit flag can pass it through
#: `cpu_flags=("--cpu-only",)` - the argument exists so that the choice is visible in the caller
#: rather than guessed here.
CPU_ONLY_FLAG = "--cpu-only"

HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"

#: A denoised matrix for a solid-tissue library is tens of megabytes. Anything below this is a
#: stub, a truncated write or an empty container. It is an argument on every function that uses
#: it because a deliberately tiny test fixture is a legitimate input.
MIN_OUTPUT_BYTES = 1_048_576

#: Barcodes whose posterior cell probability reaches this are cells. CellBender's own
#: `_cell_barcodes.csv` is written at 0.5; the value is an argument so that the two routes can be
#: compared against each other rather than one being assumed to match the other.
CELL_PROBABILITY_THRESHOLD = 0.5

#: How many trailing epochs `convergence_indicator()` measures the tail movement over.
LEARNING_CURVE_WINDOW = 20

#: Log lines that indicate the checkpoint/posterior failure mode described in the module
#: docstring. The list is a HEURISTIC over messages seen in the wild and absence of a match is
#: not evidence that the run succeeded - the file checks in `verify_products()` are the
#: authority, and this only makes the resulting message specific.
_CHECKPOINT_TROUBLE = re.compile(
    r"(unable to save checkpoint|failed to save checkpoint|could not save checkpoint|"
    r"cannot save checkpoint|no checkpoint found|checkpoint file not found|"
    r"skipping posterior|posterior .{0,40}(not|cannot|unable)|"
    r"checkpoint.{0,40}(unable|fail|error|not saved))",
    re.I,
)

_NBCONVERT_TROUBLE = re.compile(
    r"(nbconvert.{0,80}(not found|no such file|error|traceback)|"
    r"(command not found|not found).{0,40}nbconvert|jupyter: command not found)",
    re.I,
)

#: `[epoch 042]  average training loss: 1234.5678`, with `average test loss:` on test epochs.
_EPOCH_LINE = re.compile(
    r"\[epoch\s+(\d+)\]\s*(.*)$", re.I)
_TRAIN_LOSS = re.compile(r"average\s+training\s+loss:\s*(-?[\d.]+(?:[eE][-+]?\d+)?)", re.I)
_TEST_LOSS = re.compile(r"average\s+test\s+loss:\s*(-?[\d.]+(?:[eE][-+]?\d+)?)", re.I)
_VERSION_TOKEN = re.compile(r"\d+\.\d+(?:\.\d+)?(?:[.\-+][0-9A-Za-z.\-+]+)?")
_LOG_VERSION = re.compile(r"cellbender[^\n]{0,80}?v?(\d+\.\d+(?:\.\d+)?)", re.I)

#: Dataset basenames CellBender may use for each quantity. Every lookup is by EXACT basename
#: against this list and an unrecognised layout raises rather than being guessed at, because a
#: latent read from the wrong array is a number that looks measured and is not.
_CELL_PROB_NAMES = ("latent_cell_probability", "cell_probability")
_LATENT_INDEX_NAMES = ("barcode_indices_for_latents", "barcodes_analyzed_inds")
_TRAIN_ELBO_NAMES = ("learning_curve_train_elbo", "train_elbo")
_TEST_ELBO_NAMES = ("learning_curve_test_elbo", "test_elbo")
_TRAIN_EPOCH_NAMES = ("learning_curve_train_epoch", "train_epoch")
_TEST_EPOCH_NAMES = ("learning_curve_test_epoch", "test_epoch")
_CONVERGENCE_NAMES = ("convergence_indicator",)

#: Field names in a CellBender metrics table that carry an overall removal fraction. Matched
#: exactly; no fuzzy matching, because a field that merely sounds right is how a run gets audited
#: against the wrong number.
_METRICS_FRACTION_REMOVED = ("fraction_counts_removed", "fraction_removed",
                             "overall_fraction_counts_removed")

_MATRIX_KEYS = frozenset({"data", "indices", "indptr", "shape", "barcodes"})

SUMMARY_COLUMNS = ("sample", "fraction_removed_overall", "genes_fully_removed",
                   "cells_called", "droplets_total_raw", "droplets_in_output",
                   "total_counts_raw_in_scope", "total_counts_denoised_in_scope",
                   "genes_total", "scope", "cell_call_source")

PER_GENE_COLUMNS = ("sample", "symbol", "gene_id", "fraction_removed",
                    "raw_detection_frac", "denoised_detection_frac",
                    "raw_counts", "denoised_counts")


# ---------------------------------------------------------------------------- unknown is not a value

#: Type NAMES of the missing-value scalars that are neither None nor a float. Matched by name so
#: that the check costs nothing, and works, when pandas and numpy are not installed - this module
#: must stay importable with no third-party package at module scope.
_MISSING_TYPE_NAMES = frozenset({"NAType", "NaTType", "MaskedConstant"})


def is_missing(value) -> bool:
    """True when a value carries no information, in every shape one actually arrives in.

    The ONE predicate this module asks that question with. `None` is the shape everyone remembers
    and the rarest in practice: a metrics table read by pandas hands back `float('nan')` through
    the numpy dtypes, `pandas.NA` through the nullable or pyarrow-backed ones, `pandas.NaT` for a
    parsed date, and `numpy.ma.masked` through a masked array. None of those is `None`, only the
    first is a `float`, and each survives an `is not None` guard, then compares False against
    every threshold, then reaches `int()` or a CSV cell. Downstream that is indistinguishable
    from a value that was measured and did not exceed the cut - it reads as a PASS.

    Four routes, cheapest first:

      * identity against `None`;
      * blank and whitespace-only text, in `str` and in `bytes`;
      * the type NAME, which catches `pandas.NA`, `pandas.NaT` and `numpy.ma.masked` without
        importing anything;
      * `value != value`, which catches `float('nan')`, `numpy.float64('nan')`,
        `numpy.float32('nan')` - not a `float` subclass - and `numpy.datetime64('NaT')`.

    pandas is then consulted, but only if it is ALREADY in `sys.modules`: a value cannot be a
    pandas scalar in a process that never imported pandas, so looking there is both sufficient
    and free, whereas importing pandas inside a predicate that runs per gene would cost a second
    of start-up to a CLI that does not otherwise need it.

    For a value that may be a numpy boolean the rule here is `bool(x)` AFTER `is_missing(x)` has
    been checked, never `x is True`: `numpy.bool_(True) is True` is False, so identity reads a
    genuinely flagged row as unflagged.
    """
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip()
    if isinstance(value, bytes):
        return not value.decode("utf-8", "replace").strip()
    if type(value).__name__ in _MISSING_TYPE_NAMES:
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

    `numpy.bool_(True) is True` is False and `bool(pandas.NA)` raises, so the rule is: refuse
    unknown, then `bool()`.
    """
    if is_missing(value):
        raise TaskFailure(
            f"{name} is {value!r}, which is not True or False. A flag that was never set is not "
            f"the same as one set to False, and reading it as False here would record a decision "
            f"nobody made.")
    return bool(value)


# ---------------------------------------------------------------------------- command building


def _fmt_num(value) -> str:
    """Render a number for a command line without turning 0.0 into an unfamiliar string."""
    if isinstance(value, bool):
        raise TaskFailure(f"a boolean is not a numeric CellBender argument: {value!r}")
    if isinstance(value, int):
        return str(value)
    f = float(value)
    if f == int(f) and abs(f) < 1e16:
        return str(int(f))
    return repr(f)


def _seq(name: str, value) -> tuple:
    """Reject a bare string where a sequence of arguments is expected.

    `extra_args="--debug"` iterates into eleven single-character arguments and CellBender then
    fails on `-`, which reads as a CellBender problem rather than a caller's.
    """
    if value is None:
        return ()
    if isinstance(value, (str, bytes, Path)):
        raise TaskFailure(
            f"{name} must be a sequence of separate arguments, not {value!r}. A string is "
            f"iterated one character at a time and produces a nonsensical command line.")
    return tuple(str(v) for v in value)


def build_command(
    exe: str | Path,
    input_path: str | Path,
    output_h5: str | Path,
    *,
    device: str = "cuda",
    fpr: float | None = DEFAULT_FPR,
    learning_rate: float | None = None,
    epochs: int | None = None,
    expected_cells: int | None = None,
    total_droplets_included: int | None = None,
    low_count_threshold: int | None = None,
    checkpoint_mins: float | None = None,
    cpu_threads: int | None = None,
    posterior_batch_size: int | None = None,
    exclude_feature_types: Sequence[str] = (),
    cpu_flags: Sequence[str] = (),
    extra_args: Sequence[str] = (),
) -> list:
    """The `cellbender remove-background` argument vector, as a list, touching no filesystem.

    Pure so that the command a cohort will run can be reviewed, diffed and unit-tested without a
    GPU, a matrix or an installed CellBender. Every optional argument is OMITTED when it is None,
    so the default vector is exactly the reference invocation and any departure from the package
    defaults is visible in the emitted list rather than implied by it.

        cellbender remove-background --cuda --input RAW --output OUT.h5 --fpr 0

    `device` is "cuda" or "cpu". CellBender 0.3.2 has no CPU flag: the CPU is what it uses when
    `--cuda` is absent. See CPU_ONLY_FLAG.
    """
    if device not in ("cuda", "cpu"):
        raise TaskFailure(
            f"device must be 'cuda' or 'cpu', got {device!r}. There is no default: running a "
            f"cohort partly on one and partly on the other makes the denoising a technical "
            f"property that varies across the design.")
    cpu_flags = _seq("cpu_flags", cpu_flags)
    if device == "cuda" and cpu_flags:
        raise TaskFailure(
            f"device='cuda' was requested together with cpu_flags={cpu_flags!r}. These "
            f"contradict each other; one of them is not what the caller meant.")

    cmd = [str(exe), SUBCOMMAND]
    if device == "cuda":
        cmd.append("--cuda")
    else:
        cmd.extend(cpu_flags)
    cmd += ["--input", str(input_path), "--output", str(output_h5)]
    if fpr is not None:
        cmd += ["--fpr", _fmt_num(fpr)]
    if learning_rate is not None:
        cmd += ["--learning-rate", _fmt_num(learning_rate)]
    if epochs is not None:
        cmd += ["--epochs", _fmt_num(epochs)]
    if expected_cells is not None:
        cmd += ["--expected-cells", _fmt_num(expected_cells)]
    if total_droplets_included is not None:
        cmd += ["--total-droplets-included", _fmt_num(total_droplets_included)]
    if low_count_threshold is not None:
        cmd += ["--low-count-threshold", _fmt_num(low_count_threshold)]
    if checkpoint_mins is not None:
        cmd += ["--checkpoint-mins", _fmt_num(checkpoint_mins)]
    if cpu_threads is not None:
        cmd += ["--cpu-threads", _fmt_num(cpu_threads)]
    if posterior_batch_size is not None:
        cmd += ["--posterior-batch-size", _fmt_num(posterior_batch_size)]
    for ft in _seq("exclude_feature_types", exclude_feature_types):
        cmd += ["--exclude-feature-types", ft]
    cmd += list(_seq("extra_args", extra_args))
    return cmd


def build_env(env_bin: str | Path | None,
              base_path: str | None = None,
              extra: dict | None = None,
              pathsep: str | None = None) -> dict:
    """The environment overlay for a CellBender run, with the tool's own bin FIRST on PATH.

    CellBender renders its HTML report by calling `jupyter nbconvert` by bare name and does not
    check the result, so a trimmed scheduler environment produces a run that succeeds and quietly
    has no report. Prepending - rather than appending - matters: a `jupyter` from some other
    environment on PATH is a different Python and fails in a way that reads as a notebook error.

    Returns only the variables to overlay. Both executors merge this over the inherited
    environment, so nothing here erases anything.

    `pathsep` defaults to this machine's separator, which is the right one for LocalExecutor and
    the WRONG one for an orchestrator on Windows submitting to a POSIX scheduler. Pass ":"
    explicitly in that case: a PATH joined with ";" is one unusable entry, and the symptom is
    the missing-report failure this argument exists to prevent.

    A `PATH` inside `extra` is MERGED, not written over the top. `extra` used to be applied last,
    so a caller passing `extra={"PATH": ...}` - the ordinary way to add a scratch directory -
    erased the `env_bin` prefix entirely. Nothing looked wrong afterwards: a PATH was still
    present and still plausible, and the only symptom was CellBender failing on a bare
    `jupyter nbconvert`, which is the exact failure this function exists to prevent. The caller's
    PATH is now the BASE and `env_bin` is prepended to it, so neither half is lost and the
    ordering that matters - the tool's own environment first - is the one that survives.
    """
    sep = os.pathsep if pathsep is None else pathsep
    env = {str(k): str(v) for k, v in (extra or {}).items()}
    if not is_missing(env_bin):
        # Deliberately NOT routed through pathlib. An orchestrator running on Windows turns a
        # POSIX cluster path into backslashes the moment it becomes a Path, and the resulting
        # PATH entry names a directory that exists on neither machine.
        prefix = str(env_bin)
        if "PATH" in env:
            current = env["PATH"]
        elif base_path is None:
            current = os.environ.get("PATH", "")
        else:
            current = base_path
        env["PATH"] = f"{prefix}{sep}{current}" if current else prefix
    return env


#: Basenames CellBender uses for the checkpoint it writes into the WORKING directory. A leftover
#: one makes the next run resume from it, which is a feature when the parameters are the same and
#: a silent substitution when they are not - so its presence is reported, never deleted.
CHECKPOINT_PREFIX = "ckpt"


def clear_products(products: dict, *, sample: str | None = None) -> dict:
    """Delete every declared product BEFORE the command runs, and say what was deleted.

    This is the whole of the freshness guarantee, and it is deliberately done by removal rather
    than by a timestamp. Checking after the fact that an output EXISTS cannot distinguish the file
    this run wrote from the one the last run left: CellBender's checkpoint/posterior failure exits
    0 having written nothing, so with leftovers in place every check in `verify_products()` passes
    and the previous run's matrix is reported under this run's parameters. With the path emptied
    first, existence afterwards is proof of authorship. mtimes are not used as the proof because
    they are stamped by the compute node's clock while this deletion happens on the
    orchestrator's, and a freshness gate that fires on clock skew is a gate someone switches off.

    All seven declared products go, not only the `.h5`. A stale `_report.html`, `_metrics.csv` or
    `_cell_barcodes.csv` beside a fresh matrix is the same defect one step downstream - the audit
    reads those files - and a report older than the object it describes is not out of date, it is
    wrong in the one way nobody checks.

    The run's own `.log` is included. It is CellBender's account of the run that produced the
    matrix, and keeping the previous one beside a new matrix would mean `parse_epochs` and
    `detect_checkpoint_trouble` describing a different run than `verify_products` does.

    A file that cannot be deleted stops the run here rather than after several GPU-hours.
    """
    who = f"{sample}: " if sample else ""
    removed, kept = [], []
    for role, path in products.items():
        p = Path(path)
        if not p.exists():
            kept.append(role)
            continue
        try:
            p.unlink()
        except OSError as exc:
            raise TaskFailure(
                f"{who}a product of a previous run is in the way and could not be removed: {p}\n"
                f"  ({exc})\n"
                f"  It is removed BEFORE the run, not after, because a file that survives the "
                f"command is\n  indistinguishable from one the command wrote - and CellBender's "
                f"checkpoint failure exits\n  0 having written nothing. Leaving it would mean "
                f"recording the previous run's numbers\n  under this run's parameters. Move it "
                f"aside, or write into a fresh output directory.") from None
        removed.append(str(p))
    return {"removed": removed, "already_absent": sorted(kept)}


def existing_checkpoints(workdir: str | Path) -> list:
    """Checkpoint files already in the working directory. Reported, never deleted.

    CellBender resumes from `ckpt*.tar.gz` if it finds one, which is correct when the run is a
    continuation and a silent parameter substitution when it is not: the resumed run trains on
    from the earlier learning rate and epoch budget while this run records the ones it was asked
    for. Deleting it would break the legitimate case, so the fact is put in the metrics instead.
    """
    d = Path(workdir)
    if not d.is_dir():
        return []
    return sorted(str(p) for p in d.iterdir()
                  if p.is_file() and p.name.startswith(CHECKPOINT_PREFIX))


def expected_products(output_h5: str | Path) -> dict:
    """The sibling files CellBender derives from the `--output` stem.

    Named rather than globbed, so that a file this pipeline did not ask for is not swept into a
    result set, and so the caller can see what is expected before anything runs. Which of these
    actually appear is a measurement, made in `verify_products()`.
    """
    out = Path(output_h5)
    stem = out.name[:-3] if out.name.endswith(".h5") else out.stem
    d = out.parent
    return {
        "output_h5": out,
        "filtered_h5": d / f"{stem}_filtered.h5",
        "cell_barcodes_csv": d / f"{stem}_cell_barcodes.csv",
        "report_html": d / f"{stem}_report.html",
        "metrics_csv": d / f"{stem}_metrics.csv",
        "pdf": d / f"{stem}.pdf",
        "log": d / f"{stem}.log",
    }


# ---------------------------------------------------------------------------- log parsing


def detect_checkpoint_trouble(log_text: str) -> list:
    """Log lines matching the checkpoint/posterior failure mode. The LIST is the finding."""
    return [ln.strip() for ln in (log_text or "").splitlines()
            if _CHECKPOINT_TROUBLE.search(ln)]


def detect_nbconvert_trouble(log_text: str) -> list:
    """Log lines suggesting the report step could not find `jupyter nbconvert` on PATH."""
    return [ln.strip() for ln in (log_text or "").splitlines()
            if _NBCONVERT_TROUBLE.search(ln)]


def parse_epochs(log_text: str) -> dict:
    """Training and test loss per epoch, read off the run's own stdout.

    Pure and separately testable: hand it captured text and it returns the curve. Epochs are
    reported as they were logged, in order, with no interpolation over a missing test epoch -
    CellBender evaluates the test set every few epochs and inventing the gaps would smooth over
    exactly the wobble a learning-rate assessment is looking for.
    """
    epochs, train, test, train_epochs, test_epochs = [], [], [], [], []
    for line in (log_text or "").splitlines():
        m = _EPOCH_LINE.search(line)
        if not m:
            continue
        epoch = int(m.group(1))
        rest = m.group(2)
        epochs.append(epoch)
        mt = _TRAIN_LOSS.search(rest)
        if mt:
            train.append(float(mt.group(1)))
            train_epochs.append(epoch)
        ms = _TEST_LOSS.search(rest)
        if ms:
            test.append(float(ms.group(1)))
            test_epochs.append(epoch)
    return {
        "epochs_logged": epochs,
        "train_loss": train,
        "train_epochs": train_epochs,
        "test_loss": test,
        "test_epochs": test_epochs,
    }


def parse_version_from_log(log_text: str) -> str | None:
    """The CellBender version as the run itself reported it, or None if it did not.

    Preferred over asking the executable afterwards: the banner comes from the process that
    produced the matrix, whereas `cellbender --version` describes whatever is on PATH now.
    """
    for line in (log_text or "").splitlines():
        m = _LOG_VERSION.search(line)
        if m:
            return m.group(1)
    return None


def _version_token(text: str | None) -> str | None:
    if not text:
        return None
    m = _VERSION_TOKEN.search(text)
    return m.group(0) if m else None


# ---------------------------------------------------------------------------- verification


def looks_like_hdf5(path: str | Path) -> bool:
    """Does this file begin with the HDF5 signature? Standard library only, reads 8 bytes."""
    p = Path(path)
    try:
        with p.open("rb") as fh:
            return fh.read(8) == HDF5_MAGIC
    except OSError:
        return False


def verify_products(output_h5: str | Path,
                    *,
                    sample: str | None = None,
                    log_text: str | None = None,
                    min_bytes: int = MIN_OUTPUT_BYTES,
                    check_structure: bool = True) -> dict:
    """Prove the run produced a usable matrix, or raise naming the mode that produces none.

    Called after the command returned zero, because a zero exit is precisely what the
    checkpoint/posterior failure mode gives. Four checks, cheapest first, and each failure names
    what to look at rather than reporting that something went wrong.
    """
    out = Path(output_h5)
    who = f"{sample}: " if sample else ""
    products = expected_products(out)
    present = {k: v for k, v in products.items() if v.exists()}
    ckpt_lines = detect_checkpoint_trouble(log_text or "")
    nb_lines = detect_nbconvert_trouble(log_text or "")

    def _silent_failure(detail: str) -> TaskFailure:
        msg = [
            f"{who}CellBender exited 0 but {detail}",
            f"  expected output : {out}",
            f"  files that ARE present: "
            f"{', '.join(sorted(p.name for p in present.values())) or 'none of the expected set'}",
            "  A zero exit is not success for this tool. CellBender 0.3.2 against some",
            "  pyro/torch pairs cannot write its checkpoint, and its posterior step refuses to",
            "  run without one - so it trains for the full epoch budget, logs nothing fatal,",
            "  exits 0 and writes no matrix.",
        ]
        if ckpt_lines:
            msg.append("  log lines consistent with that failure mode:")
            msg += [f"    {ln}" for ln in ckpt_lines[:10]]
        else:
            msg.append("  No checkpoint message was found in the captured log; that is not "
                       "evidence against it,")
            msg.append("  the phrase list is a heuristic. Check the run's own "
                       f"{products['log'].name} as well.")
        return TaskFailure("\n".join(msg))

    if not out.exists():
        raise _silent_failure("the output .h5 does not exist.")
    size = out.stat().st_size
    if not looks_like_hdf5(out):
        raise _silent_failure(
            f"the output .h5 ({size:,} bytes) does not begin with the HDF5 signature, so it is "
            f"truncated or is not an HDF5 file.")
    if size < min_bytes:
        raise _silent_failure(
            f"the output .h5 is {size:,} bytes, below the {min_bytes:,} byte floor for a real "
            f"denoised matrix. Raise `min_bytes` if this is a deliberately small fixture.")

    structure = {"structure_checked": False, "barcodes": None, "genes": None, "stored_values": None}
    if is_true(check_structure, "check_structure"):
        try:
            import h5py  # noqa: F401
        except ImportError:
            # h5py absent is a fact about this machine, not about the matrix. Recorded as
            # not-checked rather than allowed to read as a pass.
            structure["structure_note"] = "h5py is not importable here; structure NOT checked"
        else:
            shape, nnz = _h5_matrix_shape(out)
            structure.update({"structure_checked": True, "genes": int(shape[0]),
                              "barcodes": int(shape[1]), "stored_values": int(nnz)})
            if shape[1] <= 0 or nnz <= 0:
                raise _silent_failure(
                    f"the output .h5 contains a matrix of {shape[0]:,} genes x {shape[1]:,} "
                    f"barcodes with {nnz:,} stored values - it is an empty container.")

    return {
        "output_h5": out,
        "size_bytes": size,
        "present": present,
        "checkpoint_log_lines": ckpt_lines,
        "nbconvert_log_lines": nb_lines,
        **structure,
    }


# ---------------------------------------------------------------------------- running


def run_remove_background(
    sample: str,
    input_path: str | Path,
    output_h5: str | Path,
    *,
    exe: str | Path = "cellbender",
    env_bin: str | Path | None = None,
    device: str = "cuda",
    fpr: float | None = DEFAULT_FPR,
    learning_rate: float | None = None,
    epochs: int | None = None,
    expected_cells: int | None = None,
    total_droplets_included: int | None = None,
    low_count_threshold: int | None = None,
    checkpoint_mins: float | None = None,
    cpu_threads: int | None = None,
    posterior_batch_size: int | None = None,
    exclude_feature_types: Sequence[str] = (),
    cpu_flags: Sequence[str] = (),
    extra_args: Sequence[str] = (),
    env: dict | None = None,
    log: str | Path | None = None,
    cwd: str | Path | None = None,
    timeout_s: int | None = None,
    min_bytes: int = MIN_OUTPUT_BYTES,
    check_structure: bool = True,
    executor=None,
) -> dict:
    """Denoise one library and return only what was verified afterwards.

    The working directory defaults to the output's own parent because CellBender writes its
    checkpoint into the CURRENT directory under a fixed name. Two samples sharing a working
    directory therefore overwrite each other's checkpoints, and the symptom is a resumed run that
    silently continues from the wrong library - so one output directory per sample is a
    requirement, not a tidiness preference.

    EVERY DECLARED PRODUCT IS DELETED BEFORE THE COMMAND IS LAUNCHED (`clear_products`), so a
    product found afterwards can only have been written by this invocation. Checking that an
    output exists after the command returns proves a file is there and nothing more: the
    checkpoint/posterior failure this module was written for exits 0 having written nothing, and
    with a previous run's `.h5` in place every check below passes and the old numbers are recorded
    under the new parameters. The deletion is what makes the checks mean something; the recorded
    mtimes are corroboration, not the proof, because they come from a different clock.

    A checkpoint left in the working directory is REPORTED and not deleted - CellBender resumes
    from it, which is right for a continuation and wrong for a re-run under new parameters, and
    only the caller knows which this is. `metrics["checkpoints_present_before_run"]` is the
    warning; an empty list is the ordinary case.

    Returns `{"outputs": [...], "metrics": {...}, "versions": {...}}`. `outputs` lists only files
    that were checked to exist after the command returned; the removal metrics are NOT computed
    here, because they need the raw matrix as well - see `parse_metrics()`.
    """
    if executor is None:
        raise TaskFailure(
            f"{sample}: run_remove_background needs an executor; there is no implicit local "
            f"fallback, because where a step ran is part of what a run has to record.")
    src = Path(input_path)
    if not src.exists():
        raise TaskFailure(
            f"{sample}: the raw input for CellBender does not exist: {src}\n"
            f"  Nothing is substituted for a missing input. Check step 0's output for this "
            f"sample before re-running.")
    out = Path(output_h5)
    if out.suffix != ".h5":
        raise TaskFailure(
            f"{sample}: --output must end in .h5, got {out.name!r}. CellBender derives the "
            f"names of its report, cell-barcode list and filtered matrix from this stem, and "
            f"`expected_products()` has to be able to predict them.")
    out.parent.mkdir(parents=True, exist_ok=True)
    workdir = Path(cwd) if cwd is not None else out.parent
    log_path = Path(log) if log is not None else out.parent / f"{out.stem}.scqc.stdout.log"
    check_structure = is_true(check_structure, "check_structure")

    cmd = build_command(
        exe, src, out,
        device=device, fpr=fpr, learning_rate=learning_rate, epochs=epochs,
        expected_cells=expected_cells, total_droplets_included=total_droplets_included,
        low_count_threshold=low_count_threshold, checkpoint_mins=checkpoint_mins,
        cpu_threads=cpu_threads, posterior_batch_size=posterior_batch_size,
        exclude_feature_types=exclude_feature_types, cpu_flags=cpu_flags,
        extra_args=extra_args,
    )
    overlay = build_env(env_bin, extra=env)

    # A previous run's output must never be accepted as this one's. The declared products are
    # removed BEFORE the command starts, so anything present when it returns was written by this
    # invocation - which no check made afterwards could establish on its own.
    products = expected_products(out)
    cleared = clear_products(products, sample=sample)
    checkpoints = existing_checkpoints(workdir)
    started = time.time()

    stdout = executor.shell(cmd, log_path, overlay or None, workdir, timeout_s)

    # The run's own .log is the fuller record; the captured stream may be only a tail of it.
    log_text = stdout or ""
    own_log = products["log"]
    if own_log.exists():
        try:
            log_text += "\n" + own_log.read_text(encoding="utf-8", errors="replace")
        except OSError:
            pass

    checked = verify_products(out, sample=sample, log_text=log_text,
                              min_bytes=min_bytes, check_structure=check_structure)

    outputs = [out]
    for role in ("filtered_h5", "cell_barcodes_csv", "report_html", "metrics_csv", "pdf", "log"):
        p = products[role]
        if p.exists():
            outputs.append(p)
    if log_path.exists():
        outputs.append(log_path)
    # Corroboration, recorded rather than enforced: each of these was deleted before the command
    # ran, so its presence is already proof of authorship. A negative or implausible age means the
    # two hosts disagree about the time, which is worth seeing and is not a reason to fail a run
    # that demonstrably produced output.
    output_ages = {}
    for p in outputs:
        try:
            output_ages[Path(p).name] = round(started - Path(p).stat().st_mtime, 3)
        except OSError:
            output_ages[Path(p).name] = None

    curve = parse_epochs(log_text)
    logged = parse_version_from_log(log_text)
    asked = tool_version(exe, ("--version",))
    a, b = _version_token(logged), _version_token(asked)
    if a and b and a != b:
        raise TaskFailure(
            f"{sample}: the CellBender that ran reported version {a} in its log while "
            f"`{exe} --version` reports {b}. Those are two different installations and the "
            f"provenance record cannot name one of them. Point `exe` at the environment that "
            f"actually ran, or re-run under a single installation.")
    # ONE FORMAT, ALWAYS. `versions["cellbender"]` is a `major.minor[.patch]` token or
    # NOT_INVOKED - never a raw banner line. It used to be whichever the capturing route happened
    # to produce: the log route yields a normalised token and `--version` yields a whole first
    # line, so a cohort's provenance column held "0.3.2" for some samples and
    # "cellbender, version 0.3.2" for others, and the two do not compare equal. Whatever each
    # route actually said is kept verbatim in the metrics, so normalising discards nothing.
    token = a or b
    if token is not None:
        version = token
        version_format = "normalised token major.minor[.patch]"
    else:
        version = NOT_INVOKED
        version_format = ("no version token could be extracted from either route; the raw "
                          "strings are in version_observed")

    metrics = {
        "sample": sample,
        "command": cmd,
        "device": device,
        "workdir": str(workdir),
        "stale_products_removed": cleared["removed"],
        "products_already_absent": cleared["already_absent"],
        "checkpoints_present_before_run": checkpoints,
        "output_age_s_at_start": output_ages,
        "started_epoch_s": started,
        "freshness_proof": ("every declared product was deleted before the command ran, so each "
                            "one listed in outputs was written by this invocation; mtimes are "
                            "corroboration from the compute node's clock, not the proof"),
        "version_format": version_format,
        "version_observed": {"run log": logged, "--version": asked},
        "learning_rate": learning_rate,
        "fpr": fpr,
        "output_bytes": checked["size_bytes"],
        "genes_in_output": checked.get("genes"),
        "barcodes_in_output": checked.get("barcodes"),
        "stored_values_in_output": checked.get("stored_values"),
        "structure_checked": checked["structure_checked"],
        "report_html_written": products["report_html"].exists(),
        "filtered_h5_written": products["filtered_h5"].exists(),
        "cell_barcodes_csv_written": products["cell_barcodes_csv"].exists(),
        "metrics_csv_written": products["metrics_csv"].exists(),
        "epochs_logged": len(curve["epochs_logged"]),
        "train_loss_points": len(curve["train_loss"]),
        "test_loss_points": len(curve["test_loss"]),
        "checkpoint_log_lines": checked["checkpoint_log_lines"],
        "nbconvert_log_lines": checked["nbconvert_log_lines"],
        "version_route": "run log" if logged else ("--version" if b else "not observed"),
        "log": str(log_path),
    }
    if not metrics["report_html_written"]:
        # Not a failure: the report is a by-product and the matrix is the deliverable. It is
        # reported so that a missing report is a fact in the manifest rather than a surprise.
        metrics["report_note"] = (
            "no report HTML was written. CellBender renders it by calling `jupyter nbconvert` "
            "by bare name and does not check the result; pass env_bin so the environment's own "
            "bin is first on PATH.")
    structure_note = checked.get("structure_note")
    if not is_missing(structure_note):
        metrics["structure_note"] = structure_note

    return {"outputs": outputs, "metrics": metrics, "versions": {"cellbender": version}}


# ---------------------------------------------------------------------------- h5 reading


def _decode(arr) -> list:
    out = []
    for v in arr:
        if isinstance(v, bytes):
            out.append(v.decode("utf-8", "replace"))
        else:
            out.append(str(v))
    return out


def _matrix_group_name(f) -> str:
    """Locate the one group holding a CellRanger-format sparse matrix.

    Discovered rather than assumed: CellBender writes the v3 layout under `matrix` and can write
    the v2 layout under a genome name, and a wrong guess reads real numbers out of the wrong
    array. Anything other than exactly one candidate raises with the groups it did find.
    """
    import h5py

    found = []
    if _MATRIX_KEYS <= set(f.keys()):
        found.append("/")

    def visit(name, obj):
        if isinstance(obj, h5py.Group) and _MATRIX_KEYS <= set(obj.keys()):
            found.append(name)

    f.visititems(visit)
    if len(found) != 1:
        raise TaskFailure(
            f"expected exactly one CellRanger-format matrix group in {f.filename}, found "
            f"{len(found)}: {found or 'none'}. Groups present at the top level: "
            f"{sorted(f.keys())}. A matrix read from the wrong group is a number that looks "
            f"measured and is not.")
    return found[0]


def _collect_datasets(f) -> dict:
    """basename -> [full paths]. Used only for the latents, whose location varies by version."""
    import h5py

    found: dict = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            found.setdefault(name.rsplit("/", 1)[-1], []).append(name)

    f.visititems(visit)
    return found


def _pick_dataset(found: dict, candidates: Sequence[str], what: str) -> str | None:
    """First candidate NAME that occurs, refusing only where one name occurs twice.

    The candidate tuples are an ordered preference, so a file carrying both an old and a new
    spelling of the same quantity resolves to the preferred one rather than being rejected. What
    is not resolvable is the same basename appearing at two paths: there is no preference to
    apply, and picking one silently is how a diagnostic ends up describing the wrong array.
    """
    for c in candidates:
        hits = found.get(c, [])
        if not hits:
            continue
        if len(hits) > 1:
            raise TaskFailure(
                f"{what}: the dataset name {c!r} occurs at {len(hits)} paths ({hits}). Which "
                f"one holds the value is ambiguous.")
        return hits[0]
    return None


def _h5_matrix_shape(path: Path) -> tuple:
    import h5py

    with h5py.File(str(path), "r") as f:
        g = f[_matrix_group_name(f)]
        shape = tuple(int(x) for x in g["shape"][:])
        nnz = int(g["data"].shape[0])
    return shape, nnz


def read_h5_header(path: str | Path) -> dict:
    """Barcodes and gene identifiers of a CellRanger-format .h5, without the counts."""
    import h5py

    p = Path(path)
    with h5py.File(str(p), "r") as f:
        g = f[_matrix_group_name(f)]
        shape = tuple(int(x) for x in g["shape"][:])
        barcodes = _decode(g["barcodes"][:])
        if "features" in g and hasattr(g["features"], "keys"):
            feats = g["features"]
            names = _decode(feats["name"][:]) if "name" in feats else None
            ids = _decode(feats["id"][:]) if "id" in feats else None
        else:
            names = _decode(g["gene_names"][:]) if "gene_names" in g else None
            ids = _decode(g["genes"][:]) if "genes" in g else None
    if names is None and ids is None:
        raise TaskFailure(
            f"{p} carries no gene names and no gene ids, so its rows cannot be identified. "
            f"Per-gene removal is reported by SYMBOL and an unlabelled matrix cannot supply one.")
    if names is None:
        names = list(ids)
    if ids is None:
        ids = list(names)
    if not (len(names) == len(ids) == shape[0]):
        raise TaskFailure(
            f"{p}: {shape[0]} matrix rows but {len(ids)} gene ids and {len(names)} gene names. "
            f"The file is internally inconsistent; nothing downstream can align to it.")
    if len(barcodes) != shape[1]:
        raise TaskFailure(
            f"{p}: {shape[1]} matrix columns but {len(barcodes)} barcodes.")
    return {"path": str(p), "shape": shape, "barcodes": barcodes,
            "gene_ids": ids, "gene_names": names}


def contiguous_runs(indices: Sequence[int]) -> list:
    """Maximal [start, stop) runs of consecutive integers in a sorted sequence.

    Pure, stdlib, and separately testable. It exists so that per-gene sums over a few thousand
    selected droplets can be read from a raw matrix of a few hundred thousand without loading the
    whole thing: consecutive columns of a CSC matrix occupy one contiguous slice of `data`.
    """
    idx = [int(i) for i in indices]
    if not idx:
        return []
    for a, b in zip(idx, idx[1:]):
        if b <= a:
            raise TaskFailure(
                f"contiguous_runs needs a strictly increasing sequence; got {a} before {b}.")
    runs = []
    start = prev = idx[0]
    for c in idx[1:]:
        if c == prev + 1:
            prev = c
            continue
        runs.append((start, prev + 1))
        start = prev = c
    runs.append((start, prev + 1))
    return runs


def _gene_stats_h5(path: Path, keep_cols: Sequence[int]) -> dict:
    """Per-gene summed counts and detection counts over selected columns of a CSC .h5."""
    import h5py
    import numpy as np

    with h5py.File(str(path), "r") as f:
        g = f[_matrix_group_name(f)]
        shape = tuple(int(x) for x in g["shape"][:])
        n_genes, n_cols = shape[0], shape[1]
        indptr = g["indptr"][:]
        if len(indptr) != n_cols + 1:
            raise TaskFailure(
                f"{path}: indptr has {len(indptr)} entries for {n_cols} columns; the matrix is "
                f"not the column-major CSC layout this reader assumes.")
        sums = np.zeros(n_genes, dtype=np.float64)
        det = np.zeros(n_genes, dtype=np.int64)
        total = 0.0
        data_ds, idx_ds = g["data"], g["indices"]
        for start, stop in contiguous_runs(keep_cols):
            lo, hi = int(indptr[start]), int(indptr[stop])
            if hi <= lo:
                continue
            d = np.asarray(data_ds[lo:hi], dtype=np.float64)
            i = np.asarray(idx_ds[lo:hi], dtype=np.int64)
            if i.size and int(i.max()) >= n_genes:
                raise TaskFailure(
                    f"{path}: a stored row index {int(i.max())} is outside the declared "
                    f"{n_genes} genes. The file is corrupt or is not gene-major.")
            sums += np.bincount(i, weights=d, minlength=n_genes)
            nz = d > 0
            det += np.bincount(i[nz], minlength=n_genes)
            total += float(d.sum())
    return {"sums": sums, "detected": det, "total": total, "n_cols_used": len(keep_cols)}


def _load_dense_source(path: Path) -> dict:
    """Read a MatrixMarket directory or an .h5ad into gene-major arrays.

    Whole-matrix reads, unlike the .h5 route: neither format allows a column slice to be located
    without parsing everything. Stated rather than hidden, because on an unfiltered droplet
    matrix this is the memory-expensive path.
    """
    import numpy as np
    import scipy.sparse as sp

    p = Path(path)
    if p.is_dir():
        import scipy.io as sio

        def _find(*names):
            for n in names:
                for cand in (p / n, p / (n + ".gz")):
                    if cand.exists():
                        return cand
            return None

        mtx = _find("matrix.mtx")
        feats = _find("features.tsv", "genes.tsv")
        bcs = _find("barcodes.tsv")
        if mtx is None or feats is None or bcs is None:
            raise TaskFailure(
                f"{p} is not a MatrixMarket triple: matrix.mtx={mtx}, features/genes.tsv="
                f"{feats}, barcodes.tsv={bcs}. All three are required and none is optional.")
        m = sp.csc_matrix(sio.mmread(str(mtx)))
        ids, names = [], []
        for row in _read_tsv(feats):
            ids.append(row[0])
            names.append(row[1] if len(row) > 1 else row[0])
        barcodes = [r[0] for r in _read_tsv(bcs)]
    elif p.suffix == ".h5ad" or is_anndata_h5(p):
        import anndata

        ad = anndata.read_h5ad(str(p))
        m = sp.csc_matrix(ad.X).T.tocsc()
        barcodes = [str(b) for b in ad.obs_names]
        names = [str(v) for v in ad.var_names]
        ids = [str(v) for v in (ad.var["gene_ids"] if "gene_ids" in ad.var else ad.var_names)]
    else:
        raise TaskFailure(
            f"{p}: unsupported raw matrix format. This adapter reads a CellRanger-format .h5, a "
            f"MatrixMarket directory (matrix.mtx + features/genes.tsv + barcodes.tsv) or an "
            f".h5ad.")
    if m.shape != (len(names), len(barcodes)):
        raise TaskFailure(
            f"{p}: matrix is {m.shape} but there are {len(names)} genes and {len(barcodes)} "
            f"barcodes. Expected genes as ROWS and barcodes as COLUMNS.")
    m = m.tocsc()
    return {"shape": m.shape, "barcodes": barcodes, "gene_ids": ids, "gene_names": names,
            "indptr": np.asarray(m.indptr), "indices": np.asarray(m.indices),
            "data": np.asarray(m.data)}


def _read_tsv(path: Path) -> list:
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(str(path), "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line:
                rows.append(line.split("\t"))
    return rows


def _gene_stats_dense(src: dict, keep_cols: Sequence[int]) -> dict:
    import numpy as np

    n_genes = src["shape"][0]
    indptr, indices, data = src["indptr"], src["indices"], src["data"]
    sums = np.zeros(n_genes, dtype=np.float64)
    det = np.zeros(n_genes, dtype=np.int64)
    total = 0.0
    for start, stop in contiguous_runs(keep_cols):
        lo, hi = int(indptr[start]), int(indptr[stop])
        if hi <= lo:
            continue
        d = np.asarray(data[lo:hi], dtype=np.float64)
        i = np.asarray(indices[lo:hi], dtype=np.int64)
        sums += np.bincount(i, weights=d, minlength=n_genes)
        det += np.bincount(i[d > 0], minlength=n_genes)
        total += float(d.sum())
    return {"sums": sums, "detected": det, "total": total, "n_cols_used": len(keep_cols)}


def is_anndata_h5(path: str | Path) -> bool:
    """Is this HDF5 file an AnnData object, whatever it is called?

    THE EXTENSION IS NOT THE FORMAT. A pipeline that accepts a denoised object and stores it under
    its own name - `<sample>_ambient.h5` - has an .h5ad wearing a .h5 suffix, and every reader that
    dispatches on the suffix then opens it as a CellRanger matrix and finds nothing it recognises.
    That is exactly the failure this function exists to stop: the error it produced named a missing
    matrix group, which reads as a corrupt file rather than as a reader sent down the wrong branch.

    AnnData writes `X`, `obs` and `var` at the top level; a CellRanger matrix writes a single group
    holding `data`, `indices`, `indptr` and `barcodes`. Distinguishing them is one read of the top
    level, and no heuristic beyond that is needed.
    """
    p = Path(path)
    if p.is_dir() or not p.exists():
        return False
    try:
        import h5py
    except ImportError:
        return p.suffix == ".h5ad"
    try:
        with h5py.File(p, "r") as f:
            keys = set(f.keys())
    except (OSError, ValueError):
        return False
    return {"X", "obs", "var"} <= keys


def _header_of(path: Path) -> dict:
    if path.is_dir() or path.suffix == ".h5ad" or is_anndata_h5(path):
        src = _load_dense_source(path)
        return {"path": str(path), "shape": src["shape"], "barcodes": src["barcodes"],
                "gene_ids": src["gene_ids"], "gene_names": src["gene_names"], "_dense": src}
    return read_h5_header(path)


def _gene_stats(header: dict, keep_cols: Sequence[int]) -> dict:
    if "_dense" in header:
        return _gene_stats_dense(header["_dense"], keep_cols)
    return _gene_stats_h5(Path(header["path"]), keep_cols)


# ---------------------------------------------------------------------------- cell calls


def cell_barcodes_from_csv(path: str | Path) -> list:
    """CellBender's `_cell_barcodes.csv`: one barcode per line, header optional."""
    p = Path(path)
    rows = []
    with p.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.reader(fh):
            if row and row[0].strip():
                rows.append(row[0].strip())
    if rows and rows[0].lower() in ("barcode", "barcodes", "cell_barcode", "cell_barcodes"):
        rows = rows[1:]
    if not rows:
        raise TaskFailure(
            f"{p} lists no cell barcodes. An empty cell list is not a cell count of zero to be "
            f"carried forward; it means the posterior step produced nothing.")
    if len(set(rows)) != len(rows):
        dupes = sorted({b for b in rows if rows.count(b) > 1})[:10]
        raise TaskFailure(
            f"{p} contains duplicate barcodes ({len(rows) - len(set(rows))} duplicates, e.g. "
            f"{dupes}). Counting cells from it would over-count.")
    return rows


def _cells_from_probability(h5_path: Path, barcodes: Sequence[str], threshold: float):
    """Barcodes passing the posterior threshold, or None if the file carries no such latent.

    None and `[]` mean different things and the distinction is the point: None is "this route is
    not available in this file", `[]` is "this route ran and found no cell". Returning `[]` for
    both let the caller read a route that disagreed completely with the CSV as a route that was
    not there, and skip the cross-check in the one case where it mattered most.
    """
    import h5py
    import numpy as np

    with h5py.File(str(h5_path), "r") as f:
        found = _collect_datasets(f)
        prob_path = _pick_dataset(found, _CELL_PROB_NAMES, f"{h5_path}: cell probability")
        if prob_path is None:
            return None
        prob = np.asarray(f[prob_path][:], dtype=np.float64)
        idx_path = _pick_dataset(found, _LATENT_INDEX_NAMES, f"{h5_path}: latent barcode index")
        idx = np.asarray(f[idx_path][:], dtype=np.int64) if idx_path else None
    if idx is None:
        if len(prob) != len(barcodes):
            raise TaskFailure(
                f"{h5_path}: {len(prob)} cell probabilities for {len(barcodes)} barcodes and no "
                f"index array ({list(_LATENT_INDEX_NAMES)}) saying which barcodes they describe. "
                f"Aligning them by position would be a guess.")
        sel = np.nonzero(prob >= threshold)[0]
    else:
        if len(prob) != len(idx):
            raise TaskFailure(
                f"{h5_path}: {len(prob)} cell probabilities against {len(idx)} latent barcode "
                f"indices.")
        sel = idx[prob >= threshold]
    return [barcodes[int(i)] for i in sel]


def resolve_cell_barcodes(h5_path: str | Path,
                          barcodes: Sequence[str],
                          *,
                          cell_barcodes_csv: str | Path | None = None,
                          threshold: float = CELL_PROBABILITY_THRESHOLD) -> dict:
    """Which barcodes are cells, by two independent routes, compared against each other.

    The cell count is consequential - it is the denominator of everything step 2 gates on - so it
    is derived twice where both routes are available and a disagreement stops the run rather than
    being resolved by preferring one. Where neither route is available this raises: there is no
    fallback that assumes every barcode in the file is a cell.

    AN EMPTY RESULT IS A RESULT, not a missing route. The two were conflated by `if got:`, which
    treated a probability route returning no barcodes as a route that was not there - and so
    disabled the cross-check in precisely the case where the routes disagree most: a CSV listing
    thousands of cells beside a posterior that clears none of them was reported as the CSV,
    confirmed by nothing. The route now returns None when the file has no cell-probability latent
    and a possibly-empty list when it ran, the cross-check reads the two states apart, and an
    empty cell call is refused outright - the posterior producing nothing is not a library of
    zero cells to be carried forward.
    """
    h5_path = Path(h5_path)
    from_csv = None
    csv_path = Path(cell_barcodes_csv) if cell_barcodes_csv is not None \
        else expected_products(h5_path)["cell_barcodes_csv"]
    csv_route = "absent"
    naming_note = None
    if csv_path.exists():
        from_csv = cell_barcodes_from_csv(csv_path)
        # RECONCILED BEFORE THE TWO ROUTES ARE COMPARED, not after.
        #
        # The probability route derives its barcodes from the OBJECT, so it is already in the
        # object's naming; the CSV is in the denoiser's. Comparing them unreconciled makes the
        # cross-check fire on two descriptions of one cell call - it reported 38,285 against
        # 38,285 with every barcode differing, which is a total naming mismatch presenting as a
        # total disagreement about which cells exist. Reconciling afterwards, as this first did,
        # is too late for the check that runs in between.
        from_csv, naming_note = reconcile_barcode_naming(barcodes, from_csv)
        csv_route = f"read {len(from_csv):,} barcodes"
        if naming_note:
            csv_route += f"; {naming_note}"

    from_prob = None
    if h5_path.suffix == ".h5" and not h5_path.is_dir():
        try:
            from_prob = _cells_from_probability(h5_path, barcodes, threshold)
        except ImportError:
            from_prob = None
            prob_route = "unavailable: h5py is not importable here, so the route did not run"
        else:
            prob_route = ("unavailable: the file carries no cell-probability latent "
                          f"({list(_CELL_PROB_NAMES)})" if from_prob is None
                          else f"ran: {len(from_prob):,} barcodes at probability >= {threshold}")
    else:
        prob_route = "not attempted: the denoised object is not a single .h5 file"

    routes = {"cell_barcodes_csv": csv_route, "cell_probability": prob_route}

    if from_csv is not None and from_prob is not None:
        if len(from_csv) != len(from_prob) or set(from_csv) != set(from_prob):
            raise TaskFailure(
                f"{h5_path}: the two routes to the cell call disagree - "
                f"{csv_path.name} lists {len(from_csv):,} barcodes while thresholding the "
                f"posterior cell probability at {threshold} gives {len(from_prob):,} "
                f"({len(set(from_csv) ^ set(from_prob)):,} barcodes differ). One of them does "
                f"not describe this run; do not pick one.")
        call = {"barcodes": from_csv, "cross_checked": True, "routes": routes,
                "source": f"{csv_path.name} (confirmed against cell probability >= {threshold})"}
    elif from_csv is not None:
        call = {"barcodes": from_csv, "cross_checked": False, "routes": routes,
                "source": csv_path.name}
    elif from_prob is not None:
        call = {"barcodes": from_prob, "cross_checked": False, "routes": routes,
                "source": f"cell probability >= {threshold}"}
    else:
        raise TaskFailure(
            f"{h5_path}: the cell call could not be established. Neither {csv_path.name} nor a "
            f"cell-probability latent ({list(_CELL_PROB_NAMES)}) was found. Pass "
            f"`cell_barcodes_csv` explicitly. Treating every barcode in the file as a cell is not "
            f"offered, because on the full output that is every droplet in the library.")

    if not call["barcodes"]:
        raise TaskFailure(
            f"{h5_path}: the cell call is EMPTY by {call['source']}. That is not a library of "
            f"zero cells to be carried forward as a denominator; it means the posterior step "
            f"produced nothing, and everything computed over the 'cells' scope would divide by "
            f"it. Routes: {routes}")

    # The CSV route was reconciled where it was read, above; the probability route comes from the
    # object and needs none. This records what was done so the result says how it got there.
    if naming_note:
        call["naming"] = naming_note
    return call


def reconcile_barcode_naming(object_barcodes: Sequence[str],
                             call_barcodes: Sequence[str]) -> tuple:
    """Make a cell call written in one naming usable against an object written in another.

    Returns `(barcodes, note)`. The note is None when nothing was changed, and the barcodes come
    back unaltered.

    THE SAME LIBRARY, DESCRIBED TWICE. A denoiser writes its cell call as bare barcodes; a
    pipeline that will concatenate libraries prefixes each barcode with its sample so the combined
    object has unique names. Both are correct, and they intersect in NOTHING - so every quantity
    summed over "the cells" silently becomes a quantity summed over an empty set. It is not a
    corrupt file, a version mismatch or a wrong path, and every usual explanation fails to apply,
    which is what makes it expensive to find.

    THE TRANSFORM IS VERIFIED, NOT GUESSED. A prefix is adopted only when every transformed
    barcode is then present in the object. A partial match is refused: a rule that repairs most of
    a cell call and quietly drops the rest yields a denominator wrong by an amount nobody can see.
    Where the two already agree, or where no single transform explains the difference, this
    changes nothing and the caller's own refusal stands.
    """
    obj = {str(b) for b in object_barcodes}
    call = [str(b) for b in call_barcodes]
    if not call or not obj or any(b in obj for b in call):
        return call, None

    # One candidate, derived from one barcode and then tested against every one of them. Deriving
    # it from a single pair is safe BECAUSE it is verified exhaustively afterwards.
    first = call[0]
    for cand in obj:
        if cand.endswith(first) and len(cand) > len(first):
            prefix = cand[:len(cand) - len(first)]
            moved = [prefix + b for b in call]
            if all(b in obj for b in moved):
                return moved, (f"the cell call shared no barcode with the object, and all "
                               f"{len(moved):,} are present under the prefix {prefix!r}; adopted "
                               f"after checking every one")
            break

    # The other direction: a prefixed call against an object that does not use the prefix.
    for cut in range(1, len(first)):
        if first[cut - 1] != "_":
            continue
        head = first[:cut]
        stripped = [b[cut:] if b.startswith(head) else b for b in call]
        if all(b in obj for b in stripped):
            return stripped, (f"the cell call shared no barcode with the object, and all "
                              f"{len(stripped):,} are present with the prefix {head!r} removed; "
                              f"adopted after checking every one")
    return call, None


# ---------------------------------------------------------------------------- metrics


def parse_metrics(h5_path: str | Path,
                  raw_input_path: str | Path,
                  *,
                  cell_barcodes_csv: str | Path | None = None,
                  scope: str = "cells",
                  cell_probability_threshold: float = CELL_PROBABILITY_THRESHOLD) -> dict:
    """What the denoiser removed, measured against the matrix it was given.

    `raw_input_path` is required and has no default. The denoised object does not carry the
    counts it started from, so a removal fraction computed from it alone would have to come from
    somewhere else - and the only somewhere else is an assumption.

    `scope` is the barcode set both matrices are summed over:

        "cells"   the called cells (default). This is the population every downstream step uses,
                  and the one the cohort ranges in modules/01_ambient/ambient.py describe.
        "shared"  every barcode present in both matrices, empties included. The fraction removed
                  over this scope is dominated by the empty droplets, which are reduced to almost
                  nothing by design; it is a different quantity and must not be compared with a
                  cells-scoped range.

    Genes are aligned by identifier, never by position: a run using `--exclude-feature-types`
    writes fewer rows than it read, and aligning those by row number silently relabels the whole
    matrix. Genes present in the raw matrix and absent from the output are returned as a LIST
    (docs/PRINCIPLES.md, question 1), not as a count.
    """
    if scope not in ("cells", "shared"):
        raise TaskFailure(f"scope must be 'cells' or 'shared', got {scope!r}.")
    den_path, raw_path = Path(h5_path), Path(raw_input_path)
    for p, what in ((den_path, "denoised output"), (raw_path, "raw input")):
        if not p.exists():
            raise TaskFailure(f"the {what} does not exist: {p}")

    den = _header_of(den_path)
    raw = _header_of(raw_path)

    # ---- barcode scope
    den_pos = {b: i for i, b in enumerate(den["barcodes"])}
    raw_pos = {b: i for i, b in enumerate(raw["barcodes"])}
    if len(den_pos) != len(den["barcodes"]) or len(raw_pos) != len(raw["barcodes"]):
        raise TaskFailure(
            f"duplicate barcodes: {len(den['barcodes']) - len(den_pos)} in {den_path.name}, "
            f"{len(raw['barcodes']) - len(raw_pos)} in {raw_path.name}. Summing over a scope "
            f"defined by barcode string would double-count them.")

    if scope == "cells":
        call = resolve_cell_barcodes(den_path, den["barcodes"],
                                     cell_barcodes_csv=cell_barcodes_csv,
                                     threshold=cell_probability_threshold)
        wanted, cell_source = call["barcodes"], call["source"]
    else:
        wanted = [b for b in den["barcodes"] if b in raw_pos]
        cell_source = "not used: scope='shared'"

    # TWO MATRICES, AND THEY NEED NOT BE NAMED THE SAME WAY. This is the one place in the pipeline
    # where two barcode conventions genuinely meet: the raw matrix as the aligner wrote it, and a
    # denoised object that a later step may have prefixed with its sample so a cohort can be
    # concatenated without collisions. The scope is therefore reconciled against each side
    # SEPARATELY, because neither naming is wrong and only one of them can be the scope's.
    #
    # What has to hold is that both sides select the SAME SET of droplets. It does not have to be
    # the same ORDER: the per-gene statistics below are sums and detection counts over the chosen
    # columns, and both column lists are sorted independently anyway. The reconciliation is a
    # total mapping - adopted only where every transformed barcode is present - so the two
    # selections are the same droplets or the refusal below fires.
    wanted_den, den_naming = reconcile_barcode_naming(den["barcodes"], wanted)
    wanted_raw, raw_naming = reconcile_barcode_naming(raw["barcodes"], wanted)
    missing_den = [b for b in wanted_den if b not in den_pos]
    missing_raw = [b for b in wanted_raw if b not in raw_pos]
    if missing_den or missing_raw:
        raise TaskFailure(
            f"{len(missing_den):,} of {len(wanted):,} scope barcodes are absent from "
            f"{den_path.name} and {len(missing_raw):,} from {raw_path.name} (e.g. "
            f"{(missing_den or missing_raw)[:3]}). Barcode strings must match exactly between "
            f"the raw matrix and the denoised one; a '-1' suffix added or stripped by an "
            f"intermediate step is the usual cause, and matching on a stripped form would pair "
            f"counts from different droplets.")
    if not wanted:
        raise TaskFailure(
            f"the {scope} scope is empty for {den_path.name}: there is nothing to measure "
            f"removal over.")

    if den_naming or raw_naming:
        cell_source += (f" [naming: denoised {den_naming or 'as given'}; "
                        f"raw {raw_naming or 'as given'}]")
    den_cols = sorted(den_pos[b] for b in wanted_den)
    raw_cols = sorted(raw_pos[b] for b in wanted_raw)

    # ---- gene alignment, by identifier
    den_key = den["gene_ids"] if len(set(den["gene_ids"])) == len(den["gene_ids"]) \
        else den["gene_names"]
    raw_key = raw["gene_ids"] if len(set(raw["gene_ids"])) == len(raw["gene_ids"]) \
        else raw["gene_names"]
    raw_index = {}
    for i, k in enumerate(raw_key):
        raw_index.setdefault(k, i)
    unmatched = [k for k in den_key if k not in raw_index]
    if unmatched:
        raise TaskFailure(
            f"{len(unmatched):,} genes in {den_path.name} are not in {raw_path.name} (e.g. "
            f"{unmatched[:5]}). The denoised matrix must be a subset of the matrix it was made "
            f"from; if it is not, these are not the same pair of files.")

    den_stats = _gene_stats(den, den_cols)
    raw_stats = _gene_stats(raw, raw_cols)

    total_raw = raw_stats["total"]
    total_den = den_stats["total"]
    if total_raw <= 0:
        raise TaskFailure(
            f"the raw matrix holds no counts over the {scope} scope ({len(wanted):,} barcodes). "
            f"A removal fraction against a zero denominator is undefined, not 1.0.")
    fraction_removed_overall = 1.0 - (total_den / total_raw)
    if fraction_removed_overall < 0:
        raise TaskFailure(
            f"the denoised matrix holds MORE counts than the raw one over the {scope} scope "
            f"({total_den:,.0f} against {total_raw:,.0f}). A denoiser does not add counts, so "
            f"the two files are misaligned - most often the wrong raw matrix, or barcodes that "
            f"matched by string while describing different droplets.")

    n_scope = float(len(wanted))
    per_gene = []
    genes_fully_removed = 0
    genes_no_raw_counts = 0
    genes_gaining_counts = []
    den_sums, den_det = den_stats["sums"], den_stats["detected"]
    raw_sums, raw_det = raw_stats["sums"], raw_stats["detected"]
    for j, key in enumerate(den_key):
        i = raw_index[key]
        r = float(raw_sums[i])
        d = float(den_sums[j])
        if r <= 0:
            # No counts to remove means no removal fraction. Not 0.0, not 1.0: it was never
            # measured, so the gene is excluded from the table and counted here instead.
            genes_no_raw_counts += 1
            continue
        frac = 1.0 - (d / r)
        if d <= 0:
            genes_fully_removed += 1
        if frac < 0:
            genes_gaining_counts.append(den["gene_names"][j])
        per_gene.append({
            "symbol": den["gene_names"][j],
            "gene_id": den["gene_ids"][j],
            "fraction_removed": frac,
            "raw_detection_frac": float(raw_det[i]) / n_scope,
            "denoised_detection_frac": float(den_det[j]) / n_scope,
            "raw_counts": r,
            "denoised_counts": d,
        })

    den_key_set = set(den_key)
    absent_from_output = [raw["gene_names"][i]
                          for k, i in raw_index.items() if k not in den_key_set]

    droplets_analyzed = _droplets_analyzed(den_path)

    return {
        "fraction_removed_overall": fraction_removed_overall,
        "genes_fully_removed": genes_fully_removed,
        "cells_called": len(wanted) if scope == "cells" else None,
        "cell_call_source": cell_source,
        "scope": scope,
        "scope_barcodes": len(wanted),
        "droplets_total_raw": len(raw["barcodes"]),
        "droplets_in_output": len(den["barcodes"]),
        "droplets_analyzed": droplets_analyzed,
        "total_counts_raw_in_scope": total_raw,
        "total_counts_denoised_in_scope": total_den,
        "genes_total": len(den_key),
        "genes_with_no_raw_counts_in_scope": genes_no_raw_counts,
        "genes_gaining_counts": genes_gaining_counts,
        "genes_absent_from_output": absent_from_output,
        "per_gene": per_gene,
        "denoised_path": str(den_path),
        "raw_path": str(raw_path),
    }


def _droplets_analyzed(h5_path: Path) -> int | None:
    """How many droplets CellBender actually fitted, if the file records it."""
    if h5_path.is_dir() or h5_path.suffix != ".h5":
        return None
    try:
        import h5py
    except ImportError:
        return None
    with h5py.File(str(h5_path), "r") as f:
        found = _collect_datasets(f)
        path = _pick_dataset(found, _LATENT_INDEX_NAMES, f"{h5_path}: latent barcode index")
        if path is None:
            return None
        return int(f[path].shape[0])


# ---------------------------------------------------------------------------- tables


#: Values the audit reads as numbers. Each must be PRESENT and must carry a value; a key holding
#: None, NaN or pandas.NA is not a measurement that came out blank, it is no measurement at all.
_SUMMARY_REQUIRED = ("fraction_removed_overall", "genes_fully_removed")
_PER_GENE_REQUIRED = ("symbol", "fraction_removed", "raw_detection_frac",
                      "denoised_detection_frac")


def summary_row(sample: str, metrics: dict) -> dict:
    """One row of the per-sample table `audit_ambient.audit()` reads.

    `sample`, `fraction_removed_overall` and `genes_fully_removed` are the three that module
    names. The rest travel with them so the table answers more than the gate asks.

    THE GUARD IS ON THE VALUE, NOT ON THE KEY. `if required not in metrics` is satisfied by
    `{"fraction_removed_overall": None}`, so a dict that carried the key and no measurement
    passed, the None was copied into the audit table, and every consumer downstream saw a column
    that exists and is blank rather than a run whose removal fraction was never obtained. Both
    required names are now checked with `is_missing`, which catches None, NaN, pandas.NA,
    pandas.NaT and numpy.ma.masked alike; the optional columns are normalised to None so that no
    NaN or NA sentinel is ever written into a CSV a gate reads.
    """
    for required in _SUMMARY_REQUIRED:
        if required not in metrics:
            raise TaskFailure(
                f"{sample}: metrics carry no {required!r}. The audit reads that column by name "
                f"and a row without it would be audited as though the value were absent from "
                f"the run rather than from this dict.")
        if is_missing(metrics[required]):
            raise TaskFailure(
                f"{sample}: metrics carry {required!r} as {metrics[required]!r}, which is not a "
                f"measurement. The key being present is not the same as the value being there: "
                f"copied into the audit table it reads as a run whose removal was measured and "
                f"came out blank, rather than as one where it was never obtained. Recompute it "
                f"with parse_metrics(denoised, raw), or leave this sample out of the table and "
                f"say why.")
    row = {"sample": sample}
    for col in SUMMARY_COLUMNS[1:]:
        value = metrics.get(col)
        row[col] = None if is_missing(value) else value
    return row


def per_gene_rows(sample: str, metrics: dict) -> list:
    """The per-gene table `audit_ambient.audit()` reads, one row per gene for one sample.

    Guarded on values for the same reason as `summary_row`: `parse_metrics()` never emits a gene
    with a missing fraction - a gene with no raw counts is excluded and counted separately - so a
    missing one here came from somewhere else, and writing it into the table the audit reads
    would let a NaN be compared with `>=` and pass.
    """
    if "per_gene" not in metrics:
        raise TaskFailure(f"{sample}: metrics carry no per-gene table.")
    rows = []
    for i, g in enumerate(metrics["per_gene"]):
        for required in _PER_GENE_REQUIRED:
            if required not in g or is_missing(g[required]):
                raise TaskFailure(
                    f"{sample}: per-gene entry {i} has {required!r} = "
                    f"{g.get(required)!r}. A gene whose removal was never measured is excluded "
                    f"from this table by parse_metrics() and counted as "
                    f"genes_with_no_raw_counts_in_scope; one that reaches here with a blank "
                    f"would be compared against the audit's thresholds and pass every one.")
        row = {"sample": sample}
        for col in PER_GENE_COLUMNS[1:]:
            value = g.get(col)
            row[col] = None if is_missing(value) else value
        rows.append(row)
    return rows


def write_csv(rows: Sequence[dict], path: str | Path, columns: Sequence[str]) -> Path:
    """Write a table with the standard library, so reading it back needs nothing installed."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(columns), extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    if not p.exists():
        raise TaskFailure(f"wrote {p} and it is not on disk afterwards.")
    return p


def write_summary_csv(rows: Sequence[dict], path: str | Path) -> Path:
    return write_csv(rows, path, SUMMARY_COLUMNS)


def write_per_gene_csv(rows: Sequence[dict], path: str | Path) -> Path:
    return write_csv(rows, path, PER_GENE_COLUMNS)


# ---------------------------------------------------------------------------- learning curve


def convergence_indicator(values: Sequence[float],
                          window: int = LEARNING_CURVE_WINDOW) -> float | None:
    """How much of the curve's total movement happened in its final `window` epochs, per cent.

    A scale-free summary, which it has to be: ELBO magnitudes differ by orders of magnitude
    between libraries, so any absolute quantity compares libraries by their depth rather than by
    their fit. Expressed as a percentage of the total movement from first to last recorded epoch,
    with both terms taken in absolute value so that the sign convention of the reported loss does
    not change the answer, and so that late motion in the WRONG direction still registers as
    late motion.

    `modules/01_ambient/lr_policy.py` uses this only as an outlier statistic across a cohort and
    asserts no direction of better from it. Neither does this function: a large value means the
    curve was still moving at the end relative to how far it moved overall, which is a reason to
    look, not a verdict.

    THE WINDOW IS NOT SHRUNK TO FIT A SHORT CURVE. It used to be - `w = min(window, len - 1)` -
    and the consequence was that every curve with at most `window + 1` points reported exactly
    100.0: `w` collapsed to `len - 1`, `vals[-1-w]` became `vals[0]`, and the tail was the whole
    curve by construction. A run that logged 8 epochs therefore reported perfect convergence,
    which is the most reassuring possible answer arrived at without measuring anything, and it
    survived into a cohort statistic where it read as a well-converged sample. Two further
    reasons to require the full window rather than shrink it: a shrunken window measures a
    different quantity per sample, and `lr_policy` compares samples against each other; and the
    saturating value is 100.0, so the artefact hides among plausible numbers instead of standing
    out. A curve with `window + 1` points or fewer now returns None - not derivable, which
    `lr_policy._mad_outliers` skips - rather than a number that is true by arithmetic.

    Returns None - never 0.0 - where the curve is too short to have a tail, holds a value that is
    unknown or non-finite, or did not move at all. Unknown is not a value (docs/PRINCIPLES.md
    section 4).
    """
    import math

    w = int(window)
    if w < 1:
        return None
    vals = []
    for v in values:
        # A curve carrying a missing value is not a shorter curve; it is one this statistic
        # cannot be computed over, because the gap is where the movement would have been.
        if is_missing(v):
            return None
        try:
            vals.append(float(v))
        except (TypeError, ValueError):
            return None
    if len(vals) < 3:
        return None
    if any(math.isnan(x) or math.isinf(x) for x in vals):
        return None
    if len(vals) <= w + 1:
        # The window would span the entire curve, making tail == total and the answer 100.0
        # regardless of the data. Not derivable is the honest report.
        return None
    total = abs(vals[-1] - vals[0])
    if total == 0:
        return None
    tail = abs(vals[-1] - vals[-1 - w])
    return 100.0 * tail / total


def _learning_curve_from_h5(path: Path) -> dict:
    import h5py
    import numpy as np

    with h5py.File(str(path), "r") as f:
        found = _collect_datasets(f)

        def take(names, what):
            p = _pick_dataset(found, names, f"{path}: {what}")
            return None if p is None else np.asarray(f[p][:], dtype=np.float64).tolist()

        train = take(_TRAIN_ELBO_NAMES, "training ELBO")
        test = take(_TEST_ELBO_NAMES, "test ELBO")
        train_ep = take(_TRAIN_EPOCH_NAMES, "training epochs")
        test_ep = take(_TEST_EPOCH_NAMES, "test epochs")
        observed = take(_CONVERGENCE_NAMES, "convergence indicator")
    return {"train": train or [], "test": test or [],
            "train_epochs": train_ep or [], "test_epochs": test_ep or [],
            "observed_indicator": (observed[0] if observed else None)}


def _looks_numeric(text) -> bool:
    """Does this cell hold a number? Used only to tell a header cell from a value cell."""
    try:
        float(str(text).strip().replace(",", "").rstrip("%"))
    except (TypeError, ValueError):
        return False
    return True


def _read_metrics_csv(path: Path) -> dict:
    """CellBender's own metrics table, as key -> raw string. Two shapes are accepted.

    A one-row wide table (header + one data row) and a two-column long table (name,value) both
    occur in the wild. Any other shape is left unread rather than interpreted: a value taken from
    a column that happened to be in the right place is indistinguishable from a measured one.

    THE TWO-COLUMN CASE IS THE ONE THAT WAS WRONG. The wide-form test required more than two
    columns, so `fraction_counts_removed,n_cells` over `0.09,4218` - a perfectly ordinary
    header-plus-one-row table with two columns - fell through to the long-form branch and was
    read as `{"fraction_counts_removed": "n_cells", "0.09": "4218"}`. The removal fraction then
    appeared to be the string `n_cells`, which `parse_learning_curve` reports as a field that
    exists and cannot be parsed - a confusing error about a real number that was sitting in the
    file all along.

    Two rows of two cells is genuinely ambiguous, so it is settled by the one asymmetry the two
    shapes have: the second cell of the first row is a FIELD NAME in the wide form and a VALUE in
    the long form. If it is a number, the row is a name/value pair and the table is long; if it
    is not, the row is a header and the table is wide. Where neither reading is supported - both
    cells non-numeric in a table that could be either - the wide reading is taken for the same
    reason `csv.DictReader` would take it, and the ambiguity is noted here rather than hidden:
    every field CellBender writes into this file is numeric, so a non-numeric second cell in row
    one is a header.
    """
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        rows = [r for r in csv.reader(fh) if r]
    if not rows:
        return {}
    width = len(rows[0])
    if len(rows) == 2 and all(len(r) == width for r in rows) and width >= 2:
        if width > 2 or not _looks_numeric(rows[0][1]):
            return {str(k).strip(): str(v).strip() for k, v in zip(rows[0], rows[1])}
    if all(len(r) >= 2 for r in rows):
        return {str(r[0]).strip(): str(r[1]).strip() for r in rows}
    return {}


def parse_learning_curve(h5_or_log: str | Path,
                         *,
                         window: int = LEARNING_CURVE_WINDOW,
                         metrics_csv: str | Path | None = None) -> dict:
    """Cohort-comparable diagnostics for `lr_policy.assess_cohort()`, from one run's artifacts.

    `assess_cohort` reads `{sample: {diagnostic: value}}` and consumes the two names in
    `lr_policy.DIAGNOSTICS`: `fraction_removed` and `convergence_indicator`. Those two keys are
    the contract; everything else returned here is context for a human and is ignored by the
    policy.

    `fraction_removed` IS ONLY PRESENT WHEN AN ARTIFACT SUPPLIED IT. The denoised object does not
    carry the counts it started from, so this function usually cannot produce it and the key is
    then ABSENT rather than None - which lets a caller merge this dict over the one from
    `parse_metrics()` without a blank overwriting a measurement. The intended assembly is:

        diag = {**parse_learning_curve(log),
                "fraction_removed": parse_metrics(...)["fraction_removed_overall"]}

    `convergence_indicator` prefers a field of that exact name if the run recorded one, and
    otherwise is computed by `convergence_indicator()` from the training curve. Which route was
    used is returned as `convergence_indicator_source`, because a cohort assembled from a mixture
    of the two is not comparable with itself - the policy compares samples against each other.
    """
    p = Path(h5_or_log)
    if not p.exists():
        raise TaskFailure(f"no such artifact to read learning-curve diagnostics from: {p}")

    out: dict = {"source_path": str(p)}
    curve: list = []
    observed = None
    source = "not derivable"

    if looks_like_hdf5(p):
        got = _learning_curve_from_h5(p)
        curve = got["train"]
        observed = got["observed_indicator"]
        # The number of RECORDED training points, not an epoch number: whether CellBender
        # numbers its epochs from 0 or from 1 was not established, and adding one to a maximum
        # to convert between them would be a guess dressed as a count.
        out["epochs_trained"] = len(curve) or None
        out["last_train_epoch_recorded"] = (int(max(got["train_epochs"]))
                                            if got["train_epochs"] else None)
        out["train_elbo_final"] = curve[-1] if curve else None
        out["test_elbo_final"] = got["test"][-1] if got["test"] else None
        out["curve_route"] = "h5 datasets"
    else:
        text = p.read_text(encoding="utf-8", errors="replace")
        got = parse_epochs(text)
        curve = got["train_loss"]
        out["epochs_trained"] = (max(got["epochs_logged"]) if got["epochs_logged"] else None)
        out["train_elbo_final"] = curve[-1] if curve else None
        out["test_elbo_final"] = got["test_loss"][-1] if got["test_loss"] else None
        out["curve_route"] = "log text"
        trouble = detect_checkpoint_trouble(text)
        if trouble:
            out["checkpoint_log_lines"] = trouble
    out["curve_points"] = len(curve)

    if observed is not None:
        out["convergence_indicator"] = float(observed)
        source = "observed: a field named convergence_indicator in the artifact"
    else:
        ind = convergence_indicator(curve, window=window)
        if ind is not None:
            out["convergence_indicator"] = ind
            source = (f"computed: per cent of the total training-curve movement occurring in "
                      f"the final {int(window)} recorded epochs, over {len(curve)} points")
        else:
            # No curve, a curve no longer than the window, a flat one or a non-finite one. The
            # key is present as None so that a sample which was examined and yielded nothing is
            # distinguishable from one that was never examined; `lr_policy._mad_outliers` skips
            # None and never reads it as a value.
            out["convergence_indicator"] = None
            source = (f"not derivable: fewer than three recorded epochs, no more than "
                      f"{int(window) + 1} of them (the window would span the whole curve and the "
                      f"answer would be 100.0 by construction), a curve that did not move, or a "
                      f"loss that was non-finite or missing")
    out["convergence_indicator_source"] = source

    mcsv = Path(metrics_csv) if metrics_csv is not None else expected_products(
        p if p.suffix == ".h5" else p.with_suffix(".h5"))["metrics_csv"]
    if mcsv.exists():
        fields = _read_metrics_csv(mcsv)
        for name in _METRICS_FRACTION_REMOVED:
            if name in fields:
                try:
                    out["fraction_removed"] = float(fields[name])
                except ValueError:
                    raise TaskFailure(
                        f"{mcsv}: field {name!r} is {fields[name]!r}, which is not a number. "
                        f"It is not read as missing, because a field that exists and cannot be "
                        f"parsed is a different problem from one that is absent.") from None
                out["fraction_removed_source"] = f"{mcsv.name}:{name}"
                break
    out.setdefault("fraction_removed_source",
                   "absent: no artifact here carries the raw counts, so a removal fraction "
                   "must come from parse_metrics(denoised, raw)")
    return out
