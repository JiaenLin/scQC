"""Make the orchestrator's native stack safe to use from its worker threads.

WHY THIS EXISTS
---------------
`pipeline.run` fans tasks out across a `ThreadPoolExecutor`. Under `--executor pbs` those worker
threads do two things at once that do not mix: several of them read `.h5ad` objects through
`anndata`/`pandas`, while the others sit in `executor.shell` spawning `qsub`/`qstat` subprocesses
in a poll loop. On a multi-library cohort that combination segfaults:

    Fatal Python error: Segmentation fault
    Current thread ...:
      pandas/core/arrays/string_arrow.py, line 241 in _from_sequence
      pandas/core/dtypes/dtypes.py, line 584 in validate_categories
      pandas/core/arrays/categorical.py, line 773 in from_codes
      anndata/_io/specs/methods.py, line 1197 in read_categorical
      ...
      engine/pipeline.py in _run_one

pandas 3.0 sets `future.infer_string = True`, so every categorical's `categories` Index is built as
an `ArrowStringArray` — a pyarrow allocation, on a thread, in a process that is forking for
subprocesses on other threads.

Measured, and the reason this is a one-line setting rather than a redesign:

    pandas.set_option("mode.string_storage", "python")
    Categorical.from_codes(...).categories.array  ->  StringArray       (numpy object)
    # with "pyarrow":                             ->  ArrowStringArray  (the crashing frame)

Selecting the python storage removes pyarrow from the string path, so the frame that crashes is
never entered. It does not disable Arrow anywhere else, and it does not change which values a
column holds — only how the strings are stored in memory.

IT IS NOT A REGRESSION, AND BISECTING THE TOOL CANNOT FIND IT
------------------------------------------------------------
The crash appears on commits that had run successfully for months, because what changed was not
scQC but the INTERPRETER. `mode.string_storage` defaults to `auto`, and `auto` means *use pyarrow
if it can be imported*. An environment without pyarrow silently uses the numpy backend and never
enters the crashing frame; add pyarrow — as a transitive dependency of something else, without
touching scQC at all — and the same code takes the Arrow path and starts dying under threads.

So this setting is not a workaround. It pins the backend scQC has always effectively used, instead
of leaving it to whichever packages happen to be installed alongside it.

It is also not a resource limit: peak resident memory was under 200 MB against a multi-gigabyte
allocation, and the crash occurred at 2, 4 and 8 cpus alike.

CALL IT ONCE, EARLY, IN EVERY ENTRY POINT. Both settings are process-global and must be in force
before the first pandas import does any string construction, so `main()` calls this before it
parses arguments.
"""

from __future__ import annotations

_HARDENED = False


def harden(*, faulthandler_enabled: bool = True) -> dict:
    """Apply the process-global settings and report what was applied.

    Idempotent: safe to call from several entry points. Returns a dict so a caller can record what
    took effect rather than assuming it — a setting that silently failed to apply is exactly the
    kind that is discovered by a segfault hours into a run.
    """
    global _HARDENED
    applied = {"string_storage": None, "faulthandler": False, "already": _HARDENED}

    if faulthandler_enabled:
        # A NATIVE CRASH MUST NAME A LINE. Without this a SIGSEGV prints one line from the shell
        # and nothing from Python: no frame, no thread, no step. Crashes diagnosed without it cost
        # hours and produced the wrong answer twice; the first run that had it named the offending
        # frame immediately.
        import faulthandler

        try:
            faulthandler.enable()
            applied["faulthandler"] = True
        except Exception:                                                     # noqa: BLE001
            # Never fatal. A run that cannot install a crash handler is still a run.
            pass

    try:
        import pandas as pd

        # Only pandas versions that HAVE the option are touched; on older ones the Arrow string
        # path does not exist and there is nothing to select away from.
        pd.set_option("mode.string_storage", "python")
        applied["string_storage"] = pd.get_option("mode.string_storage")
    except Exception:                                                         # noqa: BLE001
        # A pandas that refuses the option is reported as unset rather than raising: this function
        # makes a crash less likely and must never be the reason a run does not start.
        applied["string_storage"] = None

    _HARDENED = True
    return applied
