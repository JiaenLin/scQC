"""The unit of work, and what a unit of work is allowed to say about itself.

A Task is a pure description: what to run, what it needs, what it will produce. It does not know
how it will be executed - locally, on a scheduler, or not at all because a previous run already
produced its outputs. That separation is what makes resume possible and what keeps the gates
independent of the runner.

THE RULE THIS FILE EXISTS TO ENFORCE

A task reports what happened, never what was hoped for. `TaskResult.outputs` lists files that
were checked to exist after the command returned; a task that claims an output it did not write
fails at the end of `run()`, not later when something downstream reads an empty file. The failure
mode this prevents is the one this project keeps meeting: a step that "succeeded" and produced
nothing, discovered three steps later as a confusing error about a different file.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable, Optional


class Status(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    DONE = "done"
    FAILED = "failed"
    SKIPPED = "skipped"          # inputs unchanged and outputs present
    REFUSED = "refused"          # a gate said no; NOT the same as failed
    BLOCKED = "blocked"          # an upstream task failed or refused


#: Statuses from which the run cannot continue past this task.
TERMINAL_BAD = frozenset({Status.FAILED, Status.REFUSED, Status.BLOCKED})


@dataclass(frozen=True)
class Task:
    """One step, for one sample or for the cohort.

    `key` must be stable across runs: it is what resume matches on. Encoding a timestamp or a
    temporary directory into it silently disables resume, because every key becomes new.
    """

    key: str
    step: str                              # "00_ingest", "01_ambient", ...
    fn: Callable[..., dict]                # returns {"outputs": [...], "metrics": {...}}
    sample: Optional[str] = None           # None for a cohort-level task
    inputs: tuple = ()                     # paths whose content decides staleness
    params: dict = field(default_factory=dict)
    outputs: tuple = ()                    # paths this task promises to write
    needs: tuple = ()                      # keys of tasks that must finish first
    cpus: int = 1
    memory_gb: int = 8
    walltime_h: int = 4
    gpu: bool = False

    def signature(self) -> str:
        """A hash over everything that should invalidate a cached result.

        Input FILES are hashed by (size, mtime_ns) rather than by content: hashing a 40 GB
        matrix on every resume costs more than re-running most steps. That is a deliberate
        trade and it has a real failure mode - a file rewritten with identical size and mtime
        is not noticed - so `scqc run --force` exists, and the run manifest records which
        signature each result was produced under.
        """
        h = hashlib.sha256()
        h.update(self.key.encode())
        h.update(json.dumps(self.params, sort_keys=True, default=str).encode())
        for p in sorted(str(x) for x in self.inputs):
            h.update(p.encode())
            f = Path(p)
            if f.exists():
                st = f.stat()
                h.update(f"{st.st_size}:{st.st_mtime_ns}".encode())
            else:
                h.update(b"<absent>")
        return h.hexdigest()[:16]


@dataclass
class TaskResult:
    key: str
    status: Status
    step: str = ""          # the canonical step this task belongs to; the report groups on it
    sample: str = ""        # empty for a cohort-level task
    signature: str = ""
    outputs: list = field(default_factory=list)
    metrics: dict = field(default_factory=dict)
    versions: dict = field(default_factory=dict)
    message: str = ""
    seconds: float = 0.0
    log: str = ""                          # path to the captured stdout/stderr

    @property
    def ok(self) -> bool:
        return self.status in (Status.DONE, Status.SKIPPED)

    def to_json(self) -> dict:
        d = dict(self.__dict__)
        d["status"] = self.status.value
        d["outputs"] = [str(o) for o in self.outputs]
        return d


def first_line(text: str | None, limit: int | None = None) -> str:
    """The first line of a message, or "". Never raises.

    `"".splitlines()` is `[]`, not `[""]`, so `.splitlines()[0]` raises on a message that was
    never set - and `message` defaults to "" for every DONE and SKIPPED result, because nothing
    went wrong and there was nothing to say. The report is assembled from exactly those results,
    so it could not be built for any run in which something succeeded: `IndexError: list index out
    of range`, naming neither the step nor the field, at the end of a run that had worked.
    """
    lines = (text or "").splitlines()
    first = lines[0] if lines else ""
    return first[:limit] if limit else first


class TaskFailure(RuntimeError):
    """The command did not do what it said. Distinct from a gate refusal."""


class Refusal(RuntimeError):
    """A gate refused. The run stops, and this is a correct outcome, not an error.

    Kept separate from TaskFailure because the two must be reported differently: a failure is
    something to debug, a refusal is something to read. Collapsing them trains a user to treat
    refusals as flaky infrastructure.
    """
