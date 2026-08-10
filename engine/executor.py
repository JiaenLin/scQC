"""How a task actually runs: in this process, as a subprocess, or on a scheduler.

Every executor obeys the same two rules, and both exist because their opposites produce runs that
look successful and are not:

  1. A non-zero exit is a failure, always. No "the tool prints warnings to stderr so we ignore
     the exit code" special cases. If a tool genuinely misuses exit codes, its adapter states so
     explicitly and checks something else instead.
  2. Output is captured to a file and the path is reported, whether the task succeeded or not.
     A failed step whose log was discarded cannot be diagnosed, and the usual response to that
     is to re-run the whole cohort.
"""

from __future__ import annotations

import os
import os
import re
import shutil
import threading
import shlex
import subprocess
import time
from pathlib import Path
from typing import Protocol

from .task import TaskFailure


# The resources of the task currently executing on THIS thread. Adapters call executor.shell()
# without knowing which task they are inside, and tasks run concurrently, so a shared attribute on
# the executor would be read by the wrong task. A thread-local is read by exactly the task that
# set it.
_CURRENT = threading.local()


def bind_resources(**kw) -> None:
    """Record the running task's resources for any shell() call made on this thread."""
    _CURRENT.res = {k: v for k, v in kw.items() if v is not None}


def clear_resources() -> None:
    _CURRENT.res = {}


def current_resources() -> dict:
    return dict(getattr(_CURRENT, "res", {}) or {})


#: Where PBS lives when it is not on PATH. A batch job does not inherit the login shell's
#: environment, so `qsub` resolving on the submitting host says nothing about the compute node.
_PBS_DIRS = (
    "/opt/pbs/bin", "/opt/pbs/default/bin",
    "/usr/local/pbs/bin", "/usr/pbs/bin",
    "/cm/shared/apps/pbspro/current/bin",
    "/opt/torque/bin", "/usr/local/torque/bin",
)


def _find_pbs(name: str) -> str | None:
    """An absolute path to a PBS command, or None. Never a bare name.

    A bare name is resolved by whatever PATH happens to be, which is exactly the thing that
    differs between the shell a run is launched from and the job it becomes.
    """
    found = shutil.which(name)
    if found:
        return found
    for d in _PBS_DIRS:
        cand = Path(d) / name
        if cand.is_file() and os.access(cand, os.X_OK):
            return str(cand)
    return None


class Executor(Protocol):
    def shell(self, cmd: list[str], log: Path, env: dict | None = None,
              cwd: Path | None = None, timeout_s: int | None = None) -> str:
        """Run a command, capture output to `log`, return the combined output."""
        ...


class LocalExecutor:
    """Run in a subprocess on this machine."""

    name = "local"

    def shell(self, cmd, log: Path, env=None, cwd=None, timeout_s=None) -> str:
        log = Path(log)
        log.parent.mkdir(parents=True, exist_ok=True)
        full = {**os.environ, **(env or {})}
        started = time.time()
        header = (f"$ {' '.join(shlex.quote(str(c)) for c in cmd)}\n"
                  f"# cwd={cwd or Path.cwd()}\n# started={time.strftime('%Y-%m-%dT%H:%M:%S')}\n\n")
        try:
            p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True,
                               env=full, cwd=str(cwd) if cwd else None, timeout=timeout_s,
                               encoding="utf-8", errors="replace")
        except subprocess.TimeoutExpired as e:
            log.write_text(header + f"\nTIMEOUT after {timeout_s}s\n{e.stdout or ''}"
                                    f"{e.stderr or ''}", encoding="utf-8")
            raise TaskFailure(f"timed out after {timeout_s}s: {cmd[0]}") from None
        except FileNotFoundError:
            log.write_text(header + "\nexecutable not found\n", encoding="utf-8")
            raise TaskFailure(
                f"executable not found: {cmd[0]}\n"
                f"  scQC does not fall back to a different tool or a different version - the "
                f"result would not be the one the report describes. Install it, or point the "
                f"adapter at it explicitly.") from None
        out = (p.stdout or "") + (p.stderr or "")
        log.write_text(header + out + f"\n\n# exit={p.returncode} "
                                      f"seconds={time.time() - started:.1f}\n", encoding="utf-8")
        if p.returncode != 0:
            tail = "\n".join(out.strip().splitlines()[-15:])
            raise TaskFailure(f"{Path(str(cmd[0])).name} exited {p.returncode}\n"
                              f"  log: {log}\n  last lines:\n{tail}")
        return out


class PBSExecutor:
    """Submit to PBS Pro and wait.

    Written for a cluster where the login node cannot run the work itself. It polls rather than
    blocking on `qsub -W block=true`, because a blocking qsub that loses its connection leaves a
    job running with nothing watching it.
    """

    name = "pbs"
    _JOBID = re.compile(r"^\s*(\d+\S*)\s*$", re.M)

    def __init__(self, queue: str | None = None, project: str | None = None,
                 poll_s: int = 30, cpus: int = 1, memory_gb: int = 8,
                 walltime_h: int = 4, gpu: bool = False):
        self.queue, self.project, self.poll_s = queue, project, poll_s
        self.cpus, self.memory_gb, self.walltime_h, self.gpu = cpus, memory_gb, walltime_h, gpu
        # Resolved ONCE, here, so a host without PBS is refused before any task runs rather than
        # failing identically ten times with a bare FileNotFoundError that names neither the
        # command nor the reason.
        self.qsub = _find_pbs("qsub")
        self.qstat = _find_pbs("qstat")
        if not self.qsub or not self.qstat:
            missing = ", ".join(n for n, v in (("qsub", self.qsub), ("qstat", self.qstat)) if not v)
            raise SystemExit(
                f"scqc: --executor pbs was requested but {missing} could not be found.\n"
                f"    Searched $PATH and {', '.join(_PBS_DIRS)}.\n"
                f"    A batch job does not inherit the login shell's environment, so PBS being on\n"
                f"    PATH where you submitted from does not mean it is on PATH where the\n"
                f"    orchestrator runs. Either module-load PBS inside the job script, or use\n"
                f"    --executor local.")

    def shell(self, cmd, log: Path, env=None, cwd=None, timeout_s=None) -> str:
        log = Path(log)
        log.parent.mkdir(parents=True, exist_ok=True)
        script = log.with_suffix(".pbs")
        # The RUNNING TASK's resources, not this executor's defaults. The graph declares cpus,
        # memory, walltime and gpu per task and they used to be discarded here: every job was
        # submitted as 1 cpu / 8 GB / 4 h / no GPU, so the denoiser would have run without the
        # GPU it declares and the aligner would have been killed at a third of its walltime.
        res = current_resources()
        cpus = int(res.get("cpus") or self.cpus)
        mem = int(res.get("memory_gb") or self.memory_gb)
        wall = int(res.get("walltime_h") or self.walltime_h)
        gpu = bool(res.get("gpu", self.gpu))
        sel = f"select=1:ncpus={cpus}:mem={mem}gb"
        if gpu:
            sel += ":ngpus=1"
        lines = ["#!/usr/bin/env bash", "#PBS -N scqc", f"#PBS -l {sel}",
                 f"#PBS -l walltime={wall}:00:00",
                 f"#PBS -o {log}.out", f"#PBS -e {log}.err", "#PBS -j oe"]
        if self.queue:
            lines.append(f"#PBS -q {self.queue}")
        if self.project:
            lines.append(f"#PBS -P {self.project}")
        lines.append("set -euo pipefail")
        for k, v in (env or {}).items():
            lines.append(f"export {k}={shlex.quote(str(v))}")
        if cwd:
            lines.append(f"cd {shlex.quote(str(cwd))}")
        lines.append(" ".join(shlex.quote(str(c)) for c in cmd))
        script.write_text("\n".join(lines) + "\n", encoding="utf-8")

        sub = subprocess.run([self.qsub, str(script)], capture_output=True, text=True)
        if sub.returncode != 0:
            raise TaskFailure(f"qsub failed: {(sub.stderr or sub.stdout).strip()}")
        m = self._JOBID.search(sub.stdout.strip())
        if not m:
            raise TaskFailure(f"could not parse a job id from qsub output: {sub.stdout!r}")
        job = m.group(1)

        # The deadline follows the TASK's walltime, not the executor default, or a long job
        # is abandoned by the watcher while PBS is still happily running it.
        deadline = time.time() + (timeout_s or wall * 3600 + 600)
        wait = 2.0
        while True:
            q = subprocess.run([self.qstat, "-xf", job], capture_output=True, text=True)
            text = q.stdout or ""
            state = re.search(r"job_state\s*=\s*(\w)", text)
            exit_status = re.search(r"Exit_status\s*=\s*(-?\d+)", text)
            if state and state.group(1) == "F" or exit_status:
                break
            if q.returncode != 0 and not text.strip():
                # qstat forgets finished jobs on some configurations; fall through to the log.
                break
            if time.time() > deadline:
                raise TaskFailure(f"job {job} still running past the deadline; not killed - "
                                  f"check `qstat {job}` and the log at {log}.out")
            # BACK OFF rather than a flat interval. A fixed 30 s poll makes every task cost at
            # least 30 s however fast it is, so a graph of 37 short tasks spends nineteen minutes
            # asleep. Starting at 2 s and doubling to poll_s costs a handful of extra qstat calls
            # on a long job and returns a short one almost immediately.
            time.sleep(wait)
            wait = min(wait * 2, self.poll_s)

        out = ""
        for cand in (Path(f"{log}.out"), Path(f"{log}.err")):
            if cand.exists():
                out += cand.read_text(encoding="utf-8", errors="replace")
        log.write_text(out, encoding="utf-8")
        code = int(exit_status.group(1)) if exit_status else 0
        if code != 0:
            tail = "\n".join(out.strip().splitlines()[-15:])
            raise TaskFailure(f"PBS job {job} exited {code}\n  log: {log}\n{tail}")
        return out


def make_executor(kind: str, **kw) -> Executor:
    if kind == "local":
        return LocalExecutor()
    if kind == "pbs":
        return PBSExecutor(**kw)
    raise SystemExit(f"scqc: unknown executor {kind!r} (known: local, pbs)")
