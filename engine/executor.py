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
import re
import shlex
import subprocess
import time
from pathlib import Path
from typing import Protocol

from .task import TaskFailure


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

    def shell(self, cmd, log: Path, env=None, cwd=None, timeout_s=None) -> str:
        log = Path(log)
        log.parent.mkdir(parents=True, exist_ok=True)
        script = log.with_suffix(".pbs")
        sel = f"select=1:ncpus={self.cpus}:mem={self.memory_gb}gb"
        if self.gpu:
            sel += ":ngpus=1"
        lines = ["#!/usr/bin/env bash", "#PBS -N scqc", f"#PBS -l {sel}",
                 f"#PBS -l walltime={self.walltime_h}:00:00",
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

        sub = subprocess.run(["qsub", str(script)], capture_output=True, text=True)
        if sub.returncode != 0:
            raise TaskFailure(f"qsub failed: {(sub.stderr or sub.stdout).strip()}")
        m = self._JOBID.search(sub.stdout.strip())
        if not m:
            raise TaskFailure(f"could not parse a job id from qsub output: {sub.stdout!r}")
        job = m.group(1)

        deadline = time.time() + (timeout_s or self.walltime_h * 3600 + 600)
        while True:
            q = subprocess.run(["qstat", "-xf", job], capture_output=True, text=True)
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
            time.sleep(self.poll_s)

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
