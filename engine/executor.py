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


#: How long to keep asking the server for a finished job's log before giving up, in seconds.
#:
#: It is a CACHE EXPIRY, not a transfer: the file is already written, on storage this host can see.
#: The bound that matters is the NFS client's directory attribute lifetime, `acdirmax`, whose usual
#: default is 60 s - so a 30 s budget was under it even with the directory being re-read, and a
#: task whose log landed just before the orchestrator first looked could not be found in time.
#: 120 s is two of those, and it is only ever paid on a run that is about to fail anyway.
#:
#: Overridable because `acdirmax` is a mount option and some sites set it much higher. A filesystem
#: nobody here has seen is not a reason to make a user patch the source.
_LOG_VISIBILITY_S = 120.0


def _log_visibility_s() -> float:
    """The log-visibility budget, from `SCQC_LOG_VISIBILITY_S` if it is a usable number."""
    raw = os.environ.get("SCQC_LOG_VISIBILITY_S", "").strip()
    if raw:
        try:
            v = float(raw)
            if v > 0:
                return v
        except ValueError:
            # A malformed override falls back rather than stopping the run: this value controls
            # how patient a wait is, and no setting of it is worth refusing a cohort over.
            pass
    return _LOG_VISIBILITY_S


#: Where a PBS configuration lives when the environment does not carry one.
_PBS_CONFS = (
    "/etc/pbs.conf",
    "/cm/local/apps/pbspro/var/etc/pbs.conf",
    "/opt/pbs/etc/pbs.conf",
    "/var/spool/pbs/pbs.conf",
)


def _pbs_env(qsub_path: str | None) -> dict:
    """The environment qsub needs, discovered rather than required of the caller.

    PBS_CONF_FILE, PBS_EXEC and PBS_SERVER are set by the site's environment module. A BATCH JOB
    DOES NOT INHERIT THEM, so a pipeline that submits from inside a job fails with "pbsconf error:
    pbs conf variables not found" unless every user remembers to module-load PBS in a wrapper.

    Requiring that wrapper is not a feature; it is a trap. So: keep whatever the environment
    already has, and otherwise find the conf file and derive PBS_EXEC from where qsub actually
    lives (bin/qsub -> its parent's parent). Anything still missing is left alone, because a wrong
    value is worse than an absent one - qsub's own error names what it needs.
    """
    env = {k: os.environ[k] for k in
           ("PBS_CONF_FILE", "PBS_EXEC", "PBS_SERVER", "PBS_HOME") if k in os.environ}
    if "PBS_CONF_FILE" not in env:
        for c in _PBS_CONFS:
            if Path(c).is_file():
                env["PBS_CONF_FILE"] = c
                break
    if "PBS_EXEC" not in env and qsub_path:
        exec_root = Path(qsub_path).resolve().parent.parent
        if (exec_root / "bin" / "qsub").exists():
            env["PBS_EXEC"] = str(exec_root)
    # The conf file carries the rest; parsing it is how PBS itself bootstraps.
    conf = env.get("PBS_CONF_FILE")
    if conf and Path(conf).is_file():
        for line in Path(conf).read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            k, v = k.strip(), v.strip()
            if k in ("PBS_EXEC", "PBS_SERVER", "PBS_HOME") and k not in env:
                env[k] = v
    return env


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


#: Written by the generated job script as its last line, and the authoritative exit status.
#:
#: `qstat -xf` reports Exit_status only while the job remains in history (24 h here, and not at
#: all on some configurations). A watcher that treats "I could not find out" as zero reports a
#: failed task as a clean one, which is the single worst thing this module can do. The job states
#: its own outcome inside the log, so the answer survives the scheduler forgetting the job.
_RC_MARK = "# scqc-pbs-exit="
_RC_RE = re.compile(r"^# scqc-pbs-exit=(-?\d+)\s*$", re.M)


class PBSExecutor:
    """Submit to PBS Pro and wait.

    Written for a cluster where the login node cannot run the work itself. It polls rather than
    blocking on `qsub -W block=true`, because a blocking qsub that loses its connection leaves a
    job running with nothing watching it.

    THE JOB WRITES ITS OWN LOG; PBS DOES NOT DELIVER IT. `#PBS -o` names a path *on the host that
    submitted the job*, so when the orchestrator is itself a job on one node and a task lands on
    another, PBS must copy the file between two compute nodes after the job ends. On the cluster
    this was found on that copy fails and nothing says so: `Exit_status = 0`, work correct on disk,
    and
    an output file that never appears. It failed intermittently - only for tasks that happened to
    land on a different node from the orchestrator - which reads exactly like a staging delay and
    is not one. Waiting cannot fix it, because there is nothing on the way.

    So the generated script redirects its own output to the log path on shared storage and this
    class reads that file. Nothing is copied between hosts. `#PBS -o/-e` are still set, so the
    scheduler's own noise has somewhere to go, and are no longer read.
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
        self._limits: dict | None = None
        self.qsub = _find_pbs("qsub")
        self.qstat = _find_pbs("qstat")
        # FOUND IS NOT WORKING. qsub needs PBS_CONF_FILE, PBS_EXEC and PBS_SERVER, which on a
        # module-based cluster are set by `module load pbspro` and are NOT inherited by a batch
        # job. Checking only that the binary exists passed, and then every task failed with
        # "pbsconf error: pbs conf variables not found" - one failure per task for one problem
        # with the run.
        # Discovered once, then passed to every PBS call this executor makes.
        self.pbs_env = _pbs_env(self.qsub)
        if self.qsub and self.qstat:
            probe = subprocess.run([self.qstat, "-B"], capture_output=True, text=True,
                                   env={**os.environ, **self.pbs_env})
            if probe.returncode != 0:
                detail = (probe.stderr or probe.stdout).strip()[:200]
                raise SystemExit(
                    "scqc: qsub was found at " + str(self.qsub) + " but PBS is not configured "
                    "in this environment.\n"
                    "    " + detail + "\n"
                    "    PBS_CONF_FILE, PBS_EXEC and PBS_SERVER are set by the site's module, and "
                    "a batch\n"
                    "    job does not inherit them. Add `module load pbspro` (or export "
                    "PBS_CONF_FILE)\n"
                    "    to the script that runs scqc, or use --executor local.")
        if not self.qsub or not self.qstat:
            missing = ", ".join(n for n, v in (("qsub", self.qsub), ("qstat", self.qstat)) if not v)
            raise SystemExit(
                f"scqc: --executor pbs was requested but {missing} could not be found.\n"
                f"    Searched $PATH and {', '.join(_PBS_DIRS)}.\n"
                f"    A batch job does not inherit the login shell's environment, so PBS being on\n"
                f"    PATH where you submitted from does not mean it is on PATH where the\n"
                f"    orchestrator runs. Either module-load PBS inside the job script, or use\n"
                f"    --executor local.")


    def queue_limits(self) -> dict:
        """`resources_max` for the configured queue, as PBS reports it.

        Read from the scheduler, never hardcoded: a ceiling written into this repository is wrong
        on the next cluster, and silently - the job is simply rejected at submission with a
        message about a resource nobody set.
        """
        if self._limits is not None:
            return self._limits
        self._limits = {}
        if not self.queue or not self.qstat:
            return self._limits
        q = subprocess.run([self.qstat, "-Qf", self.queue], capture_output=True, text=True,
                           env={**os.environ, **self.pbs_env})
        for m in re.finditer(r"resources_max\.(\w+)\s*=\s*(\S+)", q.stdout or ""):
            self._limits[m.group(1)] = m.group(2)
        return self._limits

    def check_resources(self, tasks) -> list:
        """Which declared resources exceed the queue, reported BEFORE anything is submitted.

        The graph declares 64 GB on two tasks; a queue capping memory at 50 GB rejects them - at
        task thirty of thirty-seven, after half an hour of work that then has to be repeated. The
        cheap moment to find that is now.
        """
        lim = self.queue_limits()
        out = []
        max_mem = lim.get("mem", "")
        max_cpu = lim.get("ncpus", "")
        mem_gb = int(re.sub(r"[^0-9]", "", max_mem) or 0) if "gb" in max_mem.lower() else 0
        cpu_n = int(re.sub(r"[^0-9]", "", max_cpu) or 0)
        for t in tasks:
            if mem_gb and int(getattr(t, "memory_gb", 0) or 0) > mem_gb:
                out.append(f"{t.key}: asks {t.memory_gb} gb, queue {self.queue} allows {max_mem}")
            if cpu_n and int(getattr(t, "cpus", 0) or 0) > cpu_n:
                out.append(f"{t.key}: asks {t.cpus} cpus, queue {self.queue} allows {max_cpu}")
        return out

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
        # Everything this script emits - the command's output, and the shell's own report of a
        # child that was killed - goes straight to the log on shared storage. `exec` rather than a
        # redirect on the command alone, because "Killed" for an out-of-memory child is printed by
        # bash, not by the child, and that message is the whole diagnosis when it happens.
        lines.append(f"exec > {shlex.quote(str(log))} 2>&1")
        for k, v in (env or {}).items():
            lines.append(f"export {k}={shlex.quote(str(v))}")
        if cwd:
            lines.append(f"cd {shlex.quote(str(cwd))}")
        lines.append("rc=0")
        lines.append(" ".join(shlex.quote(str(c)) for c in cmd) + " || rc=$?")
        # A leading newline so the marker starts a line even when the command's last line had no
        # terminator, which is what makes it findable.
        lines.append(f"printf '\\n{_RC_MARK}%s\\n' \"$rc\"")
        lines.append('exit "$rc"')
        script.write_text("\n".join(lines) + "\n", encoding="utf-8")

        # A log left by an earlier attempt would otherwise be read as this attempt's output if this
        # one dies before its first line - a previous run's success, reported for a job that failed.
        for stale in (log, Path(f"{log}.out"), Path(f"{log}.err")):
            try:
                stale.unlink()
            except FileNotFoundError:
                pass

        sub = subprocess.run([self.qsub, str(script)], capture_output=True, text=True,
                             env={**os.environ, **self.pbs_env})
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
            q = subprocess.run([self.qstat, "-xf", job], capture_output=True, text=True,
                               env={**os.environ, **self.pbs_env})
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

        # WAITING FOR THE LOG TO BECOME VISIBLE, WHICH IS NOT THE SAME AS WAITING FOR IT TO EXIST.
        #
        # The compute node wrote the log itself, on shared storage this host can see, so nothing
        # has to arrive from anywhere. What is being waited out is the CLIENT'S CACHE.
        #
        # `Path.exists()` ON ITS OWN DOES NOT RE-ASK THE SERVER. On NFS a stat() for a name the
        # client has already looked up and NOT found is answered from the negative dentry cache,
        # and that entry is revalidated against the PARENT DIRECTORY's cached attributes rather
        # than by asking about the file. A loop of thirty exists() calls is therefore thirty reads
        # of ONE cached answer - thirty seconds of sleeping, and not one question. Listing the
        # parent forces a READDIR to the server, and that is what makes a retry an actual retry.
        #
        # It needs the orchestrator to be on a different node from the task, which is the only
        # condition under which the negative entry is cached at all - and under `--executor pbs`
        # that is the NORMAL condition. It stays invisible behind a green suite because every
        # local-executor test, and every run on one machine, exercises a filesystem where
        # exists() is always the truth.
        deadline_log = time.time() + _log_visibility_s()
        while True:
            # try/except rather than a guard: the directory is on the same shared storage, and a
            # transient error reading it must not become the task's failure. Whatever it raises,
            # the answer is the same - look again.
            try:
                os.listdir(log.parent)
            except OSError:
                pass
            if log.exists() or time.time() >= deadline_log:
                break
            time.sleep(1.0)
        if not log.exists():
            # WHAT IS KNOWN, NOT A DIAGNOSIS. This message used to assert the job "died before its
            # first line" - a claim about the job, contradicted by evidence this function already
            # holds whenever qstat reported an Exit_status. A guess stated as fact sends the reader
            # hunting a prologue failure that never happened.
            said = (f"The scheduler reports Exit_status = {exit_status.group(1)}, so the job DID "
                    f"run: this is the log not being VISIBLE from this host, not the job failing "
                    f"to start."
                    if exit_status else
                    "The scheduler holds no exit status for it either, so the job may genuinely "
                    "have died before its first line - a rejected resource request, a prologue "
                    "failure, or a path the compute node cannot write.")
            raise TaskFailure(
                f"PBS job {job} finished but {log} is not visible after "
                f"{_log_visibility_s():.0f}s.\n"
                f"    {said}\n"
                f"    The parent directory was re-listed on every attempt, so this is not the "
                f"negative-dentry cache alone. If your filesystem needs longer, raise it with "
                f"SCQC_LOG_VISIBILITY_S.\n"
                f"    The scheduler's own account is in {log}.out and `qstat -xf {job}`.")

        out = log.read_text(encoding="utf-8", errors="replace")
        # The job's own word first; the scheduler's only if the job did not get to speak. Neither
        # present means the outcome is unknown, and unknown is not success.
        marks = list(_RC_RE.finditer(out))
        if marks:
            code = int(marks[-1].group(1))
            out = out[:marks[-1].start()].rstrip("\n") + ("\n" if marks[-1].start() else "")
        elif exit_status:
            code = int(exit_status.group(1))
        else:
            raise TaskFailure(
                f"PBS job {job}: the log stops before the job recorded an exit status, and qstat "
                f"no longer holds one either, so whether it succeeded is unknown.\n"
                f"  log: {log}\n"
                f"    A job killed for walltime or memory ends exactly like this. Refusing rather "
                f"than assuming zero.")
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
