"""WHY THIS SUITE EXISTS

A finished PBS job's log was reported as missing while it was on disk, complete, and ending with
this executor's own exit marker. The run stopped and the message sent the reader to hunt a prologue
failure that had not happened.

THE MECHANISM. The wait was

    for _ in range(30):
        if log.exists():
            break
        time.sleep(1.0)

`Path.exists()` is a stat(), and on NFS a stat() for a name the client has already looked up and
NOT found is answered from the negative dentry cache — revalidated against the PARENT DIRECTORY's
cached attributes, never by asking about the file. Thirty stat() calls on one path are thirty reads
of ONE cached answer. The loop slept for thirty seconds and asked nothing. Only a READDIR on the
parent sends the client to the server.

It needs the orchestrator to be on a different node from the task, which is the only condition
under which the negative entry is cached at all — and under `--executor pbs` that is the NORMAL
condition. It survived a green suite because every local-executor test, and every run on one
machine, exercises a filesystem where `exists()` is always the truth.

HOW THIS SUITE MODELS IT, without needing NFS: a fake path whose `exists()` keeps returning False
until the parent directory has been LISTED. That is the cache's one observable property, and it is
enough to tell a real retry from a sleep. Sections B and C run the real wait under that model;
section D shows the model itself separates the two loops.
"""

from __future__ import annotations

import pathlib
import sys
import types

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name if cond else f"{name}  [{detail}]")
    print(f"  {'ok    ' if cond else 'FAILED'} {name}" + (f"   {detail}" if not cond else ""))


from engine import executor as ex  # noqa: E402

# --------------------------------------------------------------------- A. the budget is a budget

check("there is a named visibility budget, not a bare range()",
      isinstance(ex._LOG_VISIBILITY_S, (int, float)) and ex._LOG_VISIBILITY_S > 0)
check("it is above the usual acdirmax default of 60 s", ex._log_visibility_s() >= 60.0,
      f"got {ex._log_visibility_s()}; a budget under acdirmax cannot outlast the cache")

import os  # noqa: E402

os.environ["SCQC_LOG_VISIBILITY_S"] = "7.5"
check("it is overridable for a slower filesystem", ex._log_visibility_s() == 7.5)
os.environ["SCQC_LOG_VISIBILITY_S"] = "not-a-number"
check("a malformed override falls back instead of stopping the run",
      ex._log_visibility_s() == ex._LOG_VISIBILITY_S)
del os.environ["SCQC_LOG_VISIBILITY_S"]

# ------------------------------------------------------- B. the wait LISTS the parent, not stats
#
# BEHAVIOUR, NOT SOURCE TEXT (harness autodev W020, B1-107). This section split executor.py on the
# words of a COMMENT - "WAITING FOR THE LOG TO BECOME VISIBLE" - and looked for `os.listdir(...)`
# in the text after it: a reworded comment failed it, and a listing that never ran passed it. Now
# the real `PBSExecutor.shell()` runs against a stand-in qsub and qstat, under a filesystem model
# in which the finished job's log stays INVISIBLE until its parent directory has been listed - the
# negative dentry cache's one observable property (section D) - so what is checked is what the
# wait does.

import contextlib  # noqa: E402
import pathlib as _pl  # noqa: E402
import tempfile  # noqa: E402

from engine.task import TaskFailure  # noqa: E402


def _executor(td, exit_line=" Exit_status = 0"):
    """A PBSExecutor whose qsub runs the script at once and whose qstat says it finished."""
    qsub = td / "qsub"
    qsub.write_text('#!/bin/bash\nbash "$1" >/dev/null 2>&1; echo "4242.fake"\n')
    qstat = td / "qstat"
    qstat.write_text(f"#!/bin/bash\necho ' job_state = F'\necho '{exit_line}'\n")
    for f in (qsub, qstat):
        f.chmod(0o755)
    e = ex.PBSExecutor.__new__(ex.PBSExecutor)
    e.queue, e.project, e.gpu_queue, e.poll_s = None, None, None, 1
    e.cpus, e.memory_gb, e.walltime_h, e.gpu = 1, 1, 1, False
    e._limits_by_queue, e._limits = {}, None
    e.qsub, e.qstat, e.pbs_env = str(qsub), str(qstat), {}
    return e


@contextlib.contextmanager
def _cached(log, visible_after_listing=True):
    """`log.exists()` answers False until `os.listdir(log.parent)` has run (if ever)."""
    state = {"listed": 0, "stats": 0}
    real_listdir, real_exists = ex.os.listdir, _pl.Path.exists

    def listdir(path="."):
        if _pl.Path(path) == log.parent:
            state["listed"] += 1
        return real_listdir(path)

    def exists(self, *a, **k):
        if self == log:
            state["stats"] += 1
            if not (visible_after_listing and state["listed"]):
                return False
        return real_exists(self, *a, **k)

    ex.os.listdir, _pl.Path.exists = listdir, exists
    try:
        yield state
    finally:
        ex.os.listdir, _pl.Path.exists = real_listdir, real_exists


os.environ["SCQC_LOG_VISIBILITY_S"] = "3"
try:
    with tempfile.TemporaryDirectory() as td:
        td = _pl.Path(td)
        log = td / "logs" / "t.log"
        e = _executor(td)
        got, err = None, None
        with _cached(log) as st:
            try:
                got = e.shell(["echo", "the job wrote this"], log)
            except TaskFailure as x:
                err = x
        check("a log hidden by the cache is read once the parent is listed",
              got is not None and "the job wrote this" in got,
              f"shell() raised: {str(err)[:160]}" if err else f"got {got!r}")
        check("...and the wait listed the parent directory to get it", st["listed"] >= 1,
              "a stat()-only retry re-reads one cached answer and asks the server nothing")

    # ------------------------------------------------- C. the failure reports, does not diagnose
    for exit_line, label in ((" Exit_status = 0", "with"), (" comment = none", "without")):
        with tempfile.TemporaryDirectory() as td:
            td = _pl.Path(td)
            log = td / "logs" / "t.log"
            e = _executor(td, exit_line)
            msg = ""
            with _cached(log, visible_after_listing=False) as st:
                try:
                    e.shell(["echo", "never seen"], log)
                except TaskFailure as x:
                    msg = str(x)
            check(f"a log never visible fails the task ({label} an Exit_status)",
                  "not visible after" in msg, msg[:200] or "shell() returned")
            check(f"...re-listing the parent on every attempt ({label} an Exit_status)",
                  st["listed"] >= 2, f"listed {st['listed']} time(s) in a 3 s budget")
            check(f"...naming the override for a slower filesystem ({label} an Exit_status)",
                  "SCQC_LOG_VISIBILITY_S" in msg, msg[:200])
            if label == "with":
                check("with an Exit_status it says the job DID run, not that it died",
                      "Exit_status = 0" in msg and "DID" in msg
                      and "may genuinely" not in msg, msg[:240])
            else:
                check("without one it says the job MAY have died, not that it did",
                      "may genuinely" in msg, msg[:240])
finally:
    os.environ.pop("SCQC_LOG_VISIBILITY_S", None)

# --------------------------------------------------- D. the cache model: exists() lies until listed


class _CachedPath:
    """A path whose `exists()` is False until `os.listdir(parent)` has been called.

    This is the negative dentry cache reduced to its one observable property.
    """

    def __init__(self):
        self.listed = False
        self.stats = 0

    def exists(self):
        self.stats += 1
        return self.listed


p = _CachedPath()
for _ in range(30):                       # the OLD loop, faithfully
    if p.exists():
        break
check("the old loop never becomes true against the cache model", not p.exists(),
      "if this passes the model is wrong, not the code")
check("...and it asked the server nothing", p.listed is False)

p2 = _CachedPath()
for _ in range(30):                       # the NEW shape
    p2.listed = True                      # <- os.listdir(parent) forces a READDIR
    if p2.exists():
        break
check("listing the parent is what makes it visible", p2.exists() is True)

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
for f in FAIL:
    print(f"  FAILED: {f}")
sys.exit(1 if FAIL else 0)
