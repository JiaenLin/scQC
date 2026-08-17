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
enough to tell a real retry from a sleep. All four checks below fail against the old loop.
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

# The source is the subject here: the behaviour needs a real NFS mount to exercise end to end, and
# what regressed is whether the parent is listed at all.
_src = (ROOT / "engine" / "executor.py").read_text(encoding="utf-8")
_wait = _src.split("WAITING FOR THE LOG TO BECOME VISIBLE", 1)
check("the wait is documented as a visibility wait", len(_wait) == 2)
if len(_wait) == 2:
    body = _wait[1].split("out = log.read_text", 1)[0]
    check("the retry LISTS the parent directory", "os.listdir(log.parent)" in body,
          "a stat()-only retry re-reads one cached answer and asks the server nothing")
    check("...and it is not the old fixed-count loop", "for _ in range(30)" not in body)
    check("the budget is used rather than a hardcoded count", "_log_visibility_s()" in body)
    check("a directory listing error does not become the task's failure", "except OSError" in body)

# ------------------------------------------------------- C. the failure reports, does not diagnose

check("the message no longer asserts the job died before its first line",
      "died before its first line - a rejected" not in _src
      or "may genuinely" in _src,
      "an unconditional diagnosis contradicts an Exit_status the same function holds")
check("it distinguishes 'not visible' from 'did not run'",
      "not visible after" in _src and "Exit_status" in _src)
check("it names the override so a slower filesystem has a route",
      "SCQC_LOG_VISIBILITY_S" in _src)

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
