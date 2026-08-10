# Gate findings must survive concurrency and must survive a skipped task.
"""Findings are the report's evidence. Losing one is losing the reason a run needed a human.

TWO WAYS THEY WERE LOST, BOTH SILENT, BOTH IN THE DIRECTION THAT MAKES A COHORT LOOK CLEANER

1. A resumed run skips completed tasks, so no gate re-ran and the in-memory list was empty -
   and an empty gate list means "the gates ran and raised nothing", which reads as a pass. A
   report rebuilt on a finished cohort came out PASS with zero findings over a run that had
   raised thirteen REVIEWs.

2. A step gates once per LIBRARY and those tasks run concurrently. Two libraries both found
   the step absent from the "already gated this run" set, and the second one's reset discarded
   what the first had recorded. Five of thirty-five findings went that way.

Neither failure raises anything. The run succeeds, the report renders, and the only symptom is
a number nobody has another copy of.
"""
from __future__ import annotations

import sys
import threading
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.pipeline import Pipeline  # noqa: E402

fails: list[str] = []
print("Gate findings - concurrency and resume")
print("=" * 74)


class _FakeState:
    """Just enough RunState to record into, with no file behind it."""

    def __init__(self, data=None):
        self.data = dict(data or {})

    def flush(self):
        pass


def _pipeline(stored=None) -> Pipeline:
    p = Pipeline.__new__(Pipeline)
    p._prov_lock = threading.Lock()
    p.state = _FakeState({"findings_by_step": stored} if stored else {})
    p._findings_by_step = dict(stored or {})
    p._gated_this_run = set()
    return p


_F = SimpleNamespace(check="a check", severity="ok", message="a message", detail=[])

# ---- 1. every library's findings survive, however they interleave.
N = 24
p = _pipeline()
threads = [threading.Thread(target=p.gate, args=("01_ambient", [_F], "ok")) for _ in range(N)]
for t in threads:
    t.start()
for t in threads:
    t.join()
got = len(p.findings)
print(f"  {'ok    ' if got == N else 'LOST  '} {N} concurrent gates of one step -> {got} findings")
if got != N:
    fails.append(f"{N - got} of {N} findings lost when one step gated from several threads. A "
                 f"step gates once per library and those tasks run concurrently.")

# ---- 2. a step that does not re-run keeps what it recorded last time.
stored = {"01_ambient": [{"step": "01_ambient", "check": "kept", "severity": "REVIEW",
                          "message": "from the previous run", "detail": []}]}
p = _pipeline(stored)
p.gate("07_apply", [_F], "ok")                       # a different step runs; 01_ambient skips
kept = [f for f in p.findings if f["step"] == "01_ambient"]
print(f"  {'ok    ' if kept else 'LOST  '} a skipped step keeps its findings -> {len(kept)} kept")
if not kept:
    fails.append("a step whose tasks were all skipped lost its findings. A resumed run then "
                 "reports no findings, which reads as a pass.")

# ---- 3. a step that DOES re-run replaces its own, rather than doubling them.
p = _pipeline(stored)
p.gate("01_ambient", [_F], "ok")
same = [f for f in p.findings if f["step"] == "01_ambient"]
ok3 = len(same) == 1 and same[0]["message"] != "from the previous run"
print(f"  {'ok    ' if ok3 else 'DRIFT '} a re-run step replaces its own -> {len(same)} finding(s)")
if not ok3:
    fails.append(f"a re-running step did not replace its previous findings ({len(same)} present). "
                 f"Stale findings accumulate and the report describes runs that no longer exist.")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("findings OK - nothing is lost to concurrency, a skip, or a re-run")
