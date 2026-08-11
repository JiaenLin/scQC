"""A finished run must be resumable more than once.

THE DEFECT THIS GUARDS. `RunState.should_skip` requires the recorded status to be DONE. The
orchestrator recorded the SKIPPED result over that DONE record, so the next run found a non-done
record and re-ran the task - and the run after that re-ran everything. A run could be resumed
exactly ONCE; the second resume silently became a full re-execution.

Silent is the operative word. It exits 0, it produces the same deliverable, and the only visible
difference is RUN where SKIP should be - in a log of seventy lines that all look alike. On the
cohort that found it, the second resume re-copied 2.5 GB of matrices, re-scored every library for
doublets, re-clustered all ten, and hit a concurrency race that failed one library, all to rebuild
one report.

Run: python tests/test_resume.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.state import RunState  # noqa: E402
from engine.task import Status, Task, TaskResult  # noqa: E402

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'PASS' if cond else 'FAIL'}  {name}"
          + (f"  --  {detail}" if detail and not cond else ""))


with tempfile.TemporaryDirectory() as tmp:
    d = Path(tmp)
    out = d / "made.csv"
    out.write_text("x", encoding="utf-8")
    task = Task(key="05_quality/libA", step="05_quality", sample="libA", fn=None,
                params={"sample": "libA"}, outputs=(str(out),))
    state = RunState(d / "state.json")

    # --- the run that did the work
    state.record(TaskResult("05_quality/libA", Status.DONE, step="05_quality", sample="libA",
                            signature=task.signature(), outputs=[str(out)],
                            metrics={"valleys": {"umi": 350}}))
    skip, why = state.should_skip(task)
    check("the first resume skips", skip, why)

    # --- the first resume. The orchestrator must NOT write the skip over the done record.
    #
    # Reproduced at the level the defect lived at: what `record()` is given decides whether the
    # next resume works, and the guard is at the call site.
    skipped = TaskResult("05_quality/libA", Status.SKIPPED, step="05_quality", sample="libA",
                         signature=task.signature(), outputs=[str(out)],
                         metrics={"valleys": {"umi": 350}},
                         message="unchanged since the last completed run")
    check("a skipped result is not DONE, so recording it would break the next resume",
          skipped.status is not Status.DONE)

    # The defect, demonstrated: record the skip and the next resume re-runs.
    state.record(skipped)
    skip2, why2 = state.should_skip(task)
    check("...recording it DOES break the next resume - the probe can fail",
          not skip2, "recording a skip left the task skippable, so this test proves nothing")
    check("...and the reason names the status", "status" in (why2 or "").lower(), why2)

    # The fix: leave the done record alone. Restored here, then resumed twice more.
    state.record(TaskResult("05_quality/libA", Status.DONE, step="05_quality", sample="libA",
                            signature=task.signature(), outputs=[str(out)],
                            metrics={"valleys": {"umi": 350}}))
    for n in (2, 3, 4):
        skip_n, why_n = state.should_skip(task)
        check(f"resume {n} still skips when the skip was not recorded", skip_n, why_n)

    print("\nwhat must still invalidate a task")
    gone = Task(key="05_quality/libA", step="05_quality", sample="libA", fn=None,
                params={"sample": "libA"}, outputs=(str(out),))
    out.unlink()
    skip_g, why_g = gone.key and state.should_skip(gone)
    check("a recorded output that is gone re-runs the task", not skip_g, why_g)
    out.write_text("x", encoding="utf-8")

    changed = Task(key="05_quality/libA", step="05_quality", sample="libA", fn=None,
                   params={"sample": "libA", "light_floor": 500}, outputs=(str(out),))
    skip_c, why_c = state.should_skip(changed)
    check("a changed parameter re-runs the task", not skip_c, why_c)

    declared = Task(key="05_quality/libA", step="05_quality", sample="libA", fn=None,
                    params={"sample": "libA"}, outputs=(str(d / "never_written.csv"),))
    state.record(TaskResult("05_quality/libA", Status.DONE, step="05_quality", sample="libA",
                            signature=declared.signature(), outputs=[]))
    skip_d, why_d = state.should_skip(declared)
    check("a task that declared outputs and recorded none re-runs", not skip_d, why_d)

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
sys.exit(1 if FAIL else 0)
