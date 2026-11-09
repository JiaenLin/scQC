# Drives one Pipeline's run manifest from several threads at once, in a temporary project; runs no task.
"""The run manifest is written by one thread at a time (found by harness autodev W020).

Running the development suite's two fixture shapes side by side, one run died with
`RuntimeError: dictionary changed size during iteration` and STATUS.json said `failed`. The run
manifest (`work/<key>/state.json`) is one dict written from two places under two different locks:
the run loop records each finished task under its own lock, and a gate records its findings under
`_prov_lock`. Each write serialises the WHOLE dict - so one thread's `json.dump` iterated the dict
while the other thread added a task or a finding to it. Tasks run concurrently, so on any run a gate
and a finishing task can meet; on that run they did, and the run was lost to bookkeeping.

Checked here: many threads recording tasks and findings at once leave no exception, a manifest that
parses, and every task and finding they wrote. On the code before the fix it raises.
"""
from __future__ import annotations

import json
import sys
import tempfile
import threading
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.executor import make_executor  # noqa: E402
from engine.pipeline import Pipeline  # noqa: E402
from engine.task import Status, TaskResult  # noqa: E402

fails: list[str] = []
N_TASKS, N_FINDINGS, ROUNDS = 400, 400, 3

for rnd in range(ROUNDS):
    with tempfile.TemporaryDirectory() as td:
        pipe = Pipeline(project=Path(td) / "proj", mode="evidence",
                        executor=make_executor("local"),
                        samples=[{"sample": "A", "assay": "scrna"}], jobs=4)
        errors: list = []
        loop_lock = threading.Lock()          # the run loop's own lock, as Pipeline.run holds it

        def tasks(offset):
            try:
                for i in range(N_TASKS):
                    r = TaskResult(f"t{offset}/{i}", Status.DONE, step="x", sample="A",
                                   metrics={"i": i, "pad": list(range(20))})
                    with loop_lock:
                        pipe.results_by_key[r.key] = r
                        pipe.state.record(r)
            except Exception as e:                                        # noqa: BLE001
                errors.append(f"recording a task: {type(e).__name__}: {e}")

        def findings(step):
            try:
                for i in range(N_FINDINGS):
                    pipe.record_findings(step, [{"step": step, "check": f"c{i}",
                                                 "severity": "ok", "message": "m", "detail": []}])
            except Exception as e:                                        # noqa: BLE001
                errors.append(f"recording a finding: {type(e).__name__}: {e}")

        threads = [threading.Thread(target=tasks, args=(k,)) for k in range(2)] + \
                  [threading.Thread(target=findings, args=(s,)) for s in ("04_doublets", "05_quality")]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        if errors:
            fails.append(f"round {rnd + 1}: {len(errors)} thread(s) died - {errors[0]}")
            continue
        data = json.loads(pipe.state.path.read_text(encoding="utf-8"))
        n_tasks = len(data.get("tasks") or {})
        n_find = sum(len(v) for v in (data.get("findings_by_step") or {}).values())
        if n_tasks != 2 * N_TASKS or n_find != 2 * N_FINDINGS:
            fails.append(f"round {rnd + 1}: the manifest holds {n_tasks} tasks and {n_find} "
                         f"findings, not {2 * N_TASKS} and {2 * N_FINDINGS}")
print(f"{ROUNDS} rounds of 2 task writers and 2 finding writers at once")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the run manifest is written by one thread at a time, and keeps everything written to it")
