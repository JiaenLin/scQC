# Drives PBSExecutor.shell() against a fake qsub and qstat in a temporary directory; submits nothing.
"""A full queue is a state, not a failure (single-cell-harness ADR-0027, Q11).

A PBS queue's per-user cap counts every job the user has in it, from this run and from any other,
so `--jobs` cannot keep under it. On 2026-10-02 two runs submitted in the same minute and two of
one run's doublet tasks were refused "qsub: would exceed queue generic's per-user limit": the run
failed on a condition that cleared minutes later. Checked here, each against the code before it:

  1  a limit refusal is waited out and the job then runs;
  2  any OTHER refusal still fails at once - waiting is for a full queue, not for a wrong job;
  3  a queue that stays full past the budget fails, and says the job itself was fine.
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from engine import executor as ex  # noqa: E402
from engine.task import TaskFailure  # noqa: E402

fails: list[str] = []
ex._QSUB_LIMIT_FIRST_WAIT_S = 0.01


def _executor(td: Path, refusals: int, message: str):
    count = td / "n"
    count.write_text("0")
    qsub = td / "qsub"
    qsub.write_text(f"""#!/bin/bash
n=$(cat {count}); echo $((n+1)) > {count}
if [ "$n" -lt {refusals} ]; then echo "{message}" >&2; exit 38; fi
bash "$1" >/dev/null 2>&1; echo "4242.fake"
""")
    qstat = td / "qstat"
    qstat.write_text("#!/bin/bash\necho ' job_state = F'\necho ' Exit_status = 0'\n")
    for f in (qsub, qstat):
        f.chmod(0o755)
    e = ex.PBSExecutor.__new__(ex.PBSExecutor)
    e.queue, e.project, e.gpu_queue, e.poll_s = None, None, None, 1
    e.cpus, e.memory_gb, e.walltime_h, e.gpu = 1, 1, 1, False
    e._limits_by_queue, e._limits = {}, None
    e.qsub, e.qstat, e.pbs_env = str(qsub), str(qstat), {}
    return e, count


LIMIT = "qsub: would exceed queue generic's per-user limit"
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    e, count = _executor(td, 2, LIMIT)
    try:
        out = e.shell(["echo", "ran"], td / "logs" / "t1.log")
        if "ran" not in out:
            fails.append(f"1: the job did not run after the queue cleared: {out!r}")
        if count.read_text().strip() != "3":
            fails.append(f"1: qsub was called {count.read_text().strip()} times, expected 3")
    except TaskFailure as err:
        fails.append(f"1: a full queue failed the task: {str(err)[:120]}")
print("1. a limit refusal is waited out, and the job then runs")

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    e, count = _executor(td, 5, "qsub: Unknown queue")
    try:
        e.shell(["echo", "ran"], td / "logs" / "t2.log")
        fails.append("2: a refusal that is not a limit was retried into a success")
    except TaskFailure as err:
        if count.read_text().strip() != "1":
            fails.append(f"2: a non-limit refusal was retried ({count.read_text().strip()} calls)")
print("2. any other refusal fails at once")

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    e, count = _executor(td, 10_000, LIMIT)
    os.environ["SCQC_QSUB_LIMIT_WAIT_S"] = "0.05"
    try:
        e.shell(["echo", "ran"], td / "logs" / "t3.log")
        fails.append("3: a queue full past the budget did not fail")
    except TaskFailure as err:
        if "SCQC_QSUB_LIMIT_WAIT_S" not in str(err):
            fails.append(f"3: the failure does not say the queue stayed full: {str(err)[:120]}")
    finally:
        os.environ.pop("SCQC_QSUB_LIMIT_WAIT_S", None)
print("3. a queue full past the budget fails, naming the budget")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("a full queue is waited for; nothing else is")
