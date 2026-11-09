# Builds the task graph for a two-library project in a temporary directory and drives PBSExecutor.shell()
# against a stand-in qsub and qstat; submits nothing and reads no matrix.
"""A task's resources can be declared, not only coded (harness autodev W020, B1-043).

Every task's cpus, memory and walltime were fixed in engine/graph.py - 07_apply 8 CPU / 48 GB / 4 h,
the denoiser 64 GB on one GPU queue - and the queue refusal told the reader to "lower the task's
declaration", which only a code edit could do. A cohort of a million cells could not be run.
Checked here:

  1  the register's CHECK: main_stage(..., tools={'resources': {'07_apply': (8, 180, 12)}}) gives
     07_apply memory_gb 180, and the PBS script the executor writes for it carries mem=180gb
  2  a task-key prefix (`07_measure`) reaches only those tasks; undeclared tasks keep their defaults
  3  `--resources` is parsed before any file is opened, and a malformed one is a usage error
  4  `--gpu-queue a,b`: a GPU task goes to the first queue whose limits admit it; the refusal for a
     task no queue admits names `--resources`
"""
from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from engine import executor as ex  # noqa: E402
from engine import graph  # noqa: E402
from engine.pipeline import Pipeline  # noqa: E402
from engine.task import Refusal, Task  # noqa: E402

fails: list[str] = []

ROWS = [{"sample": s, "platform": "10x", "species": "homo_sapiens", "reference": "homo_sapiens/ref",
         "assay": "scrna", "matrix": f"/nonexistent/{s}", "aligner_cells": f"/nonexistent/{s}_f",
         "mt_prefix": "MT-", "ribo_pattern": "^RP[SL]", "dbr": "0.05", "dbr_sd": "0.05"}
        for s in ("A_1", "B_1")]


def tasks_with(td: Path, resources):
    pipe = Pipeline(project=td / "proj", mode="apply", executor=ex.make_executor("local"),
                    samples=ROWS, jobs=1)
    tools = {"dbr": 0.05, **({"resources": resources} if resources else {})}
    return {t.key: t for t in graph.main_stage(pipe, sys.executable, tools,
                                               {s["sample"]: {"mode": "accept"} for s in ROWS})}


# ---- 1-2. the graph
with tempfile.TemporaryDirectory() as td:
    base = tasks_with(Path(td), None)
with tempfile.TemporaryDirectory() as td:
    try:
        got = tasks_with(Path(td), {"07_apply": (8, 180, 12), "07_measure": (2, 16, 1)})
    except Exception as e:                                                # noqa: BLE001
        got = {}
        fails.append(f"1: main_stage with declared resources raised {type(e).__name__}: {e}")
if got:
    a = got.get("07_apply")
    if a is None or (a.cpus, a.memory_gb, a.walltime_h) != (8, 180, 12):
        fails.append(f"1: 07_apply is {a and (a.cpus, a.memory_gb, a.walltime_h)}, not (8, 180, 12)")
    m = got.get("07_measure/A_1")
    if m is None or (m.cpus, m.memory_gb, m.walltime_h) != (2, 16, 1):
        fails.append(f"2: 07_measure/A_1 is {m and (m.cpus, m.memory_gb, m.walltime_h)}")
    moved = [k for k, t in got.items() if not k.startswith("07_")
             and (t.cpus, t.memory_gb, t.walltime_h) != (base[k].cpus, base[k].memory_gb,
                                                         base[k].walltime_h)]
    if moved:
        fails.append(f"2: undeclared tasks changed their resources: {moved}")
print("1-2. declared resources reach the tasks named, and no other")

# ---- 1. the script the executor writes for 07_apply carries the declared memory
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    qsub, qstat = td / "qsub", td / "qstat"
    qsub.write_text('#!/bin/bash\ncp "$1" "$(dirname "$1")/submitted.pbs"; bash "$1" >/dev/null 2>&1; '
                    'echo "4242.fake"\n')
    qstat.write_text("#!/bin/bash\necho ' job_state = F'\necho ' Exit_status = 0'\n")
    for f in (qsub, qstat):
        f.chmod(0o755)
    e = ex.PBSExecutor.__new__(ex.PBSExecutor)
    e.queue, e.project, e.gpu_queue, e.poll_s = None, None, None, 1
    e.gpu_queues = []
    e.cpus, e.memory_gb, e.walltime_h, e.gpu = 1, 1, 1, False
    e._limits_by_queue, e._limits = {}, None
    e.qsub, e.qstat, e.pbs_env = str(qsub), str(qstat), {}
    t = (got or {}).get("07_apply")
    if t is not None:
        ex.bind_resources(cpus=t.cpus, memory_gb=t.memory_gb, walltime_h=t.walltime_h, gpu=t.gpu)
        try:
            e.shell(["echo", "ok"], td / "logs" / "apply.log")
        finally:
            ex.clear_resources()
        script = (td / "logs" / "apply.pbs").read_text()
        if "mem=180gb" not in script or "walltime=12:00:00" not in script:
            fails.append(f"1: the PBS script does not carry the declared resources:\n{script[:300]}")
print("1. the PBS script for 07_apply carries mem=180gb")

# ---- 3. the flag
spec = importlib.util.spec_from_file_location("w020_cli", ROOT / "scqc_cli.py")
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)
parse = getattr(cli, "_resources", None)
if parse is None:
    fails.append("3: scqc run has no --resources")
else:
    if parse("07_apply=8:180:12,06_cluster=8:64:6") != {"07_apply": [8, 180, 12],
                                                        "06_cluster": [8, 64, 6]}:
        fails.append("3: --resources was not parsed as declared")
    for bad in ("07_apply=8:180", "07_apply=8:x:12", "=8:1:1", "07_apply=0:1:1"):
        try:
            parse(bad)
            fails.append(f"3: a malformed --resources {bad!r} was accepted")
        except SystemExit:
            pass
print("3. --resources is parsed before anything runs, and a malformed one refused")

# ---- 4. a comma list of GPU queues, and the refusal's words
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    qstat = td / "qstat"
    qstat.write_text('#!/bin/bash\ncase "$2" in\n  gq_small) echo " resources_max.mem = 64gb" ;;\n'
                     '  gq_big) echo " resources_max.mem = 256gb" ;;\nesac\n')
    qstat.chmod(0o755)
    e = ex.PBSExecutor.__new__(ex.PBSExecutor)
    e.queue, e.project, e.gpu_queue, e.poll_s = "cpuq", None, "gq_small,gq_big", 1
    e.gpu_queues = ["gq_small", "gq_big"]
    e._limits_by_queue, e._limits = {}, None
    e.qstat, e.pbs_env = str(qstat), {}
    if e.queue_for(True, 128, 4) != "gq_big" or e.queue_for(True, 32, 4) != "gq_small":
        fails.append(f"4: a GPU task was not sent to the first queue that admits it "
                     f"({e.queue_for(True, 128, 4)}, {e.queue_for(True, 32, 4)})")
    big = Task(key="01_ambient/A_1", step="01_ambient", fn=lambda *a: {}, cpus=4,
               memory_gb=512, walltime_h=8, gpu=True)
    over = e.check_resources([big])
    if not over:
        fails.append("4: a task no GPU queue admits was not reported")
    with tempfile.TemporaryDirectory() as td2:
        pipe = Pipeline(project=Path(td2) / "p", mode="evidence", executor=e, samples=ROWS, jobs=1)
        try:
            pipe.run([big])
            fails.append("4: the run started a task no queue admits")
        except Refusal as r:
            if "--resources" not in str(r):
                fails.append(f"4: the refusal does not name --resources: {str(r)[:200]}")
print("4. a GPU task goes to the first queue that admits it; the refusal names --resources")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("a task's resources are declared, not coded")
