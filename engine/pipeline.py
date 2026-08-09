"""The orchestrator: what runs, in what order, and where it stops.

This is the only file that knows the shape of the whole pipeline. Everything else is either a
decision module (which judges numbers) or an adapter (which produces them), and neither knows
that the other exists.

THREE PROPERTIES THIS FILE IS RESPONSIBLE FOR

  1. A refusal stops the run. The gates return findings and a verdict; something has to act on
     that, and it is here. A run that prints REFUSE and carries on has a gate in name only.

  2. Evidence mode cannot remove anything. Not "does not by default" - the apply task is not
     placed in the graph at all, so there is no code path from `--mode evidence` to a deletion.
     A flag that merely defaults to safe is one typo away from not being safe.

  3. A step that did not run is never reported as one that passed. Every task ends in exactly one
     status, blocked tasks are recorded as BLOCKED rather than omitted, and the report payload
     carries them, because a missing section reads as an absent problem.
"""

from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

from .provenance import Provenance
from .state import RunState
from .task import Refusal, Status, Task, TaskFailure, TaskResult

ROOT = Path(__file__).resolve().parents[1]

#: The decision modules, which live in directories that are not valid Python identifiers.
_STEP_MODULES = {
    "verify_raw": ROOT / "lib" / "verify_raw.py",
    "ingest": ROOT / "modules" / "00_ingest" / "ingest.py",
    "ambient": ROOT / "modules" / "01_ambient" / "ambient.py",
    "lr_policy": ROOT / "modules" / "01_ambient" / "lr_policy.py",
    "audit_ambient": ROOT / "modules" / "01_ambient" / "audit_ambient.py",
    "cellcall_gate": ROOT / "modules" / "02_cells" / "cellcall_gate.py",
    "light_floor": ROOT / "modules" / "03_light_floor" / "light_floor.py",
    "doublet": ROOT / "modules" / "04_doublets" / "doublet.py",
    "doublet_health": ROOT / "modules" / "04_doublets" / "doublet_health.py",
    "quality": ROOT / "modules" / "05_quality" / "quality.py",
    "cluster_flags": ROOT / "modules" / "06_cluster_check" / "cluster_flags.py",
    "apply": ROOT / "modules" / "07_apply" / "apply.py",
}
_loaded: dict = {}


def step_module(name: str):
    if name in _loaded:
        return _loaded[name]
    path = _STEP_MODULES[name]
    lib = str(ROOT / "lib")
    if lib not in sys.path:
        sys.path.insert(0, lib)
    spec = importlib.util.spec_from_file_location(f"scqc_step_{name}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    _loaded[name] = mod
    return mod


class Pipeline:
    """Builds the task graph for a project and runs it."""

    def __init__(self, project: Path, mode: str, executor, samples: list[dict],
                 decisions: dict | None = None, force: bool = False):
        if mode not in ("evidence", "apply"):
            raise SystemExit(f"scqc: mode must be 'evidence' or 'apply', got {mode!r}")
        self.project = Path(project)
        self.mode = mode
        self.executor = executor
        self.samples = samples
        self.decisions = decisions or {}
        self.force = force

        self.work = self.project / "work"
        self.results = self.project / "results"
        self.logs = self.project / "logs"
        for d in (self.work, self.results / "tables", self.results / "figures",
                  self.results / "reports", self.results / "objects", self.logs):
            d.mkdir(parents=True, exist_ok=True)

        self.state = RunState(self.work / "state.json")
        self.prov = Provenance(ROOT)
        self.findings: list[dict] = []       # every gate finding, for section 1 of the report
        self.results_by_key: dict = {}

    # ------------------------------------------------------------------ gates

    def gate(self, step: str, findings, verdict: str) -> None:
        """Record a gate's output and stop the run if it refused.

        REVIEW does not stop anything - it means a human must look, not that the run is wrong.
        The distinction is preserved all the way into the report rather than collapsed to
        pass/fail, because collapsing it trains a reader to ignore both.
        """
        for f in findings:
            self.findings.append({
                "step": step,
                "check": getattr(f, "check", "?"),
                "severity": getattr(f, "severity", "ok"),
                "message": getattr(f, "message", str(f)),
                "detail": list(getattr(f, "detail", []) or []),
            })
        if verdict == "REFUSE":
            worst = [f for f in findings if getattr(f, "severity", "") == "REFUSE"]
            raise Refusal(
                f"{step} refused.\n" + "\n".join(f"  {f}" for f in worst)
                + "\n\nNothing after this step has run. Fix the cause, or record a decision that\n"
                  "accepts it, and run again - the completed steps will be reused.")

    # ------------------------------------------------------------------ running

    def run(self, tasks: list[Task]) -> dict:
        by_key = {t.key: t for t in tasks}
        order = _toposort(tasks)
        self.state.begin_run({
            "mode": self.mode,
            "started": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "executor": getattr(self.executor, "name", "?"),
            "samples": [s.get("sample") for s in self.samples],
            "n_tasks": len(tasks),
        })

        stopped: str | None = None
        for key in order:
            task = by_key[key]

            bad = [n for n in task.needs
                   if self.results_by_key.get(n) is None
                   or not self.results_by_key[n].ok]
            if bad:
                r = TaskResult(key, Status.BLOCKED,
                               message=f"upstream did not complete: {', '.join(sorted(bad))}")
                self.results_by_key[key] = r
                self.state.record(r)
                continue

            if not self.force:
                skip, why = self.state.should_skip(task)
                if skip:
                    rec = self.state.get(key) or {}
                    r = TaskResult(key, Status.SKIPPED, signature=rec.get("signature", ""),
                                   outputs=rec.get("outputs", []),
                                   metrics=rec.get("metrics", {}),
                                   versions=rec.get("versions", {}),
                                   message="unchanged since the last completed run")
                    self.results_by_key[key] = r
                    self.state.record(r)
                    print(f"  SKIP    {key}")
                    continue

            print(f"  RUN     {key}")
            started = time.time()
            log = self.logs / f"{key.replace('/', '_')}.log"
            try:
                out = task.fn(task=task, pipeline=self, log=log) or {}
                produced = [str(p) for p in out.get("outputs", [])]
                missing = [p for p in produced if not Path(p).exists()]
                if missing:
                    raise TaskFailure(
                        f"reported outputs that do not exist: {missing}. A step that claims a "
                        f"file it did not write fails here, not three steps later.")
                for name, ver in (out.get("versions") or {}).items():
                    self.prov.observe(name, ver)
                r = TaskResult(key, Status.DONE, signature=task.signature(), outputs=produced,
                               metrics=out.get("metrics", {}), versions=out.get("versions", {}),
                               seconds=time.time() - started, log=str(log))
            except Refusal as e:
                r = TaskResult(key, Status.REFUSED, message=str(e),
                               seconds=time.time() - started, log=str(log))
                self.results_by_key[key] = r
                self.state.record(r)
                stopped = f"{key}: refused"
                print(f"  REFUSE  {key}")
                break
            except TaskFailure as e:
                r = TaskResult(key, Status.FAILED, message=str(e),
                               seconds=time.time() - started, log=str(log))
                self.results_by_key[key] = r
                self.state.record(r)
                stopped = f"{key}: failed"
                print(f"  FAIL    {key}")
                break
            except Exception as e:                                    # noqa: BLE001
                r = TaskResult(key, Status.FAILED,
                               message=f"{type(e).__name__}: {e}",
                               seconds=time.time() - started, log=str(log))
                self.results_by_key[key] = r
                self.state.record(r)
                stopped = f"{key}: {type(e).__name__}"
                print(f"  ERROR   {key}  {type(e).__name__}: {e}")
                break

            self.results_by_key[key] = r
            self.state.record(r)

        # Tasks after the stop never ran. They are recorded, not omitted: a report with a
        # missing section reads as a section with nothing to report.
        for key in order:
            if key not in self.results_by_key:
                r = TaskResult(key, Status.BLOCKED, message="the run stopped before this step")
                self.results_by_key[key] = r
                self.state.record(r)

        return self.payload(stopped)

    # ------------------------------------------------------------------ report payload

    def payload(self, stopped: str | None) -> dict:
        counts: dict = {}
        for r in self.results_by_key.values():
            counts[r.status.value] = counts.get(r.status.value, 0) + 1
        return {
            "mode": self.mode,
            "project": str(self.project),
            "stopped": stopped,
            "task_counts": counts,
            "tasks": {k: r.to_json() for k, r in self.results_by_key.items()},
            "findings": self.findings,
            "provenance": self.prov.snapshot(
                self.project / "decisions.yml"
                if (self.project / "decisions.yml").exists() else None),
            "samples": self.samples,
        }


def _toposort(tasks: list[Task]) -> list[str]:
    """Dependency order, refusing a cycle rather than looping or silently dropping a task."""
    remaining = {t.key: set(t.needs) for t in tasks}
    known = set(remaining)
    for key, needs in remaining.items():
        unknown = needs - known
        if unknown:
            raise SystemExit(f"scqc: task {key!r} depends on unknown task(s): {sorted(unknown)}")
    order: list[str] = []
    while remaining:
        ready = sorted(k for k, v in remaining.items() if not v)
        if not ready:
            raise SystemExit(f"scqc: the task graph has a cycle among {sorted(remaining)}")
        for k in ready:
            order.append(k)
            del remaining[k]
        for v in remaining.values():
            v.difference_update(ready)
    return order
