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

import concurrent.futures as cf
import importlib.util
import sys
import threading
import time
from pathlib import Path

from .executor import bind_resources, clear_resources
from .fs import VISIBILITY_TIMEOUT_S, await_visible
from .provenance import Provenance
from .state import RunState
from .task import Refusal, Status, Task, TaskFailure, TaskResult, first_line

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
    "audit_removal": ROOT / "modules" / "07_apply" / "audit_removal.py",
}
_loaded: dict = {}
# Tasks now run concurrently, so the cache below is read and written from several threads. Without
# this lock two workers can both miss, both exec_module, and both register in sys.modules - the
# modules are idempotent so nothing breaks visibly, which is exactly why it would never be found.
_load_lock = threading.Lock()


def step_module(name: str):
    with _load_lock:
        return _step_module_locked(name)


def _step_module_locked(name: str):
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
                 decisions: dict | None = None, force: bool = False, jobs: int = 1,
                 tools: dict | None = None):
        if mode not in ("evidence", "apply"):
            raise SystemExit(f"scqc: mode must be 'evidence' or 'apply', got {mode!r}")
        self.project = Path(project)
        self.mode = mode
        self.executor = executor
        self.samples = samples
        self.decisions = decisions or {}
        self.force = force
        # How many independent tasks may run at once. 1 reproduces the old serial behaviour
        # exactly, which is what you want when a failure needs to be read in one log.
        self.jobs = max(1, int(jobs))
        self._prov_lock = threading.Lock()
        self.started_at = time.time()

        # OUTPUTS ARE STORED UNDER A NAME DERIVED FROM WHAT PRODUCED THEM.
        #
        # Same samplesheet, same parameters, same mode -> the same directory, which is what lets a
        # re-run reuse completed work. Change a threshold and the digest changes with it, so the
        # new run writes BESIDE the old one rather than over it. Nothing is ever overwritten by a
        # run that would have produced something different, and that is a property of the layout
        # rather than a rule someone has to remember.
        #
        # `results/latest` points at the newest and is the only thing here that does not
        # accumulate, which is exactly why nothing is allowed to depend on it.
        from . import runkey

        self.tools = dict(tools or {})
        self.run_key, self.run_described = runkey.compute(
            samplesheet_rows=samples, tools=self.tools, mode=mode)
        try:
            self.results = runkey.claim(self.project / "results", self.run_key,
                                        self.run_described)
        except runkey.RunKeyMismatch as e:
            raise SystemExit(f"scqc: {e}") from None
        runkey.index(self.project / "results", self.run_key, self.run_described, note=mode)
        self.work = self.project / "work" / self.run_key
        self.logs = self.project / "logs" / self.run_key

        # SCRATCH GOES TO LOCAL DISK, not the shared filesystem.
        #
        # `work/` holds the run manifest and small artefacts and belongs beside the project, on
        # whatever storage the project is on. But the big intermediates are pure inter-process
        # handoffs - a 233 MB MatrixMarket triple written by Python and read straight back by R,
        # which nothing downstream ever opens - and on a cluster the project sits on NFS.
        #
        # Measured: ten doublet tasks writing and re-reading ~1.5 GB across NFS ran at 20-30% CPU
        # over ten minutes each; the same work is 2.2 minutes when nothing else contends. The
        # processes were waiting on the network, not computing, on a node with 4 TB of idle local
        # ext4. Concurrency made that worse rather than better, because ten tasks then contend
        # for one network mount.
        #
        # WHO HAS TO SEE THE SCRATCH DECIDES WHERE IT GOES, and only the executor knows that.
        #
        # With a local executor every task runs in this process's node, so node-local disk is both
        # correct and much faster - on NFS the 233 MB inter-process matrices cost more than the
        # computation they carry.
        #
        # With PBS each task is a job on a DIFFERENT node. TMPDIR inside a job is that job's own
        # /var/tmp, so a path written by the orchestrator is invisible to every task it submits.
        # Choosing node-local scratch under PBS produced exactly that: the orchestrator wrote its
        # params file to /var/tmp/pbs.<its own job>/... and every child job died with
        # FileNotFoundError on a path that existed, on a machine it was not running on.
        #
        # So: shared storage whenever tasks run elsewhere. SCQC_SCRATCH still overrides, because a
        # site may have a genuinely shared fast scratch, but it is not consulted for the default.
        import os as _os
        import tempfile as _tf
        shared_needed = getattr(executor, "name", "local") != "local"
        override = _os.environ.get("SCQC_SCRATCH")
        if override:
            base = override
        elif shared_needed:
            base = str(self.work)            # beside the project: every node can see it
        else:
            base = _os.environ.get("TMPDIR") or _tf.gettempdir()
        try:
            self.scratch = Path(base) / f"scqc_{self.project.name}"
            self.scratch.mkdir(parents=True, exist_ok=True)
            probe = self.scratch / ".writable"
            probe.write_text("x", encoding="utf-8")
            probe.unlink()
        except OSError:
            self.scratch = self.work
        for d in (self.work, self.results / "tables", self.results / "figures",
                  self.results / "reports", self.results / "objects", self.logs):
            d.mkdir(parents=True, exist_ok=True)

        self.state = RunState(self.work / "state.json")
        self.prov = Provenance(ROOT)
        self.results_by_key: dict = {}

        # GATE FINDINGS SURVIVE A SKIP, because a skipped task's findings are still true.
        #
        # They used to live only here, in memory, for the length of one process. A re-run skips
        # completed tasks, so no gate re-ran, so the list was empty - and an empty gate list
        # means "supplied, nothing raised". A report rebuilt on a finished run therefore came
        # out PASS with zero findings over a cohort whose original run raised thirteen REVIEWs,
        # and it looked exactly like a clean run, which is the only kind of wrong a reader
        # cannot see.
        #
        # Keyed by STEP, which is what gate() is given. A step that actually re-runs replaces
        # its own findings the first time it gates in this run; a step that is skipped keeps
        # what it recorded when it last ran.
        self._findings_by_step: dict = dict(self.state.data.get("findings_by_step") or {})
        self._gated_this_run: set = set()

    # ------------------------------------------------------------------ gates

    @property
    def findings(self) -> list:
        """Every gate finding this cohort has, in step order - re-run or restored alike.

        READ-ONLY. It builds a new list on each access, so `pipeline.findings.append(...)`
        appends to a throwaway and vanishes without error. Two callers did exactly that and
        lost four findings between them - including the one saying 21,395 delivered nuclei sat
        inside a flagged cluster. Use `record_findings()`, or `gate()` when there is a verdict.
        """
        return [f for step in sorted(self._findings_by_step)
                for f in self._findings_by_step[step]]

    def record_findings(self, step: str, findings) -> None:
        """Store findings that are not a gate decision - an observation, with nothing to stop.

        Same store and same replace-on-re-run rule as `gate()`; it simply has no verdict to
        act on. Accepts dicts or finding objects, because both are already in use and making
        each caller convert is one more place a message can be dropped.
        """
        with self._prov_lock:
            if step not in self._gated_this_run:
                self._gated_this_run.add(step)
                self._findings_by_step[step] = []
            for f in findings:
                self._findings_by_step[step].append(f if isinstance(f, dict) else {
                    "step": step,
                    "check": getattr(f, "check", "?"),
                    "severity": getattr(f, "severity", "ok"),
                    "message": getattr(f, "message", str(f)),
                    "detail": list(getattr(f, "detail", []) or []),
                })
            self.state.data["findings_by_step"] = self._findings_by_step
            self.state.flush()

    def gate(self, step: str, findings, verdict: str) -> None:
        """Record a gate's output and stop the run if it refused.

        REVIEW does not stop anything - it means a human must look, not that the run is wrong.
        The distinction is preserved all the way into the report rather than collapsed to
        pass/fail, because collapsing it trains a reader to ignore both.
        """
        # UNDER THE LOCK, because a step gates once per LIBRARY and those tasks run concurrently.
        #
        # First gate of this run for this step replaces what the step recorded last time; later
        # gates of the same step append to it. A step whose tasks were all skipped never reaches
        # here and keeps the findings it stored when it did run.
        #
        # Without the lock two libraries of the same step both find the step absent from
        # `_gated_this_run`, and the second one's reset discards what the first had just
        # recorded. It cost five of thirty-five findings on the first run after this store was
        # added - a silent partial loss, in the direction that makes a cohort look cleaner.
        self.record_findings(step, findings)
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

        # INDEPENDENT TASKS RUN CONCURRENTLY.
        #
        # This walked the topological order one task at a time, which made the wall-clock the SUM
        # of every step even though the graph already knows what does not depend on what. Ten
        # libraries of doublet scoring at ~2 minutes each is 22 minutes of a 256-core node running
        # one core. Under PBS it was worse, not better: the executor submits and polls, so serial
        # execution meant submit -> queue -> run -> poll, once per job.
        #
        # Tasks are now released in dependency WAVES: everything whose needs are satisfied starts
        # together, bounded by `jobs`. Threads rather than processes because every task body is
        # waiting on a subprocess or a scheduler, not holding the GIL.
        #
        # Order of RESULTS is unchanged - `order` still drives reporting - so a run's records read
        # the same whether it ran with one worker or thirty.
        # Checked BEFORE the first task, against the queue's own resources_max. A task asking for
        # more than the queue allows is rejected by the scheduler at submission - which on a long
        # graph means half an hour of successful work and then a failure about a limit nobody set.
        checker = getattr(self.executor, "check_resources", None)
        if checker is not None:
            over = checker(tasks)
            if over:
                raise Refusal(
                    "the queue cannot run this graph as declared:\n"
                    + "\n".join(f"    - {o}" for o in over)
                    + "\n    Lower the task's declaration, or choose a queue that allows it "
                      "(`qstat -Qf <queue>` lists resources_max). Refused here rather than at "
                      "submission, so nothing runs that would be thrown away.")

        stopped: str | None = None
        lock = threading.Lock()
        pending = list(order)
        n_workers = max(1, int(self.jobs))

        def _classify(key):
            """BLOCKED / SKIPPED decisions, made under the lock before any work starts."""
            task = by_key[key]
            bad = [n for n in task.needs
                   if self.results_by_key.get(n) is None
                   or not self.results_by_key[n].ok]
            if bad:
                return TaskResult(key, Status.BLOCKED, step=task.step, sample=(task.sample or ""),
                                  message=f"upstream did not complete: {', '.join(sorted(bad))}")
            if not self.force:
                skip, _why = self.state.should_skip(task)
                if skip:
                    rec = self.state.get(key) or {}
                    return TaskResult(key, Status.SKIPPED, step=task.step,
                                      sample=(task.sample or ""),
                                      signature=rec.get("signature", ""),
                                      outputs=rec.get("outputs", []),
                                      metrics=rec.get("metrics", {}),
                                      versions=rec.get("versions", {}),
                                      message="unchanged since the last completed run")
            return None

        def _ready():
            """Keys whose upstreams have all finished, in topological order."""
            out = []
            for key in pending:
                needs = by_key[key].needs or ()
                if all(n in self.results_by_key for n in needs):
                    out.append(key)
            return out

        while pending and stopped is None:
            wave = _ready()
            if not wave:
                break
            batch = []
            for key in wave:
                pre = _classify(key)
                if pre is not None:
                    self.results_by_key[key] = pre
                    # A SKIP IS NOT RECORDED, AND THE MANIFEST IS WHY.
                    #
                    # `should_skip` requires the recorded status to be DONE. Recording the skip
                    # overwrote that DONE with `skipped`, so the NEXT run found a non-done record
                    # and re-ran the task - and the run after that re-ran everything again. A run
                    # could therefore be resumed exactly once, and the second resume silently
                    # became a full re-execution: hours of compute, and every concurrency hazard
                    # in the pipeline exercised again, to rebuild a report.
                    #
                    # There is nothing to record in any case. The record that justified the skip
                    # is already there, and it is the one with the signature, the outputs and the
                    # metrics; the skip result is a copy of it with a worse status. A BLOCKED
                    # result IS recorded, because that is new information about this run.
                    if pre.status is Status.SKIPPED:
                        print(f"  SKIP    {key}")
                    else:
                        self.state.record(pre)
                    pending.remove(key)
                    continue
                batch.append(key)
            if not batch:
                continue
            for key in batch:
                pending.remove(key)
                print(f"  RUN     {key}")
            results = {}
            with cf.ThreadPoolExecutor(max_workers=min(n_workers, len(batch))) as pool:
                futures = {pool.submit(self._run_one, by_key[k]): k for k in batch}
                for fut in cf.as_completed(futures):
                    k = futures[fut]
                    r, note = fut.result()
                    results[k] = (r, note)
                    # RECORDED AS IT LANDS, not when the wave ends. Holding a wave's results until
                    # every task in it finished meant a run killed mid-wave lost work that had
                    # completed minutes earlier, and `state.json` showed nothing while ten tasks
                    # were plainly running. Durability and visible progress both want the write
                    # here; the ORDER a reader sees is restored below.
                    with lock:
                        self.state.record(r)
            # `results_by_key` is populated in topological order, not completion order, so the
            # manifest of a graph reads identically whatever the scheduler did. Safe to do after
            # the wave because every task in a wave is independent by construction - nothing in
            # it consulted another's result.
            for key in [k for k in order if k in results]:
                r, note = results[key]
                with lock:
                    self.results_by_key[key] = r
                if note:
                    print(note)
                if r.status in (Status.REFUSED, Status.FAILED) and stopped is None:
                    stopped = f"{key}: {r.status.value}"

        # Tasks after the stop never ran. They are recorded, not omitted: a report with a
        # missing section reads as a section with nothing to report.
        for key in order:
            if key not in self.results_by_key:
                t = by_key[key]
                r = TaskResult(key, Status.BLOCKED, step=t.step, sample=(t.sample or ""),
                               message="the run stopped before this step")
                self.results_by_key[key] = r
                self.state.record(r)

        return self.payload(stopped)

    def _run_one(self, task):
        """Execute one task body. Returns (TaskResult, line-to-print). Never raises.

        Runs on a worker thread, so it must not touch shared state: the caller records the result
        under the lock, in topological order. The only shared write here is `prov.observe`, which
        is why it is guarded - two tasks reporting the same tool version at once would otherwise
        race on a dict.
        """
        key = task.key
        started = time.time()
        log = self.logs / f"{key.replace('/', '_')}.log"
        # Bound for THIS thread so any executor.shell() an adapter makes inside this task asks the
        # scheduler for what the task declared.
        bind_resources(cpus=task.cpus, memory_gb=task.memory_gb,
                       walltime_h=task.walltime_h, gpu=task.gpu)
        try:
            out = task.fn(task=task, pipeline=self, log=log) or {}
            produced = [str(p) for p in out.get("outputs", [])]
            # Waited for, not merely checked: the task ran on another node and the project is on
            # NFS, so a file written seconds ago can be invisible here for as long as the client
            # caches the directory. See engine/fs.py - a correct run reported as a failure is
            # worse than a slow one.
            missing = await_visible(produced)
            if missing:
                raise TaskFailure(
                    f"reported outputs that do not exist: {missing}. Waited "
                    f"{VISIBILITY_TIMEOUT_S}s for them in case the filesystem was behind. A step "
                    f"that claims a file it did not write fails here, not three steps later.")
            with self._prov_lock:
                for name, ver in (out.get("versions") or {}).items():
                    self.prov.observe(name, ver)
            return TaskResult(
                key, Status.DONE, step=task.step, sample=(task.sample or ""),
                signature=task.signature(), outputs=produced,
                metrics=out.get("metrics", {}), versions=out.get("versions", {}),
                seconds=time.time() - started, log=str(log)), None
        except Refusal as e:
            return TaskResult(key, Status.REFUSED, step=task.step, sample=(task.sample or ""),
                              message=str(e), seconds=time.time() - started,
                              log=str(log)), f"  REFUSE  {key}"
        except TaskFailure as e:
            return TaskResult(key, Status.FAILED, step=task.step, sample=(task.sample or ""),
                              message=str(e), seconds=time.time() - started,
                              log=str(log)), f"  FAIL    {key}"
        except Exception as e:                                        # noqa: BLE001
            return TaskResult(key, Status.FAILED, step=task.step, sample=(task.sample or ""),
                              message=f"{type(e).__name__}: {e}",
                              seconds=time.time() - started,
                              log=str(log)), f"  ERROR   {key}  {type(e).__name__}: {e}"
        finally:
            clear_resources()

    # ------------------------------------------------------------------ report payload

    def report_payload(self, stopped: str | None) -> dict:
        """The orchestrator's state, in the shape report/build.py documents.

        A separate function from payload() on purpose. The report has its OWN schema - run,
        deliverable, gates, parameters, steps, provenance - and a key it does not recognise is
        dropped rather than rejected. Handing it payload() directly produced a document with no
        refusal in it, no stop reason, and nothing anywhere saying a translation had failed. The
        mapping is written out here so a mismatch is a visible edit rather than a silent omission.
        """
        stat = {k: r.status.value for k, r in self.results_by_key.items()}
        refused = [k for k, v in stat.items() if v == "refused"]
        failed = [k for k, v in stat.items() if v == "failed"]
        stopped_after = (sorted(refused) + sorted(failed) or [None])[0]
        return {
            "run": {"project": str(self.project), "mode": self.mode,
                    "invocation": " ".join(sys.argv),
                    "started": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                    # Wall-clock, and the CPU-time the tasks actually consumed. The ratio is the
                    # speed-up concurrency bought; reporting only the first hides whether the
                    # machine was used, and only the second hides how long a person waited.
                    "elapsed_s": round(time.time() - self.started_at, 1),
                    "task_seconds_total": round(
                        sum((r.seconds or 0) for r in self.results_by_key.values()), 1),
                    "jobs": self.jobs},
            # WHAT THE RUN PRODUCED, from the step that produced it. This said "no deliverable
            # was written; this run measured and did not apply" on every run that did not stop -
            # including an apply run that had just written 117,021 nuclei to disk, because the
            # text was chosen from `stopped_after` alone and never consulted step 7 at all. The
            # counts are what the report's first block is built from, so a run's headline figure
            # was a sentence denying the run had produced anything.
            "deliverable": self._deliverable_block(stopped, stopped_after),
            # `None`, not `[]`, when no gate ran in this process AND none was restored from a
            # previous one. The two are different claims: an empty list says the gates ran and
            # raised nothing, which the report reads as a pass, and nothing distinguishes it on
            # the page from a run where no gate was ever evaluated. A resumed run whose state
            # predates finding-persistence is exactly that case.
            "gates": (list(self.findings)
                      if (self.findings or self._gated_this_run) else None),
            "steps": self._step_records(),
            # ASSEMBLED HERE, not by the report task, because the report is built from TWO call
            # sites: the `report` task, and scqc_cli.finish() afterwards - the second exists so a
            # run that STOPPED still leaves a document. The task set this key and finish() then
            # rebuilt the payload without it and wrote over the top, so the section was absent
            # from every report while the code that produced it ran correctly every time. One
            # payload builder, both callers.
            "per_sample": self._per_sample_block(),
            # HOW MANY EACH CRITERION REMOVED, PER LIBRARY. F9 draws the cohort; the cohort is not
            # where the question lives. A criterion taking 2% of one library and 42% of another is
            # a technical gradient sitting exactly where the biology is measured, and one bar
            # cannot show it. Read from the file step 7 wrote, never recomputed here.
            "removal_breakdown": self._removal_breakdown_block(),
            # Figures, for the same reason and the same way. Assembled at the report task, they
            # were assembled correctly, printed as assembled, and then finish() rebuilt the
            # payload without them and wrote a figure-less document over the top - the identical
            # failure this comment already described, reproduced one key lower.
            **self._figures_block(),
            # ONE ROW PER PARAMETER, WITH WHO WAS ALLOWED TO SET IT. `report/build.py` has
            # validated this table since the report existed and nothing ever supplied it, so the
            # section the design document calls "the point of this report" has never rendered.
            "parameters": self._parameters_block(),
            "provenance": {**self.prov.snapshot(
                self.project / "decisions.yml"
                if (self.project / "decisions.yml").exists() else None),
                **self._input_check_block(),
                **self._newest_input_block()},
            "open_items": ([f"{k} did not run: {first_line(self.results_by_key[k].message)}"
                            for k in sorted(stat) if stat[k] == "blocked"][:20]
                           or ["none recorded"]),
        }

    #: The parameters that are true regardless of dataset. Changing one is a code change, which
    #: is what FIXED means - they are listed rather than inferred because a contract nobody wrote
    #: down is not a contract.
    FIXED_PARAMETERS = (
        ("cell caller", "the denoiser's output",
         "pipeline contract: cell selection belongs to QC thresholds and doublet detection, not "
         "to the tool chosen for denoising"),
        ("variable-gene selection", "every gene, no class excluded",
         "no mitochondrial, ribosomal or haemoglobin exclusion term exists; the flagged genes "
         "selection CHOSE are reported instead"),
        ("removal point", "step 7 only",
         "no other step drops an observation, so a removal is recoverable from the ledger"),
    )

    #: Declared parameters, as the CLI names them: (tools key, label, the flag that sets it).
    DECLARED_PARAMETERS = (
        ("dbr", "doublet rate prior (dbr)", "--dbr"),
        ("dbr_sd", "doublet prior uncertainty (dbr.sd)", "--dbr-sd"),
        ("dbr_sd_sweep", "dbr.sd sweep", "--dbr-sd-sweep"),
        ("light_floor", "light floor (UMI)", "--light-floor"),
        ("resolution", "clustering resolution", "--resolution"),
        ("seed", "seed", "--seed"),
        ("device", "device", "--cpu"),
    )

    def _parameters_block(self) -> list | None:
        """Every parameter this run applied, with its CLASS and the basis for it.

        THE CLASS IS THE CLAIM. `DERIVED` says the data produced this number and it is
        reproducible from the data alone; `ADJUDICATED` says a person chose it after seeing the
        evidence, and carries their own words. The two look identical as values on a page and
        carry completely different weight, which is the whole reason this section exists.

        BUILT FROM THE MANIFEST AND THE DECLARATION, NOT FROM A STEP'S METRICS. The class of a
        threshold is decided exactly as `steps._apply_thresholds` decides it - a decisions file
        declares it or the pipeline derived it - so it is re-read from the same two sources here
        rather than carried out of step 7. That keeps this buildable on a resumed run, where
        step 7 is skipped and restores whatever metrics it recorded under an older version.

        Returns None - never `[]` - when the run applied nothing, because the report reads an
        empty list as "there were no parameters" and a null as "nobody wrote them down".
        """
        rows: list = []
        for name, value, basis in self.FIXED_PARAMETERS:
            rows.append({"name": name, "value": value, "class": "FIXED", "basis": basis})

        for key, label, flag in self.DECLARED_PARAMETERS:
            v = self.tools.get(key)
            if v is None or v == "":
                continue
            rows.append({"name": label, "value": v, "class": "DECLARED",
                         "basis": f"supplied as {flag} before the run; no default exists for it"})

        # The thresholds step 7 applies. ADJUDICATED where the decisions file declares one - and
        # the operator's own words travel with it, because a row asserting a human decided
        # without any way to check that assertion is withheld by the report rather than printed.
        from . import steps as _s

        q = (self.decisions or {}).get("quality") or {}
        cohort = (getattr(self.results_by_key.get("05_quality"), "metrics", None) or {})
        applied = (getattr(self.results_by_key.get("07_apply"), "metrics", None) or {})
        for leaf, label, derived_key in (("umi_floor", "UMI floor", "umi_proposed"),
                                         ("gene_floor", "gene floor", "genes_proposed"),
                                         ("mito_ceiling_pct", "mitochondrial ceiling", None)):
            # `_attested`, NOT `isinstance(block, dict)`. A decisions entry with a value and no
            # `approved_by`, or none of the operator's own words, is a number somebody typed and
            # step 7 uses the DERIVED value instead. Reading it as adjudicated here would print a
            # class the run did not apply - the report describing a different filter from the one
            # that ran, which is the one failure this section exists to make impossible. The
            # decision is taken by the same function step 7 takes it with, not by a copy of it.
            block = q.get(leaf)
            if _s._attested(block) is not None:
                rows.append({"name": label, "value": block.get("value"),
                             "class": "ADJUDICATED",
                             "basis": "declared in decisions.yml, overriding what step 5 derived",
                             "verbatim": block.get("verbatim") or "",
                             "decided_by": block.get("approved_by") or "",
                             "decided_on": block.get("approved_on") or ""})
                continue
            if derived_key:
                value = cohort.get(derived_key)
                basis = (f"step 5: the density valley measured per library, proposed as one "
                         f"cohort constant and bounded")
            else:
                value = applied.get("ceilings") or "per library"
                basis = ("step 5: median + k*1.4826*MAD over the barcodes above the light floor, "
                         "with k derived from each library's Tukey fence and the whole bounded")
            if value is None:
                continue
            rows.append({"name": label, "value": value, "class": "DERIVED", "basis": basis})

        return rows or None

    def _newest_input_block(self) -> dict:
        """`{"newest_input": iso}` over the files this run READ, for the freshness comparison.

        WHAT COUNTS AS AN INPUT is the whole question, and getting it wrong in either direction
        breaks the check rather than tightening it. Everything this run wrote is excluded: an
        artifact is trivially older than its own outputs, and a check that compares a report
        against files the report caused would refuse every correct run. What is included is what
        the run read and did not write - the samplesheet, the decisions file, and each library's
        declared matrix and supplied ambient object.

        Absent files are LISTED by `newest_input_time` rather than skipped, so a time computed
        over three inputs where there were five cannot pass as the same measurement.
        """
        from report.build import newest_input_time

        paths = [self.project / "samplesheet.csv"]
        d = self.project / "decisions.yml"
        if d.exists():
            paths.append(d)
        for row in self.samples:
            for col in ("matrix", "ambient_h5"):
                v = str(row.get(col) or "").strip()
                if v:
                    paths.append(v)
        got = newest_input_time(paths)
        if not got.get("newest_input"):
            return {}
        return {"newest_input": got["newest_input"],
                "newest_input_path": got.get("newest_path"),
                "inputs_absent": got.get("absent") or [],
                "inputs_checked": got.get("n_checked")}

    def _input_check_block(self) -> dict:
        """`{"input_check": [...]}` from step 0, or `{}` so the report reports the absence.

        Merged into the provenance snapshot rather than added to it by `Provenance`, because
        `Provenance` describes the ENVIRONMENT - versions, commit, reference, clock - and knows
        nothing about tasks. What each input matrix was verified to be is a result of this run,
        and it comes from the step that verified it.

        An EMPTY dict when no library recorded a verdict, never `[]`: the report distinguishes
        `input_check: null` from `input_check: []`, and only the first is honest about a run that
        never checked.
        """
        rows = []
        for key, r in sorted(self.results_by_key.items()):
            if not key.startswith("00_ingest"):
                continue
            check = (getattr(r, "metrics", None) or {}).get("input_check")
            if isinstance(check, dict) and check.get("name"):
                rows.append(check)
        return {"input_check": rows} if rows else {}

    def _deliverable_block(self, stopped, stopped_after) -> dict:
        """What this run produced: the counts from step 7, or why there are none.

        `n_in` and `n_kept` are handed over rather than a sentence, so the report composes its
        own wording and its headline block has numbers to show. A run that stopped, or one in a
        mode that applies nothing, says so - and says which, because "nothing was removed" and
        "the removal never ran" are different facts about a cohort.
        """
        because = stopped or (
            first_line(self.results_by_key[stopped_after].message)
            if stopped_after and self.results_by_key.get(stopped_after) else None)
        block = {"stopped_after": stopped_after, "stopped_because": because,
                 "unit": "observations"}
        r = self.results_by_key.get("07_apply")
        m = (getattr(r, "metrics", None) or {}) if r is not None else {}
        n_in, n_kept = m.get("n_in"), m.get("n_delivered")
        if n_in is not None and n_kept is not None:
            block["n_in"], block["n_kept"] = n_in, n_kept
            return block
        block["text"] = (f"STOPPED at {stopped_after}" if stopped_after else
                         "no deliverable was written: step 7 recorded no counts, so this run "
                         "measured and did not apply")
        return block

    def _figures_block(self) -> dict:
        """The figure data this run's tables can support, and a reason for every figure they
        cannot. Never raises: a figure that cannot be assembled must not be able to stop the
        document, which is the only record of the run it would have illustrated."""
        from report.collect import collect as _collect

        try:
            figures, notes = _collect(self.results / "tables", samplesheet_rows=self.samples)
        except Exception as exc:                                          # noqa: BLE001
            return {"figures": {},
                    "figure_notes": {"*": f"the figure assembler failed: "
                                          f"{type(exc).__name__}: {exc}"}}
        return {"figures": figures, "figure_notes": notes}

    def _per_sample_block(self) -> dict | None:
        """Every threshold this run derived, one row per library. None if it cannot be built.

        Imported inside the function because `steps` imports this module; the same shape
        `_step_records` uses. None rather than a partial block: the report states an absent
        section as a defect, and a half-filled one would not be stated at all.
        """
        from . import steps as _s

        samples = [s.get("sample") for s in self.samples if s.get("sample")]
        if not samples:
            return None
        try:
            block, _path = _s._per_sample_thresholds(self, samples)
        except Exception:                                             # noqa: BLE001
            return None
        return block

    def _removal_breakdown_block(self) -> dict | None:
        """Per library and per criterion, how many observations left. None if step 7 wrote none.

        Read from `tables/removal_by_criterion.csv` rather than recomputed, for the reason the
        per-library threshold table gives: a report that derives its own numbers can disagree with
        the run it describes, and nothing on the page would say which was right.

        None rather than an empty block. A run that measured and did not apply has no breakdown to
        show, and a table of zeroes would say every criterion removed nothing - which is a claim
        about a removal that never happened.
        """
        import csv as _csv

        path = self.results / "tables" / "removal_by_criterion.csv"
        if not path.exists():
            return None
        try:
            with path.open(encoding="utf-8", newline="") as fh:
                rows = [r for r in _csv.DictReader(fh) if r.get("criterion")]
        except OSError:
            return None
        if not rows:
            return None

        def num(v):
            try:
                return int(str(v).strip())
            except (TypeError, ValueError):
                try:
                    return float(str(v).strip())
                except (TypeError, ValueError):
                    return None

        criteria, samples = [], []
        for r in rows:
            if r["criterion"] not in criteria:
                criteria.append(r["criterion"])
            s = r.get("sample") or ""
            if s and s != "ALL" and s not in samples:
                samples.append(s)
        return {"source": str(path), "criteria": criteria, "samples": samples,
                "rows": [{"sample": r.get("sample"), "criterion": r["criterion"],
                          "n_in": num(r.get("n_in")), "n_fired": num(r.get("n_fired")),
                          "n_sole": num(r.get("n_sole")),
                          "n_removed_any": num(r.get("n_removed_any")),
                          "pct_of_library": num(r.get("pct_of_library")),
                          "pct_sole_of_library": num(r.get("pct_sole_of_library"))}
                         for r in rows]}

    def _step_records(self) -> list:
        """One record per canonical step, aggregating every task that belongs to it.

        The report is organised by step, not by library: a reader asks what happened at the cell
        call, not what happened to library 7. Emitting one record per TASK put every step under a
        key the report did not recognise, and it correctly said so - "no record of this step at
        all" for eight steps that had in fact run. The aggregation is here rather than in the
        report because only the orchestrator knows which tasks belong to which step.

        A step's status is the worst of its tasks, so one refused library cannot be averaged away
        by nine that passed.
        """
        from . import steps as _s

        rank = {"refused": 5, "failed": 4, "blocked": 3, "running": 2, "pending": 1,
                "skipped": 0, "done": 0}
        grouped: dict = {}
        for key, r in sorted(self.results_by_key.items()):
            grouped.setdefault(r.step if hasattr(r, "step") else key.split("/", 1)[0],
                               []).append((key, r))

        out = []
        for step, items in sorted(grouped.items()):
            worst = max(items, key=lambda kr: rank.get(kr[1].status.value, 0))[1]
            what, cannot = _s.step_text(step)
            found, sources = [], []
            for key, r in items:
                sources.extend(str(o) for o in (r.outputs or []))
                for mk, mv in sorted((r.metrics or {}).items()):
                    label = f"{mk} ({r.sample})" if getattr(r, "sample", None) else mk
                    found.append({"label": label, "value": mv,
                                  "source": (list(r.outputs or []) or [None])[0]})
            statuses = sorted({r.status.value for _, r in items})
            out.append({
                "key": step,
                "status": worst.status.value,
                "what_it_does": what,
                "cannot_establish": cannot,
                "found": found,
                "sources": sorted(set(sources)),
                "tasks": [{"key": k, "status": r.status.value,
                           "message": first_line(r.message, 200)}
                          for k, r in items],
                "task_statuses": statuses,
            })
        return out

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
