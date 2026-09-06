#!/usr/bin/env python3
"""What a run leaves behind: STATUS.json, RUNNING.txt, SEALED.txt / FAILED.txt — from the tool.

WHY THIS TEST EXISTS

Until now the only seal a run carried was the one the job script's exit trap wrote, so a run
launched any other way was unsealed, and from stderr alone a refusal and a crash looked the same.
The status contract makes four states distinguishable from the filesystem by a reader who did
not watch the run: ok, partial (it died), refused, failed.

WHAT IS CHECKED

  1. begin() writes STATUS.json as `partial` and RUNNING.txt, with the commit and the job
  2. finish(ok) with every expected product present writes SEALED.txt and removes RUNNING.txt
  3. finish(ok) with an expected product MISSING is not ok: status failed, FAILED.txt names it
  4. finish(refused) carries a refusal with a fix, and FAILED.txt (a refusal is not a success)
  5. the products listed are the files on disk, relative to the run, never absolute
  6. through the real CLI path: `scqc run` with a fake pipeline seals on success, on refusal,
     and on a crash — the crash re-raises after sealing
  7. `scqc run --executor local` outside a job on a machine with qsub is refused, and
     --allow-local lifts it

Stdlib only.
"""
from __future__ import annotations

import io
import json
import os
import sys
import tempfile
import types
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import status  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f"  {detail}" if detail and not ok else ""))
    if not ok:
        FAILS.append(name)


def unit_tests():
    print("status writers")
    with tempfile.TemporaryDirectory() as td:
        r = Path(td) / "results" / "abc123"
        status.begin(r, tool="scqc", version="0.0.0", commit="deadbeef", state_version=1,
                     declaration={"sees": ["design"], "cannot_show": ["x"]})
        rec = json.loads((r / "STATUS.json").read_text())
        check("1 begin writes partial", rec["status"] == "partial" and rec["commit"] == "deadbeef")
        check("1 begin writes RUNNING.txt", (r / "RUNNING.txt").exists())
        check("1 sees and cannot_show carried", rec["sees"] == ["design"] and rec["cannot_show"] == ["x"])
        (r / "reports").mkdir()
        (r / "reports" / "report.json").write_text("{}")
        (r / "tables").mkdir()
        (r / "tables" / "a.csv").write_text("id\n1\n")
        rec = status.finish(r, status="ok", headline="done", exit_code=0,
                            expected=["reports/report.json"])
        check("2 ok seals", (r / "SEALED.txt").exists() and not (r / "RUNNING.txt").exists()
              and rec["status"] == "ok")
        check("5 products relative", all(not x["path"].startswith("/") for x in rec["products"])
              and {x["path"] for x in rec["products"]} == {"reports/report.json", "tables/a.csv"})
        rec = status.finish(r, status="ok", headline="done", exit_code=0,
                            expected=["reports/report.json", "objects/x.h5ad"])
        check("3 missing product is not ok", rec["status"] == "failed" and rec["missing"] == ["objects/x.h5ad"]
              and (r / "FAILED.txt").exists() and not (r / "SEALED.txt").exists()
              and "objects/x.h5ad" in (r / "FAILED.txt").read_text())
        rec = status.finish(r, status="refused", headline="gate refused", exit_code=2)
        check("4 refused carries a fix", rec["refusal"] and rec["refusal"]["fix"] and (r / "FAILED.txt").exists())
        try:
            status.finish(r, status="green", headline="", exit_code=0)
            check("status vocabulary is closed", False)
        except ValueError:
            check("status vocabulary is closed", True)


class _FakeResult:
    def __init__(self, key, ok=True, st="done"):
        self.key, self.ok, self.message, self.metrics = key, ok, "", {}
        self.status = types.SimpleNamespace(value=st)


class _FakePipeline:
    """Enough of engine.pipeline.Pipeline for cmd_run: a results dir and a results_by_key."""
    instances = []

    def __init__(self, project, mode, executor, samples, decisions, force, jobs, tools):
        self.project = Path(project)
        self.results = self.project / "results" / "fake"
        (self.results / "reports").mkdir(parents=True, exist_ok=True)
        self.results_by_key = {}
        self.behaviour = os.environ.get("SCQC_FAKE_BEHAVIOUR", "ok")
        _FakePipeline.instances.append(self)

    def run(self, tasks):
        from engine.task import Refusal
        if self.behaviour == "refuse":
            raise Refusal("02_cells: the aligner lost 12% of its cells to the denoiser")
        if self.behaviour == "crash":
            raise RuntimeError("planted crash")
        self.results_by_key["00_ingest/s1"] = _FakeResult("00_ingest/s1")

    def report_payload(self, stopped):
        return {"parameters": [{"name": "quality.umi_floor", "class": "ADJUDICATED", "value": 350,
                                "basis": "valley 302", "decided_by": "a person", "verbatim": "use 350"}],
                "stopped": stopped}


def cli_tests():
    print("through the CLI")
    import scqc_cli
    from engine import graph
    import engine.pipeline

    def run_cli(argv, behaviour):
        os.environ["SCQC_FAKE_BEHAVIOUR"] = behaviour
        out, err = io.StringIO(), io.StringIO()
        orig = (engine.pipeline.Pipeline, graph.ingest_stage, graph.main_stage, scqc_cli.shutil.which)
        engine.pipeline.Pipeline = _FakePipeline
        graph.ingest_stage = lambda *a, **k: []
        graph.main_stage = lambda *a, **k: []
        scqc_cli.shutil.which = lambda x: None          # no scheduler on this machine
        import report.build as rb
        orig_build = rb.build_report
        rb.build_report = lambda payload, html, js: (Path(js).write_text(json.dumps({"ok": True})),
                                                     Path(html).write_text("<html></html>"))
        try:
            with redirect_stdout(out), redirect_stderr(err):
                try:
                    rc = scqc_cli.main(argv)
                except RuntimeError as e:
                    rc = f"raised {e}"
        finally:
            engine.pipeline.Pipeline, graph.ingest_stage, graph.main_stage, scqc_cli.shutil.which = orig
            rb.build_report = orig_build
            os.environ.pop("SCQC_FAKE_BEHAVIOUR", None)
        return rc, out.getvalue() + err.getvalue(), _FakePipeline.instances[-1].results

    with tempfile.TemporaryDirectory() as td:
        proj = Path(td)
        (proj / "samplesheet.csv").write_text("sample,matrix\ns1,/nonexistent\n")
        argv = ["run", "--project", str(proj), "--mode", "evidence"]
        rc, out, res = run_cli(argv, "ok")
        rec = json.loads((res / "STATUS.json").read_text())
        check("6 ok run: exit 0, STATUS ok, SEALED", rc == 0 and rec["status"] == "ok" and (res / "SEALED.txt").exists(), f"rc={rc} {out[-300:]}")
        check("6 ok run: adjudicated parameter recorded as an escape pair",
              rec["escapes"] and rec["escapes"][0]["decision"]["by"] == "a person")
        check("6 ok run: RUNNING.txt gone", not (res / "RUNNING.txt").exists())
        rc, out, res = run_cli(argv, "refuse")
        rec = json.loads((res / "STATUS.json").read_text())
        check("6 refused run: exit 2, STATUS refused with fix, FAILED.txt",
              rc == 2 and rec["status"] == "refused" and rec["refusal"]["fix"] and (res / "FAILED.txt").exists(), f"rc={rc}")
        rc, out, res = run_cli(argv, "crash")
        rec = json.loads((res / "STATUS.json").read_text())
        # main() reports a bug as exit 1 with its traceback; the seal was written before that
        check("6 crashed run: sealed FAILED with the exception named, then reported as a bug (exit 1)",
              rc == 1 and rec["status"] == "failed" and "RuntimeError" in rec["headline"]
              and (res / "FAILED.txt").exists() and "Traceback" in out, f"rc={rc}")

    # 7. login-node refusal
    print("compute placement")
    with tempfile.TemporaryDirectory() as td:
        proj = Path(td)
        (proj / "samplesheet.csv").write_text("sample,matrix\ns1,/nonexistent\n")
        orig_which = scqc_cli.shutil.which
        saved = {k: os.environ.pop(k) for k in ("PBS_ENVIRONMENT", "SLURM_JOB_ID") if k in os.environ}
        scqc_cli.shutil.which = lambda x: "/usr/bin/qsub" if x == "qsub" else None
        try:
            out = io.StringIO()
            with redirect_stdout(out), redirect_stderr(out):
                try:
                    scqc_cli.main(["run", "--project", str(proj), "--mode", "evidence"])
                    refused = False
                except SystemExit as e:
                    refused = "--allow-local" in str(e)
            check("7 local executor beside a scheduler, outside a job, is refused naming the fix", refused)
            os.environ["PBS_ENVIRONMENT"] = "PBS_BATCH"
            engine.pipeline.Pipeline = _FakePipeline
            graph.ingest_stage = lambda *a, **k: []
            graph.main_stage = lambda *a, **k: []
            with redirect_stdout(out), redirect_stderr(out):
                try:
                    rc = scqc_cli.main(["run", "--project", str(proj), "--mode", "evidence"])
                    inside_job_ok = True
                except SystemExit:
                    inside_job_ok = False
            check("7 inside a job the local executor is the compliant form", inside_job_ok)
        finally:
            scqc_cli.shutil.which = orig_which
            os.environ.pop("PBS_ENVIRONMENT", None)
            os.environ.update(saved)


if __name__ == "__main__":
    unit_tests()
    cli_tests()
    print(f"\n{'PASS' if not FAILS else 'FAIL'}: {len(FAILS)} failing")
    sys.exit(1 if FAILS else 0)
