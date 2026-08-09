"""The run manifest: what ran, under what signature, and what it produced.

Resume is the reason this file exists, and resume is where pipelines quietly lie. The two
failure modes it is written against:

  * Re-running a step that did not need it. Expensive, annoying, harmless.
  * SKIPPING a step that did need it. Cheap, invisible, and it silently mixes outputs from two
    different parameter sets into one deliverable.

Only the second is dangerous, so every decision is biased towards re-running. A task is skipped
only when its signature matches AND every promised output still exists AND the manifest records
it as DONE. Anything else runs again.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

from .task import Status, Task, TaskResult

MANIFEST_VERSION = 1


class RunState:
    """A JSON manifest beside the results, written after every task."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.data: dict = {"version": MANIFEST_VERSION, "tasks": {}, "runs": []}
        if self.path.exists():
            try:
                loaded = json.loads(self.path.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError) as e:
                # A corrupt manifest must not be silently discarded: that would turn every
                # completed task into a pending one and re-run a cohort without saying why.
                raise SystemExit(
                    f"scqc: the run manifest at {self.path} could not be read ({e}).\n"
                    f"       Move it aside to start a fresh run, or restore it from a backup.\n"
                    f"       Refusing to overwrite it, because doing so discards the record of\n"
                    f"       what has already been computed.")
            if loaded.get("version") != MANIFEST_VERSION:
                raise SystemExit(
                    f"scqc: manifest at {self.path} is version {loaded.get('version')!r}, "
                    f"this scQC writes version {MANIFEST_VERSION}. Move it aside.")
            self.data = loaded

    # ---------------------------------------------------------------- queries

    def get(self, key: str) -> dict | None:
        return self.data["tasks"].get(key)

    def should_skip(self, task: Task) -> tuple[bool, str]:
        """Can this task be skipped? Returns (skip, reason-if-not)."""
        rec = self.get(task.key)
        if rec is None:
            return False, "never run"
        if rec.get("status") != Status.DONE.value:
            return False, f"previous status was {rec.get('status')}"
        sig = task.signature()
        if rec.get("signature") != sig:
            return False, "inputs or parameters changed"
        missing = [o for o in rec.get("outputs", []) if not Path(o).exists()]
        if missing:
            return False, f"{len(missing)} recorded output(s) no longer on disk"
        # A task that promised outputs and recorded none did not do its job, whatever it said.
        if task.outputs and not rec.get("outputs"):
            return False, "recorded no outputs although the task declares some"
        return True, ""

    # ---------------------------------------------------------------- updates

    def record(self, result: TaskResult) -> None:
        self.data["tasks"][result.key] = result.to_json()
        self.flush()

    def begin_run(self, meta: dict) -> None:
        self.data["runs"].append(meta)
        self.flush()

    def flush(self) -> None:
        """Write atomically. A half-written manifest is worse than none: it parses."""
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=str(self.path.parent), suffix=".tmp")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                json.dump(self.data, fh, indent=2, default=str)
                fh.flush()
                os.fsync(fh.fileno())
            os.replace(tmp, self.path)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise

    # ---------------------------------------------------------------- summary

    def summary(self) -> dict:
        counts: dict = {}
        for rec in self.data["tasks"].values():
            counts[rec.get("status", "?")] = counts.get(rec.get("status", "?"), 0) + 1
        return counts
