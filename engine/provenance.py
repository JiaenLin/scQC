"""What produced this result.

Every version recorded here is one this run actually obtained by asking the tool. Nothing is
inferred from a lockfile, a conda environment name, or a previous run's manifest: a recorded
version that was not observed is a fabricated provenance record, and it is worse than an absent
one because it cannot be told apart from a real one.

A tool that was never invoked has no version. It is recorded as `not invoked`, never as unknown
and never omitted, so a reader can tell "we did not run CellBender" from "we ran it and failed to
capture which one".
"""

from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
import time
from pathlib import Path

NOT_INVOKED = "not invoked"


def _run(cmd: list[str], timeout: int = 30) -> str | None:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                           encoding="utf-8", errors="replace")
    except (OSError, subprocess.SubprocessError):
        return None
    if p.returncode != 0 and not (p.stdout or p.stderr):
        return None
    return ((p.stdout or "") + (p.stderr or "")).strip() or None


def tool_version(exe: str | Path, args: tuple = ("--version",),
                 pick: int = 0) -> str:
    """Ask an executable for its version. Returns the observed string, or NOT_INVOKED."""
    exe = str(exe)
    out = _run([exe, *args])
    if out is None:
        return NOT_INVOKED
    lines = [ln.strip() for ln in out.splitlines() if ln.strip()]
    return lines[pick] if len(lines) > pick else out.strip()


def git_provenance(repo: Path) -> dict:
    """Commit and cleanliness of the scQC checkout that produced the result."""
    repo = Path(repo)

    def g(*a):
        return _run(["git", "-C", str(repo), *a])

    commit = g("rev-parse", "HEAD")
    if commit is None:
        return {"commit": "not a git checkout", "dirty": None, "describe": None}
    dirty = bool(g("status", "--porcelain"))
    return {
        "commit": commit.split()[0] if commit else None,
        "dirty": dirty,
        "describe": (g("describe", "--tags", "--always", "--dirty") or "").split()[0] or None,
        "branch": (g("rev-parse", "--abbrev-ref", "HEAD") or "").split()[0] or None,
    }


def file_hash(path: Path, algo: str = "sha256", cap_mb: int | None = 64) -> str:
    """Hash a file. Large files are hashed over head+tail+size rather than in full.

    Stated plainly because a partial hash is not a content hash: it detects edits at the ends and
    any size change, and will not detect a change confined to the middle of a multi-gigabyte
    object. Files below the cap are hashed completely and say so.
    """
    p = Path(path)
    if not p.exists():
        return "absent"
    h = hashlib.new(algo)
    size = p.stat().st_size
    cap = None if cap_mb is None else cap_mb * 1024 * 1024
    with p.open("rb") as fh:
        if cap is None or size <= cap:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
            return f"{algo}:{h.hexdigest()}"
        half = cap // 2
        h.update(fh.read(half))
        fh.seek(-half, 2)
        h.update(fh.read(half))
    h.update(str(size).encode())
    return f"{algo}:{h.hexdigest()}:partial(head+tail+size)"


def environment() -> dict:
    return {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "platform": platform.platform(),
        "host": platform.node(),
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }


class Provenance:
    """Accumulates observed versions over a run."""

    def __init__(self, repo: Path):
        self.repo = Path(repo)
        self.tools: dict = {}
        self.inputs: dict = {}

    def observe(self, name: str, version: str) -> None:
        prev = self.tools.get(name)
        if prev and prev != version and prev != NOT_INVOKED and version != NOT_INVOKED:
            # Two different versions of one tool inside a single run makes the cohort
            # incomparable with itself, which is exactly the thing provenance is for.
            raise SystemExit(
                f"scqc: {name} reported two different versions during one run:\n"
                f"       {prev!r} and then {version!r}.\n"
                f"       A cohort processed by two versions of the same tool is not one cohort.")
        self.tools[name] = version

    def record_input(self, label: str, path: Path) -> None:
        self.inputs[label] = {"path": str(path), "hash": file_hash(path)}

    def snapshot(self, decisions_file: Path | None = None,
                 tools_expected: tuple = ()) -> dict:
        """The provenance block, in the shape report/build.py documents.

        The key is `pipeline`, not `scqc`: the report reads that name, and an object whose keys
        the consumer does not recognise is silently dropped rather than rejected - the commit
        simply did not appear in the document, and nothing said so.

        `tools_expected` names every tool this run SHOULD have invoked. The report uses it to
        distinguish a tool that ran from one that never did: without the list, a tool absent from
        `tools` is indistinguishable from a tool whose version capture failed.
        """
        version = (self.repo / "VERSION")
        d = {
            "pipeline": {
                "version": version.read_text(encoding="utf-8").strip()
                if version.exists() else "unknown",
                **git_provenance(self.repo),
            },
            "environment": environment(),
            "tools": dict(sorted(self.tools.items())),
            "tools_expected": sorted(tools_expected),
            "inputs": self.inputs,
            "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        }
        for name in tools_expected:
            d["tools"].setdefault(name, NOT_INVOKED)
        if decisions_file is not None:
            d["decisions"] = {"path": str(decisions_file),
                              "hash": file_hash(Path(decisions_file), cap_mb=None)}
        return d
