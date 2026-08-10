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
import re
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


def _git_dir(repo: Path) -> Path | None:
    """The `.git` directory, following the `gitdir:` pointer a worktree leaves behind."""
    g = repo / ".git"
    if g.is_dir():
        return g
    if g.is_file():
        txt = g.read_text(encoding="utf-8", errors="replace").strip()
        if txt.startswith("gitdir:"):
            p = Path(txt.split(":", 1)[1].strip())
            return p if p.is_absolute() else (repo / p)
    return None


def _commit_from_files(repo: Path) -> dict | None:
    """HEAD's commit read out of `.git`, for a host that has no git binary.

    A COMPUTE NODE USUALLY HAS NO GIT. The orchestrator runs as a batch job, `git rev-parse`
    there fails exactly as it fails in a directory that is not a checkout, and the run then
    recorded `commit: not a git checkout` for a checkout - a provenance record that is not
    absent but WRONG, which this module's own header says is the worse of the two.

    The commit is plain text inside `.git` and needs no binary. Cleanliness genuinely does need
    one, so it stays unknown here rather than being guessed: an unmodified tree and an unchecked
    one must not record the same way.
    """
    gd = _git_dir(Path(repo))
    if gd is None or not gd.is_dir():
        return None
    head = gd / "HEAD"
    if not head.is_file():
        return None
    txt = head.read_text(encoding="utf-8", errors="replace").strip()
    if not txt.startswith("ref:"):                       # detached HEAD holds the sha directly
        sha = txt.split()[0] if txt else ""
        return {"commit": sha, "branch": None} if re.fullmatch(r"[0-9a-f]{7,40}", sha) else None
    ref = txt.split(":", 1)[1].strip()
    branch = ref.rsplit("/", 1)[-1]
    loose = gd / ref
    if loose.is_file():
        return {"commit": loose.read_text(encoding="utf-8", errors="replace").strip().split()[0],
                "branch": branch}
    packed = gd / "packed-refs"                          # a ref that has been packed away
    if packed.is_file():
        for line in packed.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.endswith(" " + ref):
                return {"commit": line.split()[0], "branch": branch}
    return None


def git_provenance(repo: Path) -> dict:
    """Commit and cleanliness of the scQC checkout that produced the result."""
    repo = Path(repo)

    def g(*a):
        return _run(["git", "-C", str(repo), *a])

    commit = g("rev-parse", "HEAD")
    if commit is None:
        read = _commit_from_files(repo)
        if read is not None:
            # Which commit ran is recorded; whether it had been edited is NOT, and says so.
            return {"commit": read["commit"], "dirty": None, "describe": None,
                    "branch": read["branch"]}
        if _git_dir(repo) is not None:
            return {"commit": "a checkout whose HEAD could not be resolved", "dirty": None,
                    "describe": None, "branch": None}
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
