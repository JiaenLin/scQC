"""What a run leaves behind so that a reader who did not watch it can tell four states apart:
ok, partial, refused, died — from the filesystem alone.

    STATUS.json   written FIRST as `partial`, rewritten LAST with the outcome and the products
    RUNNING.txt   written at start; replaced at exit by
    SEALED.txt    exit 0 and every expected product present, or
    FAILED.txt    anything else, naming what is missing

A run that dies leaves `STATUS.json` saying `partial` and `RUNNING.txt` still standing, which is
not the same fact as `refused` and not the same fact as `ok`. The job script's own EXIT trap is
the second witness; this is the first, and the one that exists when there is no job script.

Nothing here decides anything about the data. The shape follows the status contract shared with
the tools this one is orchestrated beside; the field names are the contract's, not this tool's.
"""
from __future__ import annotations

import json
import os
import socket
import time
from pathlib import Path

CONTRACT = "1.0"
STATUSES = ("ok", "partial", "refused", "failed")


def _now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _job() -> dict:
    """Which scheduler job this is, if any. Read from the environment; never invented."""
    for var, name in (("PBS_JOBID", "pbs"), ("SLURM_JOB_ID", "slurm")):
        if os.environ.get(var):
            return {"scheduler": name, "id": os.environ[var], "host": socket.gethostname()}
    return {"scheduler": None, "id": None, "host": socket.gethostname()}


def _rel(path: Path, root: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def products_of(results: Path) -> list:
    """Every file under reports/, tables/ and objects/, relative to the run. What is THERE, not
    what was intended — a status that lists a product it did not write is a lie with a schema."""
    out = []
    for sub in ("reports", "tables", "objects"):
        d = results / sub
        if d.is_dir():
            for p in sorted(d.rglob("*")):
                if p.is_file() and p.suffix in (".json", ".csv", ".html", ".h5ad", ".h5", ".txt"):
                    out.append({"path": _rel(p, results), "bytes": p.stat().st_size})
    return out


def begin(results: Path, *, tool: str, version: str, commit: str | None, state_version: int,
          declaration: dict, headline: str = "started") -> Path:
    """Write STATUS.json as `partial` and RUNNING.txt. Called before any task runs."""
    results = Path(results)
    results.mkdir(parents=True, exist_ok=True)
    rec = {
        "contract": CONTRACT, "tool": tool, "version": version, "commit": commit,
        "state_version": state_version, "status": "partial", "headline": headline,
        "started": _now(), "finished": None, "job": _job(),
        "inputs": [], "products": [], "absent": [], "refusal": None,
        "sees": list(declaration.get("sees", [])),
        "escapes": [], "cannot_show": list(declaration.get("cannot_show", [])),
        "wrapped_versions": {},
    }
    (results / "STATUS.json").write_text(json.dumps(rec, indent=1, default=str) + "\n",
                                         encoding="utf-8")
    (results / "RUNNING.txt").write_text(
        f"started={rec['started']}\njobid={rec['job']['id'] or 'none'}\nhost={rec['job']['host']}\n"
        f"commit={commit or 'unidentified'}\n", encoding="utf-8")
    return results / "STATUS.json"


def finish(results: Path, *, status: str, headline: str, exit_code: int,
           refusal: dict | None = None, escapes: list | None = None,
           wrapped_versions: dict | None = None, absent: list | None = None,
           expected: list | None = None) -> dict:
    """Rewrite STATUS.json with the outcome and replace RUNNING.txt with a seal.

    `expected` names the products the run was supposed to write, relative to the run; a seal is
    SEALED only when the exit code is 0 AND every expected product is present and non-empty.
    A green light indistinguishable from a failure is worse than no light.
    """
    if status not in STATUSES:
        raise ValueError(f"status {status!r} is not one of {STATUSES}")
    results = Path(results)
    p = results / "STATUS.json"
    try:
        rec = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        rec = {"contract": CONTRACT, "started": None}
    products = products_of(results)
    present = {x["path"] for x in products if x["bytes"] > 0}
    missing = [e for e in (expected or []) if e not in present]
    if status == "ok" and missing:
        status = "failed"
        headline = f"{headline}; missing products: {', '.join(missing)}"
    if status == "refused" and not (refusal and refusal.get("fix")):
        refusal = dict(refusal or {})
        refusal.setdefault("reason", headline)
        refusal.setdefault("fix", "read the refusal in reports/report.json; a refusal names what to change")
    rec.update({
        "status": status, "headline": headline, "finished": _now(), "exit": exit_code,
        "products": products, "refusal": refusal if status == "refused" else None,
        "escapes": list(escapes or rec.get("escapes") or []),
        "wrapped_versions": dict(wrapped_versions or rec.get("wrapped_versions") or {}),
        "absent": list(absent or []),
        "missing": missing,
    })
    p.write_text(json.dumps(rec, indent=1, default=str) + "\n", encoding="utf-8")
    sealed = status == "ok" and exit_code == 0 and not missing
    seal = results / ("SEALED.txt" if sealed else "FAILED.txt")
    lines = [f"exit={exit_code}", f"status={status}", f"jobid={(rec.get('job') or {}).get('id') or 'none'}",
             f"commit={rec.get('commit') or 'unidentified'}", f"started={rec.get('started')}",
             f"finished={rec['finished']}", f"host={socket.gethostname()}",
             "products=" + " ".join(x["path"] for x in products)]
    if missing:
        lines.append("missing=" + " ".join(missing))
    seal.write_text("\n".join(lines) + "\n", encoding="utf-8")
    other = results / ("FAILED.txt" if sealed else "SEALED.txt")
    if other.exists():
        other.unlink()
    running = results / "RUNNING.txt"
    if running.exists():
        running.unlink()
    return rec


def read(results: Path) -> dict | None:
    p = Path(results) / "STATUS.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except ValueError:
        return None
