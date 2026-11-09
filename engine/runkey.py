"""Where a run's outputs go: a directory named after what produced them.

THE RULE THIS FILE EXISTS TO ENFORCE

An output is a function of its inputs, so it is stored under a name derived from those inputs.
Same samplesheet, same parameters, same mode -> same directory, which is what lets a re-run reuse
completed work. Change a threshold and the digest changes with it, so the new run writes beside
the old one instead of over it. Nothing is ever overwritten by a run that would have produced
something different, and that is a property of the layout rather than a rule anyone has to keep.

WHAT GOES INTO THE KEY, AND WHAT DELIBERATELY DOES NOT

The key covers what the OPERATOR supplied: the samplesheet's content, the declared parameters,
the mode - and THE VERSION OF THIS PIPELINE. It does not cover anything the pipeline derives -
the count floors, the mitochondrial ceilings, the doublet calls - because those are a function of
the inputs already. Putting a derived value in the key would be circular: the directory could not
be named until the run that fills it had finished.

It covers the samplesheet's CONTENT rather than its path. Two projects pointing at the same
libraries with the same thresholds are the same computation and should land in the same place; a
samplesheet edited in place is a different one and must not.

WHY THE CODE IS IN THE KEY, ADDED 2026-08-13

It was not, and the omission defeated the rule at the top of this file in the one case that
matters most. A DERIVED threshold is a function of the inputs AND of the code that derives it, so
a change to the derivation produces a different answer from an identical key. The 2026-08-13
mitochondrial change did exactly that: same samplesheet, same parameters, same mode, a materially
different deliverable. Re-running would have resolved to the previous run's directory, found every
task complete, SKIPPED all of them, and republished the old numbers under the new code - exiting
zero, with a report that opens and reads like a correct one.

That is the precise failure this layout exists to prevent, arriving through the one door left
open. So the code is part of the description, and a code change writes BESIDE the old run exactly
as a parameter change does.

THE CODE IS KEYED BY ITS CONTENT, NOT BY ITS COMMIT (harness autodev W020, B1-106)

The key held the commit, so ANY commit - one that changed only a docstring, a README or a test -
invalidated resume for every project, and every rerun of a finished cohort started again from
nothing. The commit protected nothing the code's own bytes do not: what can change an output is
the code that runs, not the history it sits in. So the key now holds a digest of the
OUTPUT-BEARING code (`OUTPUT_BEARING` below - the engine, the decision modules, the adapters, the
report, lib/, the CLI, the environment recipes and VERSION), read from the files themselves: a
docs-only or tests-only commit keeps the key, and a change to any file that runs moves it.

It also closes the hole the commit left. Two runs from one commit with UNCOMMITTED edits shared a
key, because a modified tree had no cheap identity; a digest of the files IS that identity, so an
uncommitted edit to the code that runs now writes beside the run it would otherwise have reused.
The commit and the dirty flag are still recorded - in INPUTS.json beside the description
(`code_records`, one per run that used the directory) and in STATUS.json - because a reader asks
"which commit", but they no longer decide where the outputs go.

WHAT THIS STILL DOES NOT PROTECT AGAINST

A digest is also not a guarantee of reproducibility. Two runs with the same key can still differ
if a TOOL changed underneath them - a detector's version, a library's RNG - and tool versions are
not in the key because they are observed during the run, not declared before it. The manifest
written beside each result records the versions that were actually used, which is where a
difference of that kind shows up. The key answers "was this asked for in the same way, by the same
code", not "is this the same answer".
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import time
from pathlib import Path

#: How much of the digest appears in a path. Twelve hex characters is 48 bits: with a few
#: thousand runs of one project the chance of a collision is negligible, and a directory name a
#: person can read out loud is worth more than the next four characters.
DIGEST_CHARS = 12

MANIFEST = "INPUTS.json"


class RunKeyMismatch(RuntimeError):
    """The directory this key names already holds a run described by different inputs."""


def _canonical(value):
    """A stable JSON form. Sorted keys, and paths as text, so ordering cannot change the digest."""
    return json.dumps(value, sort_keys=True, default=str, separators=(",", ":"))


#: The code whose bytes can change what a run produces, relative to the repository root. A
#: directory means every file under it except `__pycache__/`, compiled `*.pyc`/`*.pyo` and hidden
#: files. NOT here, deliberately: docs/, tests/, setup/ and the top-level markdown - none of them
#: runs inside a pipeline task, so editing one must not cost a project its completed work.
OUTPUT_BEARING = ("engine", "modules", "adapters", "report", "lib", "scqc_cli.py", "conf/env",
                  "VERSION")


def _code_files(root: Path, entry: str) -> list:
    p = root / entry
    if p.is_file():
        return [p]
    if not p.is_dir():
        return []
    out = []
    for f in p.rglob("*"):
        parts = f.relative_to(root).parts
        if not f.is_file() or "__pycache__" in parts or any(x.startswith(".") for x in parts):
            continue
        if f.suffix in (".pyc", ".pyo"):
            continue
        out.append(f)
    return sorted(out, key=lambda f: f.relative_to(root).as_posix())


def code_digest(repo: Path | None = None) -> str:
    """A SHA-256 over the output-bearing code: each file's path and its bytes, nothing else.

    No modification times and no git: two checkouts holding the same files key the same on any
    machine, a compute node with no git binary included. An entry that is absent is part of the
    digest as an absence, so a tree missing a directory cannot key as one that has it.
    """
    root = Path(repo) if repo else Path(__file__).resolve().parents[1]
    h = hashlib.sha256()
    for entry in OUTPUT_BEARING:
        files = _code_files(root, entry)
        if not files:
            h.update(f"{entry}\0<absent>\n".encode("utf-8"))
            continue
        for f in files:
            h.update(f.relative_to(root).as_posix().encode("utf-8") + b"\0")
            h.update(hashlib.sha256(f.read_bytes()).hexdigest().encode("ascii") + b"\n")
    return h.hexdigest()


def code_identity(repo: Path | None = None) -> dict:
    """What version of this pipeline is about to run, as a KEY component: its code's digest.

    The commit is not in it (see "THE CODE IS KEYED BY ITS CONTENT" above); `code_record()` is
    where it is kept.
    """
    return {"digest": code_digest(repo)}


def code_record(repo: Path | None = None) -> dict:
    """The commit and the dirty flag, for the MANIFEST beside the key - not for the key.

    Read through `engine.provenance`, which already resolves a commit on a compute node with no
    git binary. Absence is recorded as absence, and `dirty` stays three-valued: an unchecked tree
    (None) and a verified-clean one (False) are different claims.
    """
    from .provenance import git_provenance

    g = git_provenance(Path(repo) if repo else Path(__file__).resolve().parents[1]) or {}
    commit = g.get("commit")
    ok = isinstance(commit, str) and re.fullmatch(r"[0-9a-f]{7,40}", commit or "")
    return {"commit": commit if ok else "unidentified", "dirty": g.get("dirty")}


def compute(*, samplesheet_rows, tools: dict, mode: str, extra: dict | None = None,
            code: dict | None = None) -> tuple:
    """(digest, the description it was computed from).

    The description is kept and written out, because a digest nobody can explain is a directory
    name nobody can audit: the point of storing it is that a reader can see WHY two runs went to
    different places.

    `code` defaults to this checkout's identity - the digest of its output-bearing code. It is an
    argument only so a test can pin it; passing a fixed value in production would reinstate the
    defect it exists to close.
    """
    rows = [{str(k): ("" if v is None else str(v)) for k, v in sorted(r.items())}
            for r in samplesheet_rows]
    described = {
        "mode": str(mode),
        "samples": [r.get("sample", "") for r in rows],
        "samplesheet": rows,
        # Declared parameters only, and every one of them: a flag that changes the result and is
        # not here would let two different runs share a directory.
        "parameters": {str(k): ("" if v is None else str(v)) for k, v in sorted(tools.items())},
        # A DERIVED threshold is a function of the inputs and of the code that derives it. See
        # this module's header: without this, a change to a derivation republishes the previous
        # run's numbers under the new code and exits zero.
        "code": code if code is not None else code_identity(),
        **({"extra": extra} if extra else {}),
    }
    digest = hashlib.sha256(_canonical(described).encode("utf-8")).hexdigest()[:DIGEST_CHARS]
    return digest, described


def claim(root: Path, digest: str, described: dict, record: dict | None = None) -> Path:
    """Return the directory for this key, recording what it is for. Refuses a mismatch.

    On first use the description is written into the directory. On every later use it is COMPARED,
    so a directory whose contents were produced by different inputs is never written into - the
    one way a content-addressed layout can still overwrite something is if two different runs are
    handed the same name, and this is what notices.

    `record` - the commit and dirty flag (`code_record()`) - is kept BESIDE the description, never
    in it: each run that uses the directory adds its own to `code_records` (a run repeated from the
    same commit adds nothing), so a directory reused across a docs-only commit says which commits
    wrote into it.
    """
    d = Path(root) / digest
    d.mkdir(parents=True, exist_ok=True)
    m = d / MANIFEST
    rec = dict(record) if record else None
    if m.exists():
        try:
            prev = json.loads(m.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            prev = None
        if prev is not None and prev.get("described") != described:
            raise RunKeyMismatch(
                f"{d} already holds a run described by different inputs.\n"
                f"    Its manifest and this run's description disagree, so writing here would "
                f"overwrite results produced from something else. Either the samplesheet was "
                f"edited in place under the same digest, or two descriptions have collided.\n"
                f"    Compare {m} with this run's parameters and move the old directory aside.")
        if prev is not None and rec is not None:
            seen = list(prev.get("code_records") or [])
            if not seen or {k: seen[-1].get(k) for k in ("commit", "dirty")} != rec:
                seen.append({**rec, "at": time.strftime("%Y-%m-%dT%H:%M:%S%z")})
                prev["code_records"] = seen
                m.write_text(json.dumps(prev, indent=2, default=str) + "\n", encoding="utf-8")
    else:
        m.write_text(json.dumps({"digest": digest, "described": described,
                                 **({"code_records": [{**rec, "at": time.strftime(
                                     "%Y-%m-%dT%H:%M:%S%z")}]} if rec is not None else {}),
                                 "first_written": time.strftime("%Y-%m-%dT%H:%M:%S%z")},
                                indent=2, default=str) + "\n", encoding="utf-8")
    return d


def index(root: Path, digest: str, described: dict, note: str = "",
          record: dict | None = None) -> Path:
    """Append this run to the human-readable index beside the results, and point `latest` at it.

    The index is the thing a person actually reads: a directory of digests answers "where is it"
    and not "which one do I want". One line per run, newest last, with the parameters that
    distinguish it.
    """
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    idx = root / "INDEX.tsv"
    # The `code` column joined the index when the commit entered the key. Two rows with identical
    # parameters and different digests are otherwise unexplainable from the index alone, and the
    # index is the file a person actually reads - "which of these do I want" is exactly the
    # question a code change makes hard to answer.
    header = "digest\tmode\tsamples\tcode\tparameters\tfirst_seen\tnote\n"
    if not idx.exists():
        idx.write_text(header, encoding="utf-8")
    lines = idx.read_text(encoding="utf-8").splitlines()
    if lines and lines[0] != header.rstrip("\n"):
        # An index written before the column existed. Rewritten in place rather than appended to
        # with a different shape: a TSV whose rows have two widths is one no reader parses twice.
        old = [ln.split("\t") for ln in lines[1:] if ln.strip()]
        idx.write_text(header + "".join(
            "\t".join(r[:3] + ["(pre-2026-08-13, code not keyed)"] + r[3:]) + "\n"
            for r in old), encoding="utf-8")
        lines = idx.read_text(encoding="utf-8").splitlines()
    seen = {ln.split("\t", 1)[0] for ln in lines[1:]}
    if digest not in seen:
        params = " ".join(f"{k}={v}" for k, v in sorted(
            (described.get("parameters") or {}).items()) if v != "")
        # The commit as the RECORD gives it; a description written before the code was keyed by
        # content carried it itself.
        c = record or described.get("code") or {}
        code_txt = (str(c.get("commit") or "?")[:12]
                    + ("+dirty" if c.get("dirty") else
                       "" if c.get("dirty") is False else "+unchecked"))
        with idx.open("a", encoding="utf-8") as fh:
            fh.write(f"{digest}\t{described.get('mode', '')}\t"
                     f"{len(described.get('samples') or [])}\t{code_txt}\t{params}\t"
                     f"{time.strftime('%Y-%m-%dT%H:%M:%S%z')}\t{note}\n")
    # A convenience pointer, and only that. It is rewritten every run, so it is the one thing
    # here that does not accumulate - which is why nothing is allowed to depend on it.
    link = root / "latest"
    try:
        if link.is_symlink() or link.exists():
            if link.is_symlink() or link.is_file():
                link.unlink()
        os.symlink(digest, link, target_is_directory=True)
    except (OSError, NotImplementedError, AttributeError):
        # Windows without developer mode, or a filesystem with no symlinks. A text pointer says
        # the same thing and never fails.
        (root / "latest.txt").write_text(digest + "\n", encoding="utf-8")
    return idx
