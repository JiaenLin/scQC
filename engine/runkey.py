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
open. So the pipeline's commit is part of the description, and a code change now writes BESIDE the
old run exactly as a parameter change does.

The cost is real and is accepted: any commit invalidates resume for every project, including a
commit that changed only a docstring. Re-doing work that would have produced the same answer is
waste; reusing work that would not is a wrong result nobody can see.

WHAT THIS STILL DOES NOT PROTECT AGAINST

Two runs from the same commit with UNCOMMITTED edits share a key - `dirty` is recorded but a
modified tree has no identity to hash cheaply, so during development the guarantee is weaker than
it looks. Commit before a run whose output anyone will rely on.

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


def code_identity(repo: Path | None = None) -> dict:
    """What version of this pipeline is about to run, as a key component.

    Read through `engine.provenance`, which already resolves a commit on a compute node with no
    git binary. Absence is recorded as absence: a checkout whose commit cannot be read must not
    key the same as one that can, or an unidentifiable code state would silently share a
    directory with an identified one.
    """
    from .provenance import git_provenance

    g = git_provenance(Path(repo) if repo else Path(__file__).resolve().parents[1]) or {}
    commit = g.get("commit")
    ok = isinstance(commit, str) and re.fullmatch(r"[0-9a-f]{7,40}", commit or "")
    return {"commit": commit if ok else "unidentified",
            # True, False or None - and None is not False. An unchecked tree and a verified-clean
            # one are different claims and must not produce the same directory.
            "dirty": g.get("dirty")}


def compute(*, samplesheet_rows, tools: dict, mode: str, extra: dict | None = None,
            code: dict | None = None) -> tuple:
    """(digest, the description it was computed from).

    The description is kept and written out, because a digest nobody can explain is a directory
    name nobody can audit: the point of storing it is that a reader can see WHY two runs went to
    different places.

    `code` defaults to this checkout's identity. It is an argument only so a test can pin it;
    passing a fixed value in production would reinstate the defect it exists to close.
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


def claim(root: Path, digest: str, described: dict) -> Path:
    """Return the directory for this key, recording what it is for. Refuses a mismatch.

    On first use the description is written into the directory. On every later use it is COMPARED,
    so a directory whose contents were produced by different inputs is never written into - the
    one way a content-addressed layout can still overwrite something is if two different runs are
    handed the same name, and this is what notices.
    """
    d = Path(root) / digest
    d.mkdir(parents=True, exist_ok=True)
    m = d / MANIFEST
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
    else:
        m.write_text(json.dumps({"digest": digest, "described": described,
                                 "first_written": time.strftime("%Y-%m-%dT%H:%M:%S%z")},
                                indent=2, default=str) + "\n", encoding="utf-8")
    return d


def index(root: Path, digest: str, described: dict, note: str = "") -> Path:
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
        c = described.get("code") or {}
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
