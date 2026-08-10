"""Where a run's outputs go: a directory named after what produced them.

THE RULE THIS FILE EXISTS TO ENFORCE

An output is a function of its inputs, so it is stored under a name derived from those inputs.
Same samplesheet, same parameters, same mode -> same directory, which is what lets a re-run reuse
completed work. Change a threshold and the digest changes with it, so the new run writes beside
the old one instead of over it. Nothing is ever overwritten by a run that would have produced
something different, and that is a property of the layout rather than a rule anyone has to keep.

WHAT GOES INTO THE KEY, AND WHAT DELIBERATELY DOES NOT

The key covers what the OPERATOR supplied: the samplesheet's content, the declared parameters,
and the mode. It does not cover anything the pipeline derives - the count floors, the
mitochondrial ceilings, the doublet calls - because those are a function of the inputs already.
Putting a derived value in the key would be circular: the directory could not be named until the
run that fills it had finished.

It covers the samplesheet's CONTENT rather than its path. Two projects pointing at the same
libraries with the same thresholds are the same computation and should land in the same place; a
samplesheet edited in place is a different one and must not.

WHAT THIS DOES NOT PROTECT AGAINST

A digest is not a guarantee of reproducibility. Two runs with the same key can still differ if a
TOOL changed underneath them - a detector's version, a library's RNG - and tool versions are not
in the key because they are observed during the run, not declared before it. The manifest written
beside each result records the versions that were actually used, which is where a difference of
that kind shows up. The key answers "was this asked for in the same way", not "is this the same
answer".
"""

from __future__ import annotations

import hashlib
import json
import os
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


def compute(*, samplesheet_rows, tools: dict, mode: str, extra: dict | None = None) -> tuple:
    """(digest, the description it was computed from).

    The description is kept and written out, because a digest nobody can explain is a directory
    name nobody can audit: the point of storing it is that a reader can see WHY two runs went to
    different places.
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
    if not idx.exists():
        idx.write_text("digest\tmode\tsamples\tparameters\tfirst_seen\tnote\n", encoding="utf-8")
    seen = {ln.split("\t", 1)[0] for ln in idx.read_text(encoding="utf-8").splitlines()[1:]}
    if digest not in seen:
        params = " ".join(f"{k}={v}" for k, v in sorted(
            (described.get("parameters") or {}).items()) if v != "")
        with idx.open("a", encoding="utf-8") as fh:
            fh.write(f"{digest}\t{described.get('mode', '')}\t"
                     f"{len(described.get('samples') or [])}\t{params}\t"
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
