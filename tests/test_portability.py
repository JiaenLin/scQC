#!/usr/bin/env python3
"""The leak guard: no site, host, project or cohort identifier anywhere in this repository.

WHY THIS TEST EXISTS

This tool was calibrated on one cohort and is separated from it in code and custody. A cohort
name in a default, a sample name in a help string, a home path in a job script — each is how a
value fitted to one dataset becomes a default for every dataset, and each has happened at least
once here: the calibration cohort's identity was taken back out on 2026-08-10, and a test
fixture carried the study name until a later commit found it.

WHAT IS CHECKED

Every text file in the tree, including tests/, setup/, conf/, docs/ and this file's siblings —
the directories where such strings land last and are looked for least. Two kinds of pattern:

  1. SHAPES of site leakage, spelled here because they are not names: a user home path, a
     site home path, a login-node hostname, a scheduler head-node or job id, an e-mail address,
     a dated scratch directory.
  2. TERMS supplied from OUTSIDE the repository — the file or files named by
     $SCQC_FORBIDDEN_TERMS, several joined as PATH joins them (`os.pathsep`) — so that this guard
     never spells the cohort or the site it guards against. Unset, the scan proves less, and says
     so rather than passing quietly. SET, every file it names must exist and yield at least one
     term, or the guard FAILS: it read one path only, so a ':'-joined pair of lists scanned 0
     terms and printed PASS - a guard proving nothing in the words of one that proved everything.

This file exempts only itself. Stdlib only.
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_EXT = {".py", ".R", ".r", ".sh", ".pbs", ".md", ".yml", ".yaml", ".toml", ".cfg", ".txt",
            ".json", ".csv", ".tsv", ".cff", ".template"}
SKIP_DIRS = {".git", "__pycache__", "_data", ".egg-info"}
SHAPES = [
    (r"(?<![\w/])/(?:Users|home)/[A-Za-z][\w.-]*", "a user home path"),
    (r"/data/[A-Za-z][\w.-]*/home/", "a site home path"),
    (r"\blogin-\d{2}-\d{2}\b", "a login-node hostname"),
    (r"\bhn-\d{2}-\d{2}\b", "a scheduler head-node name"),
    (r"\b\d{6}\.hn-\d{2}-\d{2}\b", "a scheduler job id"),
    (r"[\w.+-]+@[\w-]+\.(?:edu|com|org|sg|ac\.uk)\b", "an e-mail address"),
    (r"scratch/\d{8}__", "a dated scratch directory"),
]


def files():
    for p in ROOT.rglob("*"):
        if any(part in SKIP_DIRS or part.endswith(".egg-info") for part in p.relative_to(ROOT).parts):
            continue
        if p.is_file() and (p.suffix in TEXT_EXT or p.name in ("VERSION", "HEAD.txt")) \
                and p.resolve() != Path(__file__).resolve():
            yield p


class TermSourceError(ValueError):
    """A term list was named and gave nothing to search for."""


def terms() -> list:
    """Every term of every list $SCQC_FORBIDDEN_TERMS names; [] only when it is not set.

    A list that is named and missing, or that holds no term, raises: a typo in the path, or a
    value carrying two lists where one was read, must not quietly shrink what the guard looks for.
    """
    spec = os.environ.get("SCQC_FORBIDDEN_TERMS")
    if spec is None:
        return []
    files = [Path(x) for x in spec.split(os.pathsep) if x.strip()]
    if not files:
        raise TermSourceError("SCQC_FORBIDDEN_TERMS is set and names no file")
    out, bad = [], []
    for f in files:
        if not f.is_file():
            bad.append(f"{f}: no such file")
            continue
        got = [l.strip() for l in f.read_text(encoding="utf-8").splitlines()
               if l.strip() and not l.strip().startswith("#")]
        if not got:
            bad.append(f"{f}: 0 terms")
        out.extend(got)
    if bad:
        raise TermSourceError("; ".join(bad))
    return list(dict.fromkeys(out))


# Author metadata is attribution, not leakage: an e-mail in the citation file names who to cite.
ATTRIBUTION = {"CITATION.cff", "pyproject.toml"}


def main() -> int:
    pats = [(re.compile(p), why) for p, why in SHAPES]
    try:
        site = terms()
    except TermSourceError as e:
        print(f"FAIL: the term source cannot be read as given - {e}. Set SCQC_FORBIDDEN_TERMS to "
              f"one or more files of terms, joined by {os.pathsep!r}; a scan of 0 terms from a "
              f"named source proves nothing and is not reported as a pass.")
        return 1
    pats += [(re.compile(re.escape(t), re.I), f"site term {t!r}") for t in site]
    hits = []
    for p in files():
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for i, line in enumerate(text.splitlines(), 1):
            for rx, why in pats:
                if why == "an e-mail address" and p.name in ATTRIBUTION:
                    continue
                if rx.search(line):
                    hits.append(f"{p.relative_to(ROOT)}:{i}: {why}: {line.strip()[:90]}")
    print(f"scanned {sum(1 for _ in files())} files; {len(site)} site term(s) "
          f"{'from ' + os.environ['SCQC_FORBIDDEN_TERMS'] if site else '(none supplied: set SCQC_FORBIDDEN_TERMS to prove more)'}")
    for h in hits:
        print("  LEAK " + h)
    if hits:
        print(f"FAIL: {len(hits)} leak(s). A portable tool carries no site or cohort identifier; "
              f"move it to the project or generalise it.")
        return 1
    print("PASS: no site or cohort identifier in the tree")
    return 0


if __name__ == "__main__":
    sys.exit(main())
