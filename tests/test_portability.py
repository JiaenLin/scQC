#!/usr/bin/env python3
"""The leak guard: no site, host, project or cohort identifier anywhere in this repository.

WHY THIS TEST EXISTS

This tool was calibrated on one cohort and is separated from it in code and custody. A cohort
name in a default, a sample name in a help string, a home path in a job script — each is how a
value fitted to one dataset becomes a default for every dataset, and each has happened at least
once here: the calibration cohort's identity was taken back out on 2026-08-10, and a test
fixture carried the study name until a later commit found it.

WHAT IS CHECKED

Every file git lists - tracked, and new files git does not ignore - READ BY ITS CONTENT, whatever
its name or extension: a file is text when its first 8 KiB are UTF-8 with no NUL byte, and a binary
file is skipped and counted, so a scan that skipped everything cannot look like one that found
nothing. Where there is no git, the tree is walked. This file is read like any other. Three kinds of
pattern, each read the way the harness's `sch conform` reads it (S1), so that the two agree:

  1. The GENERIC shapes of leakage, spelled here because they are no site's: a user home path, an
     e-mail address (an author's, in the citation file or the package metadata, is attribution), a
     dated scratch directory.
  2. The SITE's shapes - its hostnames, its job ids, its queues, its filesystem roots - from the
     file or files named by $SCQC_SITE_SHAPES, one per line, `<regex> :: <what it is> ::
     <literal>[,<literal>]`. A site's naming is the site's fact, so it is supplied from outside
     and never spelled here; a line that is not a shape, or a named file that is missing, FAILS.
  3. TERMS supplied from OUTSIDE the repository - the file or files named by
     $SCQC_FORBIDDEN_TERMS, several joined as PATH joins them (`os.pathsep`) - so that this guard
     never spells the cohort or the site it guards against. A line is matched as its prefix says:
         word:T    T as a whole word, in any case; a number is never found inside a longer number
                   (`4,321` is not in `14,321` or `4,321.5`); spaces in T match any run of spaces;
         exact:T   the same, in T's own case only;
         ...T*     after either prefix, a word that BEGINS with T;
         T         (no prefix) the text anywhere, in any case.
     Until 2026-11-11 this guard read every line as a literal substring, so a generated list's
     `word:` and `exact:` lines - nearly all of it - could match nothing, and the guard printed
     PASS over a tree in which the harness's scan, given the same lists, found hits
     (harness autodev W048).
     Unset, the scan proves less, and says so rather than passing quietly. SET, every file it
     names must exist and yield at least one term, or the guard FAILS: it read one path only, so a
     ':'-joined pair of lists scanned 0 terms and printed PASS - a guard proving nothing in the
     words of one that proved everything.

No file is exempt and no scope list excuses anything here. Stdlib only.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
#: Pruned only when there is no git to say what is the repository's: what a walk finds installed
#: beside the source - an environment, a build, a cache.
SKIP_DIRS = {".git", "__pycache__", ".venv", "venv", "node_modules", "build", "dist", "_data"}
#: The shapes that are no site's. The literal is a speed gate only, as in `sch conform`.
GENERIC = [
    (r"(?<![\w/])/(?:Users|home)/[A-Za-z][\w.-]*", "a user home path", ("/Users/", "/home/")),
    (r"[\w.+-]+@[\w-]+\.(?:edu|com|org|sg|ac\.uk)\b", "an e-mail address", ("@",)),
    (r"scratch/\d{8}__", "a dated scratch directory", ("scratch/",)),
]
SNIFF = 8192
TOKEN = re.compile(r"[a-z0-9]+")


class TermSourceError(ValueError):
    """A term or shape list was named and gave nothing to search for, or cannot be read."""


def _named(var: str) -> list | None:
    spec = os.environ.get(var)
    if spec is None:
        return None
    files = [Path(x.strip()) for x in spec.split(os.pathsep) if x.strip()]
    if not files:
        raise TermSourceError(f"{var} is set and names no file")
    bad = [f"{f}: no such file" for f in files if not f.is_file()]
    if bad:
        raise TermSourceError("; ".join(bad))
    return files


class Term:
    """One line of a word list, read as `sch conform` reads it (its `Term`)."""

    def __init__(self, line: str):
        self.line = line = line.strip()
        mode, text = "substring", line
        for mark, m in (("word:", "word"), ("exact:", "exact")):
            if line.startswith(mark):
                mode, text = m, line[len(mark):].strip()
                break
        prefix = mode != "substring" and len(text) > 1 and text.endswith("*")
        if prefix:
            text = text[:-1].rstrip()
        if mode != "substring" and text == "*":
            text = ""
        self.text, self.mode, self.prefix, self.low = text, mode, prefix, text.lower()
        self.rx = self._regex() if text else None
        self.pieces = tuple(self.low.split()) or (self.low,)
        whole = mode != "substring" and not prefix
        self.tokens = tuple(TOKEN.findall(self.low)) if whole else ()

    def _regex(self):
        text = self.text
        if self.mode == "substring":
            return re.compile(re.escape(text), re.I)
        body = r"\s+".join(re.escape(w) for w in text.split())
        left = right = ""
        if text[:1].isalnum():
            left = r"(?<![A-Za-z0-9])"
        if text[:1].isdigit():
            left += r"(?<![0-9][.,])"
        if not self.prefix and text[-1:].isalnum():
            right = r"(?![A-Za-z0-9])"
        if not self.prefix and text[-1:].isdigit():
            right += r"(?![.,][0-9])"
        return re.compile(left + body + right, re.A | (re.I if self.mode == "word" else 0))

    def admits(self, low: str, tokens: set) -> bool:
        """Can this file hold a match? Never no where it does; a yes is then asked line by line."""
        if self.mode == "substring":
            return self.low in low
        if self.tokens:
            return all(t in tokens for t in self.tokens)
        return all(p in low for p in self.pieces)

    @property
    def label(self) -> str:
        if self.mode == "substring":
            return f"site term {self.text!r}"
        return f"site term {self.text!r} ({self.mode}{', prefix' if self.prefix else ''})"


def terms() -> list:
    """Every term of every list $SCQC_FORBIDDEN_TERMS names; [] only when it is not set.

    A list that is named and missing, or that holds no term, raises: a typo in the path, or a
    value carrying two lists where one was read, must not quietly shrink what the guard looks for.
    """
    files = _named("SCQC_FORBIDDEN_TERMS")
    if files is None:
        return []
    out, bad = [], []
    for f in files:
        got = []
        for n, raw in enumerate(f.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            if not Term(s).text:
                bad.append(f"{f}:{n} names no term: {s!r}")
                continue
            got.append(s)
        if not got:
            bad.append(f"{f}: 0 terms")
        out.extend(got)
    if bad:
        raise TermSourceError("; ".join(bad))
    return [Term(s) for s in dict.fromkeys(out)]


def site_shapes() -> list:
    """(regex, what it is, literals) of every file $SCQC_SITE_SHAPES names; [] when it is not set.
    A line that is not `<regex> :: <what it is>` or not a regular expression raises: a site shape
    that silently does not load is a check that silently does not run."""
    files = _named("SCQC_SITE_SHAPES")
    if files is None:
        return []
    out, bad = [], []
    for f in files:
        for n, raw in enumerate(f.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            parts = [x.strip() for x in s.split("::")]
            if len(parts) < 2 or not parts[0] or not parts[1]:
                bad.append(f"{f}:{n} is not `<regex> :: <what it is>`")
                continue
            try:
                re.compile(parts[0])
            except re.error as e:
                bad.append(f"{f}:{n} is not a regular expression: {e}")
                continue
            lits = tuple(x for x in (parts[2].split(",") if len(parts) > 2 else []) if x)
            out.append((parts[0], parts[1], lits))
    if bad:
        raise TermSourceError("; ".join(bad))
    return out


def _git_files(root: Path):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    try:
        r = subprocess.run(["git", "-C", str(root), "ls-files", "-z", "--cached", "--others",
                            "--exclude-standard"], capture_output=True, timeout=60, env=env)
    except Exception:                                                     # noqa: BLE001
        return None
    if r.returncode != 0:
        return None
    names = sorted({n for n in r.stdout.decode("utf-8", "surrogateescape").split("\0") if n})
    return names or None


def files(root: Path = ROOT) -> list:
    """Every regular file the scan considers: git's list, else a walk; symbolic links are not
    followed (their target is not the repository's text)."""
    names = _git_files(root)
    if names is None:
        names = []
        for d, dirs, fs in os.walk(root, followlinks=False):
            dirs[:] = [x for x in dirs if x not in SKIP_DIRS and not x.endswith(".egg-info")]
            names += [os.path.relpath(os.path.join(d, f), root) for f in fs]
    out = []
    for n in sorted(set(names)):
        p = root / n
        if not p.is_symlink() and p.is_file():
            out.append(p)
    return out


def _binary(p: Path) -> bool:
    try:
        with open(p, "rb") as fh:
            head = fh.read(SNIFF)
    except OSError:
        return True
    if b"\0" in head:
        return True
    try:
        head.decode("utf-8")
    except UnicodeDecodeError as e:
        return not (len(head) == SNIFF and e.start >= SNIFF - 3)
    return False


def scan(root: Path, words: list, shapes: list) -> tuple:
    """(hits, files read, binary files skipped)."""
    compiled = [(re.compile(p), why, lits or ("",)) for p, why, lits in shapes]
    hits, n_read, n_binary = [], 0, 0
    for p in files(root):
        if _binary(p):
            n_binary += 1
            continue
        n_read += 1
        text = p.read_text(encoding="utf-8", errors="replace")
        low = text.lower()
        toks = set(TOKEN.findall(low))
        present = [t for t in words if t.admits(low, toks)]
        gated = [(rx, why) for rx, why, lits in compiled
                 if any(x in text for x in lits) and rx.search(text)]
        if not present and not gated:
            continue
        rel = p.relative_to(root).as_posix()
        # Attribution, not leakage: an address in these names who to cite.
        attribution = p.name in ("CITATION.cff", "pyproject.toml") or rel == ".mailmap"
        for i, line in enumerate(text.splitlines(), 1):
            for rx, why in gated:
                if attribution and why == "an e-mail address":
                    continue
                if rx.search(line):
                    hits.append(f"{rel}:{i}: {why}: {line.strip()[:90]}")
            for t in present:
                if t.rx.search(line):
                    hits.append(f"{rel}:{i}: {t.label}: {line.strip()[:90]}")
    return hits, n_read, n_binary


def main() -> int:
    try:
        site = terms()
        shapes = site_shapes()
    except TermSourceError as e:
        print(f"FAIL: the term source cannot be read as given - {e}. Set SCQC_FORBIDDEN_TERMS to "
              f"one or more files of terms, and SCQC_SITE_SHAPES to the site's shapes, each joined "
              f"by {os.pathsep!r}; a scan of 0 terms from a named source proves nothing and is not "
              f"reported as a pass.")
        return 1
    hits, n_read, n_binary = scan(ROOT, site, GENERIC + shapes)
    print(f"scanned {n_read} files by content ({n_binary} binary skipped); {len(site)} site "
          f"term(s) "
          f"{'from ' + os.environ['SCQC_FORBIDDEN_TERMS'] if site else '(none supplied: set SCQC_FORBIDDEN_TERMS to prove more)'}"
          f"; {len(shapes)} site shape(s) "
          f"{'from ' + os.environ['SCQC_SITE_SHAPES'] if shapes else '(none supplied: set SCQC_SITE_SHAPES to prove more)'}")
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
