# Runs conf/env/build_reference.sh under bash with a stand-in celescope, curl and wget in a temporary
# directory; downloads nothing, builds nothing, and needs only the standard library.
"""The reference build names the reference it builds: no species, assembly or release is assumed.

`conf/env/build_reference.sh` defaulted SPECIES, ASSEMBLY and ENSEMBL_RELEASE to one calibration cohort's
reference - an organism, its assembly and an Ensembl release - so a project that set none of them was given
that cohort's genome without a word (harness gap C31). A default equal to one cohort's value is the overfit a
leak guard cannot see: it spells nothing a word list would hold. The three are required now, with no default,
and the refusal names all three.

  1  none set: refused before celescope is looked for, a directory is made or anything is downloaded, and
     the message names SPECIES, ASSEMBLY and ENSEMBL_RELEASE
  2  one or two set: refused the same way, naming the ones missing
  3  a release that is not a whole number: refused the same way
  4  all three set: not refused, and the download asked for is built from exactly the three values given -
     so the script is not refusing everything (the stand-in download fails, which is where it stops)
  5  `--help` needs none of them
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "conf" / "env" / "build_reference.sh"
NAMES = ("SPECIES", "ASSEMBLY", "ENSEMBL_RELEASE")
# Values no cohort has: the test must not teach the script, or this repository, a real reference.
GIVEN = {"SPECIES": "testus_exemplaris", "ASSEMBLY": "TST9", "ENSEMBL_RELEASE": "7"}
fails: list[str] = []


def run(values: dict, *extra: str):
    """Run the script with exactly `values` in its environment; return (exit, output, calls, out_made)."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        stand = td / "bin"
        stand.mkdir()
        log = td / "calls.log"
        for tool in ("curl", "wget", "celescope"):
            p = stand / tool
            # A stand-in that records its call; the downloaders fail, as a download that found nothing does.
            p.write_text(f'#!/bin/sh\necho "{tool} $*" >> "{log}"\n'
                         f'exit {0 if tool == "celescope" else 22}\n', encoding="utf-8")
            p.chmod(0o755)
        out = td / "reference"
        env = {"PATH": os.pathsep.join([str(stand), "/usr/bin", "/bin"]), "HOME": str(td), "LC_ALL": "C"}
        env.update(values)
        argv = ["bash", str(SCRIPT)] + (list(extra) if extra else
                                        ["--out", str(out), "--celescope", str(stand / "celescope")])
        r = subprocess.run(argv, capture_output=True, text=True, env=env, cwd=str(td), timeout=120)
        calls = log.read_text(encoding="utf-8") if log.exists() else ""
        return r.returncode, r.stdout + r.stderr, calls, out.exists()


def refused(case: str, values: dict, missing: tuple):
    rc, said, calls, made = run(values)
    named = [n for n in NAMES if n not in said]
    if rc == 0 or named or any(m not in said for m in missing) or calls or made:
        fails.append(f"{case}: exit {rc}; names absent from the message {named}; calls made {calls.split()[:2]}; "
                     f"--out created {made}: {said.strip()[-240:]!r}")
    else:
        print(f"  {case}: refused (exit {rc}), naming all three; nothing looked for, made or downloaded")


print("conf/env/build_reference.sh - the reference is named by the caller, never assumed")
print("=" * 74)

refused("1 none set", {}, NAMES)
refused("2a SPECIES only", {"SPECIES": GIVEN["SPECIES"]}, ("ASSEMBLY", "ENSEMBL_RELEASE"))
refused("2b SPECIES and ASSEMBLY", {k: GIVEN[k] for k in ("SPECIES", "ASSEMBLY")}, ("ENSEMBL_RELEASE",))
refused("2c ENSEMBL_RELEASE only", {"ENSEMBL_RELEASE": GIVEN["ENSEMBL_RELEASE"]}, ("SPECIES", "ASSEMBLY"))
refused("3 a release that is not a whole number", dict(GIVEN, ENSEMBL_RELEASE="latest"), ("latest",))

rc, said, calls, made = run(GIVEN)
want = ("release-7/", "testus_exemplaris/", "Testus_exemplaris.TST9.7.gtf.gz")
asked = [ln for ln in calls.splitlines() if ln.startswith(("curl", "wget"))]
if not asked or not all(w in asked[0] for w in want):
    fails.append(f"4 all three set: the download asked for was {asked[:1]}, not one built from "
                 f"{GIVEN} (exit {rc}): {said.strip()[-240:]!r}")
elif any(f"Not set: {n}" in said for n in NAMES):
    fails.append(f"4 all three set: still refused for want of them: {said.strip()[-200:]!r}")
else:
    print(f"  4 all three set: not refused; the download asked for is built from them: {asked[0][:110]}")

rc, said, calls, made = run({}, "--help")
if rc != 0 or not all(n in said for n in NAMES):
    fails.append(f"5 --help with nothing set: exit {rc}, or the three are not described: {said[-200:]!r}")
else:
    print("  5 --help needs none of them, and describes all three")

print("=" * 74)
if fails:
    print(f"FAILED - {len(fails)}:")
    for f in fails:
        print("  -", f)
    raise SystemExit(1)
print("the reference build names its reference; nothing is assumed")
