"""Waiting for a file that already exists, on a filesystem that has not admitted it yet.

THIS IS NOT A RETRY LOOP AROUND A FLAKY WRITE.

Every task here may run on a different machine from the orchestrator, and the project lives on
NFS. An NFS client caches directory attributes - 30 to 60 seconds by default, and nothing in the
mount options on this cluster changes that - so a file CREATED on one node is genuinely invisible
to `stat()` on another until the cache expires. The write finished, the data is on the server, and
the reader is looking at a directory it last read before the file was there.

What that produces is the worst-shaped failure this pipeline has: a task that ran correctly, exited
0 and wrote its output is reported as having written nothing. It is indistinguishable from a real
failure, and it only happens when the file did not ALREADY exist from an earlier attempt - so it
survived every resumed run and appeared on the first genuinely clean one, where two of ten
libraries failed with their 22,284 and 23,584 line outputs sitting correctly on disk.

Listing the parent is what gives the client the chance to notice; a bare `stat()` of the same
name can keep answering from a cached negative lookup. Neither forces anything before the
attribute cache expires, which is why this waits rather than merely re-checking.

A file that never appears still fails, and says how long it was waited for, because "not yet
visible" and "never written" must not end up reading the same way.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

#: Comfortably past a Linux NFS client's `acdirmax`, which defaults to 60 s and is not overridden
#: on this cluster. Not a guess at how long a tool takes: by the time this is called the tool has
#: already exited, so the only thing being waited for is the filesystem catching up.
VISIBILITY_TIMEOUT_S = 90


def await_visible(paths, timeout_s: float = VISIBILITY_TIMEOUT_S) -> list:
    """Wait for every path to become visible. Returns those that never did, as strings.

    An empty list means all of them arrived. Returns immediately when they are already there,
    which is the usual case - nothing pays for this except a run that would otherwise have
    reported a false failure.
    """
    remaining = [Path(p) for p in paths]
    remaining = [p for p in remaining if not p.exists()]
    if not remaining:
        return []
    deadline = time.time() + max(0.0, float(timeout_s))
    wait = 0.5
    while True:
        for parent in {p.parent for p in remaining}:
            try:
                os.listdir(parent)
            except OSError:
                pass
        remaining = [p for p in remaining if not p.exists()]
        if not remaining or time.time() >= deadline:
            return [str(p) for p in remaining]
        time.sleep(wait)
        wait = min(wait * 2, 5.0)
