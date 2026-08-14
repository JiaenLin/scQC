"""What the delivered object says about itself, in `uns["scqc"]`.

A filtered object leaves this pipeline carrying `cluster_FLAG` - step 6's verdict on the cluster
each nucleus sits in. Nothing on the object said what that column MEANT, so a downstream tool had
two bad options: ignore it, and annotate nuclei this pipeline flagged as technical; or guess from
the column name, which is the tool deciding what is technical on scQC's behalf.

This is the third option. The object declares it:

    uns["scqc"] = {
      "schema":       "scqc/provenance@1",
      "flag_column":  "cluster_FLAG",     # "" when step 6 produced nothing for this object
      "flag_meaning": "the cluster this nucleus sits in was FLAGGED by step 6 ...",
      "n_flagged":    3873,
      "n_obs":        109140,
      "flag_digest":  "08e09dc853c538d5",
      ...
    }

A consumer that finds this knows the flag is scQC's, what it means, and - through the digest -
whether the column it is about to act on is still the one scQC wrote. **It remains the consumer's
decision what to do about it.** This module declares; it does not instruct.

THE DIGEST IS A CROSS-TOOL CONTRACT

`flag_digest()` here must agree byte for byte with `scanno/exclude.py::flag_digest()`, or the
verification is worse than useless: it would fail on correct data and teach whoever met it to
pass whatever flag disables the check. The two are deliberately NOT shared as code - one repo
importing the other would make each unusable without the other, and neither depends on the other
today. They are held together instead by a KNOWN-ANSWER VECTOR asserted in both test suites:

    [True, False, True, True, False]  ->  see tests/test_declaration.py

If either implementation drifts, its own suite fails at home rather than the pair failing in
someone's pipeline.

The mask is coerced NA -> False before hashing, because that is what a consumer acting on the
flag does: `cluster_FLAG` is three-valued and "this barcode was never examined" is not "this
barcode was flagged". Hashing the three-valued column instead would make the digest depend on a
choice the consumer has not made yet.
"""
from __future__ import annotations

import hashlib
from datetime import datetime, timezone

#: Bump only for a CHANGE OF MEANING. A consumer keys on this, and a reader five years from now
#: needs to know whether the fields mean what its code assumes.
SCHEMA = "scqc/provenance@1"

#: The uns key. Deliberately not `scqc_apply`, which already exists on the merged object and
#: describes the MERGE - libraries in, nuclei kept - rather than what the flag column means.
KEY = "scqc"

FLAG_COLUMN = "cluster_FLAG"
WATCH_COLUMN = "cluster_WATCH"

#: Written onto the object in full rather than left to documentation, because the object is what
#: a consumer has in its hands and the documentation is what it does not.
FLAG_MEANING = (
    "step 6 flagged the CLUSTER this nucleus sits in, on depth, mitochondrial content, marker "
    "informativeness and doublet fraction. scQC removed nothing on this basis: the flag is "
    "evidence for a decision downstream, not a decision. It cannot distinguish a damaged "
    "population from a cell type this pipeline does not establish the identity of."
)
WATCH_MEANING = (
    "the cluster met the marker-informativeness criterion alone, which is not a QC failure on "
    "its own - a mitochondrial-marker cluster at normal coverage is a mitochondria-rich cell "
    "type, not damage."
)


def flag_digest(flagged) -> str:
    """A short fingerprint of the exact mask. Contract: identical to scAnno's `flag_digest`.

    Hashes the packed bits together with the length, so a mask of a different size cannot
    collide with one of this size. A count cannot do this job: two different masks of the same
    size agree on every number in a summary table.
    """
    import numpy as np

    m = np.asarray(flagged, dtype=bool)
    h = hashlib.sha256()
    h.update(str(m.size).encode())
    h.update(np.packbits(m).tobytes())
    return h.hexdigest()[:16]


def _as_bool(series):
    """NA -> False, which is the coercion a consumer acting on the flag performs.

    Routed through pandas' nullable `boolean` first so the three-valued column becomes two-valued
    exactly once, in the open. `.fillna(False)` straight onto an object column downcasts, which
    pandas deprecates and which would silently change what is hashed the release it stops.
    """
    import numpy as np
    import pandas as pd

    s = pd.Series(series)
    if s.dtype != "boolean":
        s = s.astype("boolean")
    return np.asarray(s.fillna(False).to_numpy(dtype=bool))


def build(adata, *, sample: str = "", run_key: str = "", commit: str = "", version: str = "",
          flag_column: str = FLAG_COLUMN, watch_column: str = WATCH_COLUMN) -> dict:
    """The declaration for one written object. Reads the object; asserts nothing about it.

    Where the object carries no flag column the declaration is still written, with
    `flag_column` empty and `n_flagged` -1. That is deliberate: "step 6 produced nothing for this
    object" and "this object never came from scQC" are different facts, and a consumer that finds
    no declaration at all cannot tell them apart.
    """
    n_obs = int(adata.n_obs)
    has_flag = flag_column and flag_column in adata.obs
    mask = _as_bool(adata.obs[flag_column]) if has_flag else None
    n_examined = (int((~adata.obs[flag_column].isna()).sum()) if has_flag else -1)

    out = {
        "schema": SCHEMA,
        "tool": "scQC",
        "version": str(version or ""),
        "commit": str(commit or ""),
        "run_key": str(run_key or ""),
        "sample": str(sample or ""),
        "n_obs": n_obs,
        "written": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "flag_column": str(flag_column) if has_flag else "",
        "flag_meaning": FLAG_MEANING if has_flag else "",
        "n_flagged": int(mask.sum()) if has_flag else -1,
        # Three-valued: a barcode step 6 never labelled is neither flagged nor unflagged. The
        # count of those it DID examine is carried so a consumer can see the difference rather
        # than reading `n_obs - n_flagged` as "unflagged".
        "n_examined": n_examined,
        "flag_digest": flag_digest(mask) if has_flag else "",
        "flag_na_as": "False",
        "watch_column": (str(watch_column) if watch_column and watch_column in adata.obs else ""),
        "watch_meaning": (WATCH_MEANING if watch_column and watch_column in adata.obs else ""),
        "removed_on_flag": 0,
        "note": ("scQC removed nothing on the flag. What to do about it is the consumer's "
                 "decision; this records what the flag is and which nuclei carry it."),
    }
    return out


def stamp(adata, **kw) -> dict:
    """Build the declaration and write it into `uns`. Returns it."""
    decl = build(adata, **kw)
    adata.uns[KEY] = decl
    return decl


def verify(adata, decl=None) -> tuple[bool, str]:
    """Does the object still carry the flag the declaration describes?

    Returns `(ok, reason)`. A consumer should refuse rather than proceed on a mismatch: a flag
    column that has been rewritten since scQC wrote it is not scQC's decision any more, and
    acting on it while citing scQC's provenance would attribute someone else's choice to this
    pipeline.
    """
    decl = decl if decl is not None else adata.uns.get(KEY)
    if not decl:
        return False, f"the object carries no uns[{KEY!r}] declaration"
    if str(decl.get("schema", "")) != SCHEMA:
        return False, (f"declaration schema is {decl.get('schema')!r}, this reader understands "
                       f"{SCHEMA!r}")
    col = str(decl.get("flag_column") or "")
    if not col:
        return True, "the declaration records no flag column, and the object carries none"
    if col not in adata.obs:
        return False, f"the declaration names {col!r} and the object has no such column"
    if int(decl.get("n_obs", -1)) != int(adata.n_obs):
        return False, (f"the declaration is for {decl.get('n_obs'):,} observations and this "
                       f"object holds {adata.n_obs:,}; it has been subset since scQC wrote it")
    got = flag_digest(_as_bool(adata.obs[col]))
    want = str(decl.get("flag_digest") or "")
    if got != want:
        return False, (f"{col!r} does not match the declaration: digest {got} against {want}. "
                       f"The column has been altered since scQC wrote it, so acting on it is not "
                       f"acting on scQC's decision")
    return True, f"{col!r} matches the declaration ({decl.get('n_flagged'):,} flagged)"
