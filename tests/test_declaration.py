"""What the delivered object says about itself, and the digest contract that makes it checkable.

`uns["scqc"]` exists so a downstream tool can act on `cluster_FLAG` knowing whose decision it is
rather than guessing from a column name. Two things have to hold for that to be worth anything:

  1. the declaration must describe the object it is stamped on, not the one upstream meant to
     write - so the digest is computed from the column as it stands at write time;
  2. `flag_digest()` must agree byte for byte with `scanno/exclude.py::flag_digest()`, or
     verification fails on correct data and teaches whoever meets it to disable the check.

The second is a contract between two repositories that do not import each other. It is held by
the KNOWN-ANSWER VECTOR below, asserted in both test suites: if either implementation drifts,
its own suite fails at home rather than the pair failing in someone's pipeline.

    python tests/test_declaration.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

fails = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        fails.append(name)


try:
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"SKIP: needs {e.name}")
    raise SystemExit(0)

from adapters import declaration as decl  # noqa: E402

# --------------------------------------------------------------------------------------------
# THE CONTRACT. Do not "fix" this number to match the code; if they disagree, one of the two
# implementations changed and the question is which. The same vector is asserted in scAnno's
# tests/test_emit.py against scanno.exclude.flag_digest.
KNOWN_MASK = [True, False, True, True, False]
KNOWN_DIGEST = "3ba679de109f5333"

print("\n1 - the digest is a function of the bits and the length, nothing else")
check("the CONTRACT vector still hashes to the agreed value",
      decl.flag_digest(np.array(KNOWN_MASK, dtype=bool)) == KNOWN_DIGEST,
      f"got {decl.flag_digest(np.array(KNOWN_MASK, dtype=bool))}, agreed {KNOWN_DIGEST} - if "
      f"this fails, THIS implementation changed; scAnno asserts the same vector")
check("stable across calls", decl.flag_digest(np.array(KNOWN_MASK)) == KNOWN_DIGEST)
check("16 hex characters", len(KNOWN_DIGEST) == 16 and all(c in "0123456789abcdef"
                                                           for c in KNOWN_DIGEST), KNOWN_DIGEST)
check("a list and an array agree", decl.flag_digest(KNOWN_MASK) == KNOWN_DIGEST)
check("a different mask of the SAME size differs",
      decl.flag_digest([True, True, True, True, False]) != KNOWN_DIGEST)
# Length is hashed separately for this: packbits pads to a byte, so masks of different lengths
# can share packed bytes. Without the length these two would collide.
check("padding cannot collide two different lengths",
      decl.flag_digest([True] * 3) != decl.flag_digest([True] * 3 + [False] * 5))
print(f"        contract vector {KNOWN_MASK} -> {KNOWN_DIGEST}")

print("\n2 - NA is hashed as False, which is what a consumer acting on the flag does")
three = pd.Series([True, None, False], dtype="object")
check("NA coerces to False", decl.flag_digest(decl._as_bool(three)) ==
      decl.flag_digest([True, False, False]))


def toy(n=40, flag=True, watch=True, na=0):
    """As step 7 writes them: `cluster_FLAG` is a pandas NULLABLE boolean, unknown is pd.NA.

    The dtype is part of the fixture rather than incidental. `adapters/scanpy_ops.py`
    ::_annotation_column types the flags that way precisely so "this cluster was not flagged"
    and "this barcode was never examined" survive into the object as different facts - and an
    object-dtype column of Python bools, which is the obvious thing to write in a test, is not
    writable by anndata at all. A fixture that cannot be written cannot test a round trip.
    """
    import anndata as ad
    import scipy.sparse as sp
    rng = np.random.default_rng(0)
    A = ad.AnnData(X=sp.csr_matrix(rng.random((n, 6)).astype("float32")))
    A.obs_names = [f"b{i}" for i in range(n)]
    if flag:
        col = [bool(i % 5 == 0) for i in range(n)]
        for i in range(na):
            col[-(i + 1)] = pd.NA
        A.obs["cluster_FLAG"] = pd.array(col, dtype="boolean")
    if watch:
        A.obs["cluster_WATCH"] = pd.array([bool(i % 7 == 0) for i in range(n)], dtype="boolean")
    return A


try:
    import anndata  # noqa: F401
    import scipy.sparse  # noqa: F401
except ImportError as e:
    print(f"\nSKIP the object checks: needs {e.name}")
    print("\n" + "=" * 64)
    if fails:
        print(f"declaration: {len(fails)} FAILED - " + ", ".join(fails))
        raise SystemExit(1)
    print("declaration OK (object checks skipped - a skip is not a pass)")
    raise SystemExit(0)

print("\n3 - the declaration describes the object it is stamped on")
A = toy()
d = decl.stamp(A, sample="lib1", run_key="abc123", commit="deadbeef", version="0.3.0")
check("it lands in uns under the documented key", decl.KEY in A.uns)
check("schema is declared", d["schema"] == decl.SCHEMA, d["schema"])
check("the flag column is named", d["flag_column"] == "cluster_FLAG")
check("n_flagged counts the flag", d["n_flagged"] == int(np.sum([i % 5 == 0 for i in range(40)])),
      str(d["n_flagged"]))
check("n_obs matches", d["n_obs"] == A.n_obs)
check("the digest is of THIS object's column",
      d["flag_digest"] == decl.flag_digest(decl._as_bool(A.obs["cluster_FLAG"])))
check("run identity is carried", (d["run_key"], d["commit"]) == ("abc123", "deadbeef"))
check("what the flag MEANS travels with it, not just its name",
      "removed nothing" in d["flag_meaning"], d["flag_meaning"][:60])
check("and it records that scQC removed nothing on it", d["removed_on_flag"] == 0)

print("\n4 - verify() accepts the object it was written for")
ok, why = decl.verify(A)
check("verified", ok, why)

print("\n5 - verify() REFUSES a flag that has been altered since")
B = toy()
decl.stamp(B, sample="lib1")
B.obs["cluster_FLAG"] = pd.array([True] * B.n_obs, dtype="boolean")
ok, why = decl.verify(B)
check("altered flag is refused", not ok, why)
check("...and the reason names the digest mismatch", "digest" in why)

print("\n6 - verify() REFUSES an object that has been subset since")
C = toy()
decl.stamp(C, sample="lib1")
ok, why = decl.verify(C[:10].copy())
check("subset object is refused", not ok, why)
check("...and says so in those terms", "subset" in why, why)

print("\n7 - an object with no flag is declared as HAVING no flag, not left silent")
D = toy(flag=False, watch=False)
d2 = decl.stamp(D, sample="lib2")
check("flag_column is empty", d2["flag_column"] == "")
check("n_flagged is -1, not 0 - 'none produced' is not 'none flagged'", d2["n_flagged"] == -1)
check("the declaration is still written", decl.KEY in D.uns)
ok, why = decl.verify(D)
check("and it verifies", ok, why)

print("\n8 - three-valued flags: examined is counted apart from flagged")
E = toy(na=6)
d3 = decl.stamp(E, sample="lib3")
check("n_examined excludes the NA rows", d3["n_examined"] == E.n_obs - 6, str(d3["n_examined"]))
check("n_obs - n_flagged is NOT the unflagged count, and the record lets a reader see that",
      d3["n_examined"] < d3["n_obs"])
check("the coercion is declared", d3["flag_na_as"] == "False")

print("\n9 - no declaration at all is refused rather than treated as 'no flag'")
F = toy()
ok, why = decl.verify(F)
check("undeclared object is refused", not ok, why)

print("\n10 - it survives a round trip through h5ad, which is the point")
import tempfile  # noqa: E402

with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp) / "one.h5ad"
    G = toy()
    decl.stamp(G, sample="lib1", run_key="k", commit="c", version="0.3.0")
    G.write_h5ad(p)
    back = anndata.read_h5ad(p)
    check("the declaration is readable off disk", decl.KEY in back.uns)
    ok, why = decl.verify(back)
    check("and still verifies against the object it came with", ok, why)
    check("the digest survived unchanged",
          str(back.uns[decl.KEY]["flag_digest"]) == G.uns[decl.KEY]["flag_digest"])

print("\n" + "=" * 64)
if fails:
    print(f"declaration: {len(fails)} FAILED - " + ", ".join(fails))
    raise SystemExit(1)
print("declaration OK - the object says what its flag means, and the digest proves it is intact")
