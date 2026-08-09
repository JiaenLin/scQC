# Contributing

Issues and pull requests are welcome. The expectations below are not style preferences; each one
exists because its absence produces a defect that a passing test suite cannot see.

## Before you open a pull request

```bash
./bin/scqc selftest        # ten unit suites + the adversarial suite; non-zero exit if any failed
```

Every suite but one runs on the standard library alone; `tests/test_audit_ambient.py` needs pandas
and a `COHORT_DIR`, and reports **SKIP** without them. Read a SKIP as what it is: that suite has
not been checked, which is not the same as passing.

Do not weaken a test to make it pass. If a test fails because the code is wrong, the code is what
changes.

## 1. Anything that removes observations goes through step 7

Steps 0–6 measure, score, flag and refuse. **Step 7 is the only step permitted to drop an
observation**, it requires a recorded approval carrying the operator's own words for that specific
action, and there is no force flag. A pull request that adds a removal anywhere else will be
declined regardless of how well it is tested.

A pull request that adds or changes a removal answers the removal checklist from
[docs/PRINCIPLES.md](docs/PRINCIPLES.md) in its description:

1. What exactly does it remove — the actual list, not the category?
2. Does the removed signal exist anywhere else in the data?
3. Is the removal differential across the design? Compute it per level.
4. Is there a reversible version, and why is it not the one being proposed?
5. If it must be irreversible, is what was removed still recoverable?

## 2. A new gate needs a materiality bound

A gate that fires on correct behaviour gets switched off, which is worse than not having it. A
ratio between two near-zero rates is arithmetic, not evidence: a 5.42× differential computed from
541 observations against 76 is not a finding. State the bound, state the value you chose, and say
what data the value was chosen against.

The same applies in the other direction. A gate that has never fired on any data available is a
hypothesis about the threshold, not a tested one — say so in the module rather than letting the
number read as calibrated.

## 3. A parameter carries a class

Every parameter is FIXED, DERIVED, DECLARED or ADJUDICATED — see the table in
[README.md](README.md). A new parameter without a class, or an ADJUDICATED parameter with a
default, is a defect. DERIVED means the procedure is fixed and the value is computed per dataset;
if you find yourself adding a command-line argument for a value the data can produce, it belongs in
a derivation instead.

## 4. Unknown is not a value

Never let a missing measurement read as zero, false or passing. `or 0`, `or 1e-12` and a bare
comparison against `None` are the three shapes this bug takes. Never-examined is its own category
in the tables, in the gates and in the figures.

## 5. Claims about the pipeline are checked against the tree

Documentation that describes behaviour the code does not have is the same defect class as a stale
report: it reads exactly like a correct document. If you add a claim to README.md or to `docs/`,
either point at the code that implements it or mark it as specification. The Status table in
README.md is the honest inventory and is expected to be updated by the pull request that changes
what is true.

## Style

- **Explain why a rule exists.** A threshold with no recorded justification is a bug, not a choice.
- **Conclusions belong in Markdown, not in code comments.** Comments explain the code.
- **LF line endings.** Enforced by `.gitattributes`; these files run on Linux clusters.
- **Module docstrings state what the module does not do.** A reviewer needs to establish quickly
  that a step cannot silently cost them observations.
