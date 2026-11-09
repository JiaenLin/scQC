"""Reading decisions.yml, and what makes an entry in it count.

THE FILE IS OPTIONAL. Apply mode, the default, applies the DERIVED and DECLARED thresholds and
records who decided each one; nothing waits for a person. A decisions file is how a person
OVERRIDES one of them, and it is an INPUT, not a command line: it is versioned beside the results,
so two runs with the same data and the same decisions file are the same run, while two runs with
the same command line need not be, because a command line is not kept.

WHAT MAKES AN OVERRIDE COUNT

An ADJUDICATED value needs three things: a `value`, an `approved_by`, and the approver's own words
in `verbatim`. A value missing either attestation is not an override - `engine/steps.py::_attested`
reads it as absent, the derived or declared value is applied instead, and the report says so. The
whole weight of an adjudicated parameter is that a named person looked at the evidence and decided;
a value with no words attached is indistinguishable from a default someone forgot to change.

`check_verbatim()` is the stricter reading of the words - not empty, not a bare acknowledgement
("ok", "yes", "fine"), not too short to be a decision - and `validate()` applies it to every
required entry at once, listing every problem in one pass. Where approval words are given for the
removal itself (`apply.verbatim`), step 7 matches them against the exact action text they were
given for (`action_string()`): words spoken about a different action are not consent to this one.
"""

from __future__ import annotations

from pathlib import Path

#: Words that are agreement with a question rather than a decision about data.
_EMPTY_CONSENT = {
    "ok", "okay", "yes", "y", "sure", "fine", "agreed", "agree", "approved", "confirm",
    "confirmed", "done", "good", "lgtm", "proceed", "go", "accept", "accepted", "n/a", "na",
    "-", "tbd", "todo", "default", "as discussed", "see above",
}

_MIN_VERBATIM_WORDS = 3


class DecisionError(SystemExit):
    """Raised as a hard stop; there is no override path."""


def load(path: Path) -> dict:
    """Parse decisions.yml. PyYAML is imported here so the rest of scQC stays stdlib-only."""
    path = Path(path)
    if not path.exists():
        raise DecisionError(
            f"scqc: no decisions file at {path}.\n"
            f"       --decisions names a file that does not exist. A decisions file is optional:\n"
            f"       without one, apply mode applies the derived and declared thresholds and\n"
            f"       records them as such. Write one only to override a threshold, with its\n"
            f"       value, approved_by and verbatim.")
    try:
        import yaml
    except ImportError:
        raise DecisionError(
            "scqc: reading decisions.yml needs PyYAML.\n"
            "       pip install pyyaml   (or: pip install 'scqc[run]')\n"
            "       Only the run command needs it; the gates and CLI are stdlib-only.") from None
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as e:
        raise DecisionError(f"scqc: {path} is not valid YAML: {e}") from None
    if not isinstance(data, dict):
        raise DecisionError(f"scqc: {path} did not parse to a mapping (got {type(data).__name__})")
    return data


def _blank(v) -> bool:
    return v is None or (isinstance(v, str) and not v.strip())


def check_verbatim(text, where: str) -> list[str]:
    """Is this a decision, or an acknowledgement? Returns problems."""
    if _blank(text):
        return [f"{where}: `verbatim` is empty. An adjudicated value carries the operator's own "
                f"words about THIS decision; without them it cannot be told from a default."]
    s = str(text).strip()
    if s.strip(" .!").lower() in _EMPTY_CONSENT:
        return [f"{where}: `verbatim` is {s!r}, which is agreement with a question rather than a "
                f"decision about the data. Say what was decided and against what evidence."]
    if len(s.split()) < _MIN_VERBATIM_WORDS:
        return [f"{where}: `verbatim` is {s!r} - too short to be a decision. "
                f"At least {_MIN_VERBATIM_WORDS} words."]
    return []


def validate(data: dict, required: dict) -> dict:
    """Check every required adjudicated value is present and properly attested.

    `required` maps a dotted path to a human description, e.g.
    {"quality.umi_floor": "the UMI floor"}. Every one must resolve to a block with a non-blank
    `value`, an `approved_by`, and a `verbatim` that is a decision.

    Returns the resolved {path: value}. Raises DecisionError listing EVERY problem at once -
    reporting them one per run turns a five-minute edit into five rounds.
    """
    problems: list[str] = []
    resolved: dict = {}

    for dotted, what in required.items():
        node = data
        for part in dotted.split("."):
            node = node.get(part) if isinstance(node, dict) else None
            if node is None:
                break
        if node is None:
            problems.append(f"{dotted} ({what}) is missing entirely.")
            continue
        if not isinstance(node, dict):
            # A bare scalar is a value with no attestation. Accept it only for non-adjudicated
            # entries, which never reach this function.
            problems.append(f"{dotted} ({what}) is a bare value. It needs `value`, "
                            f"`approved_by` and `verbatim`.")
            continue
        val = node.get("value")
        if _blank(val):
            problems.append(f"{dotted} ({what}): `value` is empty.")
        if _blank(node.get("approved_by")):
            problems.append(f"{dotted} ({what}): `approved_by` is empty - an approval nobody is "
                            f"named on cannot be questioned later.")
        problems.extend(check_verbatim(node.get("verbatim"), dotted))
        if not _blank(val):
            resolved[dotted] = val

    if problems:
        raise DecisionError(
            "scqc: the decisions file is not complete enough to apply.\n\n  "
            + "\n  ".join(problems)
            + "\n\nEvery adjudicated value needs a number, a person, and that person's own words "
              "about\nthis decision. There is no flag to skip this: it guards the only step that "
              "deletes data.")
    return resolved


def action_string(resolved: dict) -> str:
    """The exact text an approval authorises.

    Built from the resolved values, so changing any threshold changes the string and invalidates
    the approval. That is deliberate: an approval is consent to a specific action, and the same
    words against different numbers are a different decision.
    """
    parts = [f"{k}={resolved[k]}" for k in sorted(resolved)]
    return "scQC filter: " + ", ".join(parts)
