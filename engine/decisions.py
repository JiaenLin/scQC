"""Reading decisions.yml, and refusing to proceed on one that is incomplete.

A decisions file is an INPUT, not a command line. That is the point of the two-phase design: the
adjudicated values are chosen against the evidence report, written down in the operator's own
words, and versioned beside the results. Two runs with the same data and the same decisions file
are the same run; two runs with the same command line are not, because a command line is not kept.

WHAT THIS REFUSES, AND WHY EACH REFUSAL IS NOT PEDANTRY

  * An ADJUDICATED value with no `verbatim` text. The whole weight of an adjudicated parameter is
    that a person looked at the evidence and decided; a value with no words attached is
    indistinguishable from a default someone forgot to change.
  * A `verbatim` that is only the value ("350") or a bare acknowledgement ("ok", "yes", "fine").
    Those are consent to be asked, not a decision about the data.
  * An adjudicated block whose `approved_by` is absent. An approval nobody is named on cannot be
    questioned later.

None of these can be waived by a flag. There is deliberately no --force: a gate with an override
is a gate that gets overridden, and this one guards the only step that deletes data.
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
            f"       Run `scqc run --mode evidence` first, read the report, then copy\n"
            f"       decisions.template.yml to decisions.yml and record your choices.\n"
            f"       Apply mode has no defaults for the adjudicated values, on purpose.")
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
