"""Write resolved tool versions as a Nextflow ``versions.yml`` fragment.

The data-channel steps run one tool family each and know the versions they
resolved (via ``check_binaries``/``preflight``). They emit those as YAML lines
indented to slot under a process key, so the calling Nextflow module can compose

    "PROCESS":
        repgenr: 2.0.0
        skder: 1.3.0

and the underlying bioinformatics tools -- not just repgenr -- end up in
provenance even though the data-channel path keeps no shared ``repgenr.yaml``.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from pathlib import Path

# A value a YAML 1.1 loader would not keep as a string: a date or timestamp
# (the GTDB API query date) or one holding ': '. Such values are written
# double-quoted; tool versions such as 1.4.6 stay plain.
_NEEDS_QUOTES = re.compile(r"^\d{4}-\d{2}-\d{2}|: |^[\s'\"{\[&*!|>%@`#]")


def _scalar(value: str) -> str:
    return json.dumps(value) if _NEEDS_QUOTES.search(value) else value


def merge_stage_versions(stage_versions: Mapping[str, Mapping[str, str]]) -> dict[str, str]:
    """Combine per-stage ``tool -> version`` maps into one mapping.

    A tool recorded with one version by every stage that ran it keeps its
    plain name. A tool recorded with different versions (a pinned image in one
    stage, the host binary in another) gets one ``"<tool> (<stage>)"`` entry
    per stage, so no value is lost and the keys stay unique.
    """
    by_tool: dict[str, dict[str, str]] = {}
    for stage, versions in stage_versions.items():
        for tool, version in versions.items():
            by_tool.setdefault(tool, {})[stage] = version
    merged: dict[str, str] = {}
    for tool, per_stage in by_tool.items():
        if len(set(per_stage.values())) == 1:
            merged[tool] = next(iter(per_stage.values()))
        else:
            for stage, version in per_stage.items():
                merged[f"{tool} ({stage})"] = version
    return merged


def write_versions_fragment(path: str | Path, versions: dict[str, str]) -> None:
    """Write ``versions`` as 4-space-indented ``tool: version`` lines.

    An empty mapping writes an empty file (the module still records repgenr).
    """
    lines = [f"    {tool}: {_scalar(str(ver))}" for tool, ver in sorted(versions.items())]
    text = "\n".join(lines) + "\n" if lines else ""
    # The steps write this before their output directory exists (a Nextflow
    # task writes it into the task cwd, a manual run often into -o).
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text, encoding="utf-8")
