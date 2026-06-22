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

from pathlib import Path


def write_versions_fragment(path: str | Path, versions: dict[str, str]) -> None:
    """Write ``versions`` as 4-space-indented ``tool: version`` lines.

    An empty mapping writes an empty file (the module still records repgenr).
    """
    lines = [f"    {tool}: {ver}" for tool, ver in sorted(versions.items())]
    text = "\n".join(lines) + "\n" if lines else ""
    Path(path).write_text(text)
