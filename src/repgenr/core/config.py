"""Workdir configuration and provenance (``repgenr.yaml``).

A single YAML file at the workdir root records, per stage, the tool chosen, the
parameters used and the resolved tool versions. This replaces the old scattered
``*_parameters.txt`` files and the ``str(dict)`` / ``pickle`` state blobs.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

CONFIG_FILENAME = "repgenr.yaml"
SCHEMA_VERSION = 1


@dataclass
class StageRecord:
    """Provenance for one completed (or in-progress) stage."""

    tool: str | None = None
    params: dict[str, Any] = field(default_factory=dict)
    tool_versions: dict[str, str] = field(default_factory=dict)
    completed: str | None = None  # ISO timestamp, set by caller
    fingerprint: str | None = None  # hash of the stage invocation, for resume

    def to_dict(self) -> dict[str, Any]:
        return {
            "tool": self.tool,
            "params": self.params,
            "tool_versions": self.tool_versions,
            "completed": self.completed,
            "fingerprint": self.fingerprint,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> StageRecord:
        return cls(
            tool=data.get("tool"),
            params=dict(data.get("params") or {}),
            tool_versions=dict(data.get("tool_versions") or {}),
            completed=data.get("completed"),
            fingerprint=data.get("fingerprint"),
        )


@dataclass
class Config:
    """In-memory view of ``repgenr.yaml``."""

    schema_version: int = SCHEMA_VERSION
    repgenr_version: str = ""
    stages: dict[str, StageRecord] = field(default_factory=dict)

    @classmethod
    def load(cls, workdir: str | os.PathLike[str]) -> Config:
        path = Path(workdir) / CONFIG_FILENAME
        if not path.exists():
            from .. import __version__

            return cls(repgenr_version=__version__)
        with open(path) as fo:
            data = yaml.safe_load(fo) or {}
        stages = {
            name: StageRecord.from_dict(rec or {})
            for name, rec in (data.get("stages") or {}).items()
        }
        return cls(
            schema_version=data.get("schema_version", SCHEMA_VERSION),
            repgenr_version=data.get("repgenr_version", ""),
            stages=stages,
        )

    def save(self, workdir: str | os.PathLike[str]) -> Path:
        path = Path(workdir) / CONFIG_FILENAME
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "schema_version": self.schema_version,
            "repgenr_version": self.repgenr_version,
            "stages": {name: rec.to_dict() for name, rec in self.stages.items()},
        }
        with open(path, "w") as fo:
            yaml.safe_dump(data, fo, sort_keys=False, default_flow_style=False)
        return path

    def record_stage(
        self,
        name: str,
        *,
        tool: str | None = None,
        params: dict[str, Any] | None = None,
        tool_versions: dict[str, str] | None = None,
        completed: str | None = None,
    ) -> StageRecord:
        record = StageRecord(
            tool=tool,
            params=dict(params or {}),
            tool_versions=dict(tool_versions or {}),
            completed=completed,
        )
        self.stages[name] = record
        return record
