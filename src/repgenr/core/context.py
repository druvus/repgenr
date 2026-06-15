"""Per-run execution context.

A :class:`WorkdirContext` bundles the working directory, the canonical path
layout, the logger, the loaded config and (lazily) the SQLite manifest. Stages
receive a context plus their parsed parameters; they never re-derive paths or
re-open the logger themselves.
"""

from __future__ import annotations

import logging
import os
from functools import cached_property
from pathlib import Path

from .config import Config
from .logging import configure_logging
from .manifest import Manifest


class WorkdirContext:
    """Resolved paths and shared services for one RepGenR invocation."""

    def __init__(
        self,
        workdir: str | os.PathLike[str],
        *,
        logger: logging.Logger | None = None,
        create: bool = False,
    ):
        self.workdir = Path(workdir).resolve()
        if create:
            self.workdir.mkdir(parents=True, exist_ok=True)
        self.logger = logger or configure_logging(self.workdir)
        self.config = Config.load(self.workdir)

    # -- canonical path layout (see plan "Canonical data contracts") ----------
    @property
    def genomes_dir(self) -> Path:
        return self.workdir / "genomes"

    @property
    def outgroup_dir(self) -> Path:
        return self.workdir / "outgroup"

    @property
    def derep_dir(self) -> Path:
        return self.workdir / "derep"

    @property
    def representatives_dir(self) -> Path:
        return self.derep_dir / "representatives"

    @property
    def align_dir(self) -> Path:
        return self.workdir / "align"

    @property
    def snp_dir(self) -> Path:
        return self.workdir / "snp"

    @property
    def tree_dir(self) -> Path:
        return self.workdir / "tree"

    @property
    def scratch_dir(self) -> Path:
        return self.workdir / "scratch"

    # -- services -------------------------------------------------------------
    @cached_property
    def manifest(self) -> Manifest:
        return Manifest.open(self.workdir)

    def save_config(self) -> None:
        self.config.save(self.workdir)

    def require_dir(self, path: Path, hint: str) -> Path:
        from .errors import WorkdirError

        if not path.exists():
            raise WorkdirError(f"Expected directory not found: {path}\n{hint}")
        return path
