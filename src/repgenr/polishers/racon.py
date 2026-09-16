"""racon: consensus polishing from read-to-draft overlaps (minimap2), for ONT
and PacBio CLR assemblies; each round maps the reads to the previous consensus."""

from __future__ import annotations

import logging
from pathlib import Path

from ..assemblers.base import ReadSet
from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.plugins import ToolCapabilities
from .base import Polisher, PolishParams, PolishResult, read_dirs

# HiFi reads are already accurate; polishing them with racon does more harm than good.
_HIFI_INSTRUMENTS = ("sequel ii", "revio")

# minimap2 runs in its own pinned image; the racon BioContainer holds racon only.
_MINIMAP2 = ToolCapabilities(
    name="minimap2",
    container="quay.io/biocontainers/minimap2:2.28--h577a1d6_4",
    conda=("bioconda::minimap2",),
    required_binaries=(BinarySpec("minimap2", version_args=("--version",), min_version="2.17"),),
)


class RaconPolisher(Polisher):
    capabilities = ToolCapabilities(
        name="racon",
        container="quay.io/biocontainers/racon:1.5.0--hdcf5f25_5",
        conda=("bioconda::racon", "bioconda::minimap2"),
        required_binaries=(
            BinarySpec("racon", version_args=("--version",), min_version="1.4"),
            *_MINIMAP2.required_binaries,
        ),
    )
    read_types = frozenset({"OXFORD_NANOPORE", "PACBIO_SMRT"})

    def accepts(self, reads: ReadSet) -> bool:
        model = reads.instrument_model.lower()
        return reads.platform in self.read_types and not any(
            tag in model for tag in _HIFI_INSTRUMENTS
        )

    def polish(
        self,
        reads: ReadSet,
        draft: Path,
        out_dir: Path,
        params: PolishParams,
        logger: logging.Logger,
    ) -> PolishResult:
        out_dir.mkdir(parents=True, exist_ok=True)
        preset = "map-ont" if reads.platform == "OXFORD_NANOPORE" else "map-pb"
        query = reads.files[0]
        current = draft
        rounds = max(1, params.rounds)
        for i in range(1, rounds + 1):
            overlaps = out_dir / f"round{i}.paf"
            consensus = out_dir / f"round{i}.fasta"
            run_tool(
                _MINIMAP2,
                [
                    "minimap2",
                    "-x",
                    preset,
                    "-t",
                    str(params.threads),
                    "-o",
                    overlaps,
                    current,
                    query,
                ],
                logger=logger,
                log_prefix="minimap2",
                cwd=out_dir,
                extra_mounts=read_dirs(reads, current),
            )
            run_tool(
                self.capabilities,
                ["racon", "-t", str(params.threads), query, overlaps, current],
                logger=logger,
                log_prefix="racon",
                cwd=out_dir,
                stdout_path=consensus,
                extra_mounts=read_dirs(reads, current),
            )
            current = consensus
        return PolishResult(contigs=current, rounds=rounds, tool_stats={"preset": preset})
