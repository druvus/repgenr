"""medaka: neural-network consensus polishing of ONT assemblies.

``medaka_consensus`` needs the basecaller model. Reads basecalled with Dorado
carry it in their FASTQ headers (``basecall_model_version_id``) and medaka
resolves it from there; for older reads the model must be given with
``--tool-arg model=...`` or medaka's own default is used. The bacterial
methylation-aware model set (``--bacteria``) is on by default.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from ..assemblers.base import ReadSet
from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import UserInputError
from ..core.plugins import ToolCapabilities
from .base import Polisher, PolishParams, PolishResult, read_dirs

_MODEL_TAG = "basecall_model_version_id="


class MedakaPolisher(Polisher):
    capabilities = ToolCapabilities(
        name="medaka",
        container="quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0",
        conda=("bioconda::medaka",),
        required_binaries=(BinarySpec("medaka", version_args=("--version",), min_version="2.0"),),
        default_params={"bacteria": "true"},
        accepted_extras=frozenset({"model", "bacteria"}),
        ignored_params=frozenset({"rounds"}),
    )
    read_types = frozenset({"OXFORD_NANOPORE"})

    def polish(
        self,
        reads: ReadSet,
        draft: Path,
        out_dir: Path,
        params: PolishParams,
        logger: logging.Logger,
    ) -> PolishResult:
        if len(reads.files) != 1:
            raise UserInputError(
                f"medaka takes one read file; run {reads.run_accession} has {len(reads.files)}."
            )
        out_dir.mkdir(parents=True, exist_ok=True)
        result_dir = out_dir / "medaka_out"
        cmd: list[str | Path] = [
            "medaka_consensus",
            "-i",
            reads.files[0],
            "-d",
            draft,
            "-o",
            result_dir,
            "-t",
            str(params.threads),
        ]
        if params.extra.get("bacteria", self.capabilities.default_params["bacteria"]) != "false":
            cmd.append("--bacteria")
        model = params.extra.get("model")
        if model:
            cmd += ["-m", model]
            source = "--tool-arg"
        elif _header_names_model(reads.files[0]):
            source = "read headers"
        else:
            source = "medaka default"
            logger.warning(
                "%s: the reads name no basecaller model and none was given "
                "(--tool-arg model=...); medaka uses its default model",
                reads.run_accession,
            )
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            log_prefix="medaka",
            cwd=out_dir,
            extra_mounts=read_dirs(reads, draft),
        )
        return PolishResult(
            contigs=result_dir / "consensus.fasta",
            rounds=1,
            tool_stats={"model": model or "", "model_source": source},
        )


def _header_names_model(fastq: Path) -> bool:
    """Whether the first record's header carries a Dorado basecaller model tag."""
    opener = gzip.open if str(fastq).endswith(".gz") else open
    try:
        with opener(fastq, "rt", encoding="utf-8", errors="replace") as fo:
            return _MODEL_TAG in fo.readline()
    except OSError:
        return False
