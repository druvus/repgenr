"""medaka: neural-network consensus polishing of ONT assemblies.

``medaka_consensus`` needs the basecaller model. Reads basecalled with Dorado
carry it in their FASTQ headers (``basecall_model_version_id``) and medaka
resolves it from there, choosing the bacterial variant with ``--bacteria``.
Reads mirrored through SRA have their headers rewritten and carry no such
tag, so for them the model is ``--tool-arg model=...`` or, by default, ONT's
bacterial methylation-aware model for R10.4.1 400 bps chemistry, the common
case for public bacterial ONT runs since 2023; the marker records which of
the three applied. ``--tool-arg bacteria=false`` uses medaka's default
model instead.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from ..assemblers.base import ReadSet
from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import ToolExecutionError
from ..core.plugins import ToolCapabilities
from .base import Polisher, PolishParams, PolishResult, one_read_file, read_dirs

_MODEL_TAG = "basecall_model_version_id="
# ONT's bacterial methylation-aware model (R10.4.1, 400 bps), assumed when the
# reads name no basecaller and no model is given.
BACTERIAL_MODEL = "r1041_e82_400bps_bacterial_methylation"
_INCOMPATIBLE = "--bacteria was specified but input data was not compatible"


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
        out_dir.mkdir(parents=True, exist_ok=True)
        reads = one_read_file(reads, out_dir)
        result_dir = out_dir / "medaka_out"
        bacteria = (
            params.extra.get("bacteria", self.capabilities.default_params["bacteria"]) != "false"
        )
        model = params.extra.get("model") or ""
        auto = False
        if model:
            source = "--tool-arg"
        elif _header_names_model(reads.files[0]):
            source, auto = "read headers", True
        elif bacteria:
            model, source = BACTERIAL_MODEL, "assumed bacterial R10.4.1 model"
            logger.warning(
                "%s: the reads name no basecaller model (SRA rewrites FASTQ headers); "
                "using %s. Give the basecaller's model with --tool-arg model=... if known.",
                reads.run_accession,
                model,
            )
        else:
            source = "medaka default"
            logger.warning(
                "%s: the reads name no basecaller model and none was given "
                "(--tool-arg model=...); medaka uses its default model",
                reads.run_accession,
            )

        def command(with_bacteria: bool) -> list[str | Path]:
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
            if model:
                cmd += ["-m", model]
            if with_bacteria:
                cmd.append("--bacteria")
            return cmd

        # With a header-named basecaller, --bacteria lets medaka pick the
        # bacterial variant; it refuses when none exists for that basecaller,
        # in which case the plain auto-selected model is used.
        with_bacteria = auto and bacteria
        try:
            self._run(command(with_bacteria), reads, draft, out_dir, logger)
        except ToolExecutionError as exc:
            if not (with_bacteria and exc.output and _INCOMPATIBLE in exc.output):
                raise
            logger.warning(
                "%s: medaka has no bacterial model for these reads' basecaller; "
                "polishing with the general model",
                reads.run_accession,
            )
            with_bacteria = False
            self._run(command(False), reads, draft, out_dir, logger)
        return PolishResult(
            contigs=result_dir / "consensus.fasta",
            rounds=1,
            tool_stats={
                "model": model,
                "model_source": source,
                "bacteria": with_bacteria or bool(model == BACTERIAL_MODEL),
            },
        )

    def _run(
        self,
        cmd: list[str | Path],
        reads: ReadSet,
        draft: Path,
        out_dir: Path,
        logger: logging.Logger,
    ) -> None:
        run_tool(
            self.capabilities,
            cmd,
            logger=logger,
            log_prefix="medaka",
            cwd=out_dir,
            extra_mounts=read_dirs(reads, draft),
        )


def _header_names_model(fastq: Path) -> bool:
    """Whether the first record's header carries a Dorado basecaller model tag."""
    opener = gzip.open if str(fastq).endswith(".gz") else open
    try:
        with opener(fastq, "rt", encoding="utf-8", errors="replace") as fo:
            return _MODEL_TAG in fo.readline()
    except OSError:
        return False
