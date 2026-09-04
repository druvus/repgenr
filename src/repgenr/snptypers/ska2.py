"""ska2 (split k-mer analysis) reference-free SNP typer.

Every genome is compared through split k-mers, so no genome is privileged as
the reference and reference-private errors do not bias distances. The output
is a variable-site alignment; there is no positional whole-genome alignment,
so recombination masking is not available with this typer (the snptype stage
refuses ``--mask`` for it with a message naming the alternatives).

Command surface checked against ska 0.5.1 (bioconda ``ska2``): ``ska build -f
<name TAB path list> -k K --threads N -o <prefix>`` writes ``<prefix>.skf``;
``ska align --min-freq F --filter no-ambig-or-const --threads N -o <out>
<prefix>.skf`` writes the alignment.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_int
from .base import SnpParams, SnpResult, SnpTyper


class Ska2Typer(SnpTyper):
    capabilities = ToolCapabilities(
        name="ska2",
        conda=("bioconda::ska2",),
        required_binaries=(BinarySpec("ska", version_args=("--version",), min_version="0.3"),),
        # Split k-mer building is linear in the input; the alignment step reads
        # the whole .skf into memory, which is the practical ceiling.
        recommended_max_genomes=5000,
        default_params={"ksize": 31, "min_freq": 0.9},
        accepted_extras=frozenset({"ksize", "min_freq"}),
    )
    requires_reference = False

    def call(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: SnpParams,
        logger: logging.Logger,
    ) -> SnpResult:
        # ``reference`` is accepted for interface parity and deliberately unused:
        # every genome, including any would-be reference, is an ordinary sample.
        genomes = list(genomes)
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = parse_extra_int(params.extra, "ksize", self.capabilities.default_params["ksize"])
        min_freq = str(params.extra.get("min_freq", self.capabilities.default_params["min_freq"]))

        # ska's input list is two tab-separated columns: sample name, FASTA path.
        # Paths are absolute but not symlink-resolved, matching write_fofn, so
        # the container backend can bind the same directories it sees here.
        listing = out_dir / "genomes.tsv"
        listing.write_text(
            "".join(f"{g.stem}\t{os.path.abspath(g)}\n" for g in genomes), encoding="utf-8"
        )
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})

        prefix = out_dir / "split_kmers"
        run_tool(
            self.capabilities,
            [
                "ska",
                "build",
                "-f",
                listing,
                "-k",
                str(ksize),
                "--threads",
                str(params.threads),
                "-o",
                prefix,
            ],
            logger=logger,
            cwd=out_dir,
            log_prefix="ska-build",
            extra_mounts=genome_dirs,
        )
        skf = Path(str(prefix) + ".skf")
        if not skf.exists():
            raise WorkdirError("ska build did not produce a split k-mer (.skf) file")

        core = out_dir / "core_snp.fasta"
        run_tool(
            self.capabilities,
            [
                "ska",
                "align",
                "--min-freq",
                min_freq,
                "--filter",
                "no-ambig-or-const",
                "--threads",
                str(params.threads),
                "-o",
                core,
                skf,
            ],
            logger=logger,
            cwd=out_dir,
            log_prefix="ska-align",
        )
        if not core.exists() or core.stat().st_size == 0:
            raise WorkdirError(
                "ska align produced no variable-site alignment: the genomes may be "
                "identical, or too divergent for split k-mers at this k (try a smaller "
                "--tool-arg ksize, or lower min_freq)."
            )
        return SnpResult(core_snp_fasta=core, masked=False, full_alignment=None)
