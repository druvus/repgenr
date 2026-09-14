"""sourmash gather against a GTDB sketch database, resolved to a lineage with
``sourmash tax genome`` and the matching lineages CSV."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_chain, run_tool
from ..core.errors import UserInputError
from ..core.executors import parallel_map
from ..core.plugins import ToolCapabilities
from .base import Classification, Classifier, ClassifyParams, db_version


class SourmashClassifier(Classifier):
    capabilities = ToolCapabilities(
        name="sourmash",
        container="quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0",
        conda=("bioconda::sourmash",),
        required_binaries=(BinarySpec("sourmash", version_args=("--version",), min_version="4.0"),),
        default_params={"ksize": 31, "scaled": 1000, "threshold_bp": 50000},
        accepted_extras=frozenset({"ksize", "scaled", "threshold_bp"}),
    )

    def classify(
        self,
        genomes: list[Path],
        out_dir: Path,
        params: ClassifyParams,
        logger: logging.Logger,
    ) -> dict[str, Classification]:
        if params.lineages is None:
            raise UserInputError(
                "The sourmash classifier needs the lineages CSV published with the GTDB "
                "sketch (--gtdb-lineages or REPGENR_GTDB_LINEAGES)."
            )
        defaults = self.capabilities.default_params
        ksize = int(params.extra.get("ksize", defaults["ksize"]))
        scaled = int(params.extra.get("scaled", defaults["scaled"]))
        threshold = int(params.extra.get("threshold_bp", defaults["threshold_bp"]))
        out_dir.mkdir(parents=True, exist_ok=True)
        version = db_version(params.db)
        mounts = [
            str(Path(params.db).resolve().parent),
            str(Path(params.lineages).resolve().parent),
            *sorted({str(g.resolve().parent) for g in genomes}),
        ]

        # One sketch-and-gather chain per genome. A gather is single-threaded
        # and takes tens of seconds against a GTDB-sized sketch, so the chains
        # run side by side within the thread budget; the lineages are then
        # resolved for every gather in one tax call.
        def gather(genome: Path) -> Path:
            work = out_dir / genome.stem
            work.mkdir(parents=True, exist_ok=True)
            sig = work / "query.sig"
            gather_csv = work / "gather.csv"
            steps: list[tuple[str, list[str | Path]]] = [
                (
                    "sourmash",
                    [
                        "sourmash",
                        "sketch",
                        "dna",
                        "-p",
                        f"k={ksize},scaled={scaled}",
                        "--name",
                        genome.stem,
                        "-o",
                        sig,
                        genome,
                    ],
                ),
                (
                    "sourmash",
                    [
                        "sourmash",
                        "gather",
                        sig,
                        params.db,
                        "-k",
                        str(ksize),
                        "--threshold-bp",
                        str(threshold),
                        "-o",
                        gather_csv,
                    ],
                ),
            ]
            run_chain(self.capabilities, steps, logger=logger, extra_mounts=mounts)
            return gather_csv

        workers = max(1, min(params.threads, len(genomes)))
        gathers = parallel_map(gather, genomes, workers, logger=logger)
        # A gather with no match writes only a header; tax genome still lists it.
        base = out_dir / "tax"
        run_tool(
            self.capabilities,
            [
                "sourmash",
                "tax",
                "genome",
                "--gather-csv",
                *gathers,
                "--taxonomy-csv",
                params.lineages,
                "--output-base",
                base,
                "--force",
            ],
            logger=logger,
            log_prefix="sourmash",
            extra_mounts=mounts,
        )
        by_query = _parse_classifications(Path(str(base) + ".classifications.csv"))
        results: dict[str, Classification] = {}
        for genome in genomes:
            hit = by_query.get(genome.stem)
            if hit is not None:
                taxonomy, rank, score = hit
                results[genome.name] = Classification(taxonomy, rank, score, version)
            else:
                logger.warning("%s: no GTDB match above the threshold", genome.name)
        return results


def _parse_classifications(path: Path) -> dict[str, tuple[str, str, float | None]]:
    """Query name -> (lineage, rank, fraction) from a ``tax genome`` table."""
    out: dict[str, tuple[str, str, float | None]] = {}
    if not path.exists():
        return out
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo):
            if rec.get("status") in ("match", "nomatch", "below_threshold") and rec.get("lineage"):
                score = rec.get("fraction")
                out[rec["query_name"]] = (
                    rec["lineage"],
                    rec.get("rank", ""),
                    float(score) if score else None,
                )
    return out
