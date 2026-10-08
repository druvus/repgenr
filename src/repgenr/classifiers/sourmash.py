"""sourmash gather against a GTDB sketch database, resolved to a lineage with
``sourmash tax genome`` and the matching lineages CSV.

A genome given in ``ClassifyParams.sketches`` is not sketched here: gather
reads its ``.sig.zip`` and selects the signature at the requested k-mer size
with ``-k``. Each query signature is named by the genome's record name.
"""

from __future__ import annotations

import csv
import logging
from collections.abc import Mapping
from pathlib import Path

from ..core.containers import run_chain, run_tool
from ..core.contracts import record_name
from ..core.errors import UserInputError
from ..core.executors import parallel_map
from ..core.plugins import parse_extra_int
from ..core.sourmash import sourmash_capabilities
from .base import Classification, Classifier, ClassifyParams, db_version

# Peak resident memory of one gather against the GTDB rs226 representatives
# sketch (k=31): about 0.6 GB, measured natively and in the container.
GATHER_GB = 0.6


def gather_workers(n_genomes: int, threads: int, memory_gb: float | None) -> int:
    """Concurrent gathers: bounded by the genomes, the threads and the memory budget."""
    workers = max(1, min(threads, n_genomes))
    if memory_gb is not None:
        workers = min(workers, max(1, int(memory_gb / GATHER_GB + 1e-9)))
    return workers


class SourmashClassifier(Classifier):
    capabilities = sourmash_capabilities(
        default_params={"ksize": 31, "scaled": 1000, "threshold_bp": 50000},
        accepted_extras=frozenset({"ksize", "scaled", "threshold_bp"}),
    )
    needs_lineages = True

    def sketch_request(self, extra: Mapping[str, object]) -> tuple[int, int]:
        defaults = self.capabilities.default_params
        return (
            parse_extra_int(extra, "ksize", defaults["ksize"]),
            parse_extra_int(extra, "scaled", defaults["scaled"]),
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
        given = {
            g: Path(params.sketches[g]) for g in genomes if params.sketches and g in params.sketches
        }
        if given:
            mounts += sorted({str(p.resolve().parent) for p in given.values()})
            logger.info(
                "sourmash: %d of %d queries read from the genome sketches", len(given), len(genomes)
            )

        def gather(genome: Path) -> Path:
            name = record_name(genome)
            work = out_dir / name
            work.mkdir(parents=True, exist_ok=True)
            gather_csv = work / "gather.csv"
            steps: list[tuple[str, list[str | Path]]] = []
            sig = given.get(genome)
            if sig is None:
                sig = work / "query.sig"
                steps.append(
                    (
                        "sourmash",
                        [
                            "sourmash",
                            "sketch",
                            "dna",
                            "-p",
                            f"k={ksize},scaled={scaled}",
                            "--name",
                            name,
                            "-o",
                            sig,
                            genome,
                        ],
                    )
                )
            steps.append(
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
                )
            )
            run_chain(self.capabilities, steps, logger=logger, extra_mounts=mounts)
            return gather_csv

        workers = gather_workers(len(genomes), params.threads, params.memory_gb)
        logger.info(
            "sourmash: %d gathers, %d at a time (threads %d, memory budget %s at about "
            "%.1f GB per gather)",
            len(genomes),
            workers,
            params.threads,
            "none" if params.memory_gb is None else f"{params.memory_gb:g} GB",
            GATHER_GB,
        )
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
            hit = by_query.get(record_name(genome))
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
