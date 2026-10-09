"""sourmash gather against a GTDB sketch database, resolved to a lineage with
``sourmash tax genome`` and the matching lineages CSV.

A genome given in ``ClassifyParams.sketches`` is not sketched here: gather
reads its ``.sig.zip`` and selects the signature at the requested k-mer size
with ``-k``. Each query signature is named by the genome's record name.
"""

from __future__ import annotations

import csv
import logging
import os
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from ..core.containers import run_chain, run_tool
from ..core.contracts import record_name
from ..core.errors import UserInputError
from ..core.executors import parallel_map
from ..core.plugins import ToolCapabilities, parse_extra_int
from ..core.process import write_fofn
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


GATHER_KSIZE = 31
GATHER_THRESHOLD_BP = 50000


def gather_command(
    query: Path | str,
    db: Path | str,
    out_csv: Path,
    *,
    ksize: int = GATHER_KSIZE,
    threshold_bp: int = GATHER_THRESHOLD_BP,
) -> list[str | Path]:
    """The one ``sourmash gather`` call of a query signature against a reference
    database; shared by the classifier and the reads screen of assemble."""
    return [
        "sourmash",
        "gather",
        query,
        db,
        "-k",
        str(ksize),
        "--threshold-bp",
        str(threshold_bp),
        "-o",
        out_csv,
    ]


def tax_genome(
    caps: ToolCapabilities,
    gather_csvs: Sequence[Path],
    lineages: Path | str,
    base: Path,
    *,
    logger: logging.Logger,
    mounts: Sequence[str] = (),
    containment_threshold: float | None = None,
) -> Path:
    """Resolve gather results to lineages with ``sourmash tax genome``.

    Returns the classifications CSV. ``containment_threshold`` None keeps
    sourmash's default (0.1); the reads screen passes 0, since most k-mers of
    a reads sketch are sequencing errors that match nothing, so the unweighted
    fraction of a pure isolate stays well below 0.1.
    """
    extra: list[str | Path] = []
    if containment_threshold is not None:
        extra = ["--containment-threshold", f"{containment_threshold:g}"]
    # One gather CSV per genome: pass them in a list file (--from-file), not
    # on argv, where thousands of paths can exceed ARG_MAX. The list is
    # scratch, kept out of the working directory (the temp directory is bound
    # by the container backends).
    csv_dirs = sorted({os.path.dirname(os.path.abspath(c)) for c in gather_csvs})
    with tempfile.TemporaryDirectory(prefix="repgenr_tax_") as scratch:
        csv_list = write_fofn(gather_csvs, Path(scratch) / "gather_csvs.txt")
        run_tool(
            caps,
            [
                "sourmash",
                "tax",
                "genome",
                "--from-file",
                csv_list,
                "--taxonomy-csv",
                lineages,
                "--output-base",
                base,
                "--force",
                *extra,
            ],
            logger=logger,
            log_prefix="sourmash",
            extra_mounts=[*mounts, *csv_dirs],
        )
    return Path(str(base) + ".classifications.csv")


@dataclass(frozen=True)
class TaxRow:
    """One query of a ``tax genome`` table."""

    lineage: str
    rank: str
    fraction: float | None
    # The abundance-weighted fraction of the query at that rank (sourmash 4.4+).
    f_weighted: float | None


def read_tax_genome(path: Path) -> dict[str, TaxRow]:
    """Query name -> its row of a ``tax genome`` classifications CSV."""
    out: dict[str, TaxRow] = {}
    if not path.exists():
        return out
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo):
            if rec.get("status") in ("match", "nomatch", "below_threshold") and rec.get("lineage"):
                out[rec["query_name"]] = TaxRow(
                    rec["lineage"],
                    rec.get("rank", ""),
                    _opt_float(rec.get("fraction")),
                    _opt_float(rec.get("f_weighted_at_rank")),
                )
    return out


def read_gather_top(path: Path) -> tuple[str, float | None] | None:
    """The first (largest) match of a gather CSV: its name and f_unique_weighted.

    None when the gather found nothing above its threshold.
    """
    if not path.exists():
        return None
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo):
            return rec.get("name", ""), _opt_float(rec.get("f_unique_weighted"))
    return None


def _opt_float(value: str | None) -> float | None:
    try:
        return float(value) if value else None
    except ValueError:
        return None


class SourmashClassifier(Classifier):
    capabilities = sourmash_capabilities(
        default_params={
            "ksize": GATHER_KSIZE,
            "scaled": 1000,
            "threshold_bp": GATHER_THRESHOLD_BP,
        },
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
                    gather_command(sig, params.db, gather_csv, ksize=ksize, threshold_bp=threshold),
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
        table = tax_genome(
            self.capabilities,
            gathers,
            params.lineages,
            out_dir / "tax",
            logger=logger,
            mounts=mounts,
        )
        by_query = _parse_classifications(table)
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
    return {q: (r.lineage, r.rank, r.fraction) for q, r in read_tax_genome(path).items()}
