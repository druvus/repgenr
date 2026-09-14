"""sourmash gather against a GTDB sketch database, resolved to a lineage with
``sourmash tax genome`` and the matching lineages CSV."""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_chain
from ..core.errors import UserInputError
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
        results: dict[str, Classification] = {}
        for genome in genomes:
            work = out_dir / genome.stem
            work.mkdir(parents=True, exist_ok=True)
            sig = work / "query.sig"
            gather = work / "gather.csv"
            base = work / "tax"
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
                        gather,
                    ],
                ),
                (
                    "sourmash",
                    [
                        "sourmash",
                        "tax",
                        "genome",
                        "--gather-csv",
                        gather,
                        "--taxonomy-csv",
                        params.lineages,
                        "--output-base",
                        base,
                        "--force",
                    ],
                ),
            ]
            run_chain(
                self.capabilities,
                steps,
                logger=logger,
                extra_mounts=[
                    str(Path(params.db).resolve().parent),
                    str(Path(params.lineages).resolve().parent),
                    str(genome.resolve().parent),
                ],
            )
            hit = _parse_classification(Path(str(base) + ".classifications.csv"))
            if hit is not None:
                taxonomy, rank, score = hit
                results[genome.name] = Classification(taxonomy, rank, score, version)
            else:
                logger.warning("%s: no GTDB match above the threshold", genome.name)
        return results


def _parse_classification(path: Path) -> tuple[str, str, float | None] | None:
    if not path.exists():
        return None
    with open(path, encoding="utf-8", newline="") as fo:
        for rec in csv.DictReader(fo):
            if rec.get("status") in ("match", "nomatch", "below_threshold") and rec.get("lineage"):
                score = rec.get("fraction")
                return rec["lineage"], rec.get("rank", ""), float(score) if score else None
    return None
