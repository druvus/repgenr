"""Force the sourmash DENSE dereplication path for one benchmark cell.

There is no CLI switch to disable the branchwater sparse path, so this driver
monkeypatches the availability probe and calls the adapter directly. Run under
``/usr/bin/time -l`` by the harness; child sourmash processes count toward the
measured RSS.

Usage: python -m benchmarks.dense_driver <genomes.fofn> <out_dir> <threads>
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from repgenr.core.contracts import write_clusters, write_genome_status
from repgenr.dereplicators import sourmash as sourmash_mod
from repgenr.dereplicators.base import DerepParams


def main() -> None:
    fofn, out_dir, threads = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
    genomes = [Path(line) for line in fofn.read_text(encoding="utf-8").splitlines() if line]

    logging.basicConfig(level=logging.INFO)
    sourmash_mod._branchwater_available = lambda caps, logger: False  # type: ignore[assignment]

    adapter = sourmash_mod.SourmashDereplicator()
    result = adapter.dereplicate(
        genomes, out_dir, DerepParams(threads=threads), logging.getLogger("dense")
    )
    write_clusters(out_dir / "clusters.tsv", result.clusters)
    write_genome_status(out_dir / "genome_status.tsv", result.genome_status)
    print(f"dense: {len(result.representatives)} representatives")


if __name__ == "__main__":
    main()
