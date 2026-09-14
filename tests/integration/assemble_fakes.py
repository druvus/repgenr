"""Fakes shared by the assemble stage and step tests: an assembler that writes
one contig, a classifier returning canned GTDB lineages, a CheckM2 stand-in
and a reads.tsv row builder over tiny local FASTQ files."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

from repgenr.assemblers.base import Assembler, AssemblyResult
from repgenr.assemblers.base import registry as asm_registry
from repgenr.classifiers.base import Classification, Classifier, db_version
from repgenr.classifiers.base import registry as cls_registry
from repgenr.core.contracts import ReadRow
from repgenr.core.errors import ToolExecutionError
from repgenr.core.plugins import ToolCapabilities

GENOME = "ACGT" * 300  # 1200 bp, above the default contig floor


class FakeAssembler(Assembler):
    """Writes one contig from the reads it was given; records every call."""

    capabilities = ToolCapabilities(name="fakeasm")
    read_types = frozenset({"ILLUMINA", "OXFORD_NANOPORE"})
    calls: list[str] = []
    fail_runs: frozenset[str] = frozenset()

    def preflight(self) -> dict[str, str]:
        return {"fakeasm": "1.0"}

    def assemble(self, reads, out_dir, params, logger) -> AssemblyResult:  # noqa: ANN001
        type(self).calls.append(reads.run_accession)
        if reads.run_accession in type(self).fail_runs:
            raise ToolExecutionError(["fakeasm"], 1, "boom")
        assert all(f.exists() for f in reads.files)
        out_dir.mkdir(parents=True, exist_ok=True)
        contigs = out_dir / "raw.fa"
        contigs.write_text(f">node1\n{GENOME}\n>tiny\nACGT\n", encoding="utf-8")
        return AssemblyResult(contigs=contigs, tool_stats={"threads": params.threads})


def register_fake_assembler():
    asm_registry._load()
    asm_registry.register("fakeasm", FakeAssembler, replace=True)
    FakeAssembler.calls = []
    FakeAssembler.fail_runs = frozenset()


def unregister_fake_assembler():
    asm_registry._classes.pop("fakeasm", None)


class FakeClassifier(Classifier):
    """Registered classifier returning a canned GTDB lineage per genome (by file name)."""

    capabilities = ToolCapabilities(name="fakecls")
    lineages: dict[str, str] = {}

    def preflight(self) -> dict[str, str]:
        return {"fakecls": "1.0"}

    def classify(self, genomes, out_dir, params, logger):  # noqa: ANN001
        return {
            g.name: Classification(
                taxonomy=type(self).lineages[g.name],
                rank="species",
                score=0.9,
                db_version=db_version(params.db),
            )
            for g in genomes
            if g.name in type(self).lineages
        }


def register_fake_classifier():
    cls_registry._load()
    cls_registry.register("fakecls", FakeClassifier, replace=True)
    FakeClassifier.lineages = {}


def unregister_fake_classifier():
    cls_registry._classes.pop("fakecls", None)


def fake_checkm2(quality: dict[str, tuple[float, float]]):
    """A ``run_checkm2`` stand-in keyed by the genome file name."""

    def run_checkm2(genomes, out_dir, *, db, threads, logger):
        return {g.name: quality[g.name] for g in genomes if g.name in quality}

    return run_checkm2


def fastq(path: Path) -> tuple[str, str, int]:
    """Write a tiny gzipped FASTQ; return (local url, md5, bytes)."""
    data = gzip.compress(b"@r1\nACGT\n+\nIIII\n")
    path.write_bytes(data)
    return str(path), hashlib.md5(data).hexdigest(), len(data)


def read_row(tmp_path: Path, run: str, platform: str = "ILLUMINA", layout: str = "PAIRED", **over):
    files = [
        fastq(tmp_path / f"{run}_{i}.fastq.gz") for i in ((1, 2) if layout == "PAIRED" else (1,))
    ]
    base = dict(
        run_accession=run,
        biosample=f"SAM{run}",
        bioproject="PRJ1",
        organism="Francisella tularensis",
        taxid="263",
        platform=platform,
        instrument_model="MiSeq",
        layout=layout,
        bases=1_200_000,
        read_count=1000,
        family="Francisellaceae",
        genus="Francisella",
        species="tularensis",
        fastq_urls=tuple(f[0] for f in files),
        fastq_md5=tuple(f[1] for f in files),
        fastq_bytes=tuple(f[2] for f in files),
    )
    base.update(over)
    return ReadRow(**base)
