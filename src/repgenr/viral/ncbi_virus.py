"""NCBI Virus data source via the ``datasets`` CLI.

A single ``datasets download virus genome taxon`` call returns both the sequences
(``genomic.fna``) and structured per-sequence metadata (``data_report.jsonl``):
accession, ranked taxonomy, length, completeness, host, segment, isolate. This
replaces the BV-BRC FTP download + the separate NCBI Entrez taxonomy step with
one maintained dependency RepGenR already uses, and yields stable accessions so
viral genomes can adopt the canonical ``Family_Genus_species_Accession`` naming.
"""

from __future__ import annotations

import json
import logging
import re
import shutil
import zipfile
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities

DATASETS_CAPS = ToolCapabilities(
    name="datasets",
    required_binaries=(BinarySpec("datasets", version_args=("--version",)),),
    conda=("conda-forge::ncbi-datasets-cli",),
)

_ANONYMOUS = "ANONYMOUS"  # NCBI's value for non-segmented virus sequences


@dataclass
class VirusRecord:
    accession: str
    taxid: str
    organism: str
    family: str
    genus: str
    species: str
    length: int
    completeness: str
    segment: str
    isolate: str


def _sanitize(name: str) -> str:
    """Make a taxonomy name a single safe token (no spaces/underscores), so the
    canonical filename round-trips through ``parse_genome_filename``."""
    token = name.strip().replace(" ", "-").replace("_", "-")
    token = re.sub(r"[^A-Za-z0-9.-]", "", token)
    return token or "NA"


def _classify(lineage_names: list[str], organism: str) -> tuple[str, str, str]:
    """Derive (family, genus, species) from an NCBI Virus lineage.

    The report's lineage entries carry names but no usable rank labels, so we use
    ICTV name suffixes: family ends in ``-viridae``; a genus is a single-word name
    ending in ``-virus`` (species names are multi-word, e.g.
    'Orthomarburgvirus marburgense'); species is the organism (the leaf).
    """
    species = organism or (lineage_names[-1] if lineage_names else "")
    family = next((n for n in lineage_names if n.lower().endswith("viridae")), "")
    genus = ""
    for n in lineage_names:
        if " " not in n and n.lower().endswith("virus") and n != species:
            genus = n  # keep the last (closest to the leaf)
    if not genus and len(lineage_names) >= 2 and " " not in lineage_names[-2]:
        genus = lineage_names[-2]
    return family, genus, species


def parse_report(lines: Iterable[str]) -> list[VirusRecord]:
    """Parse a ``data_report.jsonl`` (or ``summary`` JSONL) stream into records.

    Handles both field casings: the download report uses ``organismName``/
    ``taxId``; the summary report uses ``organism_name``/``tax_id``.
    """
    records: list[VirusRecord] = []
    for raw in lines:
        raw = raw.strip()
        if not raw:
            continue
        row = json.loads(raw)
        virus = row.get("virus", {}) or {}
        organism = virus.get("organismName") or virus.get("organism_name") or ""
        taxid = virus.get("taxId") or virus.get("tax_id") or ""
        lineage_names = [e.get("name", "") for e in (virus.get("lineage") or []) if e.get("name")]
        family, genus, species = _classify(lineage_names, organism)
        records.append(
            VirusRecord(
                accession=row.get("accession", "") or "",
                taxid=str(taxid),
                organism=organism,
                family=_sanitize(family),
                genus=_sanitize(genus),
                species=_sanitize(species or organism),
                length=int(row.get("length", 0) or 0),
                completeness=(row.get("completeness", "") or "").upper(),
                segment=(row.get("segment", "") or _ANONYMOUS),
                isolate=((row.get("isolate") or {}).get("name", "") or ""),
            )
        )
    return records


def fetch(
    target: str,
    out_dir: Path,
    *,
    complete_only: bool = False,
    host: str | None = None,
    released_after: str | None = None,
    logger: logging.Logger,
    runner: Callable[..., int] = run_tool,
) -> list[VirusRecord]:
    """Download an NCBI Virus package for ``target`` and return its records.

    Writes the package sequences to ``out_dir/download.fa`` (headers
    ``>accession ...``). ``runner`` is injectable for tests.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    zip_path = out_dir / "ncbi_virus.zip"
    extract = out_dir / "ncbi_virus_pkg"
    if extract.exists():
        shutil.rmtree(extract)

    cmd = [
        "datasets", "download", "virus", "genome", "taxon", target,
        "--filename", str(zip_path), "--no-progressbar",
    ]
    if complete_only:
        cmd.append("--complete-only")
    if host:
        cmd += ["--host", host]
    if released_after:
        cmd += ["--released-after", released_after]
    runner(DATASETS_CAPS, cmd, logger=logger, log_prefix="datasets")

    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(extract)
    report = extract / "ncbi_dataset" / "data" / "data_report.jsonl"
    fna = extract / "ncbi_dataset" / "data" / "genomic.fna"
    if not report.exists() or not fna.exists():
        raise WorkdirError(
            f"NCBI Virus returned no genomes for taxon '{target}'. Check the name "
            "or loosen the filters (--complete-only/--host/--released-after)."
        )
    records = parse_report(report.read_text().splitlines())
    shutil.copyfile(fna, out_dir / "download.fa")
    zip_path.unlink(missing_ok=True)
    shutil.rmtree(extract, ignore_errors=True)
    logger.info("NCBI Virus: %d sequences for taxon '%s'", len(records), target)
    return records


def write_records(path: Path, records: list[VirusRecord]) -> None:
    path.write_text(json.dumps([asdict(r) for r in records]))


def read_records(path: Path) -> list[VirusRecord]:
    return [VirusRecord(**row) for row in json.loads(path.read_text())]
