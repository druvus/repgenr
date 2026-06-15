"""metadata stage: select a taxon's genomes from GTDB.

Ports ``metadata.py``: download the GTDB metadata table, select accessions whose
GTDB taxonomy matches the target down to the requested level, determine an
outgroup, and record the selection in the SQLite manifest plus ``repgenr.yaml``
provenance. Replaces the old ``str(dict)`` state files.
"""

from __future__ import annotations

import gzip
import shutil
import tarfile
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.errors import UserInputError, WorkdirError
from ..core.manifest import GenomeRecord

TAXONOMY = ("domain", "phylum", "class", "family", "genus", "species")


@dataclass
class MetadataParams:
    release: str
    version: str  # bac120 | ar53
    dataset: str  # all | rep
    level: str  # family | genus | species
    target_family: str | None = None
    target_genus: str | None = None
    target_species: str | None = None
    outgroup_accession: str | None = None
    metadata_path: str | None = None
    nodownload: bool = False
    limit: int | None = None


def run(ctx: WorkdirContext, params: MetadataParams) -> int:
    logger = ctx.logger
    _validate(params)

    metadata_file = _obtain_metadata(ctx, params, logger)
    accessions = _parse_metadata(metadata_file, params, logger)

    target_levels = _target_levels(accessions, params)
    if not target_levels:
        raise UserInputError(
            f"Target not found in GTDB: family={params.target_family} "
            f"genus={params.target_genus} species={params.target_species} at level {params.level}"
        )

    selected = _select(accessions, target_levels, params.limit)
    logger.info("Selected %d genomes", len(selected))

    outgroup_acc, outgroup_data = _pick_outgroup(accessions, selected, target_levels, params)
    logger.info("Outgroup: %s", outgroup_acc)

    _populate_manifest(ctx, selected, outgroup_acc, outgroup_data)
    _write_outgroup_file(ctx, outgroup_acc)

    ctx.config.record_stage(
        "metadata",
        params={
            "release": params.release, "version": params.version, "dataset": params.dataset,
            "level": params.level, "target_family": params.target_family,
            "target_genus": params.target_genus, "target_species": params.target_species,
            "selected_count": len(selected), "outgroup": outgroup_acc,
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(selected)


def _validate(params: MetadataParams) -> None:
    if "." not in params.release:
        raise UserInputError("Release must look like '207.0' (major.minor).")
    if not (params.target_genus or params.target_family):
        raise UserInputError("Supply --target-genus or --target-family.")
    if params.level == "species" and not params.target_species:
        raise UserInputError("Level 'species' needs --target-species.")
    if params.level == "genus" and not params.target_genus:
        raise UserInputError("Level 'genus' needs --target-genus.")
    if params.level == "family" and not params.target_family:
        raise UserInputError("Level 'family' needs --target-family.")


def _obtain_metadata(ctx: WorkdirContext, params: MetadataParams, logger) -> Path:
    if params.metadata_path and Path(params.metadata_path).exists():
        logger.info("Using provided metadata: %s", params.metadata_path)
        return Path(params.metadata_path)

    ctx.workdir.mkdir(parents=True, exist_ok=True)
    major = int(float(params.release))
    base = (
        f"https://data.gtdb.ecogenomic.org/releases/release{major}/"
        f"{params.release}/{params.version}_metadata_r{major}"
    )
    for ext in (".tar.gz", ".tsv.gz"):
        url = base + ext
        dest = ctx.workdir / Path(url).name
        if params.nodownload and dest.exists():
            logger.info("Using previously downloaded %s", dest.name)
            return dest
        try:
            logger.info("Downloading %s", url)
            urllib.request.urlretrieve(url, dest)
            return dest
        except Exception as exc:  # try the next naming scheme
            logger.warning("Download failed (%s); trying next target", exc)
    raise WorkdirError("Could not download GTDB metadata; check release/version.")


def _open_metadata(path: Path, workdir: Path):
    if path.name.endswith(".tar.gz"):
        tsv_gz = workdir / path.name.replace(".tar.gz", ".tsv.gz")
        # Stream the member straight out of the tarball into a gzipped TSV.
        # Avoid extract-to-directory + rmtree, which fails on exFAT/NTFS volumes
        # (shutil.rmtree's dir_fd traversal is unsupported there).
        with tarfile.open(path) as tar:
            member = next((m for m in tar.getmembers() if m.name.endswith(".tsv")), None)
            if member is None:
                raise WorkdirError("No .tsv inside GTDB tarball")
            source = tar.extractfile(member)
            if source is None:
                raise WorkdirError("Could not read .tsv from GTDB tarball")
            with source, gzip.open(tsv_gz, "wb") as fo:
                shutil.copyfileobj(source, fo)
        path = tsv_gz
    return gzip.open(path, "rt")


def _parse_metadata(path: Path, params: MetadataParams, logger) -> dict[str, dict]:
    logger.info("Parsing GTDB metadata")
    accessions: dict[str, dict] = {}
    with _open_metadata(path, path.parent) as fo:
        header = fo.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fo:
            fields = line.rstrip("\n").split("\t")
            acc_raw = fields[idx["accession"]]
            accession = acc_raw.replace("GB_", "").replace("RS_", "")
            rep = fields[idx["gtdb_genome_representative"]]
            is_rep = acc_raw == rep

            if params.dataset == "rep" and not is_rep:
                if not (params.outgroup_accession and accession == params.outgroup_accession):
                    continue

            tax = _parse_taxonomy(fields[idx["gtdb_taxonomy"]])
            accessions[accession] = {
                "accession": accession,
                "accession_ncbi": fields[idx.get("ncbi_genbank_assembly_accession", 0)],
                "tax": tax,
                "is_rep": is_rep,
            }
    logger.info("Parsed %d accessions", len(accessions))
    return accessions


def _parse_taxonomy(raw: str) -> dict[str, str]:
    tax = {level: "" for level in TAXONOMY}
    for chunk in raw.split(";"):
        for level in TAXONOMY:
            key = level[0] + "__"
            if chunk.startswith(key):
                tax[level] = chunk[len(key) :]
    # normalize species/genus/family like the legacy code
    tax["species"] = tax["species"].replace(tax["genus"], "").replace(" ", "").replace("_", "-")
    for level in ("genus", "family"):
        tax[level] = tax[level].replace(" ", "").replace("_", "-")
    return tax


def _target_levels(accessions: dict[str, dict], params: MetadataParams) -> dict[str, str]:
    def norm(value: str | None) -> str | None:
        return value.lower().replace("_", "-") if value else None

    tf, tg, ts = norm(params.target_family), norm(params.target_genus), norm(params.target_species)
    for data in accessions.values():
        tax = data["tax"]
        matched = False
        if ts and tg:
            matched = tax["species"].lower() == ts and tax["genus"].lower() == tg
        elif tg:
            matched = tax["genus"].lower() == tg
        elif tf:
            matched = tax["family"].lower() == tf
        if matched:
            levels: dict[str, str] = {}
            for level in TAXONOMY:
                levels[level] = tax[level]
                if level == params.level:
                    break
            return levels
    return {}


def _select(accessions: dict[str, dict], target_levels: dict[str, str], limit: int | None):
    selected = {}
    for acc, data in accessions.items():
        if all(data["tax"][lvl] == val for lvl, val in target_levels.items()):
            selected[acc] = data
            if limit and len(selected) >= limit:
                break
    return selected


def _pick_outgroup(accessions, selected, target_levels, params):
    if params.outgroup_accession:
        if params.outgroup_accession not in accessions:
            raise UserInputError(
                f"Outgroup accession {params.outgroup_accession} not in GTDB metadata."
            )
        return params.outgroup_accession, accessions[params.outgroup_accession]

    # one level above the selection level
    levels = list(target_levels)
    upper = levels[-2] if len(levels) >= 2 else levels[-1]
    upper_val = next(iter(selected.values()))["tax"][upper]
    for acc, data in accessions.items():
        if acc in selected:
            continue
        if data["tax"][upper] == upper_val and data["is_rep"]:
            return acc, data
    raise WorkdirError("Could not determine an outgroup; specify --outgroup-accession.")


def _populate_manifest(ctx, selected, outgroup_acc, outgroup_data) -> None:
    manifest = ctx.manifest
    records = []
    for acc, data in selected.items():
        tax = data["tax"]
        records.append(
            GenomeRecord(
                accession=acc, source="gtdb",
                family=tax["family"], genus=tax["genus"], species=tax["species"],
            )
        )
    tax = outgroup_data["tax"]
    records.append(
        GenomeRecord(
            accession=outgroup_acc, source="gtdb", is_outgroup=True,
            family=tax["family"], genus=tax["genus"], species=tax["species"],
        )
    )
    manifest.upsert_many(records)


def _write_outgroup_file(ctx, outgroup_acc) -> None:
    (ctx.workdir / "outgroup_accession.txt").write_text(outgroup_acc + "\n")
