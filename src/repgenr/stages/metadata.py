"""metadata stage: select a taxon's genomes from GTDB.

Two data sources:

* ``tsv`` (default): download the full GTDB metadata table and parse it. Robust
  and release-pinned, but downloads the whole table.
* ``api``: query the GTDB API (https://gtdb-api.ecogenomic.org), fetching only
  the target taxon's genomes. Much smaller transfer; uses the API's current
  release (``--release``/``--version`` are advisory for this source).

Either way the selection is recorded in the SQLite manifest plus ``repgenr.yaml``
provenance.
"""

from __future__ import annotations

import gzip
import shutil
import tarfile
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import requests

from ..core.context import WorkdirContext
from ..core.errors import UserInputError, WorkdirError
from ..core.manifest import GenomeRecord

TAXONOMY = ("domain", "phylum", "class", "family", "genus", "species")
GTDB_API_BASE = "https://gtdb-api.ecogenomic.org"
# Single-letter GTDB rank prefixes used by the API and taxon strings.
_RANK_PREFIX = {"family": "f", "genus": "g", "species": "s"}


@dataclass
class MetadataParams:
    dataset: str  # all | rep
    level: str  # family | genus | species
    release: str | None = None  # required for tsv source
    version: str | None = None  # bac120 | ar53; required for tsv source
    source: str = "tsv"  # tsv | api
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

    if params.source == "api":
        selected, outgroup = _select_via_api(params, logger)
    else:
        selected, outgroup = _select_via_tsv(ctx, params, logger)

    logger.info("Selected %d genomes; outgroup: %s", len(selected), outgroup.accession)
    _populate_manifest(ctx, selected, outgroup)
    _write_outgroup_file(ctx, outgroup.accession)

    ctx.config.record_stage(
        "metadata",
        params={
            "source": params.source,
            "release": params.release, "version": params.version, "dataset": params.dataset,
            "level": params.level, "target_family": params.target_family,
            "target_genus": params.target_genus, "target_species": params.target_species,
            "selected_count": len(selected), "outgroup": outgroup.accession,
        },
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(selected)


def _validate(params: MetadataParams) -> None:
    if params.source not in ("tsv", "api"):
        raise UserInputError("--source must be 'tsv' or 'api'.")
    if params.source == "tsv":
        if not params.release or "." not in params.release:
            raise UserInputError("tsv source needs --release like '207.0' (major.minor).")
        if not params.version:
            raise UserInputError("tsv source needs --version (bac120 or ar53).")
    if not (params.target_genus or params.target_family):
        raise UserInputError("Supply --target-genus or --target-family.")
    if params.level == "species" and not params.target_species:
        raise UserInputError("Level 'species' needs --target-species.")
    if params.level == "genus" and not params.target_genus:
        raise UserInputError("Level 'genus' needs --target-genus.")
    if params.level == "family" and not params.target_family:
        raise UserInputError("Level 'family' needs --target-family.")


def _select_via_tsv(
    ctx: WorkdirContext, params: MetadataParams, logger
) -> tuple[list[GenomeRecord], GenomeRecord]:
    metadata_file = _obtain_metadata(ctx, params, logger)
    accessions = _parse_metadata(metadata_file, params, logger)

    target_levels = _target_levels(accessions, params)
    if not target_levels:
        raise UserInputError(
            f"Target not found in GTDB: family={params.target_family} "
            f"genus={params.target_genus} species={params.target_species} at level {params.level}"
        )
    selected = _select(accessions, target_levels, params.limit)
    outgroup_acc, outgroup_data = _pick_outgroup(accessions, selected, target_levels, params)

    records = [_record_from_tax(acc, data["tax"]) for acc, data in selected.items()]
    outgroup = _record_from_tax(outgroup_acc, outgroup_data["tax"], is_outgroup=True)
    return records, outgroup


def _record_from_tax(accession: str, tax: dict, is_outgroup: bool = False) -> GenomeRecord:
    return GenomeRecord(
        accession=accession, source="gtdb", is_outgroup=is_outgroup,
        family=tax["family"], genus=tax["genus"], species=tax["species"],
    )


def _obtain_metadata(ctx: WorkdirContext, params: MetadataParams, logger) -> Path:
    if params.metadata_path and Path(params.metadata_path).exists():
        logger.info("Using provided metadata: %s", params.metadata_path)
        return Path(params.metadata_path)

    ctx.workdir.mkdir(parents=True, exist_ok=True)
    assert params.release is not None  # _validate guarantees this for the tsv source
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


def _populate_manifest(ctx, selected: list[GenomeRecord], outgroup: GenomeRecord) -> None:
    ctx.manifest.upsert_many([*selected, outgroup])


def _write_outgroup_file(ctx, outgroup_acc) -> None:
    (ctx.workdir / "outgroup_accession.txt").write_text(outgroup_acc + "\n")


# --- GTDB API source --------------------------------------------------------

def _api_get(path: str, params: dict | None = None) -> dict:
    url = f"{GTDB_API_BASE}{path}"
    try:
        resp = requests.get(url, params=params, timeout=120)
        resp.raise_for_status()
        return resp.json()
    except requests.RequestException as exc:
        raise WorkdirError(f"GTDB API request failed: {url} ({exc})") from exc


def _target_taxon(params: MetadataParams) -> str:
    """Build the GTDB taxon string for the selection level (e.g. g__Francisella)."""
    prefix = _RANK_PREFIX[params.level]
    if params.level == "family":
        name = params.target_family
    elif params.level == "genus":
        name = params.target_genus
    else:  # species
        name = f"{params.target_genus} {params.target_species}"
    assert name is not None  # guaranteed by _validate for the chosen level
    return f"{prefix}__{_capitalize_taxon(name)}"


def _capitalize_taxon(name: str) -> str:
    """GTDB names capitalize the first word only (e.g. 'francisella tularensis')."""
    parts = name.strip().split()
    if not parts:
        return name
    parts[0] = parts[0].capitalize()
    return " ".join(parts)


def _normalize_api_tax(row: dict) -> dict:
    """Normalize an API row's gtdb* fields to the manifest tax dict.

    Strips the rank prefix and applies the same cleanup as the TSV path:
    species has the genus removed, spaces dropped, underscores -> hyphens.
    """
    def strip(value: str) -> str:
        return value.split("__", 1)[1] if "__" in value else value

    genus = strip(row.get("gtdbGenus", ""))
    species = strip(row.get("gtdbSpecies", ""))
    family = strip(row.get("gtdbFamily", ""))
    species = species.replace(genus, "").replace(" ", "").replace("_", "-")
    genus = genus.replace(" ", "").replace("_", "-")
    family = family.replace(" ", "").replace("_", "-")
    return {"family": family, "genus": genus, "species": species,
            "gtdbFamily": row.get("gtdbFamily", ""), "gtdbGenus": row.get("gtdbGenus", "")}


def _select_via_api(
    params: MetadataParams, logger
) -> tuple[list[GenomeRecord], GenomeRecord]:
    taxon = _target_taxon(params)
    logger.info("Querying GTDB API for genomes in %s", taxon)
    sp_reps = params.dataset == "rep"
    rows = _api_genomes_detail(taxon, sp_reps, logger)
    if not rows:
        raise UserInputError(f"GTDB API returned no genomes for {taxon}. Check the target name.")

    if params.limit:
        rows = rows[: params.limit]

    records: list[GenomeRecord] = []
    selected_acc: set[str] = set()
    for row in rows:
        acc = row["gid"]
        records.append(_record_from_tax(acc, _normalize_api_tax(row)))
        selected_acc.add(acc)

    outgroup = _select_outgroup_via_api(params, rows, selected_acc, logger)
    return records, outgroup


def _api_genomes_detail(taxon: str, sp_reps_only: bool, logger) -> list[dict]:
    encoded = urllib.parse.quote(taxon, safe="")
    data = _api_get(f"/taxon/{encoded}/genomes-detail", {"sp_reps_only": str(sp_reps_only).lower()})
    return data.get("rows", [])


def _select_outgroup_via_api(
    params: MetadataParams, rows: list[dict], selected_acc: set[str], logger
) -> GenomeRecord:
    # Explicit outgroup: fetch its card for taxonomy.
    if params.outgroup_accession:
        card = _api_get(f"/genome/{urllib.parse.quote(params.outgroup_accession, safe='')}/card")
        tax_row = card.get("metadataTaxonomy", {})
        return _record_from_tax(
            params.outgroup_accession, _normalize_api_tax(tax_row), is_outgroup=True
        )

    # Otherwise pick a representative from the parent taxon (one rank up) that is
    # not part of the selection. Use the same rank order as the TSV path.
    parent_rank = {"species": "genus", "genus": "family", "family": "class"}[params.level]
    parent_field = {"genus": "gtdbGenus", "family": "gtdbFamily", "class": "gtdbClass"}[parent_rank]
    parent_taxon = rows[0].get(parent_field)
    if not parent_taxon:
        raise WorkdirError("Could not determine a parent taxon for outgroup selection.")

    logger.info("Selecting outgroup from parent taxon %s", parent_taxon)
    parent_rows = _api_genomes_detail(parent_taxon, sp_reps_only=True, logger=logger)
    for row in parent_rows:
        if row["gid"] in selected_acc:
            continue
        if not row.get("gtdbIsRep", False):
            continue
        # must be outside the selected sub-taxon
        if row.get(_rank_field(params.level)) == rows[0].get(_rank_field(params.level)):
            continue
        return _record_from_tax(row["gid"], _normalize_api_tax(row), is_outgroup=True)
    raise WorkdirError(
        "Could not determine an outgroup via the API; specify --outgroup-accession."
    )


def _rank_field(level: str) -> str:
    return {"family": "gtdbFamily", "genus": "gtdbGenus", "species": "gtdbSpecies"}[level]
