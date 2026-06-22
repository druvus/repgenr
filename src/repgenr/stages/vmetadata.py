"""viral metadata stage (BV-BRC + NCBI Entrez).

Ports ``vmetadata.py``: download a virus group's FASTA from the BV-BRC FTP,
derive per-taxid sequence-length statistics from the FASTA headers, enrich the
taxonomy via NCBI Entrez, and write the metadata tables plus a taxname->data
JSON (replacing the old pickle) under ``virus_download_wd/``. The viral genome
stage consumes these files.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import UTC, datetime
from ftplib import FTP, FTP_TLS
from pathlib import Path
from statistics import mean, median

from Bio import SeqIO

from ..core.context import WorkdirContext
from ..core.errors import UserInputError, WorkdirError
from ..viral.entrez import TAXNAMES_ORDERED, get_taxon_data_from_entrez

BVBRC_FTP = "ftp.bvbrc.org"
BVBRC_FTP_DIR = "viruses"


class _ReuseFTP_TLS(FTP_TLS):
    """FTPS client that reuses the control channel's TLS session for data.

    BV-BRC now requires SSL/TLS on the control channel and TLS session reuse on
    the data channel (a common vsftpd ``require_ssl_reuse`` setup); plain FTP and
    a vanilla FTP_TLS both fail ("550 SSL/TLS required" / "425 Unable to build
    data connection").
    """

    def ntransfercmd(self, cmd, rest=None):
        conn, size = FTP.ntransfercmd(self, cmd, rest)
        if self._prot_p:
            conn = self.context.wrap_socket(
                conn, server_hostname=self.host, session=self.sock.session
            )
        return conn, size


def _bvbrc_connect(timeout: int = 120) -> _ReuseFTP_TLS:
    ftp = _ReuseFTP_TLS(BVBRC_FTP, timeout=timeout)
    ftp.login(user="anonymous", passwd="anonymous")
    ftp.prot_p()
    return ftp


@dataclass
class VmetadataParams:
    target: str | None = None
    filter: str = "complete genome"
    list_targets: bool = False
    source: str = "ncbi_virus"  # ncbi_virus | bvbrc
    host: str | None = None  # ncbi_virus: restrict to a host species
    complete_only: bool = False  # ncbi_virus: only COMPLETE sequences
    released_after: str | None = None  # ncbi_virus: MM/DD/YYYY


def run(ctx: WorkdirContext, params: VmetadataParams) -> int:
    logger = ctx.logger
    if params.list_targets:
        _list_targets(logger)
        return 0
    if not params.target:
        raise UserInputError("Supply --target (e.g. adenoviridae) or --list.")

    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True, exist_ok=True)
    if params.source == "ncbi_virus":
        return _run_ncbi_virus(ctx, params, download_wd, logger)
    return _run_bvbrc(ctx, params, download_wd, logger)


def _run_ncbi_virus(ctx, params, download_wd, logger) -> int:
    from ..core.errors import MissingBinaryError
    from ..core.plugins import preflight
    from ..viral import ncbi_virus

    assert params.target is not None
    records = ncbi_virus.fetch(
        params.target, download_wd,
        complete_only=params.complete_only, host=params.host,
        released_after=params.released_after, logger=logger,
    )
    if not records:
        raise WorkdirError(f"NCBI Virus returned no genomes for '{params.target}'.")
    ncbi_virus.write_records(download_wd / "virus_records.json", records)
    _write_base_from_records(download_wd / "metadata_base.tsv", records)
    (ctx.workdir / "virus_metadata_base.tsv").write_text(
        (download_wd / "metadata_base.tsv").read_text()
    )

    # Best-effort provenance: fetch already used datasets, so this normally
    # resolves; never fail the stage just to record a version.
    try:
        tool_versions = preflight(ncbi_virus.DATASETS_CAPS)
    except MissingBinaryError:
        tool_versions = {}
    ctx.config.record_stage(
        "vmetadata",
        params={
            "source": "ncbi_virus", "target": params.target,
            "complete_only": params.complete_only, "host": params.host,
            "sequences": len(records),
        },
        tool_versions=tool_versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Viral metadata (NCBI Virus): %d sequences under %s", len(records), download_wd)
    return len(records)


def _write_base_from_records(path: Path, records) -> None:
    """taxid-keyed length stats (human-readable), grouped from NCBI Virus records."""
    from statistics import mean as _mean
    from statistics import median as _median

    by_taxid: dict[str, dict] = {}
    for r in records:
        g = by_taxid.setdefault(r.taxid, {"lens": [], "species": r.species, "organism": r.organism})
        g["lens"].append(r.length)
    header = ["taxid", "name", "num", "seq_min", "seq_max", "seq_med", "seq_mean", "description"]
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for taxid, g in sorted(by_taxid.items(), key=lambda x: len(x[1]["lens"]), reverse=True):
            lens = g["lens"]
            row = [
                taxid, g["species"], len(lens), min(lens), max(lens),
                int(_median(lens)), int(_mean(lens)), g["organism"],
            ]
            fo.write("\t".join(map(str, row)) + "\n")


def _run_bvbrc(ctx, params, download_wd, logger) -> int:
    target = params.target.lower()
    download_fa = download_wd / "download.fa"

    if not (download_fa.exists() and download_fa.stat().st_size > 0):
        _download_group(target, download_fa, logger)
    else:
        logger.info("Group FASTA already present; reusing %s", download_fa)

    base, all_taxids, taxid_bvbrc = _parse_fasta(download_fa, params.filter, logger)

    logger.info("Enriching %d taxids via NCBI Entrez", len(all_taxids))
    ncbi_data, missing, _alts = get_taxon_data_from_entrez(all_taxids, logger)
    taxnames_data = _taxnames_data(ncbi_data, missing, taxid_bvbrc)

    _write_base(download_wd / "metadata_base.tsv", base)
    _write_ncbi(download_wd / "metadata_ncbi.tsv", ncbi_data, base)
    serializable = {
        name: {**data, "datasets": sorted(data["datasets"])}
        for name, data in taxnames_data.items()
    }
    with open(download_wd / "metadata_ncbi_taxnames_data.json", "w") as fo:
        json.dump(serializable, fo)

    # copy human-readable tables to the workdir root with a virus_ prefix
    for name in ("metadata_base.tsv", "metadata_ncbi.tsv"):
        (ctx.workdir / f"virus_{name}").write_text((download_wd / name).read_text())

    ctx.config.record_stage(
        "vmetadata",
        params={"target": target, "filter": params.filter, "taxids": len(all_taxids)},
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Viral metadata written under %s", download_wd)
    return len(all_taxids)


def _list_targets(logger) -> None:
    with _bvbrc_connect() as ftp:
        ftp.cwd(BVBRC_FTP_DIR)
        targets = sorted(f.replace(".fna", "") for f in ftp.nlst() if f.endswith(".fna"))
    logger.info("Available targets:\n%s", "\n".join(targets))


def _download_group(target: str, dest: Path, logger) -> None:
    capitalized = target[0].upper() + target[1:]
    logger.info("Downloading %s from BV-BRC (FTPS)", capitalized)
    with _bvbrc_connect() as ftp:
        remote = f"{BVBRC_FTP_DIR}/{capitalized}.fna"
        ftp.sendcmd("TYPE I")
        try:
            remote_size = ftp.size(remote)
        except Exception as exc:
            raise WorkdirError(
                f"Could not find virus group '{capitalized}' at BV-BRC. "
                "Check the name (try --list), or download all with --target viruses."
            ) from exc
        # Binary transfer with an exact size check: ASCII mode + a 1000-byte
        # tolerance previously let a silently-truncated FASTA pass as complete.
        with open(dest, "wb") as fo:
            ftp.retrbinary("RETR " + remote, fo.write)
    local_size = dest.stat().st_size if dest.exists() else 0
    if remote_size and local_size != remote_size:
        dest.unlink(missing_ok=True)
        raise WorkdirError(
            f"Incomplete BV-BRC download: got {local_size} of {remote_size} bytes."
        )
    logger.info("Download finished (%d bytes)", local_size)


def _parse_fasta(path: Path, tag: str, logger):
    base: dict[str, dict] = {}
    all_taxids: set[str] = set()
    taxid_bvbrc: dict[str, set[str]] = {}
    counts: dict[str, int] = {}
    lens: dict[str, list[int]] = {}
    descriptions: dict[str, dict] = {}

    for entry in SeqIO.parse(str(path), "fasta"):
        bvbrc_id = entry.description.split("| ")[1].split("]")[0]
        taxid = bvbrc_id.split(".")[0]
        all_taxids.add(taxid)
        taxid_bvbrc.setdefault(taxid, set()).add(bvbrc_id)
        if tag not in entry.description:
            continue
        seq_len = len(entry.seq)
        counts[taxid] = counts.get(taxid, 0) + 1
        lens.setdefault(taxid, []).append(seq_len)
        if taxid not in descriptions:
            species = entry.description.split("[")[1].split(" | ")[0]
            full = entry.description.replace(entry.name, "").lstrip().split("[")[0].strip()
            descriptions[taxid] = {"species": species, "full": full}

    for taxid, count in counts.items():
        taxid_lens = lens[taxid]
        base[taxid] = {
            "num": count,
            "min": min(taxid_lens), "max": max(taxid_lens),
            "mean": int(mean(taxid_lens)), "median": int(median(taxid_lens)),
            "desc": descriptions[taxid],
        }
    logger.info("Found %d taxids; %d passed the '%s' filter", len(all_taxids), len(base), tag)
    return base, all_taxids, taxid_bvbrc


def _taxnames_data(ncbi_data, missing, taxid_bvbrc) -> dict[str, dict]:
    taxnames: dict[str, dict] = {}
    for taxid, data in ncbi_data.items():
        if taxid in missing:
            continue
        for level, level_data in data["taxdata"].items():
            name = level_data["name"]
            if name is None:
                continue
            entry = taxnames.setdefault(
                name, {"taxid": level_data["taxid"], "level": level, "datasets": set()}
            )
            for bvbrc_id in taxid_bvbrc.get(taxid, set()):
                entry["datasets"].add(bvbrc_id)
    return taxnames


def _write_base(path: Path, base: dict) -> None:
    header = ["taxid", "name", "num", "seq_min", "seq_max", "seq_med", "seq_mean", "description"]
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        for taxid, data in sorted(base.items(), key=lambda x: x[1]["num"], reverse=True):
            row = [
                taxid, data["desc"]["species"], data["num"], data["min"],
                data["max"], data["median"], data["mean"], data["desc"]["full"],
            ]
            fo.write("\t".join(map(str, row)) + "\n")


def _write_ncbi(path: Path, ncbi_data: dict, base: dict) -> None:
    header = ["taxid", "name", "num_with_tag"] + TAXNAMES_ORDERED + [
        f"{x}_taxid" for x in TAXNAMES_ORDERED
    ]
    with open(path, "w") as fo:
        fo.write("\t".join(header) + "\n")
        ordered = sorted(
            ncbi_data.items(),
            key=lambda x: base[x[0]]["num"] if x[0] in base else 0,
            reverse=True,
        )
        for taxid, data in ordered:
            num = base[taxid]["num"] if taxid in base else None
            row = [taxid, data["name"], num]
            row += [data["taxdata"][name]["name"] for name in TAXNAMES_ORDERED]
            row += [data["taxdata"][name]["taxid"] for name in TAXNAMES_ORDERED]
            fo.write("\t".join("" if v is None else str(v) for v in row) + "\n")
