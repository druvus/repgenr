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
from ftplib import FTP
from pathlib import Path
from statistics import mean, median

from Bio import SeqIO

from ..core.context import WorkdirContext
from ..core.errors import UserInputError, WorkdirError
from ..viral.entrez import TAXNAMES_ORDERED, get_taxon_data_from_entrez

BVBRC_FTP = "ftp.bvbrc.org"
BVBRC_FTP_DIR = "viruses"


@dataclass
class VmetadataParams:
    target: str | None = None
    filter: str = "complete genome"
    list_targets: bool = False


def run(ctx: WorkdirContext, params: VmetadataParams) -> int:
    logger = ctx.logger
    if params.list_targets:
        _list_targets(logger)
        return 0
    if not params.target:
        raise UserInputError("Supply --target (e.g. adenoviridae) or --list.")

    target = params.target.lower()
    download_wd = ctx.workdir / "virus_download_wd"
    download_wd.mkdir(parents=True, exist_ok=True)
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
    with FTP(BVBRC_FTP) as ftp:
        ftp.login(user="anonymous", passwd="")
        ftp.cwd(BVBRC_FTP_DIR)
        targets = sorted(f.replace(".fna", "") for f in ftp.nlst() if f.endswith(".fna"))
    logger.info("Available targets:\n%s", "\n".join(targets))


def _download_group(target: str, dest: Path, logger) -> None:
    capitalized = target[0].upper() + target[1:]
    logger.info("Downloading %s from BV-BRC FTP", capitalized)
    with FTP(BVBRC_FTP) as ftp:
        ftp.login(user="anonymous", passwd="")
        remote = f"{BVBRC_FTP_DIR}/{capitalized}.fna"
        ftp.sendcmd("TYPE I")
        try:
            remote_size = ftp.size(remote)
        except Exception as exc:
            raise WorkdirError(
                f"Could not find virus group '{capitalized}' at BV-BRC. "
                "Check the name (try --list), or download all with --target viruses."
            ) from exc
        ftp.sendcmd("TYPE A")
        with open(dest, "w") as fo:
            ftp.retrlines("RETR " + remote, lambda line: fo.write(line + "\n"))
    local_size = dest.stat().st_size if dest.exists() else 0
    if remote_size and abs(remote_size - local_size) >= 1000:
        raise WorkdirError(
            f"Downloaded size {local_size} differs from remote {remote_size} by >1000 bytes."
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
