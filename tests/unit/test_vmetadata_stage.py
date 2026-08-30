"""End-to-end tests of the vmetadata stage and the vgenome dispatch (WP2-5).

The NCBI Virus path runs through public ``run()`` against a faked datasets CLI
that writes a real zip package (data_report.jsonl + genomic.fna); the BV-BRC
path runs against a canned group FASTA and canned Entrez taxon data. Network
and FTP never happen.
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.stages import vgenome, vmetadata
from repgenr.stages.vgenome import VgenomeParams
from repgenr.stages.vmetadata import VmetadataParams
from repgenr.viral import ncbi_virus
from repgenr.viral.entrez import TAXNAMES_ORDERED


def _report_line(accession: str, organism: str, family: str, length: int) -> str:
    return json.dumps({
        "accession": accession,
        "length": length,
        "completeness": "complete",
        "segment": "",
        "isolate": {"name": f"iso-{accession}"},
        "virus": {
            "organismName": organism,
            "taxId": 10508,
            "lineage": [{"name": "Viruses"}, {"name": family}, {"name": "Mastadenovirus"}],
        },
    })


@pytest.fixture()
def fake_datasets(monkeypatch):
    """Fake run_tool_with_retries writing a real NCBI Virus zip package."""
    calls: list[list[str]] = []

    def fake(caps, cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        zip_path = Path(cmd[cmd.index("--filename") + 1])
        report = "\n".join([
            _report_line("NC_001.1", "Human adenovirus 1", "Adenoviridae", 34000),
            _report_line("NC_002.1", "Human adenovirus 2", "Adenoviridae", 35000),
        ])
        fna = ">NC_001.1 Human adenovirus 1\nACGT\n>NC_002.1 Human adenovirus 2\nACGT\n"
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("ncbi_dataset/data/data_report.jsonl", report)
            zf.writestr("ncbi_dataset/data/genomic.fna", fna)
        return 0

    monkeypatch.setattr(ncbi_virus, "run_tool_with_retries", fake)
    monkeypatch.setattr(vmetadata, "preflight", lambda caps: {"datasets": "16.0"},
                        raising=False)
    return calls


def test_ncbi_virus_end_to_end(tmp_path, fake_datasets, monkeypatch) -> None:
    import repgenr.core.plugins as plugins

    monkeypatch.setattr(plugins, "check_binaries", lambda specs: {})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    count = vmetadata.run(ctx, VmetadataParams(target="adenoviridae"))
    assert count == 2

    wd = ctx.workdir / "virus_download_wd"
    records = json.loads((wd / "virus_records.json").read_text(encoding="utf-8"))
    assert {r["accession"] for r in records} == {"NC_001.1", "NC_002.1"}
    assert (wd / "download.fa").exists()
    assert (ctx.workdir / "virus_metadata_base.tsv").exists()
    assert "vmetadata" in ctx.config.stages
    # filters forwarded to the datasets CLI
    flat = [tok for cmd in fake_datasets for tok in cmd]
    assert "adenoviridae" in flat


def test_ncbi_virus_filters_forwarded(tmp_path, fake_datasets, monkeypatch) -> None:
    import repgenr.core.plugins as plugins

    monkeypatch.setattr(plugins, "check_binaries", lambda specs: {})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = VmetadataParams(target="adenoviridae", complete_only=True,
                             host="homo sapiens", released_after="01/31/2024")
    vmetadata.run(ctx, params)
    flat = [tok for cmd in fake_datasets for tok in cmd]
    for token in ("--complete-only", "--host", "homo sapiens",
                  "--released-after", "01/31/2024"):
        assert token in flat


def test_ncbi_virus_empty_package_raises(tmp_path, monkeypatch) -> None:
    def empty(caps, cmd, *, logger, **kw):
        cmd = [str(c) for c in cmd]
        zip_path = Path(cmd[cmd.index("--filename") + 1])
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("ncbi_dataset/data/README.md", "empty")
        return 0

    monkeypatch.setattr(ncbi_virus, "run_tool_with_retries", empty)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="no genomes"):
        vmetadata.run(ctx, VmetadataParams(target="nosuchvirus"))


def test_missing_target_raises(tmp_path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="--target"):
        vmetadata.run(ctx, VmetadataParams(target=None))


# --- BV-BRC source ------------------------------------------------------------


_BVBRC_FASTA = (
    ">acc1 Strain one, complete genome [Adenovirus A | 111.1]\nACGTACGT\n"
    ">acc2 Strain two, complete genome [Adenovirus A | 111.2]\nACGTACGTAC\n"
    ">acc3 Strain three, partial sequence [Adenovirus B | 222.1]\nACGT\n"
)


def _taxdata(name: str, taxid: str) -> dict:
    levels = {level: {"name": None, "taxid": None} for level in TAXNAMES_ORDERED}
    levels["family"] = {"name": "Adenoviridae", "taxid": "10508"}
    levels["species"] = {"name": name, "taxid": taxid}
    return {"name": name, "taxdata": levels}


def test_bvbrc_end_to_end(tmp_path, monkeypatch) -> None:
    def fake_download_group(target, dest, logger):
        dest.write_text(_BVBRC_FASTA, encoding="utf-8")

    def fake_entrez(taxids, logger):
        data = {t: _taxdata(f"Species {t}", t) for t in taxids}
        return data, set(), {}

    monkeypatch.setattr(vmetadata, "_download_group", fake_download_group)
    monkeypatch.setattr(vmetadata, "get_taxon_data_from_entrez", fake_entrez)

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = VmetadataParams(target="adenoviridae", source="bvbrc")
    count = vmetadata.run(ctx, params)
    assert count == 2  # two distinct taxids seen

    wd = ctx.workdir / "virus_download_wd"
    base = (wd / "metadata_base.tsv").read_text(encoding="utf-8").splitlines()
    # only taxid 111 passes the "complete genome" filter, with 2 sequences
    assert len(base) == 2
    assert base[1].startswith("111\tAdenovirus A\t2\t8\t10")
    assert (wd / "metadata_ncbi.tsv").exists()
    assert (wd / "metadata_ncbi_taxnames_data.json").exists()
    assert (ctx.workdir / "virus_metadata_base.tsv").exists()
    assert "vmetadata" in ctx.config.stages


# --- vgenome dispatch ---------------------------------------------------------


def test_vgenome_missing_metadata_raises(tmp_path) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(WorkdirError, match="vmetadata"):
        vgenome.run(ctx, VgenomeParams())


def test_vgenome_dispatches_to_records_backend(tmp_path, monkeypatch) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    wd = ctx.workdir / "virus_download_wd"
    wd.mkdir(parents=True)
    (wd / "download.fa").write_text(">a\nACGT\n", encoding="utf-8")
    (wd / "virus_records.json").write_text("[]", encoding="utf-8")

    import repgenr.viral.selection as selection

    called = {}
    monkeypatch.setattr(
        selection, "run_records",
        lambda ctx, params, dwd, fasta, records, logger: called.setdefault("n", 7) or 7,
    )
    assert vgenome.run(ctx, VgenomeParams()) == 7
    assert called["n"] == 7


def test_vgenome_dispatches_to_bvbrc_backend(tmp_path, monkeypatch) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    wd = ctx.workdir / "virus_download_wd"
    wd.mkdir(parents=True)
    (wd / "download.fa").write_text(">a\nACGT\n", encoding="utf-8")
    (wd / "metadata_base.tsv").write_text("taxid\n", encoding="utf-8")
    (wd / "metadata_ncbi.tsv").write_text("taxid\n", encoding="utf-8")

    import repgenr.viral.bvbrc as bvbrc

    monkeypatch.setattr(
        bvbrc, "run_select",
        lambda ctx, params, fasta, base, ncbi, logger: 5,
    )
    assert vgenome.run(ctx, VgenomeParams()) == 5
