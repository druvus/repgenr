"""NCBI Virus source: report parsing, vmetadata path, and vgenome records path."""

from __future__ import annotations

import json
from pathlib import Path

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import read_selection
from repgenr.viral import ncbi_virus
from repgenr.viral.ncbi_virus import VirusRecord, parse_report, write_records


def _report_line(acc, family, genus, species, length, segment="ANONYMOUS"):
    return json.dumps({
        "accession": acc, "length": length, "completeness": "COMPLETE", "segment": segment,
        "virus": {
            "organism_name": species, "tax_id": 11676,
            "lineage": [
                {"name": family, "tax_id": 1, "rank": "family"},
                {"name": genus, "tax_id": 2, "rank": "genus"},
                {"name": species, "tax_id": 3, "rank": "species"},
            ],
        },
    })


def test_parse_report_extracts_taxonomy_and_sanitizes() -> None:
    line = _report_line("NC_001802.1", "Retroviridae", "Lentivirus",
                        "Human immunodeficiency virus 1", 9181)
    (rec,) = parse_report([line])
    assert rec.accession == "NC_001802.1"
    assert rec.family == "Retroviridae" and rec.genus == "Lentivirus"
    assert rec.species == "Human-immunodeficiency-virus-1"  # spaces -> hyphen, round-trippable
    assert rec.length == 9181 and rec.completeness == "COMPLETE" and rec.segment == "ANONYMOUS"


def test_parse_download_report_schema_camelcase_and_no_ranks() -> None:
    # The download package's data_report.jsonl uses camelCase and has no rank
    # labels -- taxonomy must come from ICTV name suffixes.
    line = json.dumps({
        "accession": "AY358025.2", "length": 19111, "completeness": "COMPLETE",
        "isolate": {"name": "M/S.Africa/1975"},
        "virus": {
            "organismName": "Orthomarburgvirus marburgense", "taxId": 3052505,
            "lineage": [
                {"name": "Viruses"}, {"name": "Mononegavirales"},
                {"name": "Filoviridae"}, {"name": "Orthomarburgvirus"},
                {"name": "Orthomarburgvirus marburgense"},
            ],
        },
    })
    (rec,) = parse_report([line])
    assert rec.family == "Filoviridae"
    assert rec.genus == "Orthomarburgvirus"
    assert rec.species == "Orthomarburgvirus-marburgense"
    assert rec.isolate == "M/S.Africa/1975"


def _fake_records() -> list[VirusRecord]:
    return [
        VirusRecord("NC_001802.1", "11676", "HIV-1", "Retroviridae", "Lentivirus",
                    "Human-immunodeficiency-virus-1", 9181, "COMPLETE", "ANONYMOUS", ""),
        VirusRecord("AF033819.3", "11676", "HIV-1", "Retroviridae", "Lentivirus",
                    "Human-immunodeficiency-virus-1", 9100, "COMPLETE", "ANONYMOUS", ""),
    ]


def test_vmetadata_ncbi_virus_path(workdir: Path, monkeypatch) -> None:
    from repgenr.stages.vmetadata import VmetadataParams
    from repgenr.stages.vmetadata import run as vmetadata_run

    def fake_fetch(target, out_dir, **kw):
        recs = _fake_records()
        (out_dir / "download.fa").write_text(
            "".join(f">{r.accession} x\nACGT\n" for r in recs)
        )
        return recs

    monkeypatch.setattr(ncbi_virus, "fetch", fake_fetch)
    ctx = WorkdirContext(workdir, create=True)
    n = vmetadata_run(ctx, VmetadataParams(target="lentivirus", source="ncbi_virus"))
    assert n == 2
    dl = workdir / "virus_download_wd"
    assert (dl / "virus_records.json").exists()
    assert (dl / "metadata_base.tsv").exists()


def test_vgenome_records_path_canonical_names_and_selection(workdir: Path) -> None:
    from repgenr.stages.vgenome import VgenomeParams
    from repgenr.stages.vgenome import run as vgenome_run

    dl = workdir / "virus_download_wd"
    dl.mkdir(parents=True)
    recs = _fake_records()
    (dl / "download.fa").write_text("".join(f">{r.accession} desc\nACGTACGT\n" for r in recs))
    write_records(dl / "virus_records.json", recs)

    ctx = WorkdirContext(workdir, create=True)
    params = VgenomeParams(target_genus="lentivirus", length_all=True, no_outgroup=True)
    n = vgenome_run(ctx, params)
    assert n == 2
    names = sorted(p.name for p in ctx.genomes_dir.iterdir())
    assert names == [
        "Retroviridae_Lentivirus_Human-immunodeficiency-virus-1_AF033819.3.fasta",
        "Retroviridae_Lentivirus_Human-immunodeficiency-virus-1_NC_001802.1.fasta",
    ]
    rows = read_selection(workdir / "selection.tsv")
    assert {r.accession for r in rows} == {"NC_001802.1", "AF033819.3"}
    assert all(not r.is_outgroup for r in rows)


def test_vgenome_records_no_taxonomy_match(workdir: Path) -> None:
    import pytest

    from repgenr.core.errors import UserInputError
    from repgenr.stages.vgenome import VgenomeParams
    from repgenr.stages.vgenome import run as vgenome_run

    dl = workdir / "virus_download_wd"
    dl.mkdir(parents=True)
    recs = _fake_records()
    (dl / "download.fa").write_text("".join(f">{r.accession} d\nACGT\n" for r in recs))
    write_records(dl / "virus_records.json", recs)
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError):
        vgenome_run(ctx, VgenomeParams(target_genus="nonexistentgenus", no_outgroup=True))
