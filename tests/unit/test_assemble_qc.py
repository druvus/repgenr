"""CheckM2 quality assessment and the classifier family (sourmash gather)."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.classifiers.base import ClassifyParams, registry
from repgenr.stages.assemble_qc import CHECKM2_CAPS, parse_checkm2_report, run_checkm2

_LOG = logging.getLogger("test")

_REPORT = (
    "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\n"
    "Fam_Gen_sp_SRR1\t99.12\t0.35\tNeural Network (Specific Model)\t11\n"
    "Fam_Gen_sp_SRR2\t41.00\t12.50\tGradient Boost (General Model)\t11\n"
)


def test_parse_checkm2_report(tmp_path: Path) -> None:
    report = tmp_path / "quality_report.tsv"
    report.write_text(_REPORT, encoding="utf-8")
    assert parse_checkm2_report(report) == {
        "Fam_Gen_sp_SRR1": (99.12, 0.35),
        "Fam_Gen_sp_SRR2": (41.0, 12.5),
    }


def test_run_checkm2_passes_db_threads_and_reads_the_report(tmp_path: Path, monkeypatch) -> None:
    import repgenr.stages.assemble_qc as qc

    seen = {}

    def fake_run_tool(caps, command, *, logger, **kwargs):
        cmd = [str(c) for c in command]
        seen["cmd"] = cmd
        seen["mounts"] = [str(m) for m in kwargs.get("extra_mounts", ())]
        out = Path(cmd[cmd.index("--output-directory") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "quality_report.tsv").write_text(_REPORT, encoding="utf-8")
        return 0

    monkeypatch.setattr(qc, "run_tool", fake_run_tool)
    genomes = [tmp_path / "Fam_Gen_sp_SRR1.fasta", tmp_path / "Fam_Gen_sp_SRR2.fasta"]
    for g in genomes:
        g.write_text(">a\nACGT\n", encoding="utf-8")
    quality = run_checkm2(genomes, tmp_path / "qc", db=tmp_path / "db.dmnd", threads=3, logger=_LOG)
    assert quality["Fam_Gen_sp_SRR1.fasta"] == (99.12, 0.35)
    cmd = seen["cmd"]
    assert cmd[:2] == ["checkm2", "predict"]
    assert cmd[cmd.index("--threads") + 1] == "3"
    assert cmd[cmd.index("--database_path") + 1] == str(tmp_path / "db.dmnd")
    assert "-x" in cmd and cmd[cmd.index("-x") + 1] == "fasta"
    assert str(tmp_path / "db.dmnd") in seen["mounts"] or str(tmp_path) in seen["mounts"]
    assert CHECKM2_CAPS.container.startswith("quay.io/biocontainers/checkm2:")


_GATHER_HEADER = "intersect_bp,f_orig_query,f_match,f_unique_to_query,name,query_name\n"
_CLASSIFICATION = (
    "query_name,status,rank,fraction,lineage,query_md5,query_filename\n"
    "Fam_Gen_sp_SRR1,match,species,0.94,"
    "d__Bacteria;p__Bacillota_B;c__Mycoplasmoidia;o__Mycoplasmoidales;f__Mycoplasmoidaceae;"
    "g__Mycoplasmoides;s__Mycoplasmoides genitalium,abc,x.fasta\n"
)


def test_sourmash_classifier_chain_and_parse(tmp_path: Path, monkeypatch) -> None:
    import repgenr.classifiers.sourmash as sm

    recorded: list[list[str]] = []

    def fake_run_chain(caps, steps, *, logger, **kwargs):
        for _prefix, argv in steps:
            cmd = [str(c) for c in argv]
            recorded.append(cmd)
            if cmd[:3] == ["sourmash", "tax", "genome"]:
                base = Path(cmd[cmd.index("--output-base") + 1])
                base.parent.mkdir(parents=True, exist_ok=True)
                Path(str(base) + ".classifications.csv").write_text(
                    _CLASSIFICATION, encoding="utf-8"
                )
            elif cmd[:2] == ["sourmash", "gather"]:
                Path(cmd[cmd.index("-o") + 1]).write_text(_GATHER_HEADER, encoding="utf-8")
        return 0

    monkeypatch.setattr(sm, "run_chain", fake_run_chain)
    genome = tmp_path / "Fam_Gen_sp_SRR1.fasta"
    genome.write_text(">a\nACGT\n", encoding="utf-8")
    db, lineages = tmp_path / "gtdb-rs226-reps.k31-sc10k.sig.zip", tmp_path / "lineages.csv"
    db.write_bytes(b"zip"), lineages.write_text("ident,superkingdom\n", encoding="utf-8")
    result = registry.create("sourmash").classify(
        [genome], tmp_path / "cls", ClassifyParams(db=db, lineages=lineages, threads=2), _LOG
    )
    cls = result["Fam_Gen_sp_SRR1.fasta"]
    assert cls.taxonomy.endswith("s__Mycoplasmoides genitalium")
    assert cls.rank == "species" and cls.score == 0.94
    assert cls.db_version == "gtdb-rs226-reps.k31-sc10k"
    flat = [tok for cmd in recorded for tok in cmd]
    assert "sketch" in flat and "k=31,scaled=1000" in flat
    assert "gather" in flat and "--threshold-bp" in flat and str(db) in flat
    assert "--taxonomy-csv" in flat and str(lineages) in flat


def test_classifier_family_is_registered_and_listed() -> None:
    from typer.testing import CliRunner

    from repgenr.cli.main import app

    assert "sourmash" in registry.names()
    result = CliRunner().invoke(app, ["list-tools"])
    assert result.exit_code == 0
    assert any(
        ln.startswith("classifiers:") and "sourmash" in ln for ln in result.output.splitlines()
    )
