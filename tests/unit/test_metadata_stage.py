"""End-to-end tests of the metadata stage over a fixture GTDB TSV slice (WP2-4).

The full GTDB table download is replaced by a small gzipped TSV; the API source
runs against canned ``_api_get`` responses. Both paths run through the public
``run()`` and are asserted on the published contract: selection.tsv, the
outgroup file, the manifest, and the recorded config stage.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.contracts import SELECTION_TSV
from repgenr.core.errors import UserInputError, WorkdirError
from repgenr.stages import metadata
from repgenr.stages.metadata import MetadataParams

_TAX_PREFIX = "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Francisellales"


def _tax(family: str, genus: str, species: str) -> str:
    return f"{_TAX_PREFIX};f__{family};g__{genus};s__{genus} {species}"


# accession, representative accession, taxonomy, ncbi accession
_ROWS = [
    ("GCF_000001.1", "GCF_000001.1", _tax("Francisellaceae", "Francisella", "tularensis")),
    ("GCF_000002.1", "GCF_000001.1", _tax("Francisellaceae", "Francisella", "tularensis")),
    ("GCF_000003.1", "GCF_000003.1", _tax("Francisellaceae", "Francisella", "tularensis")),
    # same genus, different species: the natural outgroup for a species run
    ("GCF_000010.1", "GCF_000010.1", _tax("Francisellaceae", "Francisella", "philomiragia")),
    # different family entirely
    ("GCF_000020.1", "GCF_000020.1", _tax("Piscirickettsiaceae", "Piscirickettsia", "salmonis")),
]


@pytest.fixture()
def gtdb_tsv(tmp_path) -> Path:
    path = tmp_path / "bac120_metadata_r232.tsv.gz"
    header = "accession\tgtdb_genome_representative\tgtdb_taxonomy\tncbi_genbank_assembly_accession"
    lines = [header]
    for acc, rep, tax in _ROWS:
        lines.append(f"RS_{acc}\tRS_{rep}\t{tax}\t{acc}")
    with gzip.open(path, "wt", encoding="utf-8") as fo:
        fo.write("\n".join(lines) + "\n")
    return path


def _params(gtdb_tsv: Path, **overrides) -> MetadataParams:
    defaults = dict(
        dataset="all", level="species", release="232.0", version="bac120",
        target_genus="Francisella", target_species="tularensis",
        metadata_path=str(gtdb_tsv),
    )
    defaults.update(overrides)
    return MetadataParams(**defaults)


def _read_selection(workdir: Path) -> list[dict]:
    with open(workdir / SELECTION_TSV, encoding="utf-8", newline="") as fo:
        return list(csv.DictReader(fo, delimiter="\t"))


def test_tsv_species_selection_end_to_end(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    count = metadata.run(ctx, _params(gtdb_tsv))
    assert count == 3

    rows = _read_selection(ctx.workdir)
    by_acc = {r["accession"]: r for r in rows}
    assert set(by_acc) == {
        "GCF_000001.1", "GCF_000002.1", "GCF_000003.1", "GCF_000010.1",
    }
    # the same-genus representative is the outgroup
    outgroups = [r for r in rows if r["is_outgroup"] in ("1", "True", "true")]
    assert [r["accession"] for r in outgroups] == ["GCF_000010.1"]
    assert (ctx.workdir / "outgroup_accession.txt").read_text(
        encoding="utf-8"
    ).strip() == "GCF_000010.1"
    assert "metadata" in ctx.config.stages


def test_tsv_rep_dataset_selects_representatives_only(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    count = metadata.run(ctx, _params(gtdb_tsv, dataset="rep"))
    # GCF_000002.1 is not its own representative
    assert count == 2
    accs = {r["accession"] for r in _read_selection(ctx.workdir)}
    assert "GCF_000002.1" not in accs


def test_tsv_limit_caps_selection(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    assert metadata.run(ctx, _params(gtdb_tsv, limit=1)) == 1


def test_tsv_unknown_target_raises(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="not found"):
        metadata.run(ctx, _params(gtdb_tsv, target_genus="Yersinia", target_species="pestis"))


def test_tsv_explicit_outgroup_accession(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(gtdb_tsv, outgroup_accession="GCF_000020.1"))
    assert (ctx.workdir / "outgroup_accession.txt").read_text(
        encoding="utf-8"
    ).strip() == "GCF_000020.1"


def test_tsv_unknown_outgroup_accession_raises(tmp_path, gtdb_tsv) -> None:
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    with pytest.raises(UserInputError, match="Outgroup accession"):
        metadata.run(ctx, _params(gtdb_tsv, outgroup_accession="GCF_999999.9"))


def test_obtain_metadata_downloads_and_verifies(tmp_path, gtdb_tsv, monkeypatch) -> None:
    """Without --metadata-path the stage downloads the modern .tsv.gz layout."""
    calls: dict[str, list] = {"download": [], "verify": []}

    def fake_download(url, dest, *, logger, **kw):
        calls["download"].append(url)
        Path(dest).write_bytes(gtdb_tsv.read_bytes())

    def fake_verify(dest, manifest_url, *, logger, **kw):
        calls["verify"].append(manifest_url)

    monkeypatch.setattr(metadata.http, "download", fake_download)
    monkeypatch.setattr(metadata.http, "verify_md5_manifest", fake_verify)

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    count = metadata.run(ctx, _params(gtdb_tsv, metadata_path=None))
    assert count == 3
    assert calls["download"] and calls["download"][0].endswith(
        "232.0/bac120_metadata_r232.tsv.gz"
    )
    assert calls["verify"], "MD5 manifest verification must run on the download"


def test_obtain_metadata_falls_back_to_tarball_layout(tmp_path, gtdb_tsv, monkeypatch) -> None:
    def fake_download(url, dest, *, logger, **kw):
        if url.endswith(".tsv.gz"):
            raise WorkdirError("404")
        # legacy layout: a .tar.gz containing the .tsv
        import tarfile

        tsv = tmp_path / "member.tsv"
        with gzip.open(gtdb_tsv, "rb") as fo:
            tsv.write_bytes(fo.read())
        with tarfile.open(dest, "w:gz") as tar:
            tar.add(tsv, arcname="bac120_metadata_r232.tsv")

    monkeypatch.setattr(metadata.http, "download", fake_download)
    monkeypatch.setattr(metadata.http, "verify_md5_manifest", lambda *a, **k: None)

    ctx = WorkdirContext(tmp_path / "wd", create=True)
    assert metadata.run(ctx, _params(gtdb_tsv, metadata_path=None)) == 3


# --- API source ---------------------------------------------------------------


def _api_row(gid: str, genus: str, species: str, family="f__Francisellaceae", is_rep=True):
    return {
        "gid": gid,
        "gtdbFamily": family,
        "gtdbGenus": f"g__{genus}",
        "gtdbSpecies": f"s__{genus} {species}",
        "gtdbIsRep": is_rep,
    }


def test_api_species_selection_end_to_end(tmp_path, monkeypatch) -> None:
    selection_rows = [
        _api_row("GCF_000001.1", "Francisella", "tularensis"),
        _api_row("GCF_000002.1", "Francisella", "tularensis", is_rep=False),
    ]
    parent_rows = selection_rows + [
        _api_row("GCF_000010.1", "Francisella", "philomiragia"),
    ]

    def fake_api_get(path, params=None):
        if "/taxon/" in path and path.endswith("/genomes-detail"):
            if "tularensis" in path:
                return {"rows": selection_rows}
            return {"rows": parent_rows}
        raise AssertionError(f"unexpected API path {path}")

    monkeypatch.setattr(metadata, "_api_get", fake_api_get)
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(
        dataset="all", level="species", source="api",
        target_genus="Francisella", target_species="tularensis",
    )
    count = metadata.run(ctx, params)
    assert count == 2
    rows = _read_selection(ctx.workdir)
    outgroups = [r["accession"] for r in rows if r["is_outgroup"] in ("1", "True", "true")]
    assert outgroups == ["GCF_000010.1"]


def test_api_no_rows_raises(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(metadata, "_api_get", lambda path, params=None: {"rows": []})
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    params = MetadataParams(
        dataset="all", level="genus", source="api", target_genus="Francisella"
    )
    with pytest.raises(UserInputError, match="no genomes"):
        metadata.run(ctx, params)


def test_narrower_reselection_shrinks_manifest(tmp_path, gtdb_tsv) -> None:
    """Re-running metadata with a narrower selection removes de-selected rows."""
    ctx = WorkdirContext(tmp_path / "wd", create=True)
    metadata.run(ctx, _params(gtdb_tsv))  # 3 tularensis + outgroup
    assert len(ctx.manifest.all_genomes(include_outgroup=True)) == 4

    metadata.run(ctx, _params(gtdb_tsv, limit=1))
    accs = {g.accession for g in ctx.manifest.all_genomes(include_outgroup=True)}
    assert len(accs) == 2  # one selected + one outgroup, older rows gone
