"""`repgenr versions` surfaces recorded tool versions from repgenr.yaml."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from repgenr.cli.main import app
from repgenr.core.config import Config

_runner = CliRunner()


def _workdir_with_versions(tmp_path: Path) -> Path:
    cfg = Config()
    cfg.record_stage("vmetadata", tool_versions={"datasets": "16.0.0"}, completed="t")
    cfg.record_stage("vgenome", tool_versions={"mashtree": "1.4.6"}, completed="t")
    cfg.save(tmp_path)
    return tmp_path


def test_versions_stdout(tmp_path: Path) -> None:
    wd = _workdir_with_versions(tmp_path)
    result = _runner.invoke(app, ["versions", "-wd", str(wd)])
    assert result.exit_code == 0
    assert "datasets: 16.0.0" in result.stdout
    assert "mashtree: 1.4.6" in result.stdout


def test_versions_fragment_file(tmp_path: Path) -> None:
    wd = _workdir_with_versions(tmp_path)
    out = tmp_path / "frag.yml"
    result = _runner.invoke(app, ["versions", "-wd", str(wd), "--versions-out", str(out)])
    assert result.exit_code == 0
    # 4-space-indented, sorted -> slots under a process key in versions.yml
    assert out.read_text() == "    datasets: 16.0.0\n    mashtree: 1.4.6\n"


def test_versions_empty_workdir(tmp_path: Path) -> None:
    Config().save(tmp_path)  # no stages recorded
    out = tmp_path / "frag.yml"
    result = _runner.invoke(app, ["versions", "-wd", str(tmp_path), "--versions-out", str(out)])
    assert result.exit_code == 0
    assert out.read_text() == ""  # nothing recorded -> empty fragment


def test_versions_missing_workdir_exits_3(tmp_path: Path) -> None:
    """A path with no repgenr.yaml is a wrong -wd, not an empty run: exit 3 and
    write no fragment, so a Nextflow module cannot publish an empty versions.yml."""
    out = tmp_path / "frag.yml"
    result = _runner.invoke(
        app, ["versions", "-wd", str(tmp_path / "missing"), "--versions-out", str(out)]
    )
    assert result.exit_code == 3, result.output
    assert "repgenr.yaml" in result.output
    assert not out.exists()


def _metadata_workdir(tmp_path: Path, **params) -> Path:
    tmp_path.mkdir(exist_ok=True)
    cfg = Config()
    cfg.record_stage("metadata", tool=params.pop("tool"), params=params, completed="t")
    cfg.save(tmp_path)
    return tmp_path


def test_versions_names_the_gtdb_release_of_the_table(tmp_path: Path) -> None:
    wd = _metadata_workdir(
        tmp_path, tool="gtdb-table", source="tsv", release="232.0", api_query_date=None
    )
    out = tmp_path / "frag.yml"
    result = _runner.invoke(app, ["versions", "-wd", str(wd), "--versions-out", str(out)])
    assert result.exit_code == 0, result.output
    assert out.read_text() == "    gtdb_release: 232.0\n"


def test_versions_names_the_gtdb_api_query_date(tmp_path: Path) -> None:
    """The API path recorded release: null and nothing else to date the taxonomy."""
    wd = _metadata_workdir(
        tmp_path,
        tool="gtdb-api",
        source="api",
        release=None,
        api_query_date="2026-10-08T07:30:12+00:00",
    )
    result = _runner.invoke(app, ["versions", "-wd", str(wd)])
    assert result.exit_code == 0, result.output
    assert "gtdb_api_query_date: 2026-10-08T07:30:12+00:00" in result.stdout
    assert "gtdb_release" not in result.stdout
    # Quoted in the fragment, so a YAML 1.1 loader keeps the string and does
    # not turn it into a timestamp.
    out = tmp_path / "frag.yml"
    _runner.invoke(app, ["versions", "-wd", str(wd), "--versions-out", str(out)])
    assert out.read_text() == '    gtdb_api_query_date: "2026-10-08T07:30:12+00:00"\n'
    import yaml

    loaded = yaml.safe_load(out.read_text())
    assert loaded == {"gtdb_api_query_date": "2026-10-08T07:30:12+00:00"}


def test_status_shows_the_gtdb_release_or_api_query_date(tmp_path: Path) -> None:
    table = _metadata_workdir(
        tmp_path / "t", tool="gtdb-table", source="tsv", release="232.0", api_query_date=None
    )
    api = _metadata_workdir(
        tmp_path / "a",
        tool="gtdb-api",
        source="api",
        release=None,
        api_query_date="2026-10-08T07:30:12+00:00",
    )
    out_table = _runner.invoke(app, ["status", "-wd", str(table)]).output
    out_api = _runner.invoke(app, ["status", "-wd", str(api)]).output
    assert "metadata [gtdb-table]  t  (GTDB release 232.0)" in out_table
    assert "metadata [gtdb-api]  t  (GTDB API queried 2026-10-08T07:30:12+00:00)" in out_api


def test_a_tool_with_different_versions_per_stage_keeps_each(tmp_path: Path) -> None:
    # dereplicate ran sourmash in its image, phylo ran the host sourmash: the
    # merged output kept only the phylo value and lost the dereplicate one.
    cfg = Config()
    image = "quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0"
    cfg.record_stage("dereplicate", tool_versions={"sourmash": image}, completed="t")
    cfg.record_stage(
        "phylo", tool_versions={"sourmash": "4.9.4", "iqtree": "3.1.2"}, completed="t"
    )
    cfg.record_stage("glance", tool_versions={"iqtree": "3.1.2"}, completed="t")
    cfg.save(tmp_path)
    result = _runner.invoke(app, ["versions", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert result.stdout.splitlines() == [
        "iqtree: 3.1.2",
        f"sourmash (dereplicate): {image}",
        "sourmash (phylo): 4.9.4",
    ]
    out = tmp_path / "frag.yml"
    _runner.invoke(app, ["versions", "-wd", str(tmp_path), "--versions-out", str(out)])
    import yaml

    assert yaml.safe_load("p:\n" + out.read_text())["p"] == {
        "iqtree": "3.1.2",
        "sourmash (dereplicate)": image,
        "sourmash (phylo)": "4.9.4",
    }


def test_versions_of_an_unfinished_stage_are_flagged(tmp_path: Path) -> None:
    cfg = Config()
    cfg.record_stage("dereplicate", tool="skder", tool_versions={"skder": "1.3.6"})
    cfg.save(tmp_path)
    result = _runner.invoke(app, ["versions", "-wd", str(tmp_path)])
    assert result.exit_code == 0
    assert result.stdout == "skder: 1.3.6\n"
    assert "dereplicate did not finish" in result.stderr


def test_a_fragment_value_that_is_not_plain_yaml_is_quoted(tmp_path: Path) -> None:
    # An unnumbered version line is recorded verbatim and may hold ": ".
    from repgenr.core.versions import write_versions_fragment

    out = tmp_path / "frag.yml"
    write_versions_fragment(out, {"tool": "build: abc #1", "other": "1.2.3"})
    import yaml

    assert yaml.safe_load("p:\n" + out.read_text())["p"] == {
        "tool": "build: abc #1",
        "other": "1.2.3",
    }
    assert "    other: 1.2.3\n" in out.read_text()
