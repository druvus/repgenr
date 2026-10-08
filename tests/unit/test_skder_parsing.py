"""Offline tests for skDER output parsing + edge-based membership."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.dereplicators.skder import _parse_skder_output

_LOG = logging.getLogger("test")


def _make_skder_out(tmp_path: Path, edges: str) -> Path:
    out = tmp_path / "skder_out"
    reps = out / "Dereplicated_Representative_Genomes"
    reps.mkdir(parents=True)
    for name in ("repA.fasta", "repB.fasta"):
        (reps / name).write_text(">x\nACGT\n")
    (out / "Skani_Triangle_Edge_Output.txt").write_text(edges)
    return out


def test_membership_from_edges(tmp_path: Path) -> None:
    # memX is closest (highest ANI) to repA; memY only similar to repB
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/repA.fasta\t/g/memX.fasta\t99.5\t90\t88\trA\tmX\n"
        "/g/repB.fasta\t/g/memX.fasta\t97.0\t80\t80\trB\tmX\n"
        "/g/memY.fasta\t/g/repB.fasta\t99.2\t85\t85\tmY\trB\n"
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "memX", "memY")]

    result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=_LOG)

    assert {p.name for p in result.representatives} == {"repA.fasta", "repB.fasta"}
    # memX assigned to repA (99.5 > 97.0); memY assigned to repB
    assert "memX.fasta" in result.clusters["repA.fasta"]
    assert "memY.fasta" in result.clusters["repB.fasta"]
    assert result.genome_status["memX.fasta"] == "contained"
    assert result.genome_status["repA.fasta"] == "representative"


def test_edges_below_cutoff_not_clustered(tmp_path: Path) -> None:
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/repA.fasta\t/g/memX.fasta\t90.0\t90\t90\trA\tmX\n"  # ANI below 99 cutoff
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "memX")]
    result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=_LOG)
    # memX not assigned to any cluster, but still marked contained
    assert result.clusters["repA.fasta"] == []
    assert result.genome_status["memX.fasta"] == "contained"


def test_skder_warns_when_argv_may_overflow(tmp_path, monkeypatch, caplog) -> None:
    import logging

    import pytest

    from repgenr.dereplicators import skder as skder_mod
    from repgenr.dereplicators.base import DerepParams
    from repgenr.dereplicators.skder import SkderDereplicator

    monkeypatch.setattr(skder_mod, "_ARGV_WARN_GENOMES", 3)

    def stop(caps, cmd, **kwargs):
        raise RuntimeError("stop before running the tool")

    monkeypatch.setattr(skder_mod, "run_tool", stop)
    genomes = [tmp_path / f"g{i}.fasta" for i in range(4)]
    params = DerepParams(primary_ani=0.9, secondary_ani=0.99, threads=1)
    with caplog.at_level(logging.WARNING), pytest.raises(RuntimeError):
        SkderDereplicator().dereplicate(
            genomes, tmp_path / "out", params, logging.getLogger("test")
        )
    assert "--process-size" in caplog.text


def test_partial_genome_joins_by_its_own_aligned_fraction(tmp_path: Path) -> None:
    # A 40 percent fragment of a repA-like genome: 99 percent of the fragment
    # aligns, 40 percent of repA does. skDER covers the fragment by repA, so the
    # cutoff applies to the member's aligned fraction, not the smaller of the two.
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/frag.fasta\t/g/repA.fasta\t99.6\t99.3\t39.7\tf\trA\n"
        "/g/repB.fasta\t/g/frag.fasta\t95.0\t39.7\t99.2\trB\tf\n"
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "frag")]
    result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=_LOG)
    assert result.clusters["repA.fasta"] == ["frag.fasta"]
    assert result.genome_status["frag.fasta"] == "contained"


def test_representative_aligned_fraction_does_not_admit_a_member(tmp_path: Path) -> None:
    # The member's own aligned fraction (30) is below the cutoff; only the
    # representative's side (95) passes. The edge does not assign the member by
    # the cutoff rule, so it is placed under its closest representative at or
    # above the ANI cutoff, with a warning naming it.
    edges = (
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "/g/repA.fasta\t/g/memX.fasta\t99.5\t95\t30\trA\tmX\n"
    )
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path(f"/g/{n}.fasta") for n in ("repA", "repB", "memX")]
    log = logging.getLogger("test.skder.nearest")
    records: list[logging.LogRecord] = []
    handler = logging.Handler()
    handler.emit = records.append  # type: ignore[method-assign]
    log.addHandler(handler)
    try:
        result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=log)
    finally:
        log.removeHandler(handler)
    assert result.clusters["repA.fasta"] == ["memX.fasta"]
    assert any("memX.fasta" in r.getMessage() for r in records if r.levelno == logging.WARNING)


def test_skder_refuses_ani_below_80_percent_before_running(tmp_path, monkeypatch) -> None:
    import pytest

    from repgenr.core.errors import UserInputError
    from repgenr.dereplicators import skder as skder_mod
    from repgenr.dereplicators.base import DerepParams
    from repgenr.dereplicators.skder import SkderDereplicator

    def must_not_run(caps, cmd, **kwargs):
        raise AssertionError("skDER was started")

    monkeypatch.setattr(skder_mod, "run_tool", must_not_run)
    params = DerepParams(secondary_ani=0.78, threads=1)
    with pytest.raises(UserInputError, match="80 percent"):
        SkderDereplicator().dereplicate([tmp_path / "g.fasta"], tmp_path / "out", params, _LOG)


def test_many_placements_are_one_line_capped_on_the_console(tmp_path: Path) -> None:
    # Twelve members pass only the ANI cutoff: one warning, whose console form
    # names the first few and counts the rest, while the full text lists all.
    rows = "".join(f"/g/repA.fasta\t/g/mem{i:02d}.fasta\t99.5\t95\t30\trA\tm\n" for i in range(12))
    edges = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tR\tQ\n" + rows
    out = _make_skder_out(tmp_path, edges)
    genomes = [Path("/g/repA.fasta"), Path("/g/repB.fasta")]
    genomes += [Path(f"/g/mem{i:02d}.fasta") for i in range(12)]
    log = logging.getLogger("test.skder.many")
    records: list[logging.LogRecord] = []
    handler = logging.Handler()
    handler.emit = records.append  # type: ignore[method-assign]
    log.addHandler(handler)
    try:
        result = _parse_skder_output(out, genomes, ani_cutoff=99.0, af_cutoff=50.0, logger=log)
    finally:
        log.removeHandler(handler)
    assert len(result.clusters["repA.fasta"]) == 12
    warnings = [r for r in records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert all(f"mem{i:02d}.fasta" in warnings[0].getMessage() for i in range(12))
    console = warnings[0].repgenr_console  # type: ignore[attr-defined]
    assert "and 7 more" in console and "mem11.fasta" not in console
