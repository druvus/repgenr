"""derep_unpack stage test (no external tools)."""

from __future__ import annotations

from pathlib import Path

from repgenr.core.context import WorkdirContext
from repgenr.stages.derep_unpack import DerepUnpackParams, run


def test_unpack_clusters(workdir: Path) -> None:
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for name in ("a.fasta", "b.fasta", "c.fasta"):
        (gdir / name).write_text(">x\nACGT\n")

    derep = workdir / "derep"
    derep.mkdir()
    (derep / "clusters.tsv").write_text(
        "representative\tmember\na.fasta\ta.fasta\na.fasta\tb.fasta\nc.fasta\tc.fasta\n"
    )

    ctx = WorkdirContext(workdir, create=True)
    unpack = run(ctx, DerepUnpackParams())

    a_dir = unpack / "a"
    assert sorted(p.name for p in a_dir.iterdir()) == ["a.fasta", "b.fasta"]
    # cluster c has only the representative -> still emitted with representative
    assert (unpack / "c" / "c.fasta").exists()


def test_unpack_no_representant(workdir: Path) -> None:
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for name in ("a.fasta", "b.fasta"):
        (gdir / name).write_text(">x\nACGT\n")
    derep = workdir / "derep"
    derep.mkdir()
    (derep / "clusters.tsv").write_text(
        "representative\tmember\na.fasta\ta.fasta\na.fasta\tb.fasta\n"
    )
    ctx = WorkdirContext(workdir, create=True)
    unpack = run(ctx, DerepUnpackParams(no_representant=True))
    assert [p.name for p in (unpack / "a").iterdir()] == ["b.fasta"]


def _unpack_with_missing(workdir: Path, n_missing: int, caplog) -> list[str]:
    import logging

    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    (gdir / "a.fasta").write_text(">x\nACGT\n")
    missing = [f"m{i:02d}.fasta" for i in range(n_missing)]
    rows = "".join(f"a.fasta\t{m}\n" for m in missing)
    derep = workdir / "derep"
    derep.mkdir()
    (derep / "clusters.tsv").write_text(f"representative\tmember\na.fasta\ta.fasta\n{rows}")
    ctx = WorkdirContext(workdir, create=True)
    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.WARNING):
        run(ctx, DerepUnpackParams())
    return [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]


def test_a_member_missing_from_genomes_is_named(workdir: Path, caplog) -> None:
    warnings = _unpack_with_missing(workdir, 2, caplog)
    assert len(warnings) == 2
    assert "m00.fasta" in warnings[0] and "genomes" in warnings[0]
    assert "m01.fasta" in warnings[1]


def test_many_missing_members_are_listed_on_one_line(workdir: Path, caplog) -> None:
    warnings = _unpack_with_missing(workdir, 11, caplog)
    assert len(warnings) == 1
    assert "11" in warnings[0] and "m00.fasta" in warnings[0] and "m10.fasta" in warnings[0]


def test_representatives_sharing_a_stem_get_distinct_directories(workdir: Path) -> None:
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for name in (
        "dup.fasta",
        "dup.fna",
        "m.fasta",
        "a.b.c.fasta",
        "z.fasta.gz",
        "Up.fa",
        "up.fna",
        "n.fna.gz",
        "f.fa.gz",
    ):
        (gdir / name).write_text(f">{name}\nACGT\n")
    derep = workdir / "derep"
    derep.mkdir()
    (derep / "clusters.tsv").write_text(
        "representative\tmember\n"
        "dup.fasta\tdup.fasta\ndup.fasta\tm.fasta\n"
        "dup.fna\tdup.fna\n"
        "a.b.c.fasta\ta.b.c.fasta\n"
        "z.fasta.gz\tz.fasta.gz\n"
        "Up.fa\tUp.fa\n"
        "up.fna\tup.fna\n"
        "n.fna.gz\tn.fna.gz\n"
        "f.fa.gz\tf.fa.gz\n"
    )
    ctx = WorkdirContext(workdir, create=True)
    unpack = run(ctx, DerepUnpackParams())

    # Names that collide (also ignoring case) fall back to the full file name;
    # others drop the genome extension, the gzip suffixes included.
    assert sorted(p.name for p in unpack.iterdir()) == [
        "Up.fa",
        "a.b.c",
        "dup.fasta",
        "dup.fna",
        "f",
        "n",
        "up.fna",
        "z",
    ]
    assert sorted(p.name for p in (unpack / "dup.fasta").iterdir()) == ["dup.fasta", "m.fasta"]
    assert [p.name for p in (unpack / "dup.fna").iterdir()] == ["dup.fna"]


def test_a_failed_run_keeps_the_previous_unpacked_tree(workdir: Path, monkeypatch) -> None:
    import pytest

    from repgenr.stages import derep_unpack

    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for name in ("a.fasta", "b.fasta", "c.fasta"):
        (gdir / name).write_text(">x\nACGT\n")
    derep = workdir / "derep"
    derep.mkdir()
    (derep / "clusters.tsv").write_text(
        "representative\tmember\na.fasta\ta.fasta\na.fasta\tb.fasta\nc.fasta\tc.fasta\n"
    )
    ctx = WorkdirContext(workdir, create=True)
    unpack = run(ctx, DerepUnpackParams())
    before = sorted(str(p.relative_to(unpack)) for p in unpack.rglob("*"))

    calls = {"n": 0}
    real = derep_unpack.link_or_copy

    def fail_on_second(src, dst):  # noqa: ANN001
        calls["n"] += 1
        if calls["n"] == 2:
            raise OSError("disk full")
        return real(src, dst)

    monkeypatch.setattr(derep_unpack, "link_or_copy", fail_on_second)
    with pytest.raises(OSError):
        run(ctx, DerepUnpackParams())

    assert sorted(str(p.relative_to(unpack)) for p in unpack.rglob("*")) == before
    assert sorted(p.name for p in derep.iterdir()) == ["clusters.tsv", "unpacked"]


def _three_genome_workdir(workdir: Path) -> WorkdirContext:
    gdir = workdir / "genomes"
    gdir.mkdir(parents=True)
    for name in ("a.fasta", "b.fasta", "c.fasta"):
        (gdir / name).write_text(">x\nACGT\n")
    derep = workdir / "derep"
    derep.mkdir()
    (derep / "clusters.tsv").write_text(
        "representative\tmember\na.fasta\ta.fasta\na.fasta\tb.fasta\nc.fasta\tc.fasta\n"
    )
    return WorkdirContext(workdir, create=True)


def _info_messages(ctx: WorkdirContext, caplog, params: DerepUnpackParams) -> list[str]:
    import logging

    ctx.logger.addHandler(caplog.handler)
    with caplog.at_level(logging.INFO):
        run(ctx, params)
    return [r.getMessage() for r in caplog.records]


def test_summary_counts_directories_written_not_clusters_read(workdir: Path, caplog) -> None:
    ctx = _three_genome_workdir(workdir)
    messages = _info_messages(ctx, caplog, DerepUnpackParams(no_representant=True))
    # Cluster c holds only its representative, so no directory is written for it.
    assert any("1 genome files into 1 cluster directories" in m for m in messages)


def test_copying_instead_of_linking_is_logged(workdir: Path, caplog, monkeypatch) -> None:
    import os

    def no_link(src, dst):  # noqa: ANN001
        raise OSError("hard links not supported")

    monkeypatch.setattr(os, "link", no_link)
    ctx = _three_genome_workdir(workdir)
    messages = _info_messages(ctx, caplog, DerepUnpackParams())
    assert sum("copying files" in m for m in messages) == 1
    assert (ctx.derep_dir / "unpacked" / "a" / "b.fasta").read_text() == ">x\nACGT\n"


def test_linking_logs_no_copy_notice(workdir: Path, caplog) -> None:
    ctx = _three_genome_workdir(workdir)
    messages = _info_messages(ctx, caplog, DerepUnpackParams())
    assert not any("copying files" in m for m in messages)
    assert any("3 genome files into 2 cluster directories" in m for m in messages)
