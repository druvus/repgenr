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
