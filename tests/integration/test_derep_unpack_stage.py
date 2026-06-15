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
        "representative\tmember\n"
        "a.fasta\ta.fasta\n"
        "a.fasta\tb.fasta\n"
        "c.fasta\tc.fasta\n"
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
