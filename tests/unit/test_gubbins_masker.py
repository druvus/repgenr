"""Gubbins is invoked on the whole-genome alignment with a thread budget."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.maskers import gubbins as mod
from repgenr.maskers.base import MaskParams


def test_gubbins_argv(tmp_path: Path, monkeypatch) -> None:
    calls: list[list] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        out_prefix = str(kw["cwd"] / "gubbins") + ".filtered_polymorphic_sites.fasta"
        Path(out_prefix).write_text(">a\nA\n", encoding="utf-8")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda: True)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n", encoding="utf-8")
    masker = mod.GubbinsMasker()
    out = masker.mask(full, tmp_path / "gub", MaskParams(threads=4), logging.getLogger("t"))
    assert out.exists()
    argv = calls[0]
    assert argv[0] == "run_gubbins.py"
    assert argv[argv.index("--threads") + 1] == "4"
    assert "--tree-builder" not in argv, "Gubbins' own default stands when RAxML can run"
    assert argv[-1] == str(full)


def test_extras_select_gubbins_tree_builders(tmp_path: Path, monkeypatch) -> None:
    calls: list[list] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        Path(str(kw["cwd"] / "gubbins") + ".filtered_polymorphic_sites.fasta").write_text(
            ">a\nA\n", encoding="utf-8"
        )

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda: False)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n", encoding="utf-8")
    params = MaskParams(
        threads=8,
        extra={
            "gubbins_tree_builder": "fasttree",
            "gubbins_first_tree_builder": "rapidnj",
            "gubbins_args": "--min-snps 5 --iterations 3",
        },
    )
    mod.GubbinsMasker().mask(full, tmp_path / "gub", params, logging.getLogger("t"))
    argv = calls[0]
    assert argv[argv.index("--tree-builder") + 1] == "fasttree"
    assert argv[argv.index("--first-tree-builder") + 1] == "rapidnj"
    assert argv[argv.index("--min-snps") + 1] == "5"
    assert argv[argv.index("--iterations") + 1] == "3"
    assert argv[argv.index("--threads") + 1] == "8", "a chosen builder keeps the thread budget"
    assert argv[-2] == "--prefix" or argv[-3] == "--prefix"


def test_resolve_tree_builder_falls_back_without_threaded_raxml(monkeypatch, caplog) -> None:
    """Gubbins exits when asked for threads without a PTHREADS RAxML build."""
    logger = logging.getLogger("t")
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda: False)

    monkeypatch.setattr(mod.shutil, "which", lambda name: "/bin/x" if name == "iqtree2" else None)
    with caplog.at_level(logging.WARNING):
        assert mod.resolve_tree_builder(None, 8, logger, on_host=True) == ("iqtree", 8)
    assert "iqtree" in caplog.text

    monkeypatch.setattr(mod.shutil, "which", lambda name: None)
    assert mod.resolve_tree_builder(None, 8, logger, on_host=True) == (None, 1)

    # One thread, a container run, or an explicit choice never trigger it.
    assert mod.resolve_tree_builder(None, 1, logger, on_host=True) == (None, 1)
    assert mod.resolve_tree_builder(None, 8, logger, on_host=False) == (None, 8)
    assert mod.resolve_tree_builder("raxmlng", 8, logger, on_host=True) == ("raxmlng", 8)
    monkeypatch.setattr(mod, "multithreaded_raxml_available", lambda: True)
    assert mod.resolve_tree_builder(None, 8, logger, on_host=True) == (None, 8)


def test_sanitise_replaces_iupac_with_n(tmp_path, caplog) -> None:
    """IUPAC codes from a diploid-style consensus become N before Gubbins."""
    import logging

    from repgenr.maskers.gubbins import sanitise_alignment

    src = tmp_path / "aln.fasta"
    src.write_text(">a\nACGTRYACGT\n>b\nACGTACGT-N\n", encoding="utf-8")
    with caplog.at_level(logging.WARNING):
        out = sanitise_alignment(src, tmp_path / "clean.fasta", logging.getLogger("t"))
    assert out.read_text(encoding="utf-8") == ">a\nACGTNNACGT\n>b\nACGTACGT-N\n"
    assert "2 ambiguous base(s)" in caplog.text


def test_sanitise_keeps_a_clean_alignment(tmp_path) -> None:
    import logging

    from repgenr.maskers.gubbins import sanitise_alignment

    src = tmp_path / "aln.fasta"
    src.write_text(">a\nACGT\n", encoding="utf-8")
    assert sanitise_alignment(src, tmp_path / "clean.fasta", logging.getLogger("t")) == src
    assert not (tmp_path / "clean.fasta").exists()


def test_exclude_runs_gubbins_on_the_ingroup_and_masks_everyone(
    tmp_path: Path, monkeypatch
) -> None:
    """D-10: the outgroup stays out of the scan but keeps its place in the output."""
    calls: list[list] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        calls.append([str(a) for a in argv])
        prefix = str(kw["cwd"] / "gubbins")
        # Gubbins predicts a recombinant block in `a` over columns 2-3.
        Path(prefix + ".recombination_predictions.gff").write_text(
            "##gff-version 3\n"
            'SEQUENCE\tGUBBINS\tCDS\t2\t3\t0.0\t.\t.\tnode="a";taxa="  a";snp_count="2"\n',
            encoding="utf-8",
        )

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGTA\n>b\nATTTA\n>c\nATTTA\n>og\nGGGGG\n", encoding="utf-8")
    out = mod.GubbinsMasker().mask(
        full,
        tmp_path / "gub",
        MaskParams(threads=2, exclude=frozenset({"og"})),
        logging.getLogger("t"),
    )
    argv = calls[0]
    ingroup = mod.read_fasta(Path(argv[-1]))
    assert set(ingroup) == {"a", "b", "c"}, "the outgroup is not handed to Gubbins"
    masked = mod.read_fasta(out)
    assert set(masked) == {"a", "b", "c", "og"}, "the outgroup is in the masked alignment"
    # Column 1 (A/G) and column 5 (A/G) vary; columns 2-3 of `a` are N so only
    # b/c/og decide those: T vs G still varies. Column 4: T/T/T/G varies.
    assert masked["og"] == "GGGGG"
    assert masked["a"].startswith("A") and "N" in masked["a"]


def test_gff_regions_and_polymorphic_sites(tmp_path: Path) -> None:
    gff = tmp_path / "p.gff"
    gff.write_text(
        'SEQ\tGUBBINS\tCDS\t1\t2\t0\t.\t.\tnode="x";taxa="  a  b";snp_count="1"\n'
        'SEQ\tGUBBINS\tCDS\t4\t4\t0\t.\t.\tnode="y";taxa="  b";snp_count="1"\n',
        encoding="utf-8",
    )
    regions = mod.read_recombination_gff(gff)
    assert regions == {"a": [(1, 2)], "b": [(1, 2), (4, 4)]}
    masked = mod.apply_masks({"a": "ACGT", "b": "ACGT", "c": "TCGA"}, regions)
    assert masked == {"a": "NNGT", "b": "NNGN", "c": "TCGA"}
    sites = mod.polymorphic_sites(masked)
    assert sites == {"a": "T", "b": "N", "c": "A"}, "only column 4 still varies among ACGT"
