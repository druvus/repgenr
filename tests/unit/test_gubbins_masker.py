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
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n", encoding="utf-8")
    masker = mod.GubbinsMasker()
    out = masker.mask(full, tmp_path / "gub", MaskParams(threads=4), logging.getLogger("t"))
    assert out.exists()
    argv = calls[0]
    assert argv[0] == "run_gubbins.py"
    assert argv[argv.index("--threads") + 1] == "4"
    assert argv[-1] == str(full)


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
