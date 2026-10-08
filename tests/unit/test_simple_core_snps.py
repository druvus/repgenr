"""Tests for the simple SNP typer's core-SNP reduction."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from repgenr.snptypers.simple import _write_core_snps


def test_core_snp_reduction(tmp_path: Path) -> None:
    consensuses = {
        "ref": "ACGTACGT",
        "s1": "ACGAACGT",  # differs at col 3
        "s2": "ACGTACGA",  # differs at col 7
    }
    core = tmp_path / "core.fasta"
    matrix = tmp_path / "dist.tsv"
    n = _write_core_snps(consensuses, core, matrix)
    assert n == 2  # columns 3 and 7 are variable

    records = _read_fasta(core)
    assert records["ref"] == "TT"  # ref bases at the two variable columns
    assert records["s1"] == "AT"
    assert records["s2"] == "TA"

    # distance matrix: ref vs s1 = 1, ref vs s2 = 1, s1 vs s2 = 2
    lines = matrix.read_text().splitlines()
    assert lines[0].split("\t")[1:] == ["ref", "s1", "s2"]


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records


def test_thread_split_prefers_concurrent_genomes() -> None:
    """The budget buys workers first; threads only once genomes run out."""
    from repgenr.snptypers.simple import _split_threads

    assert _split_threads(8, 68) == (8, 1)
    assert _split_threads(8, 2) == (2, 4)
    assert _split_threads(1, 68) == (1, 1)
    assert _split_threads(8, 0) == (1, 1)
    assert _split_threads(0, 4) == (1, 1)


def test_per_genome_chain_is_threaded_and_clears_its_scratch(tmp_path: Path, monkeypatch) -> None:
    """Every tool gets the thread count, and the intermediates go once read."""
    from repgenr.snptypers import simple as mod

    calls: list[list[str]] = []

    def fake_run_tool(caps, cmd, **kw):  # noqa: ANN001
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        tool = cmd[0] if cmd[0] != "bcftools" else f"bcftools {cmd[1]}"
        out = kw.get("stdout_path")
        if out is None and "-o" in cmd:
            out = cmd[cmd.index("-o") + 1]
        if tool == "minimap2":
            Path(out).write_text(_sam(("g1", 0, "ref", 1, "4M")), encoding="utf-8")
        elif tool == "bcftools consensus":
            Path(out).write_text(">ref\nACGT\n", encoding="utf-8")
        elif tool == "bcftools index":
            Path(cmd[-1] + ".csi").write_text("", encoding="utf-8")
        elif tool == "samtools index":
            Path(cmd[-1] + ".bai").write_text("", encoding="utf-8")
        elif out is not None:
            Path(out).write_text("", encoding="utf-8")

    def fake_run_chain(caps, steps, *, logger, **kwargs):
        for prefix, command in steps:
            fake_run_tool(caps, command, log_prefix=prefix)

    monkeypatch.setattr(mod, "run_chain", fake_run_chain)
    work = tmp_path / "per_genome"
    work.mkdir()
    ref = tmp_path / "reference.fasta"
    ref.write_text(">ref\nACGT\n", encoding="utf-8")
    genome = tmp_path / "g1.fasta"
    genome.write_text(">g1\nACGA\n", encoding="utf-8")

    from repgenr.snptypers.base import SnpParams

    seq = mod._call_one(genome, ref, work, 4, SnpParams(), logging.getLogger("t"))
    assert seq == "ACGT"
    for cmd in calls:
        if cmd[0] == "minimap2":
            assert cmd[cmd.index("-t") + 1] == "4"
        if cmd[0] == "samtools":
            assert cmd[cmd.index("-@") + 1] == "4"
    pileup = [c for c in calls if c[:2] == ["bcftools", "mpileup"]][0]
    assert "-Ob" in pileup, "the pileup is written compressed"
    assert list(work.iterdir()) == [], "nothing is left behind once the consensus is read"


def test_columns_are_found_across_block_boundaries(tmp_path: Path, monkeypatch) -> None:
    """The scan reads the alignment in column blocks; a site must not fall between two."""
    from repgenr.snptypers import simple as mod

    monkeypatch.setattr(mod, "_COLUMN_BLOCK", 4)
    consensuses = {
        "ref": "AAAAAAAAAAAA",
        "s1": "AAATAAAATAAA",  # last column of block 1, first of block 3
    }
    n = mod._write_core_snps(consensuses, tmp_path / "core.fasta", tmp_path / "dist.tsv")
    assert n == 2
    assert _read_fasta(tmp_path / "core.fasta")["s1"] == "TT"
    assert (tmp_path / "dist.tsv").read_text().splitlines()[1].split("\t")[1:] == ["0", "2"]


def test_n_against_one_base_is_not_a_variable_column(tmp_path: Path) -> None:
    """A column with N (or a gap) in some genomes and one base in the rest is not variable."""
    from repgenr.snptypers.simple import _write_core_snps

    consensuses = {"ref": "ACGT", "s1": "ANGT", "s2": "A-GN"}
    assert _write_core_snps(consensuses, tmp_path / "c.fasta", tmp_path / "d.tsv") == 0
    # No core site, so no pair shares one: every distance is NA.
    lines = (tmp_path / "d.tsv").read_text().splitlines()
    assert [line.split("\t")[1:] for line in lines[1:]] == [["NA"] * 3] * 3


def test_a_real_snp_still_counts_beside_missing_data(tmp_path: Path) -> None:
    """Two bases in a column make it variable; N and lower case are handled."""
    from repgenr.snptypers.simple import _write_core_snps

    consensuses = {"ref": "ACGT", "s1": "ATGN", "s2": "aNgt"}
    n = _write_core_snps(consensuses, tmp_path / "core.fasta", tmp_path / "dist.tsv")
    assert n == 1
    assert _read_fasta(tmp_path / "core.fasta") == {"ref": "C", "s1": "T", "s2": "N"}


def test_distances_count_only_sites_where_both_genomes_have_a_base(tmp_path: Path) -> None:
    from repgenr.snptypers.simple import _write_core_snps

    consensuses = {
        "ref": "AAAA",
        "s1": "TTTA",
        "s2": "NTTT",  # no base where s1 differs from ref at column 0
    }
    n = _write_core_snps(consensuses, tmp_path / "core.fasta", tmp_path / "dist.tsv")
    assert n == 4
    rows = {
        line.split("\t")[0]: line.split("\t")[1:]
        for line in (tmp_path / "dist.tsv").read_text().splitlines()[1:]
    }
    assert rows["ref"] == ["0", "3", "3"]
    assert rows["s1"] == ["3", "0", "1"]  # column 0 is not shared with s2
    assert rows["s2"] == ["3", "1", "0"]


def test_ragged_consensuses_are_truncated_to_the_shortest(tmp_path: Path) -> None:
    from repgenr.snptypers.simple import _write_core_snps

    consensuses = {"ref": "ACGTACGT", "s1": "ACGA", "s2": "ACGT"}
    n = _write_core_snps(consensuses, tmp_path / "core.fasta", tmp_path / "dist.tsv")
    assert n == 1, "only the columns every genome has are compared"


def test_no_variable_sites_message_names_the_cause(tmp_path: Path, monkeypatch) -> None:
    """Identical genomes give a WorkdirError that names the cause and the options."""
    import pytest

    from repgenr.core.errors import WorkdirError
    from repgenr.snptypers import simple as mod
    from repgenr.snptypers.base import SnpParams

    ref = tmp_path / "refgenome.fasta"
    ref.write_text(">c1\nACGTACGT\n")
    genomes = [ref]
    for name in ("g1", "g2"):
        p = tmp_path / f"{name}.fasta"
        p.write_text(">c1\nACGTACGT\n")
        genomes.append(p)

    monkeypatch.setattr(mod, "run_tool", lambda *a, **k: None)
    monkeypatch.setattr(mod, "_call_one", lambda genome, *a, **k: "ACGTACGT")

    with pytest.raises(WorkdirError) as exc:
        mod.SimpleSnpTyper().call(
            genomes, ref, tmp_path / "out", SnpParams(threads=1), logging.getLogger("t")
        )
    msg = str(exc.value)
    assert "No variable sites" in msg
    assert "3 genomes" in msg
    assert "refgenome" in msg
    assert "closer reference" in msg
    assert "alignment-free" in msg


def _sam(*records: tuple[str, int, str, int, str]) -> str:
    """SAM text from (query, flag, reference, 1-based position, CIGAR) tuples."""
    lines = ["@HD\tVN:1.6"]
    for qname, flag, rname, pos, cigar in records:
        lines.append(f"{qname}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*")
    return "\n".join(lines) + "\n"


def test_covered_positions_follow_primary_and_supplementary_cigar_spans(tmp_path: Path) -> None:
    from repgenr.snptypers.simple import _covered_positions

    sam = tmp_path / "g.sam"
    sam.write_text(
        _sam(
            ("q", 0, "c1", 2, "2S3M1I2D1M5H"),  # 2..4 and 7 aligned; 5..6 deleted
            ("q", 2048, "c2", 1, "2="),  # supplementary on the second contig: 1..2
            ("q", 256, "c2", 3, "3M"),  # secondary: ignored
            ("q", 4, "*", 0, "*"),  # unmapped: ignored
            ("q", 16, "c2", 5, "1X"),  # reverse strand primary: 5
        ),
        encoding="utf-8",
    )
    covered = _covered_positions(sam, [("c1", 8), ("c2", 6)])
    assert "".join("1" if c else "0" for c in covered) == "01110010" + "110010"


def _fake_genome_chain(monkeypatch, sam_text: str, consensus: str, calls: list) -> None:
    from repgenr.snptypers import simple as mod

    def fake_run_chain(caps, steps, *, logger, **kwargs):
        for _prefix, command in steps:
            cmd = [str(c) for c in command]
            calls.append(cmd)
            out = cmd[cmd.index("-o") + 1] if "-o" in cmd else None
            if cmd[0] == "minimap2":
                Path(out).write_text(sam_text, encoding="utf-8")
            elif cmd[:2] == ["bcftools", "consensus"]:
                Path(out).write_text(f">c1\n{consensus}\n", encoding="utf-8")
            elif out is not None:
                Path(out).write_text("", encoding="utf-8")

    monkeypatch.setattr(mod, "run_chain", fake_run_chain)


def test_a_genome_missing_a_region_has_n_there_and_no_snps_against_its_copy(
    tmp_path: Path, monkeypatch
) -> None:
    """A copy with a region removed differs from the full genome at no site.

    Without masking the copy would keep the reference base over the removed
    region, and every SNP the full genome has there would separate the two.
    """
    from repgenr.snptypers import simple as mod
    from repgenr.snptypers.base import SnpParams

    ref = tmp_path / "reference.fasta"
    ref.write_text(">c1\nAAAAAAAAAA\n", encoding="utf-8")
    contigs = [("c1", 10)]
    full_consensus = "ATAAAAATAA"  # SNPs at columns 1 and 7
    log = logging.getLogger("t")
    work = tmp_path / "work"
    work.mkdir()

    _fake_genome_chain(monkeypatch, _sam(("full", 0, "c1", 1, "10M")), full_consensus, [])
    full = mod._call_one(tmp_path / "full.fasta", ref, work, 1, SnpParams(), log, contigs=contigs)
    # The copy lacks reference positions 6..10, so its alignment ends at 5 and
    # the consensus keeps the reference base there.
    _fake_genome_chain(monkeypatch, _sam(("cut", 0, "c1", 1, "5M")), "ATAAAAAAAA", [])
    cut = mod._call_one(tmp_path / "cut.fasta", ref, work, 1, SnpParams(), log, contigs=contigs)
    assert full == full_consensus
    assert cut == "ATAAANNNNN"

    n = mod._write_core_snps(
        {"reference": "AAAAAAAAAA", "full": full, "cut": cut},
        tmp_path / "core.fasta",
        tmp_path / "dist.tsv",
    )
    assert n == 2, "column 7 stays variable (reference against full); N adds no column"
    rows = {
        line.split("\t")[0]: line.split("\t")[1:]
        for line in (tmp_path / "dist.tsv").read_text().splitlines()[1:]
    }
    assert rows["full"] == ["2", "0", "0"]
    assert rows["cut"] == ["1", "0", "0"]


def test_assemblies_are_mapped_with_the_asm20_preset_unless_overridden(
    tmp_path: Path, monkeypatch
) -> None:
    from repgenr.snptypers import simple as mod
    from repgenr.snptypers.base import SnpParams

    ref = tmp_path / "reference.fasta"
    ref.write_text(">c1\nACGT\n", encoding="utf-8")
    work = tmp_path / "work"
    work.mkdir()
    log = logging.getLogger("t")

    def minimap2_args(params: SnpParams) -> list[str]:
        calls: list[list[str]] = []
        _fake_genome_chain(monkeypatch, _sam(("g", 0, "c1", 1, "4M")), "ACGT", calls)
        mod._call_one(tmp_path / "g.fasta", ref, work, 1, params, log, contigs=[("c1", 4)])
        return next(c for c in calls if c[0] == "minimap2")

    default = minimap2_args(SnpParams())
    assert default[default.index("-x") + 1] == "asm20"
    custom = minimap2_args(SnpParams(extra={"preset": "asm5"}))
    assert custom[custom.index("-x") + 1] == "asm5"
    assert "-x" not in minimap2_args(SnpParams(extra={"preset": "none"}))
    assert "preset" in mod.SimpleSnpTyper.capabilities.accepted_extras


def test_reference_contigs_are_read_in_file_order(tmp_path: Path) -> None:
    from repgenr.snptypers.simple import _reference_contigs

    ref = tmp_path / "reference.fasta"
    ref.write_text(">c2 desc\nAC\nGT\n>c1\nA\n", encoding="utf-8")
    assert _reference_contigs(ref) == [("c2", 4), ("c1", 1)]


def test_gzipped_reference_and_query_genomes_are_read(tmp_path: Path, monkeypatch) -> None:
    """A .fasta.gz reference is decompressed for faidx; minimap2 gets the query as it is."""
    import gzip

    from repgenr.snptypers import simple as mod
    from repgenr.snptypers.base import SnpParams

    ref = tmp_path / "refgenome.fasta.gz"
    with gzip.open(ref, "wt", encoding="utf-8") as fh:
        fh.write(">c1\nACGTACGT\n")
    query = tmp_path / "g1.fasta.gz"
    with gzip.open(query, "wt", encoding="utf-8") as fh:
        fh.write(">c1\nACGAACGT\n")

    calls: list[list[str]] = []
    _fake_genome_chain(monkeypatch, _sam(("c1", 0, "c1", 1, "8M")), "ACGAACGT", calls)
    monkeypatch.setattr(mod, "run_tool", lambda *a, **k: None)

    out = tmp_path / "out"
    result = mod.SimpleSnpTyper().call(
        [ref, query], ref, out, SnpParams(threads=1), logging.getLogger("t")
    )
    assert (out / "reference.fasta").read_text(encoding="utf-8") == ">c1\nACGTACGT\n"
    minimap2 = next(c for c in calls if c[0] == "minimap2")
    assert minimap2[-1] == str(query.resolve())
    # Records are named without the FASTA suffix and .gz, as elsewhere.
    assert _read_fasta(result.core_snp_fasta) == {"refgenome": "T", "g1": "A"}


def test_a_deletion_in_the_genome_is_n_in_its_consensus(tmp_path: Path, monkeypatch) -> None:
    from repgenr.snptypers import simple as mod
    from repgenr.snptypers.base import SnpParams

    ref = tmp_path / "reference.fasta"
    ref.write_text(">c1\nACGTACGTAC\n", encoding="utf-8")
    work = tmp_path / "work"
    work.mkdir()
    _fake_genome_chain(monkeypatch, _sam(("g", 0, "c1", 1, "3M4D3M")), "ACGTACGTAC", [])
    seq = mod._call_one(
        tmp_path / "g.fasta",
        ref,
        work,
        1,
        SnpParams(),
        logging.getLogger("t"),
        contigs=[("c1", 10)],
    )
    assert seq == "ACGNNNNTAC"


def test_pairs_without_a_shared_site_have_no_distance(tmp_path: Path) -> None:
    from repgenr.snptypers.simple import _write_core_snps

    consensuses = {"ref": "AAAA", "s1": "TTNN", "s2": "NNTT"}
    assert _write_core_snps(consensuses, tmp_path / "core.fasta", tmp_path / "dist.tsv") == 4
    rows = {
        line.split("\t")[0]: line.split("\t")[1:]
        for line in (tmp_path / "dist.tsv").read_text().splitlines()[1:]
    }
    assert rows["s1"] == ["2", "0", "NA"]
    assert rows["s2"] == ["2", "NA", "0"]


def _typer_with_consensuses(tmp_path: Path, monkeypatch, called: dict[str, str]):
    from repgenr.snptypers import simple as mod

    ref = tmp_path / "ref.fasta"
    ref.write_text(">c1\nAAAAAAAAAA\n", encoding="utf-8")
    genomes = [ref]
    for name in called:
        p = tmp_path / f"{name}.fasta"
        p.write_text(">c1\nAAAAAAAAAA\n", encoding="utf-8")
        genomes.append(p)
    monkeypatch.setattr(mod, "run_tool", lambda *a, **k: None)
    monkeypatch.setattr(mod, "_call_one", lambda genome, *a, **k: called[genome.stem])
    return ref, genomes


def test_a_genome_without_a_base_at_any_core_site_is_refused(tmp_path: Path, monkeypatch) -> None:
    """An unaligned genome would be all N; tree builders refuse such a sequence."""
    from repgenr.core.errors import WorkdirError
    from repgenr.snptypers import simple as mod
    from repgenr.snptypers.base import SnpParams

    ref, genomes = _typer_with_consensuses(
        tmp_path, monkeypatch, {"g1": "TAAAAAAAAA", "outg": "NNNNNNNNNN"}
    )
    with pytest.raises(WorkdirError) as exc:
        mod.SimpleSnpTyper().call(
            genomes, ref, tmp_path / "out", SnpParams(threads=1), logging.getLogger("t")
        )
    assert exc.value.exit_code == 3
    assert "outg" in str(exc.value)
    assert "preset=none" in str(exc.value)
    assert "alignment-free" in str(exc.value)


def test_coverage_is_logged_and_low_coverage_warned(tmp_path: Path, monkeypatch, caplog) -> None:
    from repgenr.snptypers import simple as mod
    from repgenr.snptypers.base import SnpParams

    ref, genomes = _typer_with_consensuses(
        tmp_path, monkeypatch, {"g1": "TAAAAAAAAA", "g2": "ATNNNNNNNN"}
    )
    with caplog.at_level(logging.INFO):
        mod.SimpleSnpTyper().call(
            genomes, ref, tmp_path / "out", SnpParams(threads=1), logging.getLogger("t")
        )
    text = caplog.text
    assert "g1: 100.0% of the reference covered" in text
    assert "g2: 20.0% of the reference covered" in text
    assert "minimum 20.0%" in text
    warned = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warned) == 1 and "g2 covers only 20.0%" in warned[0].getMessage()


def test_an_unknown_preset_is_refused_before_any_work(tmp_path: Path, monkeypatch) -> None:
    from repgenr.core.errors import UserInputError
    from repgenr.snptypers import simple as mod
    from repgenr.snptypers.base import SnpParams

    ref, genomes = _typer_with_consensuses(tmp_path, monkeypatch, {"g1": "TAAAAAAAAA"})
    with pytest.raises(UserInputError) as exc:
        mod.SimpleSnpTyper().call(
            genomes,
            ref,
            tmp_path / "out",
            SnpParams(threads=1, extra={"preset": "asm30"}),
            logging.getLogger("t"),
        )
    assert exc.value.exit_code == 2
    assert "asm20" in str(exc.value)
    assert not (tmp_path / "out").exists()
