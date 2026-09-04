"""Tests for the benchmark cell matrix and runner command construction."""

from __future__ import annotations

from benchmarks.cells import Cell, all_cells, chunked_cells, derep_cells, tiers
from benchmarks.run_bench import _RSS_RE, _command


def test_cell_ids_unique():
    ids = [c.id for c in all_cells()]
    assert len(ids) == len(set(ids))


def test_cell_n_from_set_name_and_subset():
    assert Cell(id="x", kind="derep_step", tool="skder", set_name="clonal_5000_random").n == 5000
    assert (
        Cell(
            id="y",
            kind="tree_step",
            tool="fasttree",
            set_name="balanced_1000_clustered",
            subset=100,
        ).n
        == 100
    )


def test_no_sourmash_treebuilder_cell_at_5000():
    # The sourmash tree builder is a pure-Python O(n^3) NJ; the matrix caps it
    # at n=1000 and extrapolates instead of scheduling n=5000.
    for cell in all_cells():
        if cell.kind == "tree_step" and cell.tool == "sourmash":
            assert cell.n <= 1000


def test_chunked_cells_are_stage_runs_with_process_size():
    cells = chunked_cells()
    assert cells, "chunked matrix must not be empty"
    for cell in cells:
        assert cell.kind == "derep_stage"
        assert "--process-size" in cell.extra_args
    orders = {c.set_name.rsplit("_", 1)[1] for c in cells}
    assert orders == {"clustered", "random"}  # the B2 order-sensitivity pair


def test_tiers_partition_all_cells():
    cells = all_cells()
    parts = tiers(cells)
    combined = parts["smoke"] + parts["mid"] + parts["heavy"]
    assert sorted(c.id for c in combined) == sorted(c.id for c in cells)
    assert all(c.n <= 100 for c in parts["smoke"])
    assert all(c.n <= 1000 for c in parts["mid"])


def test_derep_matrix_covers_tools_sizes_scenarios():
    cells = derep_cells()
    tools = {c.tool for c in cells if c.kind == "derep_step"}
    assert tools == {"skder", "galah", "sourmash"}
    assert any(c.kind == "derep_dense" for c in cells)
    sizes = {c.n for c in cells}
    assert sizes == {100, 1000, 5000}


def test_command_derep_step_argv(tmp_path):
    set_dir = tmp_path / "set"
    set_dir.mkdir()
    for name in ("b.fasta", "a.fasta", "._a.fasta"):  # AppleDouble must be skipped
        (set_dir / name).write_text(">x\nACGT\n", encoding="utf-8")
    work = tmp_path / "work"
    work.mkdir()
    cell = Cell(id="t", kind="derep_step", tool="skder", set_name="clonal_100_clustered")
    argv, out = _command(cell, set_dir, work)
    assert "dereplicate-chunk" in argv
    assert argv[0].endswith("repgenr")  # console script, not python -m
    assert "--genomes-fofn" in argv
    assert argv[argv.index("--tool") + 1] == "skder"
    fofn = work / "genomes.fofn"
    lines = fofn.read_text(encoding="utf-8").splitlines()
    assert [line.rsplit("/", 1)[1] for line in lines] == ["a.fasta", "b.fasta"]  # sorted
    assert out == work / "derep_out"


def test_command_tree_step_subset_copies_first_n(tmp_path):
    set_dir = tmp_path / "set"
    set_dir.mkdir()
    for i in range(4):
        (set_dir / f"g{i}.fasta").write_text(">x\nACGT\n", encoding="utf-8")
    work = tmp_path / "work"
    work.mkdir()
    cell = Cell(
        id="t", kind="tree_step", tool="fasttree", set_name="balanced_1000_clustered", subset=2
    )
    argv, _out = _command(cell, set_dir, work)
    assert "phylo-build" in argv
    assert "--no-outgroup" in argv
    copied = sorted(p.name for p in (work / "genomes_subset").glob("*.fasta"))
    assert copied == ["g0.fasta", "g1.fasta"]


def test_rss_regex_parses_time_l_output():
    stderr = "        7.51 real ...\n            123456789  maximum resident set size\n"
    match = _RSS_RE.search(stderr)
    assert match is not None
    assert int(match.group(1)) == 123456789
