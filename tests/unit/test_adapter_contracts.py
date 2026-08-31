"""Contract tests over every registered dereplicator and treebuilder adapter.

Each adapter runs against a faked ``run_tool`` that records argv and writes
canned tool outputs, so the public contract is pinned without the real
binaries: every input genome ends up in exactly one cluster, normalized
parameters land in argv, and tree builders return an existing Newick file.
"""

from __future__ import annotations

import logging
import shutil
import sys
from pathlib import Path

import pytest

from repgenr.dereplicators.base import DerepParams, DerepResult
from repgenr.dereplicators.base import registry as derep_registry
from repgenr.treebuilders.base import InputKind, TreeParams
from repgenr.treebuilders.base import registry as tree_registry

_LOG = logging.getLogger("test")

_GENOMES = ["g1.fasta", "g2.fasta", "g3.fasta", "g4.fasta"]
# canned clustering: g1 represents {g2}, g3 represents {g4}
_REP_OF = {"g1.fasta": "g1.fasta", "g2.fasta": "g1.fasta",
           "g3.fasta": "g3.fasta", "g4.fasta": "g3.fasta"}
_NEWICK = "(g1:0.1,g2:0.1,(g3:0.1,g4:0.1):0.1);\n"


@pytest.fixture()
def genomes(tmp_path) -> list[Path]:
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    out = []
    for name in _GENOMES:
        p = gdir / name
        p.write_text(f">{Path(name).stem}\nACGTACGTACGT\n")
        out.append(p)
    return out


def _flag_value(cmd: list[str], flag: str) -> str:
    return cmd[cmd.index(flag) + 1]


def _read_fofn(path: str) -> list[Path]:
    return [Path(line) for line in Path(path).read_text(encoding="utf-8").splitlines() if line]


# --- canned tool output writers -----------------------------------------------


def _fake_skder(cmd: list[str]) -> None:
    result_dir = Path(_flag_value(cmd, "-o"))
    i = cmd.index("-g") + 1
    genome_args: list[Path] = []
    while i < len(cmd) and not cmd[i].startswith("-"):
        genome_args.append(Path(cmd[i]))
        i += 1
    rep_dir = result_dir / "Dereplicated_Representative_Genomes"
    rep_dir.mkdir(parents=True)
    by_name = {p.name: p for p in genome_args}
    for rep in ("g1.fasta", "g3.fasta"):
        shutil.copy2(by_name[rep], rep_dir / rep)
    edges = result_dir / "Skani_Triangle_Edge_Output.txt"
    lines = ["#query\treference\tANI\tAF_query\tAF_reference"]
    for member, rep in _REP_OF.items():
        if member != rep:
            lines.append(f"{by_name[member]}\t{by_name[rep]}\t99.5\t95.0\t95.0")
    edges.write_text("\n".join(lines) + "\n")


def _fake_drep(cmd: list[str]) -> None:
    drep_wd = Path(cmd[2])
    staged = _read_fofn(_flag_value(cmd, "-g"))
    by_name = {p.name: p for p in staged}
    rep_dir = drep_wd / "dereplicated_genomes"
    rep_dir.mkdir(parents=True)
    for rep in ("g1.fasta", "g3.fasta"):
        shutil.copy2(by_name[rep], rep_dir / rep)
    tables = drep_wd / "data_tables"
    tables.mkdir()
    rows = ["genome,secondary_cluster"]
    cluster_id = {"g1.fasta": "1_1", "g3.fasta": "2_1"}
    for member, rep in _REP_OF.items():
        rows.append(f"{member},{cluster_id[rep]}")
    (tables / "Cdb.csv").write_text("\n".join(rows) + "\n")


def _fake_galah(cmd: list[str]) -> None:
    clusters_file = Path(_flag_value(cmd, "--output-cluster-definition"))
    listed = _read_fofn(_flag_value(cmd, "--genome-fasta-list"))
    by_name = {p.name: p for p in listed}
    rows = [f"{by_name[rep]}\t{by_name[member]}" for member, rep in _REP_OF.items()]
    clusters_file.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _fake_sourmash_sketch(cmd: list[str]) -> None:
    outdir = Path(_flag_value(cmd, "--outdir"))
    for p in _read_fofn(_flag_value(cmd, "--from-file")):
        (outdir / f"{p.name}.sig").write_text("{}")


def _fake_sourmash_compare(cmd: list[str]) -> None:
    csv_path = Path(_flag_value(cmd, "--csv"))
    sigs = _read_fofn(_flag_value(cmd, "--from-file"))
    labels = [p.name.removesuffix(".sig") for p in sigs]

    def sim(a: str, b: str) -> float:
        if a == b:
            return 1.0
        return 0.995 if _REP_OF[a] == _REP_OF[b] else 0.2

    rows = [",".join(labels)]
    for a in labels:
        rows.append(",".join(f"{sim(a, b):.4f}" for b in labels))
    csv_path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _make_fake_run_tool(recorded: list[list[str]]):
    def fake_run_tool(caps, command, *, logger, stdout_path=None, **kwargs):
        cmd = [str(part) for part in command]
        recorded.append(cmd)
        if cmd[:3] == ["sourmash", "scripts", "pairwise"]:
            return 1  # no branchwater plugin -> dense path
        if cmd[0] == "skder":
            _fake_skder(cmd)
        elif cmd[0] == "dRep" and cmd[1] == "dereplicate":
            _fake_drep(cmd)
        elif cmd[0] == "galah":
            _fake_galah(cmd)
        elif cmd[:2] == ["sourmash", "sketch"]:
            _fake_sourmash_sketch(cmd)
        elif cmd[:2] == ["sourmash", "compare"]:
            _fake_sourmash_compare(cmd)
        elif cmd[0] in ("mashtree", "FastTree", "fasttree"):
            Path(stdout_path).write_text(_NEWICK, encoding="utf-8")
        elif cmd[0] == "iqtree":
            msa = Path(_flag_value(cmd, "-s"))
            Path(str(msa) + ".treefile").write_text(_NEWICK, encoding="utf-8")
        elif cmd[0] == "raxml-ng":
            prefix = _flag_value(cmd, "--prefix")
            Path(prefix + ".raxml.bestTree").write_text(_NEWICK, encoding="utf-8")
        else:
            # A third-party adapter's unknown argv must not crash the whole
            # suite; its own (skipped) parametrization covers it.
            return 0
        return 0

    return fake_run_tool


@pytest.fixture()
def recorded(monkeypatch) -> list[list[str]]:
    """Patch run_tool in every adapter module; return the recorded argvs."""
    calls: list[list[str]] = []
    fake = _make_fake_run_tool(calls)
    for reg in (derep_registry, tree_registry):
        for name in reg.names():
            module = sys.modules[reg.get(name).__module__]
            monkeypatch.setattr(module, "run_tool", fake, raising=False)
    import repgenr.dereplicators.sourmash as sm

    sm._BRANCHWATER_CACHE.clear()
    return calls


# --- dereplicator contract ----------------------------------------------------


def _assert_partition(result: DerepResult) -> None:
    covered: list[str] = []
    for rep, members in result.clusters.items():
        covered.append(rep)
        covered.extend(members)
    assert sorted(covered) == sorted(_GENOMES), "every genome in exactly one cluster"
    assert {p.name for p in result.representatives} == set(result.clusters)
    assert set(result.genome_status) == set(_GENOMES)


# tool -> tokens that must land in some recorded argv for these params
_DEREP_PARAM_TOKENS = {
    "skder": ["-i", "99", "-c", "7"],
    "drep": ["-sa", "0.99", "--processors", "7"],
    "galah": ["--ani", "99", "--threads", "7"],
    "sourmash": ["k=31,scaled=1000"],
}


@pytest.mark.parametrize("tool", sorted(_DEREP_PARAM_TOKENS))
def test_dereplicator_contract(tool, genomes, recorded, tmp_path) -> None:
    if tool not in derep_registry.names():
        pytest.skip(f"{tool} not registered")
    adapter = derep_registry.create(tool)
    params = DerepParams(primary_ani=0.90, secondary_ani=0.99,
                         aligned_fraction=0.50, threads=7)
    result = adapter.dereplicate(genomes, tmp_path / "out", params, _LOG)

    _assert_partition(result)
    assert len(result.clusters) == 2
    for rep in result.representatives:
        assert rep.exists(), f"representative {rep} must exist on disk"

    flat = [tok for cmd in recorded for tok in cmd]
    for token in _DEREP_PARAM_TOKENS[tool]:
        assert token in flat, f"{token!r} missing from recorded argv for {tool}"


def test_galah_empty_clusters_yields_empty_result(genomes, recorded, tmp_path, monkeypatch) -> None:
    """Pinned as-is: an empty galah clusters.tsv produces an empty DerepResult."""
    import repgenr.dereplicators.galah as galah_mod

    def empty_galah(caps, command, *, logger, **kwargs):
        cmd = [str(c) for c in command]
        Path(_flag_value(cmd, "--output-cluster-definition")).write_text("", encoding="utf-8")
        return 0

    monkeypatch.setattr(galah_mod, "run_tool", empty_galah)
    result = galah_mod.GalahDereplicator().dereplicate(
        genomes, tmp_path / "out", DerepParams(), _LOG
    )
    assert result.representatives == []
    assert result.clusters == {}
    assert result.genome_status == {}


# --- treebuilder contract -----------------------------------------------------


_TREE_PARAM_TOKENS = {
    "mashtree": ["--numcpus", "5"],
    "fasttree": ["-nt", "-gtr"],
    "iqtree": ["--threads-max", "5"],
    "raxmlng": ["--threads", "auto{5}"],
    "sourmash": ["k=31,scaled=1000"],
}


@pytest.fixture()
def msa(tmp_path) -> Path:
    path = tmp_path / "input_msa.fasta"
    path.write_text("".join(f">{Path(n).stem}\nACGTACGT\n" for n in _GENOMES))
    return path


@pytest.mark.parametrize("tool", sorted(_TREE_PARAM_TOKENS))
def test_treebuilder_contract(tool, genomes, msa, recorded, tmp_path) -> None:
    if tool not in tree_registry.names():
        pytest.skip(f"{tool} not registered")
    builder = tree_registry.create(tool)
    source = genomes if builder.input_kind is InputKind.GENOMES else msa
    tree = builder.build(source, tmp_path / "tree_out", TreeParams(threads=5), _LOG)

    assert tree.exists()
    content = tree.read_text()
    assert content.rstrip().endswith(";"), "output must be a Newick tree"
    for name in _GENOMES:
        assert Path(name).stem in content

    flat = [tok for cmd in recorded for tok in cmd]
    for token in _TREE_PARAM_TOKENS[tool]:
        assert token in flat, f"{token!r} missing from recorded argv for {tool}"


def test_every_registered_dereplicator_has_contract_coverage() -> None:
    """A new in-tree adapter must add itself to the contract parametrization."""
    assert set(derep_registry.names()) >= set(_DEREP_PARAM_TOKENS)
    builtin = {"drep", "galah", "skder", "sourmash"}
    assert builtin & set(derep_registry.names()) <= set(_DEREP_PARAM_TOKENS)


def test_every_registered_treebuilder_has_contract_coverage() -> None:
    builtin = {"fasttree", "iqtree", "mashtree", "raxmlng", "sourmash"}
    assert builtin & set(tree_registry.names()) <= set(_TREE_PARAM_TOKENS)
