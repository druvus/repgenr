"""sourmash sketch caching: reuse signatures across --target-reps iterations.

The cache directory may be shared between sequential --target-reps iterations
(same genome set) and, under chunked dereplication, between parallel chunk
workers (disjoint genome sets). Reuse must therefore be keyed to the requested
genome set, never to a bare file count: a chunk must only ever compare its own
genomes' signatures, and a cached zip from a different genome set must not be
mistaken for this one's.
"""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.dereplicators import sourmash
from repgenr.dereplicators.base import DerepParams
from repgenr.dereplicators.sourmash import SourmashDereplicator

_LOG = logging.getLogger("test")


def _genomes(tmp_path: Path, names: tuple[str, ...] = ("g1.fasta", "g2.fasta")) -> list[Path]:
    out = []
    for name in names:
        p = tmp_path / name
        p.write_text(">x\nACGT\n")
        out.append(p)
    return out


def _fake_run_tool(calls: list[tuple[str, list[str]]]):
    """Record (subcommand, inputs) and fabricate each tool's output files.

    sketch: writes one ``<stem>.sig`` per path listed in the ``--from-file``
    fofn. compare: records the signature fofn and writes an all-similar matrix.
    manysketch: records the genome column of its CSV and creates the ``-o`` zip.
    pairwise: writes a header-only edge list (all genomes become singletons).
    """

    def run_tool(caps, cmd, **kwargs):
        parts = [str(c) for c in cmd]
        if "manysketch" in parts:
            csv_path = Path(parts[parts.index("manysketch") + 1])
            rows = csv_path.read_text().splitlines()[1:]
            calls.append(("manysketch", [r.split(",")[1] for r in rows]))
            Path(parts[parts.index("-o") + 1]).write_text("zip")
        elif "pairwise" in parts:
            calls.append(("pairwise", [parts[parts.index("pairwise") + 1]]))
            Path(parts[parts.index("-o") + 1]).write_text("query_name,match_name,jaccard\n")
        elif "sketch" in parts:
            fofn = Path(parts[parts.index("--from-file") + 1])
            paths = fofn.read_text().splitlines()
            calls.append(("sketch", paths))
            outdir = Path(parts[parts.index("--outdir") + 1])
            outdir.mkdir(parents=True, exist_ok=True)
            for p in paths:
                (outdir / f"{Path(p).stem}.sig").write_text("sig")
        elif "compare" in parts:
            fofn = Path(parts[parts.index("--from-file") + 1])
            sigs = fofn.read_text().splitlines()
            calls.append(("compare", sigs))
            labels = [Path(s).stem for s in sigs]
            n = len(labels)
            body = "\n".join(",".join("1.0" for _ in range(n)) for _ in range(n))
            Path(parts[parts.index("--csv") + 1]).write_text(",".join(labels) + "\n" + body + "\n")
        return 0

    return run_tool


def _names(calls: list[tuple[str, list[str]]], kind: str) -> list[list[str]]:
    return [[Path(p).name for p in inputs] for c, inputs in calls if c == kind]


# --- dense path ---------------------------------------------------------------


def test_dense_sketches_when_cache_empty(tmp_path: Path, monkeypatch) -> None:
    genomes = _genomes(tmp_path)
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))
    cache = tmp_path / "sketches"

    clusters, _ = SourmashDereplicator()._dense_dereplicate(
        genomes, tmp_path / "out0", 31, 1000, 0.9, _LOG, sketch_cache=cache
    )
    assert _names(calls, "sketch") == [["g1.fasta", "g2.fasta"]]
    assert len(clusters) == 1  # identical genomes collapse to one representative


def test_dense_reuses_cached_signatures(tmp_path: Path, monkeypatch) -> None:
    genomes = _genomes(tmp_path)
    cache = tmp_path / "sketches"
    cache.mkdir()
    for g in genomes:  # pre-populate as a prior iteration would have
        (cache / f"{g.stem}.sig").write_text("sig")

    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    SourmashDereplicator()._dense_dereplicate(
        genomes, tmp_path / "out1", 31, 1000, 0.9, _LOG, sketch_cache=cache
    )
    assert _names(calls, "sketch") == []  # reused the cache, only compared
    assert len(_names(calls, "compare")) == 1


def test_dense_subset_of_cache_compares_only_subset(tmp_path: Path, monkeypatch) -> None:
    # Cache holds signatures for a *superset* (another chunk's genomes included);
    # this chunk must compare exactly its own two, not everything in the dir.
    all_genomes = _genomes(tmp_path, ("g1.fasta", "g2.fasta", "g3.fasta"))
    cache = tmp_path / "sketches"
    cache.mkdir()
    for g in all_genomes:
        (cache / f"{g.stem}.sig").write_text("sig")

    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    clusters, status = SourmashDereplicator()._dense_dereplicate(
        all_genomes[:2], tmp_path / "out", 31, 1000, 0.9, _LOG, sketch_cache=cache
    )
    assert _names(calls, "compare") == [["g1.sig", "g2.sig"]]
    assert "g3.fasta" not in status


def test_dense_foreign_cache_entries_do_not_satisfy_reuse(tmp_path: Path, monkeypatch) -> None:
    # Two cached signatures exist, but for *different* genomes: a bare count
    # check would wrongly skip sketching.
    genomes = _genomes(tmp_path)
    cache = tmp_path / "sketches"
    cache.mkdir()
    (cache / "h1.sig").write_text("sig")
    (cache / "h2.sig").write_text("sig")

    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    SourmashDereplicator()._dense_dereplicate(
        genomes, tmp_path / "out", 31, 1000, 0.9, _LOG, sketch_cache=cache
    )
    assert _names(calls, "sketch") == [["g1.fasta", "g2.fasta"]]
    assert _names(calls, "compare") == [["g1.sig", "g2.sig"]]


def test_dense_sketches_only_missing_genomes(tmp_path: Path, monkeypatch) -> None:
    genomes = _genomes(tmp_path)
    cache = tmp_path / "sketches"
    cache.mkdir()
    (cache / "g1.sig").write_text("sig")  # g2 not cached yet

    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    SourmashDereplicator()._dense_dereplicate(
        genomes, tmp_path / "out", 31, 1000, 0.9, _LOG, sketch_cache=cache
    )
    assert _names(calls, "sketch") == [["g2.fasta"]]
    assert _names(calls, "compare") == [["g1.sig", "g2.sig"]]


# --- sparse path --------------------------------------------------------------


def _sparse(genomes: list[Path], out_dir: Path, cache: Path):
    out_dir.mkdir(parents=True, exist_ok=True)  # done by dereplicate() in production
    params = DerepParams(primary_ani=0.99, secondary_ani=0.99, threads=1)
    return SourmashDereplicator()._sparse_dereplicate(
        genomes, out_dir, 31, 1000, 0.99, params, _LOG, sketch_cache=cache
    )


def test_sparse_cache_reused_for_same_genome_set(tmp_path: Path, monkeypatch) -> None:
    genomes = _genomes(tmp_path)
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    _sparse(genomes, tmp_path / "iter0", cache)
    _sparse(genomes, tmp_path / "iter1", cache)
    # same genome set: sketched once, second iteration reused the zip
    assert len(_names(calls, "manysketch")) == 1
    assert len(_names(calls, "pairwise")) == 2


def test_sparse_cache_keyed_by_genome_set(tmp_path: Path, monkeypatch) -> None:
    # Two disjoint chunks sharing one cache dir must not reuse each other's zip.
    chunk_a = _genomes(tmp_path, ("g1.fasta", "g2.fasta"))
    chunk_b = _genomes(tmp_path, ("g3.fasta", "g4.fasta"))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    clusters_a, _ = _sparse(chunk_a, tmp_path / "chunk_a", cache)
    clusters_b, _ = _sparse(chunk_b, tmp_path / "chunk_b", cache)

    assert _names(calls, "manysketch") == [
        ["g1.fasta", "g2.fasta"],
        ["g3.fasta", "g4.fasta"],
    ]
    assert set(clusters_a) == {"g1.fasta", "g2.fasta"}
    assert set(clusters_b) == {"g3.fasta", "g4.fasta"}
