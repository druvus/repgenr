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
import os
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
            rows = csv_path.read_text(encoding="utf-8").splitlines()[1:]
            calls.append(("manysketch", [r.split(",")[1] for r in rows]))
            Path(parts[parts.index("-o") + 1]).write_text("zip", encoding="utf-8")
        elif "sig" in parts and "cat" in parts:
            start = parts.index("cat") + 1
            zips = []
            for tok in parts[start:]:
                if tok.startswith("-"):
                    break
                zips.append(tok)
            pick = parts[parts.index("--picklist") + 1].split(":")[0]
            names = Path(pick).read_text(encoding="utf-8").splitlines()[1:]
            calls.append(("sigcat", zips))
            calls.append(("picklist", names))
            Path(parts[parts.index("-o") + 1]).write_text("zip", encoding="utf-8")
        elif "pairwise" in parts:
            calls.append(("pairwise", [parts[parts.index("pairwise") + 1]]))
            Path(parts[parts.index("-o") + 1]).write_text(
                "query_name,match_name,jaccard\n", encoding="utf-8"
            )
        elif "sketch" in parts:
            fofn = Path(parts[parts.index("--from-file") + 1])
            paths = fofn.read_text(encoding="utf-8").splitlines()
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
            Path(parts[parts.index("--csv") + 1]).write_text(
                ",".join(labels) + "\n" + body + "\n", encoding="utf-8"
            )
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


# --- sparse path: merge-level collections from cached chunk zips ---------------


def test_union_of_chunk_genomes_is_collected_not_resketched(tmp_path: Path, monkeypatch) -> None:
    # Chunked --target-reps: the merge level dereplicates a union of chunk
    # representatives that changes with the threshold. Its signatures are
    # already in the chunk zips, so the union zip is assembled with sig cat.
    chunk_a = _genomes(tmp_path, ("g1.fasta", "g2.fasta"))
    chunk_b = _genomes(tmp_path, ("g3.fasta", "g4.fasta"))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    _sparse(chunk_a, tmp_path / "chunk_a", cache)
    _sparse(chunk_b, tmp_path / "chunk_b", cache)
    calls.clear()
    clusters, _ = _sparse([chunk_a[0], chunk_b[1]], tmp_path / "merge", cache)

    assert _names(calls, "manysketch") == []
    assert len(_names(calls, "sigcat")) == 1
    assert len(_names(calls, "sigcat")[0]) == 2  # one zip per chunk
    assert _names(calls, "picklist") == [["g1", "g4"]]
    assert set(clusters) == {"g1.fasta", "g4.fasta"}
    # pairwise runs on the assembled zip, keyed by the union's own digest
    digest = sourmash._genome_set_digest([chunk_a[0], chunk_b[1]], 31, 1000)
    assert _names(calls, "pairwise") == [[f"signatures-{digest}.zip"]]


def test_union_collection_is_reused_for_the_same_union(tmp_path: Path, monkeypatch) -> None:
    chunk = _genomes(tmp_path, ("g1.fasta", "g2.fasta", "g3.fasta"))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    _sparse(chunk, tmp_path / "chunk", cache)
    _sparse(chunk[:2], tmp_path / "merge0", cache)
    _sparse(chunk[:2], tmp_path / "merge1", cache)
    assert len(_names(calls, "manysketch")) == 1
    assert len(_names(calls, "sigcat")) == 1


def test_genome_outside_the_cached_zips_is_sketched(tmp_path: Path, monkeypatch) -> None:
    chunk = _genomes(tmp_path, ("g1.fasta", "g2.fasta"))
    extra = _genomes(tmp_path, ("g9.fasta",))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    _sparse(chunk, tmp_path / "chunk", cache)
    _sparse([chunk[0], *extra], tmp_path / "merge", cache)
    assert _names(calls, "sigcat") == []
    assert _names(calls, "manysketch") == [["g1.fasta", "g2.fasta"], ["g1.fasta", "g9.fasta"]]


def test_cached_zips_with_other_sketch_parameters_are_not_collected(
    tmp_path: Path, monkeypatch
) -> None:
    chunk = _genomes(tmp_path, ("g1.fasta", "g2.fasta", "g3.fasta"))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    _sparse(chunk, tmp_path / "chunk", cache)  # k=31
    params = DerepParams(secondary_ani=0.99, threads=1)
    (tmp_path / "k21").mkdir()
    SourmashDereplicator()._sparse_dereplicate(
        chunk[:2], tmp_path / "k21", 21, 1000, 0.99, params, _LOG, sketch_cache=cache
    )
    assert _names(calls, "sigcat") == []
    assert len(_names(calls, "manysketch")) == 2


def test_collected_zips_are_not_offered_again(tmp_path: Path, monkeypatch) -> None:
    # An assembled union zip overlaps the chunk zips; offering it as a source
    # would let sig cat emit a genome twice. Only sketched zips are sources.
    chunk_a = _genomes(tmp_path, ("g1.fasta", "g2.fasta", "g3.fasta"))
    chunk_b = _genomes(tmp_path, ("g4.fasta", "g5.fasta", "g6.fasta"))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    _sparse(chunk_a, tmp_path / "a", cache)
    _sparse(chunk_b, tmp_path / "b", cache)
    _sparse([*chunk_a, *chunk_b[:2]], tmp_path / "u0", cache)
    calls.clear()
    _sparse([chunk_a[0], *chunk_b], tmp_path / "u1", cache)
    zips = _names(calls, "sigcat")[0]
    digests = {
        sourmash._genome_set_digest(chunk_a, 31, 1000),
        sourmash._genome_set_digest(chunk_b, 31, 1000),
    }
    assert sorted(zips) == sorted(f"signatures-{d}.zip" for d in digests)


def test_failed_collection_falls_back_to_sketching(tmp_path: Path, monkeypatch) -> None:
    chunk = _genomes(tmp_path, ("g1.fasta", "g2.fasta", "g3.fasta"))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    fake = _fake_run_tool(calls)

    def run_tool(caps, cmd, **kwargs):
        parts = [str(c) for c in cmd]
        if "sig" in parts and "cat" in parts:
            raise sourmash.ToolExecutionError(parts, 2, "unknown option", tool="sourmash")
        return fake(caps, cmd, **kwargs)

    monkeypatch.setattr(sourmash, "run_tool", run_tool)
    _sparse(chunk, tmp_path / "chunk", cache)
    clusters, _ = _sparse(chunk[:2], tmp_path / "merge", cache)
    assert _names(calls, "manysketch") == [
        ["g1.fasta", "g2.fasta", "g3.fasta"],
        ["g1.fasta", "g2.fasta"],
    ]
    assert set(clusters) == {"g1.fasta", "g2.fasta"}


def test_chunked_target_reps_sketches_each_chunk_once(tmp_path: Path, monkeypatch) -> None:
    from repgenr.stages.dereplicate import DereplicateParams, _search_target_reps

    genomes = _genomes(tmp_path, tuple(f"g{i}.fasta" for i in range(6)))
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))
    monkeypatch.setattr(sourmash, "_branchwater_available", lambda caps, logger: True)
    params = DereplicateParams(tool="sourmash", process_size=3, num_processes=1, threads=1)
    derep_params = DerepParams(threads=1, extra={"ksize": 31, "scaled": 1000})

    _search_target_reps(
        SourmashDereplicator(), genomes, tmp_path / "scratch", derep_params, params, 2, _LOG
    )
    # two chunks sketched once; every search step collects its merge-level union
    assert len(_names(calls, "manysketch")) == 2
    assert len(_names(calls, "pairwise")) > 3


# --- cache keys follow the genome file, not only its name ------------------------


def _bump(path: Path, text: str) -> None:
    """Replace a genome under the same name with different content and a later mtime."""
    st = path.stat()
    path.write_text(text, encoding="utf-8")
    os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))


def test_replaced_genome_is_not_collected_from_the_old_zip(tmp_path: Path, monkeypatch) -> None:
    chunk = _genomes(tmp_path, ("g1.fasta", "g2.fasta", "g3.fasta"))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    _sparse(chunk, tmp_path / "chunk", cache)
    _bump(chunk[0], ">x\nACGTACGTACGT\n")
    _sparse(chunk[:2], tmp_path / "merge", cache)
    assert _names(calls, "sigcat") == []
    assert _names(calls, "manysketch")[-1] == ["g1.fasta", "g2.fasta"]


def test_replaced_genome_changes_the_set_digest(tmp_path: Path) -> None:
    genomes = _genomes(tmp_path)
    before = sourmash._genome_set_digest(genomes, 31, 1000)
    _bump(genomes[1], ">x\nACGTAC\n")
    assert sourmash._genome_set_digest(genomes, 31, 1000) != before


def test_dense_signature_older_than_its_genome_is_resketched(tmp_path: Path, monkeypatch) -> None:
    genomes = _genomes(tmp_path)
    cache = tmp_path / "sketches"
    cache.mkdir()
    for g in genomes:
        (cache / f"{g.stem}.sig").write_text("sig")
    # g1.sig was made from an earlier g1.fasta: older than the current file
    old = genomes[0].stat().st_mtime_ns - 10 * 10**9
    os.utime(cache / "g1.sig", ns=(old, old))
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    SourmashDereplicator()._dense_dereplicate(
        genomes, tmp_path / "out", 31, 1000, 0.9, _LOG, sketch_cache=cache
    )
    assert _names(calls, "sketch") == [["g1.fasta"]]


def test_names_with_commas_stay_one_csv_field(tmp_path: Path, monkeypatch) -> None:
    import csv

    chunk = _genomes(tmp_path, ("a,b.fasta", "c.fasta", "d.fasta"))
    cache = tmp_path / "sketches"
    calls: list[tuple[str, list[str]]] = []
    monkeypatch.setattr(sourmash, "run_tool", _fake_run_tool(calls))

    _sparse(chunk, tmp_path / "chunk", cache)
    _sparse(chunk[:2], tmp_path / "merge", cache)
    with open(tmp_path / "chunk" / "manysketch.csv", encoding="utf-8", newline="") as fo:
        rows = list(csv.reader(fo))
    assert [r[0] for r in rows[1:]] == ["a,b", "c", "d"]
    assert all(len(r) == 3 for r in rows)
    with open(tmp_path / "merge" / "picklist.csv", encoding="utf-8", newline="") as fo:
        assert list(csv.reader(fo)) == [["name"], ["a,b"], ["c"]]
