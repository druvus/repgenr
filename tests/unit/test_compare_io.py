"""compare_io: the pairwise CSV and dendrogram writers behind glance's sourmash backend."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless: no display in CI

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from repgenr.dereplicators.compare_io import (  # noqa: E402
    average_linkage,
    leaf_order,
    write_dendrogram,
    write_pairwise_csv,
)


def _two_groups() -> tuple[list[str], np.ndarray]:
    # a,b close; c,d close; the groups far apart.
    names = ["a.fasta", "b.fasta", "c.fasta", "d.fasta"]
    sim = np.array(
        [
            [1.00, 0.99, 0.90, 0.91],
            [0.99, 1.00, 0.92, 0.90],
            [0.90, 0.92, 1.00, 0.98],
            [0.91, 0.90, 0.98, 1.00],
        ]
    )
    return names, sim


def test_pairwise_csv_has_one_row_per_unordered_pair(tmp_path: Path) -> None:
    names, sim = _two_groups()
    out = write_pairwise_csv(names, sim, tmp_path / "sub" / "pairs.csv")
    with open(out, encoding="utf-8", newline="") as fo:
        rows = list(csv.DictReader(fo))
    assert list(rows[0]) == ["genome1", "genome2", "similarity"]
    assert len(rows) == len(names) * (len(names) - 1) // 2
    pairs = {(r["genome1"], r["genome2"]): float(r["similarity"]) for r in rows}
    assert pairs[("a.fasta", "b.fasta")] == 0.99
    assert pairs[("c.fasta", "d.fasta")] == 0.98
    assert all(g1 != g2 for g1, g2 in pairs)  # no self-comparisons


def test_pairwise_csv_is_read_by_glance_pair_counter(tmp_path: Path) -> None:
    from repgenr.stages.glance import _pair_similarities

    names, sim = _two_groups()
    out = write_pairwise_csv(names, sim, tmp_path / "pairs.csv")
    assert sorted(_pair_similarities(out, 0.0, 1.0)) == sorted(
        sim[np.triu_indices(4, k=1)].tolist()
    )


def test_writers_reject_mismatched_shapes(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="square"):
        write_pairwise_csv(["a"], np.zeros((1, 2)), tmp_path / "x.csv")
    with pytest.raises(ValueError, match="names"):
        write_pairwise_csv(["a", "b", "c"], np.eye(2), tmp_path / "x.csv")
    with pytest.raises(ValueError, match="at least two"):
        average_linkage(np.zeros((1, 1)))


def test_average_linkage_on_two_groups() -> None:
    names, sim = _two_groups()
    z = average_linkage(1.0 - sim)
    assert z.shape == (3, 4)
    # First the closest pair (a,b) at 0.01, then (c,d) at 0.02, then the
    # groups at the mean of the four between-group distances.
    assert z[0].tolist() == pytest.approx([0, 1, 0.01, 2])
    assert z[1].tolist() == pytest.approx([2, 3, 0.02, 2])
    assert z[2].tolist() == pytest.approx([4, 5, 1 - np.mean([0.90, 0.91, 0.92, 0.90]), 4])
    assert leaf_order(z) == [0, 1, 2, 3]


@pytest.mark.parametrize("n", [2, 3, 7, 40, 150])
def test_average_linkage_matches_scipy(n: int) -> None:
    hierarchy = pytest.importorskip("scipy.cluster.hierarchy")
    distance = pytest.importorskip("scipy.spatial.distance")
    rng = np.random.default_rng(n)
    dist = distance.squareform(distance.pdist(rng.random((n, 4))))
    ours = average_linkage(dist)
    ref = hierarchy.linkage(distance.squareform(dist, checks=False), method="average")
    assert np.allclose(ours[:, 2], ref[:, 2])
    assert np.allclose(ours[:, 3], ref[:, 3])
    # Same tree: every merge joins the same two leaf sets.
    assert _merged_sets(ours) == _merged_sets(ref)
    assert sorted(leaf_order(ours)) == list(range(n))


def _merged_sets(z: np.ndarray) -> list[frozenset[int]]:
    n = z.shape[0] + 1
    members: dict[int, frozenset[int]] = {i: frozenset([i]) for i in range(n)}
    out = []
    for k, (a, b, _h, _s) in enumerate(z):
        members[n + k] = members[int(a)] | members[int(b)]
        out.append(members[n + k])
    return sorted(out, key=sorted)


def test_average_linkage_with_ties_and_identical_genomes() -> None:
    # Clonal sets have many pairs at exactly the same distance (often 0).
    dist = np.zeros((5, 5))
    z = average_linkage(dist)
    assert z[:, 2].tolist() == [0.0] * 4
    assert z[-1, 3] == 5
    assert sorted(leaf_order(z)) == list(range(5))


def test_dendrogram_pdf_keeps_close_genomes_adjacent(tmp_path: Path) -> None:
    from matplotlib.figure import Figure

    names, sim = _two_groups()
    seen: dict[str, str] = {}
    real_savefig = Figure.savefig

    def _capture(self, fname, *args, **kwargs):
        seen["xlabel"] = self.axes[0].get_xlabel()
        return real_savefig(self, fname, *args, **kwargs)

    pdf = tmp_path / "d.pdf"
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Figure, "savefig", _capture)
        leaves = write_dendrogram(names, sim, pdf, measure="ANI")
    assert pdf.read_bytes().startswith(b"%PDF")
    assert sorted(leaves) == sorted(names)
    pos = {name: i for i, name in enumerate(leaves)}
    assert abs(pos["a.fasta"] - pos["b.fasta"]) == 1
    assert abs(pos["c.fasta"] - pos["d.fasta"]) == 1
    assert seen["xlabel"].startswith("1 - ANI")


def test_dendrogram_of_many_genomes_caps_the_page_height(tmp_path: Path) -> None:
    from matplotlib.figure import Figure

    n = 2000
    rng = np.random.default_rng(0)
    sim = 1.0 - rng.random((n, n)) * 0.05
    sim = (sim + sim.T) / 2
    np.fill_diagonal(sim, 1.0)
    heights: list[float] = []

    def _capture(self, fname, *args, **kwargs):
        heights.append(self.get_size_inches()[1])

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Figure, "savefig", _capture)
        leaves = write_dendrogram([f"g{i}" for i in range(n)], sim, tmp_path / "d.pdf")
    assert len(leaves) == n
    assert heights == [200.0]
