"""Regression tests pinning known bias mechanisms (scaling/bias audit).

These document behavior found by the 2026-08 audit rather than assert it is
desirable: the --target-reps search on a clonal (step-function) dataset, and
the taxonomy-reduce keeper preferring the most over-represented genotype.
"""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    DerepParams,
    DerepResult,
)
from repgenr.stages.dereplicate import (
    DereplicateParams,
    _reduce_by_taxonomy,
    _search_target_reps,
)

_LOGGER = logging.getLogger("bias-regressions")


class _StepFunctionDerep:
    """A clonal dataset's rep-count curve: 1 rep below the clone ANI, n at/above.

    Real clone blocks sit at ~99.99 percent ANI, so representative count jumps
    discontinuously as --secondary-ani crosses that value; every count between
    1 and n is unreachable by threshold choice.
    """

    def __init__(self, n: int, clone_ani: float = 0.9995) -> None:
        self.n = n
        self.clone_ani = clone_ani
        self.calls: list[float] = []

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:  # noqa: ANN001
        genomes = list(genomes)
        self.calls.append(params.secondary_ani)
        keep = len(genomes) if params.secondary_ani >= self.clone_ani else 1
        reps = genomes[:keep]
        leftover = [g.name for g in genomes[keep:]]
        clusters: dict[str, list[str]] = {r.name: [] for r in reps}
        clusters[reps[0].name] = leftover
        status = {r.name: STATUS_REPRESENTATIVE for r in reps}
        for name in leftover:
            status[name] = STATUS_CONTAINED
        return DerepResult(representatives=reps, clusters=clusters, genome_status=status)


def _genomes(tmp_path: Path, n: int) -> list[Path]:
    out = []
    for i in range(n):
        p = tmp_path / f"g{i:03d}.fasta"
        p.write_text(">x\nACGT\n", encoding="utf-8")
        out.append(p)
    return out


def test_target_reps_step_function_cannot_reach_intermediate_target(tmp_path):
    # Target 20 of 40 clonal genomes: only 1 or 40 are achievable, and the
    # search silently returns the closest count instead of warning.
    adapter = _StepFunctionDerep(n=40)
    genomes = _genomes(tmp_path, 40)
    result = _search_target_reps(
        adapter, genomes, tmp_path / "scratch", DerepParams(),
        DereplicateParams(tool="mock"), target=20, logger=_LOGGER,
    )
    achieved = len(result.representatives)
    assert achieved in (1, 40)
    assert achieved != 20  # the requested target is unreachable on a step function
    assert len(adapter.calls) <= 12  # bounded by _MAX_TARGET_ITERS


def test_target_reps_exact_target_short_circuits(tmp_path):
    # A reachable target (the step's top) is found and the search stops early.
    adapter = _StepFunctionDerep(n=40)
    genomes = _genomes(tmp_path, 40)
    result = _search_target_reps(
        adapter, genomes, tmp_path / "scratch", DerepParams(),
        DereplicateParams(tool="mock"), target=40, logger=_LOGGER,
    )
    assert len(result.representatives) == 40


class _FakeManifest:
    def __init__(self, species_of: dict[str, str]) -> None:
        self._records = [
            SimpleNamespace(filename=name, species=sp, genus=sp.split(" ")[0])
            for name, sp in species_of.items()
        ]

    def all_genomes(self, include_outgroup: bool = True):  # noqa: ANN001
        return self._records


def test_taxonomy_reduce_keeper_is_the_overrepresented_genotype(tmp_path):
    # Two same-species ANI representatives: a clone-block rep containing 30
    # genomes and a diverse rep containing 2. The keeper rule (largest existing
    # cluster) always keeps the over-represented genotype and folds the diverse
    # lineage under it.
    clone_rep = tmp_path / "clone_rep.fasta"
    diverse_rep = tmp_path / "diverse_rep.fasta"
    for p in (clone_rep, diverse_rep):
        p.write_text(">x\nACGT\n", encoding="utf-8")
    clone_members = [f"clone_{i:02d}.fasta" for i in range(30)]
    diverse_members = ["div_a.fasta", "div_b.fasta"]
    result = DerepResult(
        representatives=[clone_rep, diverse_rep],
        clusters={"clone_rep.fasta": clone_members, "diverse_rep.fasta": diverse_members},
        genome_status={
            **{n: STATUS_CONTAINED for n in clone_members + diverse_members},
            "clone_rep.fasta": STATUS_REPRESENTATIVE,
            "diverse_rep.fasta": STATUS_REPRESENTATIVE,
        },
    )
    species = {
        name: "Benchgen benchsp"
        for name in ["clone_rep.fasta", "diverse_rep.fasta", *clone_members, *diverse_members]
    }
    ctx = SimpleNamespace(manifest=_FakeManifest(species))

    reduced = _reduce_by_taxonomy(ctx, result, "species", "tool", _LOGGER)

    assert [r.name for r in reduced.representatives] == ["clone_rep.fasta"]
    assert "diverse_rep.fasta" in reduced.clusters["clone_rep.fasta"]
    assert reduced.genome_status["diverse_rep.fasta"] == STATUS_CONTAINED


def test_taxonomy_reduce_keeps_unannotated_genomes(tmp_path):
    # A representative with no manifest taxon must stay its own group.
    rep = tmp_path / "unknown_rep.fasta"
    rep.write_text(">x\nACGT\n", encoding="utf-8")
    result = DerepResult(
        representatives=[rep],
        clusters={"unknown_rep.fasta": []},
        genome_status={"unknown_rep.fasta": STATUS_REPRESENTATIVE},
    )
    ctx = SimpleNamespace(manifest=_FakeManifest({}))
    reduced = _reduce_by_taxonomy(ctx, result, "species", "tool", _LOGGER)
    assert [r.name for r in reduced.representatives] == ["unknown_rep.fasta"]
