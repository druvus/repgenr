"""Parameter hygiene from the CLI audit (PR-D)."""

from __future__ import annotations

import dataclasses
import importlib
import json
from pathlib import Path

from benchmarks.genomegen import generate_set
from repgenr.cli.param_builders import tree2tax_params
from repgenr.core.contracts import MSA_FASTA
from repgenr.stages.glance import GlanceParams
from repgenr.stages.tree2tax import Tree2taxParams


def test_tree2tax_params_has_no_dead_all_genomes_field() -> None:
    assert "all_genomes" not in {f.name for f in dataclasses.fields(Tree2taxParams)}
    assert not hasattr(tree2tax_params(), "all_genomes")


def test_glance_threads_default_matches_cli() -> None:
    from repgenr.cli.base import DEFAULT_THREADS

    assert GlanceParams().threads == DEFAULT_THREADS


def test_aligners_write_the_contract_msa_name() -> None:
    for mod in ("cactus", "sibeliaz", "progressivemauve"):
        source = Path(importlib.import_module(f"repgenr.aligners.{mod}").__file__).read_text(
            encoding="utf-8"
        )
        assert "out_dir / MSA_FASTA" in source, mod
    assert MSA_FASTA == "msa.fasta"


def test_genomegen_small_n_yields_exactly_n_genomes(tmp_path: Path) -> None:
    for scenario in ("balanced", "clonal", "mixed"):
        for n in (1, 2, 3, 5, 8):
            out = tmp_path / f"{scenario}_{n}"
            truth = generate_set(out, scenario=scenario, n=n, seed=1, genome_length=2000)
            assert len(truth["clusters"]) == n
            assert len(list(out.glob("*.fasta"))) == n
            assert all(v for v in json.loads((out / "truth.json").read_text())["clusters"])


def test_run_bench_storage_from_env(monkeypatch, tmp_path: Path) -> None:
    import benchmarks.run_bench as rb

    monkeypatch.setenv("REPGENR_BENCH_STORAGE", str(tmp_path))
    rb = importlib.reload(rb)
    assert rb.STORAGE == tmp_path
    monkeypatch.delenv("REPGENR_BENCH_STORAGE")
    rb = importlib.reload(rb)
    assert rb.STORAGE == rb.DEFAULT_STORAGE
