# SWOT Remediation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the verified weaknesses from the 2026-09-01 code audit in four independently shippable phases: silent-failure fixes, quality-aware representative selection, a reference-free SNP typer, and release hygiene.

**Architecture:** Every fix lands at the contract layer (stage or shared helper), never inside one adapter, so all tools in a family benefit at once. New behaviour is opt-out by flag with the safer choice as default; provenance records the choice. Each phase is one pull request with green CI before the next starts.

**Tech Stack:** Python 3.12, Typer CLI, pytest with fake adapters registered on the real registries, ruff 0.15.18, mypy, Nextflow DSL2 with nf-test 0.9.5.

**Spec:** The audit report (artifact `https://claude.ai/code/artifact/ae44ec08-749a-41ab-8d10-08c2b868037f`) and its companions `docs/swot-derep.md`, `docs/swot-phylo.md`, `docs/swot-viral.md`, `docs/scaling-audit.md`.

## Global Constraints

- Python floor `>=3.12`; CI matrix 3.12 and 3.13.
- Ruff is pinned to `0.15.18`; run `ruff check src/ tests/` and `mypy src/repgenr` before every commit.
- Coverage gate is `--cov-fail-under=85`; do not lower it.
- No Unicode in any file under `nextflow/`.
- Modest scientific language in docs, comments and log messages.
- External tools are never executed in Python tests; mock `run_tool` or register a fake adapter with the `register_tool` fixture from `tests/conftest.py`.
- Every deliverable file is written through `atomic_path` or `atomic_replace` from `repgenr.core.contracts`.
- Commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.
- Run the quality gates from the `repgenr_dev` conda env: `/Users/andreassjodin/miniforge3/envs/repgenr_dev/bin/python -m pytest -q`.
- Assumption: the next release is `3.0.0` because the Unreleased changelog section already contains breaking changes and the project follows SemVer. If the maintainer prefers `2.1.0`, change the number in Task 13 only.

## Phases and pull requests

| Phase | Tasks | Branch | Ships |
|---|---|---|---|
| A. Silent-failure fixes | 1 to 6 | `fix/silent-failures` | completeness assertion, extras gating, masker input, read-only doctor, dead flank path, Nextflow versions |
| B. Quality-aware keeper | 7 to 9 | `feat/quality-keeper` | GTDB quality columns in manifest and selection, keeper rescoring, taxonomy reduce keeper |
| C. Reference-free typer | 10 | `feat/ska2-snptyper` | ska2 adapter |
| D. Release and hygiene | 11 to 14 | `chore/release-3.0.0` | version test, wiki cleanup, CI format gate and caching, changelog and tag |

Not in this plan (separate plans after brainstorming): nf-core meta maps and `ext.args` in Nextflow, quality-aware `--limit` sampling, tree2tax support-based collapse, benchmark path parameterisation.

---

## Phase A: Silent-failure fixes

### Task 1: Completeness assertion after every dereplication

**Files:**
- Modify: `src/repgenr/dereplicators/base.py` (add `check_result_complete`)
- Modify: `src/repgenr/stages/dereplicate.py:469` (`_write_contract`)
- Modify: `src/repgenr/stages/derep_steps.py:97,143` (chunk and merge)
- Test: `tests/unit/test_derep_result_complete.py`

**Interfaces:**
- Produces: `check_result_complete(result: DerepResult, genome_names: Collection[str]) -> None` raising `WorkdirError`.

- [ ] **Step 1: Write the failing tests**

```python
"""Every input genome must leave dereplication with a status and a home."""

from __future__ import annotations

from pathlib import Path

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    DerepResult,
    check_result_complete,
)


def _ok() -> DerepResult:
    return DerepResult(
        representatives=[Path("a.fasta")],
        clusters={"a.fasta": ["b.fasta"]},
        genome_status={"a.fasta": STATUS_REPRESENTATIVE, "b.fasta": STATUS_CONTAINED},
    )


def test_complete_result_passes() -> None:
    check_result_complete(_ok(), ["a.fasta", "b.fasta"])


def test_missing_status_raises() -> None:
    with pytest.raises(WorkdirError, match="c.fasta"):
        check_result_complete(_ok(), ["a.fasta", "b.fasta", "c.fasta"])


def test_contained_without_cluster_raises() -> None:
    r = _ok()
    r.clusters = {"a.fasta": []}
    with pytest.raises(WorkdirError, match="no representative"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_representative_without_status_raises() -> None:
    r = _ok()
    r.genome_status = {"b.fasta": STATUS_CONTAINED}
    with pytest.raises(WorkdirError, match="a.fasta"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_unknown_status_value_raises() -> None:
    r = _ok()
    r.genome_status["b.fasta"] = "weird"
    with pytest.raises(WorkdirError, match="weird"):
        check_result_complete(r, ["a.fasta", "b.fasta"])


def test_fail_qc_needs_no_cluster() -> None:
    r = DerepResult(
        representatives=[Path("a.fasta")],
        clusters={"a.fasta": []},
        genome_status={"a.fasta": STATUS_REPRESENTATIVE, "b.fasta": "fail_qc"},
    )
    check_result_complete(r, ["a.fasta", "b.fasta"])
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/unit/test_derep_result_complete.py -q`
Expected: ImportError on `check_result_complete`.

- [ ] **Step 3: Implement in `src/repgenr/dereplicators/base.py`**

Add after the `DerepResult` dataclass:

```python
_VALID_STATUS = frozenset({STATUS_REPRESENTATIVE, STATUS_CONTAINED, STATUS_FAIL_QC})


def check_result_complete(result: DerepResult, genome_names: Collection[str]) -> None:
    """Refuse a result that silently drops genomes.

    Every input genome must carry a known status; every representative must be
    marked as such; every contained genome must belong to exactly one cluster.
    Adapters that parse tool tables can otherwise return an empty membership
    (e.g. after an upstream column-layout change) and the stage would complete
    with genomes missing from every deliverable.
    """
    from ..core.errors import WorkdirError

    names = set(genome_names)
    status = result.genome_status
    missing = sorted(names - status.keys())
    if missing:
        raise WorkdirError(
            f"Dereplication left {len(missing)} of {len(names)} genome(s) without a "
            f"status (e.g. {', '.join(missing[:3])}). The adapter output is incomplete."
        )
    bad = sorted(f"{g}={s}" for g, s in status.items() if s not in _VALID_STATUS)
    if bad:
        raise WorkdirError(f"Unknown dereplication status value(s): {', '.join(bad[:3])}")

    rep_names = {p.name for p in result.representatives}
    unmarked = sorted(r for r in rep_names if status.get(r) != STATUS_REPRESENTATIVE)
    if unmarked:
        raise WorkdirError(
            f"{len(unmarked)} representative(s) are not marked as such in genome_status "
            f"(e.g. {', '.join(unmarked[:3])})."
        )

    home: dict[str, int] = {}
    for rep, members in result.clusters.items():
        for m in members:
            if m != rep:
                home[m] = home.get(m, 0) + 1
    orphans = sorted(
        g for g, s in status.items() if s == STATUS_CONTAINED and home.get(g, 0) == 0
    )
    if orphans:
        raise WorkdirError(
            f"{len(orphans)} contained genome(s) have no representative "
            f"(e.g. {', '.join(orphans[:3])}). The adapter returned an empty cluster table."
        )
    doubled = sorted(g for g, n in home.items() if n > 1)
    if doubled:
        raise WorkdirError(
            f"{len(doubled)} genome(s) appear in more than one cluster "
            f"(e.g. {', '.join(doubled[:3])})."
        )
```

Add `Collection` to the `collections.abc` import at the top of the file.

- [ ] **Step 4: Wire into the three contract writers**

In `src/repgenr/stages/dereplicate.py` `run()`, immediately before `_write_contract(ctx, result)`:

```python
    check_result_complete(result, [g.name for g in genomes])
```

Import `check_result_complete` from `..dereplicators.base` at the top of the module. Note `_reduce_by_taxonomy` runs before this line, so the reduced result is what gets checked.

In `src/repgenr/stages/derep_steps.py` `dereplicate_chunk()`, before `_write_step_contract(...)`:

```python
    check_result_complete(result, [g.name for g in params.genomes])
```

In `dereplicate_merge()`, after the composed result is built and before its `_write_step_contract`, check against every stage-1 genome:

```python
    all_names = [
        name
        for r in stage1
        for rep, members in r.clusters.items()
        for name in (rep, *members)
    ]
    check_result_complete(composed, all_names)
```

Use the variable name that holds the `_compose_two_stage` return value in that function.

- [ ] **Step 5: Run the full suite**

Run: `pytest -q`
Expected: all pass. If a fake adapter in an existing test returns an incomplete result (e.g. `_load_chunk` returns `genome_status={}` but that is never validated directly), fix the fake, not the check.

- [ ] **Step 6: Commit**

```bash
git add src/repgenr/dereplicators/base.py src/repgenr/stages/dereplicate.py src/repgenr/stages/derep_steps.py tests/unit/test_derep_result_complete.py
git commit -m "fix(derep): refuse results that drop genomes

Every input genome must leave dereplication with a known status and,
if contained, exactly one representative. Turns the skder empty-edge
and sourmash unmatched-label cases into hard failures."
```

---

### Task 2: Gate stage extras by `accepted_extras` across all families

**Files:**
- Modify: `src/repgenr/core/plugins.py` (add `warn_unconsumed_extras`)
- Modify: `src/repgenr/cli/base.py` (add `gated_extra`)
- Modify: `src/repgenr/cli/cmd_run.py:32-46` (replace `_virus_extra` body with `gated_extra`)
- Modify: `src/repgenr/cli/cmd_bacterial.py:122-125`, `src/repgenr/cli/cmd_steps.py:93-96,246-249`
- Modify: `src/repgenr/stages/dereplicate.py:94-100`, `src/repgenr/stages/derep_steps.py:91-97,123-129`, `src/repgenr/stages/phylo.py:110-114,299-315`, `src/repgenr/stages/snptype.py:100-105`
- Modify: `src/repgenr/aligners/progressivemauve.py`, `src/repgenr/aligners/sibeliaz.py`, `src/repgenr/treebuilders/mashtree.py` (declare `accepted_extras`)
- Test: `tests/unit/test_extras_gating.py`

**Interfaces:**
- Produces: `warn_unconsumed_extras(caps: ToolCapabilities, extra: Mapping[str, object], logger, *, family: str) -> None`
- Produces: `gated_extra(registry, tool: str, key: str, value: object) -> dict` in `cli/base.py`
- Consumes: `ToolCapabilities.accepted_extras` (existing)

- [ ] **Step 1: Write the failing tests**

```python
"""Extras reach adapters only when they read them; unread keys are named."""

from __future__ import annotations

import logging

from repgenr.cli.base import gated_extra
from repgenr.core.plugins import Registry, ToolCapabilities, warn_unconsumed_extras
from repgenr.dereplicators.base import Dereplicator


class _Reads(Dereplicator):
    capabilities = ToolCapabilities(name="reads", accepted_extras=frozenset({"virus"}))

    def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
        raise NotImplementedError


class _Ignores(Dereplicator):
    capabilities = ToolCapabilities(name="ignores")

    def dereplicate(self, genomes, out_dir, params, logger):  # noqa: ANN001
        raise NotImplementedError


def _registry(register_tool):
    reg: Registry = Registry("repgenr.dereplicators")
    register_tool(reg, "reads", _Reads)
    register_tool(reg, "ignores", _Ignores)
    return reg


def test_gated_extra_injects_for_reader(register_tool) -> None:
    reg = _registry(register_tool)
    assert gated_extra(reg, "reads", "virus", True) == {"virus": True}


def test_gated_extra_skips_for_ignorer(register_tool) -> None:
    reg = _registry(register_tool)
    assert gated_extra(reg, "ignores", "virus", True) == {}


def test_gated_extra_passes_through_for_auto(register_tool) -> None:
    reg = _registry(register_tool)
    assert gated_extra(reg, "auto", "virus", True) == {"virus": True}


def test_warn_names_unread_keys(caplog) -> None:
    caps = ToolCapabilities(name="t", accepted_extras=frozenset({"a"}))
    with caplog.at_level(logging.WARNING):
        warn_unconsumed_extras(caps, {"a": 1, "b": 2, "c": 3}, logging.getLogger("x"), family="Aligner")
    msgs = [r.message for r in caplog.records]
    assert any("Aligner 't' ignores" in m and "b, c" in m for m in msgs)


def test_warn_silent_when_all_read(caplog) -> None:
    caps = ToolCapabilities(name="t", accepted_extras=frozenset({"a"}))
    with caplog.at_level(logging.WARNING):
        warn_unconsumed_extras(caps, {"a": 1}, logging.getLogger("x"), family="Aligner")
    assert not caplog.records
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/unit/test_extras_gating.py -q`
Expected: ImportError.

- [ ] **Step 3: Add the two helpers**

In `src/repgenr/core/plugins.py`, after `ToolCapabilities`:

```python
def warn_unconsumed_extras(
    caps: ToolCapabilities,
    extra: Mapping[str, object],
    logger: logging.Logger,
    *,
    family: str,
) -> None:
    """Warn once, by name, about extra keys the adapter never reads.

    Keys that are also in ``default_params`` are the adapter's own defaults and
    are not reported.
    """
    unread = sorted(set(extra) - caps.accepted_extras - set(caps.default_params))
    if unread:
        logger.warning(
            "%s '%s' ignores extra parameter(s): %s", family, caps.name, ", ".join(unread)
        )
```

In `src/repgenr/cli/base.py`, near `_parse_key_values`:

```python
def gated_extra(registry, tool: str, key: str, value: object) -> dict:
    """Return ``{key: value}`` only when ``tool`` reads that extra.

    Injecting a key a tool ignores would change the resume fingerprint without
    changing the result. ``auto`` passes the key through; the stage warns after
    it has picked a concrete tool.
    """
    if tool != "auto":
        caps = registry.get(tool).capabilities
        if key not in caps.accepted_extras:
            return {}
    return {key: value}
```

- [ ] **Step 4: Use them at every injection site**

`cli/cmd_run.py`: replace the body of `_virus_extra` with

```python
    if not viral:
        return {}
    from ..dereplicators.base import registry as _derep_registry

    return gated_extra(_derep_registry, derep_tool, "virus", True)
```

and import `gated_extra` from `.base`.

`cli/cmd_bacterial.py` `dereplicate`: replace `**({"virus": True} if virus else {})` with `**(gated_extra(_derep_registry, tool, "virus", True) if virus else {})`, importing the registry inside `build()` as `from ..dereplicators.base import registry as _derep_registry`.

`cli/cmd_steps.py`: same replacement at both sites; `_derep_registry` is already imported there.

`stages/dereplicate.py:94-100`: replace the inline `unconsumed` block with `warn_unconsumed_extras(caps, params.extra or {}, logger, family="Dereplicator")`.

`stages/derep_steps.py`: in both `dereplicate_chunk` and `dereplicate_merge`, after `adapter = registry.create(params.tool)` add

```python
    caps = adapter.capabilities
    warn_unconsumed_extras(caps, params.extra or {}, logger, family="Dereplicator")
```

and build `DerepParams(... extra={**caps.default_params, **(params.extra or {})})` so the step path merges defaults the same way the in-process path does.

`stages/phylo.py`: define at module level

```python
# Keys the phylo stage reads itself; never forwarded to an adapter.
_STAGE_EXTRA_KEYS = frozenset({"mask"})


def _adapter_extra(extra: dict) -> dict:
    return {k: v for k, v in extra.items() if k not in _STAGE_EXTRA_KEYS}
```

At line 110-114 (tree builder) use `extra={**builder.capabilities.default_params, **_adapter_extra(params.extra)}` and call `warn_unconsumed_extras(builder.capabilities, _adapter_extra(params.extra), logger, family="Tree builder")` before it. At line 299-301 (aligner) use `extra=_adapter_extra(params.extra)` and warn with `family="Aligner"`. At line 305-315 (snptype branch) add `extra=_adapter_extra(params.extra)` to the `SnptypeParams(...)` call so aligner-arg tuning reaches the typer.

`stages/snptype.py` after `typer = snp_registry.create(params.tool)`: `warn_unconsumed_extras(typer.capabilities, params.extra, logger, family="SNP typer")`.

- [ ] **Step 5: Declare what adapters actually read**

- `aligners/progressivemauve.py`: `accepted_extras=frozenset({"seed_weight"})`
- `aligners/sibeliaz.py`: `accepted_extras=frozenset({"filtermemory", "kmer", "abundance", "bubble"})`
- `treebuilders/mashtree.py`: `accepted_extras=frozenset({"genomesize"})`

Verify with `grep -n "extra" src/repgenr/aligners/*.py src/repgenr/treebuilders/*.py src/repgenr/snptypers/*.py` that no other key is read.

- [ ] **Step 6: Run the gates**

Run: `pytest -q && ruff check src/ tests/ && mypy src/repgenr`
Expected: pass. `tests/unit/test_cli_wiring.py` may assert the old `{"virus": True}` in extras for skder; update that assertion to expect `{}` and add the same test for dRep expecting `{"virus": True}`.

- [ ] **Step 7: Update `docs/adding-tools.md:61-64`**

Change the sentence claiming only the dereplicate stage warns to: "Every stage warns, by name, about extra keys the chosen adapter does not declare in `accepted_extras`, and the CLI does not inject stage-level flags (such as `--virus`) into the extras of a tool that does not read them."

- [ ] **Step 8: Commit**

```bash
git add -A src/repgenr tests/unit/test_extras_gating.py tests/unit/test_cli_wiring.py docs/adding-tools.md
git commit -m "fix: gate extras by accepted_extras in every family

Stage-level flags no longer alter fingerprints of tools that ignore
them; phylo, snptype and the Nextflow step path now warn about unread
keys the way dereplicate already did."
```

---

### Task 3: Give maskers a whole-genome alignment

**Files:**
- Modify: `src/repgenr/snptypers/base.py` (`SnpResult.full_alignment`, drop `SnpParams.mask`)
- Modify: `src/repgenr/maskers/base.py` (`MaskParams`, new signature)
- Modify: `src/repgenr/maskers/gubbins.py`
- Modify: `src/repgenr/snptypers/snippy.py`, `src/repgenr/snptypers/parsnp.py`, `src/repgenr/snptypers/simple.py`
- Modify: `src/repgenr/stages/snptype.py:100-121`
- Test: `tests/integration/test_snptype_stage.py`, `tests/unit/test_gubbins_masker.py`

**Interfaces:**
- Produces: `SnpResult.full_alignment: Path | None` (whole-genome alignment in reference coordinates, or None when the typer cannot provide one).
- Produces: `MaskParams(threads: int = 16)`; `Masker.mask(full_alignment: Path, out_dir: Path, params: MaskParams, logger) -> Path` returning the masked variable-site FASTA.

- [ ] **Step 1: Write the failing tests**

Append to `tests/integration/test_snptype_stage.py`:

```python
from repgenr.core.errors import UserInputError
from repgenr.maskers.base import MaskParams, Masker
from repgenr.maskers.base import registry as masker_registry


class _FullTyper(_FakeTyper):
    def call(self, genomes, reference, out_dir, params, logger) -> SnpResult:  # noqa: ANN001
        core = out_dir / "core.fasta"
        full = out_dir / "full.fasta"
        core.write_text("".join(f">{g.stem}\nACGT\n" for g in genomes))
        full.write_text("".join(f">{g.stem}\nACGTACGTACGT\n" for g in genomes))
        return SnpResult(core_snp_fasta=core, full_alignment=full)


class _RecordingMasker(Masker):
    capabilities = ToolCapabilities(name="fakemask")
    seen: dict = {}

    def preflight(self) -> dict[str, str]:
        return {"fakemask": "1"}

    def mask(self, full_alignment, out_dir, params, logger):  # noqa: ANN001
        _RecordingMasker.seen = {"input": full_alignment, "threads": params.threads}
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / "masked.fasta"
        out.write_text(">a\nAC\n")
        return out


@pytest.fixture
def fake_masker(register_tool):
    register_tool(masker_registry, "fakemask", _RecordingMasker)
    yield


def test_masker_receives_full_alignment(workdir, genome_files, register_tool, fake_masker):
    register_tool(registry, "fulltyper", _FullTyper)
    ctx = WorkdirContext(workdir, create=True)
    run(ctx, SnptypeParams(tool="fulltyper", all_genomes=True, mask="fakemask", threads=3))
    assert _RecordingMasker.seen["input"].name == "full.fasta"
    assert _RecordingMasker.seen["threads"] == 3
    assert (ctx.snp_dir / CORE_SNP_FASTA).read_text() == ">a\nAC\n"


def test_mask_refused_without_full_alignment(workdir, genome_files, fake_typer, fake_masker):
    ctx = WorkdirContext(workdir, create=True)
    with pytest.raises(UserInputError, match="whole-genome alignment"):
        run(ctx, SnptypeParams(tool="faketyper", all_genomes=True, mask="fakemask"))
```

Create `tests/unit/test_gubbins_masker.py`:

```python
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
        Path(str(kw["cwd"] / "gubbins") + ".filtered_polymorphic_sites.fasta").write_text(">a\nA\n")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    full = tmp_path / "full.fasta"
    full.write_text(">a\nACGT\n")
    out = mod.GubbinsMasker().mask(full, tmp_path / "gub", MaskParams(threads=4), logging.getLogger("t"))
    assert out.exists()
    argv = calls[0]
    assert argv[0] == "run_gubbins.py"
    assert argv[argv.index("--threads") + 1] == "4"
    assert argv[-1] == str(full)
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/integration/test_snptype_stage.py tests/unit/test_gubbins_masker.py -q`
Expected: TypeError on `full_alignment`, ImportError on `MaskParams`.

- [ ] **Step 3: Change the contracts**

`snptypers/base.py`:

```python
@dataclass
class SnpParams:
    threads: int = 16
    reference: Path | None = None
    extra: dict = field(default_factory=dict)


@dataclass
class SnpResult:
    core_snp_fasta: Path
    vcf: Path | None = None
    snp_distance_matrix: Path | None = None
    masked: bool = False
    # Whole-genome alignment in reference coordinates (one record per genome,
    # equal lengths). Required input for recombination maskers; None when the
    # typer only produces variable sites.
    full_alignment: Path | None = None
```

`maskers/base.py`:

```python
@dataclass
class MaskParams:
    threads: int = 16


class Masker(ABC):
    capabilities: ToolCapabilities

    def preflight(self) -> dict[str, str]:
        return preflight(self.capabilities)

    @abstractmethod
    def mask(
        self,
        full_alignment: Path,
        out_dir: Path,
        params: MaskParams,
        logger: logging.Logger,
    ) -> Path:
        """Mask recombinant regions of ``full_alignment``; return the masked
        variable-site FASTA that replaces the typer's core-SNP alignment."""
        raise NotImplementedError
```

Update the module docstring: a masker takes the whole-genome alignment, not the SNP alignment. Add `from dataclasses import dataclass` import.

`maskers/gubbins.py` `mask` signature per the interface; argv becomes `["run_gubbins.py", "--threads", str(params.threads), "--prefix", prefix, full_alignment]`.

- [ ] **Step 4: Make the typers provide the full alignment**

`snptypers/snippy.py`: after `core_aln` check add

```python
        full_aln = Path(str(core_prefix) + ".full.aln")
        full: Path | None = None
        if full_aln.exists():
            full = out_dir / "full_alignment.fasta"
            full.write_text(full_aln.read_text(encoding="utf-8"), encoding="utf-8")
        ...
        return SnpResult(core_snp_fasta=core_fasta, masked=False, full_alignment=full)
```

`snptypers/parsnp.py`: after the `-S` harvesttools call, add a second call `["harvesttools", "-i", ggr, "-M", full_fasta]` with `full_fasta = out_dir / "full_alignment.fasta"` and return it as `full_alignment`.

`snptypers/simple.py`: the `consensuses` dict already holds full-length reference-coordinate sequences. Write them with

```python
        full_fasta = out_dir / "full_alignment.fasta"
        with open(full_fasta, "w", encoding="utf-8") as fo:
            for name, seq in consensuses.items():
                fo.write(f">{name}\n{seq}\n")
```

before `_write_core_snps`, and return `full_alignment=full_fasta`.

- [ ] **Step 5: Rewire the stage**

In `stages/snptype.py` `snptype_core`, remove `mask=params.mask` from the `SnpParams(...)` call and replace the masking block with:

```python
    if params.mask not in ("none", ""):
        from ..maskers.base import MaskParams
        from ..maskers.base import registry as masker_registry

        if result.full_alignment is None:
            raise UserInputError(
                f"--mask {params.mask} needs a whole-genome alignment, which SNP typer "
                f"'{params.tool}' does not produce. Use snippy, parsnp or simple, or drop --mask."
            )
        masker = masker_registry.create(params.mask)
        versions.update(masker.preflight())
        filtered = masker.mask(
            result.full_alignment, scratch / params.mask,
            MaskParams(threads=params.threads), logger,
        )
        with atomic_path(core) as tmp:
            shutil.copy2(filtered, tmp)
        masked = True
```

Also copy the full alignment into the deliverable set when present: `snp/full_alignment.fasta` via `atomic_path`, and add `full_alignment=` to the returned `SnpResult`.

- [ ] **Step 6: Run the gates and fix fakes**

Run: `pytest -q && ruff check src/ tests/ && mypy src/repgenr`
Expected: pass. Any test constructing `SnpParams(mask=...)` must drop the argument.

- [ ] **Step 7: Document**

In `docs/usage.md`, under the snptype command, add: "Recombination masking (`--mask gubbins`) runs on the typer's whole-genome alignment and replaces the core-SNP alignment with Gubbins' filtered polymorphic sites. Typers that only emit variable sites cannot be masked." Update `docs/output.md` to list `snp/full_alignment.fasta`.

- [ ] **Step 8: Commit**

```bash
git add -A src/repgenr/snptypers src/repgenr/maskers src/repgenr/stages/snptype.py tests docs/usage.md docs/output.md
git commit -m "fix(snptype): mask on the whole-genome alignment

Gubbins now receives the full alignment (snippy core.full.aln,
harvesttools -M, or the simple typer's consensuses) and a thread
budget; masking a variable-site-only typer is refused with a clear
message. Drops the unread SnpParams.mask field."
```

---

### Task 4: Make `doctor` read-only

**Files:**
- Modify: `src/repgenr/core/manifest.py` (add `Manifest.open_readonly`)
- Modify: `src/repgenr/core/doctor.py:134,274`
- Test: `tests/unit/test_doctor.py`

**Interfaces:**
- Produces: `Manifest.open_readonly(path) -> Manifest` classmethod; raises `WorkdirError` if the file does not exist.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_doctor.py`:

```python
def test_doctor_leaves_manifest_untouched(workdir: Path) -> None:
    from repgenr.core.doctor import run_doctor
    from repgenr.core.manifest import GenomeRecord, Manifest

    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1", filename="x.fasta"))
    (workdir / "selection.tsv").write_text("accession\tfilename\nGCA_1\tx.fasta\n")
    manifest = workdir / "manifest.sqlite"
    before = (manifest.stat().st_mtime_ns, sorted(p.name for p in workdir.iterdir()))

    run_doctor(workdir)

    after = (manifest.stat().st_mtime_ns, sorted(p.name for p in workdir.iterdir()))
    assert before == after


def test_open_readonly_refuses_writes(workdir: Path) -> None:
    import sqlite3

    from repgenr.core.manifest import GenomeRecord, Manifest

    workdir.mkdir(parents=True)
    with Manifest.open(workdir) as m:
        m.upsert(GenomeRecord(accession="GCA_1"))
    ro = Manifest.open_readonly(workdir / "manifest.sqlite")
    try:
        assert [g.accession for g in ro.all_genomes()] == ["GCA_1"]
        with pytest.raises(sqlite3.OperationalError):
            ro.upsert(GenomeRecord(accession="GCA_2"))
    finally:
        ro.close()
```

Look up the actual name of the doctor entry function in `src/repgenr/core/doctor.py` (the function the CLI calls) and use it instead of `run_doctor` if it differs. Note that a WAL database opened read-only may create `-shm`; if the listing assertion fails only on that file, close the writing manifest with `PRAGMA journal_mode=DELETE` inside the test before measuring, and keep the mtime assertion.

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/unit/test_doctor.py -q -k "untouched or readonly"`
Expected: AttributeError `open_readonly`; the first test may fail on mtime.

- [ ] **Step 3: Implement**

In `manifest.py`, refactor `__init__` so the connection setup lives in a private constructor path:

```python
    def __init__(self, path: str | os.PathLike[str], *, readonly: bool = False):
        self.path = Path(path)
        if readonly:
            if not self.path.exists():
                raise WorkdirError(f"Manifest not found: {self.path}")
            self._conn = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True)
            self._conn.row_factory = sqlite3.Row
            self._conn.execute(f"PRAGMA busy_timeout={BUSY_TIMEOUT_MS}")
            self._check_version()
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        ... (existing body unchanged)

    def _check_version(self) -> None:
        version = int(self._conn.execute("PRAGMA user_version").fetchone()[0])
        if version > SCHEMA_VERSION:
            raise WorkdirError(...)  # same message as _migrate

    @classmethod
    def open_readonly(cls, path: str | os.PathLike[str]) -> Manifest:
        """Open an existing manifest without creating, migrating or journaling."""
        return cls(path, readonly=True)
```

Have `_migrate` call `_check_version` for the too-new case so the message is defined once.

In `doctor.py`, replace both `Manifest(manifest_path)` with `Manifest.open_readonly(manifest_path)`.

- [ ] **Step 4: Run the gates and commit**

Run: `pytest -q && ruff check src/ tests/ && mypy src/repgenr`

```bash
git add src/repgenr/core/manifest.py src/repgenr/core/doctor.py tests/unit/test_doctor.py
git commit -m "fix(doctor): open the manifest read-only

doctor promised no writes but the Manifest constructor created the
directory, switched to WAL and applied the schema."
```

---

### Task 5: Remove the dead flank path from the XMFA converter

**Files:**
- Modify: `src/repgenr/converters/xmfa_to_fasta.py:33-170`
- Modify: `src/repgenr/aligners/progressivemauve.py:85`
- Test: `tests/unit/test_xmfa_to_fasta.py`

- [ ] **Step 1: Update the tests to the new signature**

In `tests/unit/test_xmfa_to_fasta.py`, every call `xmfa_to_fasta(path, ref, 0, out)` becomes `xmfa_to_fasta(path, ref, out)`. Add:

```python
def test_flank_parameter_removed() -> None:
    import inspect

    from repgenr.converters.xmfa_to_fasta import xmfa_to_fasta

    assert "flank" not in inspect.signature(xmfa_to_fasta).parameters
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/unit/test_xmfa_to_fasta.py -q`
Expected: TypeError (wrong arity) and the signature assertion fails.

- [ ] **Step 3: Delete the path**

Remove the `flank` parameter and its docstring line; delete the `if flank > 0:` blocks (the gap-list maintenance inside the reference-gap loop at lines 128-133, the whole flank block at 134-152, and the `start += flank` / `end -= flank` / `seq[flank:-flank]` branches in the projection loop). `list_of_gaps` becomes unused and is removed too. Update the caller in `aligners/progressivemauve.py:85` to the three-argument form.

- [ ] **Step 4: Run the gates and commit**

Run: `pytest -q && ruff check src/ tests/ && mypy src/repgenr`

```bash
git add src/repgenr/converters/xmfa_to_fasta.py src/repgenr/aligners/progressivemauve.py tests/unit/test_xmfa_to_fasta.py
git commit -m "refactor(converters): drop the unused XMFA flank screening

The only caller passed flank=0; the untested path also scanned only the
first non-reference sequence for gaps."
```

---

### Task 6: Publish software versions from Nextflow; fix help text

**Files:**
- Modify: `nextflow/main.nf:30,41-46`
- Modify: `nextflow/modules/local/dataflow/tree2tax.nf:11`
- Modify: `nextflow/tests/local_dataflow.nf.test`
- Modify: `docs/architecture.md:32`, `docs/output.md`

- [ ] **Step 1: Add the assertion to the e2e nf-test**

In `nextflow/tests/local_dataflow.nf.test`, inside the existing `then` block, add:

```groovy
            assert path("${params.outdir}/pipeline_info/software_versions.yml").exists()
            assert path("${params.outdir}/pipeline_info/software_versions.yml").text.contains("repgenr")
```

Match the variable used by the neighbouring assertions for the output directory (`params.outdir` or `outputDir`).

- [ ] **Step 2: Run to verify failure**

Run: `nf-test test nextflow/tests/local_dataflow.nf.test --tag e2e` (needs the e2e conda env from `.github/workflows/ci.yml`; if unavailable locally, rely on CI for this task and state so in the PR).
Expected: assertion fails, file absent.

- [ ] **Step 3: Publish the versions channel**

In `nextflow/main.nf`, replace the workflow body after `paramsSummaryLog` with:

```groovy
    ch_versions = Channel.empty()
    if (params.mode == 'viral') {
        VIRAL_DATAFLOW()
        ch_versions = VIRAL_DATAFLOW.out.versions
    }
    else {
        BACTERIAL_DATAFLOW()
        ch_versions = BACTERIAL_DATAFLOW.out.versions
    }
    ch_versions
        .unique()
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
        )
```

Change the help line 30 to `--vmetadata_args / --vgenome_args '<str>'  NCBI Virus selection (viral)`.

- [ ] **Step 4: Keep versions.yml out of the outdir root**

In `tree2tax.nf` change `publishDir "${params.outdir}", mode: 'copy'` to `publishDir "${params.outdir}", mode: 'copy', pattern: '*.tsv'`.

- [ ] **Step 5: Docs**

`docs/architecture.md:32`: "Five families" and add `maskers` to the list. `docs/output.md`: add `pipeline_info/software_versions.yml` under the Nextflow output section.

- [ ] **Step 6: Verify and commit**

Run: `grep -rnP '[^\x00-\x7F]' nextflow` (expect no output) and `nf-test test --tag stub`.

```bash
git add nextflow docs/architecture.md docs/output.md
git commit -m "feat(nextflow): publish software_versions.yml

The per-process versions.yml fragments were collected and then dropped
by main.nf. Also corrects the viral help text and keeps versions.yml
out of the outdir root."
```

- [ ] **Step 7: Open the Phase A pull request**

Push `fix/silent-failures`, open a PR titled "Silent-failure fixes from the 2026-09-01 audit", wait for CI, merge on green.

---

## Phase B: Quality-aware representative selection

### Task 7: Carry GTDB quality into the manifest and selection

**Files:**
- Modify: `src/repgenr/core/manifest.py` (schema v2, `GenomeRecord.completeness/contamination`)
- Modify: `src/repgenr/core/contracts.py:145-186` (`SelectionRow`, `write_selection`, `read_selection`)
- Modify: `src/repgenr/stages/metadata.py:200-227,133,309-319`
- Test: `tests/unit/test_manifest.py`, `tests/unit/test_metadata_stage.py`, `tests/unit/test_selection_contract.py`

**Interfaces:**
- Produces: `GenomeRecord.completeness: float | None`, `GenomeRecord.contamination: float | None`; `SelectionRow` gains the same two optional fields; `selection.tsv` gains two trailing columns `completeness`, `contamination` (empty string when unknown).
- Produces: `Manifest.quality() -> dict[str, tuple[float, float]]` keyed by filename, only for genomes where both values are known.

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_manifest.py`:

```python
def test_quality_round_trip_and_migration(tmp_path):
    from repgenr.core.manifest import SCHEMA_VERSION, GenomeRecord, Manifest

    m = Manifest(tmp_path / "m.sqlite")
    m.upsert(GenomeRecord(accession="A", filename="a.fasta", completeness=99.1, contamination=0.4))
    m.upsert(GenomeRecord(accession="B", filename="b.fasta"))
    assert m.quality() == {"a.fasta": (99.1, 0.4)}
    assert SCHEMA_VERSION == 2
    m.close()


def test_v1_database_is_migrated(tmp_path):
    import sqlite3

    from repgenr.core.manifest import Manifest

    path = tmp_path / "old.sqlite"
    conn = sqlite3.connect(path)
    conn.executescript(
        "CREATE TABLE genomes (accession TEXT PRIMARY KEY, filename TEXT, source TEXT, "
        "family TEXT, genus TEXT, species TEXT, is_outgroup INTEGER DEFAULT 0, "
        "derep_status TEXT, representative TEXT); PRAGMA user_version=1;"
    )
    conn.execute("INSERT INTO genomes (accession, filename) VALUES ('A', 'a.fasta')")
    conn.commit()
    conn.close()
    m = Manifest(path)
    assert m.quality() == {}
    assert m.all_genomes()[0].completeness is None
    m.close()
```

`tests/unit/test_selection_contract.py`:

```python
def test_selection_quality_columns_round_trip(tmp_path):
    from repgenr.core.contracts import SelectionRow, read_selection, write_selection

    rows = [
        SelectionRow("A", "f", "g", "s", False, "f_g_s_A.fasta", completeness=98.5, contamination=1.0),
        SelectionRow("B", "f", "g", "s", True, "f_g_s_B.fasta"),
    ]
    write_selection(tmp_path / "s.tsv", rows)
    back = read_selection(tmp_path / "s.tsv")
    assert (back[0].completeness, back[0].contamination) == (98.5, 1.0)
    assert (back[1].completeness, back[1].contamination) == (None, None)


def test_selection_without_quality_columns_still_reads(tmp_path):
    from repgenr.core.contracts import read_selection

    p = tmp_path / "s.tsv"
    p.write_text("accession\tfamily\tgenus\tspecies\tis_outgroup\tfilename\nA\tf\tg\ts\t0\tx.fasta\n")
    assert read_selection(p)[0].completeness is None
```

`tests/unit/test_metadata_stage.py`: extend the existing GTDB TSV fixture header with `checkm2_completeness` and `checkm2_contamination` columns (fall back to `checkm_completeness`/`checkm_contamination` names in a second test) and assert the parsed accession dict carries `"completeness": 99.0, "contamination": 0.5` and that the resulting `GenomeRecord` has the values.

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/unit/test_manifest.py tests/unit/test_selection_contract.py tests/unit/test_metadata_stage.py -q`

- [ ] **Step 3: Manifest schema v2**

In `manifest.py`: `SCHEMA_VERSION = 2`; add `completeness REAL, contamination REAL` to `_SCHEMA`; add both to `_UPSERT_SQL` (columns, values and the `DO UPDATE SET` list); add fields `completeness: float | None = None` and `contamination: float | None = None` to `GenomeRecord`; extend `_record_params` and `_row_to_record` (use `row["completeness"] if "completeness" in row.keys() else None`). In `_migrate` implement the step:

```python
        if version < 2:
            cols = {r[1] for r in self._conn.execute("PRAGMA table_info(genomes)")}
            for col in ("completeness", "contamination"):
                if col not in cols:
                    self._conn.execute(f"ALTER TABLE genomes ADD COLUMN {col} REAL")
            version = 2
        self._conn.execute(f"PRAGMA user_version={SCHEMA_VERSION}")
```

Add:

```python
    def quality(self) -> dict[str, tuple[float, float]]:
        """filename -> (completeness, contamination) for genomes with both values."""
        rows = self._conn.execute(
            "SELECT filename, completeness, contamination FROM genomes "
            "WHERE filename IS NOT NULL AND completeness IS NOT NULL AND contamination IS NOT NULL"
        )
        return {r["filename"]: (float(r["completeness"]), float(r["contamination"])) for r in rows}
```

- [ ] **Step 4: Selection contract**

In `contracts.py`, add `completeness: float | None = None` and `contamination: float | None = None` to `SelectionRow`. `write_selection` header gains `completeness`, `contamination`; write `""` when None, else `f"{value:.2f}"`. `read_selection` parses with a helper `_opt_float(row.get("completeness"))` returning None for missing or empty.

- [ ] **Step 5: Metadata parsing**

In `_parse_metadata`, after computing `tax`:

```python
            completeness = _opt_column(fields, idx, "checkm2_completeness", "checkm_completeness")
            contamination = _opt_column(fields, idx, "checkm2_contamination", "checkm_contamination")
            accessions[accession] = {
                ..., "completeness": completeness, "contamination": contamination,
            }
```

with

```python
def _opt_column(fields: list[str], idx: dict[str, int], *names: str) -> float | None:
    """First present, non-empty, numeric column among ``names``; else None."""
    for name in names:
        i = idx.get(name)
        if i is None or i >= len(fields):
            continue
        raw = fields[i].strip()
        if not raw or raw.lower() in {"none", "na", "n/a"}:
            continue
        try:
            return float(raw)
        except ValueError:
            continue
    return None
```

Thread the values through `_record_from_tax(acc, data["tax"], completeness=data.get("completeness"), contamination=data.get("contamination"))` and into `SelectionRow(...)` in `_write_selection`. The API source path (`_api_row_to_tax` region) sets both to None unless the API row carries `checkm2_completeness`.

- [ ] **Step 6: Gates and commit**

Run: `pytest -q && ruff check src/ tests/ && mypy src/repgenr`

```bash
git add -A src/repgenr/core src/repgenr/stages/metadata.py tests
git commit -m "feat(metadata): carry GTDB CheckM quality into manifest and selection

Manifest schema v2 adds completeness/contamination with an ALTER TABLE
migration; selection.tsv gains two optional trailing columns."
```

---

### Task 8: Quality-aware keeper rescoring

**Files:**
- Create: `src/repgenr/stages/derep_keeper.py`
- Modify: `src/repgenr/stages/dereplicate.py` (`DereplicateParams.keeper`, call site, provenance)
- Modify: `src/repgenr/cli/cmd_bacterial.py`, `src/repgenr/cli/param_builders.py`, `src/repgenr/cli/cmd_run.py` (`--keeper`)
- Test: `tests/unit/test_derep_keeper.py`, `tests/integration/test_dereplicate_stage.py`

**Interfaces:**
- Produces: `quality_score(completeness: float, contamination: float) -> float` = `completeness - 5.0 * contamination` (dRep's default weighting, without the size terms).
- Produces: `rescore_representatives(result: DerepResult, quality: Mapping[str, tuple[float, float]], logger) -> tuple[DerepResult, int]` returning the rescored result and the number of swaps.
- Consumes: `Manifest.quality()` from Task 7; `check_result_complete` from Task 1.

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_derep_keeper.py`:

```python
"""Quality-aware keeper: the best-scoring cluster member becomes the representative."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.dereplicators.base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    DerepResult,
    check_result_complete,
)
from repgenr.stages.derep_keeper import quality_score, rescore_representatives

_LOG = logging.getLogger("keeper")


def _result() -> DerepResult:
    return DerepResult(
        representatives=[Path("/g/rep.fasta"), Path("/g/solo.fasta")],
        clusters={"rep.fasta": ["m1.fasta", "m2.fasta"], "solo.fasta": []},
        genome_status={
            "rep.fasta": STATUS_REPRESENTATIVE, "solo.fasta": STATUS_REPRESENTATIVE,
            "m1.fasta": STATUS_CONTAINED, "m2.fasta": STATUS_CONTAINED,
        },
    )


def test_score_penalises_contamination() -> None:
    assert quality_score(100.0, 0.0) > quality_score(100.0, 2.0)
    assert quality_score(95.0, 0.0) == 95.0


def test_better_member_replaces_representative() -> None:
    quality = {"rep.fasta": (90.0, 3.0), "m1.fasta": (99.0, 0.2), "m2.fasta": (95.0, 1.0)}
    out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 1
    assert sorted(p.name for p in out.representatives) == ["m1.fasta", "solo.fasta"]
    assert out.clusters["m1.fasta"] == ["m2.fasta", "rep.fasta"]
    assert out.genome_status["m1.fasta"] == STATUS_REPRESENTATIVE
    assert out.genome_status["rep.fasta"] == STATUS_CONTAINED
    assert out.representatives[0].parent == Path("/g")
    check_result_complete(out, ["rep.fasta", "solo.fasta", "m1.fasta", "m2.fasta"])


def test_unscored_member_never_wins() -> None:
    quality = {"rep.fasta": (90.0, 3.0)}
    out, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 0
    assert out.clusters == _result().clusters


def test_no_quality_keeps_adapter_choice() -> None:
    out, swaps = rescore_representatives(_result(), {}, _LOG)
    assert swaps == 0
    assert out == _result()


def test_tie_keeps_current_representative() -> None:
    quality = {"rep.fasta": (99.0, 0.0), "m1.fasta": (99.0, 0.0), "m2.fasta": (50.0, 0.0)}
    _, swaps = rescore_representatives(_result(), quality, _LOG)
    assert swaps == 0
```

Append to `tests/integration/test_dereplicate_stage.py` a test that seeds the manifest with quality for the three `genome_files` (rep worst), runs `run(ctx, DereplicateParams(tool=<fake that clusters all three under the first>, keeper="quality"))`, and asserts the representative in `derep/representatives/` is the best-quality file and `ctx.config.stages["dereplicate"].params["keeper"] == "quality"` and `params["keeper_swaps"] == 1`. Add a second test with `keeper="tool"` asserting no swap.

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/unit/test_derep_keeper.py -q`
Expected: ModuleNotFoundError.

- [ ] **Step 3: Implement `derep_keeper.py`**

```python
"""Quality-aware representative selection applied after any dereplicator.

Dereplicators pick representatives by connectivity, first-seen order or tool
defaults, which favours the most-sequenced genotype in a cluster. This step
re-picks each cluster's representative by assembly quality when the manifest
carries CheckM-style completeness and contamination, and leaves the adapter's
choice in place otherwise.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

from ..dereplicators.base import STATUS_CONTAINED, STATUS_REPRESENTATIVE, DerepResult

CONTAMINATION_WEIGHT = 5.0


def quality_score(completeness: float, contamination: float) -> float:
    """dRep-style score: completeness minus five times contamination."""
    return completeness - CONTAMINATION_WEIGHT * contamination


def rescore_representatives(
    result: DerepResult,
    quality: Mapping[str, tuple[float, float]],
    logger: logging.Logger,
) -> tuple[DerepResult, int]:
    """Return a copy of ``result`` whose representatives are the best-scoring
    cluster members, and the number of clusters whose representative changed.

    A member replaces the representative only when its score is strictly higher;
    genomes without quality never replace a scored representative, and a cluster
    with no scored genome keeps the adapter's choice.
    """
    if not quality:
        return result, 0

    rep_paths = {p.name: p for p in result.representatives}
    new_reps: list[Path] = []
    new_clusters: dict[str, list[str]] = {}
    status = dict(result.genome_status)
    swaps = 0

    def score(name: str) -> float | None:
        q = quality.get(name)
        return None if q is None else quality_score(*q)

    for rep_name, members in result.clusters.items():
        candidates = [rep_name, *members]
        current = score(rep_name)
        best_name, best = rep_name, current
        for m in members:
            s = score(m)
            if s is None:
                continue
            if best is None or s > best:
                best_name, best = m, s
        if best_name == rep_name:
            new_reps.append(rep_paths[rep_name])
            new_clusters[rep_name] = list(members)
            continue
        swaps += 1
        new_reps.append(rep_paths[rep_name].with_name(best_name))
        new_clusters[best_name] = sorted(n for n in candidates if n != best_name)
        status[best_name] = STATUS_REPRESENTATIVE
        status[rep_name] = STATUS_CONTAINED
        logger.info(
            "Keeper: %s replaces %s (score %.2f vs %s)",
            best_name, rep_name, best, "n/a" if current is None else f"{current:.2f}",
        )

    if swaps:
        logger.info("Quality-aware keeper changed %d of %d representatives", swaps, len(new_clusters))
    return DerepResult(
        representatives=sorted(new_reps),
        clusters=new_clusters,
        genome_status=status,
        genome_information=result.genome_information,
    ), swaps
```

The replacement path is `rep_paths[rep_name].with_name(best_name)`; `_write_contract` already falls back to `ctx.genomes_dir / rep.name` when that path does not exist.

- [ ] **Step 4: Wire into the stage and CLI**

`stages/dereplicate.py`: add `keeper: str = "quality"  # quality | tool` to `DereplicateParams`. In `run()`, after the adapter result and before `_reduce_by_taxonomy`:

```python
    keeper_swaps = 0
    if params.keeper == "quality":
        from .derep_keeper import rescore_representatives

        quality = _quality_lookup(ctx)
        if quality:
            result, keeper_swaps = rescore_representatives(result, quality, logger)
        else:
            logger.info("No assembly quality in the manifest; keeping the tool's representatives")
```

with `_quality_lookup(ctx)` mirroring `_taxon_lookup`'s try/except around `ctx.manifest.quality()`. Record `"keeper": params.keeper` and `"keeper_swaps": keeper_swaps` in `record_stage` params.

`cli/cmd_bacterial.py` `dereplicate` and `cli/cmd_run.py` `run`: add `keeper: str = typer.Option("quality", "--keeper", help="Representative choice per cluster: quality (CheckM score from GTDB) or tool (adapter's own).")`, validate with `_require_choice(keeper, {"quality", "tool"}, "--keeper")`, and pass through `dereplicate_params(...)` in `param_builders.py`.

- [ ] **Step 5: Gates, docs, commit**

Run: `pytest -q && ruff check src/ tests/ && mypy src/repgenr`

`docs/usage.md`: document `--keeper`. `docs/swot-derep.md` gap 2: mark closed with a pointer to `stages/derep_keeper.py`.

```bash
git add -A src/repgenr tests docs
git commit -m "feat(derep): quality-aware representative selection

After any dereplicator, each cluster's representative is re-picked by
completeness - 5 x contamination when the manifest carries GTDB CheckM
values. --keeper tool restores the adapter's choice."
```

---

### Task 9: Same keeper rule for taxonomy reduce and the Nextflow merge

**Files:**
- Modify: `src/repgenr/stages/dereplicate.py:401-448` (`_reduce_by_taxonomy`)
- Modify: `src/repgenr/stages/derep_steps.py` (`MergeParams.selection_tsv`), `src/repgenr/cli/cmd_steps.py` (`--selection-tsv` on `dereplicate-merge`)
- Modify: `nextflow/modules/local/derep_merge.nf`, `nextflow/subworkflows/local/dereplicate_scatter.nf`, `nextflow/subworkflows/local/bacterial_dataflow.nf`
- Test: `tests/integration/test_taxonomy_reduce.py`, `tests/integration/test_derep_steps.py`, `nextflow/tests/derep_merge.nf.test`

**Interfaces:**
- Consumes: `quality_score`, `rescore_representatives` from Task 8; `read_selection` with quality columns from Task 7.
- Produces: `MergeParams.selection_tsv: Path | None`; CLI `repgenr dereplicate-merge --selection-tsv PATH`; `DEREPLICATE_SCATTER` takes a second input `ch_selection`.

- [ ] **Step 1: Failing tests**

`tests/integration/test_taxonomy_reduce.py`: add a test where two representatives share a species, the one with the larger cluster has worse quality, and assert the better-quality one is kept under `keeper="quality"`; a second test with `keeper="tool"` asserts the largest-cluster rule still applies.

`tests/integration/test_derep_steps.py`: write a `selection.tsv` with quality columns for the genomes, run `dereplicate_merge(MergeParams(..., selection_tsv=path))` with the `_Halver` adapter and assert the merged representative is the best-quality member of its cluster.

`nextflow/tests/derep_merge.nf.test`: pass a selection file fixture to the process and assert the run succeeds (stub tag).

- [ ] **Step 2: Implement**

`_reduce_by_taxonomy`: the keeper line becomes

```python
        keeper = max(members, key=lambda n: (_score_or_min(n), len(result.clusters.get(n, [])), n))
```

with `_score_or_min` returning `quality_score(*quality[n])` when known, else `float("-inf")`, and `quality` fetched once via `_quality_lookup(ctx)` when `params.keeper == "quality"`, else `{}`. Pass `keeper` into `_reduce_by_taxonomy` as a parameter.

`derep_steps.py`: add `selection_tsv: Path | None = None` to `MergeParams`; after composing the merged result:

```python
    if params.selection_tsv is not None:
        quality = {
            r.filename: (r.completeness, r.contamination)
            for r in read_selection(params.selection_tsv)
            if r.completeness is not None and r.contamination is not None
        }
        composed, swaps = rescore_representatives(composed, quality, logger)
```

`cmd_steps.py` `dereplicate_merge`: add `selection_tsv: Path | None = typer.Option(None, "--selection-tsv", help="selection.tsv with quality columns; enables quality-aware representatives.")`.

Nextflow: `derep_merge.nf` gains `path selection` input (may be an empty list) and appends `--selection-tsv ${selection}` when present, following the optional-input idiom already used for `outgroup` in `tree2tax.nf`. `DEREPLICATE_SCATTER` takes `ch_selection` and forwards it; `bacterial_dataflow.nf` passes `ACQUIRE.out.selection`; `viral_dataflow.nf` passes `Channel.value([])`.

- [ ] **Step 3: Gates and commit**

Run: `pytest -q && ruff check src/ tests/ && mypy src/repgenr && nf-test test --tag stub && grep -rnP '[^\x00-\x7F]' nextflow`

```bash
git add -A src/repgenr tests nextflow
git commit -m "feat(derep): quality-aware keeper for taxonomy reduce and Nextflow merge"
```

- [ ] **Step 4: Open the Phase B pull request**

Push `feat/quality-keeper`; PR title "Quality-aware representative selection"; merge on green.

---

## Phase C: Reference-free SNP typer

### Task 10: ska2 adapter

**Files:**
- Create: `src/repgenr/snptypers/ska2.py`
- Modify: `pyproject.toml` (entry point `ska2 = "repgenr.snptypers.ska2:Ska2Typer"`)
- Modify: `CITATIONS.md`, `docs/usage.md`, `docs/swot-phylo.md`
- Test: `tests/unit/test_ska2_typer.py`

**Interfaces:**
- Produces: `Ska2Typer(SnpTyper)` with `requires_reference = False`, `accepted_extras = {"ksize", "min_freq"}`, `full_alignment = None`.
- Consumes: `run_tool`, `write_fofn` from `repgenr.core.process`, `ToolCapabilities`.

- [ ] **Step 1: Confirm the ska2 command-line surface**

Run `micromamba create -n ska2check -c bioconda -c conda-forge ska2 && micromamba run -n ska2check ska build --help && micromamba run -n ska2check ska align --help` (or check the ska2 README). The plan assumes: `ska build -f <tsv of name<TAB>path> -k 31 --threads N -o <prefix>` producing `<prefix>.skf`, and `ska align --min-freq 0.9 --filter no-ambig-or-const -o <out.aln> <prefix>.skf`. Adjust flag names in Step 3 if they differ, and record the version tested in the adapter docstring.

- [ ] **Step 2: Write the failing test**

```python
"""ska2 reference-free typer: build a split k-mer file, align variable sites."""

from __future__ import annotations

import logging
from pathlib import Path

from repgenr.snptypers import ska2 as mod
from repgenr.snptypers.base import SnpParams


def test_ska2_builds_then_aligns(tmp_path: Path, monkeypatch) -> None:
    calls: list[list[str]] = []

    def fake_run_tool(caps, argv, **kw):  # noqa: ANN001
        argv = [str(a) for a in argv]
        calls.append(argv)
        if argv[1] == "build":
            Path(argv[argv.index("-o") + 1] + ".skf").write_bytes(b"skf")
        elif argv[1] == "align":
            Path(argv[argv.index("-o") + 1]).write_text(">a\nA\n>b\nC\n")

    monkeypatch.setattr(mod, "run_tool", fake_run_tool)
    genomes = []
    for n in ("a", "b"):
        p = tmp_path / f"{n}.fasta"
        p.write_text(">x\nACGT\n")
        genomes.append(p)

    result = mod.Ska2Typer().call(
        genomes, None, tmp_path / "out", SnpParams(threads=2, extra={"ksize": 21}), logging.getLogger("t")
    )
    build, align = calls
    assert build[:2] == ["ska", "build"]
    assert build[build.index("-k") + 1] == "21"
    assert build[build.index("--threads") + 1] == "2"
    fofn = Path(build[build.index("-f") + 1])
    assert fofn.read_text().splitlines() == [f"a\t{genomes[0]}", f"b\t{genomes[1]}"]
    assert align[:2] == ["ska", "align"]
    assert result.core_snp_fasta.read_text().startswith(">a")
    assert result.full_alignment is None


def test_ska2_is_registered() -> None:
    from repgenr.snptypers.base import registry

    assert registry.get("ska2").requires_reference is False
```

- [ ] **Step 3: Implement**

```python
"""ska2 (split k-mer analysis) reference-free SNP typer.

Every genome is compared through split k-mers, so no genome is privileged as
the reference and reference-private errors do not bias distances. The output
is a variable-site alignment; there is no positional whole-genome alignment,
so recombination masking is not available with this typer.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities, parse_extra_int
from .base import SnpParams, SnpResult, SnpTyper


class Ska2Typer(SnpTyper):
    capabilities = ToolCapabilities(
        name="ska2",
        conda=("bioconda::ska2",),
        required_binaries=(BinarySpec("ska", version_args=("--version",)),),
        recommended_max_genomes=5000,
        default_params={"ksize": 31, "min_freq": 0.9},
        accepted_extras=frozenset({"ksize", "min_freq"}),
    )
    requires_reference = False

    def call(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: SnpParams,
        logger: logging.Logger,
    ) -> SnpResult:
        genomes = list(genomes)
        out_dir.mkdir(parents=True, exist_ok=True)
        ksize = parse_extra_int(params.extra, "ksize", self.capabilities.default_params["ksize"])
        min_freq = str(params.extra.get("min_freq", self.capabilities.default_params["min_freq"]))

        fofn = out_dir / "genomes.tsv"
        fofn.write_text(
            "".join(f"{g.stem}\t{g.resolve()}\n" for g in genomes), encoding="utf-8"
        )
        prefix = out_dir / "split_kmers"
        run_tool(
            self.capabilities,
            ["ska", "build", "-f", fofn, "-k", str(ksize), "--threads", str(params.threads), "-o", prefix],
            logger=logger, cwd=out_dir, log_prefix="ska-build",
            extra_mounts=sorted({g.parent for g in genomes}),
        )
        skf = Path(str(prefix) + ".skf")
        if not skf.exists():
            raise WorkdirError("ska build did not produce a .skf file")

        core = out_dir / "core_snp.fasta"
        run_tool(
            self.capabilities,
            ["ska", "align", "--min-freq", min_freq, "--filter", "no-ambig-or-const",
             "--threads", str(params.threads), "-o", core, skf],
            logger=logger, cwd=out_dir, log_prefix="ska-align",
        )
        if not core.exists() or core.stat().st_size == 0:
            raise WorkdirError("ska align produced no variable-site alignment")
        return SnpResult(core_snp_fasta=core, masked=False, full_alignment=None)
```

Check that `run_tool` accepts `extra_mounts` (it is used by `galah.py` and `cactus.py`); if the keyword differs, match it.

- [ ] **Step 4: Register and document**

`pyproject.toml` under `[project.entry-points."repgenr.snptypers"]`: `ska2 = "repgenr.snptypers.ska2:Ska2Typer"`. Reinstall with `pip install -e .` so the entry point is visible. `CITATIONS.md`: add Derelle et al. 2024, "Seamless, rapid, and accurate analyses of outbreak genomic data using split k-mer analysis", Genome Research. `docs/usage.md`: list `ska2` under snptypers with the sentence "reference-free; not compatible with `--mask`". `docs/swot-phylo.md`: mark the reference-bias gap as mitigated by an available reference-free typer.

- [ ] **Step 5: Gates and commit; open the Phase C pull request**

Run: `pip install -e . && pytest -q && ruff check src/ tests/ && mypy src/repgenr`

```bash
git add src/repgenr/snptypers/ska2.py pyproject.toml tests/unit/test_ska2_typer.py CITATIONS.md docs
git commit -m "feat(snptype): ska2 reference-free SNP typer"
```

Push `feat/ska2-snptyper`; PR title "ska2 reference-free SNP typer"; merge on green. If a real ska2 binary is available locally, run one manual `repgenr snptype --tool ska2` on the 27-leaf Francisella set and paste the site count into the PR.

---

## Phase D: Release and hygiene

### Task 11: Single source of truth for the version

**Files:**
- Create: `tests/unit/test_version_consistency.py`
- Modify: `src/repgenr/core/versions.py` (if the repgenr version there is a literal, derive it from `importlib.metadata.version("repgenr")`)

- [ ] **Step 1: Write the test**

```python
"""pyproject, the Nextflow manifest and the installed package agree on the version."""

from __future__ import annotations

import re
import tomllib
from importlib.metadata import version
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_versions_agree() -> None:
    py = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    nf = re.search(r"version\s*=\s*'([^']+)'", (ROOT / "nextflow" / "nextflow.config").read_text(encoding="utf-8"))
    assert nf is not None
    assert nf.group(1) == py
    assert version("repgenr") == py


def test_changelog_has_section_for_version() -> None:
    py = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    assert f"## [{py}]" in (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
```

- [ ] **Step 2: Run, make it pass, commit**

Run: `pytest tests/unit/test_version_consistency.py -q`. It passes today at 2.0.0 and will guard Task 13. If `core/versions.py` hard-codes `2.0.0`, replace it with `importlib.metadata.version("repgenr")`.

```bash
git add tests/unit/test_version_consistency.py src/repgenr/core/versions.py
git commit -m "test: pin pyproject, Nextflow manifest and changelog to one version"
```

---

### Task 12: Move binaries out of the tracked tree

**Files:**
- Remove from index: `wiki/mock_dataset.fasta.gz`, `wiki/images/MBDC_2025_Ohrman.pdf` (the tracked name contains a non-ASCII character; use `git ls-files wiki` to get it exactly)
- Move: `wiki/images/*.png` to `docs/images/`, `wiki/fasta_simulate_sequences.py` to `scripts/`
- Modify: `.gitignore`, `README.md`

- [ ] **Step 1: Confirm nothing references the files**

Run: `grep -rn "wiki/" README.md docs nextflow src tests benchmarks scripts` and `grep -rn "mock_dataset\|fasta_simulate" --include='*.md' --include='*.py' --include='*.nf' .`
Expected: no hits outside `wiki/` itself. If a hit exists, update the path in the referencing file as part of this task.

- [ ] **Step 2: Move and remove**

```bash
git mv wiki/images/example_francisellaceae.png wiki/images/example_francisellaceae_krakenA.png wiki/images/example_francisellaceae_krakenB.png wiki/images/example_francisellaceae_tul_expansion.png wiki/images/example_francisellaceae_tul_expansion_node_renamed.png wiki/images/examply_phylo_accurate.png docs/images/
git mv wiki/fasta_simulate_sequences.py scripts/fasta_simulate_sequences.py
git rm --cached wiki/mock_dataset.fasta.gz "wiki/images/$(git ls-files wiki/images | grep -i pdf | xargs basename)"
git rm -r --cached wiki
```

Keep the local `wiki/` directory on disk (it is already in `.gitignore` line 2). Add a line to README under "Development": "Example figures live in `docs/images/`; the synthetic-genome generator is `scripts/fasta_simulate_sequences.py`. The legacy 10 MB mock dataset is no longer tracked; regenerate a test set with `benchmarks/genomegen.py`."

- [ ] **Step 3: Verify and commit**

Run: `git ls-files wiki` (expect empty) and `pytest -q`.

```bash
git commit -am "chore: stop tracking wiki/ binaries; keep figures under docs/images"
```

History rewrite (to shrink the pack from 15.6 MB) is not part of this task: it needs a force-push and every clone to re-fetch. Offer it to the maintainer as a separate, explicitly approved step.

---

### Task 13: CI format gate and caching

**Files:**
- Modify: `.github/workflows/ci.yml`
- Create: `.git-blame-ignore-revs`
- Modify: every file `ruff format` touches under `src/`, `tests/`, `benchmarks/`

- [ ] **Step 1: Format once**

Run: `ruff format src/ tests/ benchmarks/ && ruff check src/ tests/ && pytest -q`
Expected: about 110 files rewritten, tests still green.

- [ ] **Step 2: Commit the mechanical change on its own**

```bash
git add -A src tests benchmarks
git commit -m "style: apply ruff format across the tree (mechanical)"
echo "$(git rev-parse HEAD)  # style: apply ruff format across the tree" > .git-blame-ignore-revs
git add .git-blame-ignore-revs
git commit -m "chore: ignore the format commit in git blame"
```

- [ ] **Step 3: Gate and cache**

In `.github/workflows/ci.yml`:

- Under `setup-python`, add `cache: pip`.
- After the Ruff step, add:

```yaml
      - name: Ruff format
        run: ruff format --check src/ tests/ benchmarks/
```

- Under `setup-micromamba` in `nf-test-e2e`, add `cache-environment: true`.
- In the `Install nf-test` steps, pin the installer by downloading the release tarball for 0.9.5 from GitHub instead of `wget | bash`:

```yaml
      - name: Install nf-test
        run: |
          curl -fsSL -o nf-test.tar.gz https://github.com/askimed/nf-test/releases/download/v0.9.5/nf-test-0.9.5.tar.gz
          tar -xzf nf-test.tar.gz nf-test
          sudo mv nf-test /usr/local/bin/
```

- [ ] **Step 4: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: gate ruff format, cache pip and micromamba, pin nf-test installer"
```

---

### Task 14: Changelog, version bump and tag

**Files:**
- Modify: `CHANGELOG.md`, `pyproject.toml:7`, `nextflow/nextflow.config:143`, `nextflow/nextflow.config:142` (`nextflowVersion`)

- [ ] **Step 1: Reconcile the changelog**

Rename `## [Unreleased]` to `## [3.0.0] - <today's date>` and open a new empty `## [Unreleased]` above it. Add entries for Phases A to D under Added / Changed / Fixed:

- Added: quality-aware representative selection (`--keeper`), ska2 typer, `snp/full_alignment.fasta`, `pipeline_info/software_versions.yml`, manifest schema v2 with quality columns.
- Changed: maskers receive the whole-genome alignment (breaking for third-party maskers: new `mask()` signature); `SnpParams.mask` removed; `xmfa_to_fasta` lost its `flank` argument; `selection.tsv` gains two optional columns.
- Fixed: dereplication refuses incomplete results; stage flags no longer alter fingerprints of tools that ignore them; doctor is read-only.

- [ ] **Step 2: Bump**

`pyproject.toml` version `3.0.0`; `nextflow.config` `version = '3.0.0'` and `nextflowVersion = '>=24.10.0'` (the oldest version nf-schema 2.3 supports; verify against the nf-schema changelog and use its stated floor).

- [ ] **Step 3: Gates**

Run: `pip install -e . && pytest -q && ruff check src/ tests/ && ruff format --check src/ tests/ benchmarks/ && mypy src/repgenr && nf-test test --tag stub`

- [ ] **Step 4: Commit and open the Phase D pull request**

```bash
git add CHANGELOG.md pyproject.toml nextflow/nextflow.config
git commit -m "release: 3.0.0"
```

Push `chore/release-3.0.0`; PR title "Release 3.0.0"; merge on green.

- [ ] **Step 5: Tag (requires maintainer confirmation)**

Tagging triggers `.github/workflows/release.yml`, which creates a public GitHub Release. Ask before running:

```bash
git checkout main && git pull
git tag -a v3.0.0 -m "RepGenR 3.0.0"
git push origin v3.0.0
```

Enabling PyPI Trusted Publishing (the commented job in `release.yml`) is a maintainer decision and is left out of this plan.

---

## Self-review

**Spec coverage.** Audit recommendations 1 (completeness assertion) Task 1; 2 (extras gating) Task 2; 3 (Gubbins input) Task 3; 4 (keeper scoring) Tasks 7 to 9; 5 (ska2) Task 10; 6 (release) Tasks 11, 14; 7 (hygiene: wiki, ruff format, caching, versions.yml, dead flank, read-only doctor) Tasks 4, 5, 6, 12, 13. Verified low-severity items not covered here: wrong error class for a bad `--container`, child not killed on read-loop error, non-atomic tree write in phylo, retry on exit 130, IQ-TREE bootstrap validation, Wave cache ignoring platform. These are one-line fixes each and can be folded into Phase A by the executor if time allows; otherwise they go to the backlog.

**Type consistency.** `check_result_complete(result, genome_names)` is used with the same signature in Tasks 1, 8. `MaskParams`, `SnpResult.full_alignment` are defined in Task 3 and consumed in Task 10. `Manifest.quality()` returns `dict[str, tuple[float, float]]` in Task 7 and is consumed as a `Mapping[str, tuple[float, float]]` in Tasks 8, 9. `gated_extra(registry, tool, key, value)` is used identically at all four CLI sites.

**Placeholders.** None; every step names the file, the code and the command.
