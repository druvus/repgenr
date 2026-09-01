# Adding a tool

A new tool is a self-contained adapter plus one entry-point line. The core never
imports adapters directly, so third parties can ship adapters in their own
package without modifying RepGenR: after `pip install`, the tool appears in
`repgenr list-tools`, in the registry-generated `--help` texts, and is
selectable by name (including from the Nextflow pipeline, whose schema does not
restrict tool names).

## 1. Subclass the family ABC

Each family defines an ABC and a normalized result type:

| Family | ABC | Method | Returns |
|--------|-----|--------|---------|
| Dereplicator | `repgenr.dereplicators.base.Dereplicator` | `dereplicate(genomes, out_dir, params, logger)` | `DerepResult` |
| Aligner | `repgenr.aligners.base.Aligner` | `align(genomes, reference, out_dir, params, logger)` | `AlignResult` |
| SnpTyper | `repgenr.snptypers.base.SnpTyper` | `call(genomes, reference, out_dir, params, logger)` | `SnpResult` |
| TreeBuilder | `repgenr.treebuilders.base.TreeBuilder` | `build(msa_or_genomes, out_dir, params, logger)` | Newick `Path` |
| Masker | `repgenr.maskers.base.Masker` | `mask(full_alignment, out_dir, params: MaskParams, logger)` | masked variable-site FASTA `Path` |

Per-family class attributes the stages branch on -- set them or inherit the
default deliberately:

* `TreeBuilder.input_kind`: `InputKind.GENOMES` for alignment-free builders
  (they receive genome files), else the default `InputKind.MSA_FASTA` (they
  receive one MSA path). Narrow the argument with `as_genome_list` /
  `as_msa_path` from `treebuilders.base`.
* `SnpTyper.requires_reference`: when True (default) the stage resolves a
  reference genome before calling you.
* A `Masker` reads the **whole-genome** alignment, not the SNP alignment, and
  returns the masked variable-site FASTA that replaces the typer's core-SNP
  output. A `SnpTyper` should therefore set `SnpResult.full_alignment`
  whenever it can produce one; `--mask` is refused with a clear error when the
  chosen typer leaves it None.
* Optional capability hooks: `Dereplicator.compare(...)` (all-vs-all
  similarity for `repgenr glance`; return a `CompareResult` whose CSV has
  `genome1`/`genome2`/`similarity` columns) and
  `TreeBuilder.distance_matrix(...)` (pairwise matrix for the viral outgroup
  step). Tools without an override are cleanly rejected when a user selects
  them for those features.

### ToolCapabilities

```python
capabilities = ToolCapabilities(
    name="mytool",
    required_binaries=(BinarySpec("mytool", version_args=("--version",)),),
    default_params={"mode": "fast"},          # merged under params.extra by the stage
    recommended_max_genomes=2000,             # drives auto-selection and scale warnings
    container="quay.io/biocontainers/mytool:1.2.3--0",  # or:
    conda=("bioconda::mytool",),              # Wave-resolved image when no explicit one
    accepted_extras=frozenset({"mode"}),      # extra keys the adapter actually reads
)
```

* `required_binaries` powers preflight (presence + minimum version) when
  running natively.
* **Set `container` and/or `conda`** or the tool cannot run under
  `--container`/`--wave` (the run falls back to the host with a warning).
* `recommended_max_genomes` feeds `--tool auto` (tightest fitting limit wins)
  and the over-scale warnings; unbounded (`None`) tools are picked when
  nothing bounded fits.
* `accepted_extras` names the `params.extra` keys the adapter reads. Users set
  them with `--tool-arg key=value` (dereplicate/snptype and the data-channel
  steps) or `--aligner-arg key=value` (phylo). Every stage warns, by name,
  about extra keys no adapter it runs declares in `accepted_extras` -- phylo
  drives two adapters (the MSA source and the tree builder) from one dict, so
  it warns only about keys neither of them reads -- and the CLI does not
  inject stage-level flags (such as `--virus`) into the extras of a tool that
  does not read them.
  Parse integer extras with `repgenr.core.plugins.parse_extra_int` so a bad
  value raises a clean `UserInputError`.

### Running the external tool

**Use `run_tool(self.capabilities, cmd, logger=logger, ...)` from
`repgenr.core.containers` for every external call** -- never `process.run`
directly and never `shell=True`. `run_tool` routes the command through the
active container backend (docker/singularity/Wave) using your `container` /
`conda` declaration and falls back to the host natively. Use `write_fofn`
instead of a shell glob for large input lists (declare the containing
directories via `extra_mounts=` since fofn contents are invisible to the
container backend), and pass result-on-stdout via `stdout_path=` (published
atomically only on success). `preflight()`'s return value (host versions or
the container image reference) is recorded into `repgenr.yaml` provenance and
surfaces in `repgenr versions` / the Nextflow `versions.yml` bridge.

Example dereplicator skeleton:

```python
from repgenr.core.binaries import BinarySpec
from repgenr.core.containers import run_tool
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import Dereplicator, DerepResult

class MyDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="mytool",
        required_binaries=(BinarySpec("mytool", version_args=("--version",)),),
        conda=("bioconda::mytool",),
    )

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        run_tool(self.capabilities, ["mytool", ...], logger=logger, log_prefix="mytool")
        ...  # parse the tool's output
        return DerepResult(representatives=[...], clusters={...}, genome_status={...})
```

The adapter only returns the dataclass; the stage writes the canonical contract
files. Do not write `derep/...`, `align/...` etc. from inside an adapter.

### Resume semantics

A stage's skip decision digests its declared inputs (`STAGE_INPUTS` in
`repgenr/cli/base.py`) plus its parameters and the container identity. An
adapter therefore needs no resume logic of its own -- but it must only read
the inputs the stage declares; a tool that reads undeclared files would
resume incorrectly. Defaults for stage parameters live on the params
dataclasses and are built through `repgenr/cli/param_builders.py`, which both
the manual commands and `repgenr run` share -- add new options there, not in
one entry point only.

## 2. Register the entry point

In-tree, add a line to `pyproject.toml`:

```toml
[project.entry-points."repgenr.dereplicators"]
mytool = "repgenr.dereplicators.mytool:MyDereplicator"
```

(the groups are `repgenr.dereplicators`, `.aligners`, `.snptypers`,
`.treebuilders`, `.maskers`). Third parties add the same entry point in their
own package's metadata. Programmatic registration (tests, embedders) goes
through the public API: `registry.register("mytool", MyDereplicator)` /
`registry.unregister("mytool")`, or the `register_tool` pytest fixture.

## 3. Test it

Unit-test by patching your adapter module's `run_tool` to assert the exact
argument vector and to drop canned tool output, then assert the returned
dataclass. For an in-tree adapter, add the tool to the parametrized contract
suite (`tests/unit/test_adapter_contracts.py` /
`test_snptyper_aligner_contracts.py`): one entry in the family's token dict
plus a canned-output branch in the fake runner -- a completeness assertion
fails if a built-in adapter is missing from the suite. Integration-test the
owning stage with the `register_tool` fixture. Gate any test that needs the
real binary with `@pytest.mark.requires_binary`.

## In-tree checklist for a new built-in tool

1. Adapter module + `capabilities` (with `container`/`conda`).
2. Entry point in `pyproject.toml`.
3. Contract-suite token entry + canned-output branch.
4. A row in the README adapter table and, if user-facing behavior differs,
   a note in `docs/usage.md`.

Everything else (CLI help, validation, `list-tools`, auto-selection,
containers, provenance, resume, Nextflow) picks the tool up automatically.
