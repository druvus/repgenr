# Adding a tool

A new tool is a self-contained adapter plus one entry-point line. The core never
imports adapters directly, so third parties can ship adapters in their own
package without modifying RepGenR.

## 1. Subclass the family ABC

Each family defines an ABC and a normalized result type:

| Family | ABC | Method | Returns |
|--------|-----|--------|---------|
| Dereplicator | `repgenr.dereplicators.base.Dereplicator` | `dereplicate(genomes, out_dir, params, logger)` | `DerepResult` |
| Aligner | `repgenr.aligners.base.Aligner` | `align(genomes, reference, out_dir, params, logger)` | `AlignResult` |
| SnpTyper | `repgenr.snptypers.base.SnpTyper` | `call(genomes, reference, out_dir, params, logger)` | `SnpResult` |
| TreeBuilder | `repgenr.treebuilders.base.TreeBuilder` | `build(msa_or_genomes, out_dir, params, logger)` | Newick `Path` |

Declare `capabilities: ToolCapabilities` with the binaries the tool needs (so
preflight can check presence and version), default params, and scaling hints.
Use `repgenr.core.process.run` for every external call (never `shell=True`); use
`write_fofn` instead of a shell glob for large input lists.

Example dereplicator skeleton:

```python
from repgenr.core.binaries import BinarySpec
from repgenr.core.plugins import ToolCapabilities
from repgenr.dereplicators.base import Dereplicator, DerepResult

class MyDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="mytool",
        required_binaries=(BinarySpec("mytool", version_args=("--version",)),),
        supports_native_scaling=True,
    )

    def dereplicate(self, genomes, out_dir, params, logger) -> DerepResult:
        ...  # run the tool via repgenr.core.process.run, parse its output
        return DerepResult(representatives=[...], clusters={...}, genome_status={...})
```

The adapter only returns the dataclass; the stage writes the canonical contract
files. Do not write `derep/...`, `align/...` etc. from inside an adapter.

## 2. Register the entry point

In-tree, add a line to `pyproject.toml`:

```toml
[project.entry-points."repgenr.dereplicators"]
mytool = "repgenr.dereplicators.mytool:MyDereplicator"
```

Third parties add the same entry point in their own package's metadata. After
`pip install`, `repgenr list-tools` shows the new name and it is selectable via
`--tool mytool` (or `--aligner`, `--treebuilder`, etc.).

## 3. Test it

Unit-test by patching `repgenr.core.process.run` to assert the exact argument
vector and to drop canned tool output, then assert the returned dataclass.
Integration-test the owning stage with the adapter registered into the family
registry (`registry._classes["mytool"] = MyDereplicator`). See
`tests/integration/` for examples. Gate any test that needs the real binary with
`@pytest.mark.requires_binary`.
