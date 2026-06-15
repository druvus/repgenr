# RepGenR architecture

## Overview

RepGenR is an importable Python package (`src/repgenr`, Python 3.12+). A thin
Typer CLI (`repgenr.cli.main`) dispatches each subcommand to an importable stage
function `repgenr.stages.<name>.run(ctx, params)`. Stages share a working
directory; a `WorkdirContext` resolves the canonical path layout, logging, the
`repgenr.yaml` config/provenance, and a SQLite genome manifest.

```
cli/main.py            Typer app; one command per stage
core/
  context.py           WorkdirContext: paths + services
  config.py            repgenr.yaml (per-stage tool, params, versions)
  manifest.py          SQLite genome inventory
  contracts.py         canonical inter-stage file readers/writers
  process.py           subprocess wrapper (no shell, no globs)
  binaries.py          binary presence + version preflight
  plugins.py           ToolCapabilities + entry-point Registry
stages/                metadata, genome, dereplicate, snptype, phylo, tree2tax, ...
dereplicators/         Dereplicator ABC + adapters
aligners/              Aligner ABC + adapters
snptypers/             SnpTyper ABC + adapters
treebuilders/          TreeBuilder ABC + adapters
converters/            XMFA/MAF/HAL/GFA -> MSA-FASTA
tree/                  Newick + neighbor-joining
```

## Pluggable tool families

Four families, each an ABC plus a `Registry` bound to an entry-point group:
`repgenr.dereplicators`, `repgenr.aligners`, `repgenr.snptypers`,
`repgenr.treebuilders`. In-tree and third-party adapters are discovered the same
way, so the core never imports a concrete adapter.

Adapters return normalized dataclasses (`DerepResult`, `AlignResult`,
`SnpResult`, or a Newick path); the owning stage writes the canonical contract
files. This keeps tools within a family interchangeable with no downstream
change.

## Orthogonal phylogenetics axes

The phylo stage composes three independent choices:

* genome set: dereplicated representatives or all genomes
* MSA source: a whole-genome aligner OR a SNP typer's core-SNP alignment
  (skipped for alignment-free tree builders)
* tree builder: MSA-based (iqtree/fasttree/raxmlng) or alignment-free
  (mashtree/sourmash)

Outgroup rooting is handled once in the stage, regardless of the tools chosen.

## Scaling

* Genome state lives in a SQLite manifest, not `str(dict)`/`pickle`; no repeated
  directory scans.
* `process.run` always passes argument vectors and uses file-of-filenames
  instead of shell globs, so large genome sets do not hit `ARG_MAX`.
* Dereplicators declare `supports_native_scaling`; tools that do not (dRep) are
  wrapped with two-stage chunking, while skDER/galah/sourmash run directly.
* Nextflow provides the actual parallelism and HPC/cloud execution; resource
  labels assign heavy stages (Cactus) to large nodes.

## Data contracts

Canonical files (owned by `core/contracts.py`):

```
derep/representatives/    representative genome FASTAs
derep/clusters.tsv        representative <TAB> member
derep/genome_status.tsv   genome <TAB> status
align/msa.fasta           aligner output (MSA-FASTA)
snp/core_snp.fasta        SNP typer output (core-SNP alignment)
tree/tree.nwk             Newick tree
tree2tax.tsv              child <TAB> parent (FlexTaxD)
genomes_map.tsv           accession <TAB> leaf
```
