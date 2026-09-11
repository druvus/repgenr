# RepGenR

RepGenR (Representative-Genome Repositories) builds dereplicated, representative
genome repositories for large-scale genomic studies. It selects taxa, downloads
genomes, clusters them by average nucleotide identity (ANI), computes
phylogenetic trees, and emits taxonomy files for downstream tools such as
FlexTaxD.

Version 2 is a modular, importable Python package (Python 3.12+) with four
pluggable tool families and an optional Nextflow pipeline for scaling to
thousands of genomes.

## Pipeline

```
metadata -> genome -> dereplicate -> (align | snptype) -> phylo -> tree2tax
```

Each stage is a `repgenr` subcommand and an importable function. Stages
communicate through a single working directory whose state is recorded in
`repgenr.yaml` (provenance) and a SQLite genome manifest.

## Pluggable tools

Tools are discovered through Python entry points; adding one needs no change to
the core (see [docs/developing.md](docs/developing.md)).

| Family | Built-in adapters |
|--------|-------------------|
| Dereplicators | `drep`, `skder`, `galah`, `sourmash` |
| Aligners | `progressivemauve`, `sibeliaz`, `cactus` |
| SNP typers | `simple` (samtools/bcftools), `snippy`, `parsnp`, `ska2` (reference-free) |
| Maskers | `gubbins` (recombination masking via `--mask`) |
| Tree builders | `iqtree`, `fasttree`, `raxmlng` (MSA), `mashtree`, `sourmash` (alignment-free) |

`repgenr list-tools` prints what is available in your environment.

## Installation

```bash
# Tools + Python environment (conda/mamba)
mamba env create -f environment.yml
mamba activate repgenr

# Or just the Python package (tools must be on PATH separately)
pip install .
```

Cactus is distributed separately (containers/binaries); see its documentation.

## Quick start

```bash
WD=./francisella
repgenr run -wd $WD -d rep -l genus -tg francisella --tool skder --treebuilder iqtree
repgenr status -wd $WD     # which stages are done, and what to run next
```

`run` chains `metadata -> genome -> dereplicate -> phylo -> tree2tax`; each stage
is also its own subcommand, and `repgenr ingest` starts from genomes already on
disk. `--viral` selects the NCBI Virus path. Re-running a stage is a no-op
unless its parameters or inputs changed (`--force` overrides).

The same pipeline runs as typed Nextflow data channels (no shared working
directory) for HPC and cloud:

```bash
nextflow run nextflow/main.nf -profile standard --outdir results \
    --metadata_args "-r 232.0 --gtdb-version bac120 -d rep -l genus -tg francisella" \
    --derep_tool sourmash --phylo_args "--treebuilder mashtree"
```

Nextflow 26.04 or later is required.

## Documentation

| Page | What it covers |
|------|----------------|
| [docs/usage.md](docs/usage.md) | Running the pipeline: CLI stages, local genomes, viruses, resume, representative selection, SNP typing and masking, Nextflow parameters and profiles, containers, troubleshooting. |
| [docs/cli-reference.md](docs/cli-reference.md) | Every command and option, generated from the command tree. |
| [docs/output.md](docs/output.md) | The files each stage writes. |
| [docs/developing.md](docs/developing.md) | Architecture, data contracts, and how to add a tool adapter. |
| [docs/verification.md](docs/verification.md) | Which adapters have been run against their real tools, and measured runs. |
| [docs/audit/](docs/audit/README.md) | Records: the CLI matrix and the scaling and bias audit. |

## License

MIT.
