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
the core (see `docs/adding-tools.md`).

| Family | Built-in adapters |
|--------|-------------------|
| Dereplicators | `drep`, `skder`, `galah`, `sourmash` |
| Aligners | `progressivemauve`, `sibeliaz`, `cactus` |
| SNP typers | `simple` (samtools/bcftools), `snippy`, `parsnp` (optional Gubbins masking) |
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

## Usage

```bash
WD=./francisella

# Full GTDB table (release-pinned):
repgenr metadata -wd $WD -r 232.0 -v bac120 -d rep -l genus -tg francisella
# Or query just the target taxon via the GTDB API (no full-table download):
repgenr metadata -wd $WD --source api -d rep -l genus -tg francisella

repgenr genome -wd $WD
repgenr dereplicate -wd $WD --tool skder -t 16
repgenr phylo -wd $WD --aligner progressivemauve --treebuilder iqtree
repgenr tree2tax -wd $WD --include-dereplicated
```

Alternatives:

```bash
# Scalable dereplication then an alignment-free tree
repgenr dereplicate -wd $WD --tool skder
repgenr phylo -wd $WD --treebuilder mashtree

# SNP-based phylogeny (core-SNP alignment as the MSA source)
repgenr snptype -wd $WD --tool simple
repgenr phylo -wd $WD --msa-source snptype --treebuilder iqtree --mask gubbins
```

## Nextflow

```bash
nextflow run nextflow/main.nf -profile standard --workdir $WD \
    --metadata_args "-r 232.0 -v bac120 -d rep -l genus -tg francisella"
```

Profiles: `standard` (local), `slurm`, `cloud`, `test`. Resource labels
(`process_low/medium/high`) are tuned per profile; heavy aligners such as Cactus
use `process_high`.

## Development

```bash
pip install -e ".[dev]"
ruff check src/ tests/
mypy src/repgenr
pytest -q
```

See `docs/architecture.md` for the design and `docs/adding-tools.md` for writing
a new adapter.

## License

MIT.
