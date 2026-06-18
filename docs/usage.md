# repgenr: Usage

## Introduction

RepGenR builds a representative-genome repository for a taxonomic group and the
matching FlexTaxD taxonomy. The pipeline runs five stages in order:

1. **metadata** -- select genome accessions (GTDB for bacteria/archaea, BV-BRC
   for viruses).
2. **genome** -- download and organise the selected genomes with NCBI Datasets.
3. **dereplicate** -- cluster genomes by ANI and pick representatives.
4. **phylo** -- build a phylogeny from the representatives.
5. **tree2tax** -- emit a FlexTaxD-compatible taxonomy from the tree.

The current Nextflow layer orders these stages and assigns resource labels and
profiles; the data flows through a single shared working directory
(`--workdir`). A later Phase 4 increment replaces this with typed data channels.

## Quick start

```bash
nextflow run nextflow/main.nf \
    --workdir /path/to/workdir \
    --mode bacterial \
    -profile standard
```

Run `nextflow run nextflow/main.nf --help` for the parameter summary.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--workdir` | (required) | Shared RepGenR working directory. |
| `--outdir` | `results` | Location for execution reports (and, later, published results). |
| `--mode` | `bacterial` | `bacterial` (GTDB) or `viral` (BV-BRC). |
| `--metadata_args` | see config | Arguments for the bacterial metadata stage. |
| `--genome_args` | `''` | Arguments for the genome download stage. |
| `--vmetadata_args` / `--vgenome_args` | see config | Viral metadata / genome selection arguments. |
| `--dereplicate_args` | `--tool skder` | Dereplication tool and ANI thresholds. |
| `--phylo_args` | `--treebuilder mashtree` | Aligner or tree builder for the phylogeny. |
| `--tree2tax_args` | `--include-dereplicated` | tree-to-taxonomy (FlexTaxD) arguments. |
| `--derep_process_size` | `null` | Genomes per chunk for two-stage dereplication. |
| `--derep_num_processes` | `null` | Parallel chunk workers for two-stage dereplication. |

Parameters are validated against `nextflow/nextflow_schema.json` at launch.

## Profiles

Combine an executor profile with an optional container profile, e.g.
`-profile slurm,singularity`.

- **Executors**: `standard` (local), `slurm`, `cloud` (AWS Batch).
- **Containers**: `docker`, `singularity`, `wave`. These set RepGenR's own
  adapter-level container backend (`--container ...`), which runs each external
  tool in a pinned image. RepGenR itself must be available to the Nextflow
  process.
- **`test`**: minimal resources and a small target for a quick smoke run.

## Scaling

For large dereplication inputs (10k+ genomes), enable two-stage chunking:

```bash
nextflow run nextflow/main.nf --workdir <DIR> \
    --derep_process_size 2000 --derep_num_processes 4
```

Resource labels (`process_low/medium/high`) scale memory and time with the retry
attempt, so a task killed for memory or time is resubmitted with more headroom.
Tune the label values per environment in `nextflow.config`.
