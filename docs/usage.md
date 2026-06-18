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

The Nextflow layer runs these stages as typed data channels: each stage emits its
outputs (the metadata selection, genome FASTAs, per-chunk and merged
representatives, the tree, the taxonomy) as staged files that the next stage
consumes. There is no shared working directory; results are published under
`--outdir`. Nextflow owns the fan-out (scatter-gather dereplication).

## Quick start

```bash
nextflow run nextflow/main.nf \
    --mode bacterial \
    --metadata_args '-r 232.0 -v bac120 -d rep -l genus -tg francisella' \
    --outdir results \
    -profile standard
```

Run `nextflow run nextflow/main.nf --help` for the parameter summary.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Published results and execution reports. |
| `--mode` | `bacterial` | `bacterial` (GTDB) or `viral` (BV-BRC). |
| `--metadata_args` | see config | Arguments for the bacterial metadata stage. |
| `--vmetadata_args` / `--vgenome_args` | see config | Viral metadata / genome selection arguments. |
| `--derep_tool` | `skder` | Dereplicator for the scatter-gather step. |
| `--derep_process_size` | `null` | Genomes per dereplication chunk (single chunk if unset). |
| `--derep_primary_ani` / `--derep_secondary_ani` / `--derep_aligned_fraction` | `0.90` / `0.99` / `0.50` | ANI / aligned-fraction thresholds. |
| `--phylo_args` | `--treebuilder mashtree` | Aligner or tree builder for the phylogeny. |
| `--tree2tax_args` | `--include-dereplicated` | tree-to-taxonomy (FlexTaxD) arguments. |

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

For large dereplication inputs (10k+ genomes), set a chunk size so the
dereplication scatters across tasks (one per chunk):

```bash
nextflow run nextflow/main.nf --outdir results \
    --derep_tool sourmash --derep_process_size 2000 -profile slurm
```

Resource labels (`process_low/medium/high`) scale memory and time with the retry
attempt, so a task killed for memory or time is resubmitted with more headroom.
Tune the label values per environment in `nextflow.config`.

### Scatter-gather dereplication

For horizontal scaling, the `DEREPLICATE_SCATTER` subworkflow groups genomes into
chunks of `--derep_process_size`, dereplicates each chunk as a separate task (one
per node on HPC), and dereplicates the union of the chunk representatives once
more (the two-stage reduce-tree expressed as typed data channels). It wraps the
`repgenr dereplicate-chunk` / `dereplicate-merge` CLI steps and is driven by
`--derep_tool` and the `--derep_*_ani` / `--derep_aligned_fraction` thresholds.

A standalone harness runs it on a directory of genome FASTAs (no GTDB/NCBI
front-end), useful for testing and for dereplicating a local collection:

```bash
nextflow run nextflow/tests/dereplicate_scatter.nf -c nextflow/nextflow.config \
    --genomes_dir <DIR> --derep_tool sourmash --derep_process_size 2000 \
    --outdir results -profile standard
```

Add `-stub` to exercise the wiring without running the tools.

### Pipeline structure

`nextflow/main.nf` dispatches by `--mode` to one of two data-channel subworkflows
that share the dereplication, phylo and tree2tax modules:

```
bacterial: ACQUIRE  (metadata -> genome)  -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
viral:     VACQUIRE (vmetadata -> vgenome) -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
```

`metadata` emits a portable `selection.tsv`; `genome-fetch` downloads the genomes
and emits them as a channel feeding the scatter-gather dereplication; `phylo` and
`tree2tax` run in task-local working directories and emit `tree.nwk`,
`tree2tax.tsv` and `genomes_map.tsv` to `--outdir`. Add `-stub` to any run for a
quick wiring check without external tools.
