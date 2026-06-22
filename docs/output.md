# repgenr: Output

The CLI writes all stage outputs into the shared working directory
(`--workdir`). The Nextflow data-channel pipeline instead flows results between
processes as staged channel files and publishes the key deliverables (the tree,
`tree2tax.tsv`, `genomes_map.tsv`) plus execution reports under `--outdir`
(default `results/`).

## Working directory layout

| Path | Produced by | Description |
|------|-------------|-------------|
| `manifest.sqlite` | metadata | Genome manifest (accessions, taxonomy, dereplication status). |
| `repgenr.yaml` | all stages | Provenance: tool name, parameters, resolved tool versions, completion timestamps. |
| `repgenr.log` | all stages | Run log. |
| `genomes/` | genome | Downloaded genome FASTAs, one per selected accession. |
| `outgroup/` | genome | Outgroup genome for rooting. |
| `derep/` | dereplicate | Representative genomes and per-tool intermediates. |
| `tree/` | phylo | Phylogeny (`tree.nwk`) and aligner/tree-builder intermediates. |
| `genomes_map.tsv` | tree2tax | Map from each representative to its dereplicated members. |
| `tree2tax.tsv` | tree2tax | FlexTaxD-compatible taxonomy derived from the tree. |

## Pipeline information

Under `<outdir>/pipeline_info/`, each run writes timestamped Nextflow execution
reports:

- `execution_report_*.html` -- resource usage and per-task summary.
- `execution_timeline_*.html` -- task timeline.
- `execution_trace_*.txt` -- machine-readable trace of every task.
- `pipeline_dag_*.html` -- the workflow DAG.

These are useful for diagnosing resource limits (the retry strategy scales
memory and time per attempt) and for provenance.
