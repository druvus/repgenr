# repgenr: Output

The CLI writes all stage outputs into the shared working directory
(`--workdir`). The Nextflow data-channel pipeline instead flows results between
processes as staged channel files and publishes under `--outdir` (default
`results/`): `metadata/selection.tsv`, the dereplication contract under
`dereplicate/`, the phylogeny under `phylo/` (the tree and the tree builder's
files under `phylo/tree/`, the alignment it was built from under
`phylo/align/` or `phylo/snp/`), `tree2tax.tsv` and `genomes_map.tsv` at the
top level, and the execution reports under `pipeline_info/`.

## Working directory layout

| Path | Produced by | Description |
|------|-------------|-------------|
| `manifest.sqlite` | metadata, ingest, vgenome | Genome manifest (accessions, taxonomy, CheckM quality when known, dereplication status). |
| `selection.tsv` | metadata, ingest, vgenome | The selected genomes: accession, taxonomy, outgroup flag, canonical filename and quality when known. The hand-off every later stage reads, and a resume input. |
| `outgroup_accession.txt` | metadata, ingest, vgenome | Accession of the outgroup genome, read by phylo and tree2tax for rooting. Absent when there is no outgroup. |
| `repgenr.yaml` | all stages | Provenance: tool name, parameters, resolved tool versions, completion timestamps. |
| `repgenr.log` | all stages | Run log. |
| `genomes/` | genome, ingest, vgenome | Genome FASTAs, one per selected accession. |
| `outgroup/` | genome, ingest, vgenome | Outgroup genome for rooting. |
| `missing_accessions.txt` | genome | Accessions the download did not return; the completeness guard of later stages reads it. |
| `reads.tsv` | reads | The selected sequencing runs: run, sample and study accessions, organism and taxid, the resolved family/genus/species tokens, platform, instrument, layout, bases, and the FASTQ locations, checksums and sizes ENA reports (empty when ENA holds no FASTQ mirror). |
| `virus_download_wd/` | vmetadata | Downloaded viral sequences and the metadata tables `vgenome` selects from. `virus_metadata_base.tsv` (and `virus_metadata_ncbi.tsv` on the BV-BRC path) at the workdir root are copies of those tables. |
| `derep/` | dereplicate | Representative genomes and per-tool intermediates. |
| `derep/clusters.tsv` | dereplicate | `representative<TAB>member`, one row per genome; a representative also lists itself. |
| `derep/genome_status.tsv` | dereplicate | Per-genome status: `representative`, `contained` or `fail_qc`. |
| `derep/cluster_summary.tsv` | dereplicate, cluster-summary | One row per representative: member count, species spanned and keeper quality against the members (below). |
| `snp/core_snp.fasta` | snptype | Core-SNP (variable-site) alignment; masked in place when `--mask` is set. |
| `snp/full_alignment.fasta` | snptype | Whole-genome alignment in reference coordinates, when the SNP typer produces one (snippy, parsnp, simple); required input for `--mask`. |
| `snp/snp_distance_matrix.tsv` | snptype | Pairwise SNP distances between genomes; the `simple` typer writes it, the others do not. |
| `scratch/` | snptype | The typer's per-genome intermediates. Each genome's are removed once its consensus has been read; a genome whose chain failed keeps its own. |
| `align/msa.fasta` | phylo | Whole-genome alignment from the aligner (`--msa-source aligner`). |
| `align/msa_source.json`, `snp/msa_source.json` | phylo | Stamp beside the alignment phylo built: the genome set, the source settings and the alignment's digest, so a later `phylo` that changes only the tree builder reuses it. |
| `tree/` | phylo | Phylogeny (`tree.nwk`) and the tree builder's own files (logs, bootstrap trees). |
| `segments.tsv` | vgenome | With `--group-segments`: each grouped isolate's token and its member segment accessions, one row per member. Absent otherwise. |
| `genomes_map.tsv` | tree2tax | Accession to leaf: each representative, its dereplicated members, and for a grouped viral isolate its member segment accessions. |
| `tree2tax.tsv` | tree2tax | FlexTaxD-compatible taxonomy derived from the tree. |

## Cluster summary

`derep/cluster_summary.tsv` condenses `clusters.tsv` into one row per
representative, largest cluster first. `dereplicate` writes it after every run
and `repgenr cluster-summary -wd <workdir>` regenerates it for an existing
working directory from `clusters.tsv` and the manifest, without rerunning the
dereplicator.

| Column | Meaning |
|--------|---------|
| `representative` | Keeper filename, as in `clusters.tsv`. |
| `n_members` | Genomes contained under the keeper (the keeper itself is not counted). |
| `n_species` | Distinct species across keeper and members, parsed from the canonical filenames. |
| `species` | Those species, comma-separated, the keeper's first. |
| `rep_completeness`, `rep_contamination` | CheckM values of the keeper from the manifest; blank when unknown. |
| `member_max_completeness`, `member_min_contamination` | Best values among the scored members; blank when no member is scored. |
| `best_member` | Highest-scoring genome in the cluster by completeness minus five times contamination, keeper included. Equals `representative` when the keeper is already the best; blank when nothing in the cluster is scored. |

A row whose `best_member` differs from its `representative` marks a cluster
where a member outscores the keeper. This is expected under `--keeper tool`,
and can also follow `--reduce`, which merges representatives by taxon. The
quality columns come from the manifest in the workdir CLI and from
`--selection-tsv` in the Nextflow steps; without either they stay blank and
the size and species columns still apply.

## Pipeline information

Under `<outdir>/pipeline_info/`, each run writes timestamped Nextflow execution
reports:

- `execution_report_*.html` -- resource usage and per-task summary.
- `execution_timeline_*.html` -- task timeline.
- `execution_trace_*.txt` -- machine-readable trace of every task.
- `pipeline_dag_*.html` -- the workflow DAG.
- `software_versions.yml` -- every process's resolved tool versions, collected
  and de-duplicated from each process's `versions.yml` fragment.

These are useful for diagnosing resource limits (the retry strategy scales
memory and time per attempt) and for provenance.
