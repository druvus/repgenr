# repgenr: Output

The CLI writes all stage outputs into the shared working directory
(`--workdir`). The Nextflow data-channel pipeline instead flows results between
processes as staged channel files and publishes under `--outdir` (default
`results/`): `metadata/selection.tsv` (or, in reads mode, `reads/reads.tsv`,
`reads/selection.tsv`, `reads/assembly_stats.tsv`, `reads/excused_runs.tsv` and
the `quality.tsv` and `classification.tsv` tables of the QC step), the
dereplication contract under `dereplicate/`, the phylogeny under `phylo/` (the tree and the tree builder's
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
| `ncbi_acc_download_list.txt` | genome | Accessions the download step will fetch, one per line (those without a FASTA in `genomes/`). With `--accession-list-only` it is the only output. |
| `missing_accessions.txt` | genome | Accessions the download did not return; the completeness guard of later stages reads it. |
| `reads.tsv` | reads | The selected sequencing runs: run, sample and study accessions, organism and taxid, the resolved family/genus/species tokens, platform, instrument, layout, bases, the FASTQ locations, checksums and sizes ENA reports (empty when ENA holds no FASTQ mirror), and ENA's library selection (RANDOM, MDA, PCR, ...). |
| `assemblies/<run>/` | assemble | Each assembled run's filtered contigs and its `assembly.ok` marker (assembler, version, metrics); a re-run skips runs that have one. |
| `assembly_stats.tsv` | assemble | Per-assembly metrics: assembler, contigs, total length, N50, largest contig, estimated coverage, the NCBI taxonomy used for the name, quality and classification columns once those steps run, and the polisher that corrected a long-read assembly. |
| `excused_runs.tsv` | assemble | Runs that produced no genome, with the step that gave up (`fetch`, `assemble`, ...) and the reason; the completeness guard excuses them like `missing_accessions.txt`. |
| `virus_download_wd/` | vmetadata | Downloaded viral sequences and the metadata tables `vgenome` selects from. `virus_metadata_base.tsv` (and `virus_metadata_ncbi.tsv` on the BV-BRC path) at the workdir root are copies of those tables. |
| `derep/` | dereplicate | Representative genomes and per-tool intermediates. |
| `derep/representatives/` | dereplicate | The representative genomes, one genome file per cluster; the distinct values of the first column of `clusters.tsv`. |
| `derep/unpacked/<representative>/` | derep-unpack | One directory per cluster with its member genomes (and the representative unless `--no-representant`), hard-linked from `genomes/` where the file system allows it and copied otherwise. The directory is named by the representative's file stem, or by its full file name when two representatives share a stem. Replaced on each run. |
| `derep/stock/<name>/` | derep-stock | A named, stored dereplication run written by `pack`: `clusters.tsv`, `genome_status.tsv`, `cluster_summary.tsv` and a `representatives/` directory of links to the representative genomes. `unpack` restores it. |
| `glance_clustering_dendrogram.pdf` | glance | dRep's clustering dendrogram over all genomes. |
| `glance_MASH_ANI_similarity_boxplot.png`, `glance_MASH_ANI_similarity_histogram.png` | glance | Box plot and histogram of the all-against-all Mash ANI values within `--plot-min`/`--plot-max`, one value per genome pair (self-comparisons left out); the histogram counts genome pairs. A run removes the previous glance plots and dendrogram once the comparison succeeds, so a plot with no values in range is absent rather than stale. |
| `glance_wd/` | glance | dRep working files; kept only with `--keep-files`. |
| `derep/clusters.tsv` | dereplicate | `representative<TAB>member`, one row per genome; a representative also lists itself. |
| `derep/genome_status.tsv` | dereplicate | Per-genome status: `representative`, `contained` or `fail_qc`. |
| `derep/cluster_summary.tsv` | dereplicate, cluster-summary | One row per representative: member count, species spanned and keeper quality against the members (below). |
| `snp/core_snp.fasta` | snptype | Core-SNP (variable-site) alignment; masked in place when `--mask` is set. |
| `snp/full_alignment.fasta` | snptype | Whole-genome alignment in reference coordinates, when the SNP typer produces one (snippy, parsnp, simple); required input for `--mask`. |
| `snp/snp_distance_matrix.tsv` | snptype | Pairwise SNP distances between genomes; the `simple` typer writes it, the others do not. |
| `scratch/` | genome, assemble, dereplicate, snptype, phylo | Working files, one subdirectory per stage. The files under `genome_download/` and `assemble/` are removed unless `--keep-files`. `dereplicate/` holds the dereplicator's working files (sketches, pairwise tables); it is kept after the run and cleared when `dereplicate` runs again. Under `snptype/`, each genome's intermediates are removed once its consensus has been read; a genome whose chain failed keeps its own. The `ska2` typer keeps its split k-mer file (`split_kmers.skf`) and sample list there; they are left in place and cleared when `snptype` runs again. |
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

The species columns come from the genome filenames
(`Family_genus_species_ACCESSION.fasta`), not from `selection.tsv`. A name that
does not follow that pattern gives an unreliable species. In a cluster of one
species `n_species` is 1. A larger value means the cluster joins several
species, or, for names such as the synthetic benchmark sets, that every genome
has its own species token. With no quality in the manifest (for example after
`ingest` without `--selection` columns) the four quality columns and
`best_member` are blank and the command logs that it left them blank.

To list the members of a cluster, see "Finding the members of a cluster" in
[usage.md](usage.md#inspecting-a-dereplication).

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
