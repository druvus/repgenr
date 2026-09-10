# Command reference

Generated from the command tree by `scripts/render_cli_matrix.py`; the
matrix test keeps it in sync. Global options go before the command name
(`repgenr --container docker dereplicate ...`).

## Global options

| option | default | description |
|---|---|---|
| `--version` | off | Show version and exit. |
| `--container` | `none` | Run external tools in containers: none, docker, or singularity. |
| `--container-engine` |  | Engine binary override (e.g. apptainer, podman). |
| `--container-cache` |  | Directory for Singularity .sif images / Wave cache (large; can be external). |
| `--platform` |  | Container platform, e.g. linux/amd64 for emulated BioContainers on arm64. |
| `--wave`, `--no-wave` | off | Resolve images for multi-tool adapters via the Seqera Wave CLI. |
| `--force`, `-f`, `--no-force` | off | Re-run a stage even if it already completed with the same parameters. |
| `--verbose`, `-v` | off | Verbose (DEBUG) logging. |
| `--quiet`, `-q` | off | Only warnings and errors. |

## derep-stock

Store, load, list or delete named dereplication runs.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--action` | required | list, pack, unpack or delete. |
| `--name` |  | Run name for pack/unpack/delete. |

## derep-unpack

Explode clusters into one directory per representative.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--no-representant` | off | Leave the representative out of its cluster directory. |

## dereplicate

Cluster genomes by ANI and select representatives.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--tool` | `skder` | auto, drep, galah, skder, sourmash. |
| `-pani`, `--primary-ani` | `0.9` | Primary (pre-clustering) ANI threshold in (0, 1]. |
| `-sani`, `--secondary-ani` | `0.99` | Secondary (final cluster) ANI threshold in (0, 1]. |
| `-af`, `--aligned-fraction` | `0.5` | Minimum aligned fraction in (0, 1] for a pair to be compared. |
| `-t`, `--threads` | `16` | Threads for the external tool. |
| `-s`, `--process-size` |  | Chunk size; when set and exceeded, two-stage chunking runs for any tool. |
| `-p`, `--num-processes` | `0` | Parallel stage-1 chunk workers (threads split across them). 0 = auto (~threads/4, capped by cores). |
| `--pre-primary-ani` |  | Stage-1 (intra-chunk) primary ANI; defaults to --primary-ani. |
| `--pre-secondary-ani` |  | Stage-1 (intra-chunk) secondary ANI; defaults to --secondary-ani. |
| `--reduce` | `none` | Taxonomy-aware reduction after ANI: none, species, or genus (one representative per taxon). |
| `--target-reps` | `0` | Target representative count: search --secondary-ani to land near it (0 = off; re-runs dereplication per search step). |
| `--virus` | off | Pass virus-tuned parameters to the tool. |
| `--tool-arg` |  | Tool tuning as key=value (repeatable), e.g. mode=greedy. |
| `--allow-incomplete` | off | Proceed with a warning when genomes/ is missing selected genomes. |
| `--keeper` | `quality` | Representative choice per cluster: quality (CheckM score from GTDB) or tool (adapter's own). |

## dereplicate-chunk

Dereplicate one chunk of genomes (scatter step; writes a chunk result dir).

| option | default | description |
|---|---|---|
| `--genomes-fofn` | required | File of genome FASTA paths, one per line. |
| `-o`, `--out` | required | Output directory for the chunk result. |
| `--tool` | `skder` | drep, galah, skder, sourmash. |
| `-pani`, `--primary-ani` | `0.9` | Primary (pre-clustering) ANI threshold in (0, 1]. |
| `-sani`, `--secondary-ani` | `0.99` | Secondary (final cluster) ANI threshold in (0, 1]. |
| `-af`, `--aligned-fraction` | `0.5` | Minimum aligned fraction in (0, 1] for a pair to be compared. |
| `-t`, `--threads` | `16` | Threads for the external tool. |
| `--virus` | off | Pass virus-tuned parameters to the tool. |
| `--tool-arg` |  | Tool tuning as key=value (repeatable). |
| `--selection-tsv` |  | selection.tsv with quality columns; enables quality-aware representatives. |
| `--keeper` | `quality` | Representative choice when --selection-tsv is given: quality (manifest completeness/contamination) or tool (adapter's own pick). |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

## dereplicate-merge

Dereplicate the union of chunk representatives (gather step).

| option | default | description |
|---|---|---|
| `-o`, `--out` | required | Output dir for the merged result. |
| `--chunk-dir` |  | A chunk result directory (repeatable). |
| `--chunk-fofn` |  | File listing chunk result directories, one per line. |
| `--tool` | `skder` | drep, galah, skder, sourmash. |
| `-pani`, `--primary-ani` | `0.9` | Primary (pre-clustering) ANI threshold in (0, 1]. |
| `-sani`, `--secondary-ani` | `0.99` | Secondary (final cluster) ANI threshold in (0, 1]. |
| `-af`, `--aligned-fraction` | `0.5` | Minimum aligned fraction in (0, 1] for a pair to be compared. |
| `-t`, `--threads` | `16` | Threads for the external tool. |
| `--virus` | off | Pass virus-tuned parameters to the tool. |
| `--tool-arg` |  | Tool tuning as key=value (repeatable). |
| `--selection-tsv` |  | selection.tsv with quality columns; enables quality-aware representatives. |
| `--keeper` | `quality` | Representative choice when --selection-tsv is given: quality (manifest completeness/contamination) or tool (adapter's own pick). |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

## doctor

Verify a workdir's outputs against its records (read-only health check).

`status` reports what repgenr.yaml claims; `doctor` checks the claims
against the filesystem and the manifest: interrupted stages, missing or
corrupt genomes, manifest drift, representative/cluster mismatches,
truncated deliverables, and stages whose inputs changed since completion.
Exits 1 when any failure is found.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |

## genome

Download and organize genomes selected by the metadata stage.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--accession-list-only` | off | Write the accession list and stop (no download). |
| `--keep-files` | off | Keep download and scratch intermediates. |

## genome-fetch

Download genomes listed in a selection.tsv (stateless data-channel step).

| option | default | description |
|---|---|---|
| `--selection` | required | selection.tsv from the metadata stage. |
| `-o`, `--out` | required | Output dir for downloaded genomes. |
| `--keep-files` | off | Keep download intermediates. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

## glance

Quick all-vs-all ANI overview (dRep compare dendrogram + plots).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--tool` | `drep` | drep, galah, skder, sourmash. |
| `-t`, `--threads` | `16` | Threads for the external tool. |
| `--plot-max` | `1.0` | Upper similarity bound of the values plotted. |
| `--plot-min` | `0.0` | Lower similarity bound of the values plotted. |
| `--keep-files` | off | Keep download and scratch intermediates. |

## ingest

Populate a working directory from local genomes (no download).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `--genomes-dir` | required | Directory of genome FASTA files to stage under genomes/. |
| `--selection` |  | selection.tsv (accession, taxonomy, filename, outgroup flag, quality) naming the genomes to take; default: every FASTA under --genomes-dir, taxonomy parsed from canonical Family_genus_species_ACCESSION names. |
| `--outgroup` |  | Genome to set aside as the outgroup: a filename, stem or accession under --genomes-dir, or a path to a FASTA file elsewhere. |
| `--copy` | off | Copy the files into genomes/ instead of symlinking them. |

## list-tools

List the available pluggable tools in each family.

(no options)

## metadata

Select a taxon's genomes from GTDB (full table or the GTDB API).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `-d`, `--dataset` | required | all or rep. |
| `-l`, `--level` | required | family, genus or species. |
| `--source` | `tsv` | tsv (download full table) or api (GTDB API, target only). |
| `-r`, `--release` |  | GTDB release (tsv source). |
| `--gtdb-version` |  | bac120/ar53 (tsv source). |
| `-tf`, `--target-family` |  | Restrict the selection to this family. |
| `-tg`, `--target-genus` |  | Restrict the selection to this genus. |
| `-ts`, `--target-species` |  | Restrict the selection to this species. |
| `--outgroup-accession` |  | Accession to fetch and set aside as the outgroup. |
| `--metadata-path` |  | Use this GTDB metadata table instead of downloading. |
| `--nodownload` | off | Reuse a table already present in the workdir. |
| `--limit` |  | Keep at most N genomes, round-robin over species by CheckM quality. |

## phylo

Build a phylogenetic tree from an alignment, SNP alignment, or directly.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--treebuilder` | `iqtree` | auto, fasttree, iqtree, mashtree, raxmlng, sourmash. |
| `--msa-source` | `aligner` | aligner or snptype. |
| `--aligner` | `progressivemauve` | cactus, progressivemauve, sibeliaz. |
| `--snptyper` | `simple` | SNP typer: parsnp, simple, ska2, snippy. |
| `--all-genomes` | off | Use all genomes, not reps. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `-B`, `--bootstrap` | `0` | Bootstrap replicates (>=1000). |
| `--reference` |  | Reference genome filename. |
| `--aligner-arg` |  | Aligner tuning as key=value (repeatable), e.g. kmer=15 (sibeliaz) or seed_weight=11 (progressivemauve). |
| `-t`, `--threads` | `16` | Threads for the external tool. |
| `--mask` | `none` | Recombination masking for --msa-source snptype. |
| `--allow-incomplete` | off | Proceed with a warning when the input genome set is incomplete. |

## phylo-build

Build a phylogeny from a genomes directory (stateless data-channel step).

| option | default | description |
|---|---|---|
| `--genomes-dir` | required | Directory of genome FASTA files to build the tree from. |
| `-o`, `--out` | required | Output dir (writes tree/tree.nwk). |
| `--outgroup-dir` |  | Directory holding the outgroup genome file(s). |
| `--outgroup-accession` |  | File naming the outgroup accession. |
| `--treebuilder` | `iqtree` | auto, fasttree, iqtree, mashtree, raxmlng, sourmash. |
| `--msa-source` | `aligner` | aligner or snptype. |
| `--aligner` | `progressivemauve` | cactus, progressivemauve, sibeliaz. |
| `--snptyper` | `simple` | SNP typer: parsnp, simple, ska2, snippy. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `-B`, `--bootstrap` | `0` | Bootstrap replicates (>=1000). |
| `--reference` |  | Reference genome filename. |
| `--aligner-arg` |  | Aligner tuning as key=value (repeatable). |
| `-t`, `--threads` | `16` | Threads for the external tool. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |

## run

Run the whole pipeline end to end (bacterial by default, --viral for viruses).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `--viral` | off | Run the viral chain (vmetadata -> vgenome) instead of bacterial. |
| `-d`, `--dataset` | `rep` | all or rep (bacterial). |
| `-l`, `--level` |  | family/genus/species. |
| `-tf`, `--target-family` |  | Restrict the selection to this family. |
| `-tg`, `--target-genus` |  | Restrict the selection to this genus. |
| `-ts`, `--target-species` |  | Restrict the selection to this species. |
| `-r`, `--release` |  | GTDB release (tsv source). |
| `--gtdb-version` |  | bac120/ar53. |
| `--metadata-source` | `tsv` | tsv or api. |
| `--outgroup-accession` |  | Accession to fetch and set aside as the outgroup. |
| `-t`, `--target` |  | Virus taxon (viral). |
| `--viral-source` | `ncbi_virus` | ncbi_virus or bvbrc. |
| `--group-segments` | off | Group viral segments. |
| `--tool` | `skder` | auto, drep, galah, skder, sourmash. |
| `--primary-ani` | `0.9` | Primary (pre-clustering) ANI threshold in (0, 1]. |
| `--secondary-ani` | `0.99` | Secondary (final cluster) ANI threshold in (0, 1]. |
| `--aligned-fraction` | `0.5` | Minimum aligned fraction in (0, 1] for a pair to be compared. |
| `--keeper` | `quality` | Representative choice per cluster: quality (CheckM score from GTDB) or tool (adapter's own). |
| `--treebuilder` | `iqtree` | auto, fasttree, iqtree, mashtree, raxmlng, sourmash. |
| `--msa-source` | `aligner` | aligner or snptype. |
| `--aligner` | `progressivemauve` | cactus, progressivemauve, sibeliaz. |
| `--snptyper` | `simple` | SNP typer: parsnp, simple, ska2, snippy. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `--include-dereplicated`, `--no-include-dereplicated` | on | List redundant genomes under their representative in tree2tax. |
| `--threads` | `16` | Threads for the external tool. |
| `--dry-run` | off | Print the stages and key parameters, then exit. |

## snptype

Call SNPs and build a core-SNP alignment.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--tool` | `simple` | SNP typer: parsnp, simple, ska2, snippy. |
| `--reference` |  | Reference genome filename. |
| `--all-genomes` | off | Use all genomes, not reps. |
| `--mask` | `none` | Recombination masking: none, gubbins. |
| `-t`, `--threads` | `16` | Threads for the external tool. |
| `--tool-arg` |  | Tool tuning as key=value (repeatable). |
| `--allow-incomplete` | off | Proceed with a warning when the input genome set is incomplete. |

## status

Show which pipeline stages have completed in a working directory.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |

## tree2tax

Emit FlexTaxD-compatible taxonomy relations from the tree.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--node-basename` |  | Prefix for nodes. |
| `-r`, `--root-name` | `root` | Name for the root node. |
| `--remove-outgroup` | off | Drop outgroup. |
| `--include-dereplicated`, `--no-include-dereplicated` | on | List redundant genomes under their representative. |
| `--collapse-support` |  | Merge nodes whose support is below this fraction into their parent. |
| `--collapse-length` |  | Merge nodes whose branch is shorter than this length into their parent. |

## tree2tax-relations

Emit FlexTaxD relations from a tree (stateless data-channel step).

| option | default | description |
|---|---|---|
| `--tree` | required | Rooted/unrooted tree in Newick (tree.nwk). |
| `-o`, `--out` | required | Output dir (writes tree2tax.tsv + genomes_map.tsv). |
| `--clusters` |  | derep clusters.tsv (for --include-dereplicated). |
| `--outgroup-dir` |  | Directory holding the outgroup genome file(s). |
| `--outgroup-accession` |  | File naming the outgroup accession. |
| `--node-basename` |  | Prefix for nodes. |
| `-r`, `--root-name` | `root` | Name for the root node. |
| `--remove-outgroup` | off | Drop outgroup. |
| `--include-dereplicated` | off | List redundant genomes under their representative. |
| `--versions-out` |  | Write resolved tool versions (YAML fragment) here. |
| `--collapse-support` |  | Merge nodes whose support is below this fraction into their parent. |
| `--collapse-length` |  | Merge nodes whose branch is shorter than this length into their parent. |

## versions

Print the external-tool versions recorded in a workdir's repgenr.yaml.

Lets the Nextflow bridge modules (which run a full stage in a scratch workdir)
surface the resolved tool versions into versions.yml.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `--versions-out` |  | Write a versions.yml fragment here instead of stdout. |

## vgenome

Select and organize viral genomes (virus equivalent of genome).

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory. |
| `-tg`, `--target-genus` |  | Restrict the selection to this genus. |
| `-ts`, `--target-species` |  | Restrict the selection to this species. |
| `-tse`, `--target-serotype` |  | Restrict the selection to this serotype. |
| `-tc`, `--target-custom` |  | key:value. |
| `--length-all` | off | Disable the length window (keep every length). |
| `--length-deviation` | `10` | Half-width of the length window, in percent of its center. |
| `--length-method` | `median_of_medians` | Center of the length window: median_of_medians (default) gives one vote per species, so an over-sequenced outbreak species cannot shift the window; mean averages every record and loses that defense. |
| `--length-range` |  | e.g. 25000-35000. |
| `--discard` |  | Comma-separated header tags. |
| `--no-outgroup` | off | Do not root with an outgroup. |
| `--group-segments` | off | ncbi_virus: combine an isolate's segments into one genome (segmented viruses). |
| `--outgroup-candidates-taxid-min-genomes` | `5` | Genomes a sister taxid needs to qualify as an outgroup candidate. |
| `--outgroup-treebuilder` | `mashtree` | Tree builder used for the outgroup distance matrix. |
| `--glance` | off | Print selection and stop. |
| `--print-fasta-headers` | off | Print the headers of the selected records. |
| `--ignore-duplicates` | off | bvbrc: tolerate duplicate record ids (last wins). |
| `--keep-files` | off | Keep download and scratch intermediates. |

## vmetadata

Retrieve viral metadata from NCBI Virus (default) or BV-BRC.

| option | default | description |
|---|---|---|
| `-wd`, `--workdir` | required | Working directory (created). |
| `-t`, `--target` |  | Virus taxon/group/family. |
| `--source` | `ncbi_virus` | ncbi_virus (NCBI Virus via datasets) or bvbrc. |
| `-f`, `--filter` | `complete genome` | BV-BRC header tag. |
| `--host` |  | ncbi_virus: restrict to a host species. |
| `--complete-only` | off | ncbi_virus: only COMPLETE sequences. |
| `--released-after` |  | ncbi_virus: MM/DD/YYYY. |
| `-l`, `--list` | off | List BV-BRC targets and exit. |
