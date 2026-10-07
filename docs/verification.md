# Adapter verification status

Status of each pluggable tool adapter against its real external binary. "Unit"
means covered by offline tests (mock/parse); "Live" means run end-to-end against
the installed tool on real or closely-related data.

Live verification used a small GTDB/NCBI **Francisella** set (real download) and
a closely-related **synthetic** set (for tools that need within-species
divergence). Tools were installed via conda/mamba on macOS (Apple Silicon).

## Re-running the checks: the live suite

The records below were collected by hand. `tests/live/` holds the
re-executable form: pytest tests marked `live` that run the installed
`repgenr` console script against real tools on seeded synthetic genome sets
(`benchmarks/genomegen.py`) staged with `repgenr ingest`. They are deselected
by default and never run in CI. Copy `tests/live/live.example.toml` to
`live.local.toml`, point it at the conda environments that hold each tool,
and run:

```bash
conda run -n repgenr_dev --no-capture-output \
    pytest tests/live -m live --live-config tests/live/live.local.toml -ra
```

Markers `network` and `container` select the tests that need GTDB/NCBI or
Docker; `-m "live and not network and not container"` is the offline subset.
Tests whose tool is not on PATH are skipped, not failed. See
`tests/live/README.md`.

## Adapter status

Where each built-in adapter has been run against its real tool, from the
live suite as it stands. "Native" means the tool on the host PATH of the
audit machine; "container" means through `--container docker` (a pinned
image or one minted by Wave). An adapter with neither has only the offline
contract test, which checks the argument vector against canned output.

| Family | Adapter | Native | Container | Live module |
|---|---|---|---|---|
| dereplicator | skder | yes | yes | test_dereplicators, test_container_runs |
| dereplicator | galah | yes | yes | test_dereplicators, test_container_runs (pinned image) |
| dereplicator | sourmash | yes | yes | test_dereplicators, test_aux_commands, test_container_runs (pinned image) |
| dereplicator | drep | no | yes | test_container_runs (pinned and Wave images; also `glance`) |
| aligner | sibeliaz | yes | yes | test_steps, test_container_runs |
| aligner | progressivemauve | no | yes | test_container_runs, test_nextflow (unpackaged on macOS) |
| aligner | cactus | no | yes | test_container_runs (pinned image) |
| SNP typer | simple | yes | yes | test_species_set, test_container_runs |
| SNP typer | parsnp | yes | no | test_species_set (osx-64 env) |
| SNP typer | ska2 | yes | yes | test_species_set, test_container_runs (pinned image) |
| SNP typer | snippy | no | yes | test_container_runs (pinned image) |
| masker | gubbins | yes | yes | test_species_set, test_container_runs (pinned image) |
| tree builder | iqtree | yes | yes | test_species_set, test_container_runs (pinned image) |
| tree builder | fasttree | yes | yes | test_species_set, test_container_runs (pinned image) |
| tree builder | raxmlng | yes | yes | test_species_set, test_container_runs (pinned image) |
| tree builder | mashtree | yes | yes | test_treebuilders_offline, test_species_set, test_container_runs (pinned image) |
| tree builder | sourmash | yes | yes | test_treebuilders_offline, test_container_runs (pinned image) |
| assembler | skesa | no | yes | test_container_runs (pinned image, simulated reads; a `requires_binary` test covers the host) |
| assembler | shovill | no | yes | test_container_runs (pinned image, simulated reads) |
| assembler | flye | no | yes | test_container_runs (pinned image, simulated 120 kb genome at 40x ONT-like reads, one contig) |
| classifier | sourmash | yes | yes | genome-qc and assemble on Wolbachia assemblies against the GTDB rs226 sketch, native and pinned image (2026-09-14) |
| polisher | medaka | no | yes | test_reads (SRR28800588 through the pinned image, see below); CheckM2 before/after recorded below |
| polisher | racon | no | no | offline contract test (minimap2 and racon argv, rounds, stdout capture) |
| quality | checkm2 | no | no | offline test on a canned quality report; live use needs the CheckM2 database (not on the audit machine) |

drep, progressivemauve, cactus and snippy have been verified only inside
containers. skder and SibeliaZ run in Wave-minted images only: their
BioContainer images are BusyBox-based and the GNU-only calls in their shell
wrappers fail there, which is why they carry no pinned image.

## Polishing an ONT assembly (2026-09-16)

SRR28800588 (Mycoplasmoides genitalium, MinION, 85 MB) assembled with Flye
and polished with medaka through the pinned images under amd64 emulation.
SRA rewrites FASTQ headers, so the reads named no basecaller; the adapter
used ONT's bacterial R10.4.1 model and recorded that assumption.

| Step | Result | Wall time |
|---|---|---|
| flye | 2 contigs, 593908 bp, N50 579957 | about 3 min |
| medaka (assumed bacterial model) | 2 contigs, 593932 bp | about 1 min |
| CheckM2, unpolished draft | 98.90 complete, 0.09 contamination | |
| CheckM2, polished | 98.90 complete, 0.09 contamination | 22 min for both, emulated |

On this run the Flye draft was already gene-complete at CheckM2 resolution
and polishing changed 24 bp of total length without moving the scores. The
12% contamination seen earlier on DRR351706 came from an MDA-amplified
library, which `reads --drop-selection` now excludes by default, not from a
lack of polishing. Polishing still matters for downstream SNP typing and
alignment, where indel errors are not absorbed by gene calling.

## Reads chain on a public run (2026-09-14)

`tests/live/test_reads.py` selects SRR25474756 (Mycoplasmopsis arginini,
Illumina MiSeq paired, 134 MB) from ENA, downloads and verifies both files,
and assembles them with SKESA in its pinned image under amd64 emulation on
the audit machine:

| Step | Result | Wall time |
|---|---|---|
| reads (ENA, Entrez lineage) | 1 run, labelled Metamycoplasmataceae / Mycoplasmopsis / arginini | 2 s |
| fetch (HTTPS, md5) | 134 MB | 80 s |
| skesa (4 threads, emulated) | 24 contigs, 670879 bp, N50 86129, 308x | 70 s |

The stateless steps behind the Nextflow reads mode (`assemble-run`,
`reads-gather`) are run on the same accession by
`test_reads_steps_assemble_a_public_run`; `genome-qc` has only the offline
tests with fakes, since the audit machine holds neither a CheckM2 database
nor a GTDB sketch. The Nextflow layer itself is covered by the nf-test stub
suite (`reads_select_process`, `reads_assemble_process`,
`genome_qc_process`, `reads_gather_process`, `acquire_reads`,
`reads_dataflow` and the `main.nf` reads-mode test) and, with real tools, by
`test_nextflow.py::test_main_reads_mode`, which runs `--mode reads` on
SRR25474756 (Illumina, skesa) and SRR28800588 (ONT, flye) through the pinned
images with the docker backend, then sourmash dereplication, a mashtree tree
and tree2tax (2026-09-14: two assembly tasks of 1m39s and 2m47s, the whole
run under seven minutes; two leaves in the tree). That run surfaced two
fixes: the steps resolve their paths before a containerised tool sees them,
and the container backend binds the real directory behind a symlinked
staging path.

## Results

The table below is rendered from the junit output of the last complete run
on the audit machine (Apple Silicon, Docker Desktop with Rosetta) by
`scripts/live_report.py`:

```bash
python scripts/live_report.py /path/to/junit.xml --out docs/verification.md
```

<!-- live-results:start -->
Last run 2026-09-10 13:11 UTC: 75 passed, 0 failed, 0 errors, 0 skipped, 25 min in total.

| module | test | result | seconds |
|---|---|---|---|
| test_aux_commands | test_derep_stock_round_trip | passed | 7 |
| test_aux_commands | test_derep_unpack_with_and_without_representant | passed | 4 |
| test_aux_commands | test_doctor_passes_then_fails_on_a_corrupt_genome | passed | 4 |
| test_aux_commands | test_logging_flags_and_env | passed | 12 |
| test_aux_commands | test_second_run_skips_and_force_reruns | passed | 9 |
| test_aux_commands | test_status_and_versions | passed | 4 |
| test_container_runs | test_cactus_pinned_image | passed | 129 |
| test_container_runs | test_container_engine_podman_is_reported_when_missing | passed | 1 |
| test_container_runs | test_drep_in_a_wave_container_with_virus_extra | passed | 9 |
| test_container_runs | test_glance_drep_compare | passed | 7 |
| test_container_runs | test_native_result_is_not_reused_by_a_container_run | passed | 8 |
| test_container_runs | test_progressivemauve_pinned_image | passed | 9 |
| test_container_runs | test_sibeliaz_in_a_wave_container | passed | 19 |
| test_container_runs | test_simple_typer_multi_tool_image | passed | 23 |
| test_container_runs | test_skder_in_a_wave_container_with_cache_and_env | passed | 7 |
| test_dereplicators | test_adapter_recovers_the_synthetic_partition[galah] | passed | 1 |
| test_dereplicators | test_adapter_recovers_the_synthetic_partition[skder] | passed | 5 |
| test_dereplicators | test_adapter_recovers_the_synthetic_partition[sourmash] | passed | 3 |
| test_dereplicators | test_allow_incomplete_gates_a_missing_genome | passed | 4 |
| test_dereplicators | test_keeper_quality_promotes_the_best_scored_member | passed | 7 |
| test_dereplicators | test_process_size_runs_the_chunked_path | passed | 6 |
| test_dereplicators | test_reduce_species_keeps_one_representative_per_species | passed | 3 |
| test_dereplicators | test_representatives_dir_matches_clusters | passed | 3 |
| test_dereplicators | test_secondary_ani_sweep_changes_representative_count | passed | 6 |
| test_dereplicators | test_target_reps_lands_on_the_requested_count | passed | 5 |
| test_dereplicators | test_tool_arg_reaches_the_tool_command_line | passed | 3 |
| test_dereplicators | test_virus_flag_on_auto_tool_is_reported_when_ignored | passed | 5 |
| test_ingest_flags | test_outgroup_and_copy | passed | 0 |
| test_ingest_flags | test_selection_table_drives_taxonomy_and_subset | passed | 0 |
| test_network | test_api_genus_representatives | passed | 0 |
| test_network | test_api_species_limit_and_explicit_outgroup | passed | 0 |
| test_network | test_api_target_family_widens_the_selection | passed | 13 |
| test_network | test_genome_accession_list_only_is_a_pure_query | passed | 1 |
| test_network | test_genome_fetch_step | passed | 9 |
| test_network | test_genome_keep_files_retains_the_download_scratch | passed | 9 |
| test_network | test_run_bacterial_chain_end_to_end | passed | 15 |
| test_network | test_run_dry_run_prints_the_chain_without_network | passed | 0 |
| test_network | test_run_dry_run_reports_family_and_species_targets | passed | 0 |
| test_network | test_run_viral_chain_end_to_end | passed | 21 |
| test_network | test_tsv_nodownload_and_metadata_path_reuse_the_table | passed | 32 |
| test_network | test_tsv_source_downloads_and_parses_the_release_table | passed | 0 |
| test_network | test_vgenome_bvbrc_needs_ignore_duplicates | passed | 2 |
| test_network | test_vgenome_discard_glance_headers_keep_files | passed | 2 |
| test_network | test_vgenome_selection_flags | passed | 6 |
| test_network | test_vmetadata_bvbrc_source_and_filter | passed | 5 |
| test_network | test_vmetadata_list_targets_reaches_bvbrc | passed | 2 |
| test_network | test_vmetadata_ncbi_virus_complete_only | passed | 0 |
| test_network | test_vmetadata_released_after_and_host_narrow_the_set | passed | 11 |
| test_nextflow | test_docker_profile_with_progressivemauve | passed | 17 |
| test_nextflow | test_local_dataflow_variants[skder---treebuilder sourmash] | passed | 17 |
| test_nextflow | test_local_dataflow_variants[sourmash---treebuilder mashtree] | passed | 12 |
| test_nextflow | test_main_bacterial_test_profile | passed | 18 |
| test_nextflow | test_main_viral_mode | passed | 30 |
| test_smoke | test_offline_chain_sourmash_mashtree | passed | 5 |
| test_species_set | test_fasttree_and_raxmlng_from_snptype | passed | 125 |
| test_species_set | test_gubbins_mask_changes_the_core_alignment | passed | 102 |
| test_species_set | test_iqtree_from_snptype_with_bootstrap_and_outgroup | passed | 127 |
| test_species_set | test_parsnp_typer | passed | 26 |
| test_species_set | test_phylo_mask_gubbins | passed | 124 |
| test_species_set | test_simple_typer_all_genomes_with_explicit_reference | passed | 47 |
| test_species_set | test_simple_typer_on_representatives_only | passed | 12 |
| test_species_set | test_ska2_source_with_reference_and_allow_incomplete | passed | 277 |
| test_species_set | test_ska2_typer_and_tool_arg | passed | 4 |
| test_species_set | test_snptype_allow_incomplete | passed | 34 |
| test_species_set | test_tree2tax_workdir_flags | passed | 10 |
| test_steps | test_chunk_keeper_quality_from_selection_tsv | passed | 9 |
| test_steps | test_chunk_results_carry_the_contract_and_versions | passed | 16 |
| test_steps | test_merge_by_chunk_dir_recovers_the_partition | passed | 19 |
| test_steps | test_merge_by_chunk_fofn | passed | 13 |
| test_steps | test_phylo_build_aligner_and_snp_source_variants | passed | 23 |
| test_steps | test_phylo_build_and_tree2tax_relations_with_outgroup | passed | 1 |
| test_steps | test_tree2tax_relations_collapse_flags | passed | 1 |
| test_treebuilders_offline | test_alignment_free_builder_on_representatives[mashtree] | passed | 5 |
| test_treebuilders_offline | test_alignment_free_builder_on_representatives[sourmash] | passed | 6 |
| test_treebuilders_offline | test_all_genomes_puts_every_genome_in_the_tree | passed | 4 |
<!-- live-results:end -->

What each module covers, with the flags it exercises, is in
`docs/audit/cli-matrix.md` (the `live` column) and `tests/live/README.md`.

## Genus-scale run (real data, 2026-09-10)

One `run` over every GTDB genome of a genus, on the audit machine (11 cores,
18 GB), workdir on an external volume:

```bash
repgenr run -wd work/francisella_all -d all -l genus -tg Francisella \
    --metadata-source api --tool skder --treebuilder mashtree -t 8 \
    --keeper quality --collapse-support 0.5
```

| Stage | Result | Wall time |
|---|---|---|
| metadata (API, quality cards) | 1157 genomes, 937 of them F. tularensis; outgroup GCF_003574425.1 | 4 min |
| genome | 1157 downloaded in one batch, 2.2 GB, none missing | 2.5 min |
| dereplicate (skder) | 68 representatives; largest cluster 893 genomes, 36 singletons; quality keeper replaced 11 | 5.7 min |
| phylo (mashtree, 8 threads) | 69 leaves, branch lengths 0.0001 to 0.22 (outgroup) | 1 s |
| tree2tax | 136 relations, 1158 genomes mapped | under 1 s |

12 minutes end to end, 25 CPU-minutes, 2.9 GB on disk, `doctor` clean.

## SNP phylogeny on the genus-scale set (2026-09-10)

The 68 Francisella representatives plus the outgroup, three ways, on the same
machine. `simple` maps every genome to one reference with minimap2 and calls
haploid consensus SNPs; the reference here was the alphabetically first genome,
a different species from most of the set.

| Path | Wall time | Sites | Outcome |
|---|---|---|---|
| mashtree (alignment-free) | 1 s | mash sketches | 69 leaves |
| ska2 + IQ-TREE, 8 threads | 7.4 min | 1737 variable | 32 zero-length branches, median UFBoot 75, 15 of 66 splits shared with mashtree |
| simple + Gubbins + IQ-TREE, 8 threads | 87 min, then failed | 812015 core SNP, 2054227 alignment columns | Gubbins crashed in its recombination scan |

The genus alignment is 39% variable. Gubbins is built for isolates of one
species; its scan allocates per-SNP arrays on the thread stack and dies with a
bus error on input this diverse, whatever tree builder it uses and however many
threads. Where the boundary sits, measured on subsets of the same alignment:

| Subset | Genomes | Variable columns | Gubbins |
|---|---|---|---|
| F. tularensis | 9 | 19723 (1.0%) | 5 iterations in 1.9 min, 3234 recombinant regions |
| F. philomiragia | 18 | 170342 (8.3%) | 5 iterations in 4.4 min, 8427 recombinant regions |
| tularensis + one philomiragia | 10 | 271022 (13.2%) | bus error |
| tularensis + philomiragia | 27 | 345308 (16.8%) | bus error |
| whole genus | 68 | 809060 (39.4%) | bus error |

The masker now estimates this fraction, warns above 10%, and reports the figure
if Gubbins fails. Recombination masking belongs to a within-species run; a
genus-level set is served by mashtree or an alignment-based path.

Masking changes the answer where it does run. On the nine F. tularensis
representatives, IQ-TREE with 1000 ultrafast bootstrap replicates took 27 s on
the 105597 unmasked variable sites and 21 s on the 19723 Gubbins kept. The two
topologies share no internal split (Robinson-Foulds 12 of a possible 12), and
both carry high support: median 98 unmasked, 96 masked. Note what the sites
rest on. These genomes were mapped to the alphabetically first genome of the
whole genus, another species, and 5.1% of that reference's columns vary within
a species usually described as clonal. A within-species run wants `--reference`
pointed at a genome of that species; the stage warns when it falls back.

Gubbins also needs a multi-threaded RAxML build (`raxmlHPC-PTHREADS*`) whenever
it is given more than one thread. The osx-arm64 conda package has none, so the
masker falls back to IQ-TREE, which cost 76 of the 87 minutes above on this
alignment. `--tool-arg gubbins_tree_builder=fasttree` with
`gubbins_first_tree_builder=rapidnj` did the same work in 6 minutes.

## SNP typing throughput (2026-09-11)

The `simple` typer on the 68 Francisella representatives, 8 threads, workdir on
the external volume. The stage maps each genome against the reference, then
reduces the stacked consensuses to variable columns.

| | Serial, plain VCF | Concurrent, compressed BCF |
|---|---|---|
| Mapping, calling and core-SNP reduction | 10.8 min | 2.8 min |
| Scratch left behind | 10 GB | 190 MB |
| Peak resident memory | not recorded | 0.4 GB |

The old scratch was 8.3 GB of uncompressed pileup VCF, which twice filled a
disk during this audit. What remains is the whole-genome alignment the stage
publishes for maskers.

## Command audit (2026-10-07)

Every `repgenr` subcommand was checked against its help text, argument
validation, an offline run on small local data, resume, the error surface
and the matching offline tests. The checks used the offline test suite and the
live offline subset (`-m "live and not network and not container"`), with no
network access and no Docker. Audit data lived on an external volume. The 28
commands were audited in six batches, one per help panel. In the first table a
cell is `ok` (agrees with the documentation), `fixed` (a defect was found and
corrected), `open` (a finding was left for the maintainer), `n/a` (the check
does not apply, for example resume for a stateless step or a live test that
does not exist) or `not run` (a live test exists but needs the network or a
container). The Result column takes the worst cell, in the order open, fixed,
ok, and also counts the docs agreement check.

| Command | Help | Validation | Offline run | Resume | Errors | Live subset | Result |
|---|---|---|---|---|---|---|---|
| run | ok | fixed | fixed | ok | ok | ok | fixed |
| status | ok | ok | ok | n/a | n/a | ok | fixed |
| metadata | ok | fixed | ok | fixed | fixed | not run | open |
| genome | ok | ok | ok | not run | open | not run | open |
| vmetadata | ok | ok | ok | ok | open | not run | open |
| vgenome | open | ok | ok | ok | ok | not run | open |
| ingest | ok | ok | ok | ok | ok | ok | ok |
| reads | ok | ok | ok | ok | ok | not run | ok |
| assemble | ok | ok | fixed | ok | fixed | not run | open |
| dereplicate | ok | fixed | fixed | ok | ok | ok | fixed |
| snptype | ok | ok | fixed | ok | ok | n/a | fixed |
| phylo | fixed | ok | fixed | ok | ok | ok | fixed |
| tree2tax | ok | fixed | fixed | fixed | ok | n/a | fixed |
| glance | ok | fixed | fixed | ok | open | not run | open |
| cluster-summary | ok | fixed | ok | open | ok | ok | open |
| derep-unpack | ok | ok | ok | open | ok | ok | open |
| derep-stock | ok | fixed | fixed | fixed | fixed | ok | open |
| list-tools | ok | n/a | fixed | n/a | ok | ok | fixed |
| doctor | ok | ok | ok | n/a | n/a | ok | ok |
| versions | ok | fixed | ok | n/a | n/a | ok | fixed |
| genome-fetch | ok | ok | fixed | n/a | fixed | n/a | fixed |
| dereplicate-chunk | ok | ok | ok | n/a | ok | ok | ok |
| dereplicate-merge | ok | ok | ok | n/a | ok | ok | ok |
| phylo-build | fixed | ok | ok | n/a | ok | ok | fixed |
| tree2tax-relations | ok | ok | fixed | n/a | fixed | ok | fixed |
| assemble-run | ok | ok | fixed | n/a | fixed | n/a | fixed |
| genome-qc | ok | ok | ok | n/a | ok | n/a | ok |
| reads-gather | ok | ok | fixed | n/a | ok | n/a | fixed |

Defects found and fixed, with the commit that fixed each. Several commits
change shared code (the stage runner, the contract readers and the Newick
parser), so they also apply to commands other than the one named.

| Command | Defect | Commit |
|---|---|---|
| run | `--with-snptype --snptyper bogus` passed a dry run and failed the real run with exit 5; it now exits 2 at once | 2722386 |
| run | `--with-snptype --mask gubbins` was rejected with exit 2; the mask now reaches the snptype stage | d16afbd |
| run, phylo | The `--mask` help on `run` did not mention `--with-snptype`; the bootstrap help said ">=1000" although smaller values are accepted | 3c087a7 |
| status | The "no run found" hint named only `metadata` and `vmetadata`; it now names all four starting points | 3c087a7 |
| versions | A nonexistent workdir printed nothing and exited 0; it now exits 3 and writes no fragment | 3608ad9 |
| status, doctor | A nonexistent workdir exited 0 (status printed the no-run hint, doctor a warning); both now exit 3 with the shared message, fixed after the audit at the maintainer's request | follow-up |
| list-tools | A rejected version flag (sibeliaz) was recorded as the tool's error line; it is now recorded as unknown | f1905d7 |
| metadata | The tsv source named `--version` (the global flag) instead of `--gtdb-version` when the version was missing | 5bb9f0c |
| metadata | Replacing the `--metadata-path` table was skipped on resume and kept a stale selection | 4337fab |
| vmetadata | An unreachable BV-BRC FTP server gave a traceback with exit 1; it now exits 3 with a named error | cbf839b |
| docs/usage.md | Exit code 3 was described only as workdir state; it also covers a failed remote request | 8e93630 |
| docs/output.md | `scratch/` was attributed to snptype only; it is written by five stages | 2a9f38e |
| genome-fetch | The `datasets` preflight was skipped without `--versions-out`, giving a traceback; it now exits 4 | 2762f39 |
| genome-fetch | A `selection.tsv` without the needed columns raised KeyError; it now exits 3 and names the columns | faef47b |
| assemble-run, reads-gather | A `reads.tsv` without the needed columns raised KeyError; it now exits 3 and names the columns | 6e7efc6 |
| assemble, assemble-run | A failure reason spanned several lines in `excused_runs.tsv`; each excused run now stays on one line | c091e9f |
| assemble, assemble-run | With `--assembler auto` and no assembler installed, runs were excused as `unsupported_platform`; they are now excused as `assembler_not_installed` with a warning, and the stage exits 4 only when nothing can be assembled | e7afca7, 68bef93 |
| dereplicate, tree2tax | A missing workdir gave a traceback or was created silently; the stage runner now exits 3 (it also covers glance, derep-unpack, cluster-summary, derep-stock, genome, vgenome and assemble) | c04fc22 |
| dereplicate | docs/output.md placed tool intermediates under `derep/`; they are under `scratch/` | 02f44ec |
| snptype | Switching typers left the previous typer's optional outputs (`full_alignment.fasta`, `snp_distance_matrix.tsv`, `variants.vcf`) in `snp/`; they are now dropped | 1429100 |
| phylo | Switching tree builders left the previous builder's files in `tree/`; they are now cleared | d2dfbb1 |
| phylo | The SibeliaZ macOS wrapper concatenated AppleDouble `._*` files into the alignment | a09ede2 |
| phylo | The SibeliaZ macOS wrapper left thousands of empty block temp files | 82dc1c0 |
| tree2tax | Toggling `--include-dereplicated` was logged as a changed input file | 2e01791 |
| tree2tax | A dereplicated member that is also a leaf was mapped twice in `genomes_map.tsv` | 4ffd3ee |
| tree2tax, tree2tax-relations | A malformed Newick tree gave a dendropy traceback; it now exits 3 | ddf48c6 |
| tree2tax-relations | A nonexistent `--clusters` path was skipped silently; it now exits 3 | a66b8ef |
| glance | A workdir without genomes exited 4 (dRep absent) instead of 3 | 645feaa |
| glance | The stage record held no dRep version | e003775 |
| glance | docs/usage.md said glance compares the representatives; it compares all genomes | 8cab3b7 |
| cluster-summary | A missing workdir created a manifest, or raised an OSError traceback | 50194c6 |
| derep-stock | `pack` of a workdir without dereplication outputs stored an empty run | f84e11f |
| derep-stock | `unpack` of an incomplete stored run emptied the live representatives before failing | 78d14cc |
| derep-stock | A repeat `unpack` after a new dereplication was skipped and restored nothing | fde2eb9 |
| derep-stock | `list` and `pack` on a nonexistent workdir exited 0 or created the workdir | bec185a |

Observations left for the maintainer. None changed a documented behaviour, so
they are recorded here and not fixed.

| Area | Observation |
|---|---|
| Resume | The resume fingerprint covers parameters, inputs and the environment, not outputs, so a stage whose output was deleted by hand is skipped and `--force` is needed. |
| Exit codes | When every assembly fails, `assemble` and `reads-gather` exit 3 and not 6; this is documented behaviour in the exit-code table of docs/usage.md, with the reasons in `excused_runs.tsv`. |
| Errors | A stage that fails cleanly leaves no record in `repgenr.yaml`, so `status` shows it as next and not interrupted, and `doctor` reports no failure while `tree/` holds partial files. |
| Errors | `phylo --treebuilder mashtree` on a single-representative set fails inside mashtree, and `snptype` with no variable sites exits 3; a genome-count check would give a clearer message. |
| Phylogeny | mashtree can drop degenerate genomes and exit 0, and no check compares the leaves of the tree with the input genomes. |
| Phylogeny | A `tree.nwk` with text after the final `;` is accepted by `tree2tax` but flagged as truncated by `doctor`. |
| Phylogeny | `tree2tax-relations` with an outgroup that is not a leaf logs a warning and leaves the tree unrooted, with exit 0. |
| Phylogeny | FastTree on the variable-site-only alignment of the `simple` typer gives branch lengths above one substitution per site, because there is no ascertainment correction; usage.md could say so. |
| Phylogeny | `phylo-build --msa-only` leaves `snp/` beside `msa.fasta`, and a ska2 run keeps its k-mer files in `scratch/snptype/`. |
| Records | The `tree2tax` record in `repgenr.yaml` has no tool and no versions, although output.md says every stage records its tool versions. |
| Records | `ingest` does not record `drop_foreign` in its parameters, and an outgroup row in a `--selection` is dropped when `--outgroup` names another genome. |
| Help text | `glance --tool` lists four dereplicators but only dRep supports comparison; `reads -tf/-tg/-ts` are not combined, since only the most specific is used; `tree2tax --node-basename` does not say that internal nodes otherwise get hash names. |
| Help text | The `reads --drop-selection` default is rendered as a Python list in the reference, and the no-match message of `reads` does not mention `--max-bases`. |
| reads | `--accession-file` treats only a `#` in column 1 as a comment, and a rejected invocation still creates the workdir and a log. |
| assemble | A missing CheckM2 result is kept with a warning and not excused, and `--polisher auto` with no polisher installed leaves ONT assemblies unpolished without a warning. |
| derep-unpack | A cluster member missing from `genomes/` is left out without a message, and a stored run without `cluster_summary.tsv` keeps the current summary. |
| derep-stock | Deleting an already deleted run exits 0 without naming the unknown run. |
| glance | Plots from an earlier run stay in place when no similarity falls in the plot range, and a dRep failure carries its full traceback in the error message. |
| Entry stages | `metadata --nodownload` reuses a table in the workdir that is not a declared resume input, so replacing it in place does not trigger a rerun. |
| Entry stages | The vmetadata NCBI Virus record omits `released_after` from its parameters, and the four entry records carry no tool, only tool versions. |
| Network | A BV-BRC group download writes `download.fa` in place, so an interrupted transfer can leave a partial file that the next run reuses. |
| Network | The BV-BRC path uses FTPS directly, so proxy variables do not block it, and `vmetadata --list` needs the network whatever `--source` says. |
| Network | With the network down, Entrez enrichment retries every sublist three times, about 16 minutes for 1050 taxids, before it fails with exit 3. |
| Help text | `metadata --metadata-path` does not say that `-r` and `--gtdb-version` are still required with a local table. |
| Tests | `--live-config <path>` into the main checkout from a worktree loads two conftest files and fails; `--live-config=<path>` works. |
| Environment | dRep 3.4.5 fails on exFAT volumes because macOS writes `._*` files into its cache; use an APFS workdir for glance and `dereplicate --tool drep`. |
| Environment | `status` and `doctor` on a long-running workdir (`francisella_all`) were not exercised, because that workdir was not on the audit machine. |

## Platform notes (macOS / Apple Silicon)

- Several tools lack osx-arm64 builds; some run via an osx-64 (Rosetta) conda env
  (e.g. parsnp/harvesttools). `mauve` is unpackaged on macOS entirely.
- `mashtree` pulls `perl-bio-samtools`, which pins an ancient samtools 0.1.x;
  modern samtools/bcftools for the `simple` SNP typer live in a separate env and
  are placed ahead on PATH.
- On exFAT/NTFS volumes, macOS `._*` AppleDouble files must be ignored (handled
  in the adapters/stages) and skDER must run on a local filesystem.
- amd64-only container tools run under emulation. Docker's QEMU emulation cannot
  run some SIMD-heavy binaries (e.g. Cactus's bundled `vg` hangs even on
  `vg version`). Set Docker Desktop to the **Apple Virtualization framework** with
  **"Use Rosetta for x86/amd64 emulation"** enabled; Rosetta runs these binaries
  correctly. On native amd64/Linux hosts no emulation is involved.
