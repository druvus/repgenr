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

Gubbins also needs a multi-threaded RAxML build (`raxmlHPC-PTHREADS*`) whenever
it is given more than one thread. The osx-arm64 conda package has none, so the
masker falls back to IQ-TREE, which cost 76 of the 87 minutes above on this
alignment. `--tool-arg gubbins_tree_builder=fasttree` with
`gubbins_first_tree_builder=rapidnj` did the same work in 6 minutes.

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
