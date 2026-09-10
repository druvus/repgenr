# CLI matrix

Generated from `tests/audit/cli_matrix.yaml` by `scripts/render_cli_matrix.py`;
`tests/unit/test_cli_matrix.py` keeps both in step with the command tree.

22 commands, 197 flags (197 with a live test or an n/a reason, 0 pending).

## Global flags

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--version` |  | n/a: prints the version and exits | none | n/a: prints the version and exits | n/a: unit test test_version_flag_prints_version | docs/adding-tools.md |
| `--container` |  | container.backend | choice | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--container-engine` |  | container.engine | none | params.repgenr_opts | tests/live/test_container_runs.py::test_container_engine_podman_is_reported_when_missing | docs/containers.md |
| `--container-cache` |  | container.cache_dir | none | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | docs/containers.md, docs/verification.md |
| `--platform` |  | container.platform | none | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | README.md, docs/containers.md |
| `--wave` |  | container.wave_enabled | none | params.repgenr_opts | tests/live/test_container_runs.py::test_native_result_is_not_reused_by_a_container_run | docs/adding-tools.md, docs/containers.md, docs/swot-viral.md |
| `--force` | -f | state.force | none | n/a: resume is Nextflow's -resume in the data-channel layer | tests/live/test_aux_commands.py::test_second_run_skips_and_force_reruns | README.md |
| `--verbose` | -v | state.log_level | none | params.repgenr_opts | tests/live/test_aux_commands.py::test_logging_flags_and_env | README.md |
| `--quiet` | -q | state.log_level | none | params.repgenr_opts | tests/live/test_aux_commands.py::test_logging_flags_and_env |  |

## derep-stock

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/containers.md, docs/output.md |
| `--action` |  | DerepStockParams.action | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip |  |
| `--name` |  | DerepStockParams.name | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/usage.md |

## derep-unpack

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_unpack_with_and_without_representant | docs/containers.md, docs/output.md |
| `--no-representant` |  | DerepUnpackParams.no_representant | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_unpack_with_and_without_representant |  |

## dereplicate

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--tool` |  | DereplicateParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_adapter_recovers_the_synthetic_partition | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | DereplicateParams.primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--secondary-ani` | -sani | DereplicateParams.secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_secondary_ani_sweep_changes_representative_count | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | DereplicateParams.aligned_fraction | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path |  |
| `--threads` | -t | DereplicateParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--process-size` | -s | DereplicateParams.process_size | none | params.derep_process_size | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/architecture.md, docs/scaling-audit.md, docs/swot-derep.md, docs/verification.md |
| `--num-processes` | -p | DereplicateParams.num_processes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--pre-primary-ani` |  | DereplicateParams.pre_primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--pre-secondary-ani` |  | DereplicateParams.pre_secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/verification.md |
| `--reduce` |  | DereplicateParams.reduce | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_reduce_species_keeps_one_representative_per_species | docs/scaling-audit.md, docs/swot-derep.md |
| `--target-reps` |  | DereplicateParams.target_reps | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_target_reps_lands_on_the_requested_count | docs/scaling-audit.md, docs/swot-derep.md |
| `--virus` |  | DereplicateParams.extra | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_virus_flag_on_auto_tool_is_reported_when_ignored | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | DereplicateParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_tool_arg_reaches_the_tool_command_line | docs/adding-tools.md, docs/usage.md |
| `--allow-incomplete` |  | DereplicateParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_allow_incomplete_gates_a_missing_genome |  |
| `--keeper` |  | DereplicateParams.keeper | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_keeper_quality_promotes_the_best_scored_member | README.md, docs/swot-derep.md, docs/usage.md |

## dereplicate-chunk

dispatch: `step:repgenr.stages.derep_steps.dereplicate_chunk`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--genomes-fofn` |  | ChunkParams.genomes | stage | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions |  |
| `--out` | -o | ChunkParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions |  |
| `--tool` |  | ChunkParams.tool | registry | params.derep_tool | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | ChunkParams.primary_ani | unit_interval | params.derep_primary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/verification.md |
| `--secondary-ani` | -sani | ChunkParams.secondary_ani | unit_interval | params.derep_secondary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | ChunkParams.aligned_fraction | unit_interval | params.derep_aligned_fraction | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) |  |
| `--threads` | -t | ChunkParams.threads | range | task.cpus | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/verification.md |
| `--virus` |  | ChunkParams.extra | none | params.mode | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | ChunkParams.extra | callback | n/a: not exposed; add it to ext.args in modules.config | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/adding-tools.md, docs/usage.md |
| `--selection-tsv` |  | ChunkParams.selection_tsv | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv |  |
| `--keeper` |  | ChunkParams.keeper | choice | params.derep_keeper | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | README.md, docs/swot-derep.md, docs/usage.md |
| `--versions-out` |  | ChunkParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions |  |

## dereplicate-merge

dispatch: `step:repgenr.stages.derep_steps.dereplicate_merge`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--out` | -o | MergeParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition |  |
| `--chunk-dir` |  | MergeParams.chunk_dirs | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition |  |
| `--chunk-fofn` |  | MergeParams.chunk_dirs | stage | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_fofn |  |
| `--tool` |  | MergeParams.tool | registry | params.derep_tool | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` | -pani | MergeParams.primary_ani | unit_interval | params.derep_primary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/verification.md |
| `--secondary-ani` | -sani | MergeParams.secondary_ani | unit_interval | params.derep_secondary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` | -af | MergeParams.aligned_fraction | unit_interval | params.derep_aligned_fraction | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) |  |
| `--threads` | -t | MergeParams.threads | range | task.cpus | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | docs/verification.md |
| `--virus` |  | MergeParams.extra | none | params.mode | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | README.md, docs/adding-tools.md, docs/swot-derep.md, docs/verification.md |
| `--tool-arg` |  | MergeParams.extra | callback | n/a: not exposed; add it to ext.args in modules.config | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/adding-tools.md, docs/usage.md |
| `--selection-tsv` |  | MergeParams.selection_tsv | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv |  |
| `--keeper` |  | MergeParams.keeper | choice | params.derep_keeper | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | README.md, docs/swot-derep.md, docs/usage.md |
| `--versions-out` |  | MergeParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_fofn |  |

## doctor

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_doctor_passes_then_fails_on_a_corrupt_genome | docs/containers.md, docs/output.md |

## genome

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_api_genus_representatives | docs/containers.md, docs/output.md |
| `--accession-list-only` |  | GenomeParams.accession_list_only | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_genome_accession_list_only_is_a_pure_query |  |
| `--keep-files` |  | GenomeParams.keep_files | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_genome_keep_files_retains_the_download_scratch |  |

## genome-fetch

dispatch: `step:repgenr.stages.genome_steps.genome_fetch`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--selection` |  | GenomeFetchParams.selection_tsv | stage | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step | README.md, docs/usage.md |
| `--out` | -o | GenomeFetchParams.out_dir | none | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step |  |
| `--keep-files` |  | GenomeFetchParams.keep_files | none | n/a: not exposed by the GENOME_FETCH module | tests/live/test_network.py::test_genome_fetch_step |  |
| `--versions-out` |  | GenomeFetchParams.versions_out | none | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step |  |

## glance

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | docs/containers.md, docs/output.md |
| `--tool` |  | GlanceParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--threads` | -t | GlanceParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | docs/verification.md |
| `--plot-max` |  | GlanceParams.plot_max | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: plot bounds only change the histogram; the PDF is checked for existence (test_glance_drep_compare) |  |
| `--plot-min` |  | GlanceParams.plot_min | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: plot bounds only change the histogram; the PDF is checked for existence (test_glance_drep_compare) |  |
| `--keep-files` |  | GlanceParams.keep_files | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare |  |

## ingest

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--genomes-dir` |  | IngestParams.genomes_dir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | README.md, docs/usage.md |
| `--selection` |  | IngestParams.selection | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_selection_table_drives_taxonomy_and_subset | README.md, docs/usage.md |
| `--outgroup` |  | IngestParams.outgroup | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_outgroup_and_copy | README.md, docs/usage.md |
| `--copy` |  | IngestParams.copy | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_outgroup_and_copy | README.md, docs/usage.md |

## list-tools

dispatch: `query`

(no flags)

## metadata

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_api_genus_representatives | docs/containers.md, docs/output.md |
| `--dataset` | -d | MetadataParams.dataset | choice | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup |  |
| `--level` | -l | MetadataParams.level | choice | params.metadata_args | tests/live/test_network.py::test_api_target_family_widens_the_selection |  |
| `--source` |  | MetadataParams.source | choice | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md, docs/usage.md, docs/verification.md |
| `--release` | -r | MetadataParams.release | none | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md |
| `--gtdb-version` |  | MetadataParams.version | none | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md, docs/usage.md |
| `--target-family` | -tf | MetadataParams.target_family | none | params.metadata_args | tests/live/test_network.py::test_api_target_family_widens_the_selection |  |
| `--target-genus` | -tg | MetadataParams.target_genus | none | params.metadata_args | tests/live/test_network.py::test_api_genus_representatives | README.md |
| `--target-species` | -ts | MetadataParams.target_species | none | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup |  |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup |  |
| `--metadata-path` |  | MetadataParams.metadata_path | none | params.metadata_args | tests/live/test_network.py::test_tsv_nodownload_and_metadata_path_reuse_the_table |  |
| `--nodownload` |  | MetadataParams.nodownload | none | params.metadata_args | tests/live/test_network.py::test_tsv_nodownload_and_metadata_path_reuse_the_table |  |
| `--limit` |  | MetadataParams.limit | range | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/scaling-audit.md, docs/swot-derep.md, docs/usage.md |

## phylo

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_treebuilders_offline.py::test_alignment_free_builder_on_representatives | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_iqtree_from_snptype_with_bootstrap_and_outgroup | README.md, docs/usage.md |
| `--aligner` |  | PhyloParams.aligner | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | README.md, docs/containers.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete |  |
| `--all-genomes` |  | PhyloParams.all_genomes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_treebuilders_offline.py::test_all_genomes_puts_every_genome_in_the_tree |  |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree |  |
| `--bootstrap` | -B | PhyloParams.bootstrap | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_iqtree_from_snptype_with_bootstrap_and_outgroup | docs/usage.md |
| `--reference` |  | PhyloParams.reference | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete | docs/swot-phylo.md |
| `--aligner-arg` |  | PhyloParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | docs/adding-tools.md |
| `--threads` | -t | PhyloParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | docs/verification.md |
| `--mask` |  | PhyloParams.extra | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_phylo_mask_gubbins | README.md, docs/adding-tools.md, docs/output.md, docs/usage.md, docs/verification.md |
| `--allow-incomplete` |  | PhyloParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete |  |

## phylo-build

dispatch: `step:repgenr.stages.phylo.phylo_build`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--genomes-dir` |  | PhyloBuildParams.genomes_dir | stage | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/usage.md |
| `--out` | -o | PhyloBuildParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--outgroup-dir` |  | PhyloBuildParams.outgroup_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--outgroup-accession` |  | PhyloBuildParams.outgroup_accession | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--treebuilder` |  | PhyloBuildParams.phylo.treebuilder | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloBuildParams.phylo.msa_source | choice | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | README.md, docs/usage.md |
| `--aligner` |  | PhyloBuildParams.phylo.aligner | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | README.md, docs/containers.md |
| `--snptyper` |  | PhyloBuildParams.phylo.snptyper | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants |  |
| `--no-outgroup` |  | PhyloBuildParams.phylo.no_outgroup | none | params.phylo_args | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree |  |
| `--bootstrap` | -B | PhyloBuildParams.phylo.bootstrap | range | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/usage.md |
| `--reference` |  | PhyloBuildParams.phylo.reference | none | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/swot-phylo.md |
| `--aligner-arg` |  | PhyloBuildParams.phylo.extra | callback | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/adding-tools.md |
| `--threads` | -t | PhyloBuildParams.phylo.threads | range | task.cpus | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/verification.md |
| `--versions-out` |  | PhyloBuildParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |

## run

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/containers.md, docs/output.md |
| `--viral` |  | VmetadataParams | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md |
| `--dataset` | -d | MetadataParams.dataset | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--level` | -l | MetadataParams.level | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--target-family` | -tf | MetadataParams.target_family | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets |  |
| `--target-genus` | -tg | MetadataParams.target_genus | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md |
| `--target-species` | -ts | MetadataParams.target_species | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets |  |
| `--release` | -r | MetadataParams.release | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to metadata unchanged (wiring test); the TSV path is exercised on metadata in test_network.py | README.md |
| `--gtdb-version` |  | MetadataParams.version | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to metadata unchanged (wiring test); the TSV path is exercised on metadata in test_network.py | README.md, docs/usage.md |
| `--metadata-source` |  | MetadataParams.source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--target` | -t | VmetadataParams.target | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md |
| `--viral-source` |  | VmetadataParams.source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end |  |
| `--group-segments` |  | VgenomeParams.group_segments | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to vgenome unchanged (wiring test); covered on vgenome in test_vgenome_selection_flags | README.md |
| `--tool` |  | DereplicateParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--primary-ani` |  | DereplicateParams.primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/verification.md |
| `--secondary-ani` |  | DereplicateParams.secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/scaling-audit.md, docs/usage.md, docs/verification.md |
| `--aligned-fraction` |  | DereplicateParams.aligned_fraction | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end |  |
| `--keeper` |  | DereplicateParams.keeper | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/swot-derep.md, docs/usage.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/containers.md, docs/usage.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: run forwards the phylo flags unchanged (test_species_set.py covers them on phylo) | README.md, docs/usage.md |
| `--aligner` |  | PhyloParams.aligner | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets | README.md, docs/containers.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: run forwards the phylo flags unchanged (test_species_set.py covers them on phylo) |  |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end |  |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md, docs/usage.md |
| `--threads` |  | DereplicateParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/verification.md |
| `--dry-run` |  | n/a: prints the chain and exits before any stage | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_prints_the_chain_without_network |  |

## snptype

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/containers.md, docs/output.md |
| `--tool` |  | SnptypeParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_parsnp_typer | README.md, docs/adding-tools.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/verification.md |
| `--reference` |  | SnptypeParams.reference | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/swot-phylo.md |
| `--all-genomes` |  | SnptypeParams.all_genomes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_on_representatives_only |  |
| `--mask` |  | SnptypeParams.mask | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_gubbins_mask_changes_the_core_alignment | README.md, docs/adding-tools.md, docs/output.md, docs/usage.md, docs/verification.md |
| `--threads` | -t | SnptypeParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/verification.md |
| `--tool-arg` |  | SnptypeParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_typer_and_tool_arg | docs/adding-tools.md, docs/usage.md |
| `--allow-incomplete` |  | SnptypeParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_snptype_allow_incomplete |  |

## status

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |

## tree2tax

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/containers.md, docs/output.md |
| `--node-basename` |  | Tree2taxParams.node_basename | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags |  |
| `--root-name` | -r | Tree2taxParams.root_name | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags |  |
| `--remove-outgroup` |  | Tree2taxParams.remove_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags |  |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | README.md, docs/usage.md |
| `--collapse-support` |  | Tree2taxParams.collapse_support | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/swot-phylo.md, docs/usage.md |
| `--collapse-length` |  | Tree2taxParams.collapse_length | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/swot-phylo.md, docs/usage.md |

## tree2tax-relations

dispatch: `step:repgenr.stages.tree2tax.tree2tax_relations`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--tree` |  | Tree2taxStepParams.tree | stage | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--out` | -o | Tree2taxStepParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--clusters` |  | Tree2taxStepParams.clusters | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--outgroup-dir` |  | Tree2taxStepParams.outgroup_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--outgroup-accession` |  | Tree2taxStepParams.outgroup_accession | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--node-basename` |  | Tree2taxStepParams.node_basename | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--root-name` | -r | Tree2taxStepParams.root_name | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--remove-outgroup` |  | Tree2taxStepParams.remove_outgroup | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--include-dereplicated` |  | Tree2taxStepParams.include_dereplicated | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/usage.md |
| `--versions-out` |  | Tree2taxStepParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup |  |
| `--collapse-support` |  | Tree2taxStepParams.collapse_support | range | params.tree2tax_args | tests/live/test_steps.py::test_tree2tax_relations_collapse_flags | docs/swot-phylo.md, docs/usage.md |
| `--collapse-length` |  | Tree2taxStepParams.collapse_length | range | params.tree2tax_args | tests/live/test_steps.py::test_tree2tax_relations_collapse_flags | docs/swot-phylo.md, docs/usage.md |

## versions

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | module: fixed by the process script | tests/live/test_aux_commands.py::test_status_and_versions | docs/containers.md, docs/output.md |
| `--versions-out` |  | n/a: query command, no stage parameters | none | module: fixed by the process script | tests/live/test_aux_commands.py::test_status_and_versions |  |

## vgenome

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_vgenome_selection_flags | docs/containers.md, docs/output.md |
| `--target-genus` | -tg | VgenomeParams.target_genus | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | README.md |
| `--target-species` | -ts | VgenomeParams.target_species | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--target-serotype` | -tse | VgenomeParams.target_serotype | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--target-custom` | -tc | VgenomeParams.target_custom | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--length-all` |  | VgenomeParams.length_all | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/usage.md |
| `--length-deviation` |  | VgenomeParams.length_deviation | range | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--length-method` |  | VgenomeParams.length_method | choice | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/scaling-audit.md, docs/swot-viral.md, docs/usage.md |
| `--length-range` |  | VgenomeParams.length_range | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/usage.md |
| `--discard` |  | VgenomeParams.discard | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |
| `--no-outgroup` |  | VgenomeParams.no_outgroup | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--group-segments` |  | VgenomeParams.group_segments | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | README.md |
| `--outgroup-candidates-taxid-min-genomes` |  | VgenomeParams.outgroup_candidates_taxid_min_genomes | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags |  |
| `--outgroup-treebuilder` |  | VgenomeParams.outgroup_treebuilder | registry | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |
| `--glance` |  | VgenomeParams.glance | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |
| `--print-fasta-headers` |  | VgenomeParams.print_fasta_headers | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |
| `--ignore-duplicates` |  | VgenomeParams.ignore_duplicates | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_bvbrc_needs_ignore_duplicates |  |
| `--keep-files` |  | VgenomeParams.keep_files | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files |  |

## vmetadata

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | docs/containers.md, docs/output.md |
| `--target` | -t | VmetadataParams.target | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | README.md |
| `--source` |  | VmetadataParams.source | choice | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_bvbrc_source_and_filter | README.md, docs/usage.md, docs/verification.md |
| `--filter` | -f | VmetadataParams.filter | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_bvbrc_source_and_filter |  |
| `--host` |  | VmetadataParams.host | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_released_after_and_host_narrow_the_set |  |
| `--complete-only` |  | VmetadataParams.complete_only | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only |  |
| `--released-after` |  | VmetadataParams.released_after | callback | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_released_after_and_host_narrow_the_set |  |
| `--list` | -l | VmetadataParams.list_targets | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_list_targets_reaches_bvbrc |  |

## Short-alias collisions

- `-t`: `--target`, `--threads`
- `-r`: `--release`, `--root-name`
- `-l`: `--level`, `--list`
- `-f`: `--filter`, `--force`
