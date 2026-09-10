# CLI matrix

Generated from `tests/audit/cli_matrix.yaml` by `scripts/render_cli_matrix.py`;
`tests/unit/test_cli_matrix.py` keeps both in step with the command tree.

22 commands, 197 flags (197 with a live test or an n/a reason, 0 pending).

## Global flags

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--version` |  | n/a: prints the version and exits | none | n/a: prints the version and exits | n/a: unit test test_version_flag_prints_version | docs/adding-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--container` |  | container.backend | choice | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/containers.md, docs/swot-viral.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--container-engine` |  | container.engine | none | params.repgenr_opts | tests/live/test_container_runs.py::test_container_engine_podman_is_reported_when_missing | docs/cli-reference.md, docs/containers.md, docs/audit/cli-matrix.md |
| `--container-cache` |  | container.cache_dir | none | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | docs/cli-reference.md, docs/containers.md, docs/audit/cli-matrix.md |
| `--platform` |  | container.platform | none | params.repgenr_opts | tests/live/test_container_runs.py::test_skder_in_a_wave_container_with_cache_and_env | README.md, docs/cli-reference.md, docs/containers.md, docs/audit/cli-matrix.md |
| `--wave` |  | container.wave_enabled | none | params.repgenr_opts | tests/live/test_container_runs.py::test_native_result_is_not_reused_by_a_container_run | docs/adding-tools.md, docs/cli-reference.md, docs/containers.md, docs/swot-viral.md, docs/audit/cli-matrix.md |
| `--force` | -f | state.force | none | n/a: resume is Nextflow's -resume in the data-channel layer | tests/live/test_aux_commands.py::test_second_run_skips_and_force_reruns | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--verbose` | -v | state.log_level | none | params.repgenr_opts | tests/live/test_aux_commands.py::test_logging_flags_and_env | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--quiet` | -q | state.log_level | none | params.repgenr_opts | tests/live/test_aux_commands.py::test_logging_flags_and_env | docs/cli-reference.md, docs/audit/cli-matrix.md |

## derep-stock

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--action` |  | DerepStockParams.action | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--name` |  | DerepStockParams.name | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_stock_round_trip | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |

## derep-unpack

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_unpack_with_and_without_representant | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--no-representant` |  | DerepUnpackParams.no_representant | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_derep_unpack_with_and_without_representant | docs/cli-reference.md, docs/audit/cli-matrix.md |

## dereplicate

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--tool` |  | DereplicateParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_adapter_recovers_the_synthetic_partition | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--primary-ani` | -pani | DereplicateParams.primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--secondary-ani` | -sani | DereplicateParams.secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_secondary_ani_sweep_changes_representative_count | docs/cli-reference.md, docs/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligned-fraction` | -af | DereplicateParams.aligned_fraction | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | DereplicateParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--process-size` | -s | DereplicateParams.process_size | none | params.derep_process_size | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/architecture.md, docs/cli-reference.md, docs/scaling-audit.md, docs/swot-derep.md, docs/audit/cli-matrix.md |
| `--num-processes` | -p | DereplicateParams.num_processes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--pre-primary-ani` |  | DereplicateParams.pre_primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--pre-secondary-ani` |  | DereplicateParams.pre_secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_process_size_runs_the_chunked_path | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--reduce` |  | DereplicateParams.reduce | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_reduce_species_keeps_one_representative_per_species | docs/cli-reference.md, docs/scaling-audit.md, docs/swot-derep.md, docs/audit/cli-matrix.md |
| `--target-reps` |  | DereplicateParams.target_reps | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_target_reps_lands_on_the_requested_count | docs/cli-reference.md, docs/scaling-audit.md, docs/swot-derep.md, docs/audit/cli-matrix.md |
| `--virus` |  | DereplicateParams.extra | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_virus_flag_on_auto_tool_is_reported_when_ignored | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/swot-derep.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | DereplicateParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_tool_arg_reaches_the_tool_command_line | docs/adding-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--allow-incomplete` |  | DereplicateParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_allow_incomplete_gates_a_missing_genome | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keeper` |  | DereplicateParams.keeper | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_dereplicators.py::test_keeper_quality_promotes_the_best_scored_member | README.md, docs/cli-reference.md, docs/swot-derep.md, docs/usage.md, docs/audit/cli-matrix.md |

## dereplicate-chunk

dispatch: `step:repgenr.stages.derep_steps.dereplicate_chunk`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--genomes-fofn` |  | ChunkParams.genomes | stage | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | ChunkParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--tool` |  | ChunkParams.tool | registry | params.derep_tool | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--primary-ani` | -pani | ChunkParams.primary_ani | unit_interval | params.derep_primary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--secondary-ani` | -sani | ChunkParams.secondary_ani | unit_interval | params.derep_secondary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligned-fraction` | -af | ChunkParams.aligned_fraction | unit_interval | params.derep_aligned_fraction | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | ChunkParams.threads | range | task.cpus | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--virus` |  | ChunkParams.extra | none | params.mode | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/swot-derep.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | ChunkParams.extra | callback | n/a: not exposed; add it to ext.args in modules.config | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/adding-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--selection-tsv` |  | ChunkParams.selection_tsv | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keeper` |  | ChunkParams.keeper | choice | params.derep_keeper | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | README.md, docs/cli-reference.md, docs/swot-derep.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | ChunkParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_results_carry_the_contract_and_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |

## dereplicate-merge

dispatch: `step:repgenr.stages.derep_steps.dereplicate_merge`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--out` | -o | MergeParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--chunk-dir` |  | MergeParams.chunk_dirs | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--chunk-fofn` |  | MergeParams.chunk_dirs | stage | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_fofn | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool` |  | MergeParams.tool | registry | params.derep_tool | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--primary-ani` | -pani | MergeParams.primary_ani | unit_interval | params.derep_primary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--secondary-ani` | -sani | MergeParams.secondary_ani | unit_interval | params.derep_secondary_ani | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligned-fraction` | -af | MergeParams.aligned_fraction | unit_interval | params.derep_aligned_fraction | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | MergeParams.threads | range | task.cpus | tests/live/test_steps.py::test_merge_by_chunk_dir_recovers_the_partition | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--virus` |  | MergeParams.extra | none | params.mode | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/swot-derep.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | MergeParams.extra | callback | n/a: not exposed; add it to ext.args in modules.config | n/a: shared adapter path, exercised through dereplicate (test_dereplicators.py) | docs/adding-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--selection-tsv` |  | MergeParams.selection_tsv | none | module: fixed by the process script | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keeper` |  | MergeParams.keeper | choice | params.derep_keeper | tests/live/test_steps.py::test_chunk_keeper_quality_from_selection_tsv | README.md, docs/cli-reference.md, docs/swot-derep.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | MergeParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_merge_by_chunk_fofn | docs/cli-reference.md, docs/audit/cli-matrix.md |

## doctor

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_aux_commands.py::test_doctor_passes_then_fails_on_a_corrupt_genome | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |

## genome

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_api_genus_representatives | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--accession-list-only` |  | GenomeParams.accession_list_only | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_genome_accession_list_only_is_a_pure_query | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | GenomeParams.keep_files | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_genome_keep_files_retains_the_download_scratch | docs/cli-reference.md, docs/audit/cli-matrix.md |

## genome-fetch

dispatch: `step:repgenr.stages.genome_steps.genome_fetch`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--selection` |  | GenomeFetchParams.selection_tsv | stage | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--out` | -o | GenomeFetchParams.out_dir | none | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | GenomeFetchParams.keep_files | none | n/a: not exposed by the GENOME_FETCH module | tests/live/test_network.py::test_genome_fetch_step | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | GenomeFetchParams.versions_out | none | module: fixed by the process script | tests/live/test_network.py::test_genome_fetch_step | docs/cli-reference.md, docs/audit/cli-matrix.md |

## glance

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--tool` |  | GlanceParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--threads` | -t | GlanceParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--plot-max` |  | GlanceParams.plot_max | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: plot bounds only change the histogram; the PDF is checked for existence (test_glance_drep_compare) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--plot-min` |  | GlanceParams.plot_min | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: plot bounds only change the histogram; the PDF is checked for existence (test_glance_drep_compare) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | GlanceParams.keep_files | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_glance_drep_compare | docs/cli-reference.md, docs/audit/cli-matrix.md |

## ingest

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--genomes-dir` |  | IngestParams.genomes_dir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--selection` |  | IngestParams.selection | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_selection_table_drives_taxonomy_and_subset | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--outgroup` |  | IngestParams.outgroup | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_outgroup_and_copy | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--copy` |  | IngestParams.copy | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_ingest_flags.py::test_outgroup_and_copy | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |

## list-tools

dispatch: `query`

(no flags)

## metadata

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_api_genus_representatives | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--dataset` | -d | MetadataParams.dataset | choice | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--level` | -l | MetadataParams.level | choice | params.metadata_args | tests/live/test_network.py::test_api_target_family_widens_the_selection | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--source` |  | MetadataParams.source | choice | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--release` | -r | MetadataParams.release | none | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-version` |  | MetadataParams.version | none | params.metadata_args | tests/live/test_network.py::test_tsv_source_downloads_and_parses_the_release_table | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--target-family` | -tf | MetadataParams.target_family | none | params.metadata_args | tests/live/test_network.py::test_api_target_family_widens_the_selection | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-genus` | -tg | MetadataParams.target_genus | none | params.metadata_args | tests/live/test_network.py::test_api_genus_representatives | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-species` | -ts | MetadataParams.target_species | none | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--metadata-path` |  | MetadataParams.metadata_path | none | params.metadata_args | tests/live/test_network.py::test_tsv_nodownload_and_metadata_path_reuse_the_table | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--nodownload` |  | MetadataParams.nodownload | none | params.metadata_args | tests/live/test_network.py::test_tsv_nodownload_and_metadata_path_reuse_the_table | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--limit` |  | MetadataParams.limit | range | params.metadata_args | tests/live/test_network.py::test_api_species_limit_and_explicit_outgroup | docs/cli-reference.md, docs/scaling-audit.md, docs/swot-derep.md, docs/usage.md, docs/audit/cli-matrix.md |

## phylo

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_treebuilders_offline.py::test_alignment_free_builder_on_representatives | README.md, docs/cli-reference.md, docs/containers.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_iqtree_from_snptype_with_bootstrap_and_outgroup | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligner` |  | PhyloParams.aligner | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | README.md, docs/cli-reference.md, docs/containers.md, docs/audit/cli-matrix.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--all-genomes` |  | PhyloParams.all_genomes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_treebuilders_offline.py::test_all_genomes_puts_every_genome_in_the_tree | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--bootstrap` | -B | PhyloParams.bootstrap | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_iqtree_from_snptype_with_bootstrap_and_outgroup | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--reference` |  | PhyloParams.reference | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete | docs/cli-reference.md, docs/swot-phylo.md, docs/audit/cli-matrix.md |
| `--aligner-arg` |  | PhyloParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | docs/adding-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | PhyloParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_container_runs.py::test_sibeliaz_in_a_wave_container | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--mask` |  | PhyloParams.extra | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_phylo_mask_gubbins | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/output.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--allow-incomplete` |  | PhyloParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_source_with_reference_and_allow_incomplete | docs/cli-reference.md, docs/audit/cli-matrix.md |

## phylo-build

dispatch: `step:repgenr.stages.phylo.phylo_build`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--genomes-dir` |  | PhyloBuildParams.genomes_dir | stage | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--out` | -o | PhyloBuildParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--outgroup-dir` |  | PhyloBuildParams.outgroup_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | PhyloBuildParams.outgroup_accession | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--treebuilder` |  | PhyloBuildParams.phylo.treebuilder | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/cli-reference.md, docs/containers.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--msa-source` |  | PhyloBuildParams.phylo.msa_source | choice | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligner` |  | PhyloBuildParams.phylo.aligner | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | README.md, docs/cli-reference.md, docs/containers.md, docs/audit/cli-matrix.md |
| `--snptyper` |  | PhyloBuildParams.phylo.snptyper | registry | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | PhyloBuildParams.phylo.no_outgroup | none | params.phylo_args | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--bootstrap` | -B | PhyloBuildParams.phylo.bootstrap | range | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--reference` |  | PhyloBuildParams.phylo.reference | none | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/cli-reference.md, docs/swot-phylo.md, docs/audit/cli-matrix.md |
| `--aligner-arg` |  | PhyloBuildParams.phylo.extra | callback | params.phylo_args | tests/live/test_steps.py::test_phylo_build_aligner_and_snp_source_variants | docs/adding-tools.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--threads` | -t | PhyloBuildParams.phylo.threads | range | task.cpus | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | PhyloBuildParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |

## run

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--viral` |  | VmetadataParams | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--dataset` | -d | MetadataParams.dataset | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--level` | -l | MetadataParams.level | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-family` | -tf | MetadataParams.target_family | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-genus` | -tg | MetadataParams.target_genus | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-species` | -ts | MetadataParams.target_species | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--release` | -r | MetadataParams.release | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to metadata unchanged (wiring test); the TSV path is exercised on metadata in test_network.py | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--gtdb-version` |  | MetadataParams.version | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to metadata unchanged (wiring test); the TSV path is exercised on metadata in test_network.py | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--metadata-source` |  | MetadataParams.source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | MetadataParams.outgroup_accession | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target` |  | VmetadataParams.target | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--viral-source` |  | VmetadataParams.source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--group-segments` |  | VgenomeParams.group_segments | none | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: forwarded to vgenome unchanged (wiring test); covered on vgenome in test_vgenome_selection_flags | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool` |  | DereplicateParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--primary-ani` |  | DereplicateParams.primary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--secondary-ani` |  | DereplicateParams.secondary_ani | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/scaling-audit.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligned-fraction` |  | DereplicateParams.aligned_fraction | unit_interval | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keeper` |  | DereplicateParams.keeper | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/cli-reference.md, docs/swot-derep.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--treebuilder` |  | PhyloParams.treebuilder | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | README.md, docs/cli-reference.md, docs/containers.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--msa-source` |  | PhyloParams.msa_source | choice | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: run forwards the phylo flags unchanged (test_species_set.py covers them on phylo) | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--aligner` |  | PhyloParams.aligner | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_reports_family_and_species_targets | README.md, docs/cli-reference.md, docs/containers.md, docs/audit/cli-matrix.md |
| `--snptyper` |  | PhyloParams.snptyper | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | n/a: run forwards the phylo flags unchanged (test_species_set.py covers them on phylo) | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | PhyloParams.no_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_viral_chain_end_to_end | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--threads` | -t | DereplicateParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_bacterial_chain_end_to_end | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--dry-run` |  | n/a: prints the chain and exits before any stage | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_network.py::test_run_dry_run_prints_the_chain_without_network | docs/cli-reference.md, docs/audit/cli-matrix.md |

## snptype

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--tool` |  | SnptypeParams.tool | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_parsnp_typer | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/containers.md, docs/swot-derep.md, docs/swot-phylo.md, docs/swot-viral.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--reference` |  | SnptypeParams.reference | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/cli-reference.md, docs/swot-phylo.md, docs/audit/cli-matrix.md |
| `--all-genomes` |  | SnptypeParams.all_genomes | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_on_representatives_only | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--mask` |  | SnptypeParams.mask | registry | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_gubbins_mask_changes_the_core_alignment | README.md, docs/adding-tools.md, docs/cli-reference.md, docs/output.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--threads` | -t | SnptypeParams.threads | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_simple_typer_all_genomes_with_explicit_reference | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--tool-arg` |  | SnptypeParams.extra | callback | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_ska2_typer_and_tool_arg | docs/adding-tools.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--allow-incomplete` |  | SnptypeParams.allow_incomplete | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_snptype_allow_incomplete | docs/cli-reference.md, docs/audit/cli-matrix.md |

## status

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |

## tree2tax

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--node-basename` |  | Tree2taxParams.node_basename | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--root-name` |  | Tree2taxParams.root_name | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--remove-outgroup` |  | Tree2taxParams.remove_outgroup | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--include-dereplicated` |  | Tree2taxParams.include_dereplicated | none | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_smoke.py::test_offline_chain_sourmash_mashtree | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--collapse-support` |  | Tree2taxParams.collapse_support | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/swot-phylo.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--collapse-length` |  | Tree2taxParams.collapse_length | range | n/a: workdir command; the Nextflow layer uses the stateless steps | tests/live/test_species_set.py::test_tree2tax_workdir_flags | docs/cli-reference.md, docs/swot-phylo.md, docs/usage.md, docs/audit/cli-matrix.md |

## tree2tax-relations

dispatch: `step:repgenr.stages.tree2tax.tree2tax_relations`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--tree` |  | Tree2taxStepParams.tree | stage | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--out` | -o | Tree2taxStepParams.out_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/verification.md, docs/audit/cli-matrix.md |
| `--clusters` |  | Tree2taxStepParams.clusters | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-dir` |  | Tree2taxStepParams.outgroup_dir | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-accession` |  | Tree2taxStepParams.outgroup_accession | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--node-basename` |  | Tree2taxStepParams.node_basename | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--root-name` |  | Tree2taxStepParams.root_name | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--remove-outgroup` |  | Tree2taxStepParams.remove_outgroup | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--include-dereplicated` |  | Tree2taxStepParams.include_dereplicated | none | params.tree2tax_args | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | Tree2taxStepParams.versions_out | none | module: fixed by the process script | tests/live/test_steps.py::test_phylo_build_and_tree2tax_relations_with_outgroup | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--collapse-support` |  | Tree2taxStepParams.collapse_support | range | params.tree2tax_args | tests/live/test_steps.py::test_tree2tax_relations_collapse_flags | docs/cli-reference.md, docs/swot-phylo.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--collapse-length` |  | Tree2taxStepParams.collapse_length | range | params.tree2tax_args | tests/live/test_steps.py::test_tree2tax_relations_collapse_flags | docs/cli-reference.md, docs/swot-phylo.md, docs/usage.md, docs/audit/cli-matrix.md |

## versions

dispatch: `query`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | n/a: query command, no stage parameters | none | module: fixed by the process script | tests/live/test_aux_commands.py::test_status_and_versions | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--versions-out` |  | n/a: query command, no stage parameters | none | module: fixed by the process script | tests/live/test_aux_commands.py::test_status_and_versions | docs/cli-reference.md, docs/audit/cli-matrix.md |

## vgenome

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--target-genus` | -tg | VgenomeParams.target_genus | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-species` | -ts | VgenomeParams.target_species | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-serotype` | -tse | VgenomeParams.target_serotype | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--target-custom` | -tc | VgenomeParams.target_custom | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--length-all` |  | VgenomeParams.length_all | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--length-deviation` |  | VgenomeParams.length_deviation | range | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--length-method` |  | VgenomeParams.length_method | choice | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/scaling-audit.md, docs/swot-viral.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--length-range` |  | VgenomeParams.length_range | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--discard` |  | VgenomeParams.discard | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--no-outgroup` |  | VgenomeParams.no_outgroup | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--group-segments` |  | VgenomeParams.group_segments | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-candidates-taxid-min-genomes` |  | VgenomeParams.outgroup_candidates_taxid_min_genomes | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_selection_flags | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--outgroup-treebuilder` |  | VgenomeParams.outgroup_treebuilder | registry | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--glance` |  | VgenomeParams.glance | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--print-fasta-headers` |  | VgenomeParams.print_fasta_headers | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--ignore-duplicates` |  | VgenomeParams.ignore_duplicates | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_bvbrc_needs_ignore_duplicates | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--keep-files` |  | VgenomeParams.keep_files | none | params.vgenome_args | tests/live/test_network.py::test_vgenome_discard_glance_headers_keep_files | docs/cli-reference.md, docs/audit/cli-matrix.md |

## vmetadata

dispatch: `stage`

| flag | aliases | param | validated | nextflow | live | docs |
|---|---|---|---|---|---|---|
| `--workdir` | -wd | workdir | stage | module: fixed by the process script | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | docs/cli-reference.md, docs/containers.md, docs/output.md, docs/audit/cli-matrix.md |
| `--target` |  | VmetadataParams.target | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | README.md, docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--source` |  | VmetadataParams.source | choice | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_bvbrc_source_and_filter | README.md, docs/cli-reference.md, docs/usage.md, docs/audit/cli-matrix.md |
| `--filter` |  | VmetadataParams.filter | stage | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_bvbrc_source_and_filter | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--host` |  | VmetadataParams.host | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_released_after_and_host_narrow_the_set | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--complete-only` |  | VmetadataParams.complete_only | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_ncbi_virus_complete_only | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--released-after` |  | VmetadataParams.released_after | callback | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_released_after_and_host_narrow_the_set | docs/cli-reference.md, docs/audit/cli-matrix.md |
| `--list` |  | VmetadataParams.list_targets | none | params.vmetadata_args | tests/live/test_network.py::test_vmetadata_list_targets_reaches_bvbrc | docs/cli-reference.md, docs/audit/cli-matrix.md |

## Short-alias collisions

