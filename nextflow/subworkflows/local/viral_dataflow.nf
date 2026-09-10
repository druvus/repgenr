// Viral pipeline, data-channel form:
//   VACQUIRE (vmetadata -> vgenome) -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//
// Reuses the same data-channel dereplication, phylo and tree2tax as the bacterial
// path; only the acquisition front differs (BV-BRC instead of GTDB/NCBI).

include { VACQUIRE            } from '../../modules/local/dataflow/vacquire'
include { DEREPLICATE_SCATTER } from './dereplicate_scatter'
include { PHYLO               } from '../../modules/local/dataflow/phylo'
include { TREE2TAX            } from '../../modules/local/dataflow/tree2tax'

workflow VIRAL_DATAFLOW {
    take:
    ch_meta

    main:
    def ch_versions = channel.empty()

    VACQUIRE(ch_meta)
    ch_versions = ch_versions.mix(VACQUIRE.out.versions)

    def ch_genomes = VACQUIRE.out.genomes
        .map { meta, files -> tuple(meta, files instanceof List ? files : [files]) }
    // vmetadata writes no selection.tsv and viral genomes carry no CheckM
    // quality, so the keeper input is an empty list under the run meta.
    def ch_no_selection = ch_meta.map { meta -> tuple(meta, []) }
    DEREPLICATE_SCATTER(ch_genomes, ch_no_selection)
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    def ch_outgroup = ch_genomes
        // A placeholder tuple, not a bare meta: with remainder, an unmatched
        // bare value is emitted as-is and the next closure gets one argument
        // (seen live on a viral run with --no-outgroup).
        .map { meta, _files -> tuple(meta, []) }
        .join(VACQUIRE.out.outgroup, by: 0, remainder: true)
        .map { meta, _placeholder, files ->
            def list = files == null ? [] : (files instanceof List ? files : [files])
            tuple(meta, list)
        }
    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .join(ch_outgroup, by: 0)
        .join(VACQUIRE.out.outgroup_accession, by: 0)

    PHYLO(ch_phylo_in)
    ch_versions = ch_versions.mix(PHYLO.out.versions)

    // [meta, tree] joined with [meta, reps, outgroup, accession]
    def ch_tree2tax_in = PHYLO.out.tree.join(ch_phylo_in, by: 0)

    TREE2TAX(ch_tree2tax_in)
    ch_versions = ch_versions.mix(TREE2TAX.out.versions)

    emit:
    tree        = PHYLO.out.tree
    tree2tax    = TREE2TAX.out.tree2tax
    genomes_map = TREE2TAX.out.genomes_map
    versions    = ch_versions
}
