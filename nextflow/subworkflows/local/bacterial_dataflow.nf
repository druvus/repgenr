// Bacterial pipeline, data-channel form:
//   ACQUIRE (metadata -> genome) -> DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//
// Genomes, representatives, the tree and the taxonomy flow between processes as
// staged channel files; there is no shared working directory. This is the
// data-channel replacement for the legacy done-signal BACTERIAL subworkflow.

include { ACQUIRE             } from './acquire'
include { DEREPLICATE_SCATTER } from './dereplicate_scatter'
include { PHYLO               } from '../../modules/local/dataflow/phylo'
include { TREE2TAX            } from '../../modules/local/dataflow/tree2tax'

workflow BACTERIAL_DATAFLOW {
    take:
    ch_meta

    main:
    def ch_versions = channel.empty()

    ACQUIRE(ch_meta)
    ch_versions = ch_versions.mix(ACQUIRE.out.versions)

    DEREPLICATE_SCATTER(ACQUIRE.out.genomes, ACQUIRE.out.selection)
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    def ch_phylo_in = DEREPLICATE_SCATTER.out.reps
        .join(ACQUIRE.out.outgroup, by: 0)
        .join(ACQUIRE.out.outgroup_accession, by: 0)

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
