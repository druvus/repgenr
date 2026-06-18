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
    main:
    ch_versions = Channel.empty()

    VACQUIRE()
    ch_versions = ch_versions.mix(VACQUIRE.out.versions)

    DEREPLICATE_SCATTER(VACQUIRE.out.genomes.flatten())
    ch_versions = ch_versions.mix(DEREPLICATE_SCATTER.out.versions)

    ch_reps     = DEREPLICATE_SCATTER.out.reps.map { meta, dir -> dir }
    ch_outgroup = VACQUIRE.out.outgroup.collect().ifEmpty([])
    ch_og_acc   = VACQUIRE.out.outgroup_accession

    PHYLO(ch_reps, ch_outgroup, ch_og_acc)
    ch_versions = ch_versions.mix(PHYLO.out.versions)

    TREE2TAX(PHYLO.out.tree, ch_reps, ch_outgroup, ch_og_acc)
    ch_versions = ch_versions.mix(TREE2TAX.out.versions)

    emit:
    tree        = PHYLO.out.tree
    tree2tax    = TREE2TAX.out.tree2tax
    genomes_map = TREE2TAX.out.genomes_map
    versions    = ch_versions
}
