#!/usr/bin/env nextflow

// Network-free data-channel pipeline over a local directory of genome FASTAs:
//   DEREPLICATE_SCATTER -> PHYLO -> TREE2TAX
//
// Used by the real (non-stub) end-to-end CI test with lightweight tools
// (sourmash + mashtree), so it skips the GTDB/NCBI ACQUIRE front and runs with
// no outgroup. Publishes tree2tax.tsv, genomes_map.tsv and the tree to --outdir.

nextflow.enable.dsl = 2

include { DEREPLICATE_SCATTER } from '../subworkflows/local/dereplicate_scatter'
include { PHYLO               } from '../modules/local/dataflow/phylo'
include { TREE2TAX            } from '../modules/local/dataflow/tree2tax'

params.genomes_dir = null
params.empty_accession = "${projectDir}/tests/data/empty.txt"

workflow {
    if (!params.genomes_dir) {
        error "Provide --genomes_dir <DIR> with genome FASTAs."
    }
    ch_genomes = Channel
        .fromPath("${params.genomes_dir}/*.{fasta,fa,fna,fas}")
        .filter { !it.name.startsWith('._') }

    DEREPLICATE_SCATTER(ch_genomes)

    ch_reps     = DEREPLICATE_SCATTER.out.reps.map { meta, dir -> dir }
    ch_outgroup = Channel.value([])                       // no outgroup in the test
    ch_og_acc   = Channel.fromPath(params.empty_accession)

    PHYLO(ch_reps, ch_outgroup, ch_og_acc)
    TREE2TAX(PHYLO.out.tree, ch_reps, ch_outgroup, ch_og_acc)
}
