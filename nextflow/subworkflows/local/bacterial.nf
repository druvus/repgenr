// Bacterial/archaeal pipeline: metadata -> genome -> dereplicate -> phylo -> tree2tax.
// Stages are ordered by passing a completion signal between processes; the
// shared RepGenR working directory (params.workdir) carries the data.

include { METADATA } from '../../modules/local/metadata'
include { GENOME } from '../../modules/local/genome'
include { DEREPLICATE } from '../../modules/local/dereplicate'
include { PHYLO } from '../../modules/local/phylo'
include { TREE2TAX } from '../../modules/local/tree2tax'

workflow BACTERIAL {
    main:
    METADATA()
    GENOME(METADATA.out.done)
    DEREPLICATE(GENOME.out.done)
    PHYLO(DEREPLICATE.out.done)
    TREE2TAX(PHYLO.out.done)

    emit:
    done = TREE2TAX.out.done
}
