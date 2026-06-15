// Viral pipeline: vmetadata -> vgenome -> dereplicate -> phylo -> tree2tax.
// Shares the dereplicate/phylo/tree2tax modules with the bacterial path.

include { VMETADATA } from '../../modules/local/vmetadata'
include { VGENOME } from '../../modules/local/vgenome'
include { DEREPLICATE } from '../../modules/local/dereplicate'
include { PHYLO } from '../../modules/local/phylo'
include { TREE2TAX } from '../../modules/local/tree2tax'

workflow VIRAL {
    main:
    VMETADATA()
    VGENOME(VMETADATA.out.done)
    DEREPLICATE(VGENOME.out.done)
    PHYLO(DEREPLICATE.out.done)
    TREE2TAX(PHYLO.out.done)

    emit:
    done = TREE2TAX.out.done
}
