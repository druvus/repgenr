// Emit FlexTaxD-compatible taxonomy relations from the tree.

process TREE2TAX {
    label 'process_low'
    tag "tree2tax"

    input:
    val ready

    output:
    val true, emit: done

    script:
    """
    repgenr tree2tax -wd ${params.workdir} ${params.tree2tax_args}
    """
}
