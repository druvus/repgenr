// Build a phylogenetic tree (aligner or SNP source + tree builder).
// Cactus and other heavy aligners run under the process_high label.

process PHYLO {
    label 'process_high'
    tag "phylo"

    input:
    val ready

    output:
    val true, emit: done

    script:
    """
    repgenr phylo -wd ${params.workdir} ${params.phylo_args}
    """
}
