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
    // Match threads to the reserved CPUs; a -t in phylo_args (last) overrides.
    """
    repgenr ${params.repgenr_opts} phylo -wd ${params.workdir} -t ${task.cpus} ${params.phylo_args}
    """
}
