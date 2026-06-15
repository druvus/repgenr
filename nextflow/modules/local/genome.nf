// Download and organize genomes selected by the metadata stage.

process GENOME {
    label 'process_medium'
    tag "genome"

    input:
    val ready

    output:
    val true, emit: done

    script:
    """
    repgenr ${params.repgenr_opts} genome -wd ${params.workdir} ${params.genome_args}
    """
}
