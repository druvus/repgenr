// Select and organize viral genomes.

process VGENOME {
    label 'process_medium'
    tag "vgenome"

    input:
    val ready

    output:
    val true, emit: done

    script:
    """
    repgenr ${params.repgenr_opts} vgenome -wd ${params.workdir} ${params.vgenome_args}
    """
}
