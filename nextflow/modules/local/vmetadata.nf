// Retrieve viral metadata from BV-BRC and NCBI.

process VMETADATA {
    label 'process_low'
    tag "vmetadata"

    output:
    val true, emit: done

    script:
    """
    repgenr ${params.repgenr_opts} vmetadata -wd ${params.workdir} ${params.vmetadata_args}
    """
}
