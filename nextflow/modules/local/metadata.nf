// Select a taxon's genomes from GTDB metadata.

process METADATA {
    label 'process_low'
    tag "metadata"

    output:
    val true, emit: done

    script:
    """
    repgenr ${params.repgenr_opts} metadata -wd ${params.workdir} ${params.metadata_args}
    """
}
