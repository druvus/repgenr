// Select a taxon's genomes from GTDB metadata.

process METADATA {
    label 'process_low'
    tag "metadata"

    output:
    val true, emit: done

    script:
    """
    repgenr metadata -wd ${params.workdir} ${params.metadata_args}
    """
}
