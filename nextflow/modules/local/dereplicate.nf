// Cluster genomes by ANI and select representatives.

process DEREPLICATE {
    label 'process_high'
    tag "dereplicate"

    input:
    val ready

    output:
    val true, emit: done

    script:
    """
    repgenr ${params.repgenr_opts} dereplicate -wd ${params.workdir} ${params.dereplicate_args}
    """
}
