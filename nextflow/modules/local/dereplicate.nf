// Cluster genomes by ANI and select representatives.

process DEREPLICATE {
    label 'process_high'
    tag "dereplicate"

    input:
    val ready

    output:
    val true, emit: done

    script:
    // Match the tool's thread count to the CPUs the scheduler reserved for this
    // task, and expose the two-stage chunking knobs as first-class params. A -t
    // in dereplicate_args (placed last) still wins, as an explicit override.
    def chunk = params.derep_process_size ? "--process-size ${params.derep_process_size}" : ''
    def nproc = params.derep_num_processes ? "--num-processes ${params.derep_num_processes}" : ''
    """
    repgenr ${params.repgenr_opts} dereplicate -wd ${params.workdir} -t ${task.cpus} ${chunk} ${nproc} ${params.dereplicate_args}
    """
}
