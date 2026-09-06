// Build the phylogeny (data-channel form).
//
// Calls the stateless `repgenr phylo-build` step directly on the staged channel
// files: the merged representatives directory provides the genome set, the
// outgroup files are staged into an outgroup/ directory, and tree/tree.nwk is
// emitted as a channel output. Tool flags arrive as task.ext.args from
// conf/modules.config; publishing is configured there too.

process PHYLO {
    tag "phylo"
    // debug forwards task stdout into Nextflow's own log, which is what
    // lets an nf-test assertion on process.stdout see the stub's ext.args
    // echo (the stub's .command.out is otherwise not exposed there).
    debug true
    label 'process_high'

    input:
    path reps_dir
    path outgroup, stageAs: 'outgroup/*'
    path outgroup_accession

    output:
    path "tree/tree.nwk", emit: tree
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    # Forward tool exit codes (OOM kill -> 137) so errorStrategy can retry.
    export REPGENR_PROPAGATE_TOOL_EXIT=1

    repgenr ${opts} phylo-build \\
        --genomes-dir ${reps_dir}/representatives \\
        --outgroup-dir outgroup \\
        --outgroup-accession ${outgroup_accession} \\
        -o . -t ${task.cpus} ${args} \\
        --versions-out tool_versions.yml

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    repgenr: \$(repgenr --version | sed 's/repgenr //')
END_VERSIONS
    cat tool_versions.yml >> versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo "ext.args: ${args}"
    mkdir -p tree
    names=\$(ls ${reps_dir}/representatives | sed 's/\\.[^.]*\$//' | paste -sd, -)
    echo "(\${names});" > tree/tree.nwk
    touch versions.yml
    """
}
