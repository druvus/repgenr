// Build the phylogeny (data-channel form).
//
// Calls the stateless `repgenr phylo-build` step directly on the staged channel
// files: the merged representatives directory provides the genome set, the
// outgroup files are staged into an outgroup/ directory, and tree/tree.nwk is
// emitted as a channel output. There is no shared working directory; align/SNP
// intermediates stay in task scratch and are not published.

process PHYLO {
    label 'process_high'
    tag "phylo"
    publishDir "${params.outdir}/phylo", mode: 'copy'

    input:
    path reps_dir
    path outgroup, stageAs: 'outgroup/*'
    path outgroup_accession

    output:
    path "tree/tree.nwk", emit: tree
    path "versions.yml" , emit: versions

    script:
    """
    repgenr ${params.repgenr_opts} phylo-build \\
        --genomes-dir ${reps_dir}/representatives \\
        --outgroup-dir outgroup \\
        --outgroup-accession ${outgroup_accession} \\
        -o . -t ${task.cpus} ${params.phylo_args} \\
        --versions-out tool_versions.yml

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    repgenr: \$(repgenr --version | sed 's/repgenr //')
END_VERSIONS
    cat tool_versions.yml >> versions.yml
    """

    stub:
    """
    mkdir -p tree
    names=\$(ls ${reps_dir}/representatives | sed 's/\\.[^.]*\$//' | paste -sd, -)
    echo "(\${names});" > tree/tree.nwk
    touch versions.yml
    """
}
