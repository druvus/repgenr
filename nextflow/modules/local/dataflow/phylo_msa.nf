// Build the alignment a tree builder consumes (data-channel form).
//
// The first half of `repgenr phylo-build`, run as its own task so that the
// alignment is not repeated when only the tree builder or the bootstrap
// changes: those belong to PHYLO_TREE, which caches separately. Enabled by
// params.phylo_split_msa; otherwise PHYLO does both halves in one task.

process PHYLO_MSA {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(reps_dir), path(outgroup, stageAs: 'outgroup/*'), path(outgroup_accession)

    output:
    tuple val(meta), path("msa.fasta"), emit: msa
    path "versions.yml"               , emit: versions

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
        --msa-only \\
        --versions-out tool_versions.yml

    repgenr_versions_fragment "${task.process}" tool_versions.yml
    """

    stub:
    def args = task.ext.args ?: ''
    def opts = task.ext.repgenr_opts ?: ''
    """
    echo "ext.args: ${args}"
    echo "ext.repgenr_opts: ${opts}"
    for genome in ${reps_dir}/representatives/*; do
        name=\$(basename "\$genome" | sed 's/\\.[^.]*\$//')
        echo ">\${name}"
        echo "ACGT"
    done > msa.fasta
    touch versions.yml
    """
}
