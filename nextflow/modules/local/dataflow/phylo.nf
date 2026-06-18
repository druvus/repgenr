// Build the phylogeny (data-channel form).
//
// The phylo stage drives an aligner/SNP-typer/tree-builder matrix that is tied to
// the RepGenR working-directory layout, so this module bridges: it stages the
// channel inputs (merged representatives, outgroup, selection) into a task-local
// workdir, runs the existing `repgenr phylo`, and emits tree.nwk as a channel
// output. Data flows through channels; the workdir is task scratch.

process PHYLO {
    label 'process_high'
    tag "phylo"
    publishDir "${params.outdir}/phylo", mode: 'copy'

    input:
    path reps_dir
    path outgroup
    path selection

    output:
    path "tree/tree.nwk", emit: tree
    path "versions.yml" , emit: versions

    script:
    """
    mkdir -p wd/derep/representatives wd/outgroup
    cp ${reps_dir}/representatives/* wd/derep/representatives/ 2>/dev/null || true
    [ -f ${reps_dir}/clusters.tsv ] && cp ${reps_dir}/clusters.tsv wd/derep/clusters.tsv

    for f in ${outgroup}; do
        [ -e "\$f" ] && cp "\$f" wd/outgroup/
    done
    awk -F'\\t' 'NR>1 && \$5==1 {print \$1}' ${selection} | head -1 > wd/outgroup_accession.txt

    repgenr ${params.repgenr_opts} phylo -wd wd -t ${task.cpus} ${params.phylo_args}
    mkdir -p tree
    cp wd/tree/tree.nwk tree/tree.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repgenr: \$(repgenr --version | sed 's/repgenr //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p tree
    names=\$(ls ${reps_dir}/representatives | sed 's/\\.[^.]*\$//' | paste -sd, -)
    echo "(\${names});" > tree/tree.nwk
    touch versions.yml
    """
}
