// Select a taxon's genomes (data-channel form).
//
// Runs the metadata stage in a task-local working directory and emits the
// portable selection.tsv (accession + taxonomy + outgroup flag + filename) as a
// channel output -- the hand-off the genome step consumes. No shared workdir.

process METADATA {
    label 'process_low'
    tag "metadata"
    publishDir "${params.outdir}/metadata", mode: 'copy'

    output:
    path "selection.tsv"        , emit: selection
    path "outgroup_accession.txt", emit: outgroup_accession
    path "versions.yml"         , emit: versions

    script:
    """
    repgenr ${params.repgenr_opts} metadata -wd metadata_wd ${params.metadata_args}
    cp metadata_wd/selection.tsv selection.tsv
    cp metadata_wd/outgroup_accession.txt outgroup_accession.txt

    repgenr versions -wd metadata_wd --versions-out tool_versions.yml
    cat > versions.yml <<END_VERSIONS
"${task.process}":
    repgenr: \$(repgenr --version | sed 's/repgenr //')
END_VERSIONS
    cat tool_versions.yml >> versions.yml
    """

    stub:
    """
    printf 'accession\\tfamily\\tgenus\\tspecies\\tis_outgroup\\tfilename\\n' > selection.tsv
    printf 'GCF_000001.1\\tFam\\tGen\\tsp1\\t0\\tFam_Gen_sp1_GCF_000001.1.fasta\\n' >> selection.tsv
    printf 'GCF_000002.1\\tFam\\tGen\\tsp2\\t0\\tFam_Gen_sp2_GCF_000002.1.fasta\\n' >> selection.tsv
    printf 'GCF_000009.1\\tFam\\tOut\\tgrp\\t1\\tFam_Out_grp_GCF_000009.1.fasta\\n' >> selection.tsv
    printf 'GCF_000009.1\\n' > outgroup_accession.txt
    touch versions.yml
    """
}
