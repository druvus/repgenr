// Acquire viral genomes (data-channel form): vmetadata + vgenome.
//
// The viral selection/download writes genomes directly (BV-BRC FASTA, length
// filtering, mashtree outgroup), so this module bridges: it runs vmetadata +
// vgenome in a task-local workdir and emits the genomes, the outgroup, and the
// outgroup accession as channel outputs that feed the same downstream data-channel
// dereplication/phylo/tree2tax as the bacterial path.

process VACQUIRE {
    label 'process_medium'
    tag "vacquire"
    publishDir "${params.outdir}/genomes", mode: 'copy'

    output:
    path "out/genomes/*"         , emit: genomes
    path "out/outgroup/*"        , emit: outgroup, optional: true
    path "outgroup_accession.txt", emit: outgroup_accession
    path "versions.yml"          , emit: versions

    script:
    """
    repgenr ${params.repgenr_opts} vmetadata -wd wd ${params.vmetadata_args}
    repgenr ${params.repgenr_opts} vgenome   -wd wd ${params.vgenome_args}

    mkdir -p out
    cp -r wd/genomes out/genomes
    [ -d wd/outgroup ] && cp -r wd/outgroup out/outgroup || true
    if [ -f wd/outgroup_accession.txt ]; then
        cp wd/outgroup_accession.txt outgroup_accession.txt
    else
        : > outgroup_accession.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        repgenr: \$(repgenr --version | sed 's/repgenr //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p out/genomes out/outgroup
    printf '>x\\nACGT\\n' > out/genomes/Vir_gen_sp1_iso1.fasta
    printf '>x\\nACGT\\n' > out/genomes/Vir_gen_sp2_iso2.fasta
    printf '>x\\nACGT\\n' > out/outgroup/Vir_out_grp_iso9.fasta
    printf 'iso9\\n' > outgroup_accession.txt
    touch versions.yml
    """
}
