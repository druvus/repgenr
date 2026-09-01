// Collect every process's versions.yml fragment into one published file.
//
// Fragments are de-duplicated and concatenated into
// <outdir>/pipeline_info/software_versions.yml so a run records which tool
// versions produced it.

workflow PUBLISH_VERSIONS {
    take:
    ch_versions   // channel: versions.yml fragments from every process

    main:
    ch_versions
        .unique()
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
        )
}
