nextflow.preview.types = true
process NGSBITS_SAMPLEGENDER {
    tag "$id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2be56a07ac1d5a447a10fd061be4d6144620bec00bac834f58c2bdef0330147f/data':
        'community.wave.seqera.io/library/ngs-bits:2025_09--f6ea3a4494373ed6' }"

    input:
    (id: String, bam: Path, _bai: Path, method: String): Record
    (_id2: String, fasta: Path, _fai: Path): Record

    output:
    record(id: id, tsv: file("*.tsv"))

    topic:
    tuple("${task.process}", "ngs-bits", eval("$(SampleGender --version 2>&1 | sed 's/SampleGender //')")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    def ref = fasta ? "-ref ${fasta}" : ""
    """
    SampleGender \\
        -in ${bam} \\
        -method ${method} \\
        -out ${prefix}.tsv \\
        ${ref} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"
    """
    echo "#file	gender	reads_chry	reads_chrx	ratio_chry_chrx" > ${prefix}.tsv
    echo "${id}	female	48	12423	0.0039" >> ${prefix}.tsv
    """
}
