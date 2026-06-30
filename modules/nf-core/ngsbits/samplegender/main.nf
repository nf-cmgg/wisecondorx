nextflow.enable.types = true

process NGSBITS_SAMPLEGENDER {
    tag "${input.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fb/fbf8cfd89c36e9a18a895066bb1da04b93ef585a593b0821ec7037aba6c03474/data'
        : 'community.wave.seqera.io/library/ngs-bits:2025_12--958625b0e620100a'}"

    input:
    input: NgsbitsSampleGenderInput

    output:
    record(
        id: input.id,
        tsv: file("*.tsv")
    )

    topic:
    tuple("${task.process}", 'ngsbits', eval("SampleGender --version  2>&1 | sed 's/SampleGender //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${input.id}"
    def ref = input.fasta ? "-ref ${input.fasta}" : ""
    """
    SampleGender \\
        -in ${input.input} \\
        -method ${input.method} \\
        -out ${prefix}.tsv \\
        ${ref} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${input.id}"
    """
    echo "#file	gender	reads_chry	reads_chrx	ratio_chry_chrx" > ${prefix}.tsv
    echo "${input.id}	female	48	12423	0.0039" >> ${prefix}.tsv
    """
}

record NgsbitsSampleGenderInput {
    id: String
    input: Path
    input_idx: Path
    fasta: Path
    fai: Path
    method: String // 'xy', 'hetx' or 'sry'
}
