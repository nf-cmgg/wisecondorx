nextflow.enable.types = true

process SAMTOOLS_INDEX {
    tag "${input.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    input: SamtoolsIndexInput

    output:
    record(
        id: input.id,
        input_idx: file("*.{bai,csi,crai}")
    )

    topic:
    tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${args} \\
        ${input.input}
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = input.input.getExtension() == 'cram'
        ? "crai"
        : args.contains("-c") ? "csi" : "bai"
    """
    touch ${input.input}.${extension}
    """
}

record SamtoolsIndexInput {
    id: String
    input: Path
}
