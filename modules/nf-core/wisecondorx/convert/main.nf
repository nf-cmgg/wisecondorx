nextflow.enable.types = true

process WISECONDORX_CONVERT {
    tag "${input.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/13af39819608398807612090d4b8af7dedb8db403967e71af22dbbeeb502ead1/data'
        : 'community.wave.seqera.io/library/wisecondorx:1.3.0--835c946afbce9082'}"

    input:
    input: WisecondorxConvertInput

    output:
    record(
        id: input.id,
        npz: file("*.npz")
    )

    topic:
    tuple("${task.process}", 'wisecondorx', eval("pip list |& sed -n 's/wisecondorx *//p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${input.id}"
    def reference = input.fasta ? "--reference ${input.fasta}" : ""

    """
    WisecondorX convert \\
        ${input.input} \\
        ${prefix}.npz \\
        ${reference} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${input.id}"

    """
    touch ${prefix}.npz
    """
}

record WisecondorxConvertInput {
    id: String
    input: Path
    input_idx: Path
    fasta: Path
    fai: Path
}
