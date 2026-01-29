nextflow.preview.types = true
process WISECONDORX_CONVERT {
    tag "$input.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0':
        'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0' }"

    input:
    input: WisecondorxConvertInput
    (_id2: String, fasta: Path, _fai: Path): Record

    output:
    record(id: input.id, npz: file("*.npz"))

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple("${task.process}", "wisecondorx", '1.2.9') >> 'versions'

    script:
    def args = input.args ?: ''
    def prefix = input.prefix ?: "${input.id}"
    def reference = fasta ? "--reference ${fasta}" : ""

    """
    WisecondorX convert \\
        ${input.bam} \\
        ${prefix}.npz \\
        ${reference} \\
        ${args}
    """

    stub:
    def prefix = input.prefix ?: "${input.id}"

    """
    touch ${prefix}.npz
    """
}

record WisecondorxConvertInput {
    id: String
    bam: Path
    bai: Path?
    args: String?
    prefix: String?
}