nextflow.preview.types = true
process WISECONDORX_CONVERT {
    tag "$id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0':
        'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0' }"

    input:
    (id: String, bam: Path, _bai: Path): Record
    (_id2: String, fasta: Path, _fai: Path): Record

    output:
    record(id: id, npz: file("*.npz"))

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple("${task.process}", "wisecondorx", '1.2.9') >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"
    def reference = fasta ? "--reference ${fasta}" : ""

    """
    WisecondorX convert \\
        ${bam} \\
        ${prefix}.npz \\
        ${reference} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${id}"

    """
    touch ${prefix}.npz
    """
}
