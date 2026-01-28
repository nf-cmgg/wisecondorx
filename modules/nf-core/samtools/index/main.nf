nextflow.preview.types = true
process SAMTOOLS_INDEX {
    tag "$id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    (id: String, bam: Path): Record

    output:
    record(id: id, bai: file("*.bai"), crai: file("*.crai"), csi: file("*.csi"))

    topic:
    tuple("${task.process}", "samtools", eval('samtools --version | head -1 | sed -e "s/samtools //"')) >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        $bam
    """

    stub:
    """
    touch ${bam}.bai
    touch ${bam}.crai
    touch ${bam}.csi
    """
}
