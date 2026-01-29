nextflow.preview.types = true
process SAMTOOLS_INDEX {
    tag "$input.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    input: SamtoolsIndexInput

    output:
    record(id: input.id, bai: file("*.bai"), crai: file("*.crai"), csi: file("*.csi"))

    topic:
    tuple("${task.process}", "samtools", eval('samtools --version | head -1 | sed -e "s/samtools //"')) >> 'versions'

    script:
    def args = input.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        ${input.bam.name}
    """

    stub:
    """
    touch ${input.bam.name}.bai
    touch ${input.bam.name}.crai
    touch ${input.bam.name}.csi
    """
}

record SamtoolsIndexInput {
    id: String
    bam: Path
    args: String?
}
