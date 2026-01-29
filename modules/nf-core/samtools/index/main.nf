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
        $input.bam
    """

    stub:
    """
    touch ${input.bam}.bai
    touch ${input.bam}.crai
    touch ${input.bam}.csi
    """
}

record SamtoolsIndexInput {
    id: String
    bam: Path
    args: String?
}
