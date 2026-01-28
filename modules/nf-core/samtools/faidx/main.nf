nextflow.preview.types = true
process SAMTOOLS_FAIDX {
    tag "$id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    (id: String, fasta: Path): Record

    output:
    record(id: id, fasta: file("*.{fa,fasta}"), fai: file("*.fai"), gzi: file("*.gzi"))
    
    topic:
    tuple("${task.process}", "samtools", eval('samtools --version | head -1 | sed -e "s/samtools //"')) >> 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        faidx \\
        $fasta \\
        $args
    """

    stub:
    def match = (task.ext.args =~ /-o(?:utput)?\s(.*)\s?/).findAll()
    def fastacmd = match[0] ? "touch ${match[0][1]}" : ''
    """
    ${fastacmd}
    touch ${fasta}.fai
    """
}
