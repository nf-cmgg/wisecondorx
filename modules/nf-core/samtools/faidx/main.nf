nextflow.preview.types = true
process SAMTOOLS_FAIDX {
    tag "$input.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    input: SamtoolsFaidxInput

    output:
    record(id: input.id, fasta: file("*.{fa,fasta}"), fai: file("*.fai"), gzi: file("*.gzi"))
    
    topic:
    tuple("${task.process}", "samtools", eval('samtools --version | head -1 | sed -e "s/samtools //"')) >> 'versions'

    script:
    def args = input.args ?: ''
    """
    samtools \\
        faidx \\
        $input.fasta \\
        $args
    """

    stub:
    def args = input.args ?: ''
    def match = (args =~ /-o(?:utput)?\s(.*)\s?/).findAll()
    def fastacmd = match[0] ? "touch ${match[0][1]}" : ''
    """
    ${fastacmd}
    touch ${input.fasta}.fai
    """
}

record SamtoolsFaidxInput {
    id: String
    fasta: Path
    fai: Path?
    args: String?
}
