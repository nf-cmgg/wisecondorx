nextflow.enable.types = true

process SAMTOOLS_FAIDX {
    tag "${input.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    input: SamtoolsFaidxInput

    output:
    record(
        id: input.id,
        fai: file("*.fai", optional:true),
        fasta: file("*.{fa,fasta}", optional:true),
        sizes: file("*.sizes", optional:true),
        gzi: file("*.gzi", optional:true)
    )

    topic:
    tuple("${task.process}", 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def get_sizes_command = input.get_sizes ? "cut -f 1,2 ${input.fasta}.fai > ${input.fasta}.sizes" : ''
    """
    samtools \\
        faidx \\
        ${input.fasta} \\
        ${args}

    ${get_sizes_command}
    """

    stub:
    def match = (task.ext.args =~ /-o(?:utput)?\s(.*)\s?/).findAll()
    def fastacmd = match[0] ? "touch ${match[0][1]}" : ''
    def get_sizes_command = input.get_sizes ? "touch ${input.fasta}.sizes" : ''
    """
    ${fastacmd}
    touch ${input.fasta}.fai
    if [[ "${input.fasta.extension}" == "gz" ]]; then
        touch ${input.fasta}.gzi
    fi

    ${get_sizes_command}
    """
}

record SamtoolsFaidxInput {
    id: String
    fasta: Path
    fai: Path?
    get_sizes: Boolean
}
