nextflow.preview.types = true
process WISECONDORX_NEWREF {
    tag "$input.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0':
        'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0' }"

    input:
    input: WisecondorxNewrefInput

    output:
    record(id: input.id, npz: file("*.npz"))

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple("${task.process}", "wisecondorx", '1.2.9') >> 'versions'

    script:
    def args = input.args ?: ''
    def prefix = input.prefix ?: "${input.id}"

    input.npzs.each { npz -> 
        if("${npz}" == "${prefix}.npz") error "${npz} has the same name as the output file, set prefix in module configuration to disambiguate!"
    }

    """
    WisecondorX \\
        newref \\
        *.npz \\
        ${prefix}.npz \\
        --cpus ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = input.prefix ?: "${input.id}"

    input.npzs.each { npz -> 
        if("${npz}" == "${prefix}.npz") error "${npz} has the same name as the output file, set prefix in module configuration to disambiguate!"
    }

    """
    touch ${prefix}.npz
    """
}

record WisecondorxNewrefInput {
    id: String
    npzs: List<Path>
    args: String?
    prefix: String?
}