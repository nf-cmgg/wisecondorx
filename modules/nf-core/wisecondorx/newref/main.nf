nextflow.preview.types = true
process WISECONDORX_NEWREF {
    tag "$id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wisecondorx:1.2.9--pyhdfd78af_0':
        'biocontainers/wisecondorx:1.2.9--pyhdfd78af_0' }"

    input:
    (id: String, npzs: List<Path> ) : Record

    output:
    record(id: id, npz: file("*.npz"))

    topic:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple("${task.process}", "wisecondorx", '1.2.9') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"

    npzs.each { input -> 
        if("${input}" == "${prefix}.npz") error "${input} has the same name as the output file, set prefix in module configuration to disambiguate!"
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
    def prefix = task.ext.prefix ?: "${id}"

    npzs.each { input -> 
        if("${input}" == "${prefix}.npz") error "${input} has the same name as the output file, set prefix in module configuration to disambiguate!"
    }

    """
    touch ${prefix}.npz
    """
}
