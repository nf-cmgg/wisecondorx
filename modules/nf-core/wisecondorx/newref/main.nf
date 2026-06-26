nextflow.enable.types = true

process WISECONDORX_NEWREF {
    tag "${input.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/13af39819608398807612090d4b8af7dedb8db403967e71af22dbbeeb502ead1/data'
        : 'community.wave.seqera.io/library/wisecondorx:1.3.0--835c946afbce9082'}"

    input:
    input: WisecondorxNewrefInput

    output:
    input + record(npz: file("*.npz"))

    topic:
    tuple("${task.process}", 'wisecondorx', eval("pip list |& sed -n 's/wisecondorx *//p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${input.id}"

    input.inputs.each { input_file ->
        if ("${input_file}" == "${prefix}.npz") {
            error("${input_file} has the same name as the output file, set prefix in module configuration to disambiguate!")
        }
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
    def prefix = task.ext.prefix ?: "${input.id}"

    input.inputs.each { input_file ->
        if ("${input_file}" == "${prefix}.npz") {
            error("${input_file} has the same name as the output file, set prefix in module configuration to disambiguate!")
        }
    }

    """
    touch ${prefix}.npz
    """
}

record WisecondorxNewrefInput {
    id: String
    inputs: List<Path>
    bin_size: Integer
}