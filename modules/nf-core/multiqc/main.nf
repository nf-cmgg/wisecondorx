nextflow.enable.types = true

process MULTIQC {
    tag "${input.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    input: MultiqcInput

    stage:
    stageAs input.multiqc_files, "?/*"
    stageAs input.multiqc_config, "?/*"

    output:
    input + record(
        report: file("*.html"),
        data: file("*_data"),
        plots: file("*_plots")
    )

    topic:
    tuple("${task.process}", 'multiqc', eval('multiqc --version | sed "s/.* //g"')) >> 'versions_multiqc'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = input.multiqc_config ? input.multiqc_config instanceof List ? "--config ${input.multiqc_config.join(' --config ')}" : "--config ${input.multiqc_config}" : ""
    def logo = input.multiqc_logo ? "--cl-config 'custom_logo: \"${input.multiqc_logo}\"'" : ''
    def replace = input.replace_names ? "--replace-names ${input.replace_names}" : ''
    def samples = input.sample_names ? "--sample-names ${input.sample_names}" : ''
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_plots/.stub
    touch multiqc_report.html
    """
}

record MultiqcInput {
    id: String
    multiqc_files: Set<Path>
    multiqc_config: Set<Path>
    multiqc_logo: Path?
    replace_names: Path?
    sample_names: Path?
}