nextflow.preview.types = true
process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c6c120d559d7ee04c7442b61ad7cf5a9e8970be5feefb37d68eeaa60c1034eb/data' :
        'community.wave.seqera.io/library/multiqc:1.32--d58f60e4deb769bf' }"

    input:
    input: MultiqcInput

    stage:
    stageAs "?/*", input.multiqc_files

    output:
    record(report: file("*.html"), data: file("*_data"), plots: file("*_plots"))
    // versions_multiqc: Tuple<String> = tuple("${task.process}", "multiqc", eval('multiqc --version | sed -e "s/multiqc, version //g"'))

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = input.multiqc_config ? "--config ${input.multiqc_config}" : ''
    def extra_config = input.extra_multiqc_config ? "--config ${input.extra_multiqc_config}" : ''
    def logo = input.multiqc_logo ? "--cl-config 'custom_logo: \"${input.multiqc_logo}\"'" : ''
    def replace = input.replace_names ? "--replace-names ${input.replace_names}" : ''
    def samples = input.sample_names ? "--sample-names ${input.sample_names}" : ''
    """
    multiqc \\
        --force \\
        $args \\
        $config \\
        $prefix \\
        $extra_config \\
        $logo \\
        $replace \\
        $samples \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html
    """
}

record MultiqcInput {
    multiqc_files: List<Path>
    multiqc_config: Path?
    extra_multiqc_config: Path?
    multiqc_logo: Path?
    replace_names: Path?
    sample_names: Path?
}
