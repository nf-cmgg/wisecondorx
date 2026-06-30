#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/wisecondorx
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/wisecondorx
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WISECONDORX  } from './workflows/wisecondorx'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_wisecondorx_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_wisecondorx_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_wisecondorx_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {

    // Path to comma-separated file containing information about the samples in the experiment.
    input: Path

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: String

    // The prefix to be used for the output file. If this parameter isn't given, the following format will be used: WisecondorX_DDMMYYYY
    prefix: String?

    // Email address for completion summary.
    email: String?

    // MultiQC report title. Printed as page header, used for filename if not otherwise specified.
    multiqc_title: String?

    // A comma-delimited list of bin sizes in kilobases to use for the analysis
    bin_sizes: String = '500,100,50,15,10,5'

    // Name of iGenomes reference.
    genome: String

    // Path to FASTA genome file.
    fasta: String = getGenomeAttribute('fasta')

    // Path to FASTA genome index file.
    fai: String = getGenomeAttribute('fai')

    // The path to the directory containing the files specified by the CMGG references config.
    genomes_base: String?

    // Do not load CMGG references
    genomes_ignore: Boolean

    // Directory / URL for CMGG configs
    cmgg_config_base: String = '/conf/'

    // Do not load the iGenomes reference config.
    igenomes_ignore: Boolean

    // The base path to the igenomes reference files
    igenomes_base: String?

    // Git commit id for Institutional configs.
    custom_config_version: String = 'master'

    // Base directory for Institutional configs.
    custom_config_base: String = 'https://raw.githubusercontent.com/nf-core/configs/master'

    // Institutional config name.
    config_profile_name: String?

    // Institutional config description.
    config_profile_description: String?

    // Institutional config contact information.
    config_profile_contact: String?

    // Institutional config URL link.
    config_profile_url: String?

    // Display version and exit.
    version: Boolean

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String = 'copy'

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean

    // File size limit when attaching MultiQC reports to summary emails.
    max_multiqc_email_size: String = '25.MB'

    // Do not use coloured log outputs.
    monochrome_logs: Boolean

    // Custom config file to supply to MultiQC.
    multiqc_config: Path?

    // Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file
    multiqc_logo: Path?

    // Custom MultiQC yaml file containing HTML including a methods description.
    multiqc_methods_description: Path?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: String = 'https://github.com/nf-cmgg/test-datasets/raw/refs/heads/wisecondorx/'

    // Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss.
    trace_report_suffix: String?

    // Display the help message.
    help = false

    // Display the full detailed help message.
    help_full: Boolean

    // Display hidden parameters in the help message (only works when --help or --help_full are provided).
    show_hidden: Boolean
}

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    WISECONDORX (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.fasta,
        params.fai,
        params.bin_sizes in String ? params.bin_sizes.tokenize(",") : [params.bin_sizes as String],
        params.prefix,
        params.outdir,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        WISECONDORX.out.multiqc.map { rec -> rec.report }
    )

    publish:
    multiqc        = WISECONDORX.out.multiqc
    references     = WISECONDORX.out.references
    npz            = WISECONDORX.out.npz
}

output {
    multiqc        { path "multiqc/" }
    references     { path { rec ->
        rec.npz >> "${rec.id}_${rec.bin_size}kbp.npz"
    } }
    npz            { path { rec ->
        rec.npz >> "npz/${rec.id}.npz"
    } }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
