/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { fromSamplesheet; paramsSummaryMap } from 'plugin/nf-validation'

def summary_params = paramsSummaryMap(workflow)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkPathParamList = [
    params.input,
    params.multiqc_config,
    params.fasta,
    params.fai
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'
include { NGSBITS_SAMPLEGENDER        } from '../modules/nf-core/ngsbits/samplegender/main'
include { WISECONDORX_CONVERT         } from '../modules/nf-core/wisecondorx/convert/main'
include { WISECONDORX_NEWREF          } from '../modules/nf-core/wisecondorx/newref/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CMGGWISECONDORX {

    ch_versions = Channel.empty()

    //
    // Create optional input files
    //

    ch_fasta = Channel.fromPath(params.fasta, checkIfExists:true)
        .map { [[id:params.genome ?: "fasta"], it] }
        .collect()

    if(!params.fai) {
        SAMTOOLS_FAIDX(ch_fasta)
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

        SAMTOOLS_FAIDX.out.fai
            .set { ch_fai }
    } else {
        ch_fai = Channel.fromPath(params.fai, checkIfExists:true)
            .map { [[id:params.genome ?: "fasta"], it] }
            .collect()
    }

    //
    // Validate and convert the samplesheet
    //

    Channel.fromSamplesheet("input", immutable_meta: false)
        .branch { meta, cram, crai ->
            new_meta = meta + [id:cram.baseName]
            indexed: crai
                return [ new_meta, cram, crai ]
            not_indexed: !crai
                return [ new_meta, cram ]
        }
        .set { ch_input }

    //
    // Index the non-indexed input files
    //

    SAMTOOLS_INDEX(ch_input.not_indexed)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_input.not_indexed
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_input.indexed)
        .set { ch_indexed }

    if(!params.no_metrics){

        //
        // Define the gender if it's not given
        //

        ch_indexed
            .branch { meta, cram, crai ->
                gender: meta.gender
                    [ meta, meta.gender ]
                no_gender: !meta.gender
            }
            .set { ch_ngsbits_input }

        NGSBITS_SAMPLEGENDER(
            ch_ngsbits_input.no_gender,
            ch_fasta,
            ch_fai,
            'xy'
        )
        ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions.first())

        NGSBITS_SAMPLEGENDER.out.tsv
            .map { meta, tsv ->
                gender = get_gender(tsv)
                new_meta = meta + [gender: gender]
                [ new_meta, gender ]
            }
            .mix(ch_ngsbits_input.gender)
            .set { ch_genders }

        //
        // Create a small metrics file
        //

        ch_genders
            .reduce([:]) { counts, entry ->
                meta = entry[0]
                gender = entry[1]
                counts[gender] = (counts[gender] ?: []) + meta.id
                counts
            }
            .map { genders -> create_metrics(genders)}

    }

    //
    // Convert the input files to NPZ files
    //

    WISECONDORX_CONVERT(
        ch_indexed,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(WISECONDORX_CONVERT.out.versions.first())

    //
    // Create the WisecondorX reference
    //

    // Define reference name (with timestamp) => only used when --prefix is null
    Date date = new Date()
    String dateFormat = "WisecondorX_${date.format("ddMMyyyy")}"

    WISECONDORX_CONVERT.out.npz
        .map { meta, npz ->
            new_meta = [id:params.prefix ?: dateFormat]
            [ new_meta, npz ]
        }
        .groupTuple() // All files should be present here, so no size is needed
        .combine(params.bin_sizes.split(","))
        .map { meta, npz, bin_size ->
            new_meta = meta + [bin_size:bin_size]
            [ new_meta, npz ]
        }
        .set { ch_newref_input }

    WISECONDORX_NEWREF(ch_newref_input)
    ch_versions = ch_versions.mix(WISECONDORX_NEWREF.out.versions.first())

    //
    // Dump software versions
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCmggWisecondorx.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCmggWisecondorx.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def get_gender(tsv) {
    split_tsv = tsv.splitCsv(sep:"\t", header:true, strip:true)
    return split_tsv[0].gender
}

def create_metrics(genders) {
    male = genders["male"]
    female = genders["female"]
    male_count = male.size()
    female_count = female.size()
    total_count = male_count + female_count
    ratio_male_to_female = male_count / female_count

    output = file("${params.outdir}/metrics.txt")
    output.write("Ratio male to female: ${ratio_male_to_female}\n")
    output.append("Male count: ${male_count}\n")
    output.append("Female count: ${female_count}\n")
    output.append("Total count: ${total_count}\n")
    output.append("Males: ${male.join(',')}\n")
    output.append("Females: ${female.join(',')}\n")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
