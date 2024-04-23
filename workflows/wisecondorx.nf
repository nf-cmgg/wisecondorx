/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'
include { NGSBITS_SAMPLEGENDER        } from '../modules/nf-core/ngsbits/samplegender/main'
include { WISECONDORX_CONVERT         } from '../modules/nf-core/wisecondorx/convert/main'
include { WISECONDORX_NEWREF          } from '../modules/nf-core/wisecondorx/newref/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_wisecondorx_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WISECONDORX {

    take:
    ch_samplesheet              // queue channel:   samplesheet read in from --input
    fasta                       // string:          the reference fasta file
    fai                         // string:          the index of the reference fasta file
    val_bin_sizes               // list:            a list of bin sizes to use
    no_metrics                  // boolean:         deactivate the generation of metrics
    prefix                      // string:          the prefix to be used by the output file
    outdir                      // string:          the path of the output directory
    multiqc_config              // string:          the path to the multiqc config
    multiqc_logo                // string:          the path to the multiqc logo
    multiqc_methods_description // file:            the file containing the multiqc custom method descriptions

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Create optional input files
    //

    ch_fasta = Channel.fromPath(fasta, checkIfExists:true)
        .map { [[id:"fasta"], it] }
        .collect()

    if(!fai) {
        SAMTOOLS_FAIDX(
            ch_fasta,
            [[],[]]
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

        SAMTOOLS_FAIDX.out.fai
            .set { ch_fai }
    } else {
        ch_fai = Channel.fromPath(fai, checkIfExists:true)
            .map { [[id:"fasta"], it] }
            .collect()
    }

    ch_samplesheet
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

    if(!no_metrics){

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
            .collectFile(name: "metrics_mqc.tsv")
            .set { ch_metrics }

        ch_multiqc_files = ch_multiqc_files.mix(ch_metrics)

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
            new_meta = [id:prefix ?: dateFormat]
            [ new_meta, npz ]
        }
        .groupTuple() // All files should be present here, so no size is needed
        .combine(val_bin_sizes)
        .map { meta, npz, bin_size ->
            new_meta = meta + [bin_size:bin_size]
            [ new_meta, npz ]
        }
        .set { ch_newref_input }

    WISECONDORX_NEWREF(ch_newref_input)
    ch_versions = ch_versions.mix(WISECONDORX_NEWREF.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = multiqc_config ? Channel.fromPath(multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = multiqc_logo ? Channel.fromPath(multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = multiqc_methods_description
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
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
    def List male = genders["male"]
    def List female = genders["female"]
    def Integer male_count = male.size()
    def Integer female_count = female.size()
    def Float male_to_female_ratio = male_count / female_count
    def Integer total_count = male_count + female_count

    return """# plot_type: 'table'
Male to female ratio\tMale count\tFemale count\tTotal count\tMales\tFemales
${male_to_female_ratio}\t${male_count}\t${female_count}\t${total_count}\t${male.join(",")}\t${female.join(",")}
"""
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
