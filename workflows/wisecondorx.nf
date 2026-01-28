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
    ch_samplesheet: Channel<Path>       // samplesheet read in from --input
    fasta: String                       // the reference fasta file
    fai: String                         // the index of the reference fasta file
    val_bin_sizes: List<Integer>        // a list of bin sizes to use
    prefix: String                      // the prefix to be used by the output file
    outdir: String                      // the path of the output directory
    multiqc_config: String              // the path to the multiqc config
    multiqc_logo: String                // the path to the multiqc logo
    multiqc_methods_description: Path   // the file containing the multiqc custom method descriptions

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Create optional input files
    //

    def ch_fasta = channel.fromPath(fasta, checkIfExists:true)
        .map { fasta_file -> tuple([id:"fasta"], fasta_file) }
        .collect()

    def ch_fai
    if(!fai) {
        SAMTOOLS_FAIDX(
            ch_fasta,
            [[],[]]
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

        ch_fai = SAMTOOLS_FAIDX.out.fai
    } else {
        ch_fai = channel.fromPath(fai, checkIfExists:true)
            .map { fai_file -> tuple([id:"fai"], fai_file) }
            .collect()
    }

    def ch_input = ch_samplesheet
        .branch { meta, cram, crai, npz ->
            npz: npz
                return [ meta, npz ]
            indexed: crai
                return [ meta, cram, crai ]
            not_indexed: !crai
                return [ meta, cram ]
        }

    //
    // Index the non-indexed input files
    //

    SAMTOOLS_INDEX(ch_input.not_indexed)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def ch_indexed = ch_input.not_indexed
        .join(SAMTOOLS_INDEX.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .mix(ch_input.indexed)

    //
    // Define the sex if it's not given
    //

    ch_indexed
        .branch { meta, _cram, _crai ->
            sex: meta.sex
                [ meta, meta.sex ]
            no_sex: !meta.sex
        }
        .set { ch_ngsbits_input }

    NGSBITS_SAMPLEGENDER(
        ch_ngsbits_input.no_sex,
        ch_fasta,
        ch_fai,
        'xy'
    )
    ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions.first())

    def ch_sexes = NGSBITS_SAMPLEGENDER.out.tsv
        .map { meta, tsv ->
            def sex = get_sex(tsv)
            def new_meta = meta + [sex: sex]
            [ new_meta, sex ]
        }
        .mix(ch_ngsbits_input.sex)
        .mix(ch_input.npz.map { meta, _npz -> [ meta, meta.sex ] })

    //
    // Create a small metrics file
    //

    def ch_sex_counts = ch_sexes
        .reduce([:]) { counts, entry ->
            def meta = entry[0]
            def sex = entry[1]
            counts[sex] = (counts[sex] ?: []) + meta.id
            counts
        }

    def ch_metrics = ch_sex_counts.map { sexes -> create_mqc_metrics(sexes) }
        .collectFile(name: "metrics_mqc.tsv")

    ch_multiqc_files = ch_multiqc_files.mix(ch_metrics)

    def ch_metrics_summary = ch_sex_counts
        .map { sexes ->
            def metrics = get_metrics(sexes)
            return [
                "Male/Female ratio: ${metrics.male_to_female_ratio}",
                "Male count: ${metrics.male_count}",
                "Female count: ${metrics.female_count}",
                "Total count: ${metrics.total_count}",
                "Male IDs: ${metrics.males.join(", ")}",
                "Female IDs: ${metrics.females.join(", ")}"
            ].join("\n")
        }
        .collectFile(name: "metrics_summary.txt")


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
    def Date date = new Date()
    def String dateFormat = "WisecondorX_${date.format("ddMMyyyy")}"

    def ch_newref_input = WISECONDORX_CONVERT.out.npz
        .mix(ch_input.npz)
        .map { _meta, npz ->
            def new_meta = [id:prefix ?: dateFormat]
            [ new_meta, npz ]
        }
        .groupTuple() // All files should be present here, so no size is needed
        .combine(val_bin_sizes)
        .map { meta, npz, bin_size ->
            def new_meta = meta + [bin_size:bin_size]
            [ new_meta, npz ]
        }

    WISECONDORX_NEWREF(ch_newref_input)
    ch_versions = ch_versions.mix(WISECONDORX_NEWREF.out.versions.first())

    //
    // Collate and save software versions
    //
    def ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${outdir}/pipeline_info", name: 'nf_cmgg_pipeline_software_mqc_versions.yml', sort: true, newLine: true)

    //
    // MODULE: MultiQC
    //
    def ch_multiqc_config                     = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config              = multiqc_config ? channel.fromPath(multiqc_config, checkIfExists: true) : channel.empty()
    def ch_multiqc_logo                       = multiqc_logo ? channel.fromPath(multiqc_logo, checkIfExists: true) : channel.empty()
    def summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary                   = channel.of(paramsSummaryMultiqc(summary_params))
    def ch_multiqc_custom_methods_description = multiqc_methods_description ?
                                                file(multiqc_methods_description, checkIfExists: true) :
                                                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description                = channel.of(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                          = ch_multiqc_files.mix(
                                                    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
                                                    ch_collated_versions,
                                                    ch_methods_description.collectFile(
                                                        name: 'methods_description_mqc.yaml',
                                                        sort: false
                                                    )
                                                )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report: List<Path>              = MULTIQC.out.report.toList()
    multiqc_plots: Value<Path>              = MULTIQC.out.plots
    multiqc_data: Value<Path>               = MULTIQC.out.data
    npz: Channel<Tuple<Map,Path>>           = WISECONDORX_CONVERT.out.npz
    references: Channel<Tuple<Map,Path>>    = WISECONDORX_NEWREF.out.npz
    metrics: Channel<Path>                  = ch_metrics_summary
    versions: Channel<Path>                 = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def get_sex(tsv) {
    def split_tsv = tsv.splitCsv(sep:"\t", header:true, strip:true)
    return split_tsv[0].gender
}

def create_mqc_metrics(sexes) {
    def metrics = get_metrics(sexes)

    return """# plot_type: 'table'
Male to female ratio\tMale count\tFemale count\tTotal count\tMales\tFemales
${metrics.male_to_female_ratio}\t${metrics.male_count}\t${metrics.female_count}\t${metrics.total_count}\t${metrics.males.join(",")}\t${metrics.females.join(",")}
"""
}

def get_metrics(sexes) {
    return [
        males: sexes["male"],
        females: sexes["female"],
        male_count: sexes["male"].size(),
        female_count: sexes["female"].size(),
        male_to_female_ratio: sexes["male"].size() / sexes["female"].size(),
        total_count: sexes["male"].size() + sexes["female"].size()
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
