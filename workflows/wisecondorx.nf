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

include { Input } from '../subworkflows/local/utils_nfcore_wisecondorx_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WISECONDORX {

    take:
    ch_samplesheet: Channel<Input>      // samplesheet read in from --input
    fasta: String                       // the reference fasta file
    fai: String                         // the index of the reference fasta file
    val_bin_sizes: List<Integer>        // a list of bin sizes to use
    prefix: String                      // the prefix to be used by the output file
    outdir: String                      // the path of the output directory
    multiqc_config: String              // the path to the multiqc config
    multiqc_logo: String                // the path to the multiqc logo
    multiqc_methods_description: Path   // the file containing the multiqc custom method descriptions

    main:

    def ch_versions: Channel<Path> = channel.empty()
    def ch_multiqc_files: Channel<Path> = channel.empty()

    //
    // Create optional input files
    //

    def ch_ref: Value<Reference> = channel.value(record(fasta: file(fasta), id: 'reference'))

    if(!fai) {
        SAMTOOLS_FAIDX(
            ch_ref
        )
        ch_ref = ch_ref.combine(SAMTOOLS_FAIDX.out).map { rec, out ->
            rec + record(fai: out.fai)
        }
    } else {
        ch_ref = ch_ref.map { rec -> rec + record(fai: file(fai)) }
    }

    def ch_npz: Channel<Input> = ch_samplesheet.filter { rec -> rec.npz }
    def ch_cram: Channel<Input> = ch_samplesheet.filter { rec -> !rec.npz }

    //
    // Index the non-indexed input files
    //

    def ch_index_input: Channel<Input> = ch_cram.filter { rec -> !rec.crai }
    SAMTOOLS_INDEX(ch_index_input)

    // TODO records are not supported by .join yet, update this once it is
    def ch_indexed: Channel<Input> = ch_cram.filter { rec -> rec.crai }.map { rec -> tuple(rec.id, rec)}
        .join(SAMTOOLS_INDEX.out.map { rec -> tuple(rec.id, rec) })
        .map { _id, rec1, rec2 -> rec1 + rec2 }
        .mix(ch_cram.filter { rec -> rec.crai })

    //
    // Define the sex if it's not given
    //

    def ch_no_sex: Channel<Input> = ch_indexed.filter { rec -> !rec.sex }
    NGSBITS_SAMPLEGENDER(
        ch_no_sex.map { rec -> rec + record(method:'xy')},
        ch_ref
    )

    def ch_sexes: Channel<Input> = NGSBITS_SAMPLEGENDER.out
        .map { rec ->
            def sex: String = get_sex(rec.tsv)
            record(id: rec.id, sex: sex) as Sex
        }
        .mix(ch_indexed.filter { rec -> rec.sex })
        .mix(ch_npz)

    //
    // Create a small metrics file
    //

    def ch_sex_counts: Channel<Map<String, List<String>>> = ch_sexes
        .reduce(["male":[], "female": []]) { counts: Map<String, List<String>>, rec ->
            counts[rec.sex] = counts[rec.sex] + rec.id
            counts
        }

    def ch_metrics: Value<Path> = ch_sex_counts.map { sexes -> create_mqc_metrics(sexes) }
        .collectFile(name: "metrics_mqc.tsv")
        .collect()

    ch_multiqc_files = ch_multiqc_files.mix(ch_metrics)

    def ch_metrics_summary: Value<Path> = ch_sex_counts.view()
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
        .collect()


    //
    // Convert the input files to NPZ files
    //

    def ch_wcx_npz: Channel<Input> = WISECONDORX_CONVERT(ch_indexed, ch_ref)
    def ch_all_npz: Channel<Input> = ch_wcx_npz.mix(ch_npz)

    //
    // Create the WisecondorX reference
    //

    // Define reference name (with timestamp) => only used when --prefix is null
    def Date date = new Date()
    def String dateFormat = "WisecondorX_${date.format("ddMMyyyy")}"

    def ch_newref_input: Channel<Record> = ch_all_npz
        .collect() // All files should be present here, so no size is needed
        .map { recs: List<Record> ->
            def list_npz: List<Path> = recs.collect { rec -> rec.npz }
            record(
                id: prefix ?: dateFormat,
                npzs: list_npz
            )

        }
        .combine(val_bin_sizes)
        .map { rec, bin ->
            rec + record(bin_size: bin)
        }

    def ch_refs: Channel<Record> = WISECONDORX_NEWREF(ch_newref_input)

    //
    // Collate and save software versions
    //

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'structural_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
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
    npz: Channel<Record>                    = ch_wcx_npz
    references: Channel<Record>             = ch_refs
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
    def metrics: Map<String, ?> = [
        males: sexes["male"],
        females: sexes["female"],
        male_count: sexes["male"].size(),
        female_count: sexes["female"].size(),
        total_count: sexes["male"].size() + sexes["female"].size()
    ]
    if(metrics.male_count == 0 || metrics.female_count == 0) {
        metrics['male_to_female_ratio'] = 'NA'
    } else { 
        metrics['male_to_female_ratio'] = metrics.male_count / metrics.female_count
    }
    return metrics
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RECORDS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Reference {
    id: String
    fasta: Path
    fai: Path
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Sex {
    id: String
    sex: String
}