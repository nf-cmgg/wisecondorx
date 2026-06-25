nextflow.enable.types = true

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

record Input {
    id: String
    input: Path?
    input_idx: Path?
    npz: Path?
    sex: String?
}

record InputReady {
    id: String
    input: Path?
    input_idx: Path?
    npz: Path?
    sex: String?
    fasta: Path
    fai: Path
}

workflow WISECONDORX {

    take:
    ch_samplesheet: Channel<Input>      // samplesheet read in from --input
    fasta: Path                         // the reference fasta file
    fai: Path                           // the index of the reference fasta file
    val_bin_sizes: List<Integer>        // a list of bin sizes to use
    prefix: String                      // the prefix to be used by the output file
    outdir: String                      // the path of the output directory
    multiqc_config: Path?               // the path to the multiqc config
    multiqc_logo: Path?                 // the path to the multiqc logo
    multiqc_methods_description: Path?  // the file containing the multiqc custom method descriptions

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Create optional input files
    //

    def ch_fai: Value<Record> = channel.empty()
    if(!fai) {
        SAMTOOLS_FAIDX(
            record(
                id: 'fasta',
                fasta: fasta
            )
        )
        ch_fai = SAMTOOLS_FAIDX.out
    } else {
        ch_fai = channel.of(record(id: "fasta", fai: fai))
    }

    def ch_input: Channel<InputReady> = ch_samplesheet
        .combine(ch_fai)
        .map { input_rec: Record, fai_rec: Record ->
            input_rec + record(fasta: fasta, fai: fai_rec.fai)
        }

    def ch_npz: Channel<InputReady> = ch_input.filter { rec -> rec.npz }
    def ch_crams: Channel<InputReady> = ch_input.filter { rec -> rec.input && !rec.npz }

    //
    // Index the non-indexed input files
    //

    def ch_crams_to_index = ch_crams.filter { rec -> !rec.input_idx }
    def ch_crams_with_index = ch_crams.filter { rec -> rec.input_idx }

    def index_out = SAMTOOLS_INDEX(ch_crams_to_index)
    def ch_indexed: Channel<Record> = ch_crams_with_index.mix(index_out)

    //
    // Define the sex if it's not given
    //

    def ch_indexed_with_sex: Channel<Record> = ch_indexed.filter { rec -> rec.sex }
    def ch_indexed_without_sex: Channel<Record> = ch_indexed.filter { rec -> !rec.sex }

    def ngsbits_out = NGSBITS_SAMPLEGENDER(
        ch_indexed_without_sex.map { rec -> rec + record(method: 'xy')}
    )

    def ch_sexes: Channel<Record> = ngsbits_out
        .map { rec ->
            def sex = get_sex(rec.tsv)
            rec + record(sex: sex)
        }
        .mix(ch_indexed_with_sex)
        .mix(ch_npz)

    //
    // Create a small metrics file
    //

    def ch_sex_counts = ch_sexes
        .reduce([:]) { counts: Map<String, List<String>>, rec: Record ->
            def sex = rec.sex
            counts[sex] = (counts[sex] ?: []) + rec.id
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

    def convert_out = WISECONDORX_CONVERT(
        ch_indexed
    )

    //
    // Create the WisecondorX reference
    //

    // Define reference name (with timestamp) => only used when --prefix is null
    def Date date = new Date()
    def String dateFormat = "WisecondorX_${date.format("ddMMyyyy")}"

    def ch_newref_input = convert_out
        .mix(ch_npz)
        .map { rec ->
            tuple('id', record(id: prefix ?: dateFormat, npz: rec.npz))
        }
        .groupBy() // All files should be present here, so no size is needed
        .map { key: String, items: Bag<Record> ->
            record(
                id: items[0].id,
                inputs: items*.npz
            )
        }
        .combine(channel.fromList(val_bin_sizes))
        .map { rec, bin_size ->
            rec + record(bin_size: bin_size)
        }

    def newref_out = WISECONDORX_NEWREF(ch_newref_input)

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

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_cmgg_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'nf-cmgg/wisecondorx'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    multiqc_plots  = MULTIQC.out.plots
    multiqc_data   = MULTIQC.out.data
    npz: Channel<Record>            = convert_out // WISECONDORX_CONVERT.out.npz // channel: [ val(meta), path(/path/to/npz_file.npz) ]
    references: Channel<Record>     = newref_out  // channel: [ val(meta), path(/path/to/reference.npz) ]
    metrics        = ch_metrics_summary
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def get_sex(tsv) {
    println tsv
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
