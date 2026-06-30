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
include { WISECONDORX_NEWREF; WisecondorxNewrefInput          } from '../modules/nf-core/wisecondorx/newref/main'
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

record Sample {
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
    multiqc_config: Path                // the path to the multiqc config
    multiqc_logo: Path                  // the path to the multiqc logo
    _multiqc_methods_description: Path   // the file containing the multiqc custom method descriptions

    main:

    def ch_versions: Channel<Path> = channel.empty()
    def ch_multiqc_files: Channel<Path> = channel.empty()

    //
    // Create optional input files
    //

    def ch_fai: Value<Record> = channel.empty().collect()
    if(!fai) {
        ch_fai = SAMTOOLS_FAIDX(
            record(
                id: 'fasta',
                fasta: fasta,
                get_sizes: false
            )
        )
    } else {
        ch_fai = channel.value(record(id: "fasta", fai: fai))
    }

    def ch_input: Channel<Sample> = ch_samplesheet
        .combine(ch_fai)
        .map { input_rec, fai_rec ->
            input_rec + record(fasta: fasta, fai: fai_rec.fai)
        }

    def ch_npz: Channel<Sample> = ch_input.filter { rec -> rec.npz }
    def ch_crams: Channel<Sample> = ch_input.filter { rec -> rec.input && !rec.npz }

    //
    // Index the non-indexed input files
    //

    def ch_crams_to_index: Channel<Sample> = ch_crams.filter { rec -> !rec.input_idx }
    def ch_crams_with_index: Channel<Sample> = ch_crams.filter { rec -> rec.input_idx }

    def ch_indexed: Channel<Sample> = ch_crams_with_index.mix(
        ch_crams_to_index.join(SAMTOOLS_INDEX(ch_crams_to_index), by: 'id')
    )

    //
    // Define the sex if it's not given
    //

    def ch_indexed_with_sex: Channel<Sample> = ch_indexed.filter { rec -> rec.sex }
    def ch_indexed_without_sex: Channel<Sample> = ch_indexed.filter { rec -> !rec.sex }

    def ch_sexes: Channel<Sample> = NGSBITS_SAMPLEGENDER(
            ch_indexed_without_sex
                .map { rec -> rec + record(method: 'xy') }
        ).join(ch_indexed_without_sex, by: 'id')
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
        .reduce([:]) { counts, rec ->
            def sex = rec.sex
            counts[sex] = (counts[sex] ?: []) + rec.id
            counts
        }

    def ch_metrics: Value<Path> = ch_sex_counts.map { sexes -> 
        def metrics_text: String = create_mqc_metrics(sexes)
        def f: Path = workflow.workDir.resolve("collectfiles-${workflow.sessionId}/metrics_mqc.txt")
        return create_file(f, metrics_text)
    }

    ch_multiqc_files = ch_multiqc_files.mix(ch_metrics)

    def ch_metrics_summary = ch_sex_counts
        .map { sexes ->
            def metrics = get_metrics(sexes)
            def metrics_summary: String = [
                "Male/Female ratio: ${metrics.male_to_female_ratio}",
                "Male count: ${metrics.male_count}",
                "Female count: ${metrics.female_count}",
                "Total count: ${metrics.total_count}",
                "Male IDs: ${metrics.males.join(", ")}",
                "Female IDs: ${metrics.females.join(", ")}"
            ].join("\n")
            return create_file(workflow.workDir.resolve("${params.outdir}/metrics_summary.txt"), metrics_summary)
        }

    //
    // Convert the input files to NPZ files
    //

    def convert_out = ch_indexed.join(
        WISECONDORX_CONVERT(
            ch_indexed
        ), by: 'id')

    //
    // Create the WisecondorX reference
    //

    // Define reference name (with timestamp) => only used when --prefix is null
    def day: Integer = workflow.start.getDayOfMonth()
    def month: String = workflow.start.getMonthValue() > 10 ? 
        "${workflow.start.getMonthValue()}" : 
        "0${workflow.start.getMonthValue()}"
    def year: Integer = workflow.start.getYear()
    def dateFormat: String = "WisecondorX_${day}${month}${year}"

    def ch_newref_input: Channel<WisecondorxNewrefInput> = convert_out
        .mix(ch_npz)
        .map { rec ->
            tuple('id', record(id: prefix ?: dateFormat, npz: rec.npz))
        }
        .groupBy() // All files should be present here, so no size is needed
        .map { _key, items ->
            record(
                id: items.toList()[0].id,
                inputs: items*.npz
            )
        }
        .combine(channel.fromList(val_bin_sizes))
        .map { rec, bin_size ->
            rec + record(bin_size: bin_size)
        }

    def newref_out = ch_newref_input.join(
        WISECONDORX_NEWREF(ch_newref_input), by: 'id'
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .unique()

    def topic_versions_file: Channel<Path> = topic_versions.filter { entry ->
        entry instanceof Path
    }
    def topic_versions_tuple: Channel<Tuple<String, String, String>> = topic_versions.filter { entry ->
        entry instanceof Tuple
    }

    def topic_versions_string = topic_versions_tuple
        .map { process, tool, version ->
            tuple(process.substring(process.lastIndexOf(':')+1), "  ${tool}: ${version}")
        }
        .groupBy()
        .map { process, tool_versions ->
            "${process}:\n${tool_versions.toSet().toSorted().join("\n")}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions_file))
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
    def ch_workflow_summary = channel.value(
        create_file(
            workflow.workDir.resolve("collectfiles-${workflow.sessionId}/workflow_summary_mqc.yaml"),
            paramsSummaryMultiqc(ch_summary_params)
        )
    )
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary)

    def multiqc_config_set: Set<Path> = multiqc_config ? 
        [multiqc_config].toSet() : 
        [file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)].toSet()

    def multiqc_out = MULTIQC(
        ch_multiqc_files.collect().map { files ->
            record(
                id: 'nf-cmgg/wisecondorx',
                multiqc_files: files.toSet(),
                multiqc_config: multiqc_config_set,
                multiqc_logo: multiqc_logo
            )
        }
    ).map { rec ->
        record(
            id: rec.id,
            report: rec.report,
            data: rec.data,
            plots: rec.plots
        )
    }

    emit:
    multiqc: Value<Record> = multiqc_out
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

def create_file(tmp_file: Path, content: String) -> Path {
    tmp_file.parent.mkdirs()
    tmp_file.text = content
    return tmp_file
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
