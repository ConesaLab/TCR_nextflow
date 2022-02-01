#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    mixcr_analyze; mixcr_qc;
    data_filtering; dataset_overview; correlations; overlap;
    diversity; kmers; network; ddbb; report
} from './processes.nf'

log.info """\
RNASeq - TCR Pipeline
===================================
Project parameters:
- Project Name          : ${params.pn}
- Sample List           : ${params.readfiles}
- Specie                : ${params.specie}
- Selected TCR chain    : ${params.chain}
- Output Directory      : ${params.outdir}
  """



workflow {

    // Input validation
    def valid_species = ['HomoSapiens','MusMusculus']
    is_valid_specie = params.specie in valid_species
    if (!is_valid_specie) {
        log.error "`params.specie` must be one of ${valid_species}"
        exit 1, "`params.specie` must be one of ${valid_species}"
    }

    Channel.fromPath("${params.readfiles}")
            .ifEmpty { exit 1, "File not foud: ${params.readfiles}" }
            .set { sampleInfoChannel }

    Channel.fromPath("${params.mcpas}")
            .ifEmpty { exit 1, "File not foud: ${params.mcpas}" }
            .set { mcpasChannel }
    Channel.fromPath("${params.vdjdb}")
            .ifEmpty { exit 1, "File not foud: ${params.vdjdb}" }
            .set { vdjdbChannel }

    /*
   * Create a channel that emits tuples containing three elements:
   * the SampleID, the first read-pair file and the second read-pair file
   */
    samples_channel = sampleInfoChannel
        .splitCsv(header:true)
        .map { row -> tuple (row.SampleID,
        file(row.R1),
        file(row.R2))}

    /*
    * Workflow functions
    */

    mixcr_analyze(samples_channel)

    mixcr_qc(mixcr_analyze.out.full_report_chunks.collect(), sampleInfoChannel)

    data_filtering(mixcr_analyze.out.all_clonotypes.collect(), sampleInfoChannel)

    dataset_overview(mixcr_qc.out.qc_bookdown.collect(), data_filtering.out.filt_bookdown.collect(), sampleInfoChannel)

    correlations(dataset_overview.out.overview_bookdown.collect(), sampleInfoChannel)

    overlap(data_filtering.out.filt_clones.collect(), sampleInfoChannel)

    diversity(data_filtering.out.filt_clones.collect(), sampleInfoChannel)

    kmers(data_filtering.out.filt_clones.collect(), sampleInfoChannel)

    network(data_filtering.out.filt_clones.collect(), sampleInfoChannel)

    ddbb(data_filtering.out.filt_clones.collect(), sampleInfoChannel, mcpasChannel, vdjdbChannel)

    report(
        mixcr_qc.out.qc_bookdown
        .mix(data_filtering.out.filt_bookdown)
        .mix(dataset_overview.out.overview_bookdown)
        .mix(correlations.out.corr_bookdown)
        .mix(overlap.out.overlap_bookdown)
        .mix(diversity.out.diversity_bookdown)
        .mix(kmers.out.kmers_bookdown)
        .mix(network.out.network_bookdown)
        .mix(ddbb.out.ddbb_bookdown)
        .collect(),
        file('data/scripts/10_report.Rmd')
    )

}
