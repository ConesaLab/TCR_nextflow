#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
RNASeq - TCR Pipeline
===================================
Project parameters:
- BASENAME              : ${params.bn}
- MANIFEST              : ${params.readfiles}
  """


/*
* Step 1.1. TCR quantification using MiXCR
*/
process mixcr_analyze {

    tag "$SampleID"
    label "mixcr"
    label "all_cpu"
    cache 'lenient'

    publishDir "$params.outdir/01_MiXCR",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(SampleID), file(R1), file(R2)

    output:
    path("${SampleID}.clonotypes.ALL.txt"), emit: all_clonotypes
    path("${SampleID}.report"), emit: full_report
    path("${SampleID}.report0[0-4]"), emit: full_report_chunks

    script:
    """
    mixcr analyze shotgun -t $task.cpus --species hs --starting-material rna --only-productive \
    ${R1} ${R2} ${SampleID}
    csplit -f ${SampleID}.report ${SampleID}.report '/^==/' '{*}' > mixcr_qc.log
    """
}

/*
* Step 1.2. MiXCR quality control
*/
process mixcr_qc {

    cpus 1
    label 'mhecd4tcr'

    publishDir "$params.outdir/01_MiXCR",
        mode: 'copy',
        overwrite: true

    input:
    path(report)
    path(sampleInfo)

    output:
    path("*.html"), emit: qc_report
    path(TCRanalysis_articlefigures)
    path(TCRanalysis_bookdown)

    script:
    // TODO: parametrize sampleLevels
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/data/scripts/01_mixcr_qc.Rmd', \
    params=list('workDir'=here, \
    'outputDir'='TCRanalysis_bookdown', \
    'articleDir'='TCRanalysis_articlefigures', \
    'sampleInfo'='${sampleInfo}', \
    'sampleLevels'=c('control', 'withoutMHE', 'withMHE')), \
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 2. Data filtering
*/
process data_filtering {

    label 'mhecd4tcr'

    publishDir "$params.outdir/02_DataFiltering",
        mode: 'copy',
        overwrite: true

    input:
    path(report)
    path(sampleInfo)

    output:
    path("*.html"), emit: qc_report
    path(TCRanalysis_bookdown)
    // TODO: parametrize levels
    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/data/scripts/02_datafiltering.Rmd', \
    params=list( \
        'inputDir'=here, \
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'sampleInfo'='${sampleInfo}',
        'sampleLevels'=c('control', 'withoutMHE', 'withMHE')), \
        'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}


workflow {

   Channel.fromPath("${params.readfiles}")
            .ifEmpty { exit 1, "File not foud: ${params.readfiles}" }
            .set { sampleInfoChannel }

   /*
   * Create a channel that emits tuples containing three elements:
   * the SampleID, the first read-pair file and the second read-pair file
   */
    samples_channel = sampleInfoChannel
        .splitCsv(header:true)
        .map { row -> tuple (row.SampleID,
        file(row.R1),
        file(row.R2))}

    mixcr_analyze(samples_channel)

    mixcr_qc(mixcr_analyze.out.full_report_chunks.collect(), sampleInfoChannel)

    data_filtering(mixcr_analyze.out.all_clonotypes.collect(), sampleInfoChannel)
}