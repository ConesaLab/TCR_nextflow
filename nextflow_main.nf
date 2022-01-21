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
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> Path.of(filename).getName() },
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("TCRanalysis_bookdown/*"), emit: qc_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/data/scripts/01_mixcr_qc.Rmd',
    params=list(
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'sampleInfo'='${sampleInfo}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 2. Data filtering
*/
process data_filtering {

    label 'mhecd4tcr'

    publishDir "$params.outdir/02_DataFiltering",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> Path.of(filename).getName() },
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("clones_*")
    path("TCRanalysis_bookdown/*"), emit: filt_bookdown
    
    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/data/scripts/02_datafiltering.Rmd', 
    params=list(
        'inputDir'=here, 
        'workDir'=here, 
        'outputDir'='TCRanalysis_bookdown', 
        'sampleInfo'='${sampleInfo}'), 
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 3. Dataset overview
*/
process dataset_overview {

    label 'mhecd4tcr'

    publishDir "$params.outdir/03_DatasetOverview",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> Path.of(filename).getName() },
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir1)
    path(inputDir2)
    path(sampleInfo)

    output:
    path("*.html")
    path("TCRanalysis_bookdown/*"), emit: overview_bookdown
    
    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/data/scripts/03_dataset_overview.Rmd',
    params=list(
        'inputDir1'=here,
        'inputDir2'=here,
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'sampleInfo'='${sampleInfo}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 4. Correlations
*/
process correlations {

    label 'mhecd4tcr'

    publishDir "$params.outdir/04_Correlations",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> Path.of(filename).getName() },
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("TCRanalysis_bookdown/*"), emit: corr_bookdown
    
    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/data/scripts/04_correlations.Rmd', 
    params=list(
        'inputDir'=here, 
        'workDir'=here, 
        'outputDir'='TCRanalysis_bookdown', 
        'sampleInfo'='${sampleInfo}',
        'chain'='${params.chain}'), 
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

    dataset_overview(mixcr_qc.out.qc_bookdown.collect(), data_filtering.out.filt_bookdown.collect(), sampleInfoChannel)

    correlations(dataset_overview.out.overview_bookdown.collect(), sampleInfoChannel)


}
