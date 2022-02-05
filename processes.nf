
/*
* Step 1.1. TCR quantification using MiXCR
*/
process mixcr_analyze {
    cpus "$params.cpus"
    memory "$params.memory"
    time "$params.time"

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
    mixcr analyze shotgun -t $task.cpus --species $params.specie --starting-material rna --only-productive \
    --align "-OsaveOriginalReads=true" \
    ${R1} ${R2} ${SampleID}
    csplit -f ${SampleID}.report ${SampleID}.report '/^==/' '{*}' > mixcr_qc.log
    """
}

/*
* Step 1.2. MiXCR quality control
*/
process mixcr_qc {

    label 'mhecd4tcr'

    publishDir "$params.outdir/01_MiXCR",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> filename.replaceFirst('TCRanalysis_bookdown/','')},
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
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/01_mixcr_qc.Rmd',
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
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("clones_*"), emit: filt_clones
    path("TCRanalysis_bookdown/*"), emit: filt_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/02_datafiltering.Rmd',
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
        saveAs: { filename -> filename.replaceFirst('TCRanalysis_bookdown/','')},
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir1)
    path(inputDir2)
    path(sampleInfo)

    output:
    path("*.html")
    path("03_*Rda")
    path("TCRanalysis_bookdown/*"), emit: overview_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/03_dataset_overview.Rmd',
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
        saveAs: { filename -> filename.replaceFirst('TCRanalysis_bookdown/','')},
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("04_*Rda")
    path("TCRanalysis_bookdown/*"), emit: corr_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/04_correlations.Rmd',
    params=list(
        'inputDir'=here,
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'sampleInfo'='${sampleInfo}',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 5. Overlap analysis
*/
process overlap {

    label 'mhecd4tcr'

    publishDir "$params.outdir/05_Overlap",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> filename.replaceFirst('TCRanalysis_bookdown/','')},
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("05_*Rda")
    path("TCRanalysis_bookdown/*"), emit: overlap_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/05_overlap.Rmd',
    params=list(
        'inputDir'=here,
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'sampleInfo'='${sampleInfo}',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 6. Diversity analysis
*/
process diversity {

    label 'mhecd4tcr'

    publishDir "$params.outdir/06_Diversity",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> filename.replaceFirst('TCRanalysis_bookdown/','')},
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("06_*Rda")
    path("TCRanalysis_bookdown/*"), emit: diversity_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/06_diversity.Rmd',
    params=list(
        'inputDir'=here,
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'sampleInfo'='${sampleInfo}',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 7. K-mers analysis
*/
process kmers {

    label 'mhecd4tcr'

    publishDir "$params.outdir/07_Kmers",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> filename.replaceFirst('TCRanalysis_bookdown/','')},
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("07_*Rda")
    path("TCRanalysis_bookdown/*"), emit: kmers_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/07_kmers.Rmd',
    params=list(
        'inputDir'=here,
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'sampleInfo'='${sampleInfo}',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 8. Network analysis
*/
process network {

    label 'mhecd4tcr'

    publishDir "$params.outdir/08_Network",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> filename.replaceFirst('TCRanalysis_bookdown/','')},
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)

    output:
    path("*.html")
    path("08_*Rda")
    path("TCRanalysis_bookdown/*"), emit: network_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/08_network.Rmd',
    params=list(
        'inputDir'=here,
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'sampleInfo'='${sampleInfo}',
        'chain'='${params.chain}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 9. Databases
*/
process ddbb {

    label 'mhecd4tcr'

    publishDir "$params.outdir/09_Databases",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    publishDir "$params.outdir/TCRanalysis_bookdown/",
        pattern: 'TCRanalysis_bookdown/*',
        saveAs: { filename -> filename.replaceFirst('TCRanalysis_bookdown/','')},
        mode: 'copy',
        overwrite: true

    input:
    path(inputDir)
    path(sampleInfo)
    path(mcpas)
    path(vdjdb)

    output:
    path("*.html")
    path("TCRanalysis_bookdown/*"), emit: ddbb_bookdown

    script:
    """
    Rscript -e "here<-getwd();rmarkdown::render('${projectDir}/scripts/09_ddbb.Rmd',
    params=list(
        'inputDir'=here,
        'workDir'=here,
        'outputDir'='TCRanalysis_bookdown',
        'mcpas'='${mcpas}',
        'vdjdb'='${vdjdb}',
        'sampleInfo'='${sampleInfo}',
        'chain'='${params.chain}',
        'specie'='${params.specie}'),
    'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
    """
}

/*
* Step 10. Report
*/
process report {

    label 'mhecd4tcr'

    publishDir "$params.outdir",
        pattern: '*',
        mode: 'copy',
        overwrite: true

    input:
    path(inputs)
    path(report_rmd)

    output:
    path("*.html")

    script:
    """
    Rscript -e "bookdown::render_book(
                    input='${report_rmd}',
                    params=list(
                        'inputDir'=getwd(),
                        'chain'='${params.chain}',
                        'specie'='${params.specie}'))"
    mv _main.html final-report_${params.chain}.html
    """
}