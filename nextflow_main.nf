#!/usr/bin/env nextflow

/*
 * Defines some parameters in order to specify the read pairs,
 * data information and TCR databases by using the command line options
 */
params.reads = "$baseDir/data/reads/*{1,2}.fastq.gz"
params.mcpas = "$baseDir/data/ddbb/McPAS-TRB_human.csv"
params.vdjdb = "$baseDir/data/ddbb/VDJdb-TRB_human.tsv"
params.outdir = 'results'

READFILES = params.readfiles
BASENAME  = params.bn


log.info """\
RNASeq - TCR Pipeline
===================================
Project parameters:
- BASENAME              : ${BASENAME}
- MANIFEST              : ${READFILES}
  """

/*
 * Create a channel that emits tuples containing three elements:
 * the SampleID, the first read-pair file and the second read-pair file
 */

Channel
    .fromPath( "$READFILES" )
    .splitCsv(header:true)
    .map { row -> tuple (row.SampleID,
        file(row.R1),
        file(row.R2))}
    .into { samples_channel; samples_channel_2 }


/*
* Step 1. TCR quantification using MiXCR
*/
process mixcr_analyze {

    publishDir "$params.outdir"

    input:
    set SampleID, file(R1), file(R2) from samples_channel

    output:
    file "${SampleID}.clonotypes.ALL.txt" into all_clonotypes
    file "${SampleID}.report" into full_report

    script:
    """
    mixcr analyze shotgun -t 3 --species hs --starting-material rna --only-productive \
    $baseDir/data/reads/${R1} $baseDir/data/reads/${R2} ${SampleID}
    """
}


/*
* Step 2. MiXCR quality control

* process mixcr_qc {

*      input:
*      env report from full_report

*      script:
*      """
*      for file in ${full_report}; do
*         csplit -f \$file \$file '/^==/' '{*}'
*         rm "*.report05"
*     done
*      """
*  }
*/

/*
* Step 3. Diversity Analysis
*
*process diversity {

*      input:
*      env report from full_report

*      script:
*      """
*      Rscript data/scripts/02_datafiltering_evenness.Rmd
*      """
*  }
*/
