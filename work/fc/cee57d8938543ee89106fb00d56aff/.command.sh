#!/bin/bash -ue
Rscript -e "here<-getwd();rmarkdown::render('/home/teresa/Documents/GitHub/TCR_nextflow/data/scripts/07_kmers.Rmd', 
params=list(
    'inputDir'=here, 
    'workDir'=here, 
    'outputDir'='TCRanalysis_bookdown', 
    'sampleInfo'='sampleslist.csv',
    'chain'='TRB'), 
'output_dir'= here, 'knit_root_dir'=here, quiet=TRUE)"
