# TCR_nextflow

docker pull milaboratory/mixcr

docker build -t mhecd4tcr:devel .

nextflow run nextflow_main.nf -profile docker -resume -work-dir work
