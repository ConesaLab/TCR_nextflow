# TCR_nextflow

## Get the software

To pull the official MiXCR container

```
docker pull milaboratory/mixcr
```

To pull the container with the R dependencies needed to run this workflow:

```
docker pull ssnnm/mhecd4tcr:0.1.0
```

## Run the workflow

```
nextflow run nextflow_main.nf -profile docker -resume -work-dir work
```
