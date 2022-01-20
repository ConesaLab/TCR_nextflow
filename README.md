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


### Local

Example commands:

Use the docker profile to use the containers with MiXCR and the R packages, and use cached results:

```
nextflow run nextflow_main.nf -profile docker -resume -work-dir work
```

### Cluster

If you want to use a specific path to save the singulartiy container, you need to export `NXF_SINGULARITY_CACHEDIR` before running the pipeline.
Then you run the pipeline specifying the profile you want to use for your cluster. The predefined profile are in the file `nextflow.config`.

```
export NXF_SINGULARITY_CACHEDIR=/path/to/my/singularity/images/
nextflow run nextflow_main.nf -profile singularity,slurm
```