#!/usr/bin/env nextflow

if ( params.help ) {
    log.info """
    Usage:
      nextflow run main.nf -profile conda

    Parameters:
       -profile           standard|conda|docker|singularity|slurm
    """
    exit 0
}


download.sh

create_zarr, qc_data

cluster_analysis

annotation_analysis

deg_analysis

workflow {

}

