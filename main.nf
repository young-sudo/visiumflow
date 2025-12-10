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

process CREATE_ZARR {
  output:
    val "concatenated_sdata", emit zarr_ch

  script:
    """
    python3 scripts/create_zarr.py
    """  
}

process QC_DATA {
  input:
    val zarr_ch

  output:
    val "preprocessed_adata.h5ad", emit preprocessed_ch

  script:
    """
    python3 scripts/qc_data.py
    """
}


process ANALYZE_CLUSTERS {
  input:
    val preprocessed_ch

  script:
    """
    python3 scripts/cluster_analysis.py
    """
}
process ANALYZE_ANNOTATIONS {
  input:
    val preprocessed_ch

  script:
    """
    python3 scripts/annotation_analysis.py
    """
}
process ANALYZE_DEGS {
  input:
    val preprocessed_ch

  script:
    """
    python3 scripts/deg_analysis.py
    """
}

workflow {
  spatial_data_ch = CREATE_ZARR()
  qc_data_ch = QC_DATA(spatial_data_ch)

  ANALYZE_CLUSTERS(qc_data_ch)
  ANALYZE_ANNOTATIONS(qc_data_ch)
  ANALYZE_DEGS(qc_data_ch)
}
