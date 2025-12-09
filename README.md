# Workflow for Spatial Transcriptomics datasets

<i>by Younginn Park</i>

![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Scanpy](https://img.shields.io/badge/Scanpy-ef3284?style=for-the-badge&logo=python&logoColor=white)
![Squidpy](https://img.shields.io/badge/Squidpy-949df0?style=for-the-badge&logo=python&logoColor=white)
![SpatialData](https://img.shields.io/badge/SpatialData-f1e63a?style=for-the-badge&logo=python&logoColor=white)
![PyDESeq2](https://img.shields.io/badge/PyDESeq2-2cb446?style=for-the-badge&logo=python&logoColor=white)
![Geopandas](https://img.shields.io/badge/GeoPandas-139C5A?style=for-the-badge&logo=python&logoColor=white)
![Nextflow](https://img.shields.io/badge/Nextflow-23CC85?style=for-the-badge&logo=nextflow&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-2496ED?style=for-the-badge&logo=docker&logoColor=white)

<p style="margin-top: 10px; margin-bottom: 10px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/10x_logo.png" alt="10x" width=100>
</p>

This project implements the <b>Visium HD Multi-sample Analysis</b> workflow from 10x Genomics, based on their [Tutorial](https://www.10xgenomics.com/analysis-guides/tutorial-visium-hd-multi-sample-python-colab) and [Colab Notebook](https://colab.research.google.com/github/10XGenomics/analysis_guides/blob/main/Visium_HD_multi_sample_comparison_python.ipynb#scrollTo=2WU0zKjjCKEI), transformed into a reproducible workflow that can be executed across different environments.

<b>Project under active development</b>

# Usage

## Clone repository

```bash
git clone https://github.com/young-sudo/spatialflow.git
```

## Install Nextflow

```bash
# Download Nextflow
curl -s https://get.nextflow.io | bash

# Make Nextflow executable
chmod +x nextflow

# Move Nextflow into an executable path
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/

# Confirm Nextflow is installed correctly
nextflow info
```

## Basic run

```bash
nextflow run main.nf \
  -profile conda
```

Available profiles:
- `standard` - locally without containers
- `conda` - with default mode and Conda environment
- `docker` - with Docker
- `singularity` - with Singularity/Apptainer
- `slurm` - with slurm on HPC

