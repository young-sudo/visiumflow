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

The datasets used in this Analysis Guide are publicly available from 10x Genomics:

[Human Colon Cancer and Normal Adjacent Datasets](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression/dataset-human-crc)

<p style="margin-top: 10px; margin-bottom: 30px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/samples.png" alt="10x" width=300>
<br>
<small>Sample tissues
</p>

The specific datasets used here are:
* **Visium HD, Sample P1 Colon Cancer (CRC)**: Human colon cancer tissue (Patient 1).
* **Visium HD, Sample P2 CRC**: Human colon cancer tissue (Patient 2).
* **Visium HD, Sample P3 Normal Adjacent Tissue (NAT)**: Normal adjacent human colon tissue (Patient 3).
* **Visium HD, Sample P5 NAT**: Normal adjacent human colon tissue (Patient 5).

If you are working with your own data, for each dataset, the `outs` directory will contain the cell segmentation based binned output and spatial output data. For a more detailed description of Visium HD's outputs, see the documentation on our support [site](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview).

The data was originally processed using `spaceranger count` v3.0.0. However, to generate the Space Ranger cell segmentation outputs used in this guide, the public colon cancer and normal adjacent tissue datasets were reprocessed using `spaceranger count` v4.0.1.

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

### For usage help run

```bash
nextflow run main.nf --help
```

# Methods

## The SpatialData Object and Its Components

<p style="margin-top: 10px; margin-bottom: 10px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/visium.png" alt="10x" width=300>
</p>

The `spatialdata` library is built to manage and analyze multiomic spatial datasets. It brings together multiple data types into a single, unified `SpatialData` object. These objects act as on-disk containers that utilize the Zarr file format to store various **Elements** or data types.
A `SpatialData` object created from Visium HD data typically has the following Elements:
* **Images**: CytAssist and microscopy images (e.g., H&E, fluorescence) providing spatial context. These can be accessed via `sdata.images`.
* **Shapes**: Geometric annotations such as polygons or circles representing regions of interest, cells, or spots. In Visium HD, these often represent binned regions or cell segmentations. These are accessible via `sdata.shapes`.
* **Tables**: An `AnnData` object associated with the spatial elements, typically containing gene expression data, cellular metadata, and computational results (e.g., clusters, UMAP embeddings). This Element is used for downstream analyses and is accessed via `sdata.tables`. Each `AnnData` table within `sdata.tables` has:
  * `.X`: The primary data matrix (e.g., raw counts, normalized counts).
  * `.obs`: Observation metadata (e.g., sample ID, cluster assignments).
  * `.var`: Variable metadata (e.g., gene names).
  * `.obsm`: Multi-dimensional annotations (e.g., PCA, UMAP embeddings).
  * `.layers`: Alternative representations of `.X`.

In this project, the `SpatialData` object is created from the Zarr files. Each Zarr file is read into a list of `SpatialData` objects using the `read_zarr` function before the objects are concatenated. In addition, a `sample` column is added to each `AnnData` table within the `SpatialData` object, and the `var_names_make_unique` function is used to ensure that gene names are unique.

## Clustering and Visualization

Now that we have standardized the data and performed a PCA, we will cluster and visualize the results. `Scanpy`'s `neighbors` function generates a neighbor distance matrix and a neighborhood graph, which is used by `Scanpy`'s `leiden` function to cluster the data. Finally, `Scanpy`’s `umap` function is used to visualize the results.


<p style="margin-top: 10px; margin-bottom: 10px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/umap_by_clusters.png" alt="10x" width=250>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/umap_by_sample.png" alt="10x" width=250>
<br>
<small>UMAP visualizations of cells, colored by clusters and sample</small>
</p>


### Clustering options in Scanpy

When running `Scanpy`'s `neighbors` function, the distance metric selected will impact the clustering results. Therefore, you may need to explore different metrics depending on your datasets. Common distance metrics available in `Scanpy` include:
* **Euclidean distance**: The straight-line distance between two points in multi-dimensional space. It tends to group cells with similar overall expression magnitudes. If some genes have very high expression, they can dominate the distance calculation.
* **Manhattan distance (L1 distance, City Block distance)**: The sum of the absolute differences of their coordinates. It is less sensitive to outliers than Euclidean distance. It is useful when differences in individual features are more important than overall magnitude.
* **Cosine distance/similarity**: Measures the angle between two vectors. A smaller angle (closer to 0) indicates higher similarity. It focuses on the orientation of the expression profiles, rather than their magnitude. This is particularly useful when the relative proportions of gene expression might be more informative than the absolute counts. Cell segmentation bins expressing the same genes in similar proportions will be considered close, even if one has a higher total count.
* **Correlation-based distances (e.g., Pearson, Spearman)**: These typically define distance as `1 - correlation_coefficient`. Cell segmentation bins are considered similar if their gene expression profiles are highly correlated, regardless of absolute expression values. This method identifies cell segmentation bins with similar patterns of gene activity.

<p style="margin-top: 15px; margin-bottom: 15px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/cell_dist_clusters.png" alt="10x" width=400>
<br>
<small>Cell distribution across clusters</small>
</p>


## Cell annotations

<p style="margin-top: 15px; margin-bottom: 15px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample1_clust.png" alt="10x" width=150>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample2_clust.png" alt="10x" width=150>
<br>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample3_clust.png" alt="10x" width=150>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample4_clust.png" alt="10x" width=150>
<br>
<small>Sample tissues colored by cell clusters</small>
</p>

A metric that emphasizes magnitude (like Euclidean) might connect cells based on overall transcriptional activity, while a metric emphasizing shape (like Cosine) might connect cells with similar gene expression patterns, even if their total RNA content differs.

For this analysis, we used `Scanpy`'s correlation distance metric, the clustering resolution (RES) is set to 0.8, and the number of neighbors is set to 15. These parameters will likely need to be fine-tuned for new analyses: a smaller resolution generally leads to fewer clusters, while increasing the number of neighbors will have a similar effect.

In addition, for optimal visualization, you may need to adjust the `min_dist` and `spread` parameters used in `Scanpy`'s umap function for Visium HD Gene Expression data. `min_dist` controls how tightly packed the points are in the final embedding, while `spread` determines the overall scale and the separation between clusters. These values will need to be empirically determined.

<p style="margin-top: 15px; margin-bottom: 15px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample1_anno.png" alt="10x" width=150>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample2_anno.png" alt="10x" width=150>
<br>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample3_anno.png" alt="10x" width=150>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample4_anno.png" alt="10x" width=150>
<br>
<small>Sample tissues colored by cell annotations</small>
</p>


For additional resources on plotting `SpatialData` objects, see `Spatialdata`'s Visium HD technology-focused [tutorial](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/notebooks/examples/technology_visium_hd.html).

### Marker genes

<p style="margin-top: 15px; margin-bottom: 15px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/marker_genes_canon.png" alt="10x" width=400>
<br>
<small>Canonical marker genes in clusters</small>
</p>

Using the canonical markers, we observe that Cluster 3 expresses smooth muscle markers.

Cluster annotation can often be challenging when only canonical markers are used. To assist in this process, we can use `Scanpy`’s `rank_genes_groups` function to identify marker genes for each cluster. The results can be ranked by marker score or by the log fold-change. The top-ranked genes within each cluster can then be further analyzed using tools like [Enrichr](https://maayanlab.cloud/Enrichr/) to infer the cluster's potential cell type.

<p style="margin-top: 15px; margin-bottom: 15px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/marker_genes_ranked.png" alt="10x" width=700>
<br>
<small>Marker genes ranked by Scanpy</small>
</p>


## Differential Gene Expression Analysis

Now that the clusters are annotated, we will pseudobulk the data to perform differential gene expression analysis using `DESeq2` on these aggregated counts. We will focus on the fibroblast cluster, aiming to identify genes that are differentially expressed between the `Cancer` and `Normal Adjacent` samples.

To achieve this, we use `Scanpy`’s `aggregate` function to group the `AnnData` object table by the `grouped_clusters` and `sample` in the `AnnData` object metadata. We use the `filtered_counts` layer for gene expression and sum the counts, as `DESeq2` requires raw count data as its input.

`DESeq2` requires three primary inputs: a metadata DataFrame and a counts DataFrame, along with a design formula.

* The metadata DataFrame maps each sample to its experimental conditions or covariates. In the metadata DataFrame, each row is a sample and each column is a condition or covariate.
* The counts DataFrame contains raw gene counts, with each row representing a sample and each column representing a gene.
* The design formula is the model `DESeq2` uses to estimate gene expression changes (log2 fold changes) and assess their statistical significance.

The design formula specifies which factors or covariates from the metadata DataFrame influence gene expression and how they relate to the observed counts. Written in R's formula syntax, it typically begins with a tilde `~` followed by the variables from the metadata table that you intend to include in the model. For this analysis, we use the formula `~tissue`, which will compare the fibroblasts in the `Cancer` and `Normal Adjacent` samples.

For more detailed information on `DESeq2` and setting up designs for various studies, see the official `DESeq2` R [documentation](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start) and Python [documentation](https://pydeseq2.readthedocs.io/en/stable/).


## Discussion

Upregulation of these genes in fibroblasts contributes to colon cancer progression by remodeling the extracellular matrix (ECM), secreting pro-tumorigenic factors, and fostering an immunosuppressive and pro-invasive tumor microenvironment.

### Extracellular Matrix (ECM) Remodeling and Adhesion

These genes are directly involved in building, modifying, and interacting with the ECM.

* **COL1A1, COL8A1, COL12A1:** These genes encode different types of collagen, a major structural protein of the ECM. **Cancer-associated fibroblasts (CAFs)** produce excessive amounts of collagen, leading to a stiff and dense tumor stroma (desmoplasia). This dense matrix physically promotes cancer cell invasion and migration.
* **COMP (Cartilage Oligomeric Matrix Protein):** A non-collagenous protein of the ECM that helps organize other matrix components, contributing to the structural changes in the tumor.
* **ADAM12 (ADAM Metallopeptidase Domain 12):** This gene encodes a protein with both cell adhesion and protease functions. It can cleave ECM proteins, allowing cancer cells to move through the tissue more easily.
* **ITGA11 (Integrin Subunit Alpha 11):** An integrin protein that acts as a cell surface receptor. It helps fibroblasts adhere to and remodel collagen, a critical step in tumor matrix reorganization.
* **ITGBL1 (Integrin Beta Like 1):** This protein modulates cell-ECM adhesion and signaling, influencing cell migration and survival.


### Signaling and Growth Factor Regulation

These genes are involved in signaling pathways that regulate cell proliferation, differentiation, and communication within the tumor microenvironment.

* **SFRP4 (Secreted Frizzled-Related Protein 4):** A modulator of the Wnt signaling pathway, which is often dysregulated in cancer. SFRP4 can promote tumor progression and metastasis.
* **INHBA (Inhibin Beta A):** A subunit of activin A, a growth factor that can activate fibroblasts, leading to the pro-tumorigenic desmoplastic reaction.
* **CTHRC1 (Collagen Triple Helix Repeat Containing 1):** A secreted protein that promotes cell migration and ECM remodeling. It is highly expressed in CAFs and is a known contributor to cancer cell invasion and a marker for poor prognosis.

We can now visualize the spatial expression of these genes.


<p style="margin-top: 15px; margin-bottom: 15px;" align="center">
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample1_col1a1.png" alt="10x" width=150>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample2_col1a1.png" alt="10x" width=150>
<br>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample3_col1a1.png" alt="10x" width=150>
<img src="https://raw.githubusercontent.com/young-sudo/spatialflow/main/img/sample4_col1a1.png" alt="10x" width=150>
<br>
<small>Plots for COL1A1's spatial gene expression in samples</small>
</p>

We observe that `COL1A1` is expressed by fibroblasts situated closer to the tumor in the `Cancer` samples. Though beyond the scope of this guide, to delve deeper into the specific biology of these fibroblasts, the next step would involve subsetting the fibroblast-containing cluster from the overall dataset. This isolated subset can then be re-clustered to further investigate differences and heterogeneity within the fibroblast populations.
