#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import spatialdata as spd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import utils
import matplotlib.pyplot as plt
import numpy as np


def main():
    concatenated_sdata = spd.read_zarr("concatenated_sdata")

    # aggregate by the cluster and the sample id.
    aggregated = sc.get.aggregate(concatenated_sdata["segmentation_counts"], by=["grouped_clusters","sample"], func=["sum"],layer="filtered_counts")

    #  We will focus on one cluster of interest
    aggregated_cluster_of_interest = aggregated[aggregated.obs["grouped_clusters"]== 'Fibroblast']


    # Creating metadata table
    metadata_df = aggregated_cluster_of_interest.obs["sample"].str.split("_", expand=True)
    metadata_df.columns = ["tissue", "patient"]

    #Creating counts table
    counts_df = pd.DataFrame(aggregated_cluster_of_interest.layers["sum"],index=metadata_df.index,columns=aggregated_cluster_of_interest.var_names)
    counts_df = counts_df.astype(int)


    # Determining differentially expressed genes from the aggregated data using DeSeq2
    dds = DeseqDataSet(
        counts = counts_df,
        metadata = metadata_df,
        design = "~tissue",
        refit_cooks=True
    )

    dds.deseq2()
    ds = DeseqStats(dds, contrast=["tissue", "Cancer", "Normal"])
    ds.summary()

    print(ds.results_df.sort_values(by='log2FoldChange', ascending=False).head(10))
    ds.results_df.to_csv('aggregated_diffexp.csv')

    # --- Filtering and printing the table ---
    positive_fold_change_genes = ds.results_df[ds.results_df['log2FoldChange'] > 0]

    top_10_positive_genes_table = positive_fold_change_genes.sort_values(by='padj', ascending=True).head(10)

    print("\n--- Top 10 Positively Differentially Expressed Genes (Cancer vs. Normal) ---")
    print(top_10_positive_genes_table)


    # 
    image_elements = list(concatenated_sdata.images.keys())
    shape_elements = list(concatenated_sdata.shapes.keys())

    # We are going to create a bounding box to crop the data to the capture area.
    extents = []

    for i in range(len(image_elements)):
        extent =  spd.get_extent(concatenated_sdata,elements=[shape_elements[i]],coordinate_system='downscale_to_hires')
        extents.append(extent)

    # Plotting
    if len(image_elements) != len(shape_elements):
        print("Check the spatial data to make sure that for every image there is a shape")
    else:
        for i in range(len(image_elements)):
            print("Plotting: "+ image_elements[i])
            title=image_elements[i].replace("_hires_tissue_image","")
            utils.crop0(concatenated_sdata,crs="downscale_to_hires",bbox=extents[i]).pl.render_images(image_elements[i]).pl.render_shapes(shape_elements[i],color="clusters").pl.show(coordinate_systems="downscale_to_hires", title=title)



    gene_name = "COL1A1"

    for i in range(len(image_elements)):
        print("Plotting: "+ image_elements[i])
        title=image_elements[i].replace("_hires_tissue_image","")
        utils.crop0(concatenated_sdata,crs="downscale_to_hires",bbox=extents[i]).pl.render_images(image_elements[i]).pl.render_shapes(shape_elements[i],color=gene_name).pl.show(coordinate_systems="downscale_to_hires",title=title)



    # Work directly on ds.results_df
    res = ds.results_df#.copy()

    # Remove rows with missing padj or log2FoldChange
    res = res.dropna(subset=["padj", "log2FoldChange"])

    # Compute -log10(padj)
    res["neg_log10_padj"] = -np.log10(res["padj"])

    padj_thresh = 0.05
    logfc_thresh = 1

    res["significant"] = (
        (res["padj"] < padj_thresh) &
        (abs(res["log2FoldChange"]) > logfc_thresh)
    )

    plt.figure(figsize=(8,6))

    plt.scatter(
        res["log2FoldChange"],
        res["neg_log10_padj"],
        s=10,
        alpha=0.7,
    )

    plt.scatter(
        res.loc[res["significant"], "log2FoldChange"],
        res.loc[res["significant"], "neg_log10_padj"],
        s=12,
        alpha=0.9,
    )

    plt.axhline(-np.log10(padj_thresh), linestyle="--")
    plt.axvline(logfc_thresh, linestyle="--")
    plt.axvline(-logfc_thresh, linestyle="--")

    plt.xlabel("log₂ fold change")
    plt.ylabel("–log₁₀ adjusted p-value")
    plt.title("Volcano plot: Cancer vs Normal")

    plt.tight_layout()
    plt.savefig("figures/volcano_plot.png", dpi=300)
    plt.close()



if __name__ == "__main__":
    main()
