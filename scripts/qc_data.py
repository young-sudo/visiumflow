#!/usr/bin/env python3

import utils
import os
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import gc

def qc():
    concatenated_sdata = utils.prepare_data()

    adata = concatenated_sdata["segmentation_counts"]
    # we link the AnnData Table in the SpatialData object
    # to the variable adata to make the code easier to read

    # remove low-quality cell segmentation bins
    # by visualizing total UMI distribution to estimate suitable
    # for empty or sparsely populated bins

    # Add mitochondrial gene calculation for QC
    # not used for QC here
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, percent_top=None)

    # Directory for output figures
    os.makedirs("figures", exist_ok=True)

    # 1 row, 3 columns
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))  # adjust size as needed

    # Total UMI
    sc.pl.violin(
        adata=adata,
        keys=["log1p_total_counts"],
        groupby="sample",
        stripplot=False,
        inner="box",
        show=False,
        ax=axes[0]
    )
    axes[0].set_title("Total UMI by Sample")
    axes[0].axhline(y=4, color='r', linestyle='-')
    axes[0].axhline(y=8, color='r', linestyle='-')

    # Total genes
    sc.pl.violin(
        adata=adata,
        keys=["log1p_n_genes_by_counts"],
        groupby="sample",
        stripplot=False,
        inner="box",
        show=False,
        ax=axes[1]
    )
    axes[1].set_title("Total Genes by Sample")

    # Mitochondrial genes
    sc.pl.violin(
        adata=adata,
        keys=["log1p_total_counts_mt"],
        groupby="sample",
        stripplot=False,
        inner="box",
        show=False,
        ax=axes[2]
    )
    axes[2].set_title("Mitochondrial Genes by Sample")

    # save
    plt.tight_layout()
    plt.savefig("figures/qc.png", dpi=300, bbox_inches="tight")
    plt.close()


    # cell segmentation bins with fewer than 53 counts
    # (corresponding to a log1p value of 4) and more than
    # 2,979 counts (corresponding to a log1p value of 8) are removed
            
    # Estimating the cut off
    min_counts = np.expm1(4).astype("int")
    max_counts = np.expm1(8).astype("int")

    sc.pp.filter_genes(adata, min_cells=50)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, max_counts=max_counts)

    # Visualization for QC
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    sc.pl.violin(adata, keys=["log1p_total_counts"], groupby="sample",
                stripplot=False, inner="box", show=False, ax=axes[0])
    axes[0].set_title("Total UMI by Sample")
    axes[0].axhline(4, color='r')
    axes[0].axhline(8, color='r')

    sc.pl.violin(adata, keys=["log1p_n_genes_by_counts"], groupby="sample",
                stripplot=False, inner="box", show=False, ax=axes[1])
    axes[1].set_title("Total Genes by Sample")

    sc.pl.violin(adata, keys=["log1p_total_counts_mt"], groupby="sample",
                stripplot=False, inner="box", show=False, ax=axes[2])
    axes[2].set_title("Mitochondrial Genes by Sample")

    plt.tight_layout()
    plt.savefig("figures/qc_filt.png", dpi=300, bbox_inches="tight")
    plt.close()

    # storing filtered counts
    adata.layers["filtered_counts"] = adata.X.copy()

    del max_counts, min_counts
    gc.collect()


    # Normalize and save
    sc.pp.normalize_total(adata, target_sum = None)
    sc.pp.log1p(adata)
    sc.tl.pca(adata)

    adata.write("preprocessed_adata.h5ad")


def main():
    qc()

if __name__ == "__main__":
    main()
