#!/usr/bin/env python3

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import gc
import scanpy.external as sce
import spatialdata as spd
import utils

def main():

    adata = sc.read_h5ad("filtered_sdata")

    # Elbow plot
    sc.pl.pca_variance_ratio(adata, log=True,n_pcs=50)
    plt.savefig("figures/pca_variance_ratio.png", dpi=300, bbox_inches="tight")
    plt.close()



    # neighborhood and clustering resolution
    RES = 0.5 # clustering resolution
    NEIGHBORS = 30  # number of neighbors

    MIN_DIST=0.5 #default 0.5
    SPREAD=2 #default 1

    sc.pp.neighbors(adata, n_neighbors=NEIGHBORS, use_rep="X_pca",metric="correlation")
    sc.tl.leiden(adata, flavor="igraph", key_added="clusters", resolution=RES,random_state=0)

    # To ensure that the results are reproducible we are going to reorder the clusters by size.
    adata.obs['orig_clusters'] = adata.obs['clusters']

    clusters = adata.obs['clusters'].astype(int)

    # Count cells per cluster
    cluster_sizes = clusters.value_counts().sort_values(ascending=False)

    # Create mapping: old cluster ID → new ordered ID
    cluster_order = {old: new for new, old in enumerate(cluster_sizes.index)}

    # Relabel clusters in adata
    adata.obs['clusters'] = clusters.map(cluster_order).astype(str)

    # Set random_state for reproducible UMAP
    sc.tl.umap(adata,min_dist=MIN_DIST, spread=SPREAD, random_state=0)

    # Plot UMAP
    sc.pl.umap(adata, color=["clusters"], title="UMAP by Clusters")
    plt.tight_layout()
    plt.savefig("figures/umap_by_clusters.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(adata, color=["sample"], title="UMAP by Sample")
    plt.tight_layout()
    plt.savefig("figures/umap_by_sample.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Cell distribution across clusters
    sample_names = adata.obs["sample"].unique()
    plt.imshow(pd.crosstab(adata.obs["sample"], adata.obs["clusters"]), cmap='hot', interpolation='nearest')
    plt.title("Cell Distribution Across Clusters")
    plt.xlabel("Cluster")
    plt.yticks(range(len(sample_names)), sample_names)

    plt.tight_layout()
    plt.savefig("figures/cell_dist_clusters.png", dpi=300, bbox_inches="tight")
    # plt.show()
    plt.close()

    del RES, NEIGHBORS, MIN_DIST, SPREAD
    gc.collect()



    # Correct for Batch Effects
    # Current dataset doesn't exhibit batch effects

    # neighborhood and clustering resolution
    RES = 0.5 # clustering resolution
    NEIGHBORS = 30  # number of neighbors

    MIN_DIST=0.5 #default 0.5
    SPREAD=2 #default 1

    concatenated_sdata = spd.read_zarr("concatenated_sdata")

    # Performing batch correction
    adata_harmony = concatenated_sdata["segmentation_counts"].copy()
    sce.pp.harmony_integrate(adata_harmony, key="sample", basis="X_pca",max_iter_harmony=20)

    # Copying the harmony PCA embedding results
    adata_harmony.obsm["X_pca_orig"] = adata_harmony.obsm["X_pca"]
    adata_harmony.obsm["X_pca"] = adata_harmony.obsm["X_pca_harmony"]
    adata_harmony.obs["cluster_orig"] =  adata_harmony.obs["clusters"]

    sc.pp.neighbors(adata_harmony, n_neighbors=NEIGHBORS, use_rep="X_pca",metric="correlation")
    sc.tl.leiden(adata_harmony, flavor="igraph",key_added="harmony_clusters", resolution=RES,random_state=0)
    adata_harmony.obs["clusters"] =  adata_harmony.obs["harmony_clusters"]

    # Set random_state for reproducible UMAP
    sc.tl.umap(adata_harmony,min_dist=MIN_DIST, spread=SPREAD, random_state=0)

    # Plot UMAP
    sc.pl.umap(adata_harmony, color=["clusters"], title="Harmony Corrected UMAP by Clusters")
    sc.pl.umap(adata_harmony, color=["sample"], title="Harmony Corrected UMAP by Sample")

    # Cell distribution across clusters
    sample_names = adata_harmony.obs["sample"].unique()
    plt.imshow(pd.crosstab(adata_harmony.obs["sample"], adata_harmony.obs["clusters"]), cmap='hot', interpolation='nearest')
    plt.title("Cell Distribution Across Clusters")
    plt.xlabel("Cluster")
    plt.yticks(range(len(sample_names)), sample_names)
    plt.show()
    plt.close('all')


    # If you want to use the harmony results in the analysis then overwrite
    # the AnnData table in the SpatialData object by removing the comment below.
    # concatenated_sdata["segmentation_counts"] = adata_harmony

    del adata_harmony, RES, NEIGHBORS, MIN_DIST, SPREAD
    gc.collect()


    # =================================
    # Spatial visualization of clusters
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



if __name__ == "__main__":
    main()

