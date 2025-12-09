#!/usr/bin/env python3


import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import gc
import scanpy.external as sce
import spatialdata as spd
import utils


def main():
    concatenated_sdata = spd.read_zarr("concatenated_sdata")


    # cannonical markers for annotation
    marker_genes = {
        "Fibroblasts": ["COL1A1", "MMP2", 'VIM'],
        "Goblet cells": ["FCGBP", "MUC2", "CLCA1"],
        "Enterocyte":["EPCAM","KRT8","KRT20","FABP2"],
        "Plasma B cells":['JCHAIN','IGKC','IGHG1'],
        "B cell":['MS4A1','CD74','CD19','CD22'],
        "Smooth muscle":['TAGLN','DES','MYH11'],
        "Tumor":["CEACAM6",'REG1B',"REG1A"]
    }

    # Plot dotplot for initial cluster assessment
    sc.pl.dotplot(adata = concatenated_sdata["segmentation_counts"], var_names = marker_genes, groupby="clusters", standard_scale="var")
    plt.tight_layout()
    plt.savefig("figures/marker_genes_canon.png", dpi=300, bbox_inches="tight")



    # Obtain cluster-specific marker genes
    sc.tl.rank_genes_groups(adata = concatenated_sdata["segmentation_counts"], groupby="clusters", method="wilcoxon")
    sc.pl.rank_genes_groups_dotplot(adata = concatenated_sdata["segmentation_counts"], groupby="clusters", standard_scale="var", n_genes=5)
    plt.tight_layout()
    plt.savefig("figures/marker_genes_ranked.png", dpi=300, bbox_inches="tight")

    df_marker_genes = sc.get.rank_genes_groups_df(adata = concatenated_sdata["segmentation_counts"],group = None,pval_cutoff=0.05)
    df_marker_genes.to_csv("results/marker_genes_pval.csv")




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


    # Cluster annotation mapping
    original_clusters = concatenated_sdata["segmentation_counts"].obs['clusters']

    cell_annotation = {
        '0': 'Fibroblast',
        '1': 'Tumor',
        '2': 'Plasma cell',
        '3': 'Smooth muscle',
        '4': 'Tumor',
        '5': 'B cell and T cell',
        '6': 'Enterocyte',
        '7': 'Plasma cell',
        '8': 'Goblet',
        '9': 'Enterocyte'
    }

    # Apply the mapping. This new_categories Series should have the same index as original_clusters.
    new_categories = original_clusters.astype('string').map(cell_annotation)

    # Assign to sdata_concatenate.
    concatenated_sdata["segmentation_counts"].obs["grouped_clusters"] = new_categories.astype('category')

    # Plotting with new grouped clusters
    for i in range(len(image_elements)):
    print("Plotting: "+ image_elements[i])
    title=image_elements[i].replace("_hires_tissue_image","")
    utils.crop0(concatenated_sdata,crs="downscale_to_hires",bbox=extents[i]).pl.render_images(image_elements[i]).pl.render_shapes(shape_elements[i],color="grouped_clusters").pl.show(coordinate_systems="downscale_to_hires", title=title)



    # Close up on morphology in tissue to confirm
    spd.bounding_box_query(
            concatenated_sdata,
            min_coordinate=[3235, 1500],
            max_coordinate=[4000, 1759],
            axes=("x", "y"),
            target_coordinate_system='downscale_to_hires').pl.render_images("Colon_Cancer_P2_hires_tissue_image").pl.show(coordinate_systems='downscale_to_hires', title="Colon_Cancer_P2")

    spd.bounding_box_query(
            concatenated_sdata,
            min_coordinate=[3235, 1500],
            max_coordinate=[4000, 1759],
            axes=("x", "y"),
            target_coordinate_system='downscale_to_hires').pl.render_images("Colon_Cancer_P2_hires_tissue_image").pl.render_shapes(shape_elements[1],color="grouped_clusters").pl.show(coordinate_systems='downscale_to_hires',title="Colon_Cancer_P2")



    # Save cluster results to csv
    for sample_name in concatenated_sdata["segmentation_counts"].obs['sample'].unique():
        # Filter the table for the current sample
        adata_sample = concatenated_sdata["segmentation_counts"][concatenated_sdata["segmentation_counts"].obs['sample'] == sample_name].copy()

        # Create a DataFrame
        df_output = pd.DataFrame({
            'Barcode': 'cellid_' + adata_sample.obs.index.str.split('cellid_').str[1],
            'Grouped_Annotation': adata_sample.obs['grouped_clusters']
        })

        # Save the results
        output_filename = f"{sample_name}_cell_clusters.csv"
        df_output.to_csv(output_filename, index=False)

        print(f"Saved {output_filename}")

    del adata_sample, df_output, sample_name, output_filename
    gc.collect()



if __name__ == "__main__":
    main()
