#!/usr/bin/env python3

import utils
import gc

def main():
    # Create and save Zarr files for the cell segmentation outputs.
    samples = {
        "Colon_Cancer_P1":["data/Cancer_P1_filtered_feature_cell_matrix.h5",
                        "data/Cancer_P1_tissue_hires_image.png",
                        "data/Cancer_P1_scalefactors_json.json",
                        "data/Cancer_P1_cell_segmentations.geojson",
                        "Colon_Cancer_P1"],
        "Colon_Cancer_P2":["data/Cancer_P2_filtered_feature_cell_matrix.h5",
                        "data/Cancer_P2_tissue_hires_image.png",
                        "data/Cancer_P2_scalefactors_json.json",
                        "data/Cancer_P2_cell_segmentations.geojson",
                        "Colon_Cancer_P2"],
        "Colon_Normal_P3":["data/Norm_P3_filtered_feature_cell_matrix.h5",
                        "data/Norm_P3_tissue_hires_image.png",
                        "data/Norm_P3_scalefactors_json.json",
                        "data/Norm_P3_cell_segmentations.geojson",
                        "Colon_Normal_P3"],
        "Colon_Normal_P5":["data/Norm_P5_filtered_feature_cell_matrix.h5",
                        "data/Norm_P5_tissue_hires_image.png",
                        "data/Norm_P5_scalefactors_json.json",
                        "data/Norm_P5_cell_segmentations.geojson",
                        "Colon_Normal_P5"],
    }
    print("Saving zarr files")
    for key, inputs in samples.items():
        utils.create_zarr(count_matrix_path=inputs[0],
                    image_path=inputs[1],
                    scale_factors_path=inputs[2],
                    geojson_path=inputs[3],
                    sample_name=inputs[4])

    del samples, inputs, key
    gc.collect()

if __name__ == "__main__":
    main()
