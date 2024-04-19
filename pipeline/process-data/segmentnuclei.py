#!/home/tiakju/miniconda3/envs/rapids-24.02/bin/ python3
import numpy as np
import matplotlib.pyplot as plt
import squidpy as sq
import scanpy as sc
import pandas as pd
from stardist.models import StarDist2D
import tensorflow as tf
from csbdeep.utils import normalize
from PIL import Image

# Remove decompression bomb
Image.MAX_IMAGE_PIXELS = 1000000000

input_adata = "output/v1/" + snakemake.wildcards['sample'] +".h5ad"

# Load anndata object
adata = sc.read_h5ad(input_adata)

# Get hi-res image nd-array from anndata and generate imageContainer
img = adata.uns["spatial"][args.name]['images']['hires']
img = sq.im.ImageContainer(img)

# Generate spot crops
gen = img.generate_spot_crops(adata, scale = 1, squeeze=True, return_obs=True)

# Define Stardist function
def stardist_2D_versatile_he(img, nms_thresh=None, prob_thresh=None):
    # axis_norm = (0,1)   # normalize channels independently
    axis_norm = (0, 1, 2)  # normalize channels jointly
    # Make sure to normalize the input image beforehand or supply a normalizer to the prediction function.
    # this is the default normalizer noted in StarDist examples.
    img = normalize(img, 1, 99.8, axis=axis_norm)
    model = StarDist2D.from_pretrained("2D_versatile_he")
    labels, _ = model.predict_instances(
        img, nms_thresh=nms_thresh, prob_thresh=prob_thresh
    )
    return labels

# Segmentation for each visium spot
n_spot = adata.n_obs
n_nuclei = []
barcode = []

myfun.create_folder(output_folder + "spots")

for i in range(n_spot):
    crop, id = next(gen, i)
    sq.im.segment(img=crop, layer="image", channel=None, method=stardist_2D_versatile_he, layer_added="segmented_stardist_default" ,prob_thresh=0.6, nms_thresh=None)
    n_nuclei.append(len(np.unique(crop['segmented_stardist_default'])))
    barcode.append(id)
    
    fig, axes = plt.subplots(1, 2)
    crop.show("image", ax=axes[0])
    _ = axes[0].set_title("H&E " + id)
    crop.show("segmented_stardist_default", cmap="jet", interpolation="none", ax=axes[1])
    _ = axes[1].set_title("segmentation")
    plt.savefig(output_folder + "spots/" + id + "_segmentation.png")

# Generate nuclei count dataframe
nuclei_count = pd.DataFrame(n_nuclei, columns=['n_nuclei'], index = barcode)
nuclei_count.to_csv(path_or_buf=output_files[0]['output_file_0'])

# Join anndata obs, write to adata
adata.obs = adata.obs.join(nuclei_count).fillna(0)
adata.write_h5ad(output_files[1]['output_file_1'])

# Finish Job
exit()
