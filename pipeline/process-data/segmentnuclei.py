import os
import pandas as pd
import numpy as np
import squidpy as sq
import matplotlib.pyplot as plt
import tensorflow as tf
from stardist.models import StarDist2D
from csbdeep.utils import normalize
import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

# **** Because of compatibility with squidpy I will load the raw adata. Recommend to do post-hoc addition of nuclei count to clean SpatialData object. 
# To update script for GPU implementation -> numpy problem

# For Troubleshooting
# input_sample = "donorA"
# input_folder = "data/visium/" + input_sample

input_sample = snakemake.wildcards['sample'] # replace with snakemake object
input_folder = snakemake.params['input_dir']

# Load visium object from rule, filter on tissue
adata = sq.read.visium(input_folder)
adata = adata[adata.obs['in_tissue']==1]

# Load image array, get hires, convert to imagecontainer
img = adata.uns["spatial"]
scale_factor = img[next(iter(img))]['scalefactors']["tissue_hires_scalef"]
img = img[next(iter(img))]['images']['hires']
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

spot_path = snakemake.output['spots']
if not os.path.exists(spot_path):
    os.makedirs(spot_path)
    
for i in range(n_spot):
    crop, id = next(gen, i)
    sq.im.segment(img=crop, layer="image", channel=None, method=stardist_2D_versatile_he, layer_added="segmented_stardist_default" , prob_thresh=0.6, nms_thresh=None)
    n_nuclei.append(len(np.unique(crop['segmented_stardist_default'])))
    barcode.append(id)
    fig, axes = plt.subplots(1, 2)
    crop.show("image", ax=axes[0])
    _ = axes[0].set_title("H&E " + id)
    crop.show("segmented_stardist_default", cmap="jet", interpolation="none", ax=axes[1])
    _ = axes[1].set_title("segmentation")
    plt.savefig(spot_path + "/" + id + "_segmentation.png")

# Generate nuclei count dataframe
nuclei_count = pd.DataFrame(n_nuclei, columns=['n_nuclei'], index = barcode)
nuclei_count.to_csv(path_or_buf= spot_path + "/nuclei_count.csv")

# Finish Job
exit()
