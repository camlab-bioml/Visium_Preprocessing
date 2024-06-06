import os
import glob
import skimage
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
# input_sample = "1-HT242P1H1-S1Fc1U1Z1B1_1"
# input_folder = "data/visium/" + input_sample
# spatial_folder = input_folder + "/spatial/"
# spot_path = "output/v1/data/spots/" + input_sample

input_sample = snakemake.wildcards['sample'] 
input_folder = snakemake.params['input_dir']
spatial_folder = input_folder + "/spatial/"
spot_path = snakemake.output['spots']

# Grab tiff files from input folder
types = (spatial_folder + '*.tif', spatial_folder + '*.tiff') # the tuple of file types
tiff_file = []
for i in types:
   tiff_file.extend(glob.glob(i))
tiff_file = tiff_file[0]

# Load visium object from rule, filter on tissue
adata = sq.read.visium(input_folder)
adata = adata[adata.obs['in_tissue']==1]

# Load image array, get hires, convert to imagecontainer
img = skimage.io.imread(tiff_file)
img = sq.im.ImageContainer(img, scale=1)

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
    plt.close()

# Generate nuclei count dataframe
nuclei_count = pd.DataFrame(n_nuclei, columns=['n_nuclei'], index = barcode)
nuclei_count.to_csv(path_or_buf= spot_path + "/nuclei_count.csv")

# Finish Job
exit()
