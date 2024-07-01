import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'

import cell2location
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch

torch.set_float32_matmul_precision('medium')

# Snakemake params
tumour = snakemake.wildcards["sample"]
tumour_folder = "data/visium/" + tumour 
nuclei_count_file = snakemake.input["nuclei_counts"]
sc_reference = snakemake.input["h5ad"]
stim_expr = snakemake.input["stim_expr"] 
epochs = snakemake.params["epochs"]
density_plot = snakemake.output["plot"]
abundances = snakemake.output["mat"]
adata_output = snakemake.output["h5ad"]
model_output = snakemake.output["model"]

# Load annData Object, keep tissue only
adata_vis = sc.read_visium(tumour_folder, library_id=tumour)
adata_vis.var_names_make_unique()
adata_vis = adata_vis[adata_vis.obs.in_tissue == 1]

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var.index.values]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

inf_aver = pd.read_csv(stim_expr, index_col=0)

# Make sure to find shared genegenes between visium and infer_aver
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# Get average nuclei count per spot
nuclei_count = pd.read_csv(nuclei_count_file, index_col = 0)
avg_nuclei = nuclei_count['n_nuclei'].mean().round()

# Prepare model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)
mod = cell2location.models.Cell2location(adata=adata_vis, cell_state_df=inf_aver, N_cells_per_location=int(avg_nuclei), detection_alpha=20)
mod.view_anndata_setup()

# Train
mod.train(max_epochs=epochs, batch_size=None, train_size=1, accelerator="gpu")

adata_vis = mod.export_posterior(adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True})

# Plot abundances
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

cell_types = ['B cell', 'CD4-positive, alpha-beta T cell', 'CD8-positive, alpha-beta T cell', 'blood vessel endothelial cell', 'fibroblast', 'macrophage', 'myeloid dendritic cell', 'neural cell', 'pancreatic acinar cell', 'pancreatic ductal cell', 'pancreatic epsilon cell', 'pancreatic stellate cell', 'plasma cell', 'type B pancreatic cell']

with plt.rc_context():
    sc.pl.spatial(adata_vis, cmap='magma', basis = 'spatial', color = cell_types, img_key='lowres', vmin=0, vmax='p99.2', show = False)
    plt.savefig(density_plot, bbox_inches="tight")
    plt.close()

# Save model
mod.save(model_output, overwrite = True)
adata_vis.write(adata_output)
adata_vis.obs.to_csv(abundances)

# Finish Job
exit()
