import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'

import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from types import MethodType
from scvi import REGISTRY_KEYS
from pyro import clear_param_store
import cell2location
from cell2location.models.base._pyro_mixin import PltExportMixin
import torch

torch.set_float32_matmul_precision('medium')

# Uncomment for troubleshooting
# sc_input = "data/singlecell/scRNASeq-SingleR-annotated-sce-Peng.h5ad"
# epochs = 2
# adata_output =  "output/v1/cell2loc/Peng_cell2loc_model/model_data.h5ad"
# nb_output =  "output/v1/cell2loc/Peng_cell2loc_model/stimulated_expression.csv"
# model_output = "output/v1/cell2loc/Peng_cell2loc_model/model.pt"
# train_history = "output/v1/cell2loc/Peng_cell2loc_model/train_history.png"
# train_accuracy_1 = "output/v1/cell2loc/Peng_cell2loc_model/train_accuracy_1.png"
# train_accuracy_2 = "output/v1/cell2loc/Peng_cell2loc_model/train_accuracy_2.png"

# Snakemake Params
sc_input = snakemake.input['h5ad']
epochs = snakemake.params["epochs"]
adata_output = snakemake.output['h5ad']
nb_output = snakemake.output['mat']
model_output = snakemake.output['model']
train_history = snakemake.output["history"]
train_accuracy_1 = snakemake.output["accuracy1"]
train_accuracy_2 = snakemake.output["accuracy2"]

# Solution for stacking QC plots
def plotQC_1(
    self,
    summary_name: str = "means",
    use_n_obs: int = 1000
):
    PltExportMixin.plot_QC(self, summary_name=summary_name, use_n_obs=use_n_obs)

def plotQC_2(
    self, 
    summary_name: str = "means", 
    use_n_obs: int = 1000,
    scale_average_detection: bool = True,
):
    inf_aver = self.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
    if scale_average_detection and ("detection_y_c" in list(self.samples[f"post_sample_{summary_name}"].keys())):
        inf_aver = inf_aver * self.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()
    aver = self._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
    aver = aver[self.factor_names_]

    plt.hist2d(
        np.log10(aver.values.flatten() + 1),
        np.log10(inf_aver.flatten() + 1),
        bins=50,
        norm=matplotlib.colors.LogNorm(),
    )
    plt.xlabel("Mean expression for every gene in every cluster")
    plt.ylabel("Estimated expression for every gene in every cluster")
 
# Load annData Object, set to ensembl id
adata_ref = sc.read_h5ad(sc_input)

# Select genes
selected = cell2location.utils.filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

# Set up models, train, save
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, batch_key='sample', labels_key='subset')
mod = cell2location.models.RegressionModel(adata_ref)
mod.view_anndata_setup()
mod.train(max_epochs=epochs, batch_size=2500, train_size=1) # change to 2 epoch
mod.save(model_output, overwrite=True)

# Add the new method to the instance
mod.plotQC_1 = MethodType(plotQC_1, mod)
mod.plotQC_2 = MethodType(plotQC_2, mod)

# Plots
mod.plot_history(20)
plt.savefig(train_history, bbox_inches='tight')
plt.close()

# Calculate posteior
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}) # change to 1000 samples for full run
adata_ref.write(adata_output)

# Accuracy metrics
mod.plotQC_1()
plt.savefig(train_accuracy_1, bbox_inches='tight')
plt.close()

mod.plotQC_2()
plt.savefig(train_accuracy_2, bbox_inches='tight')
plt.close()

# Export proportion matrix as pd 
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
    
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.to_csv(nb_output)

# Finish Job
exit()
