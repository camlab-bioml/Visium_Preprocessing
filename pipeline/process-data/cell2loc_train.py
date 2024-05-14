import numpy as np
import scanpy as sc
import pandas as pd
import cell2location
import matplotlib.pyplot as plt
import matplotlib as mpl

# Set up expected inputsls
sc_input = snakemake.input['h5ad'] # "data/singlecell/scRNASeq-SingleR-annotated-sce-Peng.h5ad" #

# Params
epochs = snakemake.params["epochs"]

# Set up expected outputs
adata_output = snakemake.output['h5ad'] #  "output/v1/cell2loc/Peng_cell2loc_model/model_data.h5ad"
nb_output = snakemake.output['mat'] # "output/v1/cell2loc/Peng_cell2loc_model/stimulated_expression.csv"
model_output = snakemake.output['model'] # "output/v1/cell2loc/Peng_cell2loc_model/model.pt"
train_history = snakemake.output["history"] # "output/v1/cell2loc/Peng_cell2loc_model/stimulated_expression.csv"
train_accuracy = snakemake.output["accuracy"] # "output/v1/cell2loc/Peng_cell2loc_model/train_accuracy.png"
 

# Load annData Object, set to ensembl id
adata_ref = sc.read_h5ad(sc_input)

# Select genes
selected = cell2location.utils.filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

# Set up models, train, save
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,batch_key='sample', labels_key='singler.label')
mod = cell2location.models.RegressionModel(adata_ref)
mod.view_anndata_setup()
mod.train(max_epochs=epochs, accelerator = "gpu", batch_size=2500, train_size=1) # change to 2 epoch
mod.save(model_output, overwrite=True)

# Plots
mod.plot_history()
plt.savefig(train_history)
plt.close()

# Calculate posteior
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}) # change to 1000 samples for full run
adata_ref.write(adata_output)

# Accuracy
mod.plot_QC()
plt.savefig(train_accuracy, bbox_inches='tight')
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
