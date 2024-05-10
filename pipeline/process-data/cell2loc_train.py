import numpy as np
import scanpy as sc
import pandas as pd
import cell2location

# Set up expected inputsls
sc_input = "data/singlecell/scRNASeq-SingleR-annotated-sce-Peng.h5ad"# snakemake.input['h5ad']

# Set up expected outputs
model_output = "output/v1/Peng_cell2loc_model/" # snakemake.output['dir']
nb_output =  "output/v1/Peng_cell2loc_model/stimulated_expression.csv"  # snakemake.output['mat']
adata_output = "output/v1/Peng_cell2loc_model/model.pt" # snakemake.output['model']

# Load annData Object, set to ensembl id
adata_ref =  sc.read_h5ad(sc_input)

# Select genes
selected = cell2location.utils.filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

# Set up models, train, save
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,batch_key='sample', labels_key='singler.label')
mod = cell2location.models.RegressionModel(adata_ref)
mod.view_anndata_setup()
mod.train(max_epochs=100, accelerator = "gpu", batch_size=2500, train_size=1) # change to 2 epoch
mod.save(model_output, overwrite=True)

# Calculate posteior
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}) # change to 1000 samples for full run
adata_ref.write(adata_output)

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
