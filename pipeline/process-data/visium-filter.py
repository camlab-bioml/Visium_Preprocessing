import scanpy as sc
import spatialdata as sp
import spatialdata_io
import spatialdata_plot as pl

import matplotlib.pyplot as plt

spatialdata = sp.SpatialData.read(snakemake.params['spatialdata_dir'])

## Create true/false lists to determine if spots should be removed
## if no threshold is set in config all points are kept
# Hemoglobin threshold
if snakemake.params['pct_counts_hb_thr'] != 'None':
    hb_thr = spatialdata.table.obs['pct_counts_hb'] < snakemake.params['pct_counts_hb_thr']
    hb_thr = hb_thr.to_list()
else:
    hb_thr = [True for x in range(spatialdata.table.shape[0])]

# Mitochondrial threshold
if snakemake.params['pct_counts_mt_thr'] != 'None':
    mt_thr = spatialdata.table.obs['pct_counts_mt'] < snakemake.params['pct_counts_mt_thr']
    mt_thr = mt_thr.to_list()
else:
    mt_thr = [True for x in range(spatialdata.table.shape[0])]

# Number of genes by counts
if snakemake.params['n_genes_by_counts_thr'] != 'None':
    n_genes_by_counts = spatialdata.table.obs['n_genes_by_counts'] > snakemake.params['n_genes_by_counts_thr']
    n_genes_by_counts = n_genes_by_counts.to_list()
else:
    n_genes_by_counts = [True for x in range(spatialdata.table.shape[0])]

# Filter down spatialdata object
keep = [all(elements) for elements in zip(hb_thr, mt_thr, n_genes_by_counts)]
clean_spots_table = spatialdata.table[keep,:].copy()
del spatialdata.table
spatialdata.table = clean_spots_table

# Save output
spatialdata.write(snakemake.output['filtered_spatialdata'], overwrite=True)

# Required for subsequent rules to have input
f = open(snakemake.output['finish_message'],"w+")
f.write("spatialdata filtered")