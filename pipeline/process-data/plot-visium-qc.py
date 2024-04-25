import scanpy as sc
import spatialdata as sp
import spatialdata_io
import spatialdata_plot as pl

import matplotlib.pyplot as plt

spatialdata = sp.SpatialData.read(snakemake.params['filtered_visium_spatialdata'])
sample_table = spatialdata.table[spatialdata.table.obs.region == snakemake.wildcards['sample'],:]

with plt.rc_context():
    fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(15, 4))
    sc.pl.spatial(sample_table, color = ["total_counts", "n_genes_by_counts",'pct_counts_mt', 'pct_counts_hb'], spot_size = 250)
    plt.tight_layout()
    plt.savefig(snakemake.output['qc_fig'], bbox_inches="tight")
    plt.close()