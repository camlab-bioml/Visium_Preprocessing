import scanpy as sc
import spatialdata as sp
import spatialdata_io

import matplotlib.pyplot as plt

spatialdata = sp.SpatialData.read(snakemake.params['spatialdata_dir'])

# Create violin plot
with plt.rc_context():
    sc.pl.violin(spatialdata.table, ['pct_counts_mt', 'pct_counts_hb'],
                jitter=0.4, rotation = 90, groupby = 'region', multi_panel=True)
    plt.savefig(snakemake.output['mt_hb_violin'], bbox_inches='tight')
    plt.close()