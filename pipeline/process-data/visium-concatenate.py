import scanpy as sc
import spatialdata as sp
import spatialdata_io

import matplotlib.pyplot as plt

# Read in and concatenate the spatial data objects
spatialdatas_list = []
for sample_path in snakemake.params['spatial_datas']:
    spatialdatas_list.append(sp.SpatialData.read(sample_path))

spatialdata = sp.concatenate(spatialdatas_list)
spatialdata.table.obs_names_make_unique()

# Create violin plot
with plt.rc_context():
    sc.pl.violin(spatialdata.table, ['n_genes_by_counts', 'total_counts'],
                jitter=0.4, rotation = 90, groupby = 'region', multi_panel=True)
    plt.savefig(snakemake.output['counts_violin'], bbox_inches='tight')
    plt.close()

# Save concatenated 
spatialdata.write(snakemake.output['zarr'])

# Required for subsequent rules to have input
f = open(snakemake.output['finish_message'],"w+")
f.write("spatialdata filtered")