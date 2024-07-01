import os
import anndata as ad

os.getcwd()

adata1 = ad.read_h5ad("data/singlecell/sibyl_acinarTypeB_rawcounts.h5ad")
adata2 = ad.read_h5ad("data/singlecell/sibyl_endothelial_large1_clean_rawcounts.h5ad")
adata3 = ad.read_h5ad("data/singlecell/sibyl_epithelial_rawcounts.h5ad")
adata4 = ad.read_h5ad("data/singlecell/sibyl_hemato_rawcounts.h5ad")
adata5 = ad.read_h5ad("data/singlecell/sibyl_mesenchyme_large1_cleaned_rawcounts.h5ad")

to_merge = {'acinarTypeB': adata1, 'endothelial': adata2, 'epithelial': adata3, 'hemato': adata4, 'mesenchyme': adata5}
merged = ad.concat(to_merge, label = "dataset_name", join="outer")

ad.AnnData.write_h5ad(merged, filename="data/singlecell/sibyl_merged.h5ad")

exit()