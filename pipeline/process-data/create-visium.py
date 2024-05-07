import scanpy as sc
import spatialdata_io
import spatialdata_plot
import matplotlib.pyplot as plt
import PIL.Image
PIL.Image.MAX_IMAGE_PIXELS = None

vs = spatialdata_io.visium(snakemake.params['input_dir'], dataset_id = snakemake.wildcards['sample'])
vs.table.var_names_make_unique()

vs.table.var['mt'] = vs.table.var_names.str.startswith('MT-')
vs.table.var['hb'] = vs.table.var_names.str.match('HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ')

sc.pp.calculate_qc_metrics(vs.table, qc_vars=['mt', 'hb'], percent_top=None, log1p=False, inplace=True)

#sc.pl.spatial(vs.table, color = ["total_counts", "n_genes_by_counts",'pct_counts_mt', 'pct_counts_hb'], spot_size = 300)

with plt.rc_context():
    fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(15, 4))
    vs.pl.render_shapes(color="total_counts").pl.show(ax=axs[0], title="total counts", coordinate_systems = 'global')
    vs.pl.render_shapes(color="n_genes_by_counts").pl.show(ax=axs[1], title="n_genes_by_counts", coordinate_systems = 'global')
    vs.pl.render_shapes(color="pct_counts_mt").pl.show(ax=axs[2], title="MT", coordinate_systems = 'global')
    vs.pl.render_shapes(color="pct_counts_hb").pl.show(ax=axs[3], title="HB", coordinate_systems = 'global')
    plt.tight_layout()
    plt.savefig(snakemake.output['qc_fig'], bbox_inches="tight")
    plt.close()

vs.write(snakemake.output['zarr'])
