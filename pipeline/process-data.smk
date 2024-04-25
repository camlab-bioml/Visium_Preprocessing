
process_data = {
    'processed_violin': expand(output + 'figures/qc/Filtered-{sample}-spatial-qc.png', sample = sample_ids)
}

# Conditionally add output files based on user specifications in config
if config['plot_MT_HB']:
    process_data.append(output + 'figures/qc/MT-HB-counts-violin.png')

rule create_spatialdata:
    input:
        feature_matrix = "data/output/{sample}_processed/filtered_feature_bc_matrix.h5",
        barcodes = "data/output/{sample}_processed/barcodes.tsv.gz",
        features = "data/output/{sample}_processed/features.tsv.gz",
        matrix = "data/output/{sample}_processed/matrix.mtx.gz",
        scale_factors = "data/output/{sample}_processed/spatial/scalefactors_json.json",
        hires = "data/output/{sample}_processed/spatial/tissue_hires_image.png",
        lowres = "data/output/{sample}_processed/spatial/tissue_lowres_image.png",
        tissue_pos = "data/output/{sample}_processed/spatial/tissue_positions_list.csv"
    params:
        input_dir = directory("data/output/{sample}_processed/")
    output:
        qc_fig = report(output + 'figures/qc/Raw-{sample}-spatial-qc.png', category = 'Step 1: Data QC'),
        zarr = directory(output + 'data/raw/{sample}/')
    script:
        "process-data/create-visium.py"

rule combine_anndata_and_plot_counts_violin:
    input:
        spatial_data_out = expand(output + 'figures/qc/Raw-{sample}-spatial-qc.png', sample = sample_ids)
    params:
        spatial_datas = expand(output + 'data/raw/{sample}/', sample = sample_ids)
    output:
        counts_violin = report(output + 'figures/qc/gene-counts-violin.png', category = 'Step 1: Data QC'),
        zarr = directory(output + 'data/raw/combined/'),
        finish_message = temp(output + 'messages/finish-combining-visium.txt')
    script:
        "process-data/visium-concatenate.py"

# rule plot_MT_and_HB_counts: # Needs to be tested on different dataset to ensure it works
#     input:
#         spatial_concat_out = output + 'figures/qc/gene-counts-violin.png',
#     params:
#         spatialdata_dir = output + 'data/raw/'
#     output:
#         mt_hb_violin = output + 'figures/qc/MT-HB-counts-violin.png'
#     script:
#         "process-data/plot-MT-HB-counts.py"

rule visium_filter:
    input:
        spatial_data_out = output + 'messages/finish-combining-visium.txt'
    params:
        spatialdata_dir = output + 'data/raw/combined/',
        pct_counts_hb_thr = config['filtering']['max_pct_hb'],
        pct_counts_mt_thr = config['filtering']['max_pct_mt'],
        n_genes_by_counts_thr = config['filtering']['min_n_genes_by_counts']
    output:
        filtered_spatialdata = directory(output + 'data/clean/'),
        finish_message = temp(output + 'messages/finish-visium-filter.txt')
    script:
        "process-data/visium-filter.py"

rule plot_counts_filtered:
    input:
        message = output + 'messages/finish-visium-filter.txt'
    params:
        filtered_visium_spatialdata = output + 'data/clean/'
    output:
        qc_fig = report(output + 'figures/qc/Filtered-{sample}-spatial-qc.png', category = 'Step 1: Data QC'),
    script:
        'process-data/plot-visium-qc.py'
