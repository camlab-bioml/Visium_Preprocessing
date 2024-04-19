

process_data = {
    'h5ads': directory(expand(output + 'data/{sample}/', sample = sample_ids))
}

rule create_anndata:
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
        qc_fig = output + 'figures/qc/{sample}-spatial-qc.pdf',
        h5ad = directory(output + 'data/{sample}/')
    script:
        "process-data/create-visium.py"
