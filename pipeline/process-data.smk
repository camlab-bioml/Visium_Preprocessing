process_data = {
    "processed_violin": expand(
        output + "figures/qc/Filtered-{sample}-spatial-qc.png", sample=sample_ids
    )
}

spots = {
    "spot_output": expand(output + "data/spots/{sample}/nuclei_count.csv", sample=sample_ids)
}

cell2loc_train = {
    "model": expand(output + "data/{sce}_cell2loc_model/model_adata.h5ad", sce = sce),
    "pt": expand(output + "data/{sce}_cell2loc_model", sce = sce),
    "mat": expand(output + "data/{sce}_cell2loc_model/stimulated_expression.csv", sce = sce),
    "accuracy1": expand(output + "data/{sce}_cell2loc_model/train_accuracy_1.png", sce = sce),
    "accuracy2": expand(output + "data/{sce}_cell2loc_model/train_accuracy_2.png", sce = sce),
    "history": expand(output + "data/{sce}_cell2loc_model/train_history.png", sce = sce)
}

cell2loc_predict = {
    "mat2": expand(output + "data/{sce}_cell2loc_model/predict/{sample}/celltype_abundances.csv", sample=sample_ids, sce = sce),
    "plot": expand(output + "data/{sce}_cell2loc_model/predict/{sample}/celltype_abundances.png", sample=sample_ids, sce = sce)
}

# Conditionally add output files based on user specifications in config
if config["plot_MT_HB"]:
    process_data.append(output + "figures/qc/MT-HB-counts-violin.png")


rule create_spatialdata:
    input:
        feature_matrix = "data/visium/{sample}/filtered_feature_bc_matrix.h5",
        scale_factors = "data/visium/{sample}/spatial/scalefactors_json.json",
        hires = "data/visium/{sample}/spatial/tissue_hires_image.png",
        lowres = "data/visium/{sample}/spatial/tissue_lowres_image.png",
        tissue_pos = "data/visium/{sample}/spatial/tissue_positions_list.csv",
    params:
        input_dir = directory("data/visium/{sample}/")
    threads: 1
    output:
        qc_fig = report(
            output + "figures/qc/Raw-{sample}-spatial-qc.png",
            category="Step 1: Data QC",
        ),
        zarr = directory(output + "data/raw/{sample}/"),
    script:
        "process-data/create-visium.py"

rule segment_nuclei:
    input:
        feature_matrix = "data/visium/{sample}/filtered_feature_bc_matrix.h5",
        scale_factor = "data/visium/{sample}/spatial/scalefactors_json.json",
        hires = "data/visium/{sample}/spatial/tissue_hires_image.png",
        lowres = "data/visium/{sample}/spatial/tissue_lowres_image.png",
        tissue_pos = "data/visium/{sample}/spatial/tissue_positions_list.csv"
    params:
        input_dir = directory("data/visium/{sample}/")
    threads: 1
    output:
        spots = directory(output + "data/spots/{sample}/"),
        counts = output + "data/spots/{sample}/" + "nuclei_count.csv"
    script:
        "process-data/segmentnuclei.py"

rule cell2loc_train:
    input:
        h5ad = "data/singlecell/scRNASeq-SingleR-annotated-sce-{sce}.h5ad",
    params:
        input_dir = directory("data/singlecell/"),
        epochs = config['epochs']
    output:
        h5ad = output + "data/{sce}_cell2loc_model/model_adata.h5ad",
        model = directory(output + "data/{sce}_cell2loc_model"),
        mat = output + "data/{sce}_cell2loc_model/stimulated_expression.csv",
        accuracy1 = output + "data/{sce}_cell2loc_model/train_accuracy_1.png",
        accuracy2 = output + "data/{sce}_cell2loc_model/train_accuracy_2.png",
        history = output + "data/{sce}_cell2loc_model/train_history.png"
    threads: 10
    script:
        "process-data/cell2loc_train.py"

rule cell2loc_predict:
    input:
        h5ad = "data/singlecell/scRNASeq-SingleR-annotated-sce-{sce}.h5ad",
        nuclei_counts = rules.segment_nuclei.output.counts,
        stim_expr = rules.cell2loc_train.output.mat,
        model = rules.cell2loc_train.output.model
    params:
        epochs = config['epochs']
    threads: 10
    output:
        mat = output + "data/{sce}_cell2loc_model/predict/{sample}/celltype_abundances.csv",
        plot = output + "data/{sce}_cell2loc_model/predict/{sample}/celltype_abundances.png",
        h5ad = output + "data/{sce}_cell2loc_model/predict/{sample}/{sample}_predict.h5ad",
        model = directory(output + "data/{sce}_cell2loc_model/predict/{sample}/")
    script:
        "process-data/cell2loc_predict.py"

rule combine_anndata_and_plot_counts_violin:
    input:
        spatial_data_out = expand(
            output + "figures/qc/Raw-{sample}-spatial-qc.png", sample=sample_ids
        ),
    params:
        spatial_datas = expand(output + "data/raw/{sample}/", sample=sample_ids),
    threads: 1
    output:
        counts_violin=report(
            output + "figures/qc/gene-counts-violin.png", category="Step 1: Data QC"
        ),
        zarr = directory(output + "data/raw/combined/"),
        finish_message = temp(output + "messages/finish-combining-visium.txt"),
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
        spatial_data_out=output + "messages/finish-combining-visium.txt",
    params:
        spatialdata_dir=output + "data/raw/combined/",
        pct_counts_hb_thr=config["filtering"]["max_pct_hb"],
        pct_counts_mt_thr=config["filtering"]["max_pct_mt"],
        n_genes_by_counts_thr=config["filtering"]["min_n_genes_by_counts"]
    threads: 1
    output:
        filtered_spatialdata=directory(output + "data/clean/"),
        finish_message=temp(output + "messages/finish-visium-filter.txt"),
    script:
        "process-data/visium-filter.py"


rule plot_counts_filtered:
    input:
        message=output + "messages/finish-visium-filter.txt",
    params:
        filtered_visium_spatialdata=output + "data/clean/"
    threads: 1
    output:
        qc_fig=report(
            output + "figures/qc/Filtered-{sample}-spatial-qc.png",
            category="Step 1: Data QC",
        ),
    script:
        "process-data/plot-visium-qc.py"
