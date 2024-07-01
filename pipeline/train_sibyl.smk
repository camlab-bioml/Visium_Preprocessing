cell2loc_train_sibyl = {
    "model": expand(output + "data/sibyl_cell2loc_model/model_adata.h5ad"),
    "pt": expand(output + "data/sibyl_cell2loc_model"),
    "mat": expand(output + "data/sibyl_cell2loc_model/stimulated_expression.csv"),
    "accuracy1": expand(output + "data/sibyl_cell2loc_model/train_accuracy_1.png"),
    "accuracy2": expand(output + "data/sibyl_cell2loc_model/train_accuracy_2.png"),
    "history": expand(output + "data/sibyl_cell2loc_model/train_history.png")
}

cell2loc_predict_sibyl = {
    "mat2": expand(output + "data/sibyl_cell2loc_model/predict/{sample}/celltype_abundances.csv", sample=sample_ids),
    "plot": expand(output + "data/sibyl_cell2loc_model/predict/{sample}/celltype_abundances.png", sample=sample_ids)
}

rule cell2loc_train_sibyl:
    input:
        h5ad = "data/singlecell/sibyl_merged.h5ad",
    params:
        input_dir = directory("data/singlecell/"),
        epochs = config['epochs']
    singularity: '/ddn_exa/campbell/share/containers/scvi-v2.sif'
    output:
        h5ad = output + "data/sibyl_cell2loc_model/model_adata.h5ad",
        model = directory(output + "data/sibyl_cell2loc_model"),
        mat = output + "data/sibyl_cell2loc_model/stimulated_expression.csv",
        accuracy1 = output + "data/sibyl_cell2loc_model/train_accuracy_1.png",
        accuracy2 = output + "data/sibyl_cell2loc_model/train_accuracy_2.png",
        history = output + "data/sibyl_cell2loc_model/train_history.png"
    threads: 10
    script:
        "process-data/cell2loc_train_sibyl.py"

rule cell2loc_predict_sibyl:
    input:
        h5ad = "data/singlecell/sibyl_merged.h5ad",
        nuclei_counts = rules.segment_nuclei.output.counts,
        stim_expr = rules.cell2loc_train_sibyl.output.mat,
        model = rules.cell2loc_train_sibyl.output.model
    params:
        epochs = 30000
    singularity: '/ddn_exa/campbell/share/containers/scvi-v2.sif'
    threads: 10
    output:
        mat = output + "data/sibyl_cell2loc_model/predict/{sample}/celltype_abundances.csv",
        plot = output + "data/sibyl_cell2loc_model/predict/{sample}/celltype_abundances.png",
        h5ad = output + "data/sibyl_cell2loc_model/predict/{sample}/{sample}_predict.h5ad",
        model = directory(output + "data/sibyl_cell2loc_model/predict/{sample}/")
    script:
        "process-data/cell2loc_predict_sibyl.py"