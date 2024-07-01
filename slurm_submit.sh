#!/bin/bash
#SBATCH --job-name=cell2loc_train_sdrisser
#SBATCH --cpus-per-task=24
#SBATCH --output=/ddn_exa/campbell/ttan2/visium_preprocessing/slurm/visium_preprocess_sdrissler.out
#SBATCH --error=/ddn_exa/campbell/ttan2/visium_preprocessing/slurm/visium_preprocess_sdrissler.err
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --partition=ddn
#SBATCH --qos=normal

cd /ddn_exa/campbell/ttan2/visium_preprocessing
source /home/campbell/ttan2/.local/share/virtualenvs/visium_preprocessing-UpaDixbq/bin/activate
which python
nvidia-smi
snakemake --unlock
snakemake -c 24 --touch
snakemake -c 24 --rerun-incomplete --use-singularity --singularity-args '--nv -B /ddn_exa/campbell/ttan2/visium_preprocessing'
exit
