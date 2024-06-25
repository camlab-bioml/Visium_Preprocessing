#!/bin/bash
#SBATCH --job-name=visium_preprocess
#SBATCH --cpus-per-task=10
#SBATCH --output=/ddn_exa/campbell/ttan2/visium_preprocessing/slurm/visium_preprocess.out
#SBATCH --error=/ddn_exa/campbell/ttan2/visium_preprocessing/slurm/visium_preprocess.err
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --qos=normal

cd /ddn_exa/campbell/ttan2/visium_preprocessing
source /home/campbell/ttan2/.local/share/virtualenvs/visium_preprocessing-UpaDixbq/bin/activate
which python
snakemake --unlock
snakemake -c 10 -R cell2loc_predict --rerun-incomplete --use-singularity --singularity-args '--nv -B /ddn_exa/campbell/ttan2/visium_preprocessing'
exit
