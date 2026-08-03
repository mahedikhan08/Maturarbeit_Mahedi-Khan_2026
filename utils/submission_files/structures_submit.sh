#!/bin/bash
#SBATCH --job-name=structurescorer
#SBATCH --partition=schwede-h200,schwede-a100,a100,a100-80g,l40s
#SBATCH --qos=schwede,a100-6hours,l40s-6hours
#SBATCH --gres=gpu:1
#SBATCH --output=logs/structure_output.log

# Umgebung laden
ml Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate /scicore/home/schwede/pudziu0000/mambaforge/envs/colabfold

# Offline-Modus fuer Hugging Face
export TRANSFORMERS_OFFLINE=1

# Log-Ordner erstellen falls nicht da
mkdir -p logs

# Befehl auswaehlen und ausfuehren
cd /scicore/home/schwede/khan0010/project/tea-leaves/sieve

snakemake --use-conda --profile slurm -s Snakefile structure_scores --cores 64 --configfile /scicore/home/schwede/khan0010/project/tea-leaves-workdir/cath_seq/config.yaml
