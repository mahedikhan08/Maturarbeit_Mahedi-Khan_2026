#!/bin/bash
#SBATCH --job-name=structurescorer
#SBATCH --partition=schwede-h200,schwede-a100,a100,a100-80g,l40s
#SBATCH --qos=schwede,a100-6hours,l40s-6hours
#SBATCH --gres=gpu:1
#SBATCH --output=logs/structure_output.log

# Umgebung laden
ml Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate <path to colabfold_env>

# Offline-Modus fuer Hugging Face
export TRANSFORMERS_OFFLINE=1

# Log-Ordner erstellen falls nicht da
mkdir -p logs

# Befehl auswaehlen und ausfuehren
cd <path to sieve in tea-leaves>

snakemake --use-conda --profile slurm -s Snakefile structure_scores --cores 64 --configfile <path to config.yaml>