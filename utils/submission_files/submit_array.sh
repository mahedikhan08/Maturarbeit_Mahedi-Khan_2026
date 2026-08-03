#!/bin/bash
#SBATCH --job-name=tea_leaves_array
#SBATCH --partition=schwede-h200,schwede-a100,a100,a100-80g,l40s
#SBATCH --qos=schwede,a100-1day,l40s-1day
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=82
#SBATCH --output=logs/tea_output_%A_%a.log

# Umgebung laden
ml Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate tea-leaves

# Offline-Modus fuer Hugging Face
export TRANSFORMERS_OFFLINE=1

# Log-Ordner erstellen falls nicht da
mkdir -p logs

# Befehl auswaehlen und ausfuehren
SEEDFILE="commands.cmd"
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)

echo "Starte Aufgabe $SLURM_ARRAY_TASK_ID mit Befehl: $SEED"
eval $SEED

echo "Job $SLURM_ARRAY_TASK_ID abgeschlossen am: $(date)"
