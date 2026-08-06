import os

from make_stability_jobs import find_structures as find_design_structures
from make_stability_jobs import (
    SIEVE_OUTPUT_DIR,
    TARGET_GROUP,  # <--- Importiert die neue Variable aus dem vorherigen Skript
    WORKDIR as DESIGN_WORKDIR,
    CHAIN_ID as DESIGN_CHAIN_ID,
    WEIGHTS,
    RUN_SCRIPT,
)
from make_wild_abs_job import find_wildtype_structures
from make_wild_abs_job import WILDTYPE_DIR as RAW_WILDTYPE_DIR, WORKDIR as WT_WORKDIR

#CONFIGURATION
CMD_FILE = "commands_missing.cmd"
SLURM_FILE = "submit_array_missing.sh"


def result_exists(struct_id, workdir_root):
    return os.path.isfile(os.path.join(workdir_root, struct_id, "dG_result.json"))


def main():
    weights_arg = " ".join(WEIGHTS)
    missing_cmds = []

    # 1. Fehlende Design-Strukturen finden
    print(f"Searching structures for group '{TARGET_GROUP}'...")
    design_structs = find_design_structures(SIEVE_OUTPUT_DIR, TARGET_GROUP)
    missing_design = [
        (struct_id, path) for struct_id, path in design_structs
        if not result_exists(struct_id, DESIGN_WORKDIR)
    ]
    for struct_id, struct_path in missing_design:
        out_file = f"{DESIGN_WORKDIR}/{struct_id}/dG_result.json"
        cmd = (
            f"python {RUN_SCRIPT} --input {struct_path} --chain {DESIGN_CHAIN_ID} "
            f"--output {out_file} --weights {weights_arg}\n"
        )
        missing_cmds.append(cmd)

    # 2. Fehlende Wildtyp-Strukturen finden
    wt_structs = find_wildtype_structures(RAW_WILDTYPE_DIR)
    missing_wt = [
        (struct_id, path, chain) for struct_id, path, chain in wt_structs
        if not result_exists(struct_id, WT_WORKDIR)
    ]
    for struct_id, struct_path, chain in missing_wt:
        out_file = f"{WT_WORKDIR}/{struct_id}/dG_result.json"
        cmd = (
            f"python {RUN_SCRIPT} --input {struct_path} --chain {chain} "
            f"--output {out_file} --weights {weights_arg}\n"
        )
        missing_cmds.append(cmd)

    num_missing = len(missing_cmds)
    print(f"Missing Design Structures in {TARGET_GROUP}: {len(missing_design)}")
    print(f"Missing Wildtype Structures: {len(missing_wt)}")
    print(f"Total Jobs to Repeat: {num_missing}")

    if num_missing == 0:
        print("Nothing is missing")
        return

    # 3. commands_missing.cmd schreiben
    with open(CMD_FILE, "w") as f_cmd:
        f_cmd.writelines(missing_cmds)
    print(f"-> '{CMD_FILE}' was succesfully made")

    # 4. submit_array_missing.sh schreiben
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name=stability_missing_array
#SBATCH --partition=schwede-h200,schwede-a100,a100,a100-80g,l40s
#SBATCH --qos=schwede,a100-1day,l40s-1day
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-{num_missing}
#SBATCH --output=logs/stability_missing_output_%A_%a.log

# Umgebung laden
ml Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate stability

# Offline-Modus fuer Hugging Face (Gewichte liegen schon lokal vor)
export HF_HOME=<path to hf_cache>
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

# Log-Ordner erstellen falls nicht da
mkdir -p logs

# Befehl auswaehlen und ausfuehren
SEEDFILE="{CMD_FILE}"
SEED=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p $SEEDFILE)
echo "Starte Aufgabe $SLURM_ARRAY_TASK_ID mit Befehl: $SEED"
eval $SEED
echo "Job $SLURM_ARRAY_TASK_ID abgeschlossen am: $(date)"
"""
    with open(SLURM_FILE, "w") as f_slurm:
        f_slurm.write(slurm_content)
    print(f"-> '{SLURM_FILE}' was successfully made")

    os.chmod(SLURM_FILE, 0o755)


if __name__ == "__main__":
    main()