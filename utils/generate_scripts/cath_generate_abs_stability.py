import os

# --- CONFIGURATION  ---
FASTA_FILE = "cath_sequences_selected.fasta"  # Name deiner FASTA-Datei
WORKDIR = "path to workdir"
STABILITY_WEIGHT = 0
LEN_MULTIPLIER = 200
DECAY_SCHEDULE = 2
# ------------------------------------------------

def parse_fasta(fasta_path):
    #Liest die FASTA-Datei aus und gibt eine Liste von (ID, Sequenz) zurück.
    sequences = []
    current_id = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    sequences.append((current_id, "".join(current_seq)))
                # Nimm das erste Wort nach '>' als ID (ohne Sonderzeichen/Leerzeichen)
                current_id = line[1:].split()[0].replace("/", "_").replace("|", "_")
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences.append((current_id, "".join(current_seq)))

    return sequences

def main():
    if not os.path.exists(FASTA_FILE):
        print(f"Fehler: Die Datei '{FASTA_FILE}' wurde nicht gefunden!")
        return

    # 1. FASTA parsen
    proteins = parse_fasta(FASTA_FILE)
    num_proteins = len(proteins)
    print(f"Erfolgreich {num_proteins} Proteine aus '{FASTA_FILE}' geladen.")

    # 2. commands.cmd schreiben
    with open("commands.cmd", "w") as f_cmd:
        for prot_id, seq in proteins:
            n = 5000
            step = 5000
            length = len(seq)
            while(length*LEN_MULTIPLIER >= n):
                n += step
            decay = int(n/DECAY_SCHEDULE)
            directory = f"{WORKDIR}/cath_outputs_abs_scale/dG_{STABILITY_WEIGHT}/{prot_id.replace('cath_4_4_0_', '')}_cath_4_4_0/"
            if not os.path.exists(directory):
                os.makedirs(directory)
            cmd = (
                f"tea_leaves --template {seq} -l {length} --tea_positions 0-{length-1} "
                f"--output {directory} "
                f"--temperature 0.005 --t_min 0 --cooling halving --decay {decay} "
                f"-n {n} -N 100 --abs_stability {STABILITY_WEIGHT} --log_mode simple "
                f"--aa_lm 2 --ce 1 --proposal lm --proposal_temp 1.0 --save_trajectory -s 42 \n"
            )
            f_cmd.write(cmd)
    print("-> 'commands.cmd' was succesfully made")

    # 3. submit_array.sh schreiben
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name=tea_leaves_array
#SBATCH --partition=schwede-h200,schwede-a100
#SBATCH --qos=schwede
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-1 #{num_proteins}
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
SEED=$(sed -n ${{SLURM_ARRAY_TASK_ID}}p $SEEDFILE)

echo "Starte Aufgabe $SLURM_ARRAY_TASK_ID mit Befehl: $SEED"
eval $SEED

echo "Job $SLURM_ARRAY_TASK_ID abgeschlossen am: $(date)"
"""

    with open("submit_array.sh", "w") as f_slurm:
        f_slurm.write(slurm_content)
    print("-> 'submit_array.sh' was successfully made")

    # Berechtigungen für das Shell-Skript setzen
    os.chmod("submit_array.sh", 0o755)

if __name__ == "__main__":
    main()
