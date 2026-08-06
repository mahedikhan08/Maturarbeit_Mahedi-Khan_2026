import csv
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- CONFIGURATION ---
BASE_DIR = "path to structures"
FASTA_BASE_DIR = "path to fasta files"
OUTPUT_DIR = "path for plot"
FASTA_OUTPUT_DIR = "path for fasta"
DG_CSV_PATH = "path to all_dG_.sv"
GROUPS = ["dG_0","dG_0.001","dG_0.01","dG_0.1"]
TM_THRESHOLD = 0.7 
JITTER_WIDTH = 0.15
RNG_SEED = 42
# ------------------------------------------------

L_RE = re.compile(r"\|L=(\d+)\|")


def get_protein_length(fasta_base_dir, group, cath_id):
    #Liest L= aus der ersten Kopfzeile von designs.fasta fuer dieses Protein aus FASTA_BASE_DIR
    fasta_path = os.path.join(fasta_base_dir, group, cath_id, "designs.fasta")
    if not os.path.isfile(fasta_path):
        return None
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                match = L_RE.search(line)
                if match:
                    return int(match.group(1))
    return None


def load_dg_scores(csv_path):

    #Liest die dG-Werte aus der CSV-Datei und ignoriert Wildtyp

    dg_dict = {}
    if not os.path.isfile(csv_path):
        print(f"Warnung: dG CSV-Datei '{csv_path}' nicht gefunden.")
        return dg_dict

    with open(csv_path, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["group"] == "Wildtyp":
                continue
            try:
                group = row["group"]
                design_id = row["design_id"]
                dg_val = float(row["dG"])
                dg_dict[(group, design_id)] = dg_val
            except (ValueError, KeyError):
                continue
    return dg_dict


def load_group_scores(base_dir, fasta_base_dir, group, dg_dict):
    #Liest structure_scores.json fuer eine dG_X-Gruppe ein und matcht die dG-Werte.
    scores_path = os.path.join(base_dir, group, "structure_scores.json")
    if not os.path.isfile(scores_path):
        print(f"Warnung: '{scores_path}' nicht gefunden, ueberspringe Gruppe '{group}'.")
        return []

    with open(scores_path) as f:
        data = json.load(f)

    length_cache = {}
    rows = []
    for key, values in data.items():
        idx = key.rfind("_designed_")
        if idx == -1:
            print(f"Warnung: Schluessel '{key}' passt nicht zum erwarteten Muster, ueberspringe.")
            continue
        cath_id = key[:idx]
        design_tag = key[idx + 1:]  # z.B. "designed_1"

        if cath_id not in length_cache:
            length_cache[cath_id] = get_protein_length(fasta_base_dir, group, cath_id)
        length = length_cache[cath_id]

        dg_value = dg_dict.get((group, key), None)

        rows.append({
            "design_key": key,        # Vollständiger Key (z.B. 1a0iA02_..._designed_1)
            "cath_id": cath_id,
            "design_tag": design_tag,  # z.B. designed_1
            "short_id": cath_id.split("_")[0],
            "length": length,
            "tm_model": values.get("tm_model"),
            "ptm": values.get("ptm"),
            "dg": dg_value,
        })
    return rows


def load_fasta_sequences(fasta_base_dir, group, cath_ids):
    #Liest die Aminosäure-Sequenzen aus allen relevanten designs.fasta-Dateien einer Gruppe.
    sequences = {}
    for cath_id in cath_ids:
        fasta_path = os.path.join(fasta_base_dir, group, cath_id, "designs.fasta")
        if not os.path.isfile(fasta_path):
            print(f"Warnung: FASTA '{fasta_path}' nicht gefunden.")
            continue

        current_header = None
        current_seq = []
        
        with open(fasta_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_header and current_seq:
                        tag = current_header[1:].split("|")[0]
                        sequences[(cath_id, tag)] = "".join(current_seq)
                    current_header = line
                    current_seq = []
                elif line:
                    current_seq.append(line)
            
            if current_header and current_seq:
                tag = current_header[1:].split("|")[0]
                sequences[(cath_id, tag)] = "".join(current_seq)

    return sequences


def write_best_sequences_fasta(rows, group, fasta_base_dir, fasta_output_dir, tm_threshold):
    #Filtert Designs mit TM-Score >= tm_threshold, sortiert sie zuerst alphabetisch nach Protein-Name (cath_id) und darunter aufsteigend nach dG-Wert. Schreibt sie anschliessend in eine FASTA-Datei im Zielordner.
    filtered_rows = [
        r for r in rows 
        if r["tm_model"] is not None 
        and r["tm_model"] >= tm_threshold 
        and r["dg"] is not None
    ]

    if not filtered_rows:
        print(f"Keine Sequenzen in {group} mit TM-Score >= {tm_threshold} gefunden.")
        return

    # Sortierung: 1. cath_id alphabetisch, 2. dG aufsteigend
    filtered_rows.sort(key=lambda x: (x["cath_id"], x["dg"]))

    cath_ids = {r["cath_id"] for r in filtered_rows}
    seq_dict = load_fasta_sequences(fasta_base_dir, group, cath_ids)

    out_fasta_filename = f"{group}_best_seq.fasta"
    out_fasta_path = os.path.join(fasta_output_dir, out_fasta_filename)

    written_count = 0
    with open(out_fasta_path, "w", encoding="utf-8") as out_f:
        for r in filtered_rows:
            key_pair = (r["cath_id"], r["design_tag"])
            sequence = seq_dict.get(key_pair)

            if not sequence:
                print(f"Warnung: Keine Sequenz für {r['design_key']} gefunden.")
                continue

            # Länge direkt aus der geladenen Sequenz berechnen
            seq_len = len(sequence)

            header = (
                f">{r['design_key']} "
                f"dG={r['dg']:.3f} | pTM={r['ptm']:.4f} | TM-score={r['tm_model']:.4f} | L={seq_len}"
            )
            out_f.write(f"{header}\n{sequence}\n")
            written_count += 1

    print(f"-> FASTA mit den besten Sequenzen gespeichert unter '{out_fasta_path}' ({written_count} Sequenzen).")


def plot_group(rows, group, output_dir):
    if not rows:
        return

    proteins = {}
    total_success_designs = 0
    
    for r in rows:
        proteins.setdefault(
            r["cath_id"], 
            {"short_id": r["short_id"], "length": r["length"], "success_count": 0}
        )
        if r["tm_model"] is not None and r["tm_model"] >= TM_THRESHOLD:
            proteins[r["cath_id"]]["success_count"] += 1
            total_success_designs += 1

    ordered_ids = sorted(
        proteins.keys(),
        key=lambda cid: (proteins[cid]["length"] is None, proteins[cid]["length"] or 0)
    )
    y_index = {cid: i for i, cid in enumerate(ordered_ids)}

    rng = np.random.default_rng(RNG_SEED)
    
    xs, ys = [], []
    colors_ptm, colors_dg = [], []
    
    for r in rows:
        if r["tm_model"] is None:
            continue
        base_y = y_index[r["cath_id"]]
        xs.append(r["tm_model"])
        ys.append(base_y + rng.uniform(-JITTER_WIDTH, JITTER_WIDTH))
        colors_ptm.append(r["ptm"] if r["ptm"] is not None else 0)
        colors_dg.append(r["dg"] if r["dg"] is not None else 0)

    fig_height = max(5.5, len(ordered_ids) * 0.28 + 0.5)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, fig_height), sharey=True)

    sc1 = ax1.scatter(
        xs, ys, c=colors_ptm, cmap="Blues", vmin=0, vmax=1,
        s=45, edgecolor="black", linewidth=0.4, alpha=0.85,
    )
    ax1.axvline(TM_THRESHOLD, color="gray", linestyle="--", linewidth=1)
    ax1.set_xlim(0, 1.02)
    ax1.set_xlabel("TM-score (tm_model)")
    ax1.set_title("pTM")

    sc2 = ax2.scatter(
        xs, ys, c=colors_dg, cmap="viridis",
        s=45, edgecolor="black", linewidth=0.4, alpha=0.85,
    )
    ax2.axvline(TM_THRESHOLD, color="gray", linestyle="--", linewidth=1)
    ax2.set_xlim(0, 1.02)
    ax2.set_xlabel("TM-score (tm_model)")
    ax2.set_title("dG-Wert")

    labels = []
    for cid in ordered_ids:
        length = proteins[cid]["length"]
        length_str = str(length) if length is not None else "?"
        success_count = proteins[cid]["success_count"]
        labels.append(f"[{success_count}] {proteins[cid]['short_id']} ({length_str})")

    ax1.set_yticks(range(len(ordered_ids)))
    ax1.set_yticklabels(labels, fontsize=7)
    ax1.invert_yaxis()

    ax2.yaxis.set_tick_params(labelleft=True)
    ax2.set_yticks(range(len(ordered_ids)))
    ax2.set_yticklabels(labels, fontsize=7)

    cbar1 = fig.colorbar(sc1, ax=ax1, pad=0.04, shrink=0.8)
    cbar1.set_label("pTM")
    
    cbar2 = fig.colorbar(sc2, ax=ax2, pad=0.04, shrink=0.8)
    cbar2.set_label("dG")

    plt.suptitle(f"{group}", fontsize=12, weight="bold", y=0.98)
    
    fig.text(
        0.5, 0.015, 
        f"Overall successful sequences (TM-score >= {TM_THRESHOLD}): {total_success_designs}", 
        ha="center", fontsize=10, weight="bold", 
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round,pad=0.5')
    )

    plt.tight_layout(rect=[0, 0.04, 1, 0.96])
    
    out_path = os.path.join(output_dir, f"tm_model_double_swarmplot_{group}.png")
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"-> Plot gespeichert unter '{out_path}'")
    print(f"-> Gesamtzahl erfolgreicher Designs in {group}: {total_success_designs}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(FASTA_OUTPUT_DIR, exist_ok=True)
    
    print("Lade dG-Werte aus CSV...")
    dg_dict = load_dg_scores(DG_CSV_PATH)
    
    for group in GROUPS:
        rows = load_group_scores(BASE_DIR, FASTA_BASE_DIR, group, dg_dict)
        print(f"\n=== Gruppe: {group} ({len(rows)} Eintraege) ===")
        
        plot_group(rows, group, OUTPUT_DIR)
        write_best_sequences_fasta(
            rows, group, FASTA_BASE_DIR, FASTA_OUTPUT_DIR, TM_THRESHOLD
        )


if __name__ == "__main__":
    main()