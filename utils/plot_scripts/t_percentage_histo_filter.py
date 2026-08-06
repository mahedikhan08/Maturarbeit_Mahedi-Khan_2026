import os
import re
import glob
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- CONFIGURATION ---
CATH_OUTPUTS_DIR = "path to fasta files"
SIEVE_OUTPUT_DIR = "path to structures"
OUTPUT_DIR = "output path"
BIN_EDGES = list(range(0, 101, 10))  # 0-10, 10-20, ..., 90-100
ALLOWED_GROUPS = {"dG_0","dG_0.001","dG_0.01","dG_0.1"}
# ------------------------------------------------

T_PERCENT_RE = re.compile(r"T%=([\-0-9.]+)")


def load_tm_scores(sieve_output_dir):
    """
    Liest alle structure_scores.json Dateien aus sieve_output_dir/dG_X/structure_scores.json.
    Gibt ein verschachteltes Dictionary zurück: { group: { full_design_id: tm_model_wert } }
    """
    tm_scores = {g: {} for g in ALLOWED_GROUPS}
    pattern = os.path.join(sieve_output_dir, "dG_*", "structure_scores.json")
    
    for path in glob.glob(pattern):
        # Extrahiere den Gruppen-Ordnernamen (z.B. "dG_0.1") aus dem Pfad
        group_folder = os.path.basename(os.path.dirname(path))
        if group_folder not in ALLOWED_GROUPS:
            continue
            
        try:
            with open(path) as f:
                data = json.load(f)
            for full_id, scores in data.items():
                if "tm_model" in scores:
                    tm_scores[group_folder][full_id] = scores["tm_model"]
        except (json.JSONDecodeError, KeyError, IOError) as e:
            print(f"Warnung: Konnte Scores aus '{path}' nicht laden ({e}).")
            
    return tm_scores


def load_t_percentages(base_dir, tm_scores):
    """
    Durchsucht base_dir/dG_X/<protein_dir>/designs.fasta nach den
    T%-Werten. Filtert nach tm_model > 0.7 (gruppenspezifisch).
    """
    rows = []
    skipped_filter = 0
    skipped_no_score = 0

    for dg_folder in sorted(os.listdir(base_dir)):
        dg_path = os.path.join(base_dir, dg_folder)
        if not os.path.isdir(dg_path) or dg_folder not in ALLOWED_GROUPS:
            continue

        group_scores = tm_scores.get(dg_folder, {})

        for protein_dir in sorted(os.listdir(dg_path)):
            fasta_path = os.path.join(dg_path, protein_dir, "designs.fasta")
            if not os.path.isfile(fasta_path):
                continue

            cath_id = protein_dir.split("_")[0]
            with open(fasta_path) as f:
                for line in f:
                    if not line.startswith(">"):
                        continue
                    
                    # Extrahiere die ID (z.B. "designed_1")
                    header_id = line.split("|")[0].strip(">")
                    full_design_id = f"{protein_dir}_{header_id}"

                    # Filter-Logik mit gruppenspezifischem Score
                    if full_design_id not in group_scores:
                        skipped_no_score += 1
                        continue
                    
                    if group_scores[full_design_id] <= 0.7:
                        skipped_filter += 1
                        continue

                    match = T_PERCENT_RE.search(line)
                    if not match:
                        continue
                    rows.append({
                        "group": dg_folder,
                        "cath_id": cath_id,
                        "T_percent": float(match.group(1)),
                    })

    print(f"-> Filter-Statistik: {skipped_filter} Sequenzen mit tm_model <= 0.7 entfernt.")
    if skipped_no_score > 0:
        print(f"-> Warnung: Für {skipped_no_score} Sequenzen wurde kein TM-Score in der jeweiligen Gruppe gefunden.")

    return rows


def sort_key(group_name):
    try:
        return float(group_name.replace("dG_", ""))
    except ValueError:
        return float("inf")


def main():
    if not os.path.isdir(CATH_OUTPUTS_DIR) or not os.path.isdir(SIEVE_OUTPUT_DIR):
        print("Fehler: Eingangsverzeichnisse prüfen!")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Lade TM-Scores gruppiert aus structure_scores.json...")
    tm_scores = load_tm_scores(SIEVE_OUTPUT_DIR)
    
    for g in ALLOWED_GROUPS:
        print(f"  -> {g}: {len(tm_scores[g])} Scores geladen.")

    print("\nLade und filtere FASTA-Dateien...")
    rows = load_t_percentages(CATH_OUTPUTS_DIR, tm_scores)
    print(f"Gefundene Designs nach Filter: {len(rows)}")

    if not rows:
        print("Keine T%-Werte nach Filterung übrig, breche ab.")
        return

    df = pd.DataFrame(rows)

    raw_csv_path = os.path.join(OUTPUT_DIR, "all_T_percent_filtered.csv")
    df.to_csv(raw_csv_path, index=False)

    groups = sorted(df["group"].unique(), key=sort_key)
    bin_labels = [f"{BIN_EDGES[i]}-{BIN_EDGES[i+1]}" for i in range(len(BIN_EDGES) - 1)]
    n_bins = len(bin_labels)

    counts_per_group = {}
    for group in groups:
        values = df.loc[df["group"] == group, "T_percent"]
        counts, _ = np.histogram(values, bins=BIN_EDGES)
        counts_per_group[group] = counts

    summary_df = pd.DataFrame(counts_per_group, index=bin_labels)
    summary_csv_path = os.path.join(OUTPUT_DIR, "T_percent_histogram_counts_filtered.csv")
    summary_df.to_csv(summary_csv_path)
    print("\n", summary_df)

    palette = ["red", "tan", "lime", "#7F77DD", "#D4537E"]
    color_map = {g: palette[i % len(palette)] for i, g in enumerate(groups)}

    x = np.arange(n_bins)
    n_groups = len(groups)
    bar_width = 0.8 / n_groups

    fig, ax = plt.subplots(figsize=(10, 6))
    for i, group in enumerate(groups):
        offset = (i - (n_groups - 1) / 2) * bar_width
        ax.bar(
            x + offset,
            counts_per_group[group],
            width=bar_width,
            label=group,
            color=color_map[group],
            edgecolor="black",
            linewidth=0.5,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=45, ha="right")
    ax.set_xlabel("T%")
    ax.set_ylabel("Sequences")
    ax.set_title("Distribution of T% per Weight (Filtered: tm_model > 0.7)")
    ax.legend(title="Group")
    plt.tight_layout()

    plot_path = os.path.join(OUTPUT_DIR, "T_percent_histogram_filtered.png")
    plt.savefig(plot_path, dpi=150)
    print(f"-> Histogramm gespeichert unter '{plot_path}'")


if __name__ == "__main__":
    main()