import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- CONFIGURATION ---
CATH_OUTPUTS_DIR = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/cath_outputs_abs_scale"
OUTPUT_DIR = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/stability_comparison/scale"
BIN_EDGES = list(range(0, 101, 10))  # 0-10, 10-20, ..., 90-100
# Nur diese dG_X-Ordner beruecksichtigen (andere werden ignoriert)
ALLOWED_GROUPS = {"dG_0","dG_0.001", "dG_0.01", "dG_0.1"}
# ------------------------------------------------

T_PERCENT_RE = re.compile(r"T%=([\-0-9.]+)")


def load_t_percentages(base_dir):
    """
    Durchsucht base_dir/dG_X/<protein_dir>/designs.fasta nach den
    T%-Werten in jeder FASTA-Kopfzeile. Gibt Liste von Dicts zurück.
    """
    rows = []
    for dg_folder in sorted(os.listdir(base_dir)):
        dg_path = os.path.join(base_dir, dg_folder)
        if not os.path.isdir(dg_path) or dg_folder not in ALLOWED_GROUPS:
            continue

        for protein_dir in sorted(os.listdir(dg_path)):
            fasta_path = os.path.join(dg_path, protein_dir, "designs.fasta")
            if not os.path.isfile(fasta_path):
                continue

            cath_id = protein_dir.split("_")[0]
            with open(fasta_path) as f:
                for line in f:
                    if not line.startswith(">"):
                        continue
                    match = T_PERCENT_RE.search(line)
                    if not match:
                        continue
                    rows.append({
                        "group": dg_folder,
                        "cath_id": cath_id,
                        "T_percent": float(match.group(1)),
                    })
    return rows


def sort_key(group_name):
    """dG_X aufsteigend nach Zahl sortieren."""
    try:
        return float(group_name.replace("dG_", ""))
    except ValueError:
        return float("inf")


def main():
    if not os.path.isdir(CATH_OUTPUTS_DIR):
        print(f"Fehler: '{CATH_OUTPUTS_DIR}' wurde nicht gefunden!")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    rows = load_t_percentages(CATH_OUTPUTS_DIR)
    print(f"Gefundene Designs mit T%-Wert: {len(rows)}")

    if not rows:
        print("Keine T%-Werte gefunden, breche ab.")
        return

    df = pd.DataFrame(rows)

    # Rohdaten speichern
    raw_csv_path = os.path.join(OUTPUT_DIR, "all_T_percent.csv")
    df.to_csv(raw_csv_path, index=False)
    print(f"-> Rohdaten gespeichert unter '{raw_csv_path}'")

    groups = sorted(df["group"].unique(), key=sort_key)

    # Histogramm-Bins pro Gruppe berechnen
    bin_labels = [f"{BIN_EDGES[i]}-{BIN_EDGES[i+1]}" for i in range(len(BIN_EDGES) - 1)]
    n_bins = len(bin_labels)

    counts_per_group = {}
    for group in groups:
        values = df.loc[df["group"] == group, "T_percent"]
        counts, _ = np.histogram(values, bins=BIN_EDGES)
        counts_per_group[group] = counts

    # Zusammenfassung als CSV
    summary_df = pd.DataFrame(counts_per_group, index=bin_labels)
    summary_csv_path = os.path.join(OUTPUT_DIR, "T_percent_histogram_counts.csv")
    summary_df.to_csv(summary_csv_path)
    print(f"-> Histogramm-Zaehlwerte gespeichert unter '{summary_csv_path}'")
    print(summary_df)

    # Gruppiertes Balkendiagramm, gleiche Optik wie "bars with legend"
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
    ax.set_title("Distribution of T% per Weight")
    ax.legend(title="Group")
    plt.tight_layout()

    plot_path = os.path.join(OUTPUT_DIR, "T_percent_histogram.png")
    plt.savefig(plot_path, dpi=150)
    print(f"-> Histogramm gespeichert unter '{plot_path}'")


if __name__ == "__main__":
    main()
