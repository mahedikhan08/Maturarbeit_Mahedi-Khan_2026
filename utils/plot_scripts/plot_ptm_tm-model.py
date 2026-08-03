import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# --- CONFIGURATION ---
BASE_DIR = "/scicore/home/schwede/khan0010/project/tea-leaves-workdir/cath_seq/sieve_output/scale"
DG_CSV_PATH = "/scicore/home/schwede/khan0010/project/tea-leaves-workdir/cath_seq/stability_comparison/scale/all_dG_results.csv"
OUTPUT_DIR = "/scicore/home/schwede/khan0010/project/tea-leaves-workdir/cath_seq/stability_comparison/scale"
GROUPS = ["dG_0"]
# ------------------------------------------------


def shorten_id(full_id):
    match = re.match(r"^([a-zA-Z0-9]+)_.*?_designed_(\d+)", full_id)
    if match:
        return f"{match.group(1)}__d{match.group(2)}"
    return full_id[:10]  # Fallback, falls das Muster bricht


def load_dg_lookup(csv_path):
    """
    Laedt all_dG_results.csv und baut ein Lookup (group, design_id) -> dG.
    Nur Design-Zeilen (group != 'Wildtyp') werden beruecksichtigt, da nur
    diese eine design_id passend zu den structure_scores.json Schluesseln haben.
    """
    if not os.path.isfile(csv_path):
        print(f"Fehler: '{csv_path}' nicht gefunden. Fuehre zuerst compare_stability.py aus.")
        return {}
    df = pd.read_csv(csv_path)
    df = df[df["group"] != "Wildtyp"]
    return {(row["group"], row["design_id"]): row["dG"] for _, row in df.iterrows()}


def generate_scatter_plot(group, dg_lookup):
    json_file = os.path.join(BASE_DIR, group, "structure_scores.json")
    if not os.path.exists(json_file):
        print(f"Error: {json_file} nicht gefunden.")
        return

    with open(json_file, "r") as f:
        data = json.load(f)

    parsed_data = []
    for protein_id, metrics in data.items():
        ptm = metrics.get("ptm")
        tm_model = metrics.get("tm_model")
        dg = dg_lookup.get((group, protein_id))

        if ptm is not None and tm_model is not None and ptm != "N/A" and tm_model != "N/A" and dg is not None:
            parsed_data.append({
                "Protein_ID": protein_id,
                "Label": shorten_id(protein_id),
                "pTM": float(ptm),
                "tm_model": float(tm_model),
                "dG": float(dg),
            })

    df = pd.DataFrame(parsed_data)
    print(f"{group}: {len(df)} proteins found")
    if df.empty:
        print(f"Keine gueltigen Datenpunkte fuer '{group}' zum Plotten gefunden.")
        return

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 9))

    # Scatter Plot: x=pTM, y=tm_model, Farbe=dG (aus dem unabhaengigen Predictor)
    scatter = plt.scatter(
        data=df,
        x="pTM",
        y="tm_model",
        c="dG",
        cmap="viridis",
        vmin=df["dG"].min(),
        vmax=df["dG"].max(),
        alpha=0.6,
        edgecolors="w",
        s=60,
    )

    # Titles
    plt.title(f"{group}: pTM vs. TM-score, colored after independent \u0394G", fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("pTM", fontsize=12)
    plt.ylabel("TM-score(tm_model)", fontsize=12)

    plt.axvline(x=0.7, color="gray", linestyle="--", alpha=0.5, label="pTM = 0.7")
    plt.axhline(y=0.7, color="gray", linestyle=":", alpha=0.5, label="tm_model = 0.7")
    plt.xlim(0, 1.05)
    plt.ylim(0, 1.05)

    # Color
    cbar = plt.colorbar(scatter)
    cbar.set_label("\u0394G (kcal/mol)", fontsize=10)

    plt.tight_layout()
    output_image = os.path.join(OUTPUT_DIR, f"ptm_vs_tm_model_dG_{group}.png")
    plt.savefig(output_image, dpi=300)
    plt.close()
    print(f"-> '{output_image}' got saved")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    dg_lookup = load_dg_lookup(DG_CSV_PATH)
    for group in GROUPS:
        generate_scatter_plot(group, dg_lookup)


if __name__ == "__main__":
    main()