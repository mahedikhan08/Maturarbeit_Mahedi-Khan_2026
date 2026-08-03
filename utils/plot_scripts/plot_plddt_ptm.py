import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- CONFIGURATION ---
JSON_FILE = "structure_scores.json"
OUTPUT_IMAGE = "plddt_vs_ptm_scatter.png"


def generate_scatter_plot():
    if not os.path.exists(JSON_FILE):
        print(f"Error: {JSON_FILE} nicht gefunden.")
        return

    with open(JSON_FILE, "r") as f:
        data = json.load(f)

    parsed_data = []
    for protein_id, metrics in data.items():
        plddt = metrics.get("plddt")
        ptm = metrics.get("ptm")
        
        if plddt is not None and ptm is not None and plddt != "N/A" and ptm != "N/A":
            parsed_data.append({
                "Protein_ID": protein_id,
                "pLDDT": float(plddt),
                "pTM": float(ptm)
            })

    df = pd.DataFrame(parsed_data)
    print(f" {len(df)} proteins found")

    if df.empty:
        print("Keine gültigen Proteindaten zum Plotten gefunden.")
        return

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 9)) 

    # Scatter Plot
    scatter = plt.scatter(
        data=df, 
        x="pLDDT", 
        y="pTM", 
        c="pTM", 
        cmap="viridis", 
        vmin=0.0,
        vmax=1.0,
        alpha=0.6, 
        edgecolors="w", 
        s=60
    )

    # Titles & Labels
    plt.title("Protein Design Evaluation: pLDDT vs. pTM", fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("Local Confidence (pLDDT)", fontsize=12)
    plt.ylabel("Global Topology Confidence (pTM)", fontsize=12)

    # Threshold Lines
    plt.axvline(x=70, color="gray", linestyle="--", alpha=0.5, label="pLDDT = 70")
    plt.axhline(y=0.7, color="gray", linestyle=":", alpha=0.5, label="pTM = 0.7")
    
    plt.xlim(min(25, df["pLDDT"].min() - 5), 105)
    plt.ylim(0, 1.05)

    # Colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label("pTM Score", fontsize=10)

    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, dpi=300)
    plt.close()
    
    print(f"'{OUTPUT_IMAGE}' got saved")


if __name__ == "__main__":
    generate_scatter_plot()