import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# --- CONFIGURATION ---
# CSV, die von compare_stability.py erzeugt wurde
CSV_PATH = "/scicore/home/schwede/khan0010/project/tea-leaves-workdir/cath_seq/stability_comparison/lab/all_dG_results_filtered.csv"
OUTPUT_DIR = "/scicore/home/schwede/khan0010/project/tea-leaves-workdir/cath_seq/stability_comparison/lab/per_protein_plots_scale_filter"
# ------------------------------------------------


def sort_key(group_name):
    """Wildtyp zuerst, danach dG_X aufsteigend nach Zahl."""
    if group_name == "Wildtyp":
        return (-1, 0)
    try:
        weight = float(group_name.replace("dG_", ""))
    except ValueError:
        weight = float("inf")
    return (0, weight)


def main():
    if not os.path.isfile(CSV_PATH):
        print(f"Fehler: '{CSV_PATH}' nicht gefunden. Fuehre zuerst compare_stability.py aus.")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(CSV_PATH)

    groups = sorted(df["group"].unique(), key=sort_key)
    design_groups = [g for g in groups if g != "Wildtyp"]

    palette = ["#378ADD", "#639922", "#BA7517", "#D4537E", "#7F77DD"]
    color_map = {"Wildtyp": "#888780"}
    for i, g in enumerate(design_groups):
        color_map[g] = palette[i % len(palette)]

    cath_ids = sorted(df["cath_id"].unique())
    print(f"Prüfe {len(cath_ids)} Proteine auf verfügbare Designs...")

    created_count = 0
    skipped_count = 0

    for cath_id in cath_ids:
        sub = df[df["cath_id"] == cath_id].copy()

        # SKIP-LOGIK: Wenn für dieses Protein KEIN Design vorhanden ist (sondern NUR Wildtyp), überspringen.
        unique_groups = sub["group"].unique()
        if len(unique_groups) == 1 and unique_groups[0] == "Wildtyp":
            skipped_count += 1
            continue

        sub["sort_key"] = sub["group"].apply(sort_key)
        sub = sub.sort_values(["sort_key", "design_id"]).reset_index(drop=True)

        labels = []
        for _, row in sub.iterrows():
            if row["group"] == "Wildtyp":
                labels.append("wildtype")
            else:
                design_id = str(row["design_id"])
                num = design_id.split("designed_")[-1] if "designed_" in design_id else design_id
                labels.append(f"{row['group']}\n#{num}")

        colors = [color_map[g] for g in sub["group"]]

        fig, ax = plt.subplots(figsize=(max(6, len(sub) * 0.5), 5))
        ax.bar(range(len(sub)), sub["dG"], color=colors)
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("\u0394G (kcal/mol)")
        ax.set_title(f"{cath_id}: wildtype vs. designs (0.7 < TM-Score)")
        ax.axhline(0, color="black", linewidth=0.8)

        wt_rows = sub[sub["group"] == "Wildtyp"]
        if not wt_rows.empty:
            wt_value = wt_rows["dG"].iloc[0]
            ax.axhline(wt_value, color=color_map["Wildtyp"], linestyle="--", linewidth=1, alpha=0.7)

        plt.tight_layout()
        out_path = os.path.join(OUTPUT_DIR, f"{cath_id}.png")
        plt.savefig(out_path, dpi=150)
        plt.close(fig)

        created_count += 1

    print(f"\nFertig!")
    print(f"-> {created_count} Plots wurden gespeichert unter '{OUTPUT_DIR}'")
    print(f"-> {skipped_count} Proteine wurden übersprungen (nur Wildtyp ohne Designs vorhanden).")


if __name__ == "__main__":
    main()