import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# --- CONFIGURATION ---
# CSV, die von compare_stability.py erzeugt wurde
CSV_PATH = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/stability_comparison/lab/all_dG_results_filtered.csv"
OUTPUT_DIR = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/stability_comparison/lab/violin_plots"
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

    # Gruppen ermitteln und sortieren
    groups = sorted(df["group"].unique(), key=sort_key)
    design_groups = [g for g in groups if g != "Wildtyp"]

    # Farbpalette definieren
    palette = ["#378ADD", "#639922", "#BA7517", "#D4537E", "#7F77DD"]
    color_map = {"Wildtyp": "#888780"}
    for i, g in enumerate(design_groups):
        color_map[g] = palette[i % len(palette)]

    cath_ids = sorted(df["cath_id"].unique())
    print(f"Prüfe {len(cath_ids)} Proteine auf verfügbare Designs...")

    created_count = 0
    skipped_count = 0

    # Seed für reproduzierbares "Jitter" (Streuung) der Punkte
    rng = np.random.default_rng(42)

    for cath_id in cath_ids:
        sub = df[df["cath_id"] == cath_id].copy()

        wt_rows = sub[sub["group"] == "Wildtyp"]
        design_rows = sub[sub["group"] != "Wildtyp"]

        # SKIP-LOGIK: Wenn für dieses Protein KEIN Design vorhanden ist, überspringen.
        if design_rows.empty:
            skipped_count += 1
            continue

        # Designs für die Violinplots gruppieren (sortiert nach dG_X Gewicht)
        unique_design_groups = sorted(design_rows["group"].unique(), key=sort_key)
        
        data_to_plot = []
        labels = []
        colors = []

        for g in unique_design_groups:
            g_data = design_rows[design_rows["group"] == g]["dG"].values
            if len(g_data) > 0:
                data_to_plot.append(g_data)
                labels.append(g)
                colors.append(color_map[g])

        if not data_to_plot:
            skipped_count += 1
            continue

        # Plot erstellen
        fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.5), 5))
        
        # Violinplot generieren
        parts = ax.violinplot(data_to_plot, showmeans=True, showextrema=True)

        # Farben der Violinen anpassen
        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(colors[i])
            pc.set_edgecolor('black')
            pc.set_alpha(0.6) # Etwas transparenter, damit die Punkte besser wirken

        # Farben der Linien in den Violinen anpassen
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmeans'):
            if partname in parts:
                vp = parts[partname]
                vp.set_edgecolor('black')
                vp.set_linewidth(1.2)

        # PUNKTE HINZUFÜGEN (Jitter)
        for i, g_data in enumerate(data_to_plot):
            # i+1 ist die x-Koordinate der jeweiligen Violine (1, 2, 3...)
            x_center = i + 1
            # Ein bisschen Streuung (Jitter) auf der X-Achse generieren
            x_jitter = rng.uniform(-0.08, 0.08, size=len(g_data))
            
            ax.scatter(
                x_center + x_jitter, 
                g_data, 
                color="white",         # Weiße Füllung
                edgecolor="black",     # Schwarzer Rand
                s=20,                  # Größe der Punkte
                zorder=3,              # Z-Order 3 bringt die Punkte in den Vordergrund
                alpha=0.9              # Leichte Transparenz
            )

        # X-Achse beschriften
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylabel("\u0394G (kcal/mol)")
        ax.set_title(f"{cath_id}: Distribution of dG (0.7 < TM-Score)")
        ax.axhline(0, color="black", linewidth=0.8) # Null-Linie

        # Wildtyp als Referenzlinie einzeichnen
        if not wt_rows.empty:
            wt_value = wt_rows["dG"].iloc[0]
            ax.axhline(
                wt_value, 
                color=color_map["Wildtyp"], 
                linestyle="--", 
                linewidth=2, 
                label=f"Wildtype (\u0394G = {wt_value:.2f})"
            )
            ax.legend(loc="best")

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