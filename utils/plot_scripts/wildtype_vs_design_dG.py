import glob
import json
import os

import matplotlib
matplotlib.use("Agg")  # kein Display auf dem Cluster noetig
import matplotlib.pyplot as plt
import pandas as pd

# --- CONFIGURATION ---
WILDTYPE_DIR = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/wildtype_esm3dg_predictions"
DESIGN_DIR = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/sieve_output/esm3dg_predictions_scale"
OUTPUT_DIR = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/stability_comparison/scale"
# ------------------------------------------------


def load_wildtype_results(wildtype_dir):
    """Liest alle Wildtyp dG_result.json Dateien ein. Gibt Liste von Dicts zurück."""
    rows = []
    pattern = os.path.join(wildtype_dir, "*", "dG_result.json")
    for path in sorted(glob.glob(pattern)):
        cath_id = os.path.basename(os.path.dirname(path))
        try:
            with open(path) as f:
                data = json.load(f)
            rows.append({
                "group": "Wildtyp",
                "cath_id": cath_id,
                "design_id": None,
                "dG": data["ensemble_dG"],
            })
        except (json.JSONDecodeError, KeyError) as e:
            print(f"Warnung: konnte '{path}' nicht lesen ({e}), ueberspringe.")
    return rows


def load_design_results(design_dir):
    """
    Liest alle Design dG_result.json Dateien ein. Erwartete Struktur:
    design_dir/dG_X/<protein_dir>/<design_id>/dG_result.json
    Die CATH-ID wird aus dem ersten Unterstrich-Teil von protein_dir extrahiert,
    z.B. '1a04A02_129-215_cath_4_4_0' -> '1a04A02', damit sie mit den
    Wildtyp-IDs (Dateiname ohne Endung, z.B. '1a04A02') uebereinstimmt.
    """
    rows = []
    pattern = os.path.join(design_dir, "dG_*", "*", "*", "dG_result.json")
    for path in sorted(glob.glob(pattern)):
        design_folder = os.path.dirname(path)
        design_id = os.path.basename(design_folder)
        protein_dir = os.path.basename(os.path.dirname(design_folder))
        dg_weight_folder = os.path.basename(os.path.dirname(os.path.dirname(design_folder)))
        cath_id = protein_dir.split("_")[0]
        try:
            with open(path) as f:
                data = json.load(f)
            rows.append({
                "group": dg_weight_folder,
                "cath_id": cath_id,
                "design_id": design_id,
                "dG": data["ensemble_dG"],
            })
        except (json.JSONDecodeError, KeyError) as e:
            print(f"Warnung: konnte '{path}' nicht lesen ({e}), ueberspringe.")
    return rows


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    wildtype_rows = load_wildtype_results(WILDTYPE_DIR)
    design_rows = load_design_results(DESIGN_DIR)

    print(f"Wildtyp-Ergebnisse gefunden: {len(wildtype_rows)}")
    print(f"Design-Ergebnisse gefunden: {len(design_rows)}")

    if not wildtype_rows and not design_rows:
        print("Keine Ergebnisse gefunden, breche ab.")
        return

    df = pd.DataFrame(wildtype_rows + design_rows)

    # CSV mit allen Rohdaten speichern (fuer weitere eigene Analyse)
    raw_csv_path = os.path.join(OUTPUT_DIR, "all_dG_results.csv")
    df.to_csv(raw_csv_path, index=False)
    print(f"-> Rohdaten gespeichert unter '{raw_csv_path}'")

    # Gruppen sortieren: Wildtyp zuerst, danach dG_X aufsteigend nach Zahl
    def sort_key(group_name):
        if group_name == "Wildtyp":
            return (-1, 0)
        try:
            weight = float(group_name.replace("dG_", ""))
        except ValueError:
            weight = float("inf")
        return (0, weight)

    summary = df.groupby("group")["dG"].agg(["mean", "std", "count"]).reset_index()
    summary["sort_key"] = summary["group"].apply(sort_key)
    summary = summary.sort_values("sort_key").drop(columns="sort_key")

    summary_csv_path = os.path.join(OUTPUT_DIR, "summary_dG_by_group.csv")
    summary.to_csv(summary_csv_path, index=False)
    print(f"-> Zusammenfassung gespeichert unter '{summary_csv_path}'")
    print(summary.to_string(index=False))


    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#888780" if g == "Wildtyp" else "#378ADD" for g in summary["group"]]
    ax.bar(
        summary["group"],
        summary["mean"],
        yerr=summary["std"],
        capsize=5,
        color=colors,
    )
    ax.set_ylabel("Average \u0394G (kcal/mol)")
    ax.set_xlabel("Group")
    ax.set_title("Predicted Stability: Wildtype vs. Designs")
    ax.axhline(0, color="black", linewidth=0.8)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()

    plot_path = os.path.join(OUTPUT_DIR, "wildtype_vs_designs_barplot.png")
    plt.savefig(plot_path, dpi=150)
    print(f"-> Balkendiagramm gespeichert unter '{plot_path}'")


if __name__ == "__main__":
    main()
