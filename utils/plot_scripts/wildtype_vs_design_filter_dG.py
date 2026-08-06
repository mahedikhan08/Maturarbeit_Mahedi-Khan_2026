import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")  # kein Display auf dem Cluster noetig
import matplotlib.pyplot as plt
import pandas as pd

# --- CONFIGURATION ---
WILDTYPE_DIR = "path to stab predictor wildtype output"
DESIGN_DIR = "path to stab predictor design output"
OUTPUT_DIR = "path for output"
SIEVE_BASE_DIR = "path to structures"
# ------------------------------------------------


def load_wildtype_results(wildtype_dir):
    #Liest alle Wildtyp dG_result.json Dateien ein.
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


def load_tm_scores(sieve_base_dir):
    #Liest alle structure_scores.json Dateien aus sieve_base_dir/dG_X/structure_scores.json.
    tm_scores = {}
    pattern = os.path.join(sieve_base_dir, "dG_*", "structure_scores.json")
    
    # Extrahiert den dG_X Ordnernamen aus dem Pfad
    dg_folder_re = re.compile(r"(dG_[0-9.]+)")

    for path in glob.glob(pattern):
        match = dg_folder_re.search(path)
        if not match:
            continue
        group_name = match.group(1)

        try:
            with open(path) as f:
                data = json.load(f)
            for design_id, scores in data.items():
                if "tm_model" in scores:
                    # KEY IST JETZT GRUPPENSPEZIFISCH: ('dG_0.1', 'design_1')
                    tm_scores[(group_name, design_id)] = scores["tm_model"]
        except (json.JSONDecodeError, KeyError, IOError) as e:
            print(f"Warnung: Konnte Scores aus '{path}' nicht laden ({e}).")
    return tm_scores


def load_design_results(design_dir, tm_scores):
    #Liest alle Design dG_result.json Dateien ein, filtert nach tm_model > 0.7
    rows = []
    pattern = os.path.join(design_dir, "dG_*", "*", "*", "dG_result.json")
    
    group_stats = {}
    skipped_no_score = 0

    for path in sorted(glob.glob(pattern)):
        design_folder = os.path.dirname(path)
        design_id = os.path.basename(design_folder)
        
        parts = path.split(os.sep)
        
        # 1. Die echte dG_Gruppe finden
        dg_weight_folder = None
        for part in parts:
            if part.startswith("dG_"):
                dg_weight_folder = part
                break
                
        if dg_weight_folder is None:
            dg_weight_folder = os.path.basename(os.path.dirname(os.path.dirname(design_folder)))

        # 2. Die cath_id sauber basierend auf der dG-Gruppe extrahieren
        try:
            grp_idx = parts.index(dg_weight_folder)
            protein_dir = parts[grp_idx + 1]
            cath_id = protein_dir.split("_")[0]
        except (ValueError, IndexError):
            protein_dir = os.path.basename(os.path.dirname(design_folder))
            cath_id = protein_dir.split("_")[0]
        
        # Match-Key für das gruppenspezifische Dictionary bauen
        match_key = (dg_weight_folder, design_id)

        # Filter prüfen: tm_model score holen
        if match_key not in tm_scores:
            skipped_no_score += 1
            continue
            
        if dg_weight_folder not in group_stats:
            group_stats[dg_weight_folder] = {"total": 0, "filtered": 0}
            
        group_stats[dg_weight_folder]["total"] += 1

        if tm_scores[match_key] <= 0.7:
            group_stats[dg_weight_folder]["filtered"] += 1
            continue

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
            
    filter_percentages = {}
    print("\n--- Filter-Statistiken pro Gruppe ---")
    for grp, stats in sorted(group_stats.items()):
        pct = (stats["filtered"] / stats["total"]) * 100 if stats["total"] > 0 else 0.0
        filter_percentages[grp] = pct
        print(f"Gruppe {grp}: {stats['filtered']} von {stats['total']} aussortiert ({pct:.2f}%)")
        
    if skipped_no_score > 0:
        print(f"-> Designs ohne TM-Score übersprungen: {skipped_no_score} Sequenzen.")
        
    return rows, filter_percentages


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Lade TM-Scores aus structure_scores.json...")
    tm_scores = load_tm_scores(SIEVE_BASE_DIR)
    print(f"-> Insgesamt {len(tm_scores)} TM-Scores geladen.")

    wildtype_rows = load_wildtype_results(WILDTYPE_DIR)
    design_rows, filter_percentages = load_design_results(DESIGN_DIR, tm_scores)

    print(f"\nWildtyp-Ergebnisse gefunden: {len(wildtype_rows)}")
    print(f"Design-Ergebnisse nach Filter gefunden: {len(design_rows)}")

    if not wildtype_rows and not design_rows:
        print("Keine Ergebnisse gefunden, breche ab.")
        return

    df = pd.DataFrame(wildtype_rows + design_rows)

    # CSV mit allen gefilterten Rohdaten speichern
    raw_csv_path = os.path.join(OUTPUT_DIR, "all_dG_results_filtered.csv")
    df.to_csv(raw_csv_path, index=False)
    print(f"-> Gefilterte Rohdaten gespeichert unter '{raw_csv_path}'")

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

    summary_csv_path = os.path.join(OUTPUT_DIR, "summary_dG_by_group_filtered.csv")
    summary.to_csv(summary_csv_path, index=False)
    print(f"-> Zusammenfassung gespeichert unter '{summary_csv_path}'")
    print(summary.to_string(index=False))

    # Balkendiagramm: Mittelwert +/- Standardabweichung pro Gruppe
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = ["#888780" if g == "Wildtyp" else "#378ADD" for g in summary["group"]]
    
    bars = ax.bar(
        summary["group"],
        summary["mean"],
        yerr=summary["std"].fillna(0),
        capsize=5,
        color=colors,
    )
    
    ax.set_ylabel("Average \u0394G (kcal/mol)")
    ax.set_xlabel("Group")
    ax.set_title("Predicted Stability (Filtered: tm_model > 0.7)\nWildtype vs. Designs")
    ax.axhline(0, color="black", linewidth=0.8)
    
    # Text-Labels für die Filter-Prozentzahlen über/unter die Balken schreiben
    for bar, group_name in zip(bars, summary["group"]):
        if group_name == "Wildtyp":
            continue
        
        pct = filter_percentages.get(group_name, 0.0)
        height = bar.get_height()
        
        va_dir = 'bottom' if height >= 0 else 'top'
        offset = 0.5 if height >= 0 else -0.5
        
        ax.text(
            bar.get_x() + bar.get_width() / 2, 
            height + offset, 
            f"Filtered:\n{pct:.1f}%", 
            ha='center', 
            va=va_dir, 
            fontsize=8, 
            weight='bold',
            color='#D9534F'
        )

    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()

    plot_path = os.path.join(OUTPUT_DIR, "wildtype_vs_designs_barplot_filtered.png")
    plt.savefig(plot_path, dpi=150)
    print(f"-> Balkendiagramm mit gruppenweisen Prozentangaben gespeichert unter '{plot_path}'")


if __name__ == "__main__":
    main()