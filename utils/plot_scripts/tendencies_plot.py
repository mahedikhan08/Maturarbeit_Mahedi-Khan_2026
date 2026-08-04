import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

#CONFIGURATION
CSV_PATH = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/stability_comparison/scale/filtered/all_dG_results_filtered.csv"
OUTPUT_DIR = "/scicore/home/schwede/<username>/project/tea-leaves-workdir/cath_seq/stability_comparison/scale/filtered/comparison_plots"
GROUP_LEFT = "dG_0.01" 
GROUP_RIGHT = "dG_0.1"   
#

def main():
    if not os.path.isfile(CSV_PATH):
        print(f"Error: '{CSV_PATH}' not found.")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = pd.read_csv(CSV_PATH)

    df_xy = df[df["group"].isin([GROUP_LEFT, GROUP_RIGHT])]

    if df_xy.empty:
        print(f"No data found for {GROUP_LEFT} or {GROUP_RIGHT}!")
        return

    # 2. counting designs per cath_id and group
    counts = df_xy.groupby(["cath_id", "group"]).size().unstack(fill_value=0)

    # Safety check, in case a group exists but has 0 entries
    if GROUP_LEFT not in counts.columns: counts[GROUP_LEFT] = 0
    if GROUP_RIGHT not in counts.columns: counts[GROUP_RIGHT] = 0

    # 3. filtering: keep only proteins with at least 1 successful design in either group
    counts = counts[(counts[GROUP_LEFT] > 0) | (counts[GROUP_RIGHT] > 0)].copy()
    
    print(f"Creating trend plot for {len(counts)} proteins...")

    # 4. calculating trend (between 0.0 and 1.0)
    # 0.0 = 100% GROUP_LEFT | 1.0 = 100% GROUP_RIGHT | 0.5 = equal
    counts["total"] = counts[GROUP_LEFT] + counts[GROUP_RIGHT]
    counts["ratio"] = counts[GROUP_RIGHT] / counts["total"]

    # 5. sort for the y-axis (z.B. alphabetically after cath_id)
    counts = counts.sort_index(ascending=False)
    
    cath_ids = counts.index.tolist()
    ratios = counts["ratio"].tolist()
    y_positions = range(len(cath_ids))

    y_labels = []
    for cath_id in cath_ids:
        c_left = counts.loc[cath_id, GROUP_LEFT]
        c_right = counts.loc[cath_id, GROUP_RIGHT]
        short_id = cath_id.split("_")[0]
        y_labels.append(f"[{c_left} | {c_right}] {short_id}")

    # 6. create Plot
    fig_height = max(5.5, len(cath_ids) * 0.28 + 0.5)
    fig, ax = plt.subplots(figsize=(8, fig_height))

    # creating and coloring the scatter points based on the ratio
    scatter = ax.scatter(
        ratios, 
        y_positions, 
        c=ratios, 
        cmap="coolwarm",  # blue = left, red = right
        s=60,             # size of the points
        edgecolor="black",
        linewidth=0.5,
        alpha=0.85,
        zorder=3
    )

    # 7. axis & layout
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=8)
    
    # x-axis: 0.0 = 100% GROUP_LEFT, 1.0 = 100% GROUP_RIGHT
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels([
        f"100%\n{GROUP_LEFT}", 
        "25/75", 
        "50/50\n", 
        "75/25", 
        f"100%\n{GROUP_RIGHT}"
    ], fontsize=9)

    # vertical line at 0.5 to indicate equal distribution
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=1.2, zorder=1)

    # light grid lines for better readability
    ax.xaxis.grid(True, linestyle=":", alpha=0.5, zorder=0)

    # title
    ax.set_title(
        f"Trend of the successfull designs\n{GROUP_LEFT} (left) vs {GROUP_RIGHT} (right)", 
        fontsize=12, weight="bold", pad=15
    )

    # total counts of successful designs for each group
    total_left = counts[GROUP_LEFT].sum()
    total_right = counts[GROUP_RIGHT].sum()
    fig.text(
        0.5, 0.015, 
        f"Total successful sequences: {total_left} ({GROUP_LEFT}) | {total_right} ({GROUP_RIGHT})", 
        ha="center", fontsize=10, weight="bold", 
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round,pad=0.5')
    )

    plt.tight_layout(rect=[0, 0.04, 1, 0.98])
    
    out_path = os.path.join(OUTPUT_DIR, f"tendency_plot_{GROUP_LEFT}_vs_{GROUP_RIGHT}.png")
    plt.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"-> saved trend plot to '{out_path}'")

if __name__ == "__main__":
    main()