"""
make_cm_figures.py — Generate confusion-matrix-only figures for Home page.

Draws panel (a) from fig_gpcract_perf_L2/L3.png as standalone figures,
without the threshold sweep panel or figure title.

CM values extracted from the original figures (clinical prediction CSVs
are not available on this server).
"""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

FOUT = Path(os.environ.get("SYNERGPCR_OUT", "./output/figures"))

# Style inherited from make_web_figures_v2.py
for k in ["axes.spines.top", "axes.spines.right"]:
    plt.rcParams[k] = False
plt.rcParams.update({"font.size": 11, "svg.fonttype": "none"})

# Exact values from original figures
CM_DATA = {
    "L2": {
        "tn": 102, "fp": 14, "fn": 36, "tp": 243,
        "n": 395, "bacc": 0.8751,
        "title": "G-protein Coupling (L2)",
    },
    "L3": {
        "tn": 19, "fp": 3, "fn": 12, "tp": 49,
        "n": 83, "bacc": 0.8335,
        "title": "β-arrestin Recruitment (L3)",
    },
}


def make_cm_figure(layer, fname):
    d = CM_DATA[layer]
    M = np.array([[d["tn"], d["fp"]],
                  [d["fn"], d["tp"]]])
    LAB = np.array([["TN", "FP"], ["FN", "TP"]])
    correct_mask = np.array([[1, 0], [0, 1]])
    cmap = ListedColormap(["#BBD0F0", "#2E5A9E"])

    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.imshow(correct_mask, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    tot = M.sum()
    for i in range(2):
        for j in range(2):
            txtcol = "white" if correct_mask[i, j] == 1 else "#1a2a4a"
            ax.text(j, i, f"{M[i, j]:,}\n({M[i, j] / tot * 100:.1f}%)\n{LAB[i, j]}",
                    ha="center", va="center", fontsize=12, fontweight="bold", color=txtcol)
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Antagonist", "Agonist"], fontsize=10)
    ax.set_yticklabels(["Antagonist", "Agonist"], fontsize=10, rotation=90, va="center")
    ax.set_xlabel("Predicted MoA", fontsize=11.5, fontweight="bold")
    ax.set_ylabel("True MoA (Clinical)", fontsize=11.5, fontweight="bold")
    for sp in ax.spines.values():
        sp.set_visible(True)
    # No title — will be in HTML card-header
    ax.text(0.5, -0.18,
            f"Clinical MoA BACC = {d['bacc']:.4f}  |  n = {d['n']} approved drugs\n"
            f"(Train: cell-based bioassays; Test: clinically approved drugs)",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=9, color="#333", style="italic")

    fig.savefig(FOUT / fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {fname}")


make_cm_figure("L2", "fig_gpcract_cm_L2.png")
make_cm_figure("L3", "fig_gpcract_cm_L3.png")
print("Done.")
