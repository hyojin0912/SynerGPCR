"""
make_web_figures_v2.py - STEP 5g-123fix + 5g-4.

FIG 1-3 : web figures re-rendered (centered x-tick main labels + grey "(n = ...)" lines,
          null95 line removed). FIG 4-5 : GPCRact clinical performance (confusion matrix +
          confidence threshold sweep) for L2 / L3.

Style/palette/rcParams inherited from Output/NAR_DBI/Code/webexport_step7_figures.py.
All numbers recomputed from CSV sources (final_export/ for FIG1-3; clinical prediction CSVs
for FIG4-5). FINAL_STATS asserts must pass before the run is trusted. No null95 line drawn.

Python 3.9 / matplotlib Agg / pandas 2.2.2.
"""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score

BASE = Path(os.environ.get("SYNERGPCR_OUT", "./output/figures"))
SRC = BASE / "final_export"
FOUT = SRC / "figures"
FOUT.mkdir(parents=True, exist_ok=True)

# clinical prediction CSVs that reproduce FINAL_STATS section G (verified in PART 0)
GP_ROOT = Path(os.environ.get("SYNERGPCR_OUT", "./output/figures"))
PERF_SRC = {
    "L2": (GP_ROOT / "L2_v2_clinical/ablation/withoutInactives/predictions_test.csv", 0.40),
    "L3": (GP_ROOT / "L3_v1_clinical_egnn/predictions_test.csv", 0.97),
}

# -- inherited style ---------------------------------------------------
for k in ["axes.spines.top", "axes.spines.right"]:
    plt.rcParams[k] = False
plt.rcParams.update({"font.size": 11, "svg.fonttype": "none"})
COL_BASE, COL_GP, COL_CEIL, COL_FULL = "#8C92AC", "#2E8B57", "#4169E1", "#191970"
COL_SINGLE = "#B8BCC8"
COL_N = "#888888"          # grey for "(n = ...)" annotations
LAYER_NAME = {"L1": "Binding", "L2": "G-protein", "L3": "β-arrestin", "L4": "Reporter"}

TEXTS = []


def T(s):
    TEXTS.append(str(s))
    return s


def blend(ax):
    """x in data coords, y in axes coords."""
    return mtransforms.blended_transform_factory(ax.transData, ax.transAxes)


# -- sources (final_export/) -------------------------------------------
g3 = pd.read_csv(SRC / "gpcract_3point_si.csv")
sl = pd.read_csv(SRC / "single_layer_table.csv")
syn = pd.read_csv(SRC / "synergy_stats.csv")

plotted = {}


def sig_bracket(ax, x1, x2, y, h, text):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, color="#333", zorder=4)
    ax.text((x1 + x2) / 2.0, y + h, T(text), ha="center", va="bottom",
            fontsize=9, color="#333", zorder=4)


# ============================ FIGURE 1 ================================
fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), gridspec_kw={"wspace": 0.26})
for ax, layer, tag in [(axes[0], "L2", "a"), (axes[1], "L3", "b")]:
    d = g3[g3.layer == layer].set_index("point")
    order = ["baseline", "gpcract", "ceiling"]
    cols = [COL_BASE, COL_GP, COL_CEIL]
    hatch = ["", "///", ""]
    si = [float(d.loc[p, "SI"]) for p in order]
    elo = [si[i] - float(d.loc[p, "SI_ci_lo"]) for i, p in enumerate(order)]
    ehi = [float(d.loc[p, "SI_ci_hi"]) - si[i] for i, p in enumerate(order)]
    pool_ns = [int(d.loc[p, "n_pool"]) for p in order]
    plotted[f"F1_{layer}"] = [(p, round(si[i], 4)) for i, p in enumerate(order)]
    ymax = max(float(d.loc[p, "SI_ci_hi"]) for p in order) * 1.35
    for i, p in enumerate(order):
        ax.bar(i, si[i], yerr=[[elo[i]], [ehi[i]]], width=0.64, color=cols[i], alpha=0.88,
               edgecolor="black", linewidth=1.2, hatch=hatch[i], capsize=5, ecolor="#333", zorder=2)
        ax.text(i, si[i] / 2, f"{si[i]:.2f}", ha="center", va="center", color="white",
                fontsize=11, fontweight="bold", zorder=3)
    rec = float(d.loc["gpcract", "recovery_pct"])
    gp_hi = float(d.loc["gpcract", "SI_ci_hi"])
    ax.text(1, gp_hi + ymax * 0.03, T(f"recovery {rec:.1f}%"),
            ha="center", fontsize=9, color=COL_GP, fontweight="bold")
    sig_bracket(ax, 0, 1, gp_hi + ymax * 0.12, ymax * 0.03, "p < 0.0001")
    # centered main labels + grey "(n = ...)" line, drawn via blended transform
    not_assayed = "(L2 not assayed)" if layer == "L2" else "(L3 not assayed)"
    predicted = "(L2 predicted)" if layer == "L2" else "(L3 predicted)"
    ceil_combo = "(L1+L2+L4)" if layer == "L2" else "(L1+L3+L4)"
    main_labels = [T(f"L1+L4 only\n{not_assayed}"), T(f"GPCRact\n{predicted}"),
                   T(f"Exp. ceiling\n{ceil_combo}")]
    ax.set_xticks(range(3))
    ax.set_xticklabels([""] * 3)
    ax.tick_params(axis="x", length=0)
    for xi, (ml, nv) in enumerate(zip(main_labels, pool_ns)):
        ax.text(xi, -0.06, ml, ha="center", va="top", transform=blend(ax), fontsize=10.0)
        ax.text(xi, -0.19, T(f"(n = {nv:,})"), ha="center", va="top",
                transform=blend(ax), fontsize=9, color=COL_N)
    ax.set_ylim(0, ymax)
    ax.yaxis.grid(True, ls="--", alpha=0.35, zorder=0)
    ax.set_axisbelow(True)
    title = "G-protein coupling (L2)" if layer == "L2" else "β-arrestin recruitment (L3)"
    ax.set_title(T(title), fontsize=11.5, fontweight="bold", pad=10)
    ax.text(-0.08, 1.05, tag, transform=ax.transAxes, fontsize=15, fontweight="bold")
axes[0].set_ylabel(T("Bliss Synergy Index (SI)"), fontsize=12)
fig.suptitle(T("Assay Chain Completion via AI Prediction"), fontsize=13, fontweight="bold")
fig.subplots_adjust(bottom=0.22)
fig.savefig(FOUT / "fig_gpcract_3point_v2.png", dpi=600, bbox_inches="tight")
plt.close(fig)
print("wrote fig_gpcract_3point_v2.png")

# ============================ FIGURE 2 ================================
fig, ax = plt.subplots(figsize=(6.8, 4.9))
sl2 = sl.set_index("layer").loc[["L1", "L2", "L3", "L4"]]
rate = sl2.approval_rate_pct.values
elo = rate - sl2.wilson_ci_lo_pct.values
ehi = sl2.wilson_ci_hi_pct.values - rate
ncomp = sl2.n_compounds.astype(int).values
plotted["F2"] = [(l, round(rate[i], 4)) for i, l in enumerate(["L1", "L2", "L3", "L4"])]
ax.bar(range(4), rate, width=0.62, color=COL_BASE, alpha=0.9, edgecolor="black",
       linewidth=1.2, yerr=[elo, ehi], capsize=5, ecolor="#333", zorder=2)
for i in range(4):
    ax.text(i, sl2.wilson_ci_hi_pct.values[i] + 0.02, f"{rate[i]:.3f}%", ha="center", fontsize=9)
ax.axhline(1.0, ls=":", color="#d93025", lw=1.3)
ax.text(3.45, 1.0, T(" 1% reference"), color="#d93025", va="bottom", ha="right", fontsize=9)
layer_info = [("L1", "Binding"), ("L2", "G-protein\ncoupling"),
              ("L3", "β-arrestin\nrecruitment"), ("L4", "Reporter gene")]
ax.set_xticks(range(4))
ax.set_xticklabels([""] * 4)
ax.tick_params(axis="x", length=0)
for xi, ((layer, name), nv) in enumerate(zip(layer_info, ncomp)):
    ax.text(xi, -0.06, T(f"{layer}\n{name}"), ha="center", va="top",
            transform=blend(ax), fontsize=10.0)
    ax.text(xi, -0.22, T(f"(n = {nv:,})"), ha="center", va="top",
            transform=blend(ax), fontsize=9, color=COL_N)
ax.set_ylabel(T("Clinical approval rate (%)"), fontsize=12)
ax.set_ylim(0, 1.2)
ax.yaxis.grid(True, ls="--", alpha=0.35, zorder=0)
ax.set_axisbelow(True)
ax.set_title(T("Single-layer Clinical Approval Rate"), fontsize=11.5, fontweight="bold", pad=8)
fig.subplots_adjust(bottom=0.25)
fig.savefig(FOUT / "fig_single_layer_v2.png", dpi=600, bbox_inches="tight")
plt.close(fig)
print("wrote fig_single_layer_v2.png")

# ============================ FIGURE 3 ================================
s = syn[syn.bliss_SI.notna()].copy().sort_values("bliss_SI").reset_index(drop=True)
fig, ax = plt.subplots(figsize=(12, 5.2))
colors, hatches = [], []
for cp, si in zip(s.chain_pattern, s.bliss_SI):
    if cp == "Full_chain":
        colors.append(COL_FULL)
        hatches.append("///")
    elif si == 1.0:
        colors.append(COL_SINGLE)
        hatches.append("")
    else:
        colors.append(COL_CEIL)
        hatches.append("")
elo = (s.bliss_SI - s.bliss_SI_ci_lo).clip(lower=0).values
ehi = (s.bliss_SI_ci_hi - s.bliss_SI).values
plotted["F3"] = [(cp, round(float(v), 4)) for cp, v in zip(s.chain_pattern, s.bliss_SI)]
ymax = float(s.bliss_SI_ci_hi.max()) * 1.22
for i in range(len(s)):
    cihi = s.bliss_SI_ci_hi.values[i]
    ax.bar(i, s.bliss_SI.values[i], width=0.72, color=colors[i], alpha=0.9, edgecolor="black",
           linewidth=1.0, hatch=hatches[i], yerr=[[elo[i]], [ehi[i]]], capsize=3, ecolor="#444",
           error_kw={"lw": 0.9}, zorder=2)
    ax.text(i, cihi + 0.25, f"{s.bliss_SI.values[i]:.2f}",
            ha="center", fontsize=8, color="#1a1a1a")
    ax.text(i, cihi + 0.95, T(f"(n={int(s.n_compounds.values[i]):,})"),
            ha="center", fontsize=6.5, color=COL_N)
ax.axhline(1.0, ls=":", color="#777", lw=1.0, zorder=1)
ax.set_xlim(-1.25, len(s) - 0.3)
ax.annotate(T("SI = 1\n(independence)"), xy=(-0.55, 1.0), xytext=(-1.20, 4.6),
            va="center", ha="left", fontsize=8, color="#555",
            arrowprops=dict(arrowstyle="->", color="#999", lw=0.8), zorder=5)
ax.set_xticks(range(len(s)))
ax.set_xticklabels([T(cp) for cp in s.chain_pattern], rotation=35, ha="center", fontsize=8)
ax.tick_params(axis="x", pad=8)
ax.set_ylabel(T("Bliss Synergy Index (SI)"), fontsize=12)
ax.set_ylim(0, ymax)
ax.yaxis.grid(True, ls="--", alpha=0.35, zorder=0)
ax.set_axisbelow(True)
ax.set_title(T("Assay-chain Synergy Index"), fontsize=11.5, fontweight="bold", pad=8)
ax.legend(handles=[Patch(facecolor=COL_SINGLE, edgecolor="black", label=T("single layer")),
                   Patch(facecolor=COL_CEIL, edgecolor="black", label=T("multi-layer chain")),
                   Patch(facecolor=COL_FULL, edgecolor="black", hatch="///", label=T("full chain"))],
          fontsize=8.5, loc="upper left", frameon=False)
fig.savefig(FOUT / "fig_chain_synergy_v2.png", dpi=600, bbox_inches="tight")
plt.close(fig)
print("wrote fig_chain_synergy_v2.png")


# ===================== FIG 4 / 5 : clinical performance ===============
def load_perf(layer):
    """Return (df, confusion dict, bacc, n) recomputed from the clinical test CSV."""
    path, _op = PERF_SRC[layer]
    d = pd.read_csv(path)
    d = d[d.Final_Label.isin([1, 2])].copy()          # all true-active approved pairs
    true = (d.Final_Label == 2).astype(int).values     # 1 = agonist, 0 = antagonist
    pred = (d.ago_prob.values >= d.ant_prob.values).astype(int)
    bacc = balanced_accuracy_score(true, pred)
    cm = {
        "tp": int(((true == 1) & (pred == 1)).sum()),   # agonist -> agonist
        "fn": int(((true == 1) & (pred == 0)).sum()),   # agonist -> antagonist
        "fp": int(((true == 0) & (pred == 1)).sum()),   # antagonist -> agonist
        "tn": int(((true == 0) & (pred == 0)).sum()),   # antagonist -> antagonist
    }
    d["_conf"] = d.binding_prob.values * np.maximum(d.ago_prob.values, d.ant_prob.values)
    d["_correct"] = (pred == true)
    return d, cm, bacc, len(d)


def make_perf_figure(layer, oracle_bacc, op, layer_title, fname):
    d, cm, bacc, n = load_perf(layer)
    plotted[f"PERF_{layer}"] = {"cm": cm, "bacc": round(bacc, 4), "n": n}

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.4, 4.8),
                                   gridspec_kw={"wspace": 0.34, "width_ratios": [1, 1.25]})

    # ---- Panel A : MoA confusion matrix (rows=True, cols=Pred) ----
    # layout order [Antagonist, Agonist]; diagonal = correct (dark), off-diag = error (light)
    M = np.array([[cm["tn"], cm["fp"]],
                  [cm["fn"], cm["tp"]]])
    LAB = np.array([["TN", "FP"], ["FN", "TP"]])
    correct_mask = np.array([[1, 0], [0, 1]])  # 1 = dark (correct)
    cmap = ListedColormap(["#BBD0F0", "#2E5A9E"])  # light error / dark correct
    axA.imshow(correct_mask, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    tot = M.sum()
    for i in range(2):
        for j in range(2):
            txtcol = "white" if correct_mask[i, j] == 1 else "#1a2a4a"
            axA.text(j, i, f"{M[i, j]:,}\n({M[i, j] / tot * 100:.1f}%)\n{LAB[i, j]}",
                     ha="center", va="center", fontsize=11, fontweight="bold", color=txtcol)
    axA.set_xticks([0, 1])
    axA.set_yticks([0, 1])
    axA.set_xticklabels([T("Antagonist"), T("Agonist")], fontsize=10)
    axA.set_yticklabels([T("Antagonist"), T("Agonist")], fontsize=10, rotation=90, va="center")
    axA.set_xlabel(T("Predicted MoA"), fontsize=11.5, fontweight="bold")
    axA.set_ylabel(T("True MoA (Clinical)"), fontsize=11.5, fontweight="bold")
    for sp in axA.spines.values():
        sp.set_visible(True)
    axA.set_title(T("(a) MoA confusion matrix"), fontsize=11.5, fontweight="bold", pad=8)
    axA.text(0.5, -0.14,
             T(f"Clinical MoA BACC = {bacc:.4f}  |  n = {n} approved drugs\n"
               f"(Train: cell-based bioassays; Test: clinically approved drugs)"),
             transform=axA.transAxes, ha="center", va="top",
             fontsize=9, color="#333", style="italic")

    # ---- Panel B : confidence threshold sweep ----
    thr = np.linspace(0.0, 1.0, 101)
    conf = d["_conf"].values
    corr = d["_correct"].values
    N = len(d)
    recovery, precision = [], []
    for t in thr:
        sel = conf >= t
        recovery.append(sel.sum() / N)
        precision.append(corr[sel].mean() if sel.sum() > 0 else np.nan)
    recovery = np.array(recovery)
    precision = np.array(precision)
    plotted[f"PERF_{layer}"]["op_recovery"] = round(float(recovery[np.argmin(np.abs(thr - op))]), 4)
    plotted[f"PERF_{layer}"]["op_precision"] = round(float(precision[np.argmin(np.abs(thr - op))]), 4)
    axB.plot(thr, precision, color=COL_CEIL, lw=2.0, label=T("MoA precision"))
    axB.plot(thr, recovery, color="#E8820C", lw=2.0, ls="--", label=T("Recovery rate"))
    axB.axvline(op, color=COL_GP, lw=1.6, ls=":", zorder=1)
    if op < 0.5:        # L2: text to the right of the line
        txt_x, txt_ha = op + 0.02, "left"
    else:               # L3: text to the left of the line
        txt_x, txt_ha = op - 0.02, "right"
    axB.text(txt_x, 0.62, T(f"Selected threshold\nt = {op:.2f}"),
             color=COL_GP, ha=txt_ha, va="center", fontsize=8.5, fontweight="bold")
    axB.set_xlabel(T("Confidence threshold"), fontsize=11.5)
    axB.set_ylabel(T("Rate"), fontsize=11.5)
    axB.set_xlim(0, 1.0)
    axB.set_ylim(0, 1.12)
    axB.yaxis.grid(True, ls="--", alpha=0.35, zorder=0)
    axB.set_axisbelow(True)
    axB.set_title(T(f"(b) Confidence threshold — {layer_title}"),
                  fontsize=11.5, fontweight="bold", pad=8)
    axB.legend(handles=[plt.Line2D([], [], color=COL_CEIL, lw=2.0, label=T("MoA precision")),
                        plt.Line2D([], [], color="#E8820C", lw=2.0, ls="--", label=T("Recovery rate")),
                        plt.Line2D([], [], color=COL_GP, lw=1.6, ls=":",
                                   label=T(f"Selected threshold (t = {op:.2f})"))],
               fontsize=8.5, loc="lower left", frameon=True, framealpha=0.92,
               edgecolor="#CBD5E1", bbox_to_anchor=(0.01, 0.01))

    fig.savefig(FOUT / fname, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {fname}")
    return bacc, n, cm


bacc_L2, n_L2, cm_L2 = make_perf_figure("L2", 0.875, 0.40, "G-protein Coupling (L2)",
                                        "fig_gpcract_perf_L2.png")
bacc_L3, n_L3, cm_L3 = make_perf_figure("L3", 0.834, 0.97, "β-arrestin Recruitment (L3)",
                                        "fig_gpcract_perf_L3.png")

# ============================ ORACLE ASSERTS ==========================
def gv(layer, point, col):
    return float(g3[(g3.layer == layer) & (g3.point == point)].iloc[0][col])


assert abs(gv("L2", "gpcract", "SI") - 7.74) < 0.05
assert abs(gv("L2", "ceiling", "SI") - 9.71) < 0.05
assert abs(gv("L2", "gpcract", "recovery_pct") - 79.7) < 0.5
assert abs(gv("L3", "gpcract", "SI") - 10.35) < 0.05
assert abs(gv("L3", "ceiling", "SI") - 13.01) < 0.05
assert abs(gv("L3", "gpcract", "recovery_pct") - 79.6) < 0.5
sl_rate = sl.set_index("layer").approval_rate_pct
assert abs(sl_rate["L1"] - 0.166) < 0.01
assert abs(sl_rate["L2"] - 0.211) < 0.01
assert abs(sl_rate["L3"] - 0.553) < 0.01
assert abs(sl_rate["L4"] - 0.284) < 0.01
syn_si = syn.set_index("chain_pattern").bliss_SI
assert abs(float(syn_si["Full_chain"]) - 13.45) < 0.05
assert abs(float(syn_si["L1+L2+L4"]) - 9.71) < 0.05
# FINAL_STATS section G (clinical BACC) — recomputed from CSV
assert abs(bacc_L2 - 0.875) < 0.005, f"L2 BACC {bacc_L2}"
assert abs(bacc_L3 - 0.834) < 0.005, f"L3 BACC {bacc_L3}"
print("\nAll FINAL_STATS section D + G oracle asserts PASSED.")

# ============================ SELF-CHECK =============================
print("\n" + "=" * 72)
print("SELF-CHECK: plotted values vs CSV source")
print("=" * 72)
ok = True


def chk(cond, msg):
    global ok
    ok = ok and cond
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")


for layer in ["L2", "L3"]:
    for p, si in plotted[f"F1_{layer}"]:
        row = g3[(g3.layer == layer) & (g3.point == p)].iloc[0]
        chk(round(float(row.SI), 4) == si, f"F1 {layer}/{p}: SI={si}")
for l, r in plotted["F2"]:
    row = sl[sl.layer == l].iloc[0]
    chk(round(float(row.approval_rate_pct), 4) == r, f"F2 {l}: rate={r}%")
for cp, v in plotted["F3"]:
    row = syn[syn.chain_pattern == cp].iloc[0]
    chk(round(float(row.bliss_SI), 4) == v, f"F3 {cp}: SI={v}")
for layer in ["L2", "L3"]:
    pf = plotted[f"PERF_{layer}"]
    s = pf["cm"]["tp"] + pf["cm"]["fn"] + pf["cm"]["fp"] + pf["cm"]["tn"]
    chk(s == pf["n"], f"PERF {layer}: confusion sum {s} == n {pf['n']}")
    print(f"      PERF {layer}: BACC={pf['bacc']} cm={pf['cm']} "
          f"op_recovery={pf['op_recovery']} op_precision={pf['op_precision']}")

# ============================ BANNED / STALE SCAN ====================
print("\n" + "=" * 72)
print("BANNED-TERM / STALE scan over all figure text")
print("=" * 72)
BANNED = ["permutation", "dsav", "actr", "forward", "backward"]
STALE = ["8.26", "8.96", "7.81", "7.67", "3.64", "3.699", "260938", "260,938",
         "481240", "481,240", "11428", "11,428"]
allbad = []
for t in TEXTS:
    low = t.lower()
    for b in BANNED:
        if b in low:
            allbad.append((b, t))
    for sv in STALE:
        if sv in t:
            allbad.append((sv, t))
chk(len(allbad) == 0, f"no banned term / stale value in figure text (hits={allbad})")
print(f"  scanned {len(TEXTS)} text strings")
print("\nRESULT:", "ALL PASS" if ok else "*** FAILURES ABOVE ***")
assert ok, "self-check or banned-scan failed"
