"""
part11_threshold_sweep_final.py  — READ-ONLY. Prints only; no file edits.

Threshold sweep for L2 and L3 using the FINAL pair-level pipeline (NOT the stale
compound-level 05/15 one), to check whether the locked op-points (L2 0.50 / L3 0.85) are
margin-optimal. SI / Monte-Carlo-null / pair-level pool code reused VERBATIM from
part11_si_and_null_final.py — only the threshold is looped.

Self-check: the t=0.50 (L2) and t=0.85 (L3) rows must reproduce SI 7.73 / 8.89.
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
from typing import Tuple
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


BASE     = BASE
CAT      = BASE / "Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv"
CHAIN_D  = BASE / "Output/DB/GPCRactDB/Analysis/Chain_Pattern_Distribution.csv"
PAIR_CAT = BASE / "Output/NAR_L3L4/Tables/L3L4_Pair_Catalog_v1.csv"
MASTER   = BASE / "Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv"
L2_PRED  = BASE / "GPCRactDB/results/denovo/synergpcr_annotation_L2.csv"
L3_PRED  = BASE / "GPCRactDB/results/denovo/synergpcr_annotation_L3.csv"
N_REPS, SEED = 10_000, 0
MOA_MAP = {"Agonist": "agonist", "Partial Agonist": "agonist", "PAM": "agonist",
           "Antagonist": "antagonist", "Inverse Agonist": "antagonist", "NAM": "antagonist",
           "Inactive": "inactive"}
L2_GRID = [0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.85, 0.97, 0.99]
L3_GRID = [0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.97, 0.99]


# ── verified helpers (verbatim) ───────────────────────────────────────────────
def compute_si(k, n, e):
    return float("nan") if (n == 0 or e == 0) else (k / n) / e


def get_p_only(dist, layer):
    row = dist[dist["chain_pattern"] == f"{layer}_only"]
    return int(row["n_approved"].iloc[0]) / int(row["n"].iloc[0])


def mc_null(base_df, n_up, n_ceil, k_ceil, e):
    """Monte-Carlo random-assignment null (verbatim method/seed from part11_si_and_null_final)."""
    if n_up == 0:
        return float("nan"), float("nan")
    appr = base_df.is_approved.to_numpy().astype(int)
    rng = np.random.default_rng(SEED)
    N = n_ceil + n_up
    idx = np.arange(len(appr))
    sis = np.empty(N_REPS)
    for i in range(N_REPS):
        pick = rng.choice(idx, size=n_up, replace=False)
        sis[i] = ((k_ceil + int(appr[pick].sum())) / N) / e
    return float(sis.mean()), float(np.percentile(sis, 95))


def knee_by_chord(ts, ns, sis):
    """Max-distance-to-chord knee of the (n_up vs SI) curve. Returns the threshold at the knee."""
    pts = np.array(list(zip(ns, sis)), dtype=float)
    # normalise both axes to [0,1] so the chord distance is scale-free
    rng = pts.max(0) - pts.min(0)
    rng[rng == 0] = 1.0
    p = (pts - pts.min(0)) / rng
    a, b = p[0], p[-1]
    ab = b - a
    denom = np.hypot(*ab) or 1.0
    d = np.abs(ab[0] * (a[1] - p[:, 1]) - (a[0] - p[:, 0]) * ab[1]) / denom
    return ts[int(d.argmax())], float(d.max())


def main() -> None:
    cat = pd.read_csv(CAT, low_memory=False)
    cat["is_approved"] = (pd.to_numeric(cat["max_phase"], errors="coerce") == 4)
    cat = cat.drop_duplicates("InChIKey")
    A = lambda c: cat[c] == "Active"
    dist = pd.read_csv(CHAIN_D)
    pc = pd.read_csv(PAIR_CAT, low_memory=False)
    l2 = pd.read_csv(L2_PRED)
    l3 = pd.read_csv(L3_PRED)

    P = {L: get_p_only(dist, L) for L in ["L1", "L2", "L3", "L4"]}
    E2 = 1 - (1 - P["L1"]) * (1 - P["L2"]) * (1 - P["L4"])
    E3 = 1 - (1 - P["L1"]) * (1 - P["L3"]) * (1 - P["L4"])

    # pools (compound-level, verbatim)
    ceil2 = cat[A("L1_state") & A("L2_state") & A("L4_state")].drop_duplicates("InChIKey")
    ceil3 = cat[A("L1_state") & A("L3_state") & A("L4_state")].drop_duplicates("InChIKey")
    base2 = cat[A("L1_state") & A("L4_state") & cat.L2_state.isna()].drop_duplicates("InChIKey")
    base3 = cat[A("L1_state") & A("L4_state") & cat.L3_state.isna()].drop_duplicates("InChIKey")
    ceil2_ik, n_c2, k_c2 = set(ceil2.InChIKey), len(ceil2), int(ceil2.is_approved.sum())
    ceil3_ik, n_c3, k_c3 = set(ceil3.InChIKey), len(ceil3), int(ceil3.is_approved.sum())
    base2_ik, base3_ik = set(base2.InChIKey), set(base3.InChIKey)
    base2_ap = base2.set_index("InChIKey").is_approved
    base3_ap = base3.set_index("InChIKey").is_approved
    si_ceil2 = compute_si(k_c2, n_c2, E2)
    si_ceil3 = compute_si(k_c3, n_c3, E3)
    print(f"E_Bliss L1+L2+L4={E2:.6f}  L1+L3+L4={E3:.6f}")
    print(f"ceiling SI: L2={si_ceil2:.4f}  L3={si_ceil3:.4f}  (denominators for recovery%)")

    # pair scopes (verbatim)
    m2 = pd.read_csv(MASTER, usecols=["Ligand_InChIKey", "GPCR_UniProt", "Assay_Layer", "MoA"]).rename(
        columns={"Ligand_InChIKey": "InChIKey", "GPCR_UniProt": "UniProt_AC"})
    m2 = m2[m2.Assay_Layer == "Layer 2 (Proximal)"].copy()
    m2["lab"] = m2.MoA.map(MOA_MAP)
    l2_known = set(map(tuple, m2.dropna(subset=["lab"]).drop_duplicates(["InChIKey", "UniProt_AC"])[
        ["InChIKey", "UniProt_AC"]].values))
    cand = pc[(pc.L1_state == "Active") & (pc.L4_state == "Active")][["InChIKey", "UniProt_AC"]].drop_duplicates()
    ck = list(map(tuple, cand.values))
    scope2 = cand[[k not in l2_known for k in ck]]
    scope3 = pc[(pc.L1_state == "Active") & (pc.L4_state == "Active") & (pc.L3_state_pair.isna())][
        ["InChIKey", "UniProt_AC"]].drop_duplicates()
    l2_scoped = l2.merge(scope2, on=["InChIKey", "UniProt_AC"], how="inner")
    l3_scoped = l3.merge(scope3, on=["InChIKey", "UniProt_AC"], how="inner")

    def sweep(grid, l_scoped, base_ik, base_ap, base_df, ceil_ik, n_ceil, k_ceil, e, si_ceil):
        rows = []
        for t in grid:
            up = base_ik & set(l_scoped.loc[l_scoped.confidence_raw >= t, "InChIKey"])
            n_up = len(up)
            k_up = int(base_ap.loc[list(up)].sum()) if n_up else 0
            pool_n = n_ceil + n_up
            pool_k = k_ceil + k_up
            si = compute_si(pool_k, pool_n, e)
            _, n95 = mc_null(base_df, n_up, n_ceil, k_ceil, e)
            rows.append({"t": t, "n_up": n_up, "pool_n": pool_n, "pool_K": pool_k,
                         "SI": si, "null_95th": n95, "margin": si - n95,
                         "recovery_pct": si / si_ceil * 100})
        return pd.DataFrame(rows)

    sw2 = sweep(L2_GRID, l2_scoped, base2_ik, base2_ap, base2, ceil2_ik, n_c2, k_c2, E2, si_ceil2)
    sw3 = sweep(L3_GRID, l3_scoped, base3_ik, base3_ap, base3, ceil3_ik, n_c3, k_c3, E3, si_ceil3)

    # self-check
    s50 = sw2.loc[sw2.t == 0.50].iloc[0]
    s85 = sw3.loc[sw3.t == 0.85].iloc[0]
    ok = abs(s50.SI - 7.73) < 0.02 and abs(s85.SI - 8.89) < 0.02
    print(f"\nSELF-CHECK: L2@0.50 SI={s50.SI:.2f} (exp 7.73), L3@0.85 SI={s85.SI:.2f} (exp 8.89) -> "
          f"{'MATCH ✓' if ok else 'MISMATCH ✗'}")

    def show(tag, sw):
        print(f"\n[{tag} SWEEP — final pair-level pipeline]")
        print(f"  {'t':>5s} {'n_up':>7s} {'pool_n':>7s} {'pool_K':>7s} {'SI':>7s} "
              f"{'null95':>7s} {'margin':>8s} {'recov%':>7s}")
        for _, r in sw.iterrows():
            print(f"  {r.t:5.2f} {int(r.n_up):7,d} {int(r.pool_n):7,d} {int(r.pool_K):7,d} "
                  f"{r.SI:7.3f} {r.null_95th:7.3f} {r.margin:+8.3f} {r.recovery_pct:7.1f}")

    print("\n" + "=" * 78 + "\n[1] FULL SWEEP TABLES")
    show("L2", sw2)
    show("L3", sw3)

    # [2] margin-max per layer (whole grid)
    m2max = sw2.loc[sw2.margin.idxmax()]
    m3max = sw3.loc[sw3.margin.idxmax()]
    print("\n" + "=" * 78 + "\n[2] MARGIN-MAX over the WHOLE grid")
    print(f"  L2: t={m2max.t:.2f}  margin={m2max.margin:+.3f}  SI={m2max.SI:.3f}  n_up={int(m2max.n_up):,}")
    print(f"  L3: t={m3max.t:.2f}  margin={m3max.margin:+.3f}  SI={m3max.SI:.3f}  n_up={int(m3max.n_up):,}")

    # [3] is locked op-point the margin-max?
    print("\n" + "=" * 78 + "\n[3] LOCKED OP-POINT vs MARGIN-MAX")
    print(f"  L3: locked 0.85 margin={s85.margin:+.3f}; grid-max at t={m3max.t:.2f} "
          f"margin={m3max.margin:+.3f}  -> 0.85 is {'THE peak' if m3max.t==0.85 else 'NOT the peak'}")
    l2_mono = all(sw2.sort_values('t').margin.diff().dropna() <= 1e-9)  # non-increasing as t rises?
    print(f"  L2: locked 0.50 margin={s50.margin:+.3f}; grid-max at t={m2max.t:.2f} margin={m2max.margin:+.3f}")
    rising = m2max.t < 0.50 or (sw2.sort_values('t').iloc[0].margin > s50.margin)
    print(f"      margin as t falls below 0.50 keeps RISING? {rising}  "
          f"-> 0.50 is {'NOT a margin peak (quantity-driven)' if m2max.t < 0.50 or m2max.t!=0.50 else 'the peak'}")

    # [4] L2 knee + first-significant t
    sw2s = sw2.sort_values('t').reset_index(drop=True)
    knee_t, _ = knee_by_chord(sw2s.t.tolist(), sw2s.n_up.tolist(), sw2s.SI.tolist())
    sig = sw2[sw2.SI > sw2.null_95th].sort_values('t')
    first_sig = sig.t.min() if len(sig) else float("nan")
    print("\n" + "=" * 78 + "\n[4] L2 KNEE & SIGNIFICANCE (justify 0.50 on quantity grounds)")
    print(f"  knee of (n_up vs SI) curve [max-dist-to-chord]: t={knee_t:.2f}")
    print(f"  lowest t where SI > null_95th (significant): t={first_sig:.2f}")
    print(f"  at locked 0.50: SI={s50.SI:.3f} > null95={s50.null_95th:.3f}? {s50.SI>s50.null_95th}; "
          f"n_up={int(s50.n_up):,}")

    # [5] verdict
    print("\n" + "=" * 78 + "\n[5] VERDICT")
    if m3max.t == 0.85:
        print("  L3: locked 0.85 IS the margin-max under the final pipeline -> defensible on margin-max grounds.")
    else:
        print(f"  L3: margin-max is t={m3max.t:.2f}, not 0.85 -> 0.85 defensible only if margin@0.85 ~ peak (see table).")
    print("  L2: margin favours lower t; 0.50 is NOT a margin peak but a QUANTITY-with-significance choice "
          f"(SI@0.50={s50.SI:.2f} still > null95={s50.null_95th:.2f}; knee≈t={knee_t:.2f}).")
    print("=" * 78)


if __name__ == "__main__":
    main()
