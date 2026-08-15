"""
part11_si_and_null_final.py  — READ-ONLY. Writes nothing but this script's stdout.

Bliss SI + Monte-Carlo random-assignment null at the FINAL locked op-points
(L2 conf>=0.50, L3 conf>=0.85), pair-level matching for both AI layers
(mirroring part11_reverify_pool_symmetric.py).

SI / E_Bliss / Wilson code reused VERBATIM from the verified
simulate_SI_expansion_L3.py:
  compute_si(k,n,e)=(k/n)/e ; wilson_ci via proportion_confint(method='wilson') ;
  si_from_rate(rate,e)=rate/e ; get_p_only(dist,layer)=K/N of '<layer>_only' ;
  E_Bliss(chain)=1-prod(1-P_Lk_only).

Null (locked method): from the L1+L4 baseline pool (layer experimentally NA), draw the SAME
number of compounds the AI added, 10,000 reps, recompute combined SI, take mean & 95th pct.

L2-L3 dependency: the full-chain independence SI is reported FOR REFERENCE ONLY (Bliss
independence violated, Fisher OR=9.24, p~0). No dependency-corrected SI is computed.
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
L2_T, L3_T = 0.50, 0.85
N_REPS, SEED = 10_000, 0
# catalog-consistent MoA->state (build_L3L4_pair_catalog MOA_MAP; Unknown/Binder dropped)
MOA_MAP = {"Agonist": "agonist", "Partial Agonist": "agonist", "PAM": "agonist",
           "Antagonist": "antagonist", "Inverse Agonist": "antagonist", "NAM": "antagonist",
           "Inactive": "inactive"}


# ── verified helpers (verbatim) ───────────────────────────────────────────────
def compute_si(k, n, e):
    return float("nan") if (n == 0 or e == 0) else (k / n) / e


def wilson_ci(k, n, alpha=0.05) -> Tuple[float, float]:
    return (0.0, 0.0) if n == 0 else proportion_confint(k, n, alpha=alpha, method="wilson")


def si_from_rate(rate, e):
    return float("nan") if e == 0 else rate / e


def get_p_only(dist, layer):
    row = dist[dist["chain_pattern"] == f"{layer}_only"]
    n, k = int(row["n"].iloc[0]), int(row["n_approved"].iloc[0])
    return k / n


def main() -> None:
    print("=" * 84)
    print("BLISS SI + MONTE-CARLO NULL  (final op-points L2>=0.50 / L3>=0.85, pair-level)")
    print("=" * 84)

    cat = pd.read_csv(CAT, low_memory=False)
    cat["is_approved"] = (pd.to_numeric(cat["max_phase"], errors="coerce") == 4)
    cat = cat.drop_duplicates("InChIKey")
    A = lambda c: cat[c] == "Active"
    dist = pd.read_csv(CHAIN_D)
    pc = pd.read_csv(PAIR_CAT, low_memory=False)
    l2 = pd.read_csv(L2_PRED)
    l3 = pd.read_csv(L3_PRED)

    P_L1, P_L2, P_L3, P_L4 = (get_p_only(dist, x) for x in ["L1", "L2", "L3", "L4"])
    E_L1L2L4 = 1 - (1 - P_L1) * (1 - P_L2) * (1 - P_L4)
    E_L1L3L4 = 1 - (1 - P_L1) * (1 - P_L3) * (1 - P_L4)
    E_FULL   = 1 - (1 - P_L1) * (1 - P_L2) * (1 - P_L3) * (1 - P_L4)
    print(f"P_only: L1={P_L1:.6f} L2={P_L2:.6f} L3={P_L3:.6f} L4={P_L4:.6f}")
    print(f"E_Bliss: L1+L2+L4={E_L1L2L4:.6f}  L1+L3+L4={E_L1L3L4:.6f}  full={E_FULL:.6f}")

    # ── pools (compound-level) ────────────────────────────────────────────────
    def pool(mask):
        s = cat[mask]
        return set(s.InChIKey), len(s), int(s.is_approved.sum())
    ceil_l2_ik, n_c2, k_c2 = pool(A("L1_state") & A("L2_state") & A("L4_state"))
    ceil_l3_ik, n_c3, k_c3 = pool(A("L1_state") & A("L3_state") & A("L4_state"))
    full_ik,    n_f,  k_f  = pool(A("L1_state") & A("L2_state") & A("L3_state") & A("L4_state"))
    base_l2 = cat[A("L1_state") & A("L4_state") & cat.L2_state.isna()].drop_duplicates("InChIKey")
    base_l3 = cat[A("L1_state") & A("L4_state") & cat.L3_state.isna()].drop_duplicates("InChIKey")
    base_l2_ik, n_b2, k_b2 = set(base_l2.InChIKey), len(base_l2), int(base_l2.is_approved.sum())
    base_l3_ik, n_b3, k_b3 = set(base_l3.InChIKey), len(base_l3), int(base_l3.is_approved.sum())

    # ── pair-level AI upgrades (symmetric to L3; identical to verified script) ─
    # L3 scope = catalog L1&L4 Active & L3_state_pair NA
    scope_l3 = pc[(pc.L1_state == "Active") & (pc.L4_state == "Active") & (pc.L3_state_pair.isna())][
        ["InChIKey", "UniProt_AC"]].drop_duplicates()
    hc_l3 = set(l3.merge(scope_l3, on=["InChIKey", "UniProt_AC"], how="inner")
                .query("confidence_raw >= @L3_T")["InChIKey"])
    up_l3_ik = base_l3_ik & hc_l3

    # L2 scope = catalog L1&L4 Active pairs minus pairs with experimental L2 (Active/Inactive)
    m2 = pd.read_csv(MASTER, usecols=["Ligand_InChIKey", "GPCR_UniProt", "Assay_Layer", "MoA"]).rename(
        columns={"Ligand_InChIKey": "InChIKey", "GPCR_UniProt": "UniProt_AC"})
    m2 = m2[m2.Assay_Layer == "Layer 2 (Proximal)"].copy()
    m2["lab"] = m2.MoA.map(MOA_MAP)
    l2_known = set(map(tuple, m2.dropna(subset=["lab"]).drop_duplicates(["InChIKey", "UniProt_AC"])[
        ["InChIKey", "UniProt_AC"]].values))
    cand = pc[(pc.L1_state == "Active") & (pc.L4_state == "Active")][
        ["InChIKey", "UniProt_AC"]].drop_duplicates()
    ck = list(map(tuple, cand.values))
    scope_l2 = cand[[k not in l2_known for k in ck]]
    hc_l2 = set(l2.merge(scope_l2, on=["InChIKey", "UniProt_AC"], how="inner")
                .query("confidence_raw >= @L2_T")["InChIKey"])
    up_l2_ik = base_l2_ik & hc_l2

    k_up2 = int(base_l2.set_index("InChIKey").loc[list(up_l2_ik), "is_approved"].sum())
    k_up3 = int(base_l3.set_index("InChIKey").loc[list(up_l3_ik), "is_approved"].sum())
    print(f"\nAI pair-level upgrades: up_l2={len(up_l2_ik):,} (approved {k_up2}); "
          f"up_l3={len(up_l3_ik):,} (approved {k_up3})")

    # ── 3-point panels ────────────────────────────────────────────────────────
    def bar(name, ik, e):
        n = len(ik)
        k = int(cat.set_index("InChIKey").loc[list(ik), "is_approved"].sum()) if n else 0
        rate = k / n if n else float("nan")
        lo, hi = wilson_ci(k, n)
        return {"name": name, "n": n, "k": k, "rate": rate, "e": e,
                "SI": compute_si(k, n, e), "ci_lo": si_from_rate(lo, e), "ci_hi": si_from_rate(hi, e)}

    def panel(tag, base_ik, up_ik, ceil_ik, e):
        rows = [bar(f"{tag} baseline (L1+L4, {tag[-2:]}=NA)", base_ik, e),
                bar(f"{tag} GPCRact (ceiling∪AI)", ceil_ik | up_ik, e),
                bar(f"{tag} exp ceiling", ceil_ik, e)]
        print(f"\n[{tag} PANEL]  (E_Bliss={e:.6f})")
        print(f"  {'bar':32s} {'n':>7s} {'K':>5s} {'rate%':>7s} {'SI':>8s} {'95%CI':>17s}")
        for r in rows:
            print(f"  {r['name']:32s} {r['n']:7,d} {r['k']:5d} {r['rate']*100:7.3f} "
                  f"{r['SI']:8.4f}  [{r['ci_lo']:6.3f},{r['ci_hi']:6.3f}]")
        return rows

    print("\n" + "=" * 84 + "\n[1] 3-POINT SI COMPARISON (pair-level, final op-points)")
    p2 = panel("L2", base_l2_ik, up_l2_ik, ceil_l2_ik, E_L1L2L4)
    p3 = panel("L3", base_l3_ik, up_l3_ik, ceil_l3_ik, E_L1L3L4)
    si_gp2, si_ceil2 = p2[1]["SI"], p2[2]["SI"]
    si_gp3, si_ceil3 = p3[1]["SI"], p3[2]["SI"]

    # ── [2] Monte-Carlo random-assignment null ────────────────────────────────
    def mc_null(base_df, n_up, n_ceil, k_ceil, e):
        appr = base_df.is_approved.to_numpy().astype(int)
        rng = np.random.default_rng(SEED)
        N = n_ceil + n_up
        sis = np.empty(N_REPS)
        idx = np.arange(len(appr))
        for i in range(N_REPS):
            pick = rng.choice(idx, size=n_up, replace=False)
            K = k_ceil + int(appr[pick].sum())
            sis[i] = (K / N) / e
        return float(sis.mean()), float(np.percentile(sis, 95))

    print("\n" + "=" * 84 + "\n[2] MONTE-CARLO RANDOM-ASSIGNMENT NULL  "
          f"({N_REPS:,} reps, draw n_added from L1+L4 baseline)")
    nm2, n95_2 = mc_null(base_l2, len(up_l2_ik), n_c2, k_c2, E_L1L2L4)
    nm3, n95_3 = mc_null(base_l3, len(up_l3_ik), n_c3, k_c3, E_L1L3L4)
    for tag, si, nm, n95 in [("L2", si_gp2, nm2, n95_2), ("L3", si_gp3, nm3, n95_3)]:
        margin = si - n95
        print(f"  {tag}: GPCRact_SI={si:.4f}  null_mean={nm:.4f}  null_95th={n95:.4f}  "
              f"margin={margin:+.4f}  SI>null_95th? {si > n95}")

    # ── [3] Recovery % (single-AI-layer; full chain not valid) ────────────────
    print("\n" + "=" * 84 + "\n[3] RECOVERY % (SI-based, per single-AI-layer chain)")
    rec2 = si_gp2 / si_ceil2 * 100
    rec3 = si_gp3 / si_ceil3 * 100
    print(f"  recovery_L2 = GPCRact_L2_SI / L2_ceiling_SI = {si_gp2:.4f} / {si_ceil2:.4f} = {rec2:.1f}%")
    print(f"  recovery_L3 = GPCRact_L3_SI / L3_ceiling_SI = {si_gp3:.4f} / {si_ceil3:.4f} = {rec3:.1f}%")
    print("  Reported per single-AI-layer chain (L1+L2+L4, L1+L3+L4) because the full-chain SI")
    print("  is not well-defined under L2-L3 dependence (Bliss independence violated).")

    # ── [4] Full-chain independence SI (REFERENCE ONLY) ───────────────────────
    si_full = compute_si(k_f, n_f, E_FULL)
    lo_f, hi_f = wilson_ci(k_f, n_f)
    print("\n" + "=" * 84 + "\n[4] FULL-CHAIN (L1+L2+L3+L4) INDEPENDENCE SI — REFERENCE ONLY")
    print(f"  n={n_f}, K={k_f}, rate={k_f/n_f*100:.3f}%, E_Bliss_full={E_FULL:.6f}")
    print(f"  SI = {si_full:.4f}  95%CI [{si_from_rate(lo_f,E_FULL):.3f},{si_from_rate(hi_f,E_FULL):.3f}]")
    print("  *** LIMITATION: independence assumption VIOLATED (L2-L3 dependent; Fisher OR=9.24,")
    print("      p~0). Reported for reference only — NOT a valid synergy estimate. No corrected SI. ***")

    # ── [5] Summary block for cover letter / abstract ─────────────────────────
    print("\n" + "=" * 84 + "\n[5] SUMMARY — numbers for cover letter / abstract  (tag: valid / reference-only)")
    print(f"  base L1+L2+L4 SI         = {p2[0]['SI']:.2f}   (valid)")
    print(f"  GPCRact L2@0.50 SI       = {si_gp2:.2f}   (valid)")
    print(f"  exp ceiling L1+L2+L4 SI  = {si_ceil2:.2f}   (valid)")
    print(f"  L2 null_95th             = {n95_2:.2f}   (valid; GPCRact {'>' if si_gp2>n95_2 else '<='} null)")
    print(f"  recovery_L2              = {rec2:.1f}% (valid)")
    print(f"  base L1+L4(L3=NA) SI     = {p3[0]['SI']:.2f}   (valid)")
    print(f"  GPCRact L3@0.85 SI       = {si_gp3:.2f}   (valid)")
    print(f"  exp ceiling L1+L3+L4 SI  = {si_ceil3:.2f}   (valid)")
    print(f"  L3 null_95th             = {n95_3:.2f}   (valid; GPCRact {'>' if si_gp3>n95_3 else '<='} null)")
    print(f"  recovery_L3              = {rec3:.1f}% (valid)")
    print(f"  full-chain SI            = {si_full:.2f}   (REFERENCE-ONLY; L2-L3 dependent)")
    print("=" * 84)


if __name__ == "__main__":
    main()
