"""
part11_final_op_0p40_0p97.py

FINAL op-points: L2 confidence_raw>=0.40, L3 confidence_raw>=0.97 (margin-max over the
random-assignment null). Recompute everything; reuse verified SI/null/pair-level/web-table
code VERBATIM, only thresholds change.

(1) SI 3-point + Monte-Carlo null   (from part11_si_and_null_final.py)
(2) Synergy-pool expansion          (from part11_reverify_pool_symmetric.py)   [backbone-constrained]
(3) DB web table v2                 (from part11c_assemble_webtable.py)        [NOT backbone-constrained]
    -> WRITES Output/NAR_DBI/Tables/SynerGPCR_WebTable_pairs_v2.csv (v1 kept intact)
(4) Summary block.
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
from typing import Tuple
import csv
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


B = BASE
CAT   = B / "Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv"
CHN   = B / "Output/DB/GPCRactDB/Analysis/Chain_Pattern_Distribution.csv"
PC    = B / "Output/NAR_L3L4/Tables/L3L4_Pair_Catalog_v1.csv"
MAST  = B / "Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv"
L2P   = B / "GPCRactDB/results/denovo/synergpcr_annotation_L2.csv"
L3P   = B / "GPCRactDB/results/denovo/synergpcr_annotation_L3.csv"
LOOK  = B / "Output/NAR_DBI/Tables/protein_graph_type_lookup.csv"
INFP  = B / "Output/NAR_DBI/Tables/PART11C_inference_pairs.csv"
OUTv2 = B / "Output/NAR_DBI/Tables/SynerGPCR_WebTable_pairs_v2.csv"

L2_T, L3_T = 0.40, 0.97
N_REPS, SEED = 10_000, 0
NA_TOK = {"NA", "nan", "None", "", "<NA>"}
# synergy-pool L2-NA scope rule (catalog-consistent, from symmetric script)
MOA_MAP = {"Agonist": "ag", "Partial Agonist": "ag", "PAM": "ag", "Antagonist": "an",
           "Inverse Agonist": "an", "NAM": "an", "Inactive": "inact"}
# web-table per-pair exp-state rule (verbatim from part11c_assemble_webtable.py)
ACTIVE_MOA   = {"Binder", "Agonist", "Antagonist", "Inverse Agonist", "Partial Agonist", "PAM", "NAM"}
INACTIVE_MOA = {"Inactive", "Non-binder"}
LAYER_MAP = {"Layer 1 (Binding)": "L1", "Layer 2 (Proximal)": "L2",
             "Layer 3 (Biased)": "L3", "Layer 4 (Reporter)": "L4"}


def compute_si(k, n, e): return float("nan") if (n == 0 or e == 0) else (k / n) / e
def wilson(k, n): return (0.0, 0.0) if n == 0 else proportion_confint(k, n, 0.05, "wilson")
def si_rate(r, e): return float("nan") if e == 0 else r / e
def p_only(d, L):
    r = d[d.chain_pattern == f"{L}_only"]; return int(r.n_approved.iloc[0]) / int(r.n.iloc[0])


def mc_null(base_df, n_up, n_ceil, k_ceil, e):
    if n_up == 0: return float("nan"), float("nan")
    appr = base_df.is_approved.to_numpy().astype(int)
    rng = np.random.default_rng(SEED)
    N = n_ceil + n_up; idx = np.arange(len(appr)); sis = np.empty(N_REPS)
    for i in range(N_REPS):
        sis[i] = ((k_ceil + int(appr[rng.choice(idx, n_up, replace=False)].sum())) / N) / e
    return float(sis.mean()), float(np.percentile(sis, 95))


def main() -> None:
    cat = pd.read_csv(CAT, low_memory=False).drop_duplicates("InChIKey")
    cat["is_approved"] = (pd.to_numeric(cat.max_phase, errors="coerce") == 4)
    A = lambda c: cat[c] == "Active"
    cati = cat.set_index("InChIKey")
    dist = pd.read_csv(CHN); pc = pd.read_csv(PC, low_memory=False)
    l2 = pd.read_csv(L2P); l3 = pd.read_csv(L3P)
    P = {L: p_only(dist, L) for L in ["L1", "L2", "L3", "L4"]}
    E2 = 1 - (1 - P["L1"]) * (1 - P["L2"]) * (1 - P["L4"])
    E3 = 1 - (1 - P["L1"]) * (1 - P["L3"]) * (1 - P["L4"])
    EF = 1 - (1 - P["L1"]) * (1 - P["L2"]) * (1 - P["L3"]) * (1 - P["L4"])

    def stat(ik):
        ik = list(ik)
        return (len(ik), int(cati.loc[ik].is_approved.sum())) if ik else (0, 0)

    ceil2 = set(cat.loc[A("L1_state") & A("L2_state") & A("L4_state"), "InChIKey"])
    ceil3 = set(cat.loc[A("L1_state") & A("L3_state") & A("L4_state"), "InChIKey"])
    full  = set(cat.loc[A("L1_state") & A("L2_state") & A("L3_state") & A("L4_state"), "InChIKey"])
    base2 = cat[A("L1_state") & A("L4_state") & cat.L2_state.isna()].drop_duplicates("InChIKey")
    base3 = cat[A("L1_state") & A("L4_state") & cat.L3_state.isna()].drop_duplicates("InChIKey")
    base2_ik, base3_ik = set(base2.InChIKey), set(base3.InChIKey)
    n_c2, k_c2 = stat(ceil2); n_c3, k_c3 = stat(ceil3)

    # pair scopes (synergy pool)
    m = pd.read_csv(MAST, usecols=["Ligand_InChIKey", "GPCR_UniProt", "Assay_Layer", "MoA"]).rename(
        columns={"Ligand_InChIKey": "InChIKey", "GPCR_UniProt": "UniProt_AC"})
    m2 = m[m.Assay_Layer == "Layer 2 (Proximal)"].copy(); m2["x"] = m2.MoA.map(MOA_MAP)
    l2_known = set(map(tuple, m2.dropna(subset=["x"]).drop_duplicates(["InChIKey", "UniProt_AC"])[
        ["InChIKey", "UniProt_AC"]].values))
    cand = pc[(pc.L1_state == "Active") & (pc.L4_state == "Active")][["InChIKey", "UniProt_AC"]].drop_duplicates()
    ck = list(map(tuple, cand.values))
    scope2 = cand[[k not in l2_known for k in ck]]
    scope3 = pc[(pc.L1_state == "Active") & (pc.L4_state == "Active") & (pc.L3_state_pair.isna())][
        ["InChIKey", "UniProt_AC"]].drop_duplicates()
    up_l2 = base2_ik & set(l2.merge(scope2, on=["InChIKey", "UniProt_AC"]).query("confidence_raw>=@L2_T").InChIKey)
    up_l3 = base3_ik & set(l3.merge(scope3, on=["InChIKey", "UniProt_AC"]).query("confidence_raw>=@L3_T").InChIKey)

    # ── (1) SI 3-point + null ─────────────────────────────────────────────────
    def panel(tag, base_df, base_ik, up_ik, ceil_ik, n_ceil, k_ceil, e, ceil_si):
        def bar(name, ik):
            n, k = stat(ik); lo, hi = wilson(k, n)
            return dict(name=name, n=n, k=k, rate=(k / n if n else 0), si=compute_si(k, n, e),
                        lo=si_rate(lo, e), hi=si_rate(hi, e))
        bb, bg, bc = bar("baseline", base_ik), bar("GPCRact", ceil_ik | up_ik), bar("ceiling", ceil_ik)
        nm, n95 = mc_null(base_df, len(up_ik), n_ceil, k_ceil, e)
        print(f"\n[{tag} PANEL @ t={'0.40' if tag=='L2' else '0.97'}]  E_Bliss={e:.6f}")
        for b in (bb, bg, bc):
            print(f"   {b['name']:9s} n={b['n']:6,d} K={b['k']:4d} rate={b['rate']*100:6.3f}% "
                  f"SI={b['si']:7.3f} CI[{b['lo']:6.3f},{b['hi']:6.3f}]")
        print(f"   null_mean={nm:.3f} null_95th={n95:.3f}  margin={bg['si']-n95:+.3f}  "
              f"SI>null95? {bg['si']>n95}  recovery={bg['si']/ceil_si*100:.1f}%")
        return bg["si"], n95, bg["si"] / ceil_si * 100
    si_ceil2 = compute_si(k_c2, n_c2, E2); si_ceil3 = compute_si(k_c3, n_c3, E3)
    print("=" * 88 + "\n(1) SI 3-POINT + MONTE-CARLO NULL  @ L2>=0.40 / L3>=0.97")
    si_g2, n95_2, rec2 = panel("L2", base2, base2_ik, up_l2, ceil2, n_c2, k_c2, E2, si_ceil2)
    si_g3, n95_3, rec3 = panel("L3", base3, base3_ik, up_l3, ceil3, n_c3, k_c3, E3, si_ceil3)
    print(f"\n  self-check vs sweep: L2@0.40 SI≈7.67/null≈7.17 -> got {si_g2:.2f}/{n95_2:.2f}; "
          f"L3@0.97 SI≈10.35/null≈7.36 -> got {si_g3:.2f}/{n95_3:.2f}")

    # ── (2) synergy-pool expansion (backbone-constrained) ─────────────────────
    base = ceil2; only2 = (up_l2 - up_l3) - base; only3 = (up_l3 - up_l2) - base; both = (up_l2 & up_l3) - base
    total = base | up_l2 | up_l3
    n_tot, k_tot = stat(total)
    print("\n" + "=" * 88 + "\n(2) SYNERGY-POOL EXPANSION (pair-level, backbone-constrained)  @ 0.40/0.97")
    for nm_, s in [("base(exp L1+L2+L4)", base), ("GPCRact_L2 only", only2),
                   ("GPCRact_L3 only", only3), ("GPCRact_L2&L3", both), ("TOTAL", total)]:
        n, k = stat(s); print(f"   {nm_:20s} n={n:6,d}  approved={k}")
    sb, s2, s3, sbo = stat(base)[0], stat(only2)[0], stat(only3)[0], stat(both)[0]
    assert sb + s2 + s3 + sbo == n_tot, "partition mismatch"
    print(f"   set-algebra: {sb}+{s2}+{s3}+{sbo} = {sb+s2+s3+sbo} == TOTAL {n_tot}? True")
    print(f"   AI-enabled (TOTAL-base) = {n_tot - n_c2:,};  expansion 2,580 -> {n_tot:,} (x{n_tot/n_c2:.2f})")
    known = np.zeros(len(cat), bool)
    for c in ["L1_state", "L2_state", "L3_state", "L4_state"]:
        known |= ~cat[c].astype(str).isin(NA_TOK)
    db = set(cat.loc[known, "InChIKey"])
    print(f"   carve-out 1 no-op: pool members outside DB-scale = {len(total - db)} (expect 0)")

    # ── (3) DB web table v2 (NOT backbone-constrained) ────────────────────────
    print("\n" + "=" * 88 + "\n(3) DB WEB TABLE v2 (cell-level)  @ L2 AI>=0.40 / L3 AI>=0.97")
    # carve-out 1
    l2w = l2[l2.InChIKey.isin(db)].copy(); l3w = l3[l3.InChIKey.isin(db)].copy()
    # per-pair exp states (web rule)
    mw = m.copy(); mw = mw[mw.Assay_Layer.isin(LAYER_MAP)]
    mw["layer"] = mw.Assay_Layer.map(LAYER_MAP)
    mw["st"] = np.where(mw.MoA.isin(ACTIVE_MOA), "Active", np.where(mw.MoA.isin(INACTIVE_MOA), "Inactive", "Unknown"))
    def agg(s): return "Active" if (s == "Active").any() else ("Inactive" if (s == "Inactive").any() else "NA")
    grp = mw.groupby(["InChIKey", "UniProt_AC", "layer"]).st.agg(agg).unstack("layer")
    for c in ["L1", "L2", "L3", "L4"]:
        if c not in grp: grp[c] = "NA"
    grp = grp[["L1", "L2", "L3", "L4"]].fillna("NA")
    grp.columns = ["L1_state", "L2_state", "L3_state", "L4_state"]; grp = grp.reset_index()
    bw = l2w[["InChIKey", "UniProt_AC", "ago_prob", "ant_prob", "confidence_raw"]].rename(
        columns={"ago_prob": "ag2", "ant_prob": "an2", "confidence_raw": "c2"}).merge(
        l3w[["InChIKey", "UniProt_AC", "ago_prob", "ant_prob", "confidence_raw"]].rename(
            columns={"ago_prob": "ag3", "ant_prob": "an3", "confidence_raw": "c3"}),
        on=["InChIKey", "UniProt_AC"], how="outer")
    bw = bw.merge(grp, on=["InChIKey", "UniProt_AC"], how="left")
    for c in ["L1_state", "L2_state", "L3_state", "L4_state"]:
        bw[c] = bw[c].fillna("NA")
    pgt = pd.read_csv(LOOK).set_index("UniProt_AC").protein_graph_type.to_dict()
    bw["protein_graph_type"] = bw.UniProt_AC.map(pgt)
    appr_map = cati.is_approved.to_dict()
    bw["is_approved"] = bw.InChIKey.map(appr_map).fillna(False)
    exp = {"Active", "Inactive"}
    l2e = ~bw.L2_state.isin(exp); l3e = ~bw.L3_state.isin(exp)
    bw["L2_ai_active"] = np.where(l2e, bw.c2 >= L2_T, pd.NA)
    bw["L2_confidence"] = np.where(l2e, bw.c2.round(4), pd.NA)
    bw["L2_predicted_moa"] = np.where(l2e & (bw.c2 >= L2_T), np.where(bw.ag2 >= bw.an2, "agonist", "antagonist"), pd.NA)
    bw["L3_ai_active"] = np.where(l3e, bw.c3 >= L3_T, pd.NA)
    bw["L3_confidence"] = np.where(l3e, bw.c3.round(4), pd.NA)
    bw["L3_predicted_moa"] = np.where(l3e & (bw.c3 >= L3_T), np.where(bw.ag3 >= bw.an3, "agonist", "antagonist"), pd.NA)
    bw["l1_evidence"] = np.where(bw.L1_state.isin(exp), "experimental", "predicted")
    keep = bw[l2e | l3e].copy()
    cols = ["InChIKey", "UniProt_AC", "is_approved", "protein_graph_type", "L1_state", "l1_evidence",
            "L2_state", "L2_ai_active", "L2_predicted_moa", "L2_confidence",
            "L3_state", "L3_ai_active", "L3_predicted_moa", "L3_confidence", "L4_state"]
    web = keep[cols]
    # guardrails
    v_l2 = int((web.L2_state.isin(exp) & web.L2_ai_active.notna()).sum())
    v_l3 = int((web.L3_state.isin(exp) & web.L3_ai_active.notna()).sum())
    has_l4ai = any(("L4" in c) and any(k in c for k in ("ai", "predicted", "confidence")) for c in web.columns)
    has_si = any(("SI" in c) or ("l1_ai" in c.lower()) for c in web.columns)
    assert v_l2 == 0 and v_l3 == 0 and not has_l4ai and not has_si, "guardrail violation"
    print(f"   guardrails: exp-overwrite L2={v_l2}/L3={v_l3}, L4-AI={has_l4ai}, SI/L1-AI col={has_si} (all clean)")
    a2 = web.L2_ai_active == True; a3 = web.L3_ai_active == True
    ub_comp = set(web.loc[a2 | a3, "InChIKey"]); ub_appr = web.loc[(a2 | a3) & (web.is_approved == True), "InChIKey"].nunique()
    print(f"   rows={len(web):,}  L2 AI-Active cells={int(a2.sum()):,}  L3 AI-Active cells={int(a3.sum()):,}")
    print(f"   unique compounds with >=1 AI-Active cell={len(ub_comp):,}  approved among them={ub_appr:,}")
    print(f"   (v1 @0.50/0.85 was L2=218,176 / L3=179,959 cells, 174,952 compounds, 844 approved)")

    # write v2
    header = (
        "# SynerGPCR_WebTable_pairs_v2.csv  (STEP 11C step 3, FINAL op-points L2>=0.40 / L3>=0.97)\n"
        "# v1 (0.50/0.85) kept intact. One row per (InChIKey,UniProt_AC) with >=1 emittable AI cell.\n"
        "# AI-Active: L2 confidence>=0.40, L3 confidence>=0.97 (confidence=binding_prob*max(ago,ant)).\n"
        "# Predicted L1 never a synergy field / never in SI; L4 never predicted. AI fields blank where exp exists.\n"
        "# COLOUR: green=exp Active, red=exp Inactive, grey=exp NA, blue=AI-Active (NA cells only) — frontend-applied.\n"
    )
    with open(OUTv2, "w") as fh:
        fh.write(header); web.to_csv(fh, index=False)
    print(f"   WROTE {OUTv2}  ({len(web):,} rows)")

    # repurposing decomposition (GPCR-data-bearing denominator)
    appr_all = set(cat.loc[cat.is_approved, "InChIKey"])
    inf_comp = set(pd.read_csv(INFP, usecols=["InChIKey"]).InChIKey)
    rem = appr_all & db
    B3 = rem - inf_comp; rem2 = rem & inf_comp
    web_emit_appr = set(web.loc[web.is_approved == True, "InChIKey"])  # emittable (threshold-indep)
    B1_emit = rem2 & web_emit_appr
    B2 = rem2 - web_emit_appr; B4 = appr_all - db
    # active repurposing = approved with >=1 AI-Active cell at NEW op-points
    rep_active = ub_appr
    print(f"\n   repurposing decomposition (2,420 approved):")
    print(f"     in DB-scale={len(rem):,}; emittable-in-web B1={len(B1_emit):,}; carve2 B2={len(B2):,}; "
          f"no-predictable-pair B3={len(B3):,}; outside-DB B4={len(B4):,}  sum={len(B1_emit)+len(B2)+len(B3)+len(B4):,}")
    print(f"     APPROVED with >=1 AI-Active cell @0.40/0.97 (repurposing N) = {rep_active:,}")

    # full chain (reference only)
    n_f, k_f = stat(full); si_full = compute_si(k_f, n_f, EF)

    # ── (4) SUMMARY ───────────────────────────────────────────────────────────
    def si_of(ik, e):
        n, k = stat(ik); return compute_si(k, n, e)
    si_base2, si_base3 = si_of(base2_ik, E2), si_of(base3_ik, E3)
    print("\n" + "=" * 88 + "\n(4) SUMMARY @ FINAL OP-POINTS 0.40/0.97  (tag valid / reference-only / unchanged)")
    print(f"   DB-scale universe                 = {len(db):,}            (unchanged)")
    print(f"   L2 panel: base {si_base2:.2f} -> GPCRact@0.40 {si_g2:.2f} -> ceiling {si_ceil2:.2f}; null95={n95_2:.2f}; recovery={rec2:.1f}% (valid)")
    print(f"   L3 panel: base {si_base3:.2f} -> GPCRact@0.97 {si_g3:.2f} -> ceiling {si_ceil3:.2f}; null95={n95_3:.2f}; recovery={rec3:.1f}% (valid)")
    print(f"   synergy pool 2,580 -> {n_tot:,} (approved {k_c2}->{k_tot}); AI-enabled={n_tot-n_c2:,} (valid)")
    print(f"   full-chain SI                     = {si_full:.2f}            (REFERENCE-ONLY; L2-L3 dependent)")
    print(f"   DB-loaded AI-Active cells         = L2 {int(a2.sum()):,} / L3 {int(a3.sum()):,}; compounds {len(ub_comp):,} (valid, cell-level)")
    print(f"   repurposing N (approved, AI-Active) = {rep_active:,} of 2,420 (valid)")
    print(f"   BACC (clinical, op-point-independent) = 0.875 / 0.834   (unchanged)")
    print(f"   FN Table B (approved)             = L2_FN/L1_FN 112/7 = 16x  (unchanged)")
    print("=" * 88)


if __name__ == "__main__":
    main()
