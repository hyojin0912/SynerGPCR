"""
webexport_step6_derived.py — WebExport_v2 STEP 6 (A). READ-ONLY on sources.

Two derived web tables from the SAME verified synergy-pool logic / export files (no new pool,
no SI re-derivation). Builds in memory, self-checks, writes ONLY if all pass.

  1. candidate_shortlist.csv  — the 1,357 synergy AI-enabled compounds (only via AI)
  2. repurposing_view.csv     — the 842 approved drugs gaining >=1 AI-Active prediction
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import sys
import numpy as np
import pandas as pd

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


B = BASE
OUT = B / "Output/NAR_DBI/WebExport_v2"
CAT  = B / "Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv"
PC   = B / "Output/NAR_L3L4/Tables/L3L4_Pair_Catalog_v1.csv"
MAST = B / "Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv"
L2P  = B / "GPCRactDB/results/denovo/synergpcr_annotation_L2.csv"
L3P  = B / "GPCRactDB/results/denovo/synergpcr_annotation_L3.csv"
CL   = OUT / "compound_lookup.csv"
TL   = OUT / "target_lookup.csv"
GP   = OUT / "gpcract_predictions.csv"

L2_T, L3_T = 0.40, 0.97
MOA_MAP = {"Agonist": "ag", "Partial Agonist": "ag", "PAM": "ag", "Antagonist": "an",
           "Inverse Agonist": "an", "NAM": "an", "Inactive": "inact"}


def load(p, **kw):
    if not p.exists(): sys.exit(f"STOP: missing {p}")
    df = pd.read_csv(p, **kw); print(f"  loaded {p.name}: {len(df):,} rows"); return df


print("=" * 78); print("LOADING"); print("=" * 78)
cat = load(CAT, low_memory=False).drop_duplicates("InChIKey")
cat["is_approved"] = (pd.to_numeric(cat.max_phase, errors="coerce") == 4)
A = lambda c: cat[c] == "Active"
pc = load(PC, low_memory=False)
l2 = load(L2P); l3 = load(L3P)
cl = load(CL, low_memory=False)
tl = load(TL, low_memory=False)
gp = load(GP, low_memory=False)

name_map = cl.set_index("InChIKey")["preferred_name"].to_dict()
appr_map = cat.set_index("InChIKey")["is_approved"].to_dict()
gene_map = tl.set_index("UniProt_AC")["gene_name"].to_dict()

# ── verified synergy-pool scopes (identical to part11_final_op_0p40_0p97) ──
ceil2 = set(cat.loc[A("L1_state") & A("L2_state") & A("L4_state"), "InChIKey"])     # base
base2_ik = set(cat[A("L1_state") & A("L4_state") & cat.L2_state.isna()].InChIKey)
base3_ik = set(cat[A("L1_state") & A("L4_state") & cat.L3_state.isna()].InChIKey)

m = load(MAST, low_memory=False, usecols=["Ligand_InChIKey", "GPCR_UniProt", "Assay_Layer", "MoA"]).rename(
    columns={"Ligand_InChIKey": "InChIKey", "GPCR_UniProt": "UniProt_AC"})
m2 = m[m.Assay_Layer == "Layer 2 (Proximal)"].copy(); m2["x"] = m2.MoA.map(MOA_MAP)
l2_known = set(map(tuple, m2.dropna(subset=["x"]).drop_duplicates(["InChIKey", "UniProt_AC"])[
    ["InChIKey", "UniProt_AC"]].values))
cand = pc[(pc.L1_state == "Active") & (pc.L4_state == "Active")][["InChIKey", "UniProt_AC"]].drop_duplicates()
ckk = list(map(tuple, cand.values))
scope2 = cand[[k not in l2_known for k in ckk]]
scope3 = pc[(pc.L1_state == "Active") & (pc.L4_state == "Active") & (pc.L3_state_pair.isna())][
    ["InChIKey", "UniProt_AC"]].drop_duplicates()

l2q = l2.merge(scope2, on=["InChIKey", "UniProt_AC"]).query("confidence_raw>=@L2_T")
l3q = l3.merge(scope3, on=["InChIKey", "UniProt_AC"]).query("confidence_raw>=@L3_T")
up_l2 = base2_ik & set(l2q.InChIKey)
up_l3 = base3_ik & set(l3q.InChIKey)

# best (max-confidence) qualifying pair per compound -> confidence + MoA (ago>=ant rule)
def best(qdf):
    q = qdf.sort_values("confidence_raw", ascending=False).drop_duplicates("InChIKey").copy()
    q["moa"] = np.where(q.ago_prob >= q.ant_prob, "agonist", "antagonist")
    return q.set_index("InChIKey")[["confidence_raw", "moa"]]
l2best, l3best = best(l2q), best(l3q)

# ── AI-enabled decomposition (mutually exclusive; == STEP2 1157/17/183) ──
base = ceil2
only2 = (up_l2 - up_l3) - base
only3 = (up_l3 - up_l2) - base
both  = (up_l2 & up_l3) - base
ai_enabled = only2 | only3 | both
print(f"\nAI-enabled: only2={len(only2)}  only3={len(only3)}  both={len(both)}  TOTAL={len(ai_enabled)}")

# ==========================================================================
# FILE 1 — candidate_shortlist.csv
# ==========================================================================
rows = []
for ik in ai_enabled:
    if ik in both:
        layers, chain = "L2;L3", "L1 + L2(AI) + L3(AI) + L4  (Full chain via AI)"
    elif ik in only2:
        layers, chain = "L2", "L1 + L2(AI) + L4"
    else:
        layers, chain = "L3", "L1 + L3(AI) + L4"
    l2c = l2best.loc[ik] if ik in l2best.index and layers != "L3" else None
    l3c = l3best.loc[ik] if ik in l3best.index and layers != "L2" else None
    rows.append(dict(
        InChIKey=ik, preferred_name=name_map.get(ik),
        ai_completed_chain=chain, ai_filled_layers=layers,
        L2_confidence=(round(float(l2c.confidence_raw), 4) if l2c is not None else pd.NA),
        L2_predicted_moa=(l2c.moa if l2c is not None else pd.NA),
        L3_confidence=(round(float(l3c.confidence_raw), 4) if l3c is not None else pd.NA),
        L3_predicted_moa=(l3c.moa if l3c is not None else pd.NA),
        is_approved=bool(appr_map.get(ik, False))))
cand_sl = pd.DataFrame(rows)[["InChIKey", "preferred_name", "ai_completed_chain",
    "ai_filled_layers", "L2_confidence", "L2_predicted_moa", "L3_confidence",
    "L3_predicted_moa", "is_approved"]]
print(f"\ncandidate_shortlist rows={len(cand_sl):,}  approved={int(cand_sl.is_approved.sum())}")
print(cand_sl.ai_filled_layers.value_counts().to_dict())

# ==========================================================================
# FILE 2 — repurposing_view.csv
# ==========================================================================
def ai_true(s): return (s == True) | (s.astype(str) == "True")
appr_set = set(cl.loc[cl.is_approved == True, "InChIKey"])
rv = gp[gp.InChIKey.isin(appr_set) & (ai_true(gp.L2_ai_active) | ai_true(gp.L3_ai_active))].copy()
rv["preferred_name"] = rv.InChIKey.map(name_map)
rv["gene_name"] = rv.UniProt_AC.map(gene_map)
repur = rv[["InChIKey", "preferred_name", "UniProt_AC", "gene_name",
            "L2_ai_active", "L2_predicted_moa", "L2_confidence",
            "L3_ai_active", "L3_predicted_moa", "L3_confidence"]]
print(f"\nrepurposing_view rows={len(repur):,}  unique approved drugs={repur.InChIKey.nunique():,}")

# ==========================================================================
# SELF-CHECKS
# ==========================================================================
print("\n" + "=" * 78); print("SELF-CHECKS"); print("=" * 78)
checks = []
def ck(n, c, g, e): checks.append((n, bool(c), g, e))
ck("candidate_shortlist == 1,357", len(cand_sl) == 1357, len(cand_sl), 1357)
ck("decomposition only2/only3/both == 1157/17/183",
   (len(only2), len(only3), len(both)) == (1157, 17, 183), (len(only2), len(only3), len(both)), (1157, 17, 183))
ck("candidate approved == 36 (201-165)", int(cand_sl.is_approved.sum()) == 36, int(cand_sl.is_approved.sum()), 36)
ck("repurposing_view unique approved == 842", repur.InChIKey.nunique() == 842, repur.InChIKey.nunique(), 842)
ck("repurposing all approved", set(repur.InChIKey).issubset(appr_set), len(set(repur.InChIKey) - appr_set), 0)
ck("candidate InChIKeys ⊆ compound_lookup", set(cand_sl.InChIKey).issubset(set(cl.InChIKey)),
   len(set(cand_sl.InChIKey) - set(cl.InChIKey)), 0)
allpass = all(c[1] for c in checks)
for n, ok, g, e in checks: print(f"  [{'PASS' if ok else 'FAIL'}] {n:44s} got={g} exp={e}")
if not allpass:
    print("\n*** SELF-CHECK FAILED — writing NOTHING. ***"); sys.exit(1)

print("\n" + "=" * 78); print("WRITING"); print("=" * 78)
for df, fn in [(cand_sl, "candidate_shortlist.csv"), (repur, "repurposing_view.csv")]:
    p = OUT / fn; df.to_csv(p, index=False); print(f"  wrote {p}  ({len(df):,} rows, {df.shape[1]} cols)")
print("DONE.")
