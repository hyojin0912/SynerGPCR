"""
webexport_step3_synergy.py — WebExport_v2 STEP 3. READ-ONLY on sources.

Generates THREE files, all RECALCULATED from PROVENANCE; builds in memory, runs
ALL self-checks, writes ONLY if every check passes.

  1. synergy_stats.csv        (experimental chain SI; PK chain_pattern)
  2. single_layer_table.csv   (Table A; experimental single-layer; PK layer)
  3. gpcract_3point_si.csv     (GPCRact panel @ 0.40/0.97; PK layer,point)

SI/E_Bliss/Wilson/null logic reused VERBATIM from part11_si_and_null_final.py /
part11_final_op_0p40_0p97.py. Cross-validated against synergy_stats_v3.csv.
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.proportion import proportion_confint

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


B = BASE
OUT = B / "Output/NAR_DBI/WebExport_v2"
CAT  = B / "Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv"
CHN  = B / "Output/DB/GPCRactDB/Analysis/Chain_Pattern_Distribution.csv"
PC   = B / "Output/NAR_L3L4/Tables/L3L4_Pair_Catalog_v1.csv"
MAST = B / "Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv"
L2P  = B / "GPCRactDB/results/denovo/synergpcr_annotation_L2.csv"
L3P  = B / "GPCRactDB/results/denovo/synergpcr_annotation_L3.csv"
V3   = B / "Output/NAR/Tables/synergy_stats_v3.csv"   # cross-check only

L2_T, L3_T = 0.40, 0.97
N_REPS, SEED = 10_000, 0
MOA_MAP = {"Agonist": "ag", "Partial Agonist": "ag", "PAM": "ag", "Antagonist": "an",
           "Inverse Agonist": "an", "NAM": "an", "Inactive": "inact"}


def compute_si(k, n, e): return float("nan") if (n == 0 or e == 0) else (k / n) / e
def wilson(k, n): return (0.0, 0.0) if n == 0 else proportion_confint(k, n, 0.05, "wilson")
def si_rate(r, e): return float("nan") if e == 0 else r / e

# layer-set parser for chain_pattern names
def layers_of(cp):
    if cp == "No_active_layer":
        return []
    if cp == "Full_chain":
        return ["L1", "L2", "L3", "L4"]
    if cp.endswith("_only"):
        return [cp.replace("_only", "")]
    return cp.split("+")


def load(p, **kw):
    if not p.exists():
        sys.exit(f"STOP: missing PROVENANCE {p}")
    df = pd.read_csv(p, **kw); print(f"  loaded {p.name}: {len(df):,} rows"); return df


print("=" * 80); print("LOADING PROVENANCE"); print("=" * 80)
dist = load(CHN)
cat  = load(CAT, low_memory=False).drop_duplicates("InChIKey")
cat["is_approved"] = (pd.to_numeric(cat.max_phase, errors="coerce") == 4)

# P_only per single layer (from Chain_Pattern_Distribution _only rows)
def p_only(L):
    r = dist[dist.chain_pattern == f"{L}_only"]
    return int(r.n_approved.iloc[0]) / int(r.n.iloc[0])
P = {L: p_only(L) for L in ["L1", "L2", "L3", "L4"]}
print(f"\nP_only: " + "  ".join(f"{L}={P[L]:.6f}" for L in P))

def e_bliss(layers):
    if not layers:
        return float("nan")
    prod = 1.0
    for L in layers:
        prod *= (1 - P[L])
    return 1 - prod

# ==========================================================================
# FILE 1 — synergy_stats.csv (regenerated from Chain_Pattern_Distribution)
# ==========================================================================
print("\n" + "=" * 80); print("FILE 1: synergy_stats.csv"); print("=" * 80)
rows = []
for _, r in dist.iterrows():
    cp = r.chain_pattern; n = int(r.n); k = int(r.n_approved)
    layers = layers_of(cp)
    e = e_bliss(layers)
    rate = k / n if n else float("nan")
    lo, hi = wilson(k, n)
    si = compute_si(k, n, e)
    si_lo = si_rate(lo, e); si_hi = si_rate(hi, e)
    # interpretation
    if not layers:
        interp = "baseline (no active layer)"
    elif len(layers) == 1:
        interp = "reference (single layer, SI=1)"
    elif cp == "Full_chain":
        interp = "REFERENCE-ONLY (L2-L3 dependent, Fisher OR=9.24)"
    else:
        interp = ("Synergistic" if si_lo > 1 else
                  "Sub-additive" if si_hi < 1 else "Additive")
    rows.append(dict(chain_pattern=cp, layers_in_combo=",".join(layers),
                     n_compounds=n, n_approved=k,
                     approval_rate_pct=round(rate * 100, 4),
                     wilson_ci_lo_pct=round(lo * 100, 4), wilson_ci_hi_pct=round(hi * 100, 4),
                     e_bliss=(round(e, 6) if e == e else np.nan),
                     bliss_SI=(round(si, 4) if si == si else np.nan),
                     bliss_SI_ci_lo=(round(si_lo, 4) if si_lo == si_lo else np.nan),
                     bliss_SI_ci_hi=(round(si_hi, 4) if si_hi == si_hi else np.nan),
                     bliss_interpretation=interp))
syn = pd.DataFrame(rows)
def get(cp, col): return syn.loc[syn.chain_pattern == cp, col].iloc[0]
print(syn[["chain_pattern", "n_compounds", "n_approved", "approval_rate_pct",
           "e_bliss", "bliss_SI", "bliss_interpretation"]].to_string(index=False))

# cross-validate vs synergy_stats_v3 (anchors)
v3 = load(V3)
def v3si(cp): return float(v3.loc[v3.chain_pattern == cp, "bliss_SI"].iloc[0])

# ==========================================================================
# FILE 2 — single_layer_table.csv (Table A; the 4 "_only" rows)
# ==========================================================================
print("\n" + "=" * 80); print("FILE 2: single_layer_table.csv (Table A)"); print("=" * 80)
slrows = []
for L in ["L1", "L2", "L3", "L4"]:
    r = dist[dist.chain_pattern == f"{L}_only"].iloc[0]
    n = int(r.n); k = int(r.n_approved); rate = k / n
    lo, hi = wilson(k, n); e = e_bliss([L])
    slrows.append(dict(layer=L, n_compounds=n, n_approved=k,
                       approval_rate_pct=round(rate * 100, 4),
                       wilson_ci_lo_pct=round(lo * 100, 4), wilson_ci_hi_pct=round(hi * 100, 4),
                       e_bliss=round(e, 6), bliss_SI=round(compute_si(k, n, e), 4),
                       bliss_SI_ci_lo=round(si_rate(lo, e), 4), bliss_SI_ci_hi=round(si_rate(hi, e), 4)))
sl = pd.DataFrame(slrows)
print(sl.to_string(index=False))

# ==========================================================================
# FILE 3 — gpcract_3point_si.csv (GPCRact panel @ 0.40/0.97; part11_final logic)
# ==========================================================================
print("\n" + "=" * 80); print("FILE 3: gpcract_3point_si.csv"); print("=" * 80)
A = lambda c: cat[c] == "Active"
cati = cat.set_index("InChIKey")
pc = load(PC, low_memory=False)
l2 = load(L2P); l3 = load(L3P)
E2 = e_bliss(["L1", "L2", "L4"]); E3 = e_bliss(["L1", "L3", "L4"])

def stat(ik):
    ik = list(ik)
    return (len(ik), int(cati.loc[ik].is_approved.sum())) if ik else (0, 0)

ceil2 = set(cat.loc[A("L1_state") & A("L2_state") & A("L4_state"), "InChIKey"])
ceil3 = set(cat.loc[A("L1_state") & A("L3_state") & A("L4_state"), "InChIKey"])
base2 = cat[A("L1_state") & A("L4_state") & cat.L2_state.isna()].drop_duplicates("InChIKey")
base3 = cat[A("L1_state") & A("L4_state") & cat.L3_state.isna()].drop_duplicates("InChIKey")
base2_ik, base3_ik = set(base2.InChIKey), set(base3.InChIKey)
n_c2, k_c2 = stat(ceil2); n_c3, k_c3 = stat(ceil3)

m = load(MAST, low_memory=False, usecols=["Ligand_InChIKey", "GPCR_UniProt", "Assay_Layer", "MoA"]).rename(
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

def mc_null(base_df, n_up, n_ceil, k_ceil, e):
    if n_up == 0: return float("nan"), float("nan")
    appr = base_df.is_approved.to_numpy().astype(int)
    rng = np.random.default_rng(SEED)
    N = n_ceil + n_up; idx = np.arange(len(appr)); sis = np.empty(N_REPS)
    for i in range(N_REPS):
        sis[i] = ((k_ceil + int(appr[rng.choice(idx, n_up, replace=False)].sum())) / N) / e
    return float(sis.mean()), float(np.percentile(sis, 95))

g3rows = []
def panel(layer, thr, base_ik, base_df, up_ik, ceil_ik, e):
    _nc, _kc = stat(ceil_ik)                       # stat returns (n, k)
    si_ceil = compute_si(_kc, _nc, e)              # compute_si wants (k, n, e)
    nm, n95 = mc_null(base_df, len(up_ik), _nc, _kc, e)
    for point, ik, t, n95v, recv in [
            ("baseline", base_ik, "baseline", np.nan, np.nan),
            ("gpcract", ceil_ik | up_ik, thr, n95, None),
            ("ceiling", ceil_ik, "ceiling", np.nan, np.nan)]:
        n, k = stat(ik); rate = k / n if n else float("nan")
        lo, hi = wilson(k, n); si = compute_si(k, n, e)
        rec = (si / si_ceil * 100) if point == "gpcract" else (100.0 if point == "ceiling" else np.nan)
        g3rows.append(dict(layer=layer, point=point, threshold=t, n_pool=n, n_approved=k,
                           obs_rate=round(rate, 6), e_bliss=round(e, 6), SI=round(si, 4),
                           SI_ci_lo=round(si_rate(lo, e), 4), SI_ci_hi=round(si_rate(hi, e), 4),
                           null95=(round(n95v, 4) if n95v == n95v else np.nan),
                           recovery_pct=(round(rec, 1) if rec == rec else np.nan)))
    return n95

n95_2 = panel("L2", L2_T, base2_ik, base2, up_l2, ceil2, E2)
n95_3 = panel("L3", L3_T, base3_ik, base3, up_l3, ceil3, E3)
g3 = pd.DataFrame(g3rows)
print(g3.to_string(index=False))

def g3v(layer, point, col): return g3[(g3.layer == layer) & (g3.point == point)][col].iloc[0]

# ==========================================================================
# SELF-CHECKS
# ==========================================================================
print("\n" + "=" * 80); print("SELF-CHECKS"); print("=" * 80)
checks = []
def ck(name, cond, got, exp): checks.append((name, bool(cond), got, exp))

# synergy_stats anchors (vs §3 and vs v3)
for cp, en, ek, esi in [("L1+L2+L4", 2580, 165, 9.71), ("L1+L3+L4", 331, 43, 13.01),
                        ("Full_chain", 228, 37, 13.45)]:
    ck(f"synergy {cp} n=={en}", int(get(cp, "n_compounds")) == en, int(get(cp, "n_compounds")), en)
    ck(f"synergy {cp} k=={ek}", int(get(cp, "n_approved")) == ek, int(get(cp, "n_approved")), ek)
    ck(f"synergy {cp} SI~={esi}", round(float(get(cp, "bliss_SI")), 2) == esi, round(float(get(cp, "bliss_SI")), 2), esi)
    ck(f"synergy {cp} == v3", abs(float(get(cp, "bliss_SI")) - v3si(cp)) < 0.01, round(float(get(cp,'bliss_SI')),4), round(v3si(cp),4))
ck("Full_chain label REFERENCE-ONLY", "REFERENCE-ONLY" in get("Full_chain", "bliss_interpretation"),
   get("Full_chain", "bliss_interpretation"), "contains REFERENCE-ONLY")
# all synergy rows match v3 SI
mism = []
for cp in syn.chain_pattern:
    a = float(get(cp, "bliss_SI"));
    try: bv = v3si(cp)
    except Exception: continue
    if (a == a) and (bv == bv) and abs(a - bv) >= 0.01: mism.append((cp, a, bv))
ck("ALL synergy rows match v3 SI", len(mism) == 0, mism, "[]")

# single_layer Table A: rates < 1%, SI == 1.0
ck("Table A all rates < 1%", bool((sl.approval_rate_pct < 1.0).all()), sl.approval_rate_pct.tolist(), "<1.0 each")
ck("Table A all SI == 1.0", bool((sl.bliss_SI.round(4) == 1.0).all()), sl.bliss_SI.tolist(), "1.0 each")

# gpcract 3-point
ck("L2 baseline SI ~1.58", round(float(g3v("L2","baseline","SI")),2)==1.58, round(float(g3v("L2","baseline","SI")),2),1.58)
ck("L2 gpcract SI ~7.67", round(float(g3v("L2","gpcract","SI")),2)==7.67, round(float(g3v("L2","gpcract","SI")),2),7.67)
ck("L2 ceiling SI ~9.71", round(float(g3v("L2","ceiling","SI")),2)==9.71, round(float(g3v("L2","ceiling","SI")),2),9.71)
ck("L2 null95 ~7.17", round(float(n95_2),2)==7.17, round(float(n95_2),2),7.17)
ck("L2 recovery ~79.0", round(float(g3v("L2","gpcract","recovery_pct")),1)==79.0, round(float(g3v("L2","gpcract","recovery_pct")),1),79.0)
ck("L3 baseline SI ~2.20", round(float(g3v("L3","baseline","SI")),2)==2.20, round(float(g3v("L3","baseline","SI")),2),2.20)
ck("L3 gpcract SI ~10.35", round(float(g3v("L3","gpcract","SI")),2)==10.35, round(float(g3v("L3","gpcract","SI")),2),10.35)
ck("L3 ceiling SI ~13.01", round(float(g3v("L3","ceiling","SI")),2)==13.01, round(float(g3v("L3","ceiling","SI")),2),13.01)
ck("L3 null95 ~7.36", round(float(n95_3),2)==7.36, round(float(n95_3),2),7.36)
ck("L3 recovery ~79.6", round(float(g3v("L3","gpcract","recovery_pct")),1)==79.6, round(float(g3v("L3","gpcract","recovery_pct")),1),79.6)

# §5 STALE scan — forbidden values are HEADLINE point estimates (SI / null), not CI bounds;
# Wilson CI bounds legitimately span ranges and may coincide with a forbidden number.
FORBID = {8.26, 7.23, 7.81, 13.31, 3.64, 3.67, 4.93}
POINT_COLS = {"bliss_SI", "SI", "null95", "e_bliss", "recovery_pct"}
stale_hits = []
for df in (syn, sl, g3):
    for c in df.columns:
        if c in POINT_COLS and df[c].dtype.kind in "fi":
            for v in df[c].dropna().round(2).tolist():
                if v in FORBID: stale_hits.append((c, v))
        if df[c].dtype == object and df[c].astype(str).str.contains("permutation|corrected", case=False, regex=True).any():
            stale_hits.append((c, "permutation/corrected term"))
ck("§5 STALE scan clean (point SI/null; no 8.26/7.23/7.81/13.31/3.64/3.67/4.93/permutation/corrected)",
   len(stale_hits) == 0, stale_hits, "[]")
ck("STALE guard: 7.23 not in any SI", 7.23 not in [round(float(x),2) for x in g3.SI.tolist()], [round(float(x),2) for x in g3.SI.tolist()], "no 7.23")

allpass = all(c[1] for c in checks)
for name, ok, got, exp in checks:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:46s} got={got} exp={exp}")
if not allpass:
    print("\n*** SELF-CHECK FAILED — writing NOTHING. ***"); sys.exit(1)

print("\n" + "=" * 80); print("ALL CHECKS PASS — WRITING"); print("=" * 80)
OUT.mkdir(parents=True, exist_ok=True)
for df, fn in [(syn, "synergy_stats.csv"), (sl, "single_layer_table.csv"), (g3, "gpcract_3point_si.csv")]:
    p = OUT / fn; df.to_csv(p, index=False); print(f"  wrote {p}  ({len(df):,} rows, {df.shape[1]} cols)")
print("\nDONE.")
