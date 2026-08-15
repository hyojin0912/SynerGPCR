"""STEP 2 — Task 1 target clinical enrichment. READ-ONLY sources; writes only to target_analysis/.
Hypothesis-generating; conditional approval rate with LOTO + pool baselines. No causal claims."""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import numpy as np, pandas as pd
from statsmodels.stats.proportion import proportion_confint
import statsmodels

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))

np.random.seed(0)
print("pandas", pd.__version__, "numpy", np.__version__, "statsmodels", statsmodels.__version__)

B = BASE
DB = B/"Output/DB/GPCRactDB/"; AN = B/"GPCRactDB/results/denovo/"
OUT = B/"Output/NAR_DBI/WebExport_v2/"; TA = OUT/"target_analysis"
L2_T = 0.40

# ── per-(compound,GPCR,layer) experimental state (locked: Bioassay + part11c map, any-active->Active)
ACTIVE_MOA={"Binder","Agonist","Antagonist","Inverse Agonist","Partial Agonist","PAM","NAM"}
INACTIVE_MOA={"Inactive","Non-binder"}
LMAP={"Layer 1 (Binding)":"L1","Layer 2 (Proximal)":"L2","Layer 3 (Biased)":"L3","Layer 4 (Reporter)":"L4"}
m = pd.read_csv(DB/"Bioassay_MoA_Master_Integrated.csv",
                usecols=["Ligand_InChIKey","GPCR_UniProt","Assay_Layer","MoA"]).rename(
                columns={"Ligand_InChIKey":"InChIKey","GPCR_UniProt":"UniProt_AC"})
md = m[m.Assay_Layer.isin(LMAP)].copy(); md["layer"]=md.Assay_Layer.map(LMAP)
md["st"]=np.where(md.MoA.isin(ACTIVE_MOA),"Active",np.where(md.MoA.isin(INACTIVE_MOA),"Inactive","Unknown"))
def agg(s): return "Active" if (s=="Active").any() else ("Inactive" if (s=="Inactive").any() else "NA")
grp = md.groupby(["InChIKey","UniProt_AC","layer"]).st.agg(agg).unstack("layer")
for c in ["L1","L2"]:
    if c not in grp: grp[c]="NA"
grp = grp[["L1","L2"]].fillna("NA").reset_index()

# ── universe = 600-scope ∩ Target_Indication(221)
G600 = set(md.UniProt_AC.unique())
ti = pd.read_csv(DB/"Target_Indication_Master_Integrated.csv")
U = sorted(G600 & set(ti.UniProt_AC.unique()))
print(f"universe targets = {len(U)}")
tl = pd.read_csv(OUT/"target_lookup.csv", usecols=["UniProt_AC","gene_name"])
gene = tl.set_index("UniProt_AC").gene_name.to_dict()

# ── approved sets
cat = pd.read_csv(DB/"Analysis/Compound_Chain_Catalog_v3.csv", low_memory=False,
                  usecols=["InChIKey","is_approved"]).drop_duplicates("InChIKey")
appr = set(cat.loc[cat.is_approved==True,"InChIKey"])         # max_phase==4 (primary)
di = pd.read_csv(DB/"Drug_Indication_Master_Integrated.csv")
appr_di = set(di.loc[di.Highest_Status=="Approved","InChIKey"])  # DrugInd (sensitivity)

# ── AI-active L2 (binder-gated) pairs
annL2 = pd.read_csv(AN/"synergpcr_annotation_L2.csv",
                    usecols=["InChIKey","UniProt_AC","binding_prob","confidence_raw"])
ai2 = annL2[(annL2.binding_prob>=0.5)&(annL2.confidence_raw>=L2_T)][["InChIKey","UniProt_AC"]]
ai2set = set(map(tuple, ai2.values))

# ── restrict to universe, build pair-level flags
gu = grp[grp.UniProt_AC.isin(U)].copy()
gu["L1a"] = gu.L1=="Active"; gu["L2a"] = gu.L2=="Active"
gu["aiL2"] = [ (ik,up) in ai2set for ik,up in zip(gu.InChIKey, gu.UniProt_AC) ]
gu["conc_exp"]   = gu.L1a & gu.L2a
gu["conc_aii"]   = gu.L1a & (gu.L2a | gu.aiL2)
gu["appr"]   = gu.InChIKey.isin(appr)
gu["appr_di"]= gu.InChIKey.isin(appr_di)

# per-target aggregates (compound already unique per (InChIKey,UniProt_AC))
def per_target(flag):
    sub = gu[gu[flag]]
    g = sub.groupby("UniProt_AC")
    return g.size(), g.appr.sum(), g.appr_di.sum()
nL1 = gu[gu.L1a].groupby("UniProt_AC").size()
ncE, nkE, nkE_di = per_target("conc_exp")
ncA, nkA, _      = per_target("conc_aii")
nkA_di = gu[gu.conc_aii].groupby("UniProt_AC").appr_di.sum()
apprL1 = gu[gu.L1a].groupby("UniProt_AC").appr.sum()

# pooled / LOTO baselines (stratum = concordant pairs across universe)
def baselines(nc, nk):
    tot_n = nc.sum(); tot_k = nk.sum()
    pool = tot_k/tot_n if tot_n else np.nan
    loto = {t: ((tot_k-nk.get(t,0))/(tot_n-nc.get(t,0)) if (tot_n-nc.get(t,0))>0 else np.nan) for t in U}
    return pool, loto
poolE, lotoE = baselines(ncE, nkE)
poolA, lotoA = baselines(ncA, nkA)
print(f"baseline_pool exp={poolE:.5f}  aiincl={poolA:.5f}")

def wilson(k,n): return proportion_confint(k,n,0.05,"wilson") if n>0 else (np.nan,np.nan)
rows=[]
for t in U:
    ncE_t=int(ncE.get(t,0)); nkE_t=int(nkE.get(t,0)); nL1_t=int(nL1.get(t,0))
    ncA_t=int(ncA.get(t,0)); nkA_t=int(nkA.get(t,0))
    rrate = (int(apprL1.get(t,0))/nL1_t) if nL1_t else np.nan
    enoughE = ncE_t>=5; enoughA = ncA_t>=5
    rcE = nkE_t/ncE_t if ncE_t else np.nan
    rcA = nkA_t/ncA_t if ncA_t else np.nan
    loE,hiE = wilson(nkE_t,ncE_t); loA,hiA = wilson(nkA_t,ncA_t)
    bl_lotoE = lotoE[t]; bl_lotoA = lotoA[t]
    log2E = np.log2(rcE/bl_lotoE) if (enoughE and rcE>0 and bl_lotoE and bl_lotoE>0) else np.nan
    log2A = np.log2(rcA/bl_lotoA) if (enoughA and rcA>0 and bl_lotoA and bl_lotoA>0) else np.nan
    rows.append(dict(
        UniProt=t, GPCR_name=gene.get(t,""),
        n_tested_L1_exp=nL1_t, n_concordant_exp=ncE_t, n_approved_exp=nkE_t,
        approval_rate_raw=round(rrate,6) if nL1_t else np.nan,
        approval_rate_cond_exp=round(rcE,6) if ncE_t else np.nan,
        wilson_lo_exp=round(loE,6) if ncE_t else np.nan, wilson_hi_exp=round(hiE,6) if ncE_t else np.nan,
        baseline_loto_exp=round(bl_lotoE,6) if bl_lotoE==bl_lotoE else np.nan,
        baseline_pool_exp=round(poolE,6),
        log2_enrichment_exp=round(log2E,4) if log2E==log2E else "insufficient data" if not enoughE else np.nan,
        n_concordant_aiincl=ncA_t, n_approved_aiincl=nkA_t,
        approval_rate_cond_aiincl=round(rcA,6) if ncA_t else np.nan,
        wilson_lo_aiincl=round(loA,6) if ncA_t else np.nan, wilson_hi_aiincl=round(hiA,6) if ncA_t else np.nan,
        log2_enrichment_aiincl=round(log2A,4) if log2A==log2A else "insufficient data" if not enoughA else np.nan,
        enough_n_exp=enoughE, enough_n_aiincl=enoughA,
        n_approved_exp_drugind=int(nkE_di.get(t,0)),   # sensitivity cross-check
        baseline_loto_aiincl=round(bl_lotoA,6) if bl_lotoA==bl_lotoA else np.nan,
        baseline_pool_aiincl=round(poolA,6)))
T1 = pd.DataFrame(rows)
T1.to_csv(TA/"target_clinical_enrichment.tsv", sep="\t", index=False)
print(f"\nWROTE target_clinical_enrichment.tsv  rows={len(T1)}")
print(f"  enough_n_exp True={int(T1.enough_n_exp.sum())}  aiincl True={int(T1.enough_n_aiincl.sum())}")
nume=pd.to_numeric(T1.log2_enrichment_exp,errors='coerce')
print(f"  log2_enrichment_exp computed n={nume.notna().sum()} | enriched(>0)={int((nume>0).sum())} depleted(<0)={int((nume<0).sum())}")
print("  top5 enriched (exp):")
print(T1.loc[nume.sort_values(ascending=False).head(5).index, ["UniProt","GPCR_name","n_concordant_exp","n_approved_exp","approval_rate_cond_exp","baseline_loto_exp","log2_enrichment_exp"]].to_string(index=False))

# ── B. target × disease cells (n_concordant_exp>=5)
umls = pd.read_csv(DB/"UMLS_CUI_to_Category_v2.csv")
cui2cat = umls.set_index("umls_cui").category_name.to_dict()
di["cat"] = di.UMLS_CUI.map(cui2cat)
di_appr_cat = di[(di.Highest_Status=="Approved") & di.cat.notna()][["InChIKey","cat"]].drop_duplicates()
# map: approved-by-maxphase compounds with category-d indication
appr_cat = di_appr_cat[di_appr_cat.InChIKey.isin(appr)]   # primary approved = max_phase==4
ik2cats = appr_cat.groupby("InChIKey").cat.agg(set).to_dict()
cats = sorted(umls.category_name.dropna().unique())
# concordant_exp pairs in universe
conc = gu[gu.conc_exp][["InChIKey","UniProt_AC","appr"]].copy()
cellrows=[]
# pooled baseline per disease across all targets (stratum = concordant pairs; numerator = approved w/ cat d)
def has_d(ik,d): return d in ik2cats.get(ik,set())
for d in cats:
    conc[f"num_{d}"] = [ (a and has_d(ik,d)) for ik,a in zip(conc.InChIKey, conc.appr) ]
totN = conc.groupby("UniProt_AC").size()
for d in cats:
    numt = conc.groupby("UniProt_AC")[f"num_{d}"].sum()
    tot_n=totN.sum(); tot_k=int(conc[f"num_{d}"].sum())
    pool_d = tot_k/tot_n if tot_n else np.nan
    for t in U:
        nc_t=int(totN.get(t,0))
        if nc_t<5: continue
        k_t=int(numt.get(t,0))
        loto_d = (tot_k-k_t)/(tot_n-nc_t) if (tot_n-nc_t)>0 else np.nan
        rate=k_t/nc_t; lo,hi=wilson(k_t,nc_t)
        log2=np.log2(rate/loto_d) if (rate>0 and loto_d and loto_d>0) else np.nan
        cellrows.append(dict(UniProt=t,GPCR_name=gene.get(t,""),disease_category=d,
            n_concordant_exp=nc_t,n_approved_cat=k_t,approval_rate_cond=round(rate,6),
            wilson_lo=round(lo,6),wilson_hi=round(hi,6),baseline_loto=round(loto_d,6) if loto_d==loto_d else np.nan,
            baseline_pool=round(pool_d,6),log2_enrichment=round(log2,4) if log2==log2 else np.nan))
TD = pd.DataFrame(cellrows)
TD.to_csv(TA/"target_disease_enrichment.tsv", sep="\t", index=False)
print(f"\nWROTE target_disease_enrichment.tsv  rows={len(TD)}  (cells with n_concordant_exp>=5)")
print(f"  unique targets in cells={TD.UniProt.nunique()}  categories={TD.disease_category.nunique()}")
print("DONE.")
