"""STEP 2b Part A — Task 1 co_active refinement (relabel + backend-rigor columns).
Recomputes from source to INDEPENDENTLY assert FINAL_STATS Task-1 anchors. READ-ONLY sources."""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import numpy as np, pandas as pd
from statsmodels.stats.proportion import proportion_confint
import statsmodels

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))

np.random.seed(0)
print("pandas",pd.__version__,"statsmodels",statsmodels.__version__)
B=BASE; DB=B/"Output/DB/GPCRactDB/"; AN=B/"GPCRactDB/results/denovo/"
OUT=B/"Output/NAR_DBI/WebExport_v2/"; TA=OUT/"target_analysis"; L2_T=0.40
ACTIVE_MOA={"Binder","Agonist","Antagonist","Inverse Agonist","Partial Agonist","PAM","NAM"}
INACTIVE_MOA={"Inactive","Non-binder"}
LMAP={"Layer 1 (Binding)":"L1","Layer 2 (Proximal)":"L2","Layer 3 (Biased)":"L3","Layer 4 (Reporter)":"L4"}
m=pd.read_csv(DB/"Bioassay_MoA_Master_Integrated.csv",usecols=["Ligand_InChIKey","GPCR_UniProt","Assay_Layer","MoA"]).rename(
   columns={"Ligand_InChIKey":"InChIKey","GPCR_UniProt":"UniProt_AC"})
md=m[m.Assay_Layer.isin(LMAP)].copy(); md["layer"]=md.Assay_Layer.map(LMAP)
md["st"]=np.where(md.MoA.isin(ACTIVE_MOA),"Active",np.where(md.MoA.isin(INACTIVE_MOA),"Inactive","Unknown"))
def agg(s): return "Active" if (s=="Active").any() else ("Inactive" if (s=="Inactive").any() else "NA")
grp=md.groupby(["InChIKey","UniProt_AC","layer"]).st.agg(agg).unstack("layer")
for c in ["L1","L2"]:
    if c not in grp: grp[c]="NA"
grp=grp[["L1","L2"]].fillna("NA").reset_index()
G600=set(md.UniProt_AC.unique())
ti=pd.read_csv(DB/"Target_Indication_Master_Integrated.csv")
U=sorted(G600 & set(ti.UniProt_AC.unique()))
tl=pd.read_csv(OUT/"target_lookup.csv",usecols=["UniProt_AC","gene_name"]); gene=tl.set_index("UniProt_AC").gene_name.to_dict()
cat=pd.read_csv(DB/"Analysis/Compound_Chain_Catalog_v3.csv",low_memory=False,usecols=["InChIKey","is_approved"]).drop_duplicates("InChIKey")
appr=set(cat.loc[cat.is_approved==True,"InChIKey"])
di=pd.read_csv(DB/"Drug_Indication_Master_Integrated.csv"); appr_di=set(di.loc[di.Highest_Status=="Approved","InChIKey"])
annL2=pd.read_csv(AN/"synergpcr_annotation_L2.csv",usecols=["InChIKey","UniProt_AC","binding_prob","confidence_raw"])
ai2set=set(map(tuple,annL2[(annL2.binding_prob>=0.5)&(annL2.confidence_raw>=L2_T)][["InChIKey","UniProt_AC"]].values))
gu=grp[grp.UniProt_AC.isin(U)].copy()
gu["L1a"]=gu.L1=="Active"; gu["L2a"]=gu.L2=="Active"
gu["aiL2"]=[(ik,up) in ai2set for ik,up in zip(gu.InChIKey,gu.UniProt_AC)]
gu["co_exp"]=gu.L1a&gu.L2a; gu["co_aii"]=gu.L1a&(gu.L2a|gu.aiL2)
gu["appr"]=gu.InChIKey.isin(appr); gu["appr_di"]=gu.InChIKey.isin(appr_di)

def w2(k,n): return proportion_confint(k,n,0.05,"wilson") if n>0 else (np.nan,np.nan)
def w1u(k,n): return proportion_confint(k,n,0.10,"wilson")[1] if n>0 else np.nan   # one-sided 95% upper
nL1=gu[gu.L1a].groupby("UniProt_AC").size(); apprL1=gu[gu.L1a].groupby("UniProt_AC").appr.sum()
def pt(flag):
    s=gu[gu[flag]]; g=s.groupby("UniProt_AC"); return g.size(),g.appr.sum(),g.appr_di.sum()
ncE,nkE,nkE_di=pt("co_exp"); ncA,nkA,_=pt("co_aii"); nkA_di=gu[gu.co_aii].groupby("UniProt_AC").appr_di.sum()

# ---- INDEPENDENT ANCHOR RECOMPUTE + ASSERT (FINAL_STATS H) ----
sum_co=int(ncE.sum()); sum_appr=int(nkE.sum()); base_pool=sum_appr/sum_co
sum_l1=int(nL1.sum()); sum_l1ap=int(apprL1.sum()); l1_rate=sum_l1ap/sum_l1
poolE=base_pool; poolA=int(nkA.sum())/int(ncA.sum())
def lotos(nc,nk):
    tn,tk=nc.sum(),nk.sum(); return {t:((tk-nk.get(t,0))/(tn-nc.get(t,0)) if (tn-nc.get(t,0))>0 else np.nan) for t in U}
lotoE=lotos(ncE,nkE); lotoA=lotos(ncA,nkA)
enoughE=sum(ncE.get(t,0)>=5 for t in U)
# build rows
def status(nc,nk,log2):
    if nc<5: return "insufficient_data"
    if nk==0: return "depleted (upper-bounded)"
    return "enriched" if log2>0 else "depleted"
rows=[]
for t in U:
    ncE_t=int(ncE.get(t,0)); nkE_t=int(nkE.get(t,0)); nL1_t=int(nL1.get(t,0))
    ncA_t=int(ncA.get(t,0)); nkA_t=int(nkA.get(t,0))
    rcE=nkE_t/ncE_t if ncE_t else np.nan; rcA=nkA_t/ncA_t if ncA_t else np.nan
    loE,hiE=w2(nkE_t,ncE_t); loA,hiA=w2(nkA_t,ncA_t)
    blE=lotoE[t]; blA=lotoA[t]
    log2E=np.log2(rcE/blE) if (ncE_t>=5 and rcE>0 and blE>0) else np.nan
    log2A=np.log2(rcA/blA) if (ncA_t>=5 and rcA>0 and blA>0) else np.nan
    ci_rob_E=bool(ncE_t>=5 and not np.isnan(loE) and loE>poolE)
    ci_rob_A=bool(ncA_t>=5 and not np.isnan(loA) and loA>poolA)
    high_n_no_clin=bool(ncE_t>=50 and nkE_t==0)
    rows.append(dict(UniProt=t,GPCR_name=gene.get(t,""),
        n_tested_L1_exp=nL1_t,n_co_active_exp=ncE_t,n_approved_exp=nkE_t,
        approval_rate_raw=round(int(apprL1.get(t,0))/nL1_t,6) if nL1_t else np.nan,
        approval_rate_cond_exp=round(rcE,6) if ncE_t else np.nan,
        wilson_lo_exp=round(loE,6) if ncE_t else np.nan,wilson_hi_exp=round(hiE,6) if ncE_t else np.nan,
        wilson_upper_oneside_exp=round(w1u(nkE_t,ncE_t),6) if ncE_t else np.nan,
        baseline_loto_exp=round(blE,6) if blE==blE else np.nan,baseline_pool_exp=round(poolE,6),
        log2_enrichment_exp=round(log2E,4) if log2E==log2E else np.nan,
        enrichment_status_exp=status(ncE_t,nkE_t,log2E),
        n_co_active_aiincl=ncA_t,n_approved_aiincl=nkA_t,
        approval_rate_cond_aiincl=round(rcA,6) if ncA_t else np.nan,
        wilson_lo_aiincl=round(loA,6) if ncA_t else np.nan,wilson_hi_aiincl=round(hiA,6) if ncA_t else np.nan,
        baseline_loto_aiincl=round(blA,6) if blA==blA else np.nan,baseline_pool_aiincl=round(poolA,6),
        log2_enrichment_aiincl=round(log2A,4) if log2A==log2A else np.nan,
        enrichment_status_aiincl=status(ncA_t,nkA_t,log2A),
        enough_n_exp=ncE_t>=5,enough_n_aiincl=ncA_t>=5,
        high_n_no_clinical_exp=high_n_no_clin,
        ci_robust_exp_backend=ci_rob_E,ci_robust_aiincl_backend=ci_rob_A,
        n_approved_exp_drugind=int(nkE_di.get(t,0))))
T1=pd.DataFrame(rows)
nume=pd.to_numeric(T1.log2_enrichment_exp,errors='coerce')
point_enr=int((nume>0).sum()); comp=int(nume.notna().sum())
zero_app=int(((T1.n_co_active_exp>=5)&(T1.n_approved_exp==0)).sum())
ci_rob_t=int(T1.ci_robust_exp_backend.sum()); ci_rob_a=int(T1.ci_robust_aiincl_backend.sum())
uniq_appr=gu.loc[gu.co_exp & gu.appr,"InChIKey"].nunique()

print("\n--- ANCHOR ASSERTS vs FINAL_STATS ---")
A=[("universe==211",len(U)==211,len(U)),
   ("sum co_active pairs==20145",sum_co==20145,sum_co),
   ("sum approved-co pairs==282",sum_appr==282,sum_appr),
   ("baseline pooled==1.40%",round(base_pool*100,2)==1.40,round(base_pool*100,2)),
   ("L1-active pooled==0.86%",round(l1_rate*100,2)==0.86,round(l1_rate*100,2)),
   ("enrichment 1.63x",round(base_pool/l1_rate,2)==1.63,round(base_pool/l1_rate,2)),
   ("enough_n==126",enoughE==126,enoughE),
   ("computable log2==85",comp==85,comp),
   ("n>=5 & 0-approved==41",zero_app==41,zero_app),
   ("point-enriched==50",point_enr==50,point_enr),
   ("ci_robust target==28",ci_rob_t==28,ci_rob_t),
   ("ci_robust aiincl==31 (user-confirmed aiincl-vs-own-baseline)",ci_rob_a==31,ci_rob_a)]
allok=True
for n,ok,g in A:
    print(f"  [{'PASS' if ok else 'FAIL'}] {n:32s} got={g}"); allok&=ok
print(f"  unique-compound approved (co_active) = {uniq_appr}  (vs pair-level 282)")
# high-n 0-approved spot
for gn in ["MCHR1","ADORA3","GCGR"]:
    r=T1[T1.GPCR_name==gn]
    if len(r): print(f"  {gn}: n_co_active_exp={int(r.n_co_active_exp.iloc[0])} n_approved={int(r.n_approved_exp.iloc[0])} status={r.enrichment_status_exp.iloc[0]} high_n_flag={bool(r.high_n_no_clinical_exp.iloc[0])}")
if not allok:
    print("\n*** ANCHOR MISMATCH — writing NOTHING ***"); raise SystemExit(1)
T1.to_csv(TA/"target_clinical_enrichment.tsv",sep="\t",index=False)
print(f"\nWROTE target_clinical_enrichment.tsv ({len(T1)} rows, {T1.shape[1]} cols)")

# ---- disease TSV refine ----
umls=pd.read_csv(DB/"UMLS_CUI_to_Category_v2.csv"); cui2cat=umls.set_index("umls_cui").category_name.to_dict()
di["cat"]=di.UMLS_CUI.map(cui2cat)
appr_cat=di[(di.Highest_Status=="Approved")&di.cat.notna()&di.InChIKey.isin(appr)][["InChIKey","cat"]].drop_duplicates()
ik2cats=appr_cat.groupby("InChIKey").cat.agg(set).to_dict()
cats=sorted(umls.category_name.dropna().unique())
conc=gu[gu.co_exp][["InChIKey","UniProt_AC","appr"]].copy()
def hd(ik,d): return d in ik2cats.get(ik,set())
for d in cats: conc[f"n_{d}"]=[(a and hd(ik,d)) for ik,a in zip(conc.InChIKey,conc.appr)]
totN=conc.groupby("UniProt_AC").size(); cr=[]
for d in cats:
    numt=conc.groupby("UniProt_AC")[f"n_{d}"].sum(); tn=totN.sum(); tk=int(conc[f"n_{d}"].sum()); pool_d=tk/tn
    for t in U:
        nc_t=int(totN.get(t,0))
        if nc_t<5: continue
        k_t=int(numt.get(t,0)); loto_d=(tk-k_t)/(tn-nc_t) if (tn-nc_t)>0 else np.nan
        rate=k_t/nc_t; lo,hi=w2(k_t,nc_t); log2=np.log2(rate/loto_d) if (rate>0 and loto_d>0) else np.nan
        cr.append(dict(UniProt=t,GPCR_name=gene.get(t,""),disease_category=d,n_co_active_exp=nc_t,
            n_approved_cat=k_t,approval_rate_cond=round(rate,6),wilson_lo=round(lo,6),wilson_hi=round(hi,6),
            baseline_loto=round(loto_d,6) if loto_d==loto_d else np.nan,baseline_pool=round(pool_d,6),
            log2_enrichment=round(log2,4) if log2==log2 else np.nan,
            single_approved_flag_backend=bool(k_t==1),
            ci_robust_backend=bool(not np.isnan(lo) and lo>pool_d)))
TD=pd.DataFrame(cr); ci_rob_d=int(TD.ci_robust_backend.sum())
print(f"disease cells={len(TD)}  ci_robust(disease)={ci_rob_d}  (FINAL anchor 236: {'PASS' if ci_rob_d==236 else 'FAIL'})")
if ci_rob_d!=236:
    print("*** disease ci_robust mismatch — not writing disease tsv ***")
else:
    TD.to_csv(TA/"target_disease_enrichment.tsv",sep="\t",index=False)
    print(f"WROTE target_disease_enrichment.tsv ({len(TD)} rows, {TD.shape[1]} cols)")
print("DONE.")
