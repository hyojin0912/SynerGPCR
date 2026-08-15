"""STEP 2c — high_n_no_clinical target validation + re-label. READ-ONLY sources.
Cross-checks 'does this target have an approved drug' INDEPENDENTLY of the co_active set,
adds backend columns, re-asserts FINAL_STATS anchors unchanged. Writes only target_analysis/."""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import numpy as np, pandas as pd

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))

np.random.seed(0)
print("pandas", pd.__version__)
B=BASE; DB=B/"Output/DB/GPCRactDB/"
TA=B/"Output/NAR_DBI/WebExport_v2/target_analysis"

# ── existing Task-1 table (do NOT recompute enrichment; only add columns) ──
T1=pd.read_csv(TA/"target_clinical_enrichment.tsv",sep="\t")
print(f"loaded target_clinical_enrichment.tsv: {len(T1)} rows x {T1.shape[1]} cols")

# ── independent approved-drug evidence per target ──
# (A) Drug_Indication 'Approved' drugs ⋈ Target_Indication (drug shares an indication-target link)
ti=pd.read_csv(DB/"Target_Indication_Master_Integrated.csv")   # UniProt_AC, UMLS_CUI, Source
di=pd.read_csv(DB/"Drug_Indication_Master_Integrated.csv")     # InChIKey, UMLS_CUI, Source, Highest_Status
appr_di_drugs=set(di.loc[di.Highest_Status=="Approved","InChIKey"])
# drugs approved & sharing an indication CUI with the target
appr_di_pairs=di[di.Highest_Status=="Approved"][["InChIKey","UMLS_CUI"]].merge(
    ti[["UniProt_AC","UMLS_CUI"]],on="UMLS_CUI")
tgt_has_appr_indication=appr_di_pairs.groupby("UniProt_AC").InChIKey.nunique()

# (B) Catalog max_phase==4 compounds that have ANY activity record at the target (layer-agnostic, Bioassay)
cat=pd.read_csv(DB/"Analysis/Compound_Chain_Catalog_v3.csv",low_memory=False,
                usecols=["InChIKey","is_approved"]).drop_duplicates("InChIKey")
appr_cat=set(cat.loc[cat.is_approved==True,"InChIKey"])
m=pd.read_csv(DB/"Bioassay_MoA_Master_Integrated.csv",
              usecols=["Ligand_InChIKey","GPCR_UniProt"]).rename(
              columns={"Ligand_InChIKey":"InChIKey","GPCR_UniProt":"UniProt_AC"})
m_appr=m[m.InChIKey.isin(appr_cat)]
tgt_appr_anyrecord=m_appr.groupby("UniProt_AC").InChIKey.nunique()        # any layer, any MoA
# (B') restrict to ACTIVE records only (layer-agnostic) for a stricter "target-active approved drug"
LMAP={"Layer 1 (Binding)":"L1","Layer 2 (Proximal)":"L2","Layer 3 (Biased)":"L3","Layer 4 (Reporter)":"L4"}
ACTIVE_MOA={"Binder","Agonist","Antagonist","Inverse Agonist","Partial Agonist","PAM","NAM"}
m2=pd.read_csv(DB/"Bioassay_MoA_Master_Integrated.csv",
               usecols=["Ligand_InChIKey","GPCR_UniProt","Assay_Layer","MoA"]).rename(
               columns={"Ligand_InChIKey":"InChIKey","GPCR_UniProt":"UniProt_AC"})
m2a=m2[m2.InChIKey.isin(appr_cat) & m2.MoA.isin(ACTIVE_MOA)]
tgt_appr_active=m2a.groupby("UniProt_AC").InChIKey.nunique()

# ── attach evidence ──
T1["_appr_indication"]=T1.UniProt.map(tgt_has_appr_indication).fillna(0).astype(int)
T1["_appr_anyrecord"]=T1.UniProt.map(tgt_appr_anyrecord).fillna(0).astype(int)
T1["_appr_active_anylayer"]=T1.UniProt.map(tgt_appr_active).fillna(0).astype(int)
# approved drug EXISTS for target if any independent evidence stream is positive
T1["approved_drug_exists_any"]=(T1._appr_indication>0)|(T1._appr_anyrecord>0)|(T1._appr_active_anylayer>0)
# approved exists but co_active missed it (only meaningful where co_active found 0 approved)
T1["approved_drug_not_in_coactive"]=T1.approved_drug_exists_any & (T1.n_approved_exp==0)
# clinical_status_label only for enough_n & 0-approved targets; else NA
zero_set=(T1.enough_n_exp==True)&(T1.n_approved_exp==0)
def label(r):
    if not (r.enough_n_exp and r.n_approved_exp==0): return ""
    return "approved_outside_coactive" if r.approved_drug_exists_any else "no_approved_drug"
T1["clinical_status_label"]=T1.apply(label,axis=1)

print(f"\nenough_n & 0-approved targets = {int(zero_set.sum())} (expect 41)")
z=T1[zero_set]
print(f"  approved_outside_coactive = {int((z.clinical_status_label=='approved_outside_coactive').sum())}")
print(f"  no_approved_drug          = {int((z.clinical_status_label=='no_approved_drug').sum())}")
hz=z[z.high_n_no_clinical_exp==True]
print(f"  of high_n_no_clinical(24): outside_coactive={int((hz.clinical_status_label=='approved_outside_coactive').sum())} "
      f"no_approved={int((hz.clinical_status_label=='no_approved_drug').sum())}")

# ── explicit P2RY12 / TACR1 check ──
print("\n-- explicit verification --")
for gn in ["P2RY12","TACR1","MCHR1","ADORA3","GCGR"]:
    r=T1[T1.GPCR_name==gn]
    if len(r):
        r=r.iloc[0]
        print(f"  {gn:8s} ({r.UniProt}): n_co_active={int(r.n_co_active_exp)} n_appr_coactive={int(r.n_approved_exp)} "
              f"| indic={int(r._appr_indication)} anyrec={int(r._appr_anyrecord)} active={int(r._appr_active_anylayer)} "
              f"-> {r.clinical_status_label or 'n/a'}")

# ── re-assert FINAL_STATS anchors unchanged (enrichment columns untouched) ──
print("\n-- anchor re-assert (must be unchanged) --")
nume=pd.to_numeric(T1.log2_enrichment_exp,errors='coerce')
checks=[("universe==211",len(T1)==211),
        ("sum n_co_active_exp==20145",int(T1.n_co_active_exp.sum())==20145),
        ("sum n_approved_exp==282",int(T1.n_approved_exp.sum())==282),
        ("enough_n==126",int(T1.enough_n_exp.sum())==126),
        ("point-enriched==50",int((nume>0).sum())==50),
        ("ci_robust target==28",int(T1.ci_robust_exp_backend.sum())==28),
        ("ci_robust aiincl==31",int(T1.ci_robust_aiincl_backend.sum())==31),
        ("n>=5 & 0-approved==41",int(zero_set.sum())==41)]
allok=True
for n,ok in checks: print(f"  [{'PASS' if ok else 'FAIL'}] {n}"); allok&=ok
# drop helper underscore cols except keep the 3 evidence counts as backend (rename clean)
T1=T1.rename(columns={"_appr_indication":"approved_indication_links_backend",
                      "_appr_anyrecord":"approved_anyrecord_at_target_backend",
                      "_appr_active_anylayer":"approved_active_at_target_backend"})
if not allok:
    print("\n*** ANCHOR MISMATCH — NOT writing ***"); raise SystemExit(1)
T1.to_csv(TA/"target_clinical_enrichment.tsv",sep="\t",index=False)
print(f"\nWROTE target_clinical_enrichment.tsv ({len(T1)} rows, {T1.shape[1]} cols)")
print("DONE.")
