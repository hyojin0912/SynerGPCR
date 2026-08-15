"""STEP 2d (final) — on-target APPROVED-drug reclassification, TDL-gated. READ-ONLY sources.
On-target = a target with an approved drug acting via its MoA:
  (1) DrugCentral drug.target.interaction TDL=='Tclin' targets, MoA drugs (ACTION_TYPE present), OR
  (2) ChEMBL drug_mechanism drug whose InChIKey is approved.
Modality (small molecule vs peptide/biologic): IUPHAR Type primary; SMILES heuristic fallback."""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import re
import numpy as np, pandas as pd

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))

np.random.seed(0); print("pandas", pd.__version__)
B=BASE; DB=B/"Output/DB/GPCRactDB/"; ODB=B/"Output/DB/"
TA=B/"Output/NAR_DBI/WebExport_v2/target_analysis"
T1=pd.read_csv(TA/"target_clinical_enrichment.tsv",sep="\t")
for c in ["on_target_approved_n","on_target_smallmol_n","on_target_peptide_biologic_n","on_target_sources","clinical_status_label_v2"]:
    if c in T1.columns: T1=T1.drop(columns=c)
universe=set(T1.UniProt)

# APPROVED inchikeys (small mol + peptide/biologic curated lists)
cat=pd.read_csv(DB/"Analysis/Compound_Chain_Catalog_v3.csv",low_memory=False,usecols=["InChIKey","is_approved"]).drop_duplicates("InChIKey")
appr_cat=set(cat.loc[cat.is_approved==True,"InChIKey"])
iup_drugs=pd.read_csv(ODB/"IUPHAR/NAR/IUPHAR_GPCR_Approved_Drugs.csv")
dc=pd.read_csv(ODB/"DrugCentral/NAR/DrugCentral_GPCR_Master_v2.csv")
APPROVED=appr_cat|set(iup_drugs.InChIKey.dropna())|set(dc.inchikey.dropna())

# TDL=='Tclin' targets (DrugCentral target development level = has approved on-target drug)
dti=pd.read_csv(B/"DB/DrugCentral/drug.target.interaction.tsv",sep="\t")
dti.columns=[c.strip('"') for c in dti.columns]
tclin=set(dti.loc[dti.TDL=="Tclin","ACCESSION"].dropna())
print(f"DrugCentral Tclin targets (all): {len(tclin):,}  ∩ universe: {len(tclin & universe)}")

# modality
BIO={"Peptide","Antibody"}; mod_iuphar={k:("peptide_biologic" if t in BIO else "small_molecule") for k,t in zip(iup_drugs.InChIKey,iup_drugs.Type)}
smiles_map={}
def smod(smi):
    if not isinstance(smi,str) or smi.strip()=="": return "peptide_biologic"
    if len(re.findall(r"C\(=O\)N|NC\(=O\)",smi))>=4 and len(smi)>150: return "peptide_biologic"
    return "small_molecule"
def modality(ik): return mod_iuphar[ik] if ik in mod_iuphar else smod(smiles_map.get(ik))

ontarget={}; src_flag={}
def add(up,ik,src):
    if not isinstance(up,str) or not isinstance(ik,str) or up not in universe: return
    ontarget.setdefault(up,set()).add(ik); src_flag.setdefault((up,ik),set()).add(src)
# (1) DrugCentral MoA drugs at Tclin targets
dc_t=dc[dc.ACTION_TYPE.notna() & (dc.ACTION_TYPE.astype(str).str.strip()!="") & dc.ACCESSION.isin(tclin)]
for ik,up,smi in zip(dc_t.inchikey,dc_t.ACCESSION,dc_t.smiles):
    add(up,ik,"drugcentral_Tclin");
    if isinstance(ik,str): smiles_map.setdefault(ik,smi)
# (2) ChEMBL drug_mechanism, approved-gated (curated MoA)
moa=pd.read_csv(ODB/"ChEMBL/NAR/ChEMBL_v36_MoA.csv")
t2u=pd.read_csv(ODB/"ChEMBL/NAR/ChEMBL_v36_Target_Metadata.csv").dropna(subset=["uniprot_accession"]).set_index("target_chembl_id").uniprot_accession.to_dict()
moa["uniprot"]=moa.target_chembl_id.map(t2u); moa_a=moa[moa.standard_inchi_key.isin(APPROVED)]
for ik,up,smi in zip(moa_a.standard_inchi_key,moa_a.uniprot,moa_a.canonical_smiles):
    add(up,ik,"chembl");
    if isinstance(ik,str): smiles_map.setdefault(ik,smi)
print(f"DrugCentral Tclin MoA rows={len(dc_t):,}  ChEMBL approved MoA rows={len(moa_a):,}")

rows=[]
for up in T1.UniProt:
    iks=ontarget.get(up,set()); sm=sum(modality(ik)=="small_molecule" for ik in iks); pb=len(iks)-sm
    srcs=set().union(*[src_flag.get((up,ik),set()) for ik in iks]) if iks else set()
    rows.append(dict(UniProt=up,on_target_approved_n=len(iks),on_target_smallmol_n=sm,on_target_peptide_biologic_n=pb,on_target_sources=";".join(sorted(srcs))))
T1=T1.merge(pd.DataFrame(rows).set_index("UniProt"),left_on="UniProt",right_index=True,how="left")
zero=(T1.enough_n_exp==True)&(T1.n_approved_exp==0)
def v2(r):
    if not (r.enough_n_exp and r.n_approved_exp==0): return ""
    if r.on_target_smallmol_n>0: return "approved_on_target_smallmol_outside_coactive"
    if r.on_target_peptide_biologic_n>0: return "approved_on_target_peptide_or_biologic"
    return "no_on_target_approved_drug"
T1["clinical_status_label_v2"]=T1.apply(v2,axis=1)
ren={"approved_indication_links_backend":"disease_indication_links_noise_backend",
     "approved_anyrecord_at_target_backend":"any_record_at_target_offtarget_noise_backend",
     "approved_active_at_target_backend":"active_record_at_target_offtarget_noise_backend",
     "clinical_status_label":"clinical_status_label_v1_superseded"}
T1=T1.rename(columns={k:v for k,v in ren.items() if k in T1.columns})

vc=T1[zero].clinical_status_label_v2.value_counts()
print("\n-- v2 over 41 enough_n & 0-approved --")
for k in ["approved_on_target_smallmol_outside_coactive","approved_on_target_peptide_or_biologic","no_on_target_approved_drug"]:
    print(f"  {k:48s}={int(vc.get(k,0))}")
print("\n-- explicit verification (expected: P2RY12/TACR1 smallmol; GCGR/CALCR peptide; MCHR1/ADORA3/FFAR1 no_on_target) --")
exp={"P2RY12":"smallmol","TACR1":"smallmol","GCGR":"peptide","CALCR":"peptide","MCHR1":"no","ADORA3":"no","FFAR1":"no"}
for gn in exp:
    r=T1[T1.GPCR_name==gn]
    if len(r):
        r=r.iloc[0]; lab=r.clinical_status_label_v2
        ok = (("smallmol" in lab and exp[gn]=="smallmol") or ("peptide" in lab and exp[gn]=="peptide") or ("no_on_target" in lab and exp[gn]=="no"))
        print(f"  [{'OK' if ok else 'XX'}] {gn:8s}({r.UniProt}): n={int(r.on_target_approved_n)} sm={int(r.on_target_smallmol_n)} pb={int(r.on_target_peptide_biologic_n)} src=[{r.on_target_sources}] -> {lab or 'n/a'}")
print("\n-- anchor re-assert --")
nume=pd.to_numeric(T1.log2_enrichment_exp,errors='coerce')
ck=[("universe==211",len(T1)==211),("sum n_co_active_exp==20145",int(T1.n_co_active_exp.sum())==20145),
    ("sum n_approved_exp==282",int(T1.n_approved_exp.sum())==282),("enough_n==126",int(T1.enough_n_exp.sum())==126),
    ("point-enriched==50",int((nume>0).sum())==50),("ci_robust target==28",int(T1.ci_robust_exp_backend.sum())==28),
    ("ci_robust aiincl==31",int(T1.ci_robust_aiincl_backend.sum())==31),("0-approved enough_n==41",int(zero.sum())==41)]
allok=True
for n,ok in ck: print(f"  [{'PASS' if ok else 'FAIL'}] {n}"); allok&=ok
if not allok: print("*** MISMATCH — not writing ***"); raise SystemExit(1)
T1.to_csv(TA/"target_clinical_enrichment.tsv",sep="\t",index=False)
print(f"\nWROTE target_clinical_enrichment.tsv ({len(T1)} rows, {T1.shape[1]} cols)")
print("DONE.")
