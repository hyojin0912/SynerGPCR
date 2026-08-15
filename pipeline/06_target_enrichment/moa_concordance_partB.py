"""STEP 2 Part B — MoA-direction concordance TREND + coverage (Task 2 preview). READ-ONLY.
Neutral labels only; Part A 'co_active' concept NOT mixed in here."""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import numpy as np, pandas as pd
from statsmodels.stats.proportion import proportion_confint

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))

np.random.seed(0)
B=BASE; DB=B/"Output/DB/GPCRactDB/"
print("pandas",pd.__version__)
# direction collapse
AGO={"Agonist","Partial Agonist","PAM"}; ANT={"Antagonist","Inverse Agonist","NAM"}
LMAP={"Layer 1 (Binding)":"L1","Layer 2 (Proximal)":"L2","Layer 3 (Biased)":"L3","Layer 4 (Reporter)":"L4"}
m=pd.read_csv(DB/"Bioassay_MoA_Master_Integrated.csv",usecols=["Ligand_InChIKey","GPCR_UniProt","Assay_Layer","MoA"]).rename(columns={"Ligand_InChIKey":"InChIKey","GPCR_UniProt":"UniProt_AC"})
md=m[m.Assay_Layer.isin(LMAP)].copy(); md["layer"]=md.Assay_Layer.map(LMAP)
md["dir"]=np.where(md.MoA.isin(AGO),"ago",np.where(md.MoA.isin(ANT),"ant",np.nan))
print("\nB2) direction class distribution per layer (directional records only):")
for L in ["L2","L3","L4"]:
    s=md[(md.layer==L)&md.dir.notna()]
    print(f"  {L}: ago={int((s.dir=='ago').sum()):,}  ant={int((s.dir=='ant').sum()):,}  (Inverse Agonist→ant-dir)")

# per (compound,GPCR,layer) direction: unambiguous if only one dir present
g=md[md.dir.notna()].groupby(["InChIKey","UniProt_AC","layer"]).dir.agg(lambda s:"ago" if set(s)=={"ago"} else ("ant" if set(s)=={"ant"} else "both"))
piv=g.unstack("layer")
cat=pd.read_csv(DB/"Analysis/Compound_Chain_Catalog_v3.csv",low_memory=False,usecols=["InChIKey","is_approved","human_dir"]).drop_duplicates("InChIKey")
appr=set(cat.loc[cat.is_approved==True,"InChIKey"])
def w(k,n): return proportion_confint(k,n,0.05,"wilson") if n>0 else (np.nan,np.nan)

print("\nB3/B4) pairwise both-direction coverage at same GPCR + concordant/discordant approval:")
def pair(Lx,Ly):
    if Lx not in piv or Ly not in piv: print(f"  {Lx}∩{Ly}: layer absent"); return
    d=piv[[Lx,Ly]].dropna()
    d=d[(d[Lx].isin(["ago","ant"]))&(d[Ly].isin(["ago","ant"]))]
    conc=d[d[Lx]==d[Ly]]; disc=d[d[Lx]!=d[Ly]]
    # compound-level approval (dedup; any pair in group)
    def rate(sub):
        ik=set(sub.reset_index().InChIKey); k=len(ik&appr); n=len(ik); lo,hi=w(k,n)
        return n,k,(k/n if n else np.nan),lo,hi
    nC,kC,rC,loC,hiC=rate(conc); nD,kD,rD,loD,hiD=rate(disc)
    print(f"  {Lx}∩{Ly}: both-dir (compound,GPCR) pairs={len(d):,}  concordant={len(conc):,}  discordant={len(disc):,}")
    print(f"      concordant compounds n={nC:,} approved={kC} rate={rC*100:.2f}% CI[{loC*100:.2f},{hiC*100:.2f}]" if nC else "      concordant: n=0")
    print(f"      discordant compounds n={nD:,} approved={kD} rate={rD*100:.2f}% CI[{loD*100:.2f},{hiD*100:.2f}]" if nD else "      discordant: n=0")
    if nC>=20 and nD>=20: print("      -> n sufficient (>=20 both): trend reportable")
    else: print("      -> n INSUFFICIENT (<20 in a group): report as limitation, NOT a conclusion")
pair("L2","L3"); pair("L2","L4"); pair("L3","L4")

print("\n  assay(L2)-vs-clinical(human_dir) coverage:")
# clinical dir from human_dir; assay L2 dir per compound (collapse over GPCRs: take if unambiguous)
l2c=md[(md.layer=="L2")&md.dir.notna()].groupby("InChIKey").dir.agg(lambda s:"ago" if set(s)=={"ago"} else ("ant" if set(s)=={"ant"} else "both"))
hd=cat.set_index("InChIKey").human_dir
clin_map={"Agonist":"ago","Partial Agonist":"ago","PAM":"ago","Antagonist":"ant","Inverse Agonist":"ant","NAM":"ant"}
hd2=hd.map(lambda v: clin_map.get(str(v),np.nan))
j=pd.DataFrame({"assay":l2c}).join(hd2.rename("clin")).dropna()
j=j[j.assay.isin(["ago","ant"])&j.clin.isin(["ago","ant"])]
conc=j[j.assay==j.clin]; disc=j[j.assay!=j.clin]
def rate2(ikset): k=len(ikset&appr); n=len(ikset); lo,hi=w(k,n); return n,k,(k/n if n else np.nan),lo,hi
print(f"    compounds with BOTH L2 assay-dir and clinical-dir = {len(j):,}  concordant={len(conc):,} discordant={len(disc):,}")
for nm,sub in [("concordant",conc),("discordant",disc)]:
    n,k,r,lo,hi=rate2(set(sub.index))
    print(f"      {nm}: n={n:,} approved={k} rate={r*100:.2f}% CI[{lo*100:.2f},{hi*100:.2f}]")
print("\nNOTE: L2-L3 directional discordance = biased signaling (real biology), NOT an annotation error.")
print("DONE.")
