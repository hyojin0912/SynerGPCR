"""
webexport_step2_core.py — WebExport_v2 STEP 2. READ-ONLY on sources.

Generates FOUR core export files, all stats RECALCULATED from PROVENANCE:
  1. compound_lookup.csv          (PK InChIKey)
  2. target_lookup.csv            (PK UniProt_AC)
  3. compound_chain_summary.csv   (PK InChIKey,UniProt_AC)
  4. gpcract_predictions.csv      (PK InChIKey,UniProt_AC) — replaces placeholder

Builds everything in memory, runs ALL self-checks, and writes ONLY if every
check passes. If any check fails -> raises, writes nothing.

Reference strings (preferred_name, smiles, clinical_phase) are carried from the
existing curated lookup (Output/NAR/WebDB/compound_lookup_final.csv); all COUNTS
/ stats (is_approved, in_assay_universe, GPCR universe, pair universe, AI cells)
are recomputed from the PROVENANCE sources below.
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import sys
import pandas as pd

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


B = BASE
OUT = B / "Output/NAR_DBI/WebExport_v2"

# ---- PROVENANCE sources ----
CATALOG = B / "Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv"
MASTER  = B / "Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv"
PAIRS   = B / "Output/NAR_DBI/Tables/SynerGPCR_WebTable_pairs_v2.csv"
OLDLOOK = B / "Output/NAR/WebDB/compound_lookup_final.csv"   # reference strings only
CHEMBL  = B / "Input/ChEMBL_GPCR_Info.csv"
PDBINFO = B / "Input/Human_GPCR_PDB_Info.csv"

DEF_LAYERS = {"Layer 1 (Binding)", "Layer 2 (Proximal)",
              "Layer 3 (Biased)", "Layer 4 (Reporter)"}

def load(p, **kw):
    if not p.exists():
        sys.exit(f"STOP: missing PROVENANCE file {p}")
    df = pd.read_csv(p, **kw)
    print(f"  loaded {p.name}: {len(df):,} rows")
    return df

print("="*78); print("LOADING PROVENANCE"); print("="*78)
cat = load(CATALOG, low_memory=False)
m   = load(MASTER, low_memory=False,
           usecols=["Ligand_InChIKey", "GPCR_UniProt", "Assay_Layer"])
pairs = load(PAIRS, comment="#", low_memory=False)
oldlk = load(OLDLOOK, low_memory=False)

# ---- derive experimental universes (defined layers only) ----
md = m[m.Assay_Layer.isin(DEF_LAYERS)].copy()
assay_compounds = set(md.Ligand_InChIKey.unique())          # expect 277,544
gpcr_universe   = sorted(md.GPCR_UniProt.unique())          # expect 600
exp_pairs = md[["Ligand_InChIKey", "GPCR_UniProt"]].drop_duplicates() \
              .rename(columns={"Ligand_InChIKey": "InChIKey",
                               "GPCR_UniProt": "UniProt_AC"})
print(f"\nassay-bearing compounds (>=1 defined layer): {len(assay_compounds):,}")
print(f"experimentally-covered GPCRs:                {len(gpcr_universe):,}")
print(f"defined-layer experimental pairs:            {len(exp_pairs):,}")

# --- 600-GPCR web-DB scope: exclude prediction GPCRs with no defined-layer experimental data ---
gpcr_set = set(gpcr_universe)                       # the 600 experimentally-grounded GPCRs
excluded = sorted(set(pairs.UniProt_AC.unique()) - gpcr_set)   # 14 Unknown-only prediction GPCRs
pairs = pairs[pairs.UniProt_AC.isin(gpcr_set)].copy()          # DROP the 14 from all pair files
print(f"\nEXCLUDED {len(excluded)} Unknown-only prediction GPCRs from web DB: {excluded}")
print(f"pairs_v2 after 600-scope filter: {len(pairs):,} rows (from 399,770)")
pred_pairs = set(map(tuple, pairs[["InChIKey", "UniProt_AC"]].values))
pred_gpcrs = set(pairs.UniProt_AC.unique())

# ==========================================================================
# FILE 1 — compound_lookup.csv
# ==========================================================================
print("\n" + "="*78); print("FILE 1: compound_lookup.csv"); print("="*78)
cl = cat[["InChIKey", "is_approved"]].copy()           # is_approved authoritative (bool)
ref = oldlk[["InChIKey", "preferred_name", "smiles", "clinical_phase"]]
cl = cl.merge(ref, on="InChIKey", how="left")
cl["in_assay_universe"] = cl.InChIKey.isin(assay_compounds)
cl = cl[["InChIKey", "preferred_name", "smiles", "clinical_phase",
         "is_approved", "in_assay_universe"]]
print(f"  rows={len(cl):,}  is_approved dtype={cl.is_approved.dtype}  "
      f"approved={int(cl.is_approved.sum()):,}  in_assay_universe={int(cl.in_assay_universe.sum()):,}")

# ==========================================================================
# FILE 2 — target_lookup.csv  (universe = 600 experimentally-covered GPCRs)
# ==========================================================================
print("\n" + "="*78); print("FILE 2: target_lookup.csv"); print("="*78)
tl = pd.DataFrame({"UniProt_AC": gpcr_universe})

ch = pd.read_csv(CHEMBL, low_memory=False)
ch["UniProt_AC"] = ch["UniProt Accessions"].astype(str).str.split(";").str[0].str.strip()
ch = ch.rename(columns={"Name": "protein_name", "Class": "gpcr_class_letter",
                        "Target Class": "gpcr_class_name", "Receptor Family": "gpcr_family",
                        "Organism": "organism", "ChEMBL ID": "chembl_target_id"})
ch = ch[["UniProt_AC", "protein_name", "gpcr_class_letter", "gpcr_class_name",
         "gpcr_family", "organism", "chembl_target_id"]].drop_duplicates("UniProt_AC")
pdb = pd.read_csv(PDBINFO, usecols=["Entry", "Entry Name", "Gene Names (primary)"], low_memory=False)
pdb = pdb.rename(columns={"Entry": "UniProt_AC", "Entry Name": "uniprot_entry_name",
                          "Gene Names (primary)": "gene_name"}).drop_duplicates("UniProt_AC")

tl = tl.merge(pdb[["UniProt_AC", "gene_name", "uniprot_entry_name"]], on="UniProt_AC", how="left")
tl = tl.merge(ch, on="UniProt_AC", how="left")

npairs = exp_pairs.groupby("UniProt_AC").size().rename("n_pairs_in_db")
tl = tl.merge(npairs, on="UniProt_AC", how="left")
tl["n_pairs_in_db"] = tl["n_pairs_in_db"].fillna(0).astype(int)
tl["has_gpcract_prediction"] = tl.UniProt_AC.isin(pred_gpcrs)

tl = tl[["UniProt_AC", "gene_name", "uniprot_entry_name", "protein_name",
         "gpcr_class_letter", "gpcr_class_name", "gpcr_family", "organism",
         "chembl_target_id", "n_pairs_in_db", "has_gpcract_prediction"]]
print(f"  rows(GPCRs)={len(tl):,}  has_gpcract_prediction True={int(tl.has_gpcract_prediction.sum())}  "
      f"with gene_name={int(tl.gene_name.notna().sum())}")

# ==========================================================================
# FILE 3 — compound_chain_summary.csv  (experimental backbone, ALL 600 GPCRs)
# ==========================================================================
print("\n" + "="*78); print("FILE 3: compound_chain_summary.csv"); print("="*78)
chain = cat[["InChIKey", "L1_state", "L2_state", "L3_state", "L4_state",
             "chain_pattern", "is_approved"]]
ccs = exp_pairs.merge(chain, on="InChIKey", how="left")     # experimental states (compound-level), all pairs
ccs["gpcract_prediction_available"] = [
    (ik, up) in pred_pairs for ik, up in zip(ccs.InChIKey, ccs.UniProt_AC)
]
ccs = ccs[["InChIKey", "UniProt_AC", "L1_state", "L2_state", "L3_state", "L4_state",
           "chain_pattern", "is_approved", "gpcract_prediction_available"]]
# reconcile anchors at compound level
cdedup = ccs.drop_duplicates("InChIKey")
n_full = int((cdedup.chain_pattern == "Full_chain").sum())
A = lambda c: cdedup[c] == "Active"
n_l1l2l4 = int((A("L1_state") & A("L2_state") & A("L4_state")).sum())
print(f"  rows(pairs)={len(ccs):,}  unique compounds={cdedup.InChIKey.nunique():,}")
print(f"  anchors: Full_chain compounds={n_full}  L1&L2&L4-active compounds={n_l1l2l4}  "
      f"gpcract_avail pairs={int(ccs.gpcract_prediction_available.sum()):,}")

# ==========================================================================
# FILE 4 — gpcract_predictions.csv  (from pairs_v2, AS-IS at locked op-points)
# ==========================================================================
print("\n" + "="*78); print("FILE 4: gpcract_predictions.csv"); print("="*78)
gp = pairs[["InChIKey", "UniProt_AC", "l1_evidence",
            "L2_ai_active", "L2_predicted_moa", "L2_confidence",
            "L3_ai_active", "L3_predicted_moa", "L3_confidence",
            "protein_graph_type"]].rename(columns={"l1_evidence": "L1_evidence"})
l2cells = int((gp.L2_ai_active == True).sum())
l3cells = int((gp.L3_ai_active == True).sum())
aiact = gp[(gp.L2_ai_active == True) | (gp.L3_ai_active == True)]
appr_set = set(pairs.loc[pairs.is_approved == True, "InChIKey"])
n_appr_ai = aiact[aiact.InChIKey.isin(appr_set)].InChIKey.nunique()
print(f"  rows={len(gp):,}  L2 AI cells={l2cells:,}  L3 AI cells={l3cells:,}  total={l2cells+l3cells:,}")
print(f"  unique compounds w/ >=1 AI-active={aiact.InChIKey.nunique():,}  approved among them={n_appr_ai}")

# ==========================================================================
# SELF-CHECKS (acceptance) — write nothing unless ALL pass
# ==========================================================================
print("\n" + "="*78); print("SELF-CHECKS vs CLAUDE.md §3 / §5"); print("="*78)
checks = []
def ck(name, cond, got, exp):
    checks.append((name, bool(cond), got, exp))

ck("compound_lookup PK == 481,240",            len(cl) == 481240, len(cl), 481240)
ck("compound_lookup is_approved dtype == bool", cl.is_approved.dtype == bool, str(cl.is_approved.dtype), "bool")
ck("compound_lookup approved == 2,420",        int(cl.is_approved.sum()) == 2420, int(cl.is_approved.sum()), 2420)
ck("compound_lookup in_assay_universe == 277,544", int(cl.in_assay_universe.sum()) == 277544, int(cl.in_assay_universe.sum()), 277544)
# --- target_lookup is scoped to the 600 defined-layer GPCRs (14 Unknown-only dropped) ---
ck("target GPCR universe == 600",              len(tl) == 600, len(tl), 600)
ck("target GPCR universe != stale 228",        len(tl) != 228, len(tl), "!=228")
ck("target has_gpcract_prediction == 157 (600-scope)", int(tl.has_gpcract_prediction.sum()) == 157, int(tl.has_gpcract_prediction.sum()), 157)
ck("chain_summary rows != stale 260,938",      len(ccs) != 260938, len(ccs), "!=260938")
ck("chain_summary anchor Full_chain == 228",   n_full == 228, n_full, 228)
ck("chain_summary anchor L1&L2&L4 == 2,580",   n_l1l2l4 == 2580, n_l1l2l4, 2580)
# --- AI-Active stats at the NEW 600-GPCR web-DB scope (supersede 171-GPCR figures for web/DB) ---
ck("gpcract rows == 399,520 (600-scope)",      len(gp) == 399520, len(gp), 399520)
ck("gpcract GPCRs == 157 (600-scope)",         gp.UniProt_AC.nunique() == 157, gp.UniProt_AC.nunique(), 157)
ck("gpcract L2 AI cells == 225,114 (600-scope)", l2cells == 225114, l2cells, 225114)
ck("gpcract L3 AI cells == 147,414 (600-scope)", l3cells == 147414, l3cells, 147414)
ck("gpcract total AI cells == 372,528 (600-scope)", l2cells + l3cells == 372528, l2cells + l3cells, 372528)
ck("gpcract compounds w/ AI-active == 171,821 (600-scope)", aiact.InChIKey.nunique() == 171821, aiact.InChIKey.nunique(), 171821)
ck("gpcract approved among AI-active == 842 (600-scope)", n_appr_ai == 842, n_appr_ai, 842)

# §5 STALE row-count guard already covered by 228/260,938 checks above.
allpass = all(c[1] for c in checks)
for name, ok, got, exp in checks:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:48s} got={got} exp={exp}")

if not allpass:
    print("\n*** SELF-CHECK FAILED — writing NOTHING. ***")
    sys.exit(1)

# ==========================================================================
# WRITE
# ==========================================================================
print("\n" + "="*78); print("ALL CHECKS PASS — WRITING"); print("="*78)
OUT.mkdir(parents=True, exist_ok=True)
for df, fn in [(cl, "compound_lookup.csv"), (tl, "target_lookup.csv"),
               (ccs, "compound_chain_summary.csv"), (gp, "gpcract_predictions.csv")]:
    path = OUT / fn
    df.to_csv(path, index=False)
    print(f"  wrote {path}  ({len(df):,} rows, {df.shape[1]} cols)")
print("\nDONE.")
