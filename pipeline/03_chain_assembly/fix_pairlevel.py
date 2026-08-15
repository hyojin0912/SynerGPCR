"""
fix_compound_chain_summary_pairlevel.py — WebExport_v2 hotfix. READ-ONLY on sources.

PROBLEM (reported from 229 web server):
  compound_chain_summary.csv had L1/L2/L3/L4_state filled by COMPOUND-LEVEL
  aggregation. Every (InChIKey x UniProt_AC) row of the same compound carried the
  SAME L*_state, so e.g. Aripiprazole-DRD2 inherited L* assay outcomes Aripiprazole
  had on *other* receptors. Root cause: webexport_step2_core.py line ~120 merged the
  compound-keyed Compound_Chain_Catalog_v3 states onto every pair.

FIX:
  Derive L*_state per (InChIKey x UniProt_AC) DIRECTLY from the raw assay master,
  using the canonical MoA->state mapping + "Active dominates" aggregation taken
  VERBATIM from the validated PROVENANCE script
  Output/NAR_DBI/Code/part11_final_op_0p40_0p97.py (lines 39-41, 145-153).
  chain_pattern is recomputed PAIR-LEVEL for self-consistency (display-only column;
  synergy_stats.csv derives chain_pattern from Chain_Pattern_Distribution.csv, not here).
  is_approved stays compound-level (approval is a compound property).
  gpcract_prediction_available stays pair-level (from pairs_v2), unchanged.

Self-checks must ALL pass or it writes nothing.
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import sys
import numpy as np
import pandas as pd

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


B   = BASE
OUT = B / "Output/NAR_DBI/WebExport_v2"

MASTER  = B / "Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv"
CATALOG = B / "Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv"
PAIRS   = B / "Output/NAR_DBI/Tables/SynerGPCR_WebTable_pairs_v2.csv"

# ---- canonical MoA->state (verbatim from part11_final_op_0p40_0p97.py) ----
ACTIVE_MOA   = {"Binder", "Agonist", "Antagonist", "Inverse Agonist",
                "Partial Agonist", "PAM", "NAM"}
INACTIVE_MOA = {"Inactive", "Non-binder"}
LAYER_MAP    = {"Layer 1 (Binding)": "L1", "Layer 2 (Proximal)": "L2",
                "Layer 3 (Biased)": "L3", "Layer 4 (Reporter)": "L4"}


def load(p, **kw):
    if not p.exists():
        sys.exit(f"STOP: missing PROVENANCE file {p}")
    df = pd.read_csv(p, **kw)
    print(f"  loaded {p.name}: {len(df):,} rows")
    return df


print("=" * 78); print("LOADING PROVENANCE"); print("=" * 78)
m = load(MASTER, low_memory=False,
         usecols=["Ligand_InChIKey", "GPCR_UniProt", "Assay_Layer", "MoA"]).rename(
         columns={"Ligand_InChIKey": "InChIKey", "GPCR_UniProt": "UniProt_AC"})
cat = load(CATALOG, low_memory=False).drop_duplicates("InChIKey")
pairs = load(PAIRS, comment="#", low_memory=False)

# ==========================================================================
# PAIR-LEVEL L*_state from raw assays (verbatim agg logic)
# ==========================================================================
print("\n" + "=" * 78); print("DERIVING PAIR-LEVEL L*_state"); print("=" * 78)
mw = m[m.Assay_Layer.isin(LAYER_MAP)].copy()
mw["layer"] = mw.Assay_Layer.map(LAYER_MAP)
mw["st"] = np.where(mw.MoA.isin(ACTIVE_MOA), "Active",
                    np.where(mw.MoA.isin(INACTIVE_MOA), "Inactive", "Unknown"))
agg = lambda s: "Active" if (s == "Active").any() else ("Inactive" if (s == "Inactive").any() else "NA")
grp = mw.groupby(["InChIKey", "UniProt_AC", "layer"]).st.agg(agg).unstack("layer")
for c in ["L1", "L2", "L3", "L4"]:
    if c not in grp:
        grp[c] = "NA"
grp = grp[["L1", "L2", "L3", "L4"]].fillna("NA").reset_index()
grp.columns = ["InChIKey", "UniProt_AC", "L1_state", "L2_state", "L3_state", "L4_state"]
# normalise the "NA" sentinel to true missing (matches old file's empty cells)
for c in ["L1_state", "L2_state", "L3_state", "L4_state"]:
    grp.loc[grp[c] == "NA", c] = np.nan
print(f"  pair rows={len(grp):,}  compounds={grp.InChIKey.nunique():,}  GPCRs={grp.UniProt_AC.nunique():,}")

# ---- pair-level chain_pattern (which layers are Active for THIS pair) ----
def pattern(row):
    on = [L for L in ("L1", "L2", "L3", "L4") if row[f"{L}_state"] == "Active"]
    if not on:
        return "No_active_layer"
    if len(on) == 4:
        return "Full_chain"
    return "+".join(on) if len(on) > 1 else f"{on[0]}_only"
grp["chain_pattern"] = grp.apply(pattern, axis=1)

# ---- carry is_approved (compound-level) and gpcract_prediction_available (pair-level) ----
appr = cat.set_index("InChIKey")["is_approved"].to_dict()
grp["is_approved"] = grp.InChIKey.map(appr).fillna(False)

gpcr_set = set(grp.UniProt_AC.unique())               # 600 experimentally-grounded GPCRs
pairs600 = pairs[pairs.UniProt_AC.isin(gpcr_set)]
pred_pairs = set(map(tuple, pairs600[["InChIKey", "UniProt_AC"]].values))
grp["gpcract_prediction_available"] = [
    (ik, up) in pred_pairs for ik, up in zip(grp.InChIKey, grp.UniProt_AC)
]

ccs = grp[["InChIKey", "UniProt_AC", "L1_state", "L2_state", "L3_state", "L4_state",
           "chain_pattern", "is_approved", "gpcract_prediction_available"]]

# ==========================================================================
# SELF-CHECKS
# ==========================================================================
print("\n" + "=" * 78); print("SELF-CHECKS"); print("=" * 78)
checks = []
def ck(name, cond, got, exp):
    checks.append((name, bool(cond), got, exp))

ck("pair rows == 396,169",          len(ccs) == 396169, len(ccs), 396169)
ck("unique compounds == 277,544",   ccs.InChIKey.nunique() == 277544, ccs.InChIKey.nunique(), 277544)
ck("unique GPCRs == 600",           ccs.UniProt_AC.nunique() == 600, ccs.UniProt_AC.nunique(), 600)

# bug-fix proof: states must now vary across receptors for the same compound
varies = {}
for c in ["L1_state", "L2_state", "L3_state", "L4_state"]:
    u = ccs.groupby("InChIKey")[c].nunique(dropna=False)
    varies[c] = int((u >= 2).sum())
ck("L1_state varies across receptors (>0)", varies["L1_state"] > 0, varies["L1_state"], ">0")
ck("L2_state varies across receptors (>0)", varies["L2_state"] > 0, varies["L2_state"], ">0")
ck("L3_state varies across receptors (>0)", varies["L3_state"] > 0, varies["L3_state"], ">0")
ck("L4_state varies across receptors (>0)", varies["L4_state"] > 0, varies["L4_state"], ">0")

# acceptance anchors: pair->compound aggregation must reproduce Catalog_v3 anchors
comp = ccs.groupby("InChIKey")[["L1_state", "L2_state", "L3_state", "L4_state"]].agg(
    lambda s: "Active" if (s == "Active").any() else ("Inactive" if (s == "Inactive").any() else "NA"))
A = lambda c: comp[c] == "Active"
n_full = int((A("L1_state") & A("L2_state") & A("L3_state") & A("L4_state")).sum())
n_l1l2l4 = int((A("L1_state") & A("L2_state") & A("L4_state")).sum())
ck("pair->compound Full_chain == 228",  n_full == 228, n_full, 228)
ck("pair->compound L1&L2&L4 == 2,580",  n_l1l2l4 == 2580, n_l1l2l4, 2580)

# state vocabulary sanity (only Active / Inactive / NaN)
allowed = {"Active", "Inactive"}
for c in ["L1_state", "L2_state", "L3_state", "L4_state"]:
    vals = set(ccs[c].dropna().unique())
    ck(f"{c} vocabulary subset of {{Active,Inactive}}", vals <= allowed, sorted(vals), "Active/Inactive")

allpass = all(c[1] for c in checks)
for name, ok, got, exp in checks:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:46s} got={got} exp={exp}")

if not allpass:
    print("\n*** SELF-CHECK FAILED — writing NOTHING. ***")
    sys.exit(1)

# ==========================================================================
# WRITE  (backup the buggy file first, then overwrite)
# ==========================================================================
print("\n" + "=" * 78); print("ALL CHECKS PASS — WRITING"); print("=" * 78)
target = OUT / "compound_chain_summary.csv"
backup = OUT / "compound_chain_summary.compound_level_BUGGY.csv"
if target.exists() and not backup.exists():
    target.rename(backup)
    print(f"  backed up buggy file -> {backup.name}")
ccs.to_csv(target, index=False)
print(f"  wrote {target}  ({len(ccs):,} rows, {ccs.shape[1]} cols)")
print("\nDONE.")
