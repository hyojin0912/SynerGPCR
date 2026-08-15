"""Table1: Database Summary Statistics — source DB breakdown.

Content: Source_Database | Data_Type | Assay_Layer | n_raw_records | n_unique_pairs | n_compounds | n_targets
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.
import os
from pathlib import Path
import numpy as np
import pandas as pd

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


OUT_CSV = str(BASE / "Output/NAR/Tables/Table1_db_summary.csv")

# ---------------------------------------------------------------------------
# Source 1: Bioassay (L1/L2/L4 cell assays) — from Bioassay_MoA_Master
# ---------------------------------------------------------------------------
print("Loading Bioassay_MoA_Master_Integrated...")
bio = pd.read_csv(
    str(BASE / "Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv"),
    low_memory=False
)

bioassay_layer_map = {
    "Layer 1 (Binding)":  "L1 Binding",
    "Layer 2 (Proximal)": "L2 G-protein",
    "Layer 4 (Reporter)": "L4 Reporter",
    # GLASS and IUPHAR are stored as 'Unknown' layer in the integrated file
    # but represent L1 binding assays
    "Unknown":            "L1 Binding",
}
# Override layer labels per source for known sources
source_layer_override = {
    "GLASS":  "L1 Binding (GPCR-ligand binding)",
    "IUPHAR": "L1 Binding / L2 G-protein",
}

rows = []
for src in ["ChEMBL", "BindingDB", "GLASS", "PubChem", "IUPHAR"]:
    sub = bio[bio["Source_DB"] == src]
    if src in source_layer_override:
        layer_str = source_layer_override[src]
    else:
        layers = sub["Assay_Layer"].map(bioassay_layer_map).dropna().unique()
        layer_str = " / ".join(sorted(set(layers)))
    # count unique (compound, target) pairs
    pairs = sub[["Ligand_InChIKey", "GPCR_UniProt"]].drop_duplicates()
    rows.append({
        "Source_Database": src,
        "Data_Type": "Cell/biochemical assay (LLM-parsed MoA)",
        "Assay_Layer": layer_str,
        "n_raw_records": len(sub),
        "n_unique_pairs": len(pairs),
        "n_compounds": sub["Ligand_InChIKey"].nunique(),
        "n_targets": sub["GPCR_UniProt"].nunique(),
    })
    print(f"  {src}: {len(sub):,} records, {len(pairs):,} pairs, "
          f"{sub['Ligand_InChIKey'].nunique():,} compounds, "
          f"{sub['GPCR_UniProt'].nunique()} targets")

# ---------------------------------------------------------------------------
# Source 2: Human MoA (L2 curated) — from Human_MoA_Master_Integrated
# ---------------------------------------------------------------------------
print("\nLoading Human_MoA_Master_Integrated...")
hmoa = pd.read_csv(
    str(BASE / "Output/DB/GPCRactDB/Human_MoA_Master_Integrated.csv"),
    low_memory=False
)

# Explode Source_DBs (comma-separated) to count per source
hmoa_exp = hmoa.copy()
hmoa_exp["Source_DBs"] = hmoa_exp["Source_DBs"].str.split(", ")
hmoa_exp = hmoa_exp.explode("Source_DBs")
hmoa_exp["Source_DBs"] = hmoa_exp["Source_DBs"].str.strip()

clinical_sources = ["DrugBank", "DrugCentral", "DGIdb", "GPCRdb", "IUPHAR", "TTD"]
for src in clinical_sources:
    sub = hmoa_exp[hmoa_exp["Source_DBs"] == src]
    if len(sub) == 0:
        continue
    pairs_u = sub[["Ligand_InChIKey", "GPCR_UniProt"]].drop_duplicates()
    rows.append({
        "Source_Database": src,
        "Data_Type": "Clinical MoA (curated human pharmacology)",
        "Assay_Layer": "L2 G-protein (human)",
        "n_raw_records": len(sub),
        "n_unique_pairs": len(pairs_u),
        "n_compounds": sub["Ligand_InChIKey"].nunique(),
        "n_targets": sub["GPCR_UniProt"].nunique(),
    })
    print(f"  {src}: {len(sub):,} records, {len(pairs_u):,} pairs")

# ---------------------------------------------------------------------------
# TOTAL row
# ---------------------------------------------------------------------------
# From master CSV (authoritative curated total)
master = pd.read_csv(
    str(BASE / "Output/DB/GPCRactDB/GPCRactDB_v2_master.csv"),
    low_memory=False
)
rows.append({
    "Source_Database": "TOTAL (GPCRactDB v2)",
    "Data_Type": "Curated multi-source integration",
    "Assay_Layer": "L1 + L2 + L3 + L4",
    "n_raw_records": sum(r["n_raw_records"] for r in rows),
    "n_unique_pairs": 260_938,
    "n_compounds": master["InChIKey"].nunique(),
    "n_targets": master["UniProt_AC"].nunique(),
})

df_out = pd.DataFrame(rows)
df_out.to_csv(OUT_CSV, index=False)
print("\n" + df_out.to_string())
print(f"\n=== TABLE 1 COMPLETE: {OUT_CSV} ===")
