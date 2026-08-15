"""
build_catalog_v3.py
====================
Builds Compound_Chain_Catalog_v3.csv by extending v2 with two new columns:

  max_phase              : int   (0–4) — highest clinical phase reached by the compound
  moa_concordance_status : str   ('concordant' | 'discordant' | 'NA')
                           — agreement between compound's bioassay MoA (L2/L4)
                             and its human clinical reference MoA

Sources
-------
  v2 catalog            : Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v2.csv
  Drug indications      : Output/DB/GPCRactDB/Drug_Indication_Master_Integrated.csv
  Human (clinical) MoA  : Output/DB/GPCRactDB/Human_MoA_Master_Integrated.csv
  Bioassay MoA          : Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv

Output
------
  Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.

from __future__ import annotations
import os

import sys
from pathlib import Path
from collections import Counter
from typing import Optional

import pandas as pd

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = BASE

V2_CATALOG    = ROOT / "Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v2.csv"
DRUG_IND      = ROOT / "Output/DB/GPCRactDB/Drug_Indication_Master_Integrated.csv"
HUMAN_MOA     = ROOT / "Output/DB/GPCRactDB/Human_MoA_Master_Integrated.csv"
BIOASSAY_MOA  = ROOT / "Output/DB/GPCRactDB/Bioassay_MoA_Master_Integrated.csv"
OUT_V3        = ROOT / "Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
# Map Highest_Status → numeric phase
PHASE_MAP: dict[str, int] = {
    "Approved":        4,
    "Phase 3":         3,
    "Phase 2":         2,
    "Phase 1":         1,
    "Investigational": 0,
}

# MoA direction families (case-sensitive to match source files)
MOA_AGONIST_SET: frozenset[str] = frozenset(
    {"Agonist", "Partial Agonist", "PAM"}
)
MOA_ANTAGONIST_SET: frozenset[str] = frozenset(
    {"Antagonist", "Inverse Agonist", "NAM"}
)

# Bioassay layers to use for MoA direction inference (in priority order)
BIOASSAY_LAYERS_ORDERED = [
    "Layer 2 (Proximal)",
    "Layer 4 (Reporter)",
]

# MoA values that carry no directional information
MOA_UNINFORMATIVE: frozenset[str] = frozenset(
    {"Inactive", "Binder", "Non-binder", "Unknown"}
)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def moa_direction(moa: str) -> Optional[str]:
    """Return canonical direction ('agonist' | 'antagonist' | None)."""
    if moa in MOA_AGONIST_SET:
        return "agonist"
    if moa in MOA_ANTAGONIST_SET:
        return "antagonist"
    return None


def dominant_direction(moa_series: pd.Series) -> Optional[str]:
    """Return the most frequent directional MoA in *moa_series*, or None.

    Parameters
    ----------
    moa_series:
        Series of raw MoA label strings (may contain uninformative values).

    Returns
    -------
    'agonist', 'antagonist', or None when no directional call can be made.
    """
    directions = [
        moa_direction(m)
        for m in moa_series.dropna()
        if moa_direction(m) is not None
    ]
    if not directions:
        return None
    counts = Counter(directions)
    top, top_n = counts.most_common(1)[0]
    # Require clear majority (> 50 %) to avoid ties counting as signal
    if top_n > sum(counts.values()) / 2:
        return top
    return None


# ---------------------------------------------------------------------------
# Step 0 — Column validation
# ---------------------------------------------------------------------------

def validate_columns() -> None:
    """Smoke-test all four input files for required columns."""
    requirements = {
        V2_CATALOG:   ["InChIKey", "L1_state", "L2_state", "L3_state",
                        "L4_state", "is_approved"],
        DRUG_IND:     ["InChIKey", "Highest_Status"],
        HUMAN_MOA:    ["Ligand_InChIKey", "Human_MoA", "Final_Tier"],
        BIOASSAY_MOA: ["Ligand_InChIKey", "Assay_Layer", "MoA"],
    }
    for path, required in requirements.items():
        cols = pd.read_csv(path, nrows=1).columns.tolist()
        missing = [c for c in required if c not in cols]
        if missing:
            print(f"[ERROR] {path.name} is missing columns: {missing}")
            print(f"        Found: {cols}")
            sys.exit(1)
        print(f"[OK] {path.name} — all required columns present")
    print("=== STEP 0 COMPLETE ===\n")


# ---------------------------------------------------------------------------
# Step 1 — Build max_phase lookup
# ---------------------------------------------------------------------------

def build_max_phase(drug_ind: pd.DataFrame) -> pd.Series:
    """Derive per-compound max_phase (0–4) from Drug_Indication table.

    Parameters
    ----------
    drug_ind:
        Drug_Indication_Master_Integrated DataFrame with columns
        ['InChIKey', 'Highest_Status'].

    Returns
    -------
    pd.Series indexed by InChIKey with integer max_phase values (0–4).
    Missing or unrecognised statuses map to 0.
    """
    unknown_statuses = set(drug_ind["Highest_Status"].dropna().unique()) - set(PHASE_MAP)
    if unknown_statuses:
        print(f"  [WARN] Unrecognised Highest_Status values (mapped to 0): "
              f"{sorted(unknown_statuses)}")

    drug_ind = drug_ind.copy()
    drug_ind["phase_num"] = (
        drug_ind["Highest_Status"].map(PHASE_MAP).fillna(0).astype(int)
    )
    max_phase = drug_ind.groupby("InChIKey")["phase_num"].max()
    print(f"  max_phase distribution:\n"
          f"{max_phase.value_counts().sort_index().to_string()}")
    return max_phase


# ---------------------------------------------------------------------------
# Step 2 — Build human MoA direction lookup
# ---------------------------------------------------------------------------

def build_human_moa_direction(human_moa: pd.DataFrame) -> pd.Series:
    """Return per-compound dominant human-clinical MoA direction.

    Tier 1 rows take priority over Tier 1.5 — if a compound has any
    Tier 1 entries, only those are used for the direction call.

    Parameters
    ----------
    human_moa:
        Human_MoA_Master_Integrated DataFrame.

    Returns
    -------
    pd.Series indexed by Ligand_InChIKey with values 'agonist' | 'antagonist'.
    Compounds with no clear directional call are absent from the index.
    """
    results: dict[str, str] = {}

    grouped = human_moa.groupby("Ligand_InChIKey")
    for key, grp in grouped:
        tier1 = grp[grp["Final_Tier"] == "Tier 1"]
        subset = tier1 if not tier1.empty else grp
        direction = dominant_direction(subset["Human_MoA"])
        if direction is not None:
            results[key] = direction

    series = pd.Series(results, name="human_direction")
    print(f"  Human MoA direction counts:\n"
          f"{series.value_counts().to_string()}")
    return series


# ---------------------------------------------------------------------------
# Step 3 — Build bioassay MoA direction lookup
# ---------------------------------------------------------------------------

def build_bioassay_moa_direction(bioassay: pd.DataFrame) -> pd.Series:
    """Return per-compound dominant bioassay MoA direction (L2 preferred, L4 fallback).

    Strategy
    --------
    1. Consider only directional (non-uninformative) MoA entries.
    2. For each compound, prefer Layer 2 (Proximal) entries; fall back to
       Layer 4 (Reporter) when L2 yields no directional call.

    Parameters
    ----------
    bioassay:
        Bioassay_MoA_Master_Integrated DataFrame.

    Returns
    -------
    pd.Series indexed by Ligand_InChIKey with values 'agonist' | 'antagonist'.
    Compounds with no directional bioassay call are absent.
    """
    # Filter to only the two layers of interest
    relevant = bioassay[
        bioassay["Assay_Layer"].isin(BIOASSAY_LAYERS_ORDERED)
    ].copy()

    results: dict[str, str] = {}
    grouped = relevant.groupby("Ligand_InChIKey")

    for key, grp in grouped:
        direction: Optional[str] = None
        for layer in BIOASSAY_LAYERS_ORDERED:
            layer_grp = grp[grp["Assay_Layer"] == layer]
            direction = dominant_direction(layer_grp["MoA"])
            if direction is not None:
                break
        if direction is not None:
            results[key] = direction

    series = pd.Series(results, name="bioassay_direction")
    print(f"  Bioassay MoA direction counts:\n"
          f"{series.value_counts().to_string()}")
    return series


# ---------------------------------------------------------------------------
# Step 4 — Derive moa_concordance_status
# ---------------------------------------------------------------------------

def derive_concordance(
    catalog: pd.DataFrame,
    human_direction: pd.Series,
    bioassay_direction: pd.Series,
) -> pd.Series:
    """Assign 'concordant' | 'discordant' | 'NA' per compound.

    A compound is:
      - 'concordant'  : bioassay direction == human clinical direction
      - 'discordant'  : bioassay direction != human clinical direction
      - 'NA'          : either direction is unavailable

    Parameters
    ----------
    catalog:
        v2 catalog DataFrame (index: InChIKey).
    human_direction:
        Series from build_human_moa_direction().
    bioassay_direction:
        Series from build_bioassay_moa_direction().

    Returns
    -------
    pd.Series of concordance status, aligned to catalog.index.
    """
    keys = catalog["InChIKey"]
    hd = human_direction.reindex(keys).values
    bd = bioassay_direction.reindex(keys).values

    concordance: list[str] = []
    for h, b in zip(hd, bd):
        if pd.isna(h) or pd.isna(b):
            concordance.append("NA")
        elif h == b:
            concordance.append("concordant")
        else:
            concordance.append("discordant")

    series = pd.Series(concordance, index=catalog.index, name="moa_concordance_status")
    print(f"  moa_concordance_status distribution:\n"
          f"{series.value_counts().to_string()}")
    return series


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # ---- Step 0: validate inputs ----------------------------------------
    validate_columns()

    # ---- Load data -------------------------------------------------------
    print("Loading input files …")
    catalog   = pd.read_csv(V2_CATALOG)
    drug_ind  = pd.read_csv(DRUG_IND)
    human_moa = pd.read_csv(HUMAN_MOA)
    bioassay  = pd.read_csv(BIOASSAY_MOA)
    print(f"  v2 catalog  : {catalog.shape[0]:,} rows")
    print(f"  drug_ind    : {drug_ind.shape[0]:,} rows")
    print(f"  human_moa   : {human_moa.shape[0]:,} rows")
    print(f"  bioassay    : {bioassay.shape[0]:,} rows")

    # ---- Step 1: max_phase -----------------------------------------------
    print("\n[Step 1] Building max_phase …")
    max_phase_lookup = build_max_phase(drug_ind)
    catalog["max_phase"] = (
        catalog["InChIKey"]
        .map(max_phase_lookup)
        .fillna(0)
        .astype(int)
    )
    n_with_phase = (catalog["max_phase"] > 0).sum()
    print(f"  Compounds with max_phase > 0 : {n_with_phase:,} / {len(catalog):,}")
    print("=== STEP 1 COMPLETE ===\n")

    # ---- Step 2: human MoA direction -------------------------------------
    print("[Step 2] Deriving human MoA direction …")
    human_direction = build_human_moa_direction(human_moa)
    print("=== STEP 2 COMPLETE ===\n")

    # ---- Step 3: bioassay MoA direction ----------------------------------
    print("[Step 3] Deriving bioassay MoA direction …")
    bioassay_direction = build_bioassay_moa_direction(bioassay)
    print("=== STEP 3 COMPLETE ===\n")

    # ---- Step 4: concordance status + human_dir column ------------------
    print("[Step 4] Computing moa_concordance_status and human_dir …")
    catalog["moa_concordance_status"] = derive_concordance(
        catalog, human_direction, bioassay_direction
    )
    # Expose the resolved human clinical direction so harness_core.py can
    # compare it directly against the model's predicted MoA direction.
    catalog["human_dir"] = catalog["InChIKey"].map(human_direction)
    print(f"  human_dir non-null: {catalog['human_dir'].notna().sum():,}")
    print("=== STEP 4 COMPLETE ===\n")

    # ---- Sanity checks ---------------------------------------------------
    print("[Sanity] Running validation checks …")

    # All is_approved compounds must have max_phase == 4
    approved = catalog[catalog["is_approved"]]
    bad_phase = (approved["max_phase"] != 4).sum()
    print(f"  is_approved & max_phase != 4 : {bad_phase} "
          f"(expected 0 or small — some approved drugs may lack indication data)")

    # Distribution of concordance among approved compounds
    print("  moa_concordance_status among approved drugs:")
    print(approved["moa_concordance_status"].value_counts().to_string())

    # v3 must have exactly 3 more columns than v2
    expected_new_cols = {"max_phase", "moa_concordance_status", "human_dir"}
    v2_cols = set(pd.read_csv(V2_CATALOG, nrows=0).columns)
    v3_cols = set(catalog.columns)
    added = v3_cols - v2_cols
    assert added == expected_new_cols, (
        f"Unexpected column diff: added={added}, expected={expected_new_cols}"
    )
    print(f"  New columns added (expected): {sorted(added)}")
    print("=== SANITY CHECKS PASSED ===\n")

    # ---- Save ------------------------------------------------------------
    OUT_V3.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(OUT_V3, index=False)
    print(f"[DONE] Saved {catalog.shape[0]:,} rows to:\n  {OUT_V3}")


if __name__ == "__main__":
    main()
