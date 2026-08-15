"""analyze_point3_moa_aware.py

TASK 3 — MoA-aware SI (two-step, L1+L2+L4 context)

Step A: Experimental concordance groups within the observed L1+L2+L4 pool.
Step B: Denovo-predicted concordance on the combined pool (obs + predictions).

Both use E_BLISS(L1+L2+L4) derived from Chain_Pattern_Distribution.csv.

Output:
    Output/NAR/Tables/SI_MoA_aware.csv
"""
# Data root is configurable: export SYNERGPCR_BASE=/path/to/released/tables
# (defaults to ./data). All input paths below are resolved against it.

from pathlib import Path
import os
import warnings
import numpy as np
import pandas as pd
import scipy.stats
from collections import Counter
from typing import Optional
from statsmodels.stats.proportion import proportion_confint

BASE = Path(os.environ.get("SYNERGPCR_BASE", "./data"))


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR   = str(BASE)
TABLE_DIR  = f"{BASE_DIR}/Output/NAR/Tables"
os.makedirs(TABLE_DIR, exist_ok=True)

CHAIN_DIST_PATH  = (
    f"{BASE_DIR}/Output/DB/GPCRactDB/Analysis/Chain_Pattern_Distribution.csv"
)
PERM_NULL_PATH   = (
    f"{BASE_DIR}/Output/DB/GPCRactDB/Analysis/PreTrain/Permutation_Null_SI_v3.csv"
)
CATALOG_V3_PATH  = (
    f"{BASE_DIR}/Output/DB/GPCRactDB/Analysis/Compound_Chain_Catalog_v3.csv"
)
DENOVO_PRED_PATH = (
    f"{BASE_DIR}/GPCRactDB/results/denovo/"
    "denovo_predictions_v2_clinical_withoutInactives.csv"
)

print("=" * 70)
print("TASK 3: MoA-Aware SI — Experimental + Denovo")
print("=" * 70)

# ===========================================================================
# Step 0: Column verification
# ===========================================================================
print("\n--- Step 0: Column verification ---")

cat = pd.read_csv(CATALOG_V3_PATH, keep_default_na=False)
print(f"Compound_Chain_Catalog_v3 ALL columns:\n  {list(cat.columns)}")

# Dynamic obs pool — never hardcode
_obs = (
    cat[
        (cat["L1_state"] == "Active") &
        (cat["L2_state"] == "Active") &
        (cat["L4_state"] == "Active")
    ]
    .drop_duplicates(subset="InChIKey")
    .copy()
)
_obs["max_phase"] = (
    pd.to_numeric(_obs["max_phase"], errors="coerce")
    .fillna(0).astype(int)
)
N_OBS = len(_obs)
K_OBS = int((_obs["max_phase"] == 4).sum())
print(f"Obs pool (v3, dynamic): N_OBS={N_OBS}, K_OBS={K_OBS}")
assert 2500 <= N_OBS <= 2700, f"N_OBS={N_OBS} unexpected"
assert 150  <= K_OBS <= 200,  f"K_OBS={K_OBS} unexpected"

REQUIRED = {
    "InChIKey", "L1_state", "L2_state", "L4_state",
    "max_phase", "human_dir",
}
missing_required = REQUIRED - set(cat.columns)
if missing_required:
    raise ValueError(f"STOP — missing required columns: {missing_required}")

# L2_label check (task-specified column)
if "L2_label" not in cat.columns:
    print(
        "  NOTE: 'L2_label' not found in v3 catalog. "
        "Using 'moa_concordance_status' (pre-computed 'concordant'/'discordant') "
        "and 'L2_moa' (case-normalized) as functional equivalents. "
        "Proceeding with these substitutes."
    )
    if "moa_concordance_status" not in cat.columns:
        raise ValueError(
            "STOP — neither 'L2_label' nor 'moa_concordance_status' available. "
            "Cannot determine MoA concordance."
        )

print("  All required columns present. Proceeding.")


def get_layer_rates(df: pd.DataFrame) -> dict[str, float]:
    """Return single-layer approval rates from Chain_Pattern_Distribution.

    Parameters
    ----------
    df : pd.DataFrame
        Chain_Pattern_Distribution with columns chain_pattern, n, n_approved.

    Returns
    -------
    dict[str, float]
        Keys L1 … L4; values are n_approved / n for the *_only pattern.
    """
    rates: dict[str, float] = {}
    for layer, pat in [("L1", "L1_only"), ("L2", "L2_only"),
                       ("L3", "L3_only"), ("L4", "L4_only")]:
        row = df[df["chain_pattern"] == pat]
        if row.empty:
            raise ValueError(f"Pattern '{pat}' not found in Chain_Pattern_Distribution")
        n = float(row["n"].iloc[0])
        k = float(row["n_approved"].iloc[0])
        rates[layer] = k / n
    return rates


chain_dist = pd.read_csv(CHAIN_DIST_PATH)
print(f"\nChain_Pattern_Distribution columns: {list(chain_dist.columns)}")
p = get_layer_rates(chain_dist)
for lyr, rate in p.items():
    print(f"  P_{lyr}_only = {rate:.6f}  ({100*rate:.4f}%)")

E_BLISS = 1.0 - float(
    (1 - p["L1"]) * (1 - p["L2"]) * (1 - p["L4"])
)
print(f"E_BLISS(L1+L2+L4) = {E_BLISS:.6f}")

null_df   = pd.read_csv(PERM_NULL_PATH)
print(f"\nPermutation_Null_SI_v3 columns: {list(null_df.columns)}")
if "SI_permuted" not in null_df.columns:
    raise ValueError("Missing column 'SI_permuted' in Permutation_Null_SI_v3")
null_si   = null_df["SI_permuted"].values
null_95th = float(np.percentile(null_si, 95))
print(f"Null 95th = {null_95th:.4f}")


# ---------------------------------------------------------------------------
# Helper: compute SI row dict for a group
# ---------------------------------------------------------------------------
def make_si_row(
    step: str,
    case: str,
    n: int,
    k: int,
    e_bliss: float,
    null_si_arr: np.ndarray,
    null_95th_val: float,
    note: str = "",
    combined: bool = False,
    n_obs: int = N_OBS,
    k_obs: int = K_OBS,
) -> dict:
    """Build a result row dict with SI, CI, and p-value.

    Parameters
    ----------
    step       : "A" or "B"
    case       : group label string
    n          : pool size (or n_pred for combined)
    k          : successes (or k_pred for combined)
    e_bliss    : E_BLISS denominator
    null_si_arr: permutation null SI array
    null_95th_val: null 95th percentile
    note       : free-text note
    combined   : if True, use (k_obs+k)/(n_obs+n) formula
    n_obs, k_obs: fixed obs pool for combined pool calculation
    """
    if combined:
        n_total = n_obs + n
        k_total = k_obs + k
        rate    = k_total / max(1, n_total)
    else:
        n_total = n
        k_total = k
        rate    = k_total / max(1, n_total)

    if n_total == 0 or e_bliss == 0:
        return {
            "step": step, "case": case,
            "n": n, "n_success": k,
            "success_rate_pct": None,
            "SI": None, "SI_CI_lo": None, "SI_CI_hi": None,
            "null_95th": round(null_95th_val, 4),
            "p_value": None, "is_significant": None,
            "note": "n=0, cannot compute SI",
        }

    si      = rate / e_bliss
    p_val   = float((null_si_arr >= si).sum()) / len(null_si_arr)
    is_sig  = bool(si > null_95th_val)

    ci_lo_rate, ci_hi_rate = proportion_confint(
        k_total, n_total, alpha=0.05, method="wilson"
    )
    si_lo = ci_lo_rate / e_bliss
    si_hi = ci_hi_rate / e_bliss

    if n < 20:
        warnings.warn(f"  WARNING: n={n} < 20 for '{case}' (Step {step}). SI unreliable.")
        note = (note + "; n<20, SI unreliable").lstrip("; ")

    return {
        "step":             step,
        "case":             case,
        "n":                n_total if combined else n,
        "n_success":        k_total if combined else k,
        "success_rate_pct": round(rate * 100, 4),
        "SI":               round(si, 4),
        "SI_CI_lo":         round(si_lo, 4),
        "SI_CI_hi":         round(si_hi, 4),
        "null_95th":        round(null_95th_val, 4),
        "p_value":          round(p_val, 4),
        "is_significant":   is_sig,
        "note":             note,
    }


# ===========================================================================
# Step A: Experimental concordance
# ===========================================================================
print("\n--- Step A: Experimental concordance ---")

l1l2l4 = cat[
    (cat["L1_state"] == "Active") &
    (cat["L2_state"] == "Active") &
    (cat["L4_state"] == "Active")
].drop_duplicates(subset="InChIKey").copy()

l1l2l4["max_phase"] = (
    pd.to_numeric(l1l2l4["max_phase"], errors="coerce").fillna(0).astype(int)
)

print(f"L1+L2+L4 pool (state filter, dedup): {len(l1l2l4):,}")

# Concordance groups via moa_concordance_status
# (L2_label substitute: 'concordant' ≡ L2_moa == human_dir,
#  'discordant' ≡ L2_moa != human_dir, both must be agonist/antagonist)
grp_concordant  = l1l2l4[l1l2l4["moa_concordance_status"] == "concordant"].copy()
grp_discordant  = l1l2l4[l1l2l4["moa_concordance_status"] == "discordant"].copy()
grp_directional = pd.concat([grp_concordant, grp_discordant], ignore_index=True)

rows: list[dict] = []

for label, grp in [
    ("Exp_All",        grp_directional),
    ("Exp_Concordant", grp_concordant),
    ("Exp_Discordant", grp_discordant),
]:
    n = len(grp)
    k = int((grp["max_phase"] == 4).sum())
    rate = k / max(1, n)
    si   = rate / E_BLISS
    note = "L2_label→moa_concordance_status"
    if label == "Exp_All":
        note += "; concordant+discordant"

    print(f"\n  {label}: n={n:,}, phase=4={k}, rate={rate*100:.4f}%, SI={si:.4f}")
    row = make_si_row("A", label, n, k, E_BLISS, null_si, null_95th, note=note)
    rows.append(row)
    print(f"    SI={row['SI']}, 95%CI=[{row['SI_CI_lo']}, {row['SI_CI_hi']}], "
          f"p={row['p_value']}, sig={row['is_significant']}")

# BACC-SI link
SI_conc = rows[1]["SI"]   # Exp_Concordant
SI_disc = rows[2]["SI"]   # Exp_Discordant
SI_BACC_est = 0.875 * SI_conc + 0.125 * SI_disc
SI_pred     = 7.81
match       = abs(SI_BACC_est - SI_pred) / SI_pred < 0.20
print(f"\n  BACC-SI link: estimate={SI_BACC_est:.3f}, predicted={SI_pred:.3f}, "
      f"within_20pct={match}")
print(f"    Note: SI_BACC_est uses experimental concordant/discordant approval rates.")
print(f"    SI_pred=7.81 is from combined-pool denovo analysis (t=0.99).")
print(f"    These pools differ; the BACC-SI link is indicative, not exact.")

# ===========================================================================
# Step B: Denovo predicted concordance (combined pool)
# NOTE: Step B n_pred is small (Case_B~10, Case_C~3) because human_dir
# annotation is absent for most denovo prediction targets.
# This is independent of confidence threshold.
# Step B is EXPLORATORY ONLY — not used for primary claims.
# Primary BACC-SI connection: see threshold sweep (draw_figure6.py Panel B).
# ===========================================================================
print("\n--- Step B: Denovo predicted concordance (combined pool) ---")

# Permutation targets from v3: L1+Active, L4+Active, L2 not Active/Inactive
perm_targets = cat[
    (cat["L1_state"] == "Active") &
    (cat["L4_state"] == "Active") &
    (~cat["L2_state"].isin(["Active", "Inactive"]))
].drop_duplicates(subset="InChIKey").copy()

perm_targets["max_phase"] = (
    pd.to_numeric(perm_targets["max_phase"], errors="coerce").fillna(0).astype(int)
)

print(f"Permutation targets (L1+/L4+/L2 unlabeled): {len(perm_targets):,}")

# Load denovo predictions
denovo_pred = pd.read_csv(DENOVO_PRED_PATH)
print(f"\nDenovo predictions columns: {list(denovo_pred.columns)}")
print(f"Denovo predictions shape: {denovo_pred.shape}")

for col in ("InChIKey", "is_high_confidence", "predicted_moa"):
    if col not in denovo_pred.columns:
        raise ValueError(f"Missing column '{col}' in denovo predictions")

# Aggregate to InChIKey level:
#   predicted_active = any HC row with predicted_moa != 'non-binder'
#   pred_dir = majority agonist/antagonist among HC-active rows
hc = denovo_pred[denovo_pred["is_high_confidence"] == True].copy()
print(f"HC rows: {len(hc):,}, unique IKs: {hc['InChIKey'].nunique():,}")

hc_active = hc[hc["predicted_moa"] != "non-binder"].copy()
active_iks = set(hc_active["InChIKey"].unique())
print(f"HC-active IKs: {len(active_iks):,}")


def majority_direction(group: pd.DataFrame) -> Optional[str]:
    """Return majority predicted direction among HC-active rows.

    Parameters
    ----------
    group : pd.DataFrame
        Rows for one InChIKey, all HC-active.

    Returns
    -------
    str or None
        'agonist', 'antagonist', or None if no directional calls.
    """
    dirs = [
        m for m in group["predicted_moa"]
        if m in ("agonist", "antagonist")
    ]
    if not dirs:
        return None
    cnt = Counter(dirs)
    return cnt.most_common(1)[0][0]


pred_dir_df = (
    hc_active
    .groupby("InChIKey", group_keys=False)
    .apply(majority_direction, include_groups=False)
    .reset_index()
)
pred_dir_df.columns = ["InChIKey", "pred_dir"]
print(f"pred_dir summary: {pred_dir_df['pred_dir'].value_counts().to_dict()}")

# Build per-IK summary
ik_summary = pd.DataFrame({
    "InChIKey":        list(active_iks),
    "predicted_active": True,
})
ik_summary = ik_summary.merge(pred_dir_df, on="InChIKey", how="left")

# Merge perm_targets × predictions (inner join — only IKs with predictions)
merged = perm_targets.merge(ik_summary, on="InChIKey", how="inner")
print(f"\nMerged (perm_targets ∩ predictions): {len(merged):,}")

# Active subset (Case A): predicted_active == True (all rows in merged are active)
case_A = merged.copy()

# Case B: active AND pred_dir == human_dir (human_dir not empty)
case_B = merged[
    (merged["pred_dir"].notna()) &
    (merged["pred_dir"] != "") &
    (merged["human_dir"] != "") &
    (merged["pred_dir"] == merged["human_dir"])
].copy()

# Case C: active AND pred_dir != human_dir
case_C = merged[
    (merged["pred_dir"].notna()) &
    (merged["pred_dir"] != "") &
    (merged["human_dir"] != "") &
    (merged["pred_dir"] != merged["human_dir"])
].copy()

print(f"\nCase A (any MoA active):    n={len(case_A):,}")
print(f"Case B (concordant):         n={len(case_B):,}")
print(f"Case C (discordant):         n={len(case_C):,}")

for label, grp, note in [
    ("Case_A", case_A, "denovo HC-active; combined pool (obs+pred)"),
    ("Case_B", case_B, "denovo concordant (pred_dir==human_dir); n_pred=10; exploratory only"),
    ("Case_C", case_C, "denovo discordant (pred_dir!=human_dir); n_pred=3; exploratory only"),
]:
    n_pred = len(grp)
    k_pred = int((grp["max_phase"] == 4).sum())
    si_comb = (K_OBS + k_pred) / (N_OBS + n_pred) / E_BLISS

    print(f"\n  {label}: n_pred={n_pred:,}, k_pred={k_pred}, "
          f"SI(combined)={si_comb:.4f}")

    if n_pred < 20:
        print(f"    WARNING: n_pred={n_pred} < 20. SI unreliable.")

    row = make_si_row(
        "B", label,
        n=n_pred, k=k_pred,
        e_bliss=E_BLISS,
        null_si_arr=null_si,
        null_95th_val=null_95th,
        note=note,
        combined=True,
        n_obs=N_OBS, k_obs=K_OBS,
    )
    rows.append(row)
    print(f"    n_combined={row['n']}, k_combined={row['n_success']}, "
          f"SI={row['SI']}, 95%CI=[{row['SI_CI_lo']}, {row['SI_CI_hi']}], "
          f"p={row['p_value']}, sig={row['is_significant']}")

# ===========================================================================
# Save
# ===========================================================================
print("\n--- Save ---")

COL_ORDER = [
    "step", "case", "n", "n_success", "success_rate_pct",
    "SI", "SI_CI_lo", "SI_CI_hi",
    "null_95th", "p_value", "is_significant", "note",
]
out_df = pd.DataFrame(rows)[COL_ORDER]

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)
print(out_df.to_string(index=False))

csv_path = f"{TABLE_DIR}/SI_MoA_aware.csv"
out_df.to_csv(csv_path, index=False)
print(f"\nSaved: {csv_path}  ({len(out_df)} rows)")

print("\n=== TASK 3 COMPLETE ===")
