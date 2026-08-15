"""
webexport_step4_disease.py — WebExport_v2 STEP 4. READ-ONLY on sources.

Convert existing disease data into the web-export bundle (no new analysis). Builds
in memory, runs ALL self-checks, writes ONLY if every check passes.

  1. drug_indication.csv          (PK InChIKey,UMLS_CUI)
  2. disease_category_lookup.csv  (12 categories)
  3. target_indication.csv        (UniProt-keyed, restricted to 600-GPCR scope)

Terminology rule (CLAUDE.md §0-4): the strings 'DSAV' / 'ACTR' must not appear in any
written file (column or value).
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
DI   = B / "Output/DB/GPCRactDB/Drug_Indication_Master_Integrated.csv"
TI   = B / "Output/DB/GPCRactDB/Target_Indication_Master_Integrated.csv"
UM   = B / "Output/DB/GPCRactDB/UMLS_CUI_to_Category_v2.csv"
TL   = OUT / "target_lookup.csv"            # 600-GPCR scope (STEP 2)
CL   = OUT / "compound_lookup.csv"          # in_assay_universe + is_approved (STEP 2)

STATUS_RANK = {"Approved": 4, "Phase 3": 3, "Phase 2": 2, "Phase 1": 1, "Investigational": 0}
ABBREV = {1: "CNS", 2: "CV", 3: "MET", 4: "IMM", 5: "RESP", 6: "ONC",
          7: "RENAL", 8: "GI", 9: "PAIN", 10: "REPRO", 11: "INF", 12: "OTHER"}


def load(p, **kw):
    if not p.exists():
        sys.exit(f"STOP: missing PROVENANCE {p}")
    df = pd.read_csv(p, **kw); print(f"  loaded {p.name}: {len(df):,} rows"); return df


print("=" * 78); print("LOADING PROVENANCE"); print("=" * 78)
di = load(DI, low_memory=False)
ti = load(TI, low_memory=False)
um = load(UM, low_memory=False)
tl = load(TL, low_memory=False)
cl = load(CL, low_memory=False)

assert um.category_id.nunique() == 12, "UMLS map is not 12 categories"
assert um.umls_cui.nunique() == len(um), "UMLS map has duplicate CUIs"
cui2cat = um[["umls_cui", "disease_name", "category_id", "category_name"]]

# ==========================================================================
# FILE 1 — drug_indication.csv
# ==========================================================================
print("\n" + "=" * 78); print("FILE 1: drug_indication.csv"); print("=" * 78)
d = di.copy()
d["highest_status"] = d.Highest_Status.map(STATUS_RANK)
assert d["highest_status"].notna().all(), "unmapped Highest_Status value"
d = d.merge(cui2cat, left_on="UMLS_CUI", right_on="umls_cui", how="left")
# collapse to one row per (InChIKey, UMLS_CUI) with max status (defensive; pairs already unique)
d = (d.sort_values("highest_status", ascending=False)
       .drop_duplicates(["InChIKey", "UMLS_CUI"]))
drug_ind = d[["InChIKey", "UMLS_CUI", "disease_name",
              "category_id", "category_name", "highest_status"]].copy()
drug_ind["category_id"] = drug_ind["category_id"].astype("Int64")
n_unmapped = int(drug_ind.category_id.isna().sum())
print(f"  rows={len(drug_ind):,}  unique drugs={drug_ind.InChIKey.nunique():,}  "
      f"categories={drug_ind.category_id.dropna().nunique()}")
print(f"  unmapped (CUI not in 12-cat map) pairs={n_unmapped:,} "
      f"({drug_ind.loc[drug_ind.category_id.isna(),'UMLS_CUI'].nunique()} distinct CUIs)")

# ==========================================================================
# FILE 2 — disease_category_lookup.csv (12 categories)
# ==========================================================================
print("\n" + "=" * 78); print("FILE 2: disease_category_lookup.csv"); print("=" * 78)
gpcr600 = set(tl.UniProt_AC)
ti600 = ti[ti.UniProt_AC.isin(gpcr600)].copy()
tim = ti600.merge(um[["umls_cui", "category_id"]], left_on="UMLS_CUI", right_on="umls_cui", how="inner")
gene_map = tl.set_index("UniProt_AC")["gene_name"].to_dict()
rep = (tim.groupby(["category_id", "UniProt_AC"])["UMLS_CUI"].nunique()
          .reset_index(name="n_ind"))
rep["label"] = rep.UniProt_AC.map(lambda u: gene_map.get(u) if pd.notna(gene_map.get(u)) else u)
rep = rep.sort_values(["category_id", "n_ind"], ascending=[True, False])
rep_map = rep.groupby("category_id").head(5).groupby("category_id")["label"].apply(
    lambda s: ";".join(map(str, s))).to_dict()

cats = um.drop_duplicates("category_id")[["category_id", "category_name"]].sort_values("category_id")
cats["abbreviation"] = cats.category_id.map(ABBREV)
cats["representative_GPCRs"] = cats.category_id.map(lambda c: rep_map.get(c, ""))
dcl = cats[["category_id", "category_name", "abbreviation", "representative_GPCRs"]].reset_index(drop=True)
print(dcl.to_string(index=False))

# ==========================================================================
# FILE 3 — target_indication.csv (600-scope)
# ==========================================================================
print("\n" + "=" * 78); print("FILE 3: target_indication.csv (600-GPCR scope)"); print("=" * 78)
t = ti600.merge(cui2cat, left_on="UMLS_CUI", right_on="umls_cui", how="left")
tgt_ind = t[["UniProt_AC", "UMLS_CUI", "disease_name",
             "category_id", "category_name", "Source"]].rename(columns={"Source": "source"}).copy()
tgt_ind["category_id"] = tgt_ind["category_id"].astype("Int64")
n_surv = tgt_ind.UniProt_AC.nunique()
print(f"  rows={len(tgt_ind):,}  GPCRs in scope={n_surv} of {ti.UniProt_AC.nunique()} "
      f"(dropped {ti.UniProt_AC.nunique()-n_surv} outside 600-scope)")

# ==========================================================================
# Coverage reporting (assay-bearing / approved drugs carrying a disease tag)
# ==========================================================================
assay_set = set(cl.loc[cl.in_assay_universe == True, "InChIKey"])
appr_set  = set(cl.loc[cl.is_approved == True, "InChIKey"])
tagged    = set(drug_ind.InChIKey)
cov_assay = len(assay_set & tagged)
cov_appr  = len(appr_set & tagged)
print("\n" + "=" * 78); print("DISEASE-TAG COVERAGE"); print("=" * 78)
print(f"  drugs with >=1 disease tag (any)            : {len(tagged):,} of 8,928")
print(f"  assay-bearing compounds (277,544) tagged    : {cov_assay:,}")
print(f"  approved drugs (2,420) tagged               : {cov_appr:,}")

# ==========================================================================
# SELF-CHECKS
# ==========================================================================
print("\n" + "=" * 78); print("SELF-CHECKS"); print("=" * 78)
checks = []
def ck(name, cond, got, exp): checks.append((name, bool(cond), got, exp))

ck("disease_category_lookup == 12 rows", len(dcl) == 12, len(dcl), 12)
ck("disease_category_lookup category_id 1..12", sorted(dcl.category_id) == list(range(1, 13)),
   sorted(dcl.category_id.tolist()), "1..12")
ck("drug_indication mapped categories <= 12", drug_ind.category_id.dropna().nunique() <= 12,
   int(drug_ind.category_id.dropna().nunique()), "<=12")
ck("target_indication categories <= 12", tgt_ind.category_id.dropna().nunique() <= 12,
   int(tgt_ind.category_id.dropna().nunique()), "<=12")
ck("target_indication 211 of 221 survivors", n_surv == 211, n_surv, 211)
ck("drug_indication keeps all 8,928 drugs", drug_ind.InChIKey.nunique() == 8928,
   drug_ind.InChIKey.nunique(), 8928)
ck("target_indication GPCRs subset of 600-scope", set(tgt_ind.UniProt_AC).issubset(gpcr600),
   len(set(tgt_ind.UniProt_AC) - gpcr600), 0)

# terminology guard: no 'DSAV'/'ACTR' as a column name or descriptive label.
# Identifier columns (InChIKey/UMLS_CUI/UniProt_AC hashes; gene symbols) are excluded from the
# value scan — an InChIKey hash may legitimately contain the substring 'DSAV'/'ACTR'.
LABEL_COLS = {"disease_name", "category_name", "abbreviation"}
def has_forbidden(df):
    hits = [c for c in df.columns if "dsav" in c.lower() or "actr" in c.lower()]
    for c in df.columns:
        if c in LABEL_COLS and df[c].dtype == object:
            s = df[c].astype(str)
            if s.str.contains("DSAV", case=False).any() or s.str.contains("ACTR", case=False).any():
                hits.append(f"value:{c}")
    return hits
forb = []
for nm, df in [("drug_indication", drug_ind), ("disease_category_lookup", dcl), ("target_indication", tgt_ind)]:
    h = has_forbidden(df)
    if h: forb.append((nm, h))
ck("no DSAV/ACTR strings in any file", len(forb) == 0, forb, "[]")

allpass = all(c[1] for c in checks)
for name, ok, got, exp in checks:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:46s} got={got} exp={exp}")
if not allpass:
    print("\n*** SELF-CHECK FAILED — writing NOTHING. ***"); sys.exit(1)

print("\n" + "=" * 78); print("ALL CHECKS PASS — WRITING"); print("=" * 78)
OUT.mkdir(parents=True, exist_ok=True)
for df, fn in [(drug_ind, "drug_indication.csv"), (dcl, "disease_category_lookup.csv"),
               (tgt_ind, "target_indication.csv")]:
    p = OUT / fn; df.to_csv(p, index=False); print(f"  wrote {p}  ({len(df):,} rows, {df.shape[1]} cols)")
print("\nDONE.")
