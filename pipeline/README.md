# SynerGPCR — data pipeline

This directory is the archival record of the code that produced the released SynerGPCR
dataset. It contains only the scripts and notebooks that generate a number appearing in the
manuscript, in the order the pipeline ran. It is **not** a one-command rebuild: raw source
databases are not redistributable, the LLM annotation stage is non-deterministic and calls a
paid external API, and the GPCRact neural model lives in a separate repository and requires a
GPU. What this directory does provide is a complete and auditable path from each published
number back to the code and the prompt that produced it.

All input paths resolve against a configurable root:

```bash
export SYNERGPCR_BASE=/path/to/released/tables    # defaults to ./data
```

Environment: see `requirements.txt` (Python 3.9.23, pandas 2.2.2, numpy 1.26.4,
statsmodels 0.14.5, scikit-learn 1.6.1).

---

## Reproducibility tiers

| Tier | Meaning | Stages |
|---|---|---|
| **1** | Deterministic. Re-runnable from the released tables; repeated runs give identical output. | `03_chain_assembly/` (pair chains, 396,169 records) · `04_synergy/` (Bliss SI, Wilson CI, MoA concordance) · `05_gpcract_integration/` (three-point comparison, Monte Carlo null) · `06_target_enrichment/` (pooled baseline, Wilson lower bounds) · `07_derived_views/` (disease and repurposing views) · `01_harmonisation/source_integration.py` (record counts) |
| **2** | LLM-dependent, non-deterministic. Temperature 0.0–0.1, no seed available through the API; re-running will not reproduce labels exactly. Auditable instead via the released input→output mapping tables and the verbatim prompts in `02_llm_annotation/prompts/`. | `02_llm_annotation/` (all 7 stages) · `01_harmonisation/assay_moa_integration.ipynb` and `indication_tables.ipynb`, which consume Tier-2 output |
| **3** | Not re-runnable here. | Raw source downloads (DrugBank requires account approval; ChEMBL/PubChem/BindingDB/GLASS bulk downloads) · GPCRact model training and inference (separate repository, GPU required; its predictions enter this pipeline as frozen CSVs) |

**Auditing Tier 2.** Three released tables retain the model input alongside the model output,
row by row: `ChEMBL_llm_parsed.csv` (443,285 rows, `description` column retained),
`BindingDB_llm_parsed.csv` (415,712 rows, `DESCRIPTION` retained), and
`UMLS_CUI_to_Category_v2.csv` (5,079 rows, `disease_name` retained). The PubChem equivalent
(`aid_information_llm_parsed_v2.csv`, 35,216 rows) is keyed by AID and joins back to the
assay text through `aid_information_json.csv`. Every prompt, model ID, temperature, schema
enum, retry policy and failure fallback is recorded in `02_llm_annotation/prompts/`, which is
written to stand alone without the notebooks.

---

## Claim to script

| Manuscript claim | Script | Tier |
|---|---|---|
| 1,753,921 assay + clinical MoA records across 10 source databases | `01_harmonisation/source_integration.py` | 1 |
| Source MoA harmonisation: 6-class mapping, Tier-1 > Tier-1.5 priority, conflict resolution | `01_harmonisation/assay_moa_integration.ipynb` | 2 |
| Assay layer assignment L1–L4 | `01_harmonisation/assay_moa_integration.ipynb` (`map_assay_layer`) | 2 |
| Disease terms → UMLS CUI → 12 therapeutic areas | `02_llm_annotation/indication_to_umls.ipynb`, `umls_to_category.ipynb`; `01_harmonisation/indication_tables.ipynb` | 2 |
| 396,169 compound × target pair chain records (277,544 compounds, 600 GPCRs) | `03_chain_assembly/build_pair_chains.py` → `fix_pairlevel.py` (asserts all three) | 1 |
| Full_chain n = 228 at compound level (union across targets) | `03_chain_assembly/fix_pairlevel.py`; `04_synergy/compute_synergy.py` | 1 |
| Bliss Synergy Index + Wilson 95% CI, all 15 chain patterns (L1+L2+L4 = 9.71, L1+L3+L4 = 13.01, Full_chain = 13.45) | `04_synergy/compute_synergy.py` | 1 |
| Single-layer approval-rate baselines (L1–L4) | `04_synergy/compute_synergy.py` | 1 |
| MoA concordance 123/385 (31.9%) vs 25/176 (14.2%), OR 2.84 | `04_synergy/moa_concordance.py`; concordance labels built by `04_synergy/build_catalog.py` | 1 |
| Operating point L2 ≥ 0.40 / L3 ≥ 0.97 | `05_gpcract_integration/final_operating_point.py` | 1 |
| Three-point comparison and recovery (L2 1.58 → 7.74 → 9.71, 79.7%; L3 2.20 → 10.35 → 13.01, 79.6%) | `05_gpcract_integration/final_operating_point.py`, `null_and_si.py` | 1 |
| Monte Carlo random-assignment null (null95 7.26 L2 / 7.36 L3) | `05_gpcract_integration/null_and_si.py` | 1 |
| Pooled conditional-approval baseline 1.40% over 20,145 co-active pairs across 211 targets | `06_target_enrichment/target_enrichment.py` | 1 |
| 126 analysable targets (n ≥ 5); 28 with Wilson 95% lower bound above the pooled baseline | `06_target_enrichment/target_enrichment.py` | 1 |
| High-n / zero-approved targets are a co-active coverage artefact, not absence of approved drugs | `06_target_enrichment/highn_validation.py`, `ontarget_reclass.py` | 1 |
| Disease-category views (target × 12 areas, 1,512 cells) | `07_derived_views/disease_views.py` | 1 |
| 458 cross-family repurposing records / 146 approved drugs | **producer script not retained** — see Known limitations; derivation in `reports/STEP5e*.md`, upstream tiering in `07_derived_views/repurposing_derived.py` | 1 |

---

## Known limitations

**1. `threshold_sweep_superseded.py` evaluates a superseded operating point.**
The retained sweep script scans thresholds and self-checks against L2 t = 0.50 and L3 t = 0.85
(SI 7.73 / 8.89). Those are **not** the released values. The released operating point is
L2 ≥ 0.40 / L3 ≥ 0.97, and its justification is the null comparison recorded in
`gpcract_3point_si.csv` — at the released thresholds the panel SI exceeds the Monte Carlo
null95 (L2 7.74 vs 7.26; L3 10.35 vs 7.36). The file is included because it is the only
retained record of the threshold-selection procedure; it is named `_superseded` so it is not
mistaken for the final selection. Do not quote its SI or null values.

**2. The producer script for `repurposing_view_crossfam.csv` was not retained.**
The shipped view (458 records, 146 approved drugs) was assembled by an inline script that was
never written to a file. The derivation is documented step by step in
`reports/STEP5e1_tier_canonical.md`, `STEP5e2_stat_recompute.md`, `STEP5e3_repurposing.md`,
`STEP5e3b_repurposing_final.md`, `STEP5e3c_repurposing_final.md`, `STEP5e5c_final_export.md`
and `STEP5e5e_final.md`, and all of its inputs survive. The output is released as a frozen
artefact. `07_derived_views/repurposing_derived.py` produces the *upstream* view
(`repurposing_view.csv`, 842 approved drugs), not the curated cross-family subset.

**3. Monte Carlo null: SEED = 0, 10,000 replicates, re-seeded per evaluation.**
`np.random.default_rng(SEED)` is called inside the null function, so every threshold and every
panel draws the identical replicate stream. This is a deliberate variance-reduction choice —
it makes differences between thresholds reflect the data rather than sampling noise — and it
means the per-threshold nulls are perfectly correlated by construction. Bit-level
reproducibility requires `numpy==1.26.4`; NumPy does not guarantee PCG64 stream stability
across major versions. A separate, older lineage (`analysis13_*`, not included here) used
SEED = 42 with 1,000 replicates; those values are not comparable and are not released.

**4. Tier-2 failure modes are not neutral.** Each LLM stage collapses failures to a specific
label rather than dropping the row: the assay stages to `Unknown`, the MoA stages to `binder`,
the category stage to `12 (Other / Unclassified)` for an entire failed batch of up to 30. The
size of those classes therefore mixes genuine cases with API failures and must not be read as
a biological result. Details per stage are in `02_llm_annotation/prompts/`.

**5. `01_harmonisation/source_integration.py` carries two stale hardcoded columns.**
`n_unique_pairs` is hardcoded to 260,938 and `n_targets` is read from a superseded master
table. Only `n_raw_records` (and its total, 1,753,921) is computed and current. The
authoritative pair count is 396,169 from `03_chain_assembly/`, and the authoritative target
count is 600.

---

## Directory layout

```
01_harmonisation/       source integration, MoA and indication master tables
02_llm_annotation/      Gemini annotation notebooks + prompts/ (audit trail)
03_chain_assembly/      compound × target chain records, layer states
04_synergy/             Bliss SI, Wilson CI, MoA concordance, catalog build
05_gpcract_integration/ operating point, Monte Carlo null, three-point comparison
06_target_enrichment/   co-active conditional approval, enrichment, validation
07_derived_views/       disease-category and repurposing views
reports/                step reports documenting derivations and decisions
```

`reports/` is copied from the analysis log. Two edits were applied for release: build-script
references were repointed from their scratch paths to the released paths in this tree (or
marked "build script not released"), and one absolute path was replaced with `<repo-root>`.
No numeric content was altered.
