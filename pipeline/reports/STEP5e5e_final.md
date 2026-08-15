# STEP 5e-5e — final two files + CLAUDE.md edit — **DONE (both asserts passed, exact)**

target_analysis sub-task. Sources **READ-ONLY** except the two governance files PART 3 authorizes
(`target_analysis_CLAUDE.md`, repo-root `CLAUDE.md`). Outputs only under `final_export/`. Build: inline
scripts (pandas 2.2.2, numpy 1.26.4, statsmodels 0.14.5, Python 3.9, seed 0). Both 5e-5c STOPs are now
**resolved with exact oracle reproduction — no hardcode bypass.**

The prompt's PART 1/PART 2 skeletons assumed the annotation files carry backbone/exp-state columns; they do
**not** (pure prediction tables). The canonical pool/`base2`/`ceil2`/`scope2` logic lives in
`Output/NAR_DBI/Code/part11_final_op_0p40_0p97.py` (a FINAL_STATS PROVENANCE script). I reproduced **7.74**
and **1,179** from that verified logic + the STEP0p6 binder gate (`binding_prob ≥ 0.5`), not from the skeleton.

---

## ① PART 0 — annotation file paths
- `…/GPCRactDB/results/denovo/synergpcr_annotation_L2.csv` (673,719 rows)
- `…/GPCRactDB/results/denovo/synergpcr_annotation_L3.csv` (673,719 rows)
- Schema (both): `InChIKey, UniProt_AC, binding_prob, ant_prob, ago_prob, predicted_moa, confidence_raw,
  is_high_confidence` — **prediction-only**; no backbone/exp-state. Pool logic sourced from
  `part11_final_op_0p40_0p97.py` (+ `Compound_Chain_Catalog_v3`, `L3L4_Pair_Catalog_v1`, `Bioassay_MoA_Master`).

## ② PART 1 — gpcract_3point_si.csv L2 row → **reproduced EXACTLY, written**
Pool = `ceil2` (L1&L2&L4 exp-active) ∪ `up_l2` (binder-gate AI-completion: `base2` ∩ `scope2` with
`confidence_raw ≥ 0.40 AND binding_prob ≥ 0.5`), per part11 + STEP0p6.

| quantity | recompute | §D / STEP0p6 target | assert |
|---|---:|---:|:--:|
| n_pool | **3,826** | 3,826 | ✓ (Δ≤5) |
| n_approved | **195** | 195 | ✓ (Δ≤3) |
| obs_rate | 0.050967 | — | — |
| E_Bliss | 0.006586 | 0.006586 | ✓ |
| **SI** | **7.7382** | 7.74 | ✓ (Δ<0.05) |
| Wilson SI-CI | **[6.7473, 8.8669]** | [6.747, 8.866] | ✓ |
| null95 (MC random-assignment, 10k reps, seed 0) | **7.2620** | 7.26 | ✓ (Δ<0.05) |
| recovery (vs ceiling 9.7099) | **79.7%** | 79.7% | ✓ (Δ<0.5) |

- The L2 baseline (1.5808) and ceiling (9.7099) and all L3 rows were already correct (5e-5b) → unchanged.
- ⚠ **Two write-quality fixes made (reported, not silent):** (1) the prompt skeleton wrote to a new column
  `n_appr`; I corrected it to update the existing **`n_approved`** column (was stale 198 → now 195) and removed
  the stray column. (2) I also recomputed **`SI_ci_lo/hi`** from the actual binder-gate pool (195/3,826) so the
  row is internally consistent — matches STEP0p6's [6.747, 8.866] exactly. The prompt listed neither, but a row
  carrying SI 7.74 with the old pool's CI/approved would be incoherent.
- **Why 5e-5c got 7.96:** it used `db_load_tier_retained` (population-bias-filtered) ∩ backbone. The SI panel is
  **not** re-filtered through population-bias removal (STEP5e2 §2-3); the correct pool is the binder-gate
  annotation set. **7.74 is canonical and now independently re-derived.**

## ③ PART 2 — candidate_shortlist.csv (Option A, pair-level enabling-cell) → **1,179 EXACT, written**
Enabling cells = `up_l2`/`up_l3` pairs of AI-enabled compounds (in `base`, not already in `ceil2`), filtered to
`tier == T1_anchored` in `db_load_tier_final.tsv`, grouped by compound:

| bucket | count | expected |
|---|---:|---:|
| L2-only | **1,036** | 1,036 |
| L3-only | **7** | 7 |
| both L2 & L3 | **136** | 136 |
| **TOTAL** | **1,179** | 1,179 ✓ |

- File: **1,179 rows**, 11 web cols (incl. `smiles` for ECFP4 fallback). **smiles coverage = 1,179/1,179
  (100%)**. **approved = 24** (note: differs from the old conf-only file's 36 — the canonical enabling-cell
  T1-anchored definition is narrower; 24 is the definition-consistent count, reported not forced).
- Sorted approved-first then by max confidence; `ai_filled_layers` = `L2`/`L3`/`L2+L3`.

## ④ PART 3 — CLAUDE.md edits (PART 1 passed → unblocked)
Both files had the L2 panel value on a **single line** (grep-verified, no other occurrences) → scoped edits via
the Edit tool (not blind global replace):
- **(A) `target_analysis_CLAUDE.md`** line 51: `SI 7.67` → **`SI 7.74`** (1 occurrence). `7.17`/`79.0` absent.
- **(B) repo-root `CLAUDE.md`** line 43: `1.58 → 7.67 → 9.71 (null95 7.17, recovery 79.0%)` →
  **`1.58 → 7.74 → 9.71 (null95 7.26, recovery 79.7%)`** (single occurrence each; context confirmed = the §3
  L2 acceptance line, nothing else touched). This resolves the root-vs-subtask oracle conflict flagged in 5e-5b
  in favour of the now-verified binder-gate value. ⚠ **The repo-root master was edited** — a deliberate,
  authorized, single-line change; all other §3 values (incl. the still-conf-only line 45 pool 3,937/1,357 and
  line 46 171,833) were **left untouched** (out of PART 3 scope — separate reconciliation).

## ⑤ PART 4 — final_export/ verification (16 files)
| file | status | rows (data) | oracle | ✓ |
|---|---|---:|---|:--:|
| candidate_shortlist.csv | **REGEN** | 1,179 | §D/§E 1,179 | ✓ |
| gpcract_3point_si.csv | **REGEN** | 6 | §D L2 7.74/7.26/79.7 | ✓ |
| repurposing_view_crossfam.csv | REGEN (5e-5c) | 458 / 146 drugs | §E 458 visible | ✓ |
| repurposing_view_full.csv | REGEN (5e-5c) | 1,771 / 496 drugs | §E 496 | ✓ |
| gpcract_predictions.csv | REGEN (5e-5c) | 366,797 / 104,641 ret. | §E 104,641 | ✓ |
| compound_chain_summary.csv | REGEN (5e-5c) | 396,169 / 59,888 uniq. | §E 59,888 | ✓ |
| synergy_stats / single_layer_table / compound_lookup / target_lookup / drug_indication / target_indication / disease_category_lookup / target_clinical_enrichment / target_disease_enrichment / assay_clinical_validity | KEEP | — | — | ✓ |

**Checklist:** candidate_shortlist 1,179 ✓ · gpcract_3point_si L2 SI 7.74/null95 7.26/rec 79.7% ✓ ·
crossfam 458/146 ✓ · full 496 ✓ · predictions 104,641 retained ✓ · chain_summary 59,888 ✓. `figures/` empty
(later STEP).

## ⑥ 229 schema flags (carried)
- `target_lookup` lacks `gpcr_name`/`gpcr_class` → standardized in views to `protein_name` / `gpcr_class_name`
  (`gpcr_class_letter` also available) / `gpcr_family` / `gene_name`.
- `compound_chain_summary.gpcract_prediction_available` = compound-level coverage on a pair-level table
  (row-sum 122,656 vs unique 59,888) — decide pair-vs-compound semantics for the 229 contract.
- Repurposing UI: `repurposing_view_crossfam` (curated, 146 visible drugs) vs `repurposing_view_full` (496
  approved AI-profiles) — confirm whether the public "Approved Drugs / Repurposing" tab merges or separates them.
- `gpcract_predictions.csv` is a **backend** source: it includes REMOVED rows (`is_retained=False`,
  `removed_reason=modality_target_unanchored` etc.) — these internal reasons are backend metadata, not
  user-facing labels (and are not banned terms).

## ⑦ STALE scan
New scalars {3,826 · 195 · 7.7382 · 7.2620 · 79.7 · 6.7473 · 8.8669 · 1,179 · 1,036 · 24} ∩ §1.5 STALE = NONE.
Two grep matches investigated and cleared as **false positives**: (1) `ACTR`/`DSAV` occur **only** as substrings
inside random InChIKey hashes (e.g. `ACTRVOBWPAIOHC-…`, `DSAVHAONXZQUBN-…`) — never as a column/label
(confirmed: 0 hits in non-InChIKey fields); (2) `7.67` matched the substring of `7.6786` = the experimental
L1+L2+L3 approval-rate % in the KEEP `synergy_stats.csv`, **not** the retired L2 SI. No `0.50/0.85` op-point,
no 11,428, no recovery 82.5/92.3, no corrected-4.93, no forward/backward/permutation. Null = "Monte-Carlo
random-assignment". The retired 7.67/7.17/79.0% no longer appear anywhere in `final_export/` or the two
governance files (verified: both now read 7.74/7.26/79.7%).

## ⑧ 5e-6 / 5f prep
- **5e-6 (Methods paragraph):** the artifact-removal narrative (target exclusion → modality-anchoring →
  pan-family-hub K=5) and the SI-panel definition (binder-gate, pool 3,826/195, MC random-assignment null
  10k reps) are now fully reproducible and can be written up.
- **5f (EXPORT_MANIFEST):** `final_export/` holds 16 files (10 KEEP + 6 regenerated, all oracle-verified) —
  ready to manifest. Document the column contract + the three 229 schema flags above.
- Carried open items (unchanged): §E `3,719` funnel and §F `171/13`-vs-`170/14` ambiguities (filter-invariant,
  flagged); `no_evidenced_family` hide-filter scope decision; mail revisions for 5e-4
  (171,821→59,888; 842→496; "over a thousand"→~1,180; L2 SI 7.74).

Files written: `final_export/gpcract_3point_si.csv`, `final_export/candidate_shortlist.csv`; edited
`target_analysis_CLAUDE.md` (1 line) + repo-root `CLAUDE.md` (1 line). Report:
`target_analysis/STEP5e5e_final.md`. Stopping per the one-STEP rule.
