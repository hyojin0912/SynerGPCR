# STEP 2 — Task 1: target clinical enrichment — **DONE**

target_analysis sub-task. Sources READ-ONLY; outputs only in `target_analysis/`.
Build: `step2_prelim.py (build script not released)`, `06_target_enrichment/enrichment_prelim.py`, `step2_summ.py (build script not released)`.
Env: pandas 2.2.2 / numpy 1.26.4 / **statsmodels 0.14.5**, seed=0.
**All results are hypothesis-generating** — conditional-approval *association*, never a causal/quality
claim (approval is downstream of PK/tox/commercial factors outside the target). Internal terms
(DSAV/ACTR/permutation/forward/backward) not used; banned-term scan on both TSVs = **CLEAN**.

Locked definitions (STEP 1): experimental active = Bioassay_MoA_Master + part11c Active/Inactive map,
any-active→Active (100% chain-consistent with Catalog_v3). `L1_exp_active`=L1 Active(Binder);
`L2_exp_active`=L2 Active(agonist/antagonist-class). `L2_ai_active`=binder-gated annotation L2
(`binding_prob≥0.5 & confidence_raw≥0.40`). approved=Catalog `max_phase==4` (SI-panel consistent);
Drug_Indication 'Approved' carried as sensitivity column. Minimal unit=(compound,GPCR); counts dedup
unique compounds. Universe = 600-scope ∩ Target_Indication(221) = **211 targets**.

---

## Preliminary confirmations (before compute)

### P1. Assay-bearing-approved exact (mail representative approved count)
| quantity | value |
|---|---:|
| assay-bearing (defined-layer Bioassay) | 277,544 |
| Catalog approved (max_phase==4) | 2,420 |
| **max_phase==4 ∩ assay-bearing (defined-layer)** | **932** |
| max_phase==4 ∩ any-record incl Unknown layer | 2,420 |
- All 2,420 DB-approved drugs have ≥1 record, but only **932** have a **defined-layer (L1–L4)** assay;
  the other 1,488 are Unknown-layer-only. → For the mail, the *assay-bearing approved* figure is **932**
  (same "assay-bearing" definition as the 277,544); the DB-wide approved universe is 2,420. Report both
  with their definitions; don't conflate.

### P2. Funnel 277,544 → 169,635 → 3,719 → 1,263 (binder-gate)
| stage | n |
|---|---:|
| assay-bearing compounds | 277,544 |
| AI-completion (binder-gate, 600-scope) | 169,635 |
| …of which have experimental **L1 AND L4** record (backbone present) | 3,719 |
| AI-enabled synergy pool (backbone-constrained) | **1,263** (L2 route 1,246 + L3-only 17) |

**3,719 → 1,263 decomposition** (why most backbone-record compounds don't enter the synergy pool):
- **(a) no co-localized *active* backbone:** 1,691 of 3,719 have L1 and L4 records but **not** both
  **Active at the same GPCR with L2 unmeasured** → 2,028 remain backbone-eligible (`∃ GPCR: L1&L4 Active & L2-NA`).
- **(b) backbone present but no qualifying AI prediction:** 782 of those 2,028 have **no** binder-gate
  L2 AI pass at the backbone GPCR → **1,246** reach the L2 (L1+L2+L4) pool; +17 via the L3-only route = **1,263**.
- So the reduction is dominated by **(a) missing co-localized active backbone** (not a thresholding artifact).

---

## A. `target_clinical_enrichment.tsv` — 211 rows × 23 cols

Columns: `UniProt, GPCR_name, n_tested_L1_exp, n_concordant_exp, n_approved_exp, approval_rate_raw,
approval_rate_cond_exp, wilson_lo_exp, wilson_hi_exp, baseline_loto_exp, baseline_pool_exp,
log2_enrichment_exp, n_concordant_aiincl, n_approved_aiincl, approval_rate_cond_aiincl, wilson_lo_aiincl,
wilson_hi_aiincl, log2_enrichment_aiincl, enough_n_exp, enough_n_aiincl, n_approved_exp_drugind,
baseline_loto_aiincl, baseline_pool_aiincl`.

- `concordant_exp(t)` = `L1_exp_active ∩ L2_exp_active` @ t; `concordant_aiincl(t)` =
  `L1_exp_active ∩ (L2_exp_active OR L2_ai_active)` @ t.
- `approval_rate_raw` = approved among L1-active / n_tested_L1_exp (**survivorship-biased; NOT for
  ranking** — annotated in file usage).
- `log2_enrichment_*` = `log2(approval_rate_cond / baseline_loto)` (**LOTO primary**); `baseline_pool`
  given for sanity (pool 0.01400 vs LOTO ~0.0137–0.0139 — nearly identical, large stratum).
- Wilson CI denominator = n_concordant; **n_concordant < 5 → "insufficient data"**, enrichment blank.

### Key numbers
| metric | value |
|---|---|
| baseline_pool_exp (conditional approval over all concordant pairs) | **1.40%** |
| baseline_pool_aiincl | 1.11% |
| approval_rate_raw median (L1 pool) | 0.44% |
| approval_rate_cond_exp median (enough-n) | **1.00%** |
| enough_n_exp (≥5) | **126 / 211** (85 insufficient) |
| of the 126: computable log2 enrichment | 85 (the other **41 have 0 approved concordant** → rate 0, log2 undefined) |
| enriched (log2>0) / depleted (<0) | **50 / 35**; |log2|≥1 (≥2×): 50 |

Conditioning on L1∩L2 concordance roughly **doubles** the median approval rate vs the raw L1 pool
(0.44% → 1.00%), the intended demonstration that co-active (binding + G-protein) compounds are more
clinically translatable — stated as association.

### Top conditional-enrichment targets (exp, hypothesis-generating)
| UniProt | GPCR | n_concordant | n_approved | cond rate | log2 enrichment |
|---|---|---:|---:|---:|---:|
| Q8TDU6 | GPBAR1 (TGR5) | 5 | 3 | 0.600 | 5.44 |
| P35367 | HRH1 | 11 | 6 | 0.545 | 5.31 |
| P35368 | ADRA1B | 9 | 3 | 0.333 | 4.59 |
| P35348 | ADRA1A | 23 | 7 | 0.304 | 4.48 |
| P30556 | AGTR1 | 18 | 5 | 0.278 | 4.34 |
All are well-characterised druggable GPCRs with multiple approved drugs — a sanity check that the
conditional-enrichment ranking surfaces genuinely clinically-validated targets.

### AI-inclusion preserves the signal (secondary)
Among 126 targets enough-n in both exp & aiincl (85 with both log2 computable):
**sign-preserved in 73/85 (86%)**; median log2 exp 0.564 → aiincl 0.178 (signal retained, attenuated
because AI adds concordant compounds, lowering both the rate and the baseline). → "AI inclusion
**preserves** the enrichment direction" holds as a secondary demonstration.

### Approval-definition robustness (sensitivity)
`n_approved_exp_drugind` (Drug_Indication 'Approved') = `n_approved_exp` (max_phase==4) on **every**
target (totals 282 = 282; 0 targets differ). Within the assay-bearing concordant set the two approval
definitions **coincide** → Task-1 numbers are robust to the approval definition.

---

## B. `target_disease_enrichment.tsv` — 1,512 cells (n_concordant_exp ≥ 5)

Columns: `UniProt, GPCR_name, disease_category, n_concordant_exp, n_approved_cat, approval_rate_cond,
wilson_lo, wilson_hi, baseline_loto, baseline_pool, log2_enrichment`.
**Cell definition (default, stated):** for (target t, category d), numerator = concordant-at-t compounds
that are approved (`max_phase==4`) **and** carry a Drug_Indication category-d indication; denominator =
n_concordant_exp(t). Per-disease LOTO + pool baselines; Wilson CI denom = n_concordant_exp.

- **1,512 cells** over **126 targets × 12 categories** (only n_concordant_exp ≥ 5).
- enriched (log2>0): 447 · depleted (<0): 94 · zero-numerator (no approved-in-category): 971.
- Top cells recover known target–disease biology (sanity):

| GPCR | disease_category | n_concordant | n_approved_cat | cond rate | log2 |
|---|---|---:|---:|---:|---:|
| GPBAR1 (TGR5) | Metabolic / Endocrine | 5 | 3 | 0.600 | 7.75 |
| GPBAR1 | Gastrointestinal | 5 | 3 | 0.600 | 6.87 |
| HRH1 | Inflammatory / Immune | 11 | 5 | 0.455 | 6.73 |

GPBAR1/TGR5 → Metabolic and HRH1 → Inflammatory/allergy are textbook associations, confirming the cells
answer "for this disease, which GPCR is more strongly associated with clinical success" sensibly.

---

## STALE / guardrail scan
- Computed scalars {277,544 · 932 · 169,635 · 3,719 · 1,263 · 211} ∩ §1.5 STALE = **NONE**.
- `3,843` does **not** appear; op-points 0.40/0.97 only; no `corrected 4.93`, no `0.50/0.85`.
- Banned-term scan over both output TSVs = CLEAN.

## Caveats / next-STEP note
- 85/211 targets are insufficient-n (<5 concordant) and 41 more have n≥5 but 0 approved concordant
  (rate 0, log2 undefined, reported as enough_n=True with blank enrichment) — these are sparse-data
  targets, **not** evidence of "no association"; flagged so they're not mis-ranked.
- All claims hypothesis-generating; the file annotates `approval_rate_raw` as survivorship-biased / not
  for ranking. Stopping per the one-STEP rule — awaiting your next instruction (Task 2 / Task 3).
