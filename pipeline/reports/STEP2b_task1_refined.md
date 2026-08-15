# STEP 2b — Part A: Task 1 co_active refinement — **DONE**

target_analysis sub-task. Sources READ-ONLY; outputs only in `target_analysis/`. No enrichment value
recomputed differently — this STEP **relabels** (`concordant`→`co_active`) and adds **backend-rigor
columns**, then **independently recomputes every FINAL_STATS Task-1 anchor and asserts equality**
(write-gated on all asserts). Build: `06_target_enrichment/target_enrichment.py`. Env: pandas 2.2.2 / statsmodels 0.14.5,
seed=0. All claims hypothesis-generating (association, not causation/quality).

---

## 1. Terminology: `concordant` → `co_active` (all columns/labels)
`co_active = L1_exp_active ∩ L2_exp_active` (binding + G-protein co-activity). The word "concordant" is
reserved for **MoA-direction** concordance (Part B) and is **not** used in Part A — they are different
concepts. Renamed columns: `n_co_active_exp`, `n_co_active_aiincl`, plus `enrichment_status_*`.

## 2. High-activity / low-approval targets rescued (not discarded)
`enough_n & n_approved==0` are **no longer dropped as "undefined"** — they get a one-sided Wilson upper
bound (`wilson_upper_oneside_exp`) and `enrichment_status = "depleted (upper-bounded)"`. High-n examples
flagged `high_n_no_clinical_exp=True` (`n_co_active ≥ 50 & 0 approved`) — these are **"heavily assayed
but not clinically entered"**, a genuine web/DB data feature, NOT "sparse":

| GPCR | n_co_active_exp | n_approved | status |
|---|---:|---:|---|
| MCHR1 | 760 | 0 | depleted (upper-bounded), high-n |
| ADORA3 | 601 | 0 | depleted (upper-bounded), high-n |
| GCGR | 354 | 0 | depleted (upper-bounded), high-n |

## 3. Headline enrichment = pooled, apples-to-apples (LOCKED)
**L1-active 0.86% → co_active 1.40% = 1.63×** (pooled conditional approval over all co_active pairs;
sum co_active pairs 20,145, approved-co_active pairs 282; L1-active pool rate 0.86%). The median-based
"~2× / doubling" is **secondary only** — the headline is the pooled **1.63×**.

## 4. Backend-rigor columns (NOT front-facing) + summary counts
- `ci_robust_exp_backend` = `wilson_lo_exp > exp pooled baseline (1.40%)`; `ci_robust_aiincl_backend`
  = `wilson_lo_aiincl > aiincl pooled baseline (1.11%)`.
- disease cells: `single_approved_flag_backend` = `n_approved_cat==1`; `ci_robust_backend` =
  `wilson_lo > per-disease pooled baseline`.
- **CI-robust counts: target 28 / aiincl 31 / disease 236.**

## 5. AI-inclusion narrative (secondary)
Direction **preserved in 86%** of dual-computable targets (73/85), magnitude **attenuated ~3×**
(median `log2` exp **0.564** → aiincl **0.178**). AI inclusion keeps the enrichment *direction* while
diluting its size (it adds co_active compounds, lowering both rate and baseline).

## 6. Unique-compound vs pair-level approved
Unique approved compounds among co_active (across the 211 universe targets) = **161**, distinct from the
**pair-level 282** approved-co_active instances (a drug counted once per target it is co_active at).

---

## Anchor asserts vs FINAL_STATS (write-gated) — 12/12 PASS
| anchor | expected | got |
|---|---|---|
| universe | 211 | **211** ✓ |
| sum co_active pairs | 20,145 | 20,145 ✓ |
| sum approved-co_active pairs | 282 | 282 ✓ |
| baseline pooled | 1.40% | 1.40% ✓ |
| L1-active pooled | 0.86% | 0.86% ✓ |
| pooled enrichment | 1.63× | 1.63× ✓ |
| enough_n (≥5) | 126 | 126 ✓ |
| computable log2 | 85 | 85 ✓ |
| n≥5 & 0-approved | 41 | 41 ✓ |
| point-enriched (log2>0) | 50 | 50 ✓ |
| CI-robust target | 28 | 28 ✓ |
| CI-robust aiincl | (see note) | **31** ✓ |
| CI-robust disease | 236 | 236 ✓ |

### ⚠ FINAL_STATS reconciliation needed — `AI-incl CI-robust`
FINAL_STATS §H lists **AI-incl = 20**, which I could only reproduce as the *intersection*
(exp-CI-robust ∩ aiincl-CI-robust, both measured against the **shared** 1.40% exp baseline). Per your
decision, the column uses the statistically-clean **aiincl-vs-own-baseline** definition (`wilson_lo_aiincl
> aiincl pooled 1.11%`) = **31**. → **Recommend updating FINAL_STATS §H "AI-incl 20" → "31
(aiincl-vs-own-baseline)"**, or annotating that the 20 was the intersection-vs-shared-baseline construct.
The other 11 anchors are unaffected. No headline number changes (this is a backend rigor count only).

---

## Outputs (overwrote STEP-2 versions with co_active-named, enriched-column files)
- `target_clinical_enrichment.tsv` — **211 rows × 29 cols** (co_active naming; `enrichment_status_*`,
  `wilson_upper_oneside_exp`, `high_n_no_clinical_exp`, `ci_robust_*_backend`, `n_approved_exp_drugind`).
- `target_disease_enrichment.tsv` — **1,512 rows × 13 cols** (`n_co_active_exp`, `single_approved_flag_backend`,
  `ci_robust_backend`).

## STALE / banned scan
Computed scalars {211 · 20,145 · 282 · 161 · 28 · 31 · 236} ∩ §1.5 STALE = NONE. No `3,843` (Task-1
unrelated), no 0.50/0.85, no corrected-4.93. Banned-term scan over both TSVs = CLEAN (`concordant` no
longer present as a column; reserved for Part B). Stopping per the one-STEP rule.
