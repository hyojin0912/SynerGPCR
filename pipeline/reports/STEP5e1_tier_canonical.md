# STEP 5e-1 — canonical DB-load tier output + folder cleanup — **DONE**

target_analysis sub-task. Sources READ-ONLY. **No stat / FINAL_STATS / mail changed** — this STEP only
materializes the frozen rev3 (Option B) tiering into the canonical DB-load files and cleans `tier_candidates/`.
Build: `step5e1.py (build script not released)`. Env: pandas 2.2.2, numpy 1.26.4, Python 3.9, seed=42.

Canonical source: `tier_candidates/tier_assignment_rev3.tsv` (366,801 rows; frozen Option-B policy).
Outputs written to `target_analysis/` (NOT `tier_candidates/`):
`db_load_tier_final.tsv` (RETAINED + REMOVED) and `db_load_tier_retained.tsv` (DB-load, tier ≠ REMOVED).
Columns: `inchikey | uniprot | layer | tier | removed_reason | binding_prob | confidence_raw | predicted_moa |
evidenced_families | denovo_crossfam_span`. `binding_prob/confidence_raw/predicted_moa` merged from
`synergpcr_annotation_{L2,L3}.csv` (source of truth); `evidenced_families` (defined-layer active families)
and `denovo_crossfam_span` (span excl. {GLP1R,PTH1R,CRHR2}, K-frozen) re-derived per definition — no hardcoding.

## ① Row-count check (expected vs actual, tol ±5)
| bucket | expected | actual | Δ | |
|---|---:|---:|---:|:--:|
| total | 366,801 | 366,801 | +0 | OK |
| RETAINED | 104,645 | 104,645 | +0 | OK |
| T1_anchored | 69,699 | 69,699 | +0 | OK |
| T2_within_family | 16,529 | 16,529 | +0 | OK |
| T2_focused_crossfam | 18,417 | 18,417 | +0 | OK |
| REMOVED | 262,156 | 262,156 | +0 | OK |

All exact. Annotation merge: **0 rows missing** `binding_prob` (every AI-active cell has a probe record).

## ② Sanity checks (PART 3, on db_load_tier_retained.tsv) — all PASS
1. **Unique (inchikey,uniprot) pairs per tier:** T1_anchored 69,699 cells / 63,624 pairs · T2_within_family
   16,529 / 13,451 · T2_focused_crossfam 18,417 / 14,037.
2. **Layer split:** L2 **70,765** (exp ~70,765) / L3 **33,880** (exp ~33,880) — exact.
3. **M-targets in RETAINED:** 11 of 14 carry retained cells (3,971 cells) — **all T1_anchored; non-T1 rows = 0**
   → modality policy honored (no glycoprotein/Class-B1 de-novo or L4-only cell survives). The other 3 M-targets
   (GLP2R/CALCR/VIPR2) retain nothing.
4. **REMOVED rows with empty `removed_reason` = 0** (every REMOVED cell labeled
   `modality_target_unanchored` or `panfamily_hub_crossfam`).

## ③ Deleted files (`tier_candidates/`, 13 total)
Exact-match intermediates (11): `tier_assignment_summary.tsv`, `tier_assignment_rev2.tsv`,
`T1_candidates.tsv`, `T2_candidates_sample.tsv`, `T1_candidates_rev.tsv`, `T2_candidates_sample_rev.tsv`,
`T1_candidates_rev2.tsv`, `hub_compounds_Ksweep.tsv`, `suspicious_artifact_screens.tsv`,
`removed_pairs_reasons_rev.tsv`, `removed_pairs_reasons_rev2.tsv`.

**⚠ Two name-variant mappings** — the delete list named files that exist under a slightly different name; I
deleted the clearly-intended rev1/rev2 intermediate and note the mapping for transparency:
- listed `tier_assignment_rev.tsv` → actual **`tier_assignment_summary_rev.tsv`** (rev1 full assignment) → deleted.
- listed `T2_candidates_sample_rev2.tsv` → actual **`T2_candidates_rev2.tsv`** (rev2 T2 candidates) → deleted.

Listed-but-absent (skipped, not an error): `compound_artifact_categorized.tsv`,
`removed_compounds_reasons.tsv`, `per_target_table.tsv`, `proposed_excluded_targets.tsv`,
`worst_compounds_after_target_cleanup.tsv`.

## ④ Remaining files in `tier_candidates/` (9)
**Preserve (canonical, per spec):** `tier_assignment_rev3.tsv`, `target_hub_diagnostic.tsv`,
`removed_pairs_reasons_rev3.tsv`.
**Leftover — not in the 5e-1 delete list, so NOT deleted (flagged for 5e-2 decision):**
- `tier_assignment_rev3_optA.tsv` — the Option-A variant; obsolete once Option B is confirmed final.
- `T1_candidates_rev3.tsv`, `T2_candidates_rev3.tsv` — rev3 review samples (now superseded by `db_load_tier_*`).
- `hub_compounds_detail.tsv`, `suspicious_artifact_screens_compounds.tsv` — rev2 diagnostics.
- `removed_compounds_evidenced_empty_rev.tsv` — rev1 intermediate.

I did not delete these because they were absent from the explicit delete list; recommend including them in the
5e-2 cleanup sweep (all are superseded by `db_load_tier_final.tsv` / `target_hub_diagnostic.tsv`).

## ⑤ Caveats
- **4-row binder-gate edge-band** (0.004% of retained): cells at exactly `binding_prob = 0.5` (the gate's
  `≥ 0.5` boundary) where the source `predicted_moa` field — using strict `< 0.5 → non-binder` — labels them
  `non-binder` despite passing the gate (P34972, P29275, P0DMS8, P30542; all L2, confidence_raw 0.4998–0.5000).
  This is the edge-band anticipated in STEP 0 §5-2(a). I kept the model's emitted label rather than overwrite;
  retained `predicted_moa`: agonist 62,285 / antagonist 42,356 / non-binder **4**. Negligible; flag only.
- The canonical files carry **REMOVED rows too** (`db_load_tier_final.tsv`) for provenance; the DB-load set is
  `db_load_tier_retained.tsv` (104,645 cells).

## ⑥ Next STEP (5e-2) proposal
With the canonical tier files locked, 5e-2 should: (i) re-compute the affected stats **from
`db_load_tier_retained.tsv`** (AI-completion coverage, per-tier compound/pair counts, Gini, synergy ⊂ T1) and
reconcile against FINAL_STATS §E/§F — reporting deltas before any oracle edit; (ii) confirm Option B as final
(retire `tier_assignment_rev3_optA.tsv`); (iii) finish the `tier_candidates/` sweep (the 6 leftovers in ④);
(iv) draft the paper-methods artifact-removal summary (target exclusion + modality-anchoring + pan-family-hub
K=5). No FINAL_STATS / mail edit until that recompute is reviewed.

## STALE / scope scan
Scalars {366,801 · 104,645 · 69,699 · 16,529 · 18,417 · 262,156 · 70,765 · 33,880} ∩ §1.5 STALE = NONE.
No 0.50/0.85, no 11,428, no banned tokens (DSAV/ACTR/forward/backward/permutation). **No source/FINAL_STATS/
mail modified.** Stopping per the one-STEP rule.
