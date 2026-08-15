# STEP 5e-2 — population-bias stat recompute + FINAL_STATS delta — **DONE (report only)**

target_analysis sub-task. Sources READ-ONLY. **No FINAL_STATS / mail / source edited** — recompute + delta
report only; oracle edits are a separate user-approved sub-step. Build: `step5e2.py (build script not released)` (+ verification
queries). Env: pandas 2.2.2, numpy 1.26.4, Python 3.9, seed=42. Frozen inputs: `db_load_tier_retained.tsv`
(104,645 cells / 59,889 compounds / 151 targets), tiers T1_anchored / T2_within_family / T2_focused_crossfam.

## ① PART 1 — tier_candidates/ cleanup
Deleted the 6 leftovers from 5e-1 (all existed): `tier_assignment_rev3_optA.tsv`, `T1_candidates_rev3.tsv`,
`T2_candidates_rev3.tsv`, `hub_compounds_detail.tsv`, `suspicious_artifact_screens_compounds.tsv`,
`removed_compounds_evidenced_empty_rev.tsv`. **Remaining (3 canonical only):** `tier_assignment_rev3.tsv`,
`target_hub_diagnostic.tsv`, `removed_pairs_reasons_rev3.tsv`.

## ② PART 2 — recomputed stats (source · logic · n)

**[2-1] AI-completion coverage** — source `db_load_tier_retained.tsv` (+ `Compound_Chain_Catalog_v3` max_phase).
- (a) RETAINED unique compounds = **59,889** (vs §E 169,635 pre-filter).
- (b) per-tier: T1_anchored **47,698** · T2_within_family **9,291** · T2_focused_crossfam **9,457**.
- (d) approved (max_phase==4) with ≥1 RETAINED AI-active = **496** (vs §E 841). **Verified real, not a pipeline
  error:** approved ∩ *full pre-filter universe* = **841** (reproduces oracle exactly) → the −345 is the
  genuine effect of removing approved drugs' de-novo / hub / M-target cross-family cells.
- (e) AI-enabled synergy shortlist (pre-filter) recomputes to **1,263** (reproduces §D/§E exactly). Post-filter,
  by **precise enabling-cell tier** (the synergy-enabling L2/L3 cell at the L1+L4-backbone pair): **1,179** have
  that cell as T1_anchored, **70 lost** (all enabling cells REMOVED), 14 in T2 (a `L3L4_Pair_Catalog` vs
  `Bioassay_MoA_Master` active-label reconciliation edge). → **AI-enabled synergy post-filter = 1,179** (−84).

**[2-2] Prediction-GPCR coverage / skew** — source `db_load_tier_retained.tsv`, PAIR-level (unique compound per target).
- (a) DB-loaded prediction GPCRs (unique UniProt in RETAINED) = **151** (vs §F 154; GLP2R/CALCR/VIPR2 → 0).
- (b/c) **L2:** 148 covered · top-10 share **37.3%** · Gini **0.644** (vs §F 70.6% / 0.835).
  **L3:** 137 covered · top-10 share **45.8%** · Gini **0.699** (vs §F 83.9% / 0.909).
  *Basis check:* pre-filter pair-level Gini = L2 0.823 / L3 0.905 with top-10 70.6% / 83.9% — matches the oracle
  §F (0.835/0.909; tiny gap = 167/152-covered broader scope) → the retained values are a true apples-to-apples
  flattening from hub/modality removal.

**[2-3] Synergy pool** — SI panel is experimental → AI-filter-invariant.
- experimental base = 2,580 (unchanged). pool pre-filter = 2,580+1,263 = **3,843** (reproduces §D). pool
  post-filter = 2,580+**1,179** = **3,759** (−84).
- **SI panel L2 9.71 / L3 13.01 / full-chain 13.45 = EXPERIMENTAL → UNCHANGED** (confirmed; not recomputed; no
  STALE recovery/SI/null tokens touched).

**[2-4] Funnel.**
- OLD: 277,544 → 169,635 → 3,719 → 1,263.
- NEW (binder-gate retained, reproducible path): **277,544 → 59,889 (RETAINED AI coverage) → 3,280 (RETAINED ∩
  exp L1+L4 backbone) → 1,179 (AI-enabled synergy, T1-anchored)**.
- ⚠ **`3,719` not reproducible** (see caveat ⑤) — NEW backbone step uses the reproducible 3,280; the funnel
  middle number is held for clarification, not auto-substituted.

## ③ PART 3 — FINAL_STATS delta table (`stat_delta_summary.tsv`)
| section | stat | current_oracle | recomputed | delta | action |
|---|---|---|---|---|---|
| §E | AI-completion compounds | 169,635 | **59,889** | −109,746 | UPDATE |
| §E | repurposing approved (600-scope) | 841 | **496** | −345 | UPDATE |
| §E | AI-enabled synergy shortlist | 1,263 | **1,179** | −84 | UPDATE |
| §E | funnel AI-coverage step | 169,635 | **59,889** | −109,746 | UPDATE |
| §E | funnel exp L1+L4 backbone step | 3,719 | 3,280* | flag | **HOLD-clarify** |
| §F | loaded prediction GPCRs | 154 | **151** | −3 | UPDATE |
| §F | L2 Gini (pair-level) | 0.835 | **0.644** | −0.191 | UPDATE (internal) |
| §F | L3 Gini (pair-level) | 0.909 | **0.699** | −0.210 | UPDATE (internal) |
| §F | L2 top-10 share | 70.6% | **37.3%** | −33.3pp | UPDATE (internal) |
| §F | L3 top-10 share | 83.9% | **45.8%** | −38.1pp | UPDATE (internal) |
| §D | synergy pool total | 3,843 | **3,759** | −84 | UPDATE |
| §D | AI-enabled synergy | 1,263 | **1,179** | −84 | UPDATE |
| §D | SI panel 9.71/13.01/13.45 | — | unchanged | 0 | NO ACTION |

### Draft FINAL_STATS wording (for the approval sub-step — NOT yet applied)
- §E: "AI-completion compounds (≥1 binder-gated, **population-bias-filtered (T1/T2)** high-conf L2/L3 fill) =
  **59,889**" (note prior 169,635 = pre-filter universe). "repurposing approved (600-scope) = **496**."
  "AI-enabled shortlist = **1,179** (= synergy pool AI-enabled, enabling cell T1-anchored)."
- §F: "prediction universe served = **151** loaded GPCRs (post-filter); skew after population-bias filter:
  L2 top-10 **37.3%**, Gini **0.644** (148 covered); L3 top-10 **45.8%**, Gini **0.699** (137 covered)."
  Keep the pre-filter 0.835/0.909 as the *raw* skew for context.
- §D: "synergy pool: 2,580 → **3,759** (binder-gated, **population-bias-filtered**); AI-enabled = **1,179**
  (70 removed by the filter)." SI panel line unchanged.

## ④ PART 4 — editor-mail revision points (identify only; mail edit = 5e-4)
No standalone mail file exists in `target_analysis/` scope; the figures below also recur in the webexport STEP
reports (different sub-task — not edited here). Revision points by the prompt's phrases → new values:
- "171,821 unique compounds" → **59,889** (population-bias-filtered AI-completion coverage). ⚠ large change.
- "842 approved drugs" → **496** (approved with ≥1 retained AI-active; the 600-scope pre-filter count 841/842 stays as the raw coverage figure).
- "over a thousand candidate compounds" → "**approximately 1,180 candidate compounds**" (AI-enabled synergy 1,179).
- General: any sentence quoting pre-filter AI-prediction breadth (e.g. per-target coverage / "predicted across N
  GPCRs") should carry the filtered framing; Gini / per-target skew are internal QC and absent from the mail (no check needed).

## ⑤ Caveats
- **`3,719` (funnel backbone step) not reproducible** from current sources: backbone (L1&L4 exp active) ∩
  binder-gate universe = **3,424**; ∩ conf-only universe = **3,571**; ∩ retained = 3,280 — none equal 3,719.
  This is a **pre-existing oracle-definition ambiguity** (which AI-coverage universe / backbone source the
  original funnel used), **not introduced by the population-bias filter**. Per the no-guessing rule I did NOT
  substitute it; the funnel step is marked **HOLD-clarify**. Everything else in §E/§F/§D reproduced its
  pre-filter oracle value exactly (841, 1,263, 3,843, top-10 70.6%/83.9%), so the recompute pipeline is sound.
- Synergy post-filter has two readings: **1,179** (precise — enabling cell is T1_anchored) vs 1,219
  (compound appears anywhere in RETAINED). I headline the rigorous **1,179**; the 14-cell gap is a pair-catalog
  vs MoA-master active-label reconciliation, negligible.
- approved drop 841→496 and coverage drop 169,635→59,889 are the **intended, verified** consequence of the
  filter (de-novo / hub / M-target cells removed); both pre-filter values reproduce the oracle exactly.

## ⑥ Next STEP (5e-3) proposal
With deltas quantified, 5e-3 should: (i) get user decision on the **3,719 funnel-step** definition (re-derive
vs retire the funnel middle number) before any §E edit; (ii) on approval, apply the §E/§F/§D UPDATEs to
FINAL_STATS with the draft wording above (lineage-preserving: keep pre-filter values labeled "raw"); (iii)
hand the §F internal Gini and the mail revision list to 5e-4 (mail) — no mail edit until FINAL_STATS is locked.

## STALE / scope scan
Recomputed scalars {59,889 · 496 · 1,179 · 151 · 0.644 · 0.699 · 3,759 · 3,280} ∩ §1.5 STALE = NONE. No
0.50/0.85, no 11,428, no recovery 82.5/92.3, no corrected-4.93, no banned tokens. SI panel (experimental)
confirmed unchanged. **No FINAL_STATS / mail / source modified.** Stopping per the one-STEP rule.
