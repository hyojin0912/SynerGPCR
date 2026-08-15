# STEP 5e-2 fix — non-binder edge-band removal + FINAL_STATS oracle update — **DONE**

target_analysis sub-task. Sources READ-ONLY except the three files this STEP is authorized to edit:
`db_load_tier_retained.tsv`, `db_load_tier_final.tsv`, `FINAL_STATS.md`. Build: ad-hoc verified scripts
(pandas 2.2.2, numpy 1.26.4, Python 3.9). Two fixes only; nothing else touched.

## ① FIX 1 — 4-row non-binder edge-band removed
The 4 `predicted_moa == 'non-binder'` rows (all `binding_prob = 0.5`, the gate boundary; semantically
contradictory — a binder-gated AI-active cell with no MoA) were removed from both canonical files.

| inchikey | uniprot (gene) | layer | tier (pre-fix) | binding_prob | confidence_raw |
|---|---|---|---|---:|---:|
| MTYIMKQNRSISNK-UHFFFAOYSA-N | P34972 (CNR2) | L2 | T1_anchored | 0.5 | 0.5000 |
| PAVNNTWFJALROS-UHFFFAOYSA-N | P29275 (ADORA2B) | L2 | T1_anchored | 0.5 | 0.5000 |
| SBFAOBRHMFHNPX-RSVOXEIVSA-N | P0DMS8 (ADORA3) | L2 | T1_anchored | 0.5 | 0.4998 |
| XKRDCKIVPOWNFM-QSJAHPCTSA-N | P30542 (ADORA1) | L2 | T2_within_family | 0.5 | 0.5000 |

**Post-fix distributions** (both files; final drops the same 4 RETAINED cells, REMOVED untouched):
| tier | pre-fix | post-fix | Δ |
|---|---:|---:|---:|
| T1_anchored | 69,699 | **69,696** | −3 |
| T2_within_family | 16,529 | **16,528** | −1 |
| T2_focused_crossfam | 18,417 | **18,417** | 0 |
| **RETAINED total** | 104,645 | **104,641** | −4 |
| REMOVED | 262,156 | 262,156 | 0 |
| **final total** | 366,801 | **366,797** | −4 |

Layer (retained): L2 70,765 → **70,761** (−4); L3 **33,880** (unchanged). Remaining non-binder in retained = **0**.

⚠ **One headline value shifted because of the fix — reported, not assumed.** `PAVNNTWFJALROS` had its *only*
retained cell be a non-binder cell → it dropped from the compound set entirely. So **AI-completion compounds =
59,889 → 59,888** (−1). Per the no-guessing rule I recomputed every §E/§F/§D headline on the post-fix files and
wrote the actual post-fix value (**59,888**, not the pre-fix 59,889 from the prompt). All other headlines held
exactly: repurposing approved 496 · synergy 1,179 · backbone 3,280 · loaded GPCRs 151 · Gini L2 0.644 / L3 0.699.

## ② FIX 2 — FINAL_STATS.md before/after
| section | item | before | after |
|---|---|---|---|
| §D | synergy pool total | 2,580 → 3,843 | **2,580 → 3,759** (binder-gated, population-bias-filtered) |
| §D | AI-enabled synergy | 1,263 | **1,179** (70 removed by filter) |
| §D | SI panel 9.71/13.01/13.45 | — | **unchanged** (experimental; NO ACTION) |
| §E | AI-completion compounds | 169,635 | **59,888** (pre-filter universe 169,635) |
| §E | funnel | 277,544→169,635→3,719→1,263 | **277,544→59,888→3,280→1,179** |
| §E | repurposing approved | 841 | **496** (pre-filter 841; modality-unanchored removal) |
| §E | AI-enabled shortlist | 1,263 | **1,179** (pre-filter 1,263) |
| §E | coverage≠candidates note | 169,635 / 1,263 | **59,888 / 1,179** (credible, population-bias-filtered) |
| §F | loaded prediction GPCRs | 154 | **151** (served post-filter) |
| §F | L2 skew | 70.6% / 0.835 | **37.3% / 0.644** (148 covered; pre-filter retained as raw) |
| §F | L3 skew | 83.9% / 0.909 | **45.8% / 0.699** (137 covered; pre-filter retained as raw) |

- §E note "coverage ≠ candidates" updated in place (no duplicate added) to 59,888 / 1,179.
- §E funnel: **3,719 replaced by the reproducible post-filter 3,280**, with a one-line pre/post-basis comment;
  the old `3,719→1,263` decomposition kept and relabeled "(pre-filter lineage)".
- §D guardrail header (lines 11–13) **left untouched** — it documents why the *pre-filter* 3,843 signature is
  valid (vs the 11,428/3,843 STALE pair); 3,843 still appears there only as pre-filter lineage. Flagging in case
  you want that header's "2,580 → 3,843" example refreshed to "→ 3,759" in a later pass (out of this fix's scope).

## ③ prediction-only GPCR verification (§F)
Recomputed binder-gate prediction accounting (source: `synergpcr_annotation_{L2,L3}.csv` ∩ 600-scope):
- prediction universe (any scope) = **170**; in-600-scope = 156; **prediction-only outside-600 = 14**.
- loaded (db universe, 600-scope) = 154 → **served post-filter = 151**.
- **The 3 GPCRs dropping 154→151 are CALCR (P30988), PTAFR (P25105), SMO (Q99835)** — **NOT** the prompt's
  assumed GLP2R/CALCR/VIPR2. Correction: only **CALCR** (Class-B1, M-target, modality-unanchored → 0 defined-
  anchored cells) matches; GLP2R (O95838) and VIPR2 (P41587) were **never in the loaded universe** (0 AI-active
  cells), so they could not "drop." PTAFR and SMO are non-M targets whose few cells were all de-novo/hub-removed.
- ⚠ **prediction-only outside-600 is filter-invariant** (confirmed) but binder-gate recompute gives **14, and
  universe 170** — ±1 vs the current oracle 171/13. This is a **pre-existing ambiguity** (the conf-only line
  already says 14), **not caused by the filter**. I wrote 151 served + flagged 170/14 in §F rather than
  fabricating a clean "171 = 151 + 13" decomposition (which does not close arithmetically).

## ④ STALE scan (full FINAL_STATS.md, post-edit)
No newly-introduced STALE/banned value. The only matches are legitimate guardrail/lineage context: header
lines 11–13 (the STALE-token *definition*), line 32 ("no corrected-4.93" negation), line 81 ("RETIRED legacy
31.9%/14.2%, BACC 0.90/0.85"). New headline scalars {59,888 · 496 · 1,179 · 3,759 · 3,280 · 151 · 0.644 ·
0.699} ∩ STALE = NONE. SI panel (9.71/13.01/13.45) and recovery (79.7/79.6) untouched.

## ⑤ Next STEP (5e-3) prep
- Open item: the **§E funnel `3,719`** pre-existing ambiguity and the **§F `171/13` vs `170/14`** ±1 ambiguity —
  both are flagged, filter-invariant, and need a definition decision (re-derive vs retire) if exactness matters
  for the manuscript. Not blocking the population-bias update.
- The **3 served-drop GPCR identities (CALCR/PTAFR/SMO)** correct the prompt's assumption — confirm before any
  mail/methods text quotes them.
- Mail revision points (from 5e-2): 171,821→**59,888**, 842→**496**, "over a thousand"→**~1,180** — apply in 5e-4
  after FINAL_STATS is considered locked. No mail edited here.
- Paper-methods artifact-removal summary (target exclusion → modality-anchoring → pan-family-hub K=5) still pending.

Python 3.9 / PEP8. **Files written:** `db_load_tier_retained.tsv` (104,641), `db_load_tier_final.tsv`
(366,797), `FINAL_STATS.md` (§D/§E/§F). Stopping per the one-STEP rule.
