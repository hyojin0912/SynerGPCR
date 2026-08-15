# STEP 5e-5c — stale-file regeneration + `final_export/` assembly — **PARTIAL (2 STOPs, 1 HOLD)**

target_analysis sub-task. Sources **READ-ONLY**; **no source edited** (PART 7 CLAUDE.md edit **HELD** — see ⑦).
Outputs only under `final_export/`. Build: inline scripts (pandas 2.2.2, numpy 1.26.4, Python 3.9).
Governing files re-read: `target_analysis_CLAUDE.md` (root) + `FINAL_STATS.md` (oracle).

`BASE = …/WebExport_v2` · `OUT = …/WebExport_v2/final_export`.

**Headline:** PART 0/1/4-A/4-B/5/6 completed and verified against the oracle. **PART 2 and PART 3 hit genuine
oracle mismatches → STOPPED, no file written, no value hardcoded.** PART 7 held because it depends on the
unresolved PART 2 value. 14 files now in `final_export/` (10 KEEP + 4 regenerated); two STALE files
(`gpcract_3point_si.csv`, `candidate_shortlist.csv`) deliberately **absent** pending resolution.

---

## ① PART 1 — KEEP files copied (10/10, row-counts identical to source)
| file | rows |
|---|---:|
| synergy_stats.csv | 17 |
| single_layer_table.csv | 5 |
| compound_lookup.csv | 481,241 |
| target_lookup.csv | 601 |
| drug_indication.csv | 53,563 |
| target_indication.csv | 11,545 |
| disease_category_lookup.csv | 13 |
| target_clinical_enrichment.tsv | 212 |
| target_disease_enrichment.tsv | 1,513 |
| assay_clinical_validity.tsv | 16 |

All `cp` verified byte-for-row identical (src wc -l == dst wc -l). `figures/` created empty (later STEP).

---

## ② PART 2 — gpcract_3point_si.csv L2 row → **STOP, NOT WRITTEN** ⛔
Independent recompute of the L2 operating-point SI from `db_load_tier_retained.tsv` (L2 binder-gated
AI-active) ∩ L1&L4-active backbone (9,357 compounds):

| quantity | recompute | old 06-04 file | FINAL_STATS §D |
|---|---:|---:|---:|
| pool n | 3,127 | 3,920 | (unstated) |
| approved in pool | 164 | 198 | — |
| obs_rate | 0.052446 | 0.05051 | — |
| e_bliss | 0.006586 | 0.006586 | — |
| **SI** | **7.96** | 7.67 | **7.74** |
| recovery vs ceiling 9.71 | 82.0% | 79.0% | 79.7% |

→ My independent recompute yields a **third, distinct value (7.96)** — it reproduces **neither** the old-file
7.67 (= root master §3) **nor** the oracle §D 7.74. The three values arise from **three different L2 AI-active
pools** (3,127 vs 3,920 vs §D's unstated count). Per the non-negotiable rule, **I did not write 7.74
(hardcode-bypass forbidden) and did not write my own 7.96** — the file is left out of `final_export/` until the
canonical L2 pool definition is fixed. The `assert abs(SI-7.74)<0.05` failed legitimately (Δ=0.22).
⚠ Note recovery 82.0% is *near* but **not equal** to the STALE 82.5% — it is a recompute artifact of this pool
choice, **not** the banned token; flagged only for transparency. **This is the central blocker of the STEP.**

## ③ PART 3 — candidate_shortlist.csv (1,357→?) → **STOP, NOT WRITTEN** ⛔
The PART 3 definition (L1&L4-active backbone ∩ `T1_anchored` ∩ AI-active L2/L3) yields:

| bucket | count |
|---|---:|
| L2-only | 2,065 |
| L3-only | 188 |
| L2 & L3 | 422 |
| **TOTAL** | **2,675** |

→ **2,675 ≠ 1,179** (the §D/§E AI-enabled shortlist). The discrepancy is definitional: the **1,179** figure
(reproduced exactly in my earlier enabling-cell query: L2-only 1,036 / L3-only 7 / both 136) requires the
**enabling-cell** definition — the L2/L3 layer must be **experimentally missing** (L1&L4 active backbone with
L2/L3 NA, then AI-filled) — whereas PART 3's "any T1-anchored AI-active on backbone" is broader and **double-counts
pairs whose L2/L3 already has experimental data**. I did **not** write a 2,675-row file mislabeled as the 1,179
shortlist, and did **not** force 1,179. **Decision needed (5e-5d):** the web `candidate_shortlist` should be
built from the **enabling-cell** logic (= 1,179) — confirm and I will regenerate with that definition.

## ④ PART 4 — repurposing_view × 2 → **WRITTEN ✓**
**4-A `repurposing_view_crossfam.csv`** (AI-Prediction curated tab; `web_display=='visible'`):
- **458 cells ✓** (assert pass). **Unique drugs = 146**, **not the prompt's 179.** 179 is the *total* cross-family
  drug count (`544` cells); the 33 `no_evidenced_family` drugs (86 cells) are `web_display='hidden'` (5e-3c), so
  the visible curated view is **146 drugs / 458 cells**. Reported, not forced to 179.
- `smiles` merged from `compound_lookup` — **458/458 non-null** (full ECFP4-fallback coverage for this view).
- 14 web cols incl. `is_demo`/`demo_label` (6 demo rows preserved).

**4-B `repurposing_view_full.csv`** (Drug-detail AI-profile source; approved ∩ ≥1 retained AI-active):
- **496 drugs ✓** (assert pass) / 1,771 cells.
- ⚠ **target_lookup column adaptation (flag for 229 schema):** the prompt's `gpcr_name`/`gpcr_class` **do not
  exist**. Real columns used: `protein_name` (↔ gpcr_name), `gpcr_class_name` (↔ gpcr_class), `gpcr_family`,
  `gene_name`. The 229 target table should standardize these names.

## ⑤ PART 5 — compound_chain_summary.csv flag recompute → **WRITTEN ✓**
`gpcract_prediction_available` re-set = (`InChIKey` ∈ retained-tier compounds).
- **Unique compounds flagged True = 59,888 ✓ (== §E exactly)** — assert pass at the compound level.
- Row-sum changed **303,372 → 122,656**. ⚠ The prompt's literal `assert |row_sum − 59,888| < 5%` would **fail**
  because this file is **pair-level** (396,169 rows): summing a per-row boolean counts *pairs*, not compounds.
  I validated against the intended metric (**unique compounds = 59,888**) and wrote the file. **229 schema
  caveat:** as a pair-level column, `gpcract_prediction_available=True` here means "this compound received AI
  completion *somewhere*," not "*this* (compound,target) pair has a retained prediction." If pair-level semantics
  are wanted, it must be a pair-membership join instead — flagged for the schema decision.

## ⑥ PART 6 — gpcract_predictions.csv regenerated → **WRITTEN ✓**
Rebuilt from `db_load_tier_final.tsv` (366,797 rows) with web columns + `is_retained` + `ai_active`.
- **retained (tier≠REMOVED) = 104,641 ✓** (assert pass). `ai_active` cells = 104,641 (all retained cells are
  AI-active by construction — binder-gate+threshold is how they were retained).
- Same target-column adaptation as 4-B (`gpcr_class_name`/`gpcr_family`/`gene_name`; no `gpcr_class`/`gpcr_name`).

## ⑦ PART 7 — CLAUDE.md edit → **HELD (not executed)** ✋
The instruction edits `target_analysis_CLAUDE.md` L2 `7.67→7.74` (etc.) to adopt §D. **Held** because PART 2
shows §D's `7.74` is **not independently reproducible** (my recompute = 7.96; old file = 7.67) — three values
from three pool definitions. Editing the governing file to enshrine an unverifiable scalar would violate the
"recompute-then-assert; mismatch → STOP, no guessing" rule. **No source byte changed.** Once the canonical L2
pool is fixed (PART 2 resolution), PART 7 becomes a one-line, scoped edit. Also note: the **repo-root master
`<repo-root>/CLAUDE.md` §3** still carries the full `7.67 / null95 7.17 / recovery 79.0%` triple —
the original 5e-5b conflict — which likewise should be reconciled to whatever value wins, not before.

---

## ⑧ PART 8 — `final_export/` final listing + checklist
14 files (10 KEEP + 4 regenerated) + empty `figures/`:
```
   2,729  assay_clinical_validity.tsv        (KEEP)
27,166,391  compound_chain_summary.csv        (REGEN: flag, 396,169 rows)
52,930,473  compound_lookup.csv               (KEEP)
     735  disease_category_lookup.csv         (KEEP)
4,184,192  drug_indication.csv                (KEEP)
62,656,012  gpcract_predictions.csv           (REGEN: 366,797 rows, 104,641 retained)
  127,367  repurposing_view_crossfam.csv      (REGEN: 458 cells / 146 drugs)
  466,954  repurposing_view_full.csv          (REGEN: 1,771 cells / 496 drugs)
     364  single_layer_table.csv              (KEEP)
   1,671  synergy_stats.csv                   (KEEP)
   49,332  target_clinical_enrichment.tsv     (KEEP)
  137,319  target_disease_enrichment.tsv      (KEEP)
  732,323  target_indication.csv              (KEEP)
   44,081  target_lookup.csv                  (KEEP)
```
**Verification checklist:**
- [x] repurposing_view_crossfam: **458 cells** (drugs **146**, not 179 — see 4-A) (§E)
- [x] repurposing_view_full: **496 drugs** (§E)
- [x] compound_chain_summary flag: **59,888 unique compounds** (§E)
- [x] gpcract_predictions: **104,641 retained cells** (§E)
- [ ] candidate_shortlist: **1,179** — NOT produced (PART 3 def gave 2,675; enabling-cell def needed) ⛔
- [ ] gpcract_3point_si L2 row: **7.74** — NOT produced (recompute 7.96 ≠ §D; pool def unresolved) ⛔

## ⑨ target_lookup column-gap flags (for 229 schema design)
- No `gpcr_name` → used `protein_name` (or `uniprot_entry_name`). No `gpcr_class` → used `gpcr_class_name`
  (and `gpcr_class_letter` is available). 229 target table should adopt canonical names so views join cleanly.
- `compound_chain_summary.gpcract_prediction_available` is **compound-level coverage** on a pair-level table —
  decide pair-level vs compound-level semantics for the 229 contract (see ⑤).

## ⑩ STALE scan
New scalars {458 · 146 · 496 · 1,771 · 59,888 · 104,641 · 122,656 · 366,797 · 2,675} ∩ §1.5 STALE = NONE.
The two STOP values 7.96 / 82.0% were **printed only** (never written to any file); 82.0% ≠ banned 82.5%.
No banned terms (DSAV/ACTR/forward/backward/permutation). Null always "Monte-Carlo random-assignment".

## ⑪ Next STEP (5e-5d) prep — resolve the two STOPs before manifest/Methods
1. **PART 2 (blocker):** fix the canonical L2 operating-point pool so SI is reproducible. Three candidates on
   record: 7.67 (pool 3,920/198, old file = root §3), 7.74 (§D, pool unstated), 7.96 (pool 3,127/164, this
   recompute). Pick the definition, then write `gpcract_3point_si.csv` + execute the scoped PART 7 edit
   (target_analysis_CLAUDE.md, and decide whether repo-root master CLAUDE.md §3 is updated in lockstep).
2. **PART 3:** regenerate `candidate_shortlist.csv` from the **enabling-cell** definition (= 1,179, breakdown
   L2-only 1,036 / L3-only 7 / both 136) rather than the broad backbone∩T1 (2,675).
3. **229 schema:** confirm target-column names and the pair-vs-compound semantics of
   `gpcract_prediction_available`.
4. (Carried) §E `3,719` funnel and §F `171/13`-vs-`170/14` ambiguities; `no_evidenced_family` hide-filter
   scope; paper-methods artifact-removal summary.

Report written: `target_analysis/STEP5e5c_final_export.md`. 4 files written to `final_export/`; 2 STOPs + 1
HOLD reported. No source / FINAL_STATS edited. Stopping per the one-STEP rule.
