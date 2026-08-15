# STEP 5e-3c — repurposing final: web-exclusion flag + Demo-C swap — **DONE**

target_analysis sub-task. Sources READ-ONLY; **no FINAL_STATS / source edited.** Only
`repurposing_candidates_final.tsv` overwritten (column added + demo flags re-set; **no rows deleted**).
Build: inline script (pandas 2.2.2, Python 3.9, seed=42). File now **544 rows × 19 cols**.

## ① FIX 1 — `web_display` (non-anchored compounds hidden from web)
- `internal_flag == no_evidenced_family` → `web_display = 'hidden'`; everything else `'visible'`. No rows removed.
- **Distribution:** visible **458 cells / 146 drugs** · hidden **86 cells / 33 drugs**.
- ⚠ **Drug-count deviation from the prompt:** the prompt expected 32 hidden drugs (→ 147 visible). The actual
  `no_evidenced_family` set is **33 drugs** (86 cells exact), so **visible = 179 − 33 = 146 drugs** (not 147).
  Reported, not forced to 147.
- **Guard passed:** 0 demo rows are hidden (all demos have evidenced families).
- **Hidden drugs (33):** the non-therapeutic chemicals the prompt cited are present — Toluene, Sodium acetate,
  Propylparaben, Riboflavin, Arbutin — plus Bufexamac, Pentetic acid, Quinacrine, Vidarabine, Oxibendazole,
  Sepiapterin, and 7 unnamed CHEMBL IDs.
- ⚠ **Caveat (flagged for 5e-4):** the `no_evidenced_family` filter is an *"unanchored-prediction"* proxy, not
  strictly "non-therapeutic." It **also hides real drugs** whose only issue is the absence of an experimental
  GPCR-active record — e.g. **Olmesartan, Nicorandil, Levobetaxolol, Levocabastine, Mequitazine, Miglitol,
  Stiripentol, Amifampridine, Dihydroergotoxine**. Hiding them is still defensible (their predictions are pure
  de-novo with no experimental anchor → lowest reliability), but if the web wants to surface anchored-elsewhere
  real drugs, a name-based non-therapeutic blocklist would be narrower than this anchor-based filter. No action
  taken — exactly the prompt's `no_evidenced_family` rule applied.

## ② FIX 2 — Demo C swapped: Ozanimod→CCR7 ⇒ Azilsartan→MRGPRX1
- Ozanimod→CCR7: 2 rows demoted (`is_demo=False`, `demo_label=''`).
- Azilsartan→MRGPRX1: 2 rows promoted (`is_demo=True`, `demo_label='TypeC_biased_signalling'`).
- Verified values: **L2 agonist 0.9999 / L3 antagonist 0.9965** (both very-high — fixes Ozanimod's weak L2 0.7153),
  evidenced = Angiotensin receptors, small_molecule, `repurposing_class = TYPE_NOVEL_INDICATION;TYPE_BIASED_SIGNALLING`
  (cross-domain: angiotensin→Mas-related, + biased). Azilsartan has only the MRGPRX1 pair predicted (clean).
- **Mail-ready sentence (Demo C):** "Azilsartan, approved for hypertension, carries a high-confidence agonist
  prediction at MRGPRX1 (Mas-related GPCR; DB family label 'Class A orphans') via the G-protein pathway, while
  the β-arrestin pathway shows antagonist activity — a biased-signalling hypothesis outside its
  experimentally-evidenced angiotensin-receptor family. MRGPRX1 is a sensory-neuron/mast-cell receptor
  implicated in itch, pain, and pseudo-allergic drug reactions, making this functional-selectivity profile
  relevant to sensory pharmacology. [AI-predicted hypothesis, not experimentally validated]"

## ③ Final verification (PART 3)
1. total rows **544** (no deletion) ✓
2. `web_display`: visible **458** / hidden **86** ✓
3. `is_demo=True` = **6** rows (TypeA×2 / TypeB×2 / TypeC×2) ✓; Ozanimod no longer a demo ✓
4. TypeC Azilsartan→MRGPRX1 L2 agonist 0.9999 / L3 antagonist 0.9965 ✓
5. columns **19** (18 + `web_display`) ✓
6. STALE scan: NONE.

## ④ Final demo trail (3)
| demo_label | drug | target | repurposing_class | L2_moa | L2_conf | L3_moa | L3_conf | is_clean |
|---|---|---|---|---|---|---|---|:--:|
| TypeA_polypharmacology | Cariprazine | HTR1A | TYPE_POLYPHARMACOLOGY | agonist | 0.9997 | agonist | 0.9986 | True |
| TypeB_novel_indication | Fingolimod | HCAR3 | TYPE_NOVEL_INDICATION | agonist | 0.9999 | agonist | 0.9780 | True |
| TypeC_biased_signalling | Azilsartan | MRGPRX1 | TYPE_NOVEL_INDICATION;TYPE_BIASED_SIGNALLING | agonist | 0.9999 | antagonist | 0.9965 | False* |

*Demo C `is_clean=False` by design — the L2/L3 MoA discordance IS the biased-signalling feature.
Mail sentences: A & B as locked in 5e-3b; C as above. All carry `[AI-predicted hypothesis, not experimentally validated]`.

## ⑤ Web-deployment checklist
- [x] **Web-facing 11 cols:** inchikey, drug_name, known_indication, predicted_uniprot, predicted_gpcr_name,
  predicted_gpcr_family, evidenced_families, predicted_moa, confidence_raw, layer, repurposing_class.
- [x] **Apply `web_display='visible'` filter** at load → 458 cells / 146 drugs surfaced; 86/33 hidden.
- [x] **Backend-only (hide):** compound_modality, internal_flag, prediction_reliability_note, target_note,
  is_demo, demo_label, web_display.
- [x] **`repurposing_class` human-readable badges:** POLYPHARMACOLOGY → "Polypharmacology (same therapeutic area)";
  NOVEL_INDICATION → "Novel-indication hypothesis"; BIASED_SIGNALLING → "Biased-signalling hypothesis";
  PEPTIDE_UNCERTAIN → caveat icon only.
- [x] **`known_indication`** still multi-category — display a single primary indication; raw list for filtering.
- [x] Cariprazine remains the recommended interactive web demo (coherent visible 5-HT prediction fan).

## ⑥ STALE scan
New scalars {544 · 458 · 86 · 146 · 33 · 6 · 19 · 0.9999 · 0.9965} ∩ §1.5 STALE = NONE. No banned terms.
(Substring "0.50"/"0.85" only inside confidence decimals — not the operating point.)

## ⑦ Next STEP (5e-4) prep
- Editor mail: 179 repurposing candidates (146 web-surfaced; the 33 unanchored hidden from the public tab but
  retained in the DB) + the three demo sentences (Demo C = **Azilsartan biased-signalling**, very-high both layers).
  Plus the 5e-2 mail revisions (171,821→59,888; 842→496; "over a thousand"→~1,180).
- **Decision for 5e-4:** confirm the `no_evidenced_family` hide-filter (which also hides real-but-unanchored
  drugs like Olmesartan/Nicorandil) — keep anchor-based, or switch to a narrower non-therapeutic blocklist.
- Carried flags (unchanged): §E `3,719` funnel and §F `171/13`-vs-`170/14` ambiguities; paper-methods artifact-removal summary.

File written: `repurposing_candidates_final.tsv` (544 × 19). Stopping per the one-STEP rule.
