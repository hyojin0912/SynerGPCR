# STEP 5e-3 — Task 3 repurposing rebuilt on the population-bias-filtered tier — **DONE**

target_analysis sub-task. Sources READ-ONLY. Writes ONLY to `target_analysis/`. Build: `step5e3.py (build script not released)`
(pandas 2.2.2, numpy 1.26.4, Python 3.9, seed=42). Definition (frozen): a repurposing candidate = an
**approved drug** with a **credible AI prediction in a GPCR family outside its experimentally-evidenced
families** = `approved ∩ tier == T2_focused_crossfam`. Source: `db_load_tier_retained.tsv` (post 5e-2 fix,
104,641 cells) + `Compound_Chain_Catalog_v3` (max_phase==4) + `Drug_Indication_Master` → `UMLS_CUI_to_Category`
+ `target_lookup` (gene/family). Legacy `repurposing_candidates.tsv` (pre-filter) renamed to
`repurposing_candidates_legacy.tsv` (preserved, not deleted); `repurposing_candidates_v2.tsv` is canonical.

## ① Total — approved ∩ T2_focused_crossfam
- **544 cross-family hypothesis cells / 179 unique approved drugs.**
- (context: T2_focused_crossfam overall = 18,417 cells / 9,457 compounds; the approved subset is 179 drugs.)

## ② Layer breakdown
| layer | cells | unique approved drugs |
|---|---:|---:|
| L2 (G-protein) | 456 | 179 |
| L3 (β-arrestin) | 88 | 63 |
| **total** | **544** | **179** (all 179 have ≥1 L2; 63 also carry L3) |

## ③ PART 4 — §E note added (FINAL_STATS)
Appended under the §E repurposing line (does not alter the 496 headline; it is a sub-classification):
> *of which 179 (approved ∩ T2_focused_crossfam, 544 cells) carry novel cross-family repurposing predictions
> [→ repurposing_candidates_v2.tsv; the editor-mail repurposing count]*

→ **179** is the figure the editor mail should use for "repurposing candidates" (replacing any pre-filter count).

## ④ Demo trail (3, named approved · clean cross-family · clinically-relevant target · high conf)
All carry the mandatory label **[AI-predicted hypothesis, not experimentally validated]**.

**1. Cariprazine — Dopamine → Serotonin (HTR1A, agonist, conf 0.9997, L2)**
- Approved D2/D3 partial-agonist antipsychotic (schizophrenia, bipolar). Evidenced family: Dopamine receptors.
  AI predicts agonist activity at **5-HT1A (HTR1A, P08908)** — a different family (5-hydroxytryptamine).
- *Why interesting:* crosses from its known dopaminergic class into serotonergic signalling; **5-HT1A partial
  agonism is in fact a documented part of cariprazine's real pharmacology** (Kiss et al., 2010) — so the model
  recovers a genuine, label-relevant cross-family polypharmacology — a strong validation-style example.
- **Mail sentence:** "For example, Cariprazine, currently approved for schizophrenia and bipolar disorder,
  carries a high-confidence agonist prediction at the serotonin 5-HT1A receptor (HTR1A; 5-hydroxytryptamine
  family) — outside its experimentally-evidenced dopamine-receptor family — a receptor associated with mood
  and anxiety regulation."

**2. Atosiban — Vasopressin/Oxytocin → Angiotensin (AGTR1, agonist, conf 0.9999, L2)**
- Approved tocolytic (preterm labour); oxytocin/vasopressin-receptor antagonist (peptide). Evidenced family:
  Vasopressin and oxytocin receptors. AI predicts agonist at **angiotensin AT1 (AGTR1, P30556)**.
- *Why interesting:* bridges two peptide-hormone GPCR systems (oxytocin/vasopressin ↔ angiotensin) that jointly
  govern **cardiovascular and fluid homeostasis** — a distinct, non-aminergic family pair suggesting a
  cardiovascular / uterine-perfusion repurposing axis. No direct prior repurposing report; mechanistically plausible.
- **Mail sentence:** "Atosiban, currently approved as a tocolytic for preterm labour, carries a high-confidence
  agonist prediction at the angiotensin AT1 receptor (AGTR1; angiotensin family) — distinct from its
  experimentally-evidenced vasopressin/oxytocin-receptor family — a target central to cardiovascular and fluid
  homeostasis."

**3. Pentoxyverine — Muscarinic → Opioid (OPRK1, agonist, conf 0.9997, L2)**
- Approved antitussive (cough). Evidenced families: muscarinic acetylcholine + glycoprotein hormone. AI
  predicts agonist at the **κ-opioid receptor (OPRK1, P41145)**.
- *Why interesting:* cough suppression is classically opioid-mediated (codeine, dextromethorphan); a κ-opioid
  agonist prediction proposes a concrete receptor-level mechanism for pentoxyverine's antitussive action
  outside its muscarinic evidence — a within-therapeutic-area mechanistic hypothesis with literature precedent
  for opioid involvement in the cough reflex.
- **Mail sentence:** "Pentoxyverine, an approved antitussive, carries a high-confidence agonist prediction at
  the κ-opioid receptor (OPRK1; opioid family) — outside its experimentally-evidenced muscarinic-acetylcholine
  family — consistent with the established role of opioid receptors in cough suppression."

*Additional notable hit (not in the trail — direction caveat):* **Fentanyl → ADRA2B** (α2B-adrenoceptor,
predicted **antagonist**, conf 0.9996, L2). Interesting opioid→adrenergic-analgesia adjacency, but the
antagonist direction runs opposite to α2-agonist analgesia, so it is weaker as a clean demo; logged for review.

## ③' PART 3 — files
- `repurposing_candidates_v2.tsv` (**canonical**, 544 rows / 179 drugs): `inchikey | drug_name |
  known_indication | predicted_uniprot | predicted_gpcr_name | predicted_gpcr_family | evidenced_families |
  predicted_moa | confidence_raw | layer | tier_label` (tier_label = `AI_repurposing_hypothesis`, fixed),
  sorted by confidence_raw desc.
- `repurposing_candidates_legacy.tsv` (renamed from the pre-filter `repurposing_candidates.tsv`, preserved for reference).

## ⑤ STALE / scope scan
New scalars {544 · 179 · 456 · 88 · 63} ∩ §1.5 STALE = NONE. The §E note (179/544) introduces no stale token.
⚠ A substring grep for "0.50"/"0.85" matches **confidence_raw decimal values** (e.g. 0.8594, 0.5013) in the
v2 TSV — these are legitimate model outputs, **not** the stale 0.50/0.85 operating point (no operating-point
usage present). No banned terms (DSAV/ACTR/forward/backward/permutation). FINAL_STATS edit limited to the one
§E sub-note.

## ⑥ Next STEP (5e-4) prep
- Editor mail: use **179** repurposing candidates; drop in the three demo sentences above (each with the
  [AI-predicted hypothesis] label). Apply the 5e-2 mail revisions (171,821→59,888; 842→496; "over a thousand"→
  ~1,180) in the same pass.
- `known_indication` in v2 is multi-category and noisy for some drugs (Drug_Indication maps many UMLS CUIs);
  the demo sentences use the clinically-correct primary indication, not the raw category list — flag for the
  web "repurposing" tab whether to show primary-indication-only.
- Still pending (from earlier flags): the §E `3,719` and §F `171/13`-vs-`170/14` ambiguities; the paper-methods
  artifact-removal summary.

Files written: `repurposing_candidates_v2.tsv`, `repurposing_candidates_legacy.tsv` (rename), `FINAL_STATS.md`
(§E sub-note). Stopping per the one-STEP rule.
