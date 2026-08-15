# STEP 5e-3b — repurposing validation / classification / finalization — **DONE**

target_analysis sub-task. Sources READ-ONLY; **no FINAL_STATS / source edited.** Writes only
`repurposing_candidates_final.tsv` (canonical, web-load) under `target_analysis/`; `repurposing_candidates_v2.tsv`
kept as intermediate. Build: `step5e3b.py (build script not released)` (pandas 2.2.2, numpy 1.26.4, Python 3.9, seed=42).

## ① PART 1 — input validation
- rows **544** ✓ · unique InChIKey **179** ✓ · all `tier_label == AI_repurposing_hypothesis` ✓.
- `predicted_moa ∈ {agonist, antagonist}`, **non-binder = 0** ✓. `confidence_raw` ∈ [0.4292, 1.0000].
- **evidenced_families NaN = 86 cells / 33 drugs.** Why present: T2_focused_crossfam = de-novo cross-family
  cells of *non-hub* compounds (span<5); a compound with **zero** experimental GPCR evidence is trivially
  non-hub, so its few cross-family predictions survive as "focused." These are real (binder-gated, conf-passing)
  but rest on **no experimental anchor** → tagged `internal_flag = no_evidenced_family` (kept, not deleted).

## ② PART 2 — peptide identification
- **12 cells / 3 drugs** flagged `compound_modality = peptide`: **Atosiban, Octreotide, Oxytocin** →
  `prediction_reliability_note = "peptide compound; small-molecule model prediction reliability uncertain"`.
- **Manual-review list (21 small-molecule drugs)** with a peptide-receptor *evidenced* family (Vasopressin/
  Oxytocin/Glycoprotein-hormone) but not on the peptide name-list — e.g. Formoterol, Quinidine, Pentoxyverine,
  Cyproheptadine, Mozavaptan, Phenelzine. **Reviewed → all genuine small molecules** (they simply carry an
  experimental record at a peptide-ligand receptor; not modality-mismatched) → kept `small_molecule`, no reclassification.
- **AGTR1 target_note** applied to its **31** prediction rows: *"AGTR1 endogenous ligand is a peptide
  (Angiotensin II); prediction may reflect genuine cross-peptide binding."*

## ③ PART 3 — systematic classification (`repurposing_class`, multi-label)
Domain map applied to evidenced vs predicted family (CNS / Cardiovascular / Immune / Metabolic). Same-domain →
POLYPHARMACOLOGY; different/absent-domain → NOVEL_INDICATION; per-pair L2≠L3 MoA → BIASED_SIGNALLING; peptide → PEPTIDE_UNCERTAIN.

| class | cells | drugs |
|---|---:|---:|
| TYPE_POLYPHARMACOLOGY | 274 | 95 |
| TYPE_NOVEL_INDICATION | 270 | 116 |
| TYPE_BIASED_SIGNALLING | 46 | 21 |
| TYPE_PEPTIDE_UNCERTAIN | 12 | 3 |

(Classes co-occur, so drug counts overlap.) Rows with `no_evidenced_family` (86) default to NOVEL_INDICATION
(no evidenced domain to match) **but carry the flag** — interpret as exploratory, not anchored.

**NOVEL_INDICATION, conf ≥ 0.99 — top 20** (cross-domain, highest confidence):
Tasimelteon→FFAR4 (Melatonin→free-fatty-acid, 1.000) · Sincalide→HCAR3 (1.000) · CHEMBL5896→NTSR1 (1.000) ·
Palmidrol→GPR119 (Cannabinoid→metabolic, 1.000) · Methylene-blue→GABBR1 (1.000) · **Fingolimod→HCAR3 (0.9999)** ·
Vorapaxar→GPR183 (0.9999) · Azilsartan→MRGPRX1 (0.9999) · Tauroursodeoxycholic-acid→S1PR2 (0.9999) ·
Brexanolone→PTGER2 (0.9998) · Synephrine→TAAR1 (0.9997) … (full list in file; Oxibendazole/Sepiapterin/
CHEMBL5896 carry `no_evidenced_family`).

## ④ PART 4 — demo verification
| demo | pair | repurposing_class (computed) | L2 | L3 | is_clean | conf_tier |
|---|---|---|---|---|:--:|---|
| A | Cariprazine→HTR1A | TYPE_POLYPHARMACOLOGY | agonist 0.9997 | agonist 0.9986 | **True** | very_high |
| B | Fingolimod→HCAR3 | TYPE_NOVEL_INDICATION | agonist 0.9999 | agonist 0.978 | **True** | very_high |
| C | Ozanimod→CCR7 | **TYPE_POLYPHARMACOLOGY;TYPE_BIASED_SIGNALLING** | antagonist 0.7153 | agonist 0.9999 | **False** | very_high |

**⚠ Two corrections vs the prompt's demo assumptions (data-driven):**
- **Demo C class:** the prompt labelled Ozanimod→CCR7 `TYPE_NOVEL_INDICATION`, but by the supplied domain map
  **S1P (Lysophospholipid) and CCR7 (Chemokine) are BOTH Immune_domain** → same domain → POLYPHARMACOLOGY. The
  prompt's own note ("CCR7 co-regulates lymphocyte egress — the same trafficking axis as S1P1") supports the
  same-domain reading. The **biased-signalling** label is correct (L2 antagonist vs L3 agonist). I kept it as the
  Type-C biased demo but set `repurposing_class = TYPE_POLYPHARMACOLOGY;TYPE_BIASED_SIGNALLING`. `is_clean = False`
  by design (the L2/L3 MoA discordance is the feature). ⚠ note its **L2 antagonist confidence is only 0.7153
  (moderate)** — the biased hypothesis rests on a very-high L3 agonist + moderate L2 antagonist.
- **Demo B (Fingolimod):** the prompt expected "no L3 = clean", but an **L3 prediction exists** (agonist 0.978).
  Since it is **concordant** with L2 (no MoA conflict), `is_clean = True` still holds by the formal definition.

**Same-drug high-confidence context** (for the "why trust only this one" objection):
- Cariprazine ≥0.95: HTR2B, HTR2A, HTR2C, HTR1A (all serotonin) + OPRM1 — a coherent serotonergic fan.
- Fingolimod ≥0.95: CCR7, GPR183, HCAR3 — immune-trafficking cluster.
- Ozanimod ≥0.95: CCR7 (L3 1.000), GPR183 (L3 antag 0.996), GRM5 (L2 antag 0.993).

**Mail-ready sentences** (each tagged `[AI-predicted hypothesis, not experimentally validated]`):
- **A:** "Cariprazine, approved for schizophrenia and bipolar disorder, carries a high-confidence agonist
  prediction at the serotonin 5-HT1A receptor (HTR1A; 5-hydroxytryptamine family) — outside its
  experimentally-evidenced dopamine-receptor family — recovering a serotonergic activity that is part of its
  documented clinical pharmacology (Kiss et al., 2010). [AI-predicted hypothesis, not experimentally validated]"
- **B:** "Fingolimod, approved for relapsing multiple sclerosis, carries a high-confidence agonist prediction at
  the hydroxycarboxylic-acid receptor HCAR3 (GPR109B; hydroxycarboxylic-acid family) — outside its
  experimentally-evidenced S1P-receptor family — a receptor on immune cells modulating neutrophil/macrophage
  function, a novel mechanistic hypothesis for its immunomodulation. [AI-predicted hypothesis, not experimentally validated]"
- **C:** "Ozanimod, approved for multiple sclerosis and ulcerative colitis, carries a biased-signalling
  prediction at the chemokine receptor CCR7 — β-arrestin-pathway agonism (very high confidence) with
  G-protein-pathway antagonism — within the same lymphocyte-trafficking axis as its evidenced S1P receptors, a
  functional-selectivity hypothesis. [AI-predicted hypothesis, not experimentally validated]"

## ⑤ PART 5 — final file
`repurposing_candidates_final.tsv` — 544 rows × **18 cols**, sorted NOVEL_INDICATION-first then confidence desc.
Added: `compound_modality` (12 peptide), `repurposing_class`, `internal_flag` (86 no_evidenced_family),
`prediction_reliability_note` (12), `target_note` (31 AGTR1), `is_demo` (**6** rows = 3 demos × 2 layers),
`demo_label` (TypeA/B/C). `repurposing_candidates_v2.tsv` retained.

## ⑥ Web-deployment recommendations
- **Web-facing columns (show to users):** inchikey, drug_name, known_indication, predicted_uniprot,
  predicted_gpcr_name, predicted_gpcr_family, evidenced_families, predicted_moa, confidence_raw, layer,
  repurposing_class. **Backend-only (hide):** compound_modality, internal_flag, prediction_reliability_note,
  target_note, is_demo, demo_label.
- **`repurposing_class` → human-readable badges** (not the TYPE_ codes):
  POLYPHARMACOLOGY → "Polypharmacology (same therapeutic area)"; NOVEL_INDICATION → "Novel-indication
  hypothesis"; BIASED_SIGNALLING → "Biased-signalling hypothesis"; PEPTIDE_UNCERTAIN → caveat icon only.
- **`known_indication` is multi-category/noisy** — for display prefer a single primary indication; the raw
  category list is better suited to filtering than to the headline cell.
- **Cariprazine is the best web demo:** it shows a *fan* of high-confidence, biologically-coherent serotonergic
  predictions (5-HT1A/2A/2B/2C) plus OPRM1 — i.e. multiple visible predictions that together read as credible
  polypharmacology, and its 5-HT1A activity is literature-validated, making it persuasive without overclaiming.

## ⑦ STALE scan
New scalars {544 · 179 · 86 · 274 · 270 · 46 · 12 · 95 · 116 · 21 · 33 · 31} ∩ §1.5 STALE = NONE. No banned
terms. (Substring "0.50"/"0.85" appears only inside `confidence_raw` decimals — not the stale operating point.)

## ⑧ Next STEP (5e-4) prep
- Editor mail: 179 repurposing candidates + the three demo sentences above (note Demo C is **biased-signalling**,
  not novel-indication; phrase accordingly). Apply the 5e-2 mail revisions (171,821→59,888; 842→496; "thousand"→~1,180).
- Web tab: apply the column split + human-readable badges above.
- Carried flags (unchanged): §E `3,719` funnel and §F `171/13`-vs-`170/14` ambiguities; paper-methods artifact-removal summary.

Files written: `repurposing_candidates_final.tsv`. Stopping per the one-STEP rule.
