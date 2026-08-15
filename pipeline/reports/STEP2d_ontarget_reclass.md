# STEP 2d — on-target approved-drug reclassification — **DONE**

target_analysis sub-task. Sources READ-ONLY; only `target_clinical_enrichment.tsv` extended (no
enrichment value changed). Build: `06_target_enrichment/ontarget_reclass.py`. Env: pandas 2.2.2, seed=0. FINAL_STATS Task-1
anchors independently re-asserted **unchanged** (write-gated, 8/8 PASS).

> **Why STEP 2c was wrong:** 2c's evidence (`approved_indication_links` = disease-CUI sharing;
> `approved_active_at_target` = any approved drug with an *activity record* at the target) counts
> **off-target activity and disease co-occurrence**, so it returned "approved drug exists" for **all 41**
> (0 truly missing) — structurally unable to detect a genuinely un-drugged target. This STEP replaces it
> with a true **on-target approved-drug** test and demotes the 2c columns to explicit `_noise` labels.

---

## 1. On-target approved-drug reference (co_active- and indication-independent)
A target is "on-target drugged" if an **approved** drug acts on it via its **mechanism of action**:
- **DrugCentral** — target development level **`TDL == 'Tclin'`** (= target has an approved drug with
  known MoA); take its MoA drugs (`ACTION_TYPE` present). 519 Tclin targets total, **96 in the 211 universe**.
- **ChEMBL `drug_mechanism`** (curated MoA) — drugs whose `standard_inchi_key` is in the approved set
  (Catalog max_phase==4 ∪ IUPHAR approved ∪ DrugCentral approved = 2,532 InChIKeys; spans peptides).
- DrugCentral Tclin MoA rows = 1,085; ChEMBL approved MoA rows = 1,626.
- **Modality** per drug: IUPHAR `Type` (Peptide/Antibody → peptide_biologic; else small_molecule),
  SMILES-amide heuristic as fallback.

This is exactly what distinguishes a real on-target drug (clopidogrel→P2RY12, glucagon→GCGR) from
**off-target/secondary** activity (caffeine/theophylline→ADORA3, rosiglitazone→FFAR1, which are
DrugCentral **Tchem**, not Tclin).

## 2. Columns added (40 cols total) + 2c demotion
- new: `on_target_approved_n`, `on_target_smallmol_n`, `on_target_peptide_biologic_n`,
  `on_target_sources`, `clinical_status_label_v2`.
- demoted (renamed, kept): `approved_indication_links_backend` → `disease_indication_links_noise_backend`;
  `approved_anyrecord_at_target_backend` → `any_record_at_target_offtarget_noise_backend`;
  `approved_active_at_target_backend` → `active_record_at_target_offtarget_noise_backend`;
  `clinical_status_label` → `clinical_status_label_v1_superseded`.

`clinical_status_label_v2 ∈ {approved_on_target_smallmol_outside_coactive,
approved_on_target_peptide_or_biologic, no_on_target_approved_drug}` (set only for enough_n & 0-approved).

## 3. Reclassification of the 41 enough_n & 0-approved targets — **2c's "0 missing" corrected to 29**
| label | n | members (n_co_active) |
|---|---:|---|
| `approved_on_target_smallmol_outside_coactive` (assay-coverage artifact) | **8** | CRHR1(258), ADORA2B(221), P2RY12(132), TACR1(52), CCR4(32), C5AR1(21), CASR(6), EDNRB(5) |
| `approved_on_target_peptide_or_biologic` (modality artifact) | **4** | GCGR(354), SSTR3(182), CALCR(116), SSTR4(15) |
| `no_on_target_approved_drug` (genuine non-success) | **29** | MCHR1(760), ADORA3(601), GRM5(230), CCR2(154), BDKRB1(151), NPY1R(129), FFAR1(125), CXCR3(123), NPY2R(119), NPY4R(95), UTS2R(94), CCR3(85), CCR1(85), P2RY1(73), KISS1R(72), GPR119(71), APLNR(62), NPY5R(47), NMUR2(45), BRS3(44), DRD5(34), GPR84(33), CXCR2(30), NTSR1(30), LPAR1/3, GRM4(19), FFAR2(17), LTB4R2(13) … |

→ **12 of 41 are artifacts** (8 small-molecule assay-coverage misses + 4 peptide/biologic modality
misses); **29 are genuinely without an on-target approved drug** — biologically sensible (MCHR1, GRM5/mGlu5,
CCR2, FFAR1/GPR40, KISS1R, GPR119: well-known targets with failed/no approved programs).

## 4. Explicit verification — all 7 expected, correct ✓
| GPCR | on_target n (sm/pb) | source | label | expected |
|---|---|---|---|---|
| P2RY12 | 10 (10/0) | chembl;drugcentral_Tclin | smallmol_outside_coactive | smallmol ✓ |
| TACR1 | 10 (9/1) | chembl;drugcentral_Tclin | smallmol_outside_coactive | smallmol ✓ |
| GCGR | 2 (0/2) | drugcentral_Tclin | peptide_or_biologic | peptide ✓ |
| CALCR | 2 (0/2) | chembl;drugcentral_Tclin | peptide_or_biologic | peptide ✓ |
| MCHR1 | 0 | — | no_on_target_approved_drug | no_on_target ✓ |
| ADORA3 | 0 | — | no_on_target_approved_drug | no_on_target ✓ |
| FFAR1 | 0 | — | no_on_target_approved_drug | no_on_target ✓ |

## 5. Web/DB labeling rule (applied)
The exposure label **"heavily assayed, no approved drug"** applies **only** to
`clinical_status_label_v2 == "no_on_target_approved_drug"` → **29 targets** (of which **17** are the
high-n `high_n_no_clinical_exp` subset). The other 12 high-/mid-n 0-co_active-approved targets are NOT
"no drug": **8 are small-molecule assay-coverage artifacts** (approved drug exists, missed by L1∩L2 —
prodrug/covalent/coverage, e.g. P2RY12, TACR1) and **4 are peptide/biologic modality artifacts**
(GPCRact's small-molecule assays inherently can't represent the peptide drug, e.g. GCGR, CALCR).

## 6. Anchor re-assert (unchanged) — 8/8 PASS
universe 211 · sum n_co_active_exp 20,145 · sum n_approved_exp 282 · enough_n 126 · point-enriched 50 ·
CI-robust target 28 · CI-robust aiincl 31 · 0-approved enough_n 41 — all unchanged.

## Limitations
- `Tclin` gate is from DrugCentral (single source); ChEMBL `drug_mechanism` corroborates at drug level.
  A target drugged only in a source lacking a Tclin flag could in principle be under-called (low risk for
  GPCRs; cross-checked by ChEMBL).
- Modality uses IUPHAR `Type` where available, else a SMILES-amide heuristic; rare borderline peptoids may
  be mis-binned.
- ChEMBL `drug_mechanism` is gated to approved by InChIKey match; peptide InChIKey variants are covered via
  the DrugCentral Tclin route.

Output: `target_clinical_enrichment.tsv` (211 × 40) + this report. STALE/banned scan: scalars {8·4·29·96}
∩ §1.5 STALE = NONE; no banned terms. Stopping per the one-STEP rule.
