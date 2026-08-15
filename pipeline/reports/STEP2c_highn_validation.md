# STEP 2c — high-n / 0-approved target validation + re-label — **DONE**

target_analysis sub-task. Sources READ-ONLY; only `target_clinical_enrichment.tsv` extended (no
enrichment value recomputed). Build: `06_target_enrichment/highn_validation.py`. Env: pandas 2.2.2, seed=0. FINAL_STATS Task-1
anchors independently re-asserted **unchanged** (write-gated).

> **Headline correction:** the STEP-2b label *"heavily assayed, not clinically entered"* is **WRONG for
> every high-n target**. All 41 enough_n & 0-approved targets — **including all 24 `high_n_no_clinical`** —
> actually **have an approved drug**; their 0 approved-co_active is a **co_active coverage artifact**, not a
> clinical gap. `no_approved_drug` count = **0**.

---

## 1. Independent cross-check (co_active-agnostic) of the 41 enough_n & 0-approved targets
Three independent evidence streams per target UniProt, none using the co_active set:
- **approved_indication_links** = # Drug_Indication 'Approved' drugs sharing an indication CUI with the
  target (Drug_Indication ⋈ Target_Indication).
- **approved_anyrecord_at_target** = # Catalog `max_phase==4` compounds with ANY Bioassay record at the
  target (layer-agnostic).
- **approved_active_at_target** = same, restricted to ACTIVE MoA records (any layer).
`approved_drug_exists_any` = OR of the three; `clinical_status_label` ∈
{`approved_outside_coactive`, `no_approved_drug`} assigned only to enough_n & 0-approved targets.

### Result
| group | n | approved_outside_coactive | no_approved_drug |
|---|---:|---:|---:|
| all enough_n & 0-approved | 41 | **41** | **0** |
| …of which high_n_no_clinical | 24 | **24** | **0** |

**Every** heavily-assayed 0-co_active-approved target has independent approved-drug evidence → the
"depleted" appearance is entirely a **denominator/coverage artifact** of the L1∩L2 co_active definition.

---

## 2. Columns added (backend) to `target_clinical_enrichment.tsv` (29 → 35 cols)
- `approved_indication_links_backend`, `approved_anyrecord_at_target_backend`,
  `approved_active_at_target_backend` (the three evidence counts)
- `approved_drug_exists_any` (bool), `approved_drug_not_in_coactive` (bool; True iff exists AND
  n_approved_exp==0), `clinical_status_label` (string; blank unless enough_n & 0-approved)

---

## 3. Explicit P2RY12 / TACR1 verification — both `approved_outside_coactive` ✓
| GPCR | UniProt | n_co_active | appr in co_active | indic links | anyrecord | active | label |
|---|---|---:|---:|---:|---:|---:|---|
| **P2RY12** | Q9H244 | 132 | 0 | 289 | 83 | 83 | **approved_outside_coactive** |
| **TACR1** | P25103 | 52 | 0 | 1464 | 763 | 763 | **approved_outside_coactive** |
| MCHR1 | Q99705 | 760 | 0 | 191 | 17 | 17 | approved_outside_coactive |
| ADORA3 | P0DMS8 | 601 | 0 | 1401 | 1555 | 1555 | approved_outside_coactive |
| GCGR | P47871 | 354 | 0 | 515 | 314 | 313 | approved_outside_coactive |

**Likely causes the approved drug escapes the co_active (L1∩L2) capture (reported, not asserted):**
- **P2RY12** — approved antiplatelets (clopidogrel, prasugrel) are **prodrugs** whose active metabolites
  bind **covalently/irreversibly**; reversible binding (L1) + G-protein (L2) assays on the parent
  InChIKey frequently don't register as co-active. Classic prodrug + covalent case.
- **TACR1** — NK1 antagonists (aprepitant etc.) may be recorded at one layer (e.g. L2 antagonist or a
  functional readout) but **lack a matched L1 binding record** for the exact approved InChIKey → no L1∩L2
  intersection. Assay-record missingness / antagonist-coverage gap.
- General: prodrug (inactive parent), covalent/irreversible binders, single-layer assay coverage, or
  salt/stereo InChIKey mismatch between the approved-drug key and the assay key.

---

## 4. Web/DB labeling rule (applied)
The exposure label **"heavily assayed, not clinically entered"** is restricted to
`clinical_status_label == "no_approved_drug"` targets — of which there are **ZERO**. → **None of the 24
high-n targets** should carry that label in web/db; they are surfaced instead as
**`approved_outside_coactive` (co_active coverage artifact)**. The `high_n_no_clinical_exp` column is kept
(it factually means high-n & 0 approved *within co_active*) but **must be read with
`clinical_status_label`**, which overrides any "no clinical entry" reading.

This sharpens the Task-1 caveat: a low/zero **conditional** approval rate at a target is **not** evidence
of "no approved drug" — it can be a co_active coverage miss (prodrug / covalent / assay gap). Reinforces
the hypothesis-generating framing.

---

## Anchor re-assert (unchanged) — 8/8 PASS
universe 211 · sum n_co_active_exp 20,145 · sum n_approved_exp 282 · enough_n 126 · point-enriched 50 ·
CI-robust target 28 · CI-robust aiincl 31 · n≥5 & 0-approved 41 — **all unchanged**. No enrichment number
moved; this STEP only adds labels/evidence columns.

## STALE / banned scan
New scalars {41 · 24 · 0} ∩ §1.5 STALE = NONE. No banned terms. Output:
`target_clinical_enrichment.tsv` (211 × 35) + this report. Stopping per the one-STEP rule.
