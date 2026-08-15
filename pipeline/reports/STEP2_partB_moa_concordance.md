# STEP 2 — Part B: MoA-direction concordance (coverage + trend, Task-2 preview) — **DONE**

target_analysis sub-task. Sources READ-ONLY. This is a **preview/scoping** for Task 2, kept strictly
**separate from Part A** (Part A's `co_active` = L1∩L2 *co-activity*; Part B = MoA *directional*
concordance — different concepts, never mixed). Neutral labels only (no DSAV/forward/backward).
Build: `06_target_enrichment/moa_concordance_partB.py`, `step2_partB_clin.py (build script not released)`. Env: pandas 2.2.2 / statsmodels 0.14.5, seed=0.
All trends are hypothesis-generating; n-limits stated explicitly.

---

## B1. Existing analyses on the server (definitions + numbers, NOT reinvented)
Three prior, **mutually inconsistent-in-direction** analyses exist (this inconsistency is exactly the
methodological concern raised — they are different concordance constructs on different populations):

| file (mtime) | concordance construct | concordant | discordant | direction |
|---|---|---|---|---|
| `Output/NAR/Tables/SI_MoA_aware.csv` (step A) | L2 label → `moa_concordance_status` (L2–L4 dir), SI synergy pool n=561 | 385 / 123 = **31.9%** | 176 / 25 = **14.2%** | concordant **higher** |
| `Output/DB/GPCRactDB/Analysis/L2L4_concordance_results.csv` (2026-03-31) | L2–L4 MoA dir over co-active pool | 11,713 / 240 = **2.05%** | 685 / 34 = **4.96%** (biased agonism) | concordant **lower** |
| `…/PreTrain/Full_chain_L3_Concordance_stats_v3.csv` (2026-04-06) | L2–L3 dir (full-activation/balanced/biased) | agg 195 / 32 = 16.4% | biased agg 12 / 2 = 16.7% (n=12) | **no difference, underpowered** |

→ "MoA concordance → approval" is **not one number**: the sign flips with the construct/population, and
the L2–L3 (biased-signaling) arm is underpowered (n=8–12). All three predate the 0.40/0.97 cleanup → treat
as lineage only; Task 2 must re-derive with neutral labels.

---

## B2. Direction collapse (locked for this preview) + observed class distribution
- **agonist-direction** = {Agonist, Partial Agonist, PAM}; **antagonist-direction** = {Antagonist,
  **Inverse Agonist**, NAM}. *Inverse Agonist is assigned to antagonist-direction* (negative efficacy /
  functional block); flagged because it is the one arguable case. Inactive/Unknown excluded.

| layer | agonist-dir records | antagonist-dir records |
|---|---:|---:|
| L2 G-protein | 57,122 | 52,903 |
| L3 β-arrestin | 4,154 | 2,874 |
| L4 Reporter | 25,567 | 73,411 (antagonist-heavy) |

---

## B3 / B4. Pairwise both-direction coverage + concordant-vs-discordant approval (Wilson CI)
Per `(compound, GPCR)` with an unambiguous directional call in **both** layers; approval = `max_phase==4`,
unique-compound dedup. Approval rate compared concordant (same dir) vs discordant (opposite dir):

| layer pair | both-dir pairs | concordant n / appr / rate [95% CI] | discordant n / appr / rate [95% CI] | adequacy |
|---|---:|---|---|---|
| **L2 ∩ L4** | 8,880 | 8,326 / 162 / **1.95%** [1.67, 2.27] | 116 / 7 / **6.03%** [2.95, 11.93] | both ≥20; **CIs non-overlapping** → discordant higher |
| **L2 ∩ L3** | 1,171 | 1,091 / 31 / **2.84%** [2.01, 4.00] | 32 / 2 / **6.25%** [1.73, 20.15] | discordant n=32, **CI very wide** → weak |
| **L3 ∩ L4** | 699 | 653 / 29 / 4.44% [3.11, 6.31] | 15 / 2 / 13.33% [3.74, 37.88] | discordant n=15 **<20 → INSUFFICIENT** |

**assay (L2 direction) vs clinical (`Catalog_v3.human_dir`):**
- compounds with both = **1,438**; **match rate 62.3%** (896 concordant / 542 discordant).
- concordant approval **16.52%** [14.23, 19.09] vs discordant **13.10%** [10.52, 16.20] — **CIs overlap**
  (weak positive: assay↔clinical agreement modestly tracks approval, not significant).

### Reading the trend (the key Task-2 message)
- **Within-assay cross-layer** direction *discordance* (e.g. L2 agonist but L4 antagonist) is associated
  with **higher** approval (L2∩L4: 6.03% vs 1.95%, non-overlapping). This is **biased signaling /
  functional selectivity — real biology, NOT an annotation error** (caption). Directionally divergent
  compounds are over-represented among approved drugs.
- **Assay-vs-clinical** MoA agreement instead shows a **weak positive** (concordant slightly higher,
  overlapping CIs).
- → These two "concordances" point in **opposite directions** and must be reported as **distinct
  metrics**; conflating them is precisely the prior pitfall. Adequately-powered comparisons: L2∩L4 and
  assay-vs-clinical; **L3-involving** arms (β-arrestin / biased) remain **n-limited** (discordant n=15–32)
  → report as limitation, do not conclude.

---

## Limitations / Task-2 handoff
- L3 (β-arrestin) directional discordance is the biologically interesting "biased signaling" axis but is
  **underpowered** everywhere it appears (n=8–32) → Task 2 should report coverage and flag, not force a test.
- `human_dir` clinical-direction column covers only 7,830 compounds (2,256 agonist / 5,574 antagonist);
  the assay-vs-clinical join here used L2 only — Task 2 should confirm the canonical clinical-MoA source
  and extend to L3/L4.
- No banned terms used; Part A `co_active` concept not referenced. Full Task-2 development deferred.
