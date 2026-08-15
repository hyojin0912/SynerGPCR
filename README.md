# SynerGPCR

SynerGPCR integrates GPCR in vitro functional assay data from ten public sources into per-compound, four-layer assay-chain profiles (receptor binding → G-protein coupling → β-arrestin recruitment → reporter gene activity), linked to clinical approval outcomes and completed by GPCRact predictions where experimental data is missing. The live app is available at: https://synergpcr.kaist.ac.kr

<img width="2103" height="1902" alt="Figure1" src="https://github.com/user-attachments/assets/52865cdd-b230-48c3-a608-93abf2b15273" />

*Ten public sources are harmonised to common compound and GPCR identifiers, linking 1,753,921 assay records to clinical status for 2,420 approved drugs. 
Every compound–GPCR pair is scored active, inactive or untested at four assay layers ordered along the signalling cascade, and active layers are combined across targets into a compound-level chain pattern. GPCRact completes chains missing L2 or L3 data; four entry points serve the resulting profiles.*

## Layout

- `pipeline/` — the data pipeline that produces the released tables (chain assembly, Bliss synergy index, Wilson confidence intervals, enrichment analysis, and LLM-based annotation). See `pipeline/README.md` for reproducibility tiers.
- `webapp/` — the FastAPI web application that serves the interactive site. 
- `data/examples/` — small example tables shipped with the repo so the webapp and pipeline can be run without downloading the full datasets.
  See `data/examples/README.md`, including a note on UMLS CUI usage.

## Data

The webapp reads its data directory from the `SYNERGPCR_DATA` environment
variable. Point it at `data/examples/` to run against the bundled example
data, or at a local copy of the full released tables.

The full datasets are archived on [Zenodo](https://zenodo.org/records/21943761).
GPCRact model weights are maintained in a separate repository
([hyojin0912/HJ-GPCRact](https://github.com/hyojin0912/HJ-GPCRact)) and are
not redistributed here.

## Licence

- Code: MIT — see `LICENSE`.
- Data (`data/examples/` and the full tables archived on Zenodo):
  CC BY-NC-SA 4.0 — see `LICENSE-DATA`.

This dual licensing is intended to simultaneously satisfy the ChEMBL
share-alike requirement and the DrugBank non-commercial requirement for the
underlying source data.

## Citation

A manuscript describing SynerGPCR is currently under review. Citation details
will be added here and in `CITATION.cff` upon publication.

## Contact
Hyojin Son (sonhyojin0912@gmail.com)
