# SynerGPCR

SynerGPCR is a resource and web application for exploring GPCR-relevant drug
combination synergy predictions. The live app is available at:

https://synergpcr.kaist.ac.kr

## Layout

- `pipeline/` — the data pipeline that produces the released tables (chain
  assembly, Bliss synergy index, Wilson confidence intervals, enrichment
  analysis, and LLM-based annotation). See `pipeline/README.md` for
  reproducibility tiers.
- `webapp/` — the Flask/web application that serves the interactive site.
- `data/examples/` — small example tables shipped with the repo so the
  webapp and pipeline can be run without downloading the full datasets.
  See `data/examples/README.md`, including a note on UMLS CUI usage.

## Data

The webapp reads its data directory from the `SYNERGPCR_DATA` environment
variable. Point it at `data/examples/` to run against the bundled example
data, or at a local copy of the full released tables.

The full datasets are archived on Zenodo (not bundled in this repository).
GPCRact model weights are maintained in a separate repository and are not
redistributed here.

## Licence

- Code: MIT — see `LICENSE`.
- Data (`data/examples/` and the full tables archived on Zenodo):
  CC BY-NC-SA 4.0 — see `LICENSE-DATA`.

This dual licensing is intended to simultaneously satisfy the ChEMBL
share-alike requirement and the DrugBank non-commercial requirement for the
underlying source data.

## Citation

See `CITATION.cff`.
