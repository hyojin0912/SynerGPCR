# SynerGPCR Example Data

The seven CSV files below mirror exactly the public **/download** tab on the
live site (https://synergpcr.kaist.ac.kr) — same filenames, same columns.
They let the web application and pipeline code be reviewed against real
column schemas without needing the full released tables.

- compound_lookup.csv
- target_lookup.csv
- compound_chain_summary.csv
- gpcract_predictions.csv
- candidate_shortlist.csv
- synergy_stats.csv
- single_layer_table.csv

One additional file, `gpcract_3point_si.csv`, is included alongside them —
it is read directly by the webapp (`webapp/main.py`) to render an in-app
figure and is not itself offered as a `/download` file.

Full datasets are archived on Zenodo (see the root README for the DOI link).
