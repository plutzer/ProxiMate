# Example datasets

Public datasets in SAINT input format, small enough to run the pipeline end to end.

| Directory | Source | Quant type | Notes |
| --- | --- | --- | --- |
| `TIP49_spc/` | SAINTexpress 3.6.3 demo data (Teo et al. 2014) | Spectral Counts | Gene-symbol identifiers, so BioGRID annotation matches nothing. See its README for the command. |
| `HCM_LFQ/` | Human Cell Map (Go et al. 2021, PubMed 34079125) | LFQ | `bait.txt`, `prey.txt` and the experimental design `ED.csv`. The interaction file is too large to track; regenerate it from the MaxQuant `proteinGroups.txt` in the PRIDE deposit with `run_pipeline.sh --format maxquant`. |

Each set is in the layout `run_pipeline.sh --format saint` expects: `bait.txt` (run,
bait, C/T), `prey.txt` (prey identifier plus one annotation column), `interaction.txt`
(run, bait, prey, quantity), all tab-delimited without headers.
