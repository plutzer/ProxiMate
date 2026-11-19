# Databases

ProxiMate annotates interactions against BioGRID, UniProt, the Human Protein
Atlas and CORUM.

## Which snapshot a release uses

Each release ships the BioGRID, UniProt and Human Protein Atlas files downloaded
on the day it was built. CORUM is tracked in the repository. Find the download
dates in any of:

- `/Datasets/build_info.txt` inside the image
- the release's notes on GitHub
- every `run.json`

## Sources

| File | Source |
| --- | --- |
| `BIOGRID-ALL.tab3.txt`, `BIOGRID-MV-Physical.tab3.txt` | [BioGRID](https://downloads.thebiogrid.org/BioGRID) |
| `uniprot_anns.tsv` (reviewed proteome annotations) | [UniProt](https://www.uniprot.org/uniprotkb?query=%28proteome%3AUP000005640%29&facets=reviewed%3Atrue) |
| `subcellular_location.tsv` | [Human Protein Atlas](https://www.proteinatlas.org/humanproteome/subcellular/data#locations) |
| `corum_humanComplexes.txt` | [CORUM](https://mips.helmholtz-muenchen.de/corum/download) |

To build an image with a snapshot you choose, see
[Development](development.md#pin-a-database-snapshot).

## How BioGRID interactions are matched

BioGRID records an interaction once, bait first. `biogrid_summary.csv` lists every
pair in both orientations with the evidence pooled, so:

- an interaction is `In.BioGRID` whichever partner was the bait in the published
  experiment;
- `Multivalidated` is the same for both orientations;
- only physical interactions with both partners in the dataset's organism count.

## Excluding proximity evidence

BioGRID holds proximity-labeling screens (experimental system `Proximity Label-MS`),
such as the Human Cell Map. Scoring a proximity-labeling experiment against them
counts hits as "known" on the strength of the same kind of assay.

To avoid that, tick **Exclude proximity evidence** on the Score card, or pass
`--exclude-proximity` to `run_pipeline.sh`.

- Annotation then uses `/Datasets/<organism>/biogrid_summary_no_proximity.csv`: the
  BioGRID summary with every `Proximity Label-MS` evidence line removed.
- Interactions with other evidence are kept; only their proximity entries are dropped.
- `Multivalidated` is recomputed from what is left: two publications or two
  experimental systems.
- `run.json` records the choice.
- The Quality Control tab's known-interaction and network-degree metrics use the
  same file the run was annotated against.
