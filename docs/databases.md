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

## Excluding Human Cell Map evidence

The Human Cell Map (Go et al. 2021, PubMed 34079125) is a large BioID screen
deposited in BioGRID. Scoring a proximity-labeling experiment against it counts
hits as "known" on the strength of the same kind of assay.

To avoid that, tick **Exclude Human Cell Map evidence** on the Score card, or pass
`--exclude-hcm` to `run_pipeline.sh`. Human datasets only.

- Annotation then uses `/Datasets/human/biogrid_summary_no_hcm.csv`: the human
  BioGRID summary with every evidence line from that paper removed.
- Interactions with other evidence are kept; only their HCM entries are dropped.
- `run.json` records the choice.
- The Quality Control tab's known-interaction and network-degree metrics use the
  same file the run was annotated against.
