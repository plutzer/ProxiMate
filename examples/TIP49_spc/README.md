# TIP49 spectral-count example

The SAINTexpress authors' demo dataset, shipped with SAINTexpress 3.6.3 as
`example input files/TIP49/`: BioID-style spectral counts for baits of the TIP49/INO80
complexes, with nine control runs. Prey and bait identifiers are gene symbols rather
than UniProt accessions, so BioGRID annotation matches nothing for it.

Run it as SAINT-format input with quant type "Spectral Counts":

```bash
./run_pipeline.sh --format saint examples/TIP49_spc/bait.txt \
  examples/TIP49_spc/prey.txt examples/TIP49_spc/interaction.txt \
  "Spectral Counts" output/ 2 0
```
