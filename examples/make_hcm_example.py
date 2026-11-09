"""Cut the GUI's example dataset (GUI/example_data/) from the Human Cell Map.

Reads the full HCM MaxQuant proteinGroups.txt (Go et al. 2021, PubMed 34079125) and its
experimental design, and keeps eight baits from distinct compartments plus the
BirA-FLAG and GFP control runs. Intensity, LFQ intensity and MS/MS count columns are
kept for those runs, so every quant type can be parsed.

    python examples/make_hcm_example.py <proteinGroups.txt> <ED.csv> GUI/example_data
"""
import os
import sys

import pandas as pd

BAITS = ["COIL (BirA-FLAG)", "LMNA (BirA-FLAG)", "NUP50 (BirA-FLAG)",
         "SEC61B_Nterm (BirA-FLAG)", "LAMP1 (C-BirA-FLAG)", "PEX14 (C-BirA-FLAG)",
         "COX4I1 (C-BirA-FLAG)", "RAB5A (BirA-FLAG)"]
CONTROLS = ["GFP (BirA-FLAG)", "BirAFLAG (BirA-FLAG)"]
PROTEIN_COLUMNS = ["Protein IDs", "Majority protein IDs", "Protein names", "Gene names",
                   "Sequence length", "Only identified by site", "Reverse",
                   "Potential contaminant", "id"]
QUANT_PREFIXES = ["Intensity ", "LFQ intensity ", "MS/MS count "]


def main(pg_path, ed_path, out_dir):
    ed = pd.read_csv(ed_path)
    ed = ed[ed["Bait"].isin(BAITS + CONTROLS)]
    missing = set(BAITS + CONTROLS) - set(ed["Bait"])
    if missing:
        raise ValueError(f"Baits not in {ed_path}: {sorted(missing)}")

    quant_columns = [prefix + run for run in ed["Experiment Name"] for prefix in QUANT_PREFIXES]
    pg = pd.read_csv(pg_path, sep="\t", usecols=PROTEIN_COLUMNS + quant_columns,
                     low_memory=False)[PROTEIN_COLUMNS + quant_columns]
    pg = pg[(pg[quant_columns] > 0).any(axis=1)]

    os.makedirs(out_dir, exist_ok=True)
    ed.to_csv(os.path.join(out_dir, "ED.csv"), index=False, lineterminator="\n")
    pg.to_csv(os.path.join(out_dir, "proteinGroups.txt"), sep="\t", index=False,
              lineterminator="\n")
    print(f"{len(ed)} runs, {len(pg)} protein groups written to {out_dir}")


if __name__ == "__main__":
    main(*sys.argv[1:])
