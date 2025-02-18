import pandas as pd
import numpy as np
import argparse
import os
from experimental_design import ExperimentalDesign
from protein_groups import ProteinGroups
import shutil

description = "This is the entry point to the program. It will execute the requested tasks."

# initialize the parser
parser = argparse.ArgumentParser(description=description)

# Add arguments for the annotator
parser.add_argument("--proteinGroups",
                    help="path to MaxQuant ProteinGroups.txt",
                    default=None)

parser.add_argument("--experimentalDesign",
                    help="path to experimental design file",
                    default=None)

parser.add_argument("--outputPath",
                    help="path for the output directory. If it already exists it will be overwritten.",
                    default="/srv/shiny-server/myapp/score_inputs")

# Argument for quantification type
parser.add_argument("--quantType",
                    help="SAINT: quantification type (Tntensity, Spectral Counts, LFQ)",
                    default="Intensity")

args = parser.parse_args()

# MAIN

if args.proteinGroups is not None:
    if args.experimentalDesign is not None:
        # Check to see if the output directory exists. If not, create it.
        if not os.path.exists(args.outputPath):
            os.makedirs(args.outputPath)

        # Copy the experimental design file to the output directory
        shutil.copy(args.experimentalDesign, os.path.join(args.outputPath, "ED.csv"))

        # Parse the proteinGroups.txt file
        print("Parsing Experimental Design file: " + args.experimentalDesign)

        # parse experimental design
        experimental_design = ExperimentalDesign(args.experimentalDesign)

        print("Parsing Protein Groups file: " + args.proteinGroups)
        print("Quantification type: " + args.quantType)

        # process MaxQuant proteinGroups
        protein_groups = ProteinGroups(experimental_design, args.proteinGroups,
                                    args.quantType, args.quantType)

        protein_groups.to_SAINT(args.outputPath)
        protein_groups.to_CompPASS(args.outputPath)
    else:
        print("Error: No experimental design file provided.")
else:
    print("Error: No proteinGroups file provided.")