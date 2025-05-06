import pandas as pd
import numpy as np
import argparse
import os
from experimental_design import ExperimentalDesign
from protein_groups import ProteinGroups
import shutil

def parse_ed_pg(proteinGroups, experimentalDesign, quantType, outputPath):
    """
    This is the parse call used by the GUI.

    Parse the proteinGroups.txt file and the experimental design file.
    :param proteinGroups: path to the proteinGroups.txt file
    :param experimentalDesign: path to the experimental design file
    :param quantType: quantification type (Intensity, LFQ, Spectral Counts)
    :param outputPath: path for the output directory
    """

    print("Parsing Experimental Design file: " + experimentalDesign)

    # Print the quantification type
    print("Quantification type: " + quantType)

    # Check to see if the output directory exists. If not, create it.
    if not os.path.exists(outputPath):
        print("Creating output directory: " + outputPath)
        os.makedirs(outputPath)

    # Copy the experimental design file to the output directory
    shutil.copy(experimentalDesign, f"{outputPath}/ED.csv")
    shutil.copy(proteinGroups, f"{outputPath}/proteinGroups.txt")

    # Parse experimental design
    experimental_design = ExperimentalDesign(experimentalDesign)

    num_expts = experimental_design.num_experiments
    num_ctrls = experimental_design.num_controls

    # num_baits, num_ctrls = experimental_design.get_num_baits(), experimental_design.get_num_ctrls()

    # process MaxQuant proteinGroups
    protein_groups = ProteinGroups(experimental_design, proteinGroups,
                                   quantType, quantType)

    protein_groups.to_SAINT(outputPath)
    protein_groups.to_CompPASS(outputPath)

    return num_expts, num_ctrls

def main():

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

if __name__ == "__main__":
    # main()
    num_expts, num_ctrls = parse_ed_pg("C:/Users/plutzer/Work/ProxiMate_Testing/proteinGroups.txt",
                "C:/Users/plutzer/Work/ProxiMate_Testing/ED_CUL3_NC_split_basal_1hB_3hB_IB.csv",
                "Intensity",
                "C:/Users/plutzer/Work/ProxiMate_Testing")
    
    print("Num experiments: ", num_expts)
    print("Num controls: ", num_ctrls)