#!/usr/bin/env python
import sys
import os
import argparse
import subprocess
import pandas as pd
import aft_impute_saint
import refactored_aft
import compPASS_pval
import time

SAINT_EXPRESS_INT_DIR = "/bin/SAINTexpress-int"
SAINT_EXPRESS_INT_DEFAULT_DIR = "/bin/SAINTexpress-int_default"
SAINT_EXPRESS_SPC_DIR = "/bin/SAINTexpress-spc"

def main():
    start_time = time.time()

    ########################################################################################################################
    # Command line argument parsing
    ########################################################################################################################

    description = "This is the entry point to the program. It will execute the requested tasks."

    # initialize the parser
    parser = argparse.ArgumentParser(description=description)

    # Experiment Design
    parser.add_argument("--experimentalDesign",
                        help="path to experimental design file",
                        default=None)

    # Input Files
    parser.add_argument("--scoreInputs",
                        help="path for score file inputs. If it already exists it will be overwritten.",
                        default="/srv/shiny-server/myapp/score_inputs")

    # Output Path
    parser.add_argument("--outputPath",
                        help="path for the output directory. If it already exists it will be overwritten.",
                        default="/srv/shiny-server/myapp/score_outputs")

    # SAINT arguments

    parser.add_argument("--compress-n-rep",
                        help="SAINTexpress: the number of test bait replicates used for scoring, with priotiy given to the\n" +
                            "baits with higher probability scores. If this number is greater than or equal to the number\n" +
                            "of available replicates, then the scores will use the data from all replicates. Otherwise, \n" +
                            "the highest scoring replicate scores wil lbe averaged to yield the final probability score.\n" +
                            "Default: 1000",
                        type=int,
                        default=1000)

    # CompPASS arguments

    # Argument for number of iterations to run CompPASS resampling
    parser.add_argument("--n-iterations",
                        help="CompPASS: number of iterations to run CompPASS resampling",
                        type=int,
                        default=5)

    # Argument for imputation
    parser.add_argument("--imputation",
                        help="SAINT: whether to impute missing values (1) or not (0)",
                        default="0")

    # Argument for quantification type
    parser.add_argument("--quantType",
                        help="SAINT: quantification type (intensity, spc, LFQ)",
                        default="Intensity")

    args = parser.parse_args()

    ########################################################################################################################
    # Main
    ########################################################################################################################

    # Parse the proteinGroups.txt file

    # Check to see if the output directory exists. If not, create it.
    if not os.path.exists(args.outputPath):
        os.makedirs(args.outputPath)

    # Run the imputation
    if args.imputation == "2":
        print("Running refactored AFT imputation...")
        refactored_aft.filter_impute(
            f"{args.scoreInputs}/prey.txt",
            f"{args.scoreInputs}/interaction.txt",
            f"{args.scoreInputs}/",
            args.experimentalDesign,
            impute=True)
    else:
        aft_impute_saint.filter_impute(
            f"{args.scoreInputs}/prey.txt",
            f"{args.scoreInputs}/interaction.txt",
            f"{args.scoreInputs}/",
            args.experimentalDesign,
            impute=bool(int(args.imputation)))

    ########################################################################################################################

    # Run SAINT
    print("Running SAINT...")
    if args.quantType == "Spectral Counts":
        if args.imputation == "1":
            # TODO: Imputation for spectral counts
            print("Imputation for spectral counts not yet implemented. Running SPC SAINT without imputation...")
            p = subprocess.run([SAINT_EXPRESS_SPC_DIR,
                                "-L",
                                str(args.compress_n_rep),
                                "filtered_interaction.txt",
                                "prey.txt",
                                "bait.txt"],
                                cwd=args.scoreInputs)
        else:
            # Spectral Counts
            p = subprocess.run([SAINT_EXPRESS_SPC_DIR,
                                "-L",
                                str(args.compress_n_rep),
                                "filtered_interaction.txt",
                                "prey.txt",
                                "bait.txt"],
                                cwd=args.scoreInputs)
    else:
        if args.imputation in ("1", "2"):
            p = subprocess.run([SAINT_EXPRESS_INT_DIR,
                                "-L",
                                str(args.compress_n_rep),
                                "filtered_interaction.txt",
                                "imputed_prey.txt",
                                "bait.txt"],
                                cwd=args.scoreInputs)
        else:
            # Intensity
            p = subprocess.run([SAINT_EXPRESS_INT_DEFAULT_DIR,
                                "-L",
                                str(args.compress_n_rep),
                                "filtered_interaction.txt",
                                "prey.txt",
                                "bait.txt"],
                                cwd=args.scoreInputs)

    ########################################################################################################################
    ########################################################################################################################
    # Run CompPASS
    print("Running CompPASS...")

    # Read input file
    comp_input = pd.read_csv(f"{args.outputPath}/to_CompPASS.csv", sep=',', index_col=0).astype({'Prey': str, 'Bait': str})

    compPass_result = compPASS_pval.score_compPass(comp_input, 0.98, args.n_iterations)

    compPass_result.to_csv(f"{args.outputPath}/compPASS.csv", index=False)

    ########################################################################################################################
    ########################################################################################################################

    # Merge CompPASS and SAINT outputs
    print("Merging CompPASS and SAINT outputs...")
    # Read in the SAINT output
    saint = pd.read_csv(f"{args.scoreInputs}/list.txt", sep="\t")

    # Merge the two dataframes on two columns:
    # "Bait" in SAINT corresponds to "Experiment.ID" in CompPASS
    # "Prey" in SAINT corresponds to "Prey" in CompPASS
    merged = pd.merge(saint, compPass_result, how="left", left_on=["Bait", "Prey"], right_on=["Experiment.ID", "Prey"])

    # Rename Columns
    merged.rename(columns={
        "Bait_y": "Bait.ID",
        "Prey": "Prey.ID",}, inplace=True)

    # Drop Columns
    merged.drop(columns=["Bait_x","boosted_by"], inplace=True)

    # Write the merged dataframe to a file
    merged.to_csv(f"{args.outputPath}/merged.csv", index=False)

if __name__ == "__main__":
    main()
