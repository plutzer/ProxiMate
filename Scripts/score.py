#!/usr/bin/env python
import sys
import os
import argparse
import subprocess
import pandas as pd
import aft_impute_saint
#import tab_to_JSON

print("Running parse_scoreSAINT.py...")

SAINT_EXPRESS_INT_DIR = "/bin/SAINTexpress-int"
SAINT_EXPRESS_INT_DEFAULT_DIR = "/bin/SAINTexpress-int_default"
SAINT_EXPRESS_SPC_DIR = "/bin/SAINTexpress-spc"
FILTER_INTERACTION_PATH = "C:/Users/plutzer/Repos/scoreAPMS/filter_interaction.R"

# Time this script
import time
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

# with open('/tmp/progress.txt', 'w') as f:
#     f.write("1|Running score script...")

# Check to see if the output directory exists. If not, create it.
if not os.path.exists(args.outputPath):
    os.makedirs(args.outputPath)

print("Quantification type: " + args.quantType)

# Run the imputation
aft_impute_saint.filter_impute(f"{args.scoreInputs}/prey.txt", f"{args.scoreInputs}/interaction.txt", f"{args.scoreInputs}/", args.experimentalDesign, impute=bool(int(args.imputation)))

########################################################################################################################

# with open('/tmp/progress.txt', 'w') as f:
#     f.write("70|Running SAINT...")

# Run SAINT
print("Running SAINT...")
print(args.quantType)
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
    if args.imputation == "1":
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
        

# with open('/tmp/progress.txt', 'w') as f:
#     f.write("79|Preparing files for scoring...")

########################################################################################################################
# Run CompPASS R script
print("Running CompPASS...")
p = subprocess.run(['Rscript',
                    '/neo_comppass.R',
                    args.scoreInputs,
                    args.outputPath,
                    str(args.n_iterations)],
                    cwd="/") # Output file will be called compPASS.csv and compPASS_resampled.csv

print("Calculating WD-FDR from Resampled CompPASS...")

# Read in the CompPASS output
compPASS = pd.read_csv(f"{args.outputPath}/compPASS.csv")

# Read in Resampled CompPASS output
comp_resampled = pd.read_csv(f"{args.outputPath}/compPASS_resampled.csv")

wds = comp_resampled.loc[:, comp_resampled.columns.str.contains('WD_') & ~comp_resampled.columns.str.contains('test')]
wd_list = pd.DataFrame(wds.apply(lambda x: x.dropna().tolist(), axis=1))
wd_list['WD'] = compPASS['WD']
wd_list['FDR'] = wd_list.apply(lambda x: sum(i > x['WD'] for i in x[0])/len(x[0]), axis=1)
compPASS['WDFDR'] = wd_list['FDR']

# Overwrite the original compPASS.csv file with the new one
compPASS.to_csv(f"{args.outputPath}/compPASS.csv", index=False)


########################################################################################################################
# # Run numpPASS (Under testing)
# import numpPASS

# print("Running numpPASS...")

# # Read in the CompPASS output
# try:
#     input = pd.read_csv(f'{args.scoreInputs}/to_CompPASS.csv', sep='\t', index_col=0).astype({'Prey': str, 'Bait': str})
# except:
#     print("Couldn't read input file with tab separator. Trying Commas.")
#     input = pd.read_csv(f'{args.scoreInputs}/to_CompPASS.csv', sep=',', index_col=0).astype({'Prey': str, 'Bait': str})

# result = numpPASS.score_compPass(input, 0.98, args.n_iterations)

# # Write the output to a file
# result.to_csv(f"{args.outputPath}/compPASS.csv", index=False)

# compPASS = result

########################################################################################################################
print('MADE IT HERE!!!')

# Merge output files into a single file

# Merge CompPASS and SAINT outputs for just intensity
print("Merging CompPASS and SAINT outputs for intensity...")
# Read in the SAINT output
saint = pd.read_csv(f"{args.scoreInputs}/list.txt", sep="\t")

# Print the number of rows in each dataframe
print("CompPASS rows: " + str(len(compPASS.index)))
print("SAINT rows: " + str(len(saint.index)))

# Merge the two dataframes on two columns:
# "Bait" in SAINT corresponds to "Experiment.ID" in CompPASS
# "Prey" in SAINT corresponds to "Prey" in CompPASS
merged = pd.merge(saint, compPASS, how="left", left_on=["Bait", "Prey"], right_on=["Experiment.ID", "Prey"])

# Print the number of rows in the merged dataframe
print("Merged rows: " + str(len(merged.index)))

# Rename Columns
merged.rename(columns={
    "Bait_y": "Bait.ID",
    "Prey": "Prey.ID",}, inplace=True)

# Drop Columns
merged.drop(columns=["Bait_x","boosted_by"], inplace=True)

# Write the merged dataframe to a file
merged.to_csv(f"{args.outputPath}/merged.csv", index=False)

# Stop the timer and print the elapsed time
elapsed_time = time.time() - start_time
print("Elapsed time: " + str(elapsed_time) + " seconds")


