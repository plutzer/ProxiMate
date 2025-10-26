import pandas as pd
import numpy as np
import argparse
import os
from experimental_design import ExperimentalDesign
from protein_groups import ProteinGroups
import shutil
import re
import tempfile

def validate_name(s: str, datasets):
    """
    Check that `s` is non‐empty, contains no whitespace, and only alphanumeric or underscore characters.
    Returns:
      0 if the name is valid,
      otherwise a string describing the first validation error.
    """
    if not s:
        return "Error: Dataset name cannot be empty"
    if re.search(r"\s", s):
        return "Error: Dataset name cannot contain spaces"
    # allow only letters, digits, and underscore
    if re.search(r"[^A-Za-z0-9_]", s):
        return "Error: Dataset name can only contain letters, numbers, and underscores"
    if s in datasets:
        return f"Error: Dataset name '{s}' already exists in the datasets"
    return 0

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

def parse_from_saint(bait_df, preyfile, interactionfile, outputPath):

    # bait = pd.read_csv(baitfile, sep="\t", header=None, names=["Experiment Name", "Bait", "Type", ])
    prey = pd.read_csv(preyfile, sep="\t", header=None, names=["Prey", "Prey.Name"])
    interaction = pd.read_csv(interactionfile, sep="\t", header=None, names=["Experiment.ID", "Bait", "Prey", "Spectral.Count"])
    print('Parsing SAINT files...')

    # Recreate the experimental design file
    new_ed = bait_df.copy()[["Experiment Name", "Type", "Bait", "Bait ID"]]
    # Add replicate numbers
    new_ed["Replicate"] = new_ed.groupby("Bait").cumcount() + 1
    # new_ed["Bait ID"] = "None"
    new_ed.to_csv(os.path.join(outputPath, "ED.csv"), index=False)

    # Recreate the proteinGroups file
    # Apparently, I might not need to do this...

    # Recreate the to_compass file
    compass = interaction.copy()
    compass = pd.merge(compass, new_ed.drop(columns=["Bait"]), left_on="Experiment.ID", right_on="Experiment Name", how="left")
    compass = pd.merge(compass, prey, on="Prey", how="left")

    to_compass = compass[["Bait", "Replicate", "Bait ID", "Prey", "Prey.Name", "Spectral.Count"]].rename(columns={"Bait ID": "Bait", "Bait": "Experiment.ID"})
    to_compass.to_csv(os.path.join(outputPath, "to_CompPASS.csv"), index=False)

    # Return information needed for the GUI Datasets table
    n_expts = len(new_ed["Experiment Name"].unique())
    n_ctrls = len(new_ed[new_ed["Type"] == "C"]["Experiment Name"].unique())
    return n_expts, n_ctrls

def convert_diann_to_maxquant_format(diann_file, experimental_design):
    """
    Convert DIA-NN report.pg_matrix.tsv to a MaxQuant-like proteinGroups format.

    :param diann_file: path to DIA-NN report.pg_matrix.tsv file
    :param experimental_design: ExperimentalDesign object
    :return: DataFrame in MaxQuant-like format
    """
    print("Converting DIA-NN format to MaxQuant-compatible format...")

    # Read DIA-NN file
    diann_data = pd.read_csv(diann_file, sep="\t")

    # Create a new DataFrame with MaxQuant-like columns
    mq_data = pd.DataFrame()

    # Map DIA-NN columns to MaxQuant columns
    mq_data["Majority protein IDs"] = diann_data["Protein.Group"]
    mq_data["Protein names"] = diann_data["Protein.Names"] if "Protein.Names" in diann_data.columns else ""
    mq_data["Gene names"] = diann_data["Genes"] if "Genes" in diann_data.columns else ""

    # Add fake filtering columns (DIA-NN doesn't have these, so set all to "-")
    mq_data["Reverse"] = "-"
    mq_data["Only identified by site"] = "-"
    mq_data["Potential contaminant"] = "-"

    # DIA-NN doesn't have sequence length in pg_matrix, but it's needed for spectral counts
    # For now, set to 1 (this won't affect intensity-based quantification)
    mq_data["Sequence length"] = 1

    # Copy over the intensity columns - rename them to match MaxQuant format
    # DIA-NN columns are the raw file paths, need to add "Intensity " prefix
    for col in diann_data.columns:
        # Skip the non-quantification columns
        if col in ["Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
                   "N.Sequences", "N.Proteotypic.Sequences"]:
            continue

        # Check if this column is in the experimental design
        if col in experimental_design.name2experiment:
            # Add with "Intensity " prefix to match MaxQuant format
            mq_data[f"Intensity {col}"] = diann_data[col].fillna(0)

    print(f"Converted {len(mq_data)} protein groups from DIA-NN format")

    return mq_data

def parse_diann(diannMatrix, experimentalDesign, quantType, outputPath):
    """
    Parse DIA-NN report.pg_matrix.tsv file and experimental design file.

    :param diannMatrix: path to DIA-NN report.pg_matrix.tsv file
    :param experimentalDesign: path to experimental design file
    :param quantType: quantification type (should be "Intensity" for DIA-NN)
    :param outputPath: path for the output directory
    :return: tuple of (num_experiments, num_controls)
    """
    print("Parsing DIA-NN data...")
    print("Experimental Design file: " + experimentalDesign)
    print("DIA-NN Matrix file: " + diannMatrix)
    print("Quantification type: " + quantType)

    # Check to see if the output directory exists. If not, create it.
    if not os.path.exists(outputPath):
        print("Creating output directory: " + outputPath)
        os.makedirs(outputPath)

    # Copy the experimental design file to the output directory
    shutil.copy(experimentalDesign, f"{outputPath}/ED.csv")

    # Parse experimental design
    experimental_design = ExperimentalDesign(experimentalDesign)

    num_expts = experimental_design.num_experiments
    num_ctrls = experimental_design.num_controls

    # Convert DIA-NN format to MaxQuant-like format
    mq_format_data = convert_diann_to_maxquant_format(diannMatrix, experimental_design)

    # Save converted data to a temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, newline='') as tmp_file:
        tmp_path = tmp_file.name
        mq_format_data.to_csv(tmp_file, sep="\t", index=False)

    # Copy the temporary file to the output directory as proteinGroups.txt
    shutil.copy(tmp_path, f"{outputPath}/proteinGroups.txt")

    try:
        # Process using the standard ProteinGroups class
        # For DIA-NN, we always use "Intensity" as the quantification type
        protein_groups = ProteinGroups(experimental_design, tmp_path,
                                       "Intensity", "Intensity")

        protein_groups.to_SAINT(outputPath)
        protein_groups.to_CompPASS(outputPath)
    finally:
        # Clean up temporary file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    return num_expts, num_ctrls

def main():

    description = "This is the entry point to the program. It will execute the requested tasks."

    # initialize the parser
    parser = argparse.ArgumentParser(description=description)

    # Add arguments for the annotator
    parser.add_argument("--proteinGroups",
                        help="path to MaxQuant ProteinGroups.txt",
                        default=None)

    parser.add_argument("--diannMatrix",
                        help="path to DIA-NN report.pg_matrix.tsv file",
                        default=None)

    parser.add_argument("--experimentalDesign",
                        help="path to experimental design file",
                        default=None)

    parser.add_argument("--outputPath",
                        help="path for the output directory. If it already exists it will be overwritten.",
                        default="/srv/shiny-server/myapp/score_inputs")

    # Argument for quantification type
    parser.add_argument("--quantType",
                        help="SAINT: quantification type (Intensity, Spectral Counts, LFQ)",
                        default="Intensity")

    args = parser.parse_args()

    # MAIN

    if args.diannMatrix is not None:
        # DIA-NN input mode
        if args.experimentalDesign is not None:
            parse_diann(args.diannMatrix, args.experimentalDesign, args.quantType, args.outputPath)
        else:
            print("Error: No experimental design file provided.")
    elif args.proteinGroups is not None:
        # MaxQuant input mode
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
        print("Error: No input file provided. Please specify either --proteinGroups or --diannMatrix.")

if __name__ == "__main__":
    # basedir = "C:/Users/isaac/Work/Ilah_testdata/SAINT"
    # # For testing parse from saint
    # bait_df = pd.read_csv(basedir + "/bait.txt", sep="\t", header=None, names=["Experiment Name", "Bait", "Type", "Bait ID"])
    # parse_from_saint(bait_df, basedir + "/prey.txt", basedir + "/interaction.txt", basedir)

    main()