import pandas as pd
import numpy as np
import argparse
import os
import sys
from experimental_design import ExperimentalDesign
from protein_groups import ProteinGroups
import shutil
import re
import tempfile
from ed_validation import validate_maxquant_inputs, validate_diann_inputs, validate_fragpipe_inputs
from ed_exceptions import ProxiMateError
from log_config import get_logger, add_file_handler

logger = get_logger(__name__)

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

    # VALIDATE INPUTS FIRST
    try:
        ed_df, pg_df = validate_maxquant_inputs(experimentalDesign, proteinGroups, quantType)
    except ProxiMateError:
        # Re-raise validation errors to be caught by app.py
        raise

    # Check to see if the output directory exists. If not, create it.
    if not os.path.exists(outputPath):
        logger.info("Creating output directory: %s", outputPath)
        os.makedirs(outputPath)

    add_file_handler(os.path.join(outputPath, "proximate.log"))

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
    # VALIDATE INPUTS FIRST
    try:
        ed_df, diann_df = validate_diann_inputs(experimentalDesign, diannMatrix)
    except ProxiMateError:
        # Re-raise validation errors to be caught by app.py
        raise

    # Check to see if the output directory exists. If not, create it.
    if not os.path.exists(outputPath):
        logger.info("Creating output directory: %s", outputPath)
        os.makedirs(outputPath)

    add_file_handler(os.path.join(outputPath, "proximate.log"))

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

def convert_fragpipe_to_maxquant_format(fp_file, experimental_design, quant_type):
    """
    Convert FragPipe combined_protein.tsv to a MaxQuant-like proteinGroups format.

    :param fp_file: path to FragPipe combined_protein.tsv file
    :param experimental_design: ExperimentalDesign object
    :param quant_type: quantification type (Intensity, LFQ, or Spectral Counts)
    :return: DataFrame in MaxQuant-like format
    """
    fp_data = pd.read_csv(fp_file, sep="\t")

    # Map quant_type to FragPipe column suffix and MaxQuant column prefix
    QUANT_MAP = {
        "Intensity":       (" Intensity",            "Intensity "),
        "LFQ":             (" MaxLFQ Intensity",     "LFQ intensity "),
        "Spectral Counts": (" Total Spectral Count", "MS/MS count "),
    }
    fp_suffix, mq_prefix = QUANT_MAP[quant_type]

    mq_data = pd.DataFrame()

    # Map FragPipe columns to MaxQuant columns
    mq_data["Majority protein IDs"] = fp_data["Protein ID"]
    mq_data["Protein names"] = fp_data["Description"] if "Description" in fp_data.columns else ""
    mq_data["Gene names"] = fp_data["Gene"].fillna("")

    # Filter columns — FragPipe uses contam_ prefix instead of separate columns
    mq_data["Reverse"] = "-"
    mq_data["Only identified by site"] = "-"
    mq_data["Potential contaminant"] = np.where(
        fp_data["Protein"].str.startswith("contam_"), "+", "-"
    )

    # FragPipe has actual protein length (DIA-NN lacks this and sets to 1)
    mq_data["Sequence length"] = fp_data["Protein Length"]

    # Copy and rename quantification columns
    for col in fp_data.columns:
        if col.endswith(fp_suffix):
            # Disambiguate: when in Intensity mode, skip MaxLFQ Intensity columns
            if quant_type == "Intensity" and col.endswith(" MaxLFQ Intensity"):
                continue
            sample_name = col[:-len(fp_suffix)]
            if sample_name in experimental_design.name2experiment:
                mq_data[f"{mq_prefix}{sample_name}"] = fp_data[col].fillna(0)

    return mq_data

def parse_fragpipe(fp_file, experimentalDesign, quantType, outputPath):
    """
    Parse FragPipe combined_protein.tsv file and experimental design file.

    :param fp_file: path to FragPipe combined_protein.tsv file
    :param experimentalDesign: path to experimental design file
    :param quantType: quantification type (Intensity, LFQ, or Spectral Counts)
    :param outputPath: path for the output directory
    :return: tuple of (num_experiments, num_controls)
    """
    # VALIDATE INPUTS FIRST
    try:
        ed_df, fp_df = validate_fragpipe_inputs(experimentalDesign, fp_file)
    except ProxiMateError:
        raise

    if not os.path.exists(outputPath):
        logger.info("Creating output directory: %s", outputPath)
        os.makedirs(outputPath)

    add_file_handler(os.path.join(outputPath, "proximate.log"))

    shutil.copy(experimentalDesign, f"{outputPath}/ED.csv")

    experimental_design = ExperimentalDesign(experimentalDesign)
    num_expts = experimental_design.num_experiments
    num_ctrls = experimental_design.num_controls

    mq_format_data = convert_fragpipe_to_maxquant_format(fp_file, experimental_design, quantType)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, newline='') as tmp_file:
        tmp_path = tmp_file.name
        mq_format_data.to_csv(tmp_file, sep="\t", index=False)

    shutil.copy(tmp_path, f"{outputPath}/proteinGroups.txt")

    try:
        protein_groups = ProteinGroups(experimental_design, tmp_path,
                                       quantType, quantType)
        protein_groups.to_SAINT(outputPath)
        protein_groups.to_CompPASS(outputPath)
    finally:
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

    parser.add_argument("--fragpipeFile",
                        help="path to FragPipe combined_protein.tsv file",
                        default=None)

    # SAINT format arguments
    parser.add_argument("--bait",
                        help="path to SAINT bait.txt file",
                        default=None)

    parser.add_argument("--prey",
                        help="path to SAINT prey.txt file",
                        default=None)

    parser.add_argument("--interaction",
                        help="path to SAINT interaction.txt file",
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
            logger.error("No experimental design file provided.")
    elif args.fragpipeFile is not None:
        # FragPipe input mode
        if args.experimentalDesign is not None:
            parse_fragpipe(args.fragpipeFile, args.experimentalDesign, args.quantType, args.outputPath)
        else:
            logger.error("No experimental design file provided.")
    elif args.proteinGroups is not None:
        # MaxQuant input mode
        if args.experimentalDesign is not None:
            # Check to see if the output directory exists. If not, create it.
            if not os.path.exists(args.outputPath):
                os.makedirs(args.outputPath)

            # Copy the experimental design file to the output directory
            shutil.copy(args.experimentalDesign, os.path.join(args.outputPath, "ED.csv"))

            logger.info("Parsing Experimental Design file: %s", args.experimentalDesign)

            # parse experimental design
            experimental_design = ExperimentalDesign(args.experimentalDesign)

            logger.info("Parsing Protein Groups file: %s", args.proteinGroups)
            logger.info("Quantification type: %s", args.quantType)

            # process MaxQuant proteinGroups
            protein_groups = ProteinGroups(experimental_design, args.proteinGroups,
                                        args.quantType, args.quantType)

            protein_groups.to_SAINT(args.outputPath)
            protein_groups.to_CompPASS(args.outputPath)
        else:
            logger.error("No experimental design file provided.")
    elif args.bait is not None:
        # SAINT input mode
        if args.prey is None or args.interaction is None:
            logger.error("SAINT format requires --bait, --prey, and --interaction files.")
            sys.exit(1)

        if not os.path.exists(args.outputPath):
            os.makedirs(args.outputPath)

        # Read bait.txt into a DataFrame matching GUI expectations
        bait_df = pd.read_csv(args.bait, sep="\t", header=None,
                              names=["Experiment Name", "Bait", "Type"])
        bait_df["Bait ID"] = "None"

        logger.info("Parsing SAINT input files...")
        parse_from_saint(bait_df, args.prey, args.interaction, args.outputPath)

        # Copy SAINT files to output directory (mirrors GUI behavior)
        shutil.copy(args.bait, os.path.join(args.outputPath, "bait.txt"))
        shutil.copy(args.prey, os.path.join(args.outputPath, "prey.txt"))
        shutil.copy(args.interaction, os.path.join(args.outputPath, "interaction.txt"))
    else:
        logger.error("No input file provided. Please specify --proteinGroups, --diannMatrix, --fragpipeFile, or --bait/--prey/--interaction.")

if __name__ == "__main__":
    main()