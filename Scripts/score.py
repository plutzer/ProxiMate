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
from log_config import get_logger, add_file_handler

logger = get_logger(__name__)

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

    # Check to see if the output directory exists. If not, create it.
    if not os.path.exists(args.outputPath):
        os.makedirs(args.outputPath)

    add_file_handler(os.path.join(args.outputPath, "proximate.log"))

    # Validate required input files exist
    for required_file in ["prey.txt", "interaction.txt", "bait.txt", "to_CompPASS.csv"]:
        fpath = os.path.join(args.scoreInputs, required_file)
        if not os.path.exists(fpath):
            logger.error(f"Required input file not found: {fpath}")
            sys.exit(1)

    # Run the imputation
    try:
        if args.imputation == "2":
            logger.info("Running refactored AFT imputation...")
            refactored_aft.filter_impute(
                f"{args.scoreInputs}/prey.txt",
                f"{args.scoreInputs}/interaction.txt",
                f"{args.scoreInputs}/",
                args.experimentalDesign,
                impute=True)
        else:
            logger.info("Running imputation filter (impute=%s)...", bool(int(args.imputation)))
            aft_impute_saint.filter_impute(
                f"{args.scoreInputs}/prey.txt",
                f"{args.scoreInputs}/interaction.txt",
                f"{args.scoreInputs}/",
                args.experimentalDesign,
                impute=bool(int(args.imputation)))
    except Exception:
        logger.exception("Imputation failed")
        sys.exit(1)

    # Verify imputation produced filtered_interaction.txt
    filtered_interaction = os.path.join(args.scoreInputs, "filtered_interaction.txt")
    if not os.path.exists(filtered_interaction):
        logger.error("Imputation did not produce filtered_interaction.txt")
        sys.exit(1)

    ########################################################################################################################

    # Run SAINT
    logger.info("Running SAINTexpress (quantType=%s, imputation=%s)...", args.quantType, args.imputation)
    if args.quantType == "Spectral Counts":
        if args.imputation == "1":
            logger.warning("Imputation for spectral counts not yet implemented. Running SPC SAINT without imputation...")
        saint_cmd = [SAINT_EXPRESS_SPC_DIR,
                     "-L", str(args.compress_n_rep),
                     "filtered_interaction.txt", "prey.txt", "bait.txt"]
    else:
        if args.imputation in ("1", "2"):
            saint_cmd = [SAINT_EXPRESS_INT_DIR,
                         "-L", str(args.compress_n_rep),
                         "filtered_interaction.txt", "imputed_prey.txt", "bait.txt"]
        else:
            saint_cmd = [SAINT_EXPRESS_INT_DEFAULT_DIR,
                         "-L", str(args.compress_n_rep),
                         "filtered_interaction.txt", "prey.txt", "bait.txt"]

    logger.info("SAINTexpress command: %s", " ".join(saint_cmd))
    p = subprocess.run(saint_cmd, cwd=args.scoreInputs,
                       capture_output=True, text=True)

    if p.returncode != 0:
        logger.error("SAINTexpress failed (exit code %d)", p.returncode)
        if p.stdout:
            logger.error("SAINTexpress stdout:\n%s", p.stdout)
        if p.stderr:
            logger.error("SAINTexpress stderr:\n%s", p.stderr)
        sys.exit(1)
    else:
        logger.info("SAINTexpress completed successfully")
        if p.stdout:
            logger.debug("SAINTexpress stdout:\n%s", p.stdout)

    # Verify SAINT produced list.txt
    saint_output = os.path.join(args.scoreInputs, "list.txt")
    if not os.path.exists(saint_output):
        logger.error("SAINTexpress did not produce list.txt at %s", saint_output)
        sys.exit(1)

    ########################################################################################################################
    ########################################################################################################################
    # Run CompPASS
    logger.info("Running CompPASS...")

    try:
        comp_input = pd.read_csv(f"{args.outputPath}/to_CompPASS.csv", sep=',', index_col=0).astype({'Prey': str, 'Bait': str})
        compPass_result = compPASS_pval.score_compPass(comp_input, 0.98, args.n_iterations)
        compPass_result.to_csv(f"{args.outputPath}/compPASS.csv", index=False)
        logger.info("CompPASS completed successfully")
    except Exception:
        logger.exception("CompPASS scoring failed")
        sys.exit(1)

    ########################################################################################################################
    ########################################################################################################################

    # Merge CompPASS and SAINT outputs
    logger.info("Merging CompPASS and SAINT outputs...")
    try:
        saint = pd.read_csv(f"{args.scoreInputs}/list.txt", sep="\t")
        logger.info("SAINT output: %d rows, columns: %s", len(saint), list(saint.columns))

        merged = pd.merge(saint, compPass_result, how="left", left_on=["Bait", "Prey"], right_on=["Experiment.ID", "Prey"])

        # Rename Columns
        merged.rename(columns={
            "Bait_y": "Bait.ID",
            "Prey": "Prey.ID",}, inplace=True)

        # Drop Columns
        merged.drop(columns=["Bait_x","boosted_by"], inplace=True)

        # Write the merged dataframe to a file
        merged.to_csv(f"{args.outputPath}/merged.csv", index=False)
        logger.info("Merged output written to %s/merged.csv (%d rows)", args.outputPath, len(merged))
    except Exception:
        logger.exception("Merging CompPASS and SAINT outputs failed")
        sys.exit(1)

    elapsed = time.time() - start_time
    logger.info("Scoring completed in %.1f seconds", elapsed)

if __name__ == "__main__":
    main()
