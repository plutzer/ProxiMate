#!/usr/bin/env python
import sys
import os
import argparse
import subprocess
import shutil
import json
import pandas as pd
import aft_impute_saint
import refactored_aft
import one_component_aft
import compPASS_pval
import time
from log_config import get_logger, add_file_handler
from experimental_design import ExperimentalDesign
from bfdr_pool import recompute_bfdr

logger = get_logger(__name__)

SAINT_EXPRESS_INT_DIR = "/bin/SAINTexpress-int"
SAINT_EXPRESS_INT_DEFAULT_DIR = "/bin/SAINTexpress-int_default"
SAINT_EXPRESS_SPC_DIR = "/bin/SAINTexpress-spc"


def _build_saint_cmd(compress_n_rep, quant_type, imputation, prey_filename):
    """Build the SAINTexpress command vector. Filenames are relative to cwd."""
    if quant_type == "Spectral Counts":
        if imputation == "1":
            logger.warning("Imputation for spectral counts not yet implemented. Running SPC SAINT without imputation...")
        return [SAINT_EXPRESS_SPC_DIR,
                "-L", str(compress_n_rep),
                "filtered_interaction.txt", "prey.txt", "bait.txt"]

    if imputation in ("1", "2", "3"):
        return [SAINT_EXPRESS_INT_DIR,
                "-L", str(compress_n_rep),
                "filtered_interaction.txt", prey_filename, "bait.txt"]

    return [SAINT_EXPRESS_INT_DEFAULT_DIR,
            "-L", str(compress_n_rep),
            "filtered_interaction.txt", prey_filename, "bait.txt"]


def _choose_prey_filename(quant_type, imputation):
    """SAINT reads imputed_prey.txt for intensity+imputation runs; prey.txt otherwise."""
    if quant_type != "Spectral Counts" and imputation in ("1", "2", "3"):
        return "imputed_prey.txt"
    return "prey.txt"


def _run_saint(cwd, compress_n_rep, quant_type, imputation, prey_filename):
    """Invoke SAINTexpress in `cwd`. Writes list.txt there. Exits on failure."""
    saint_cmd = _build_saint_cmd(compress_n_rep, quant_type, imputation, prey_filename)
    logger.info("SAINTexpress command (cwd=%s): %s", cwd, " ".join(saint_cmd))
    p = subprocess.run(saint_cmd, cwd=cwd, capture_output=True, text=True)
    if p.returncode != 0:
        logger.error("SAINTexpress failed (exit code %d) in %s", p.returncode, cwd)
        if p.stdout:
            logger.error("SAINTexpress stdout:\n%s", p.stdout)
        if p.stderr:
            logger.error("SAINTexpress stderr:\n%s", p.stderr)
        sys.exit(1)
    logger.info("SAINTexpress completed successfully (cwd=%s)", cwd)
    if p.stdout:
        logger.debug("SAINTexpress stdout:\n%s", p.stdout)

    list_path = os.path.join(cwd, "list.txt")
    if not os.path.exists(list_path):
        logger.error("SAINTexpress did not produce list.txt at %s", list_path)
        sys.exit(1)


def _build_group_saint_inputs(src_dir, group_dir, ed, group, use_imputed_prey):
    """Build per-group SAINT inputs by filtering bait/interaction files to the group's
    experiments and copying the global prey file(s) unchanged."""
    tests, ctrls = ed.get_experiments_for_group(group)
    # Force str — ED values round-trip as CSV so they are strings; if Experiment
    # Names look numeric, pandas will infer int64 on the files below and the
    # isin() match would silently return zero rows.
    allowed = {str(e.attributes["Experiment Name"]) for e in (tests + ctrls)}
    logger.info("Group %s: %d test + %d control experiments; allowed names: %s",
                group, len(tests), len(ctrls), sorted(allowed))

    bait_in = pd.read_csv(
        os.path.join(src_dir, "bait.txt"),
        sep="\t", header=None, names=["ExperimentName", "Bait", "Type"],
        usecols=[0, 1, 2],
        dtype={"ExperimentName": str, "Bait": str, "Type": str},
    )
    bait_out = bait_in[bait_in["ExperimentName"].isin(allowed)]
    logger.info("Group %s bait.txt: %d rows in, %d rows after filter",
                group, len(bait_in), len(bait_out))
    if bait_out.empty and not bait_in.empty:
        logger.error("Group %s bait.txt filter produced zero rows. "
                     "Sample ExperimentName values in bait.txt (dtype=%s): %s",
                     group, bait_in["ExperimentName"].dtype,
                     bait_in["ExperimentName"].head(5).tolist())
    bait_out.to_csv(
        os.path.join(group_dir, "bait.txt"),
        sep="\t", header=False, index=False,
    )

    inter_in = pd.read_csv(
        os.path.join(src_dir, "filtered_interaction.txt"),
        sep="\t", header=None,
        names=["ExperimentName", "Bait", "Prey", "Intensity"],
        usecols=[0, 1, 2, 3],
        dtype={"ExperimentName": str, "Bait": str, "Prey": str},
    )
    inter_out = inter_in[inter_in["ExperimentName"].isin(allowed)]
    logger.info("Group %s filtered_interaction.txt: %d rows in, %d rows after filter",
                group, len(inter_in), len(inter_out))
    if inter_out.empty and not inter_in.empty:
        logger.error("Group %s filtered_interaction.txt filter produced zero rows. "
                     "Sample ExperimentName values in file (dtype=%s): %s",
                     group, inter_in["ExperimentName"].dtype,
                     inter_in["ExperimentName"].head(5).tolist())
    inter_out.to_csv(
        os.path.join(group_dir, "filtered_interaction.txt"),
        sep="\t", header=False, index=False,
    )

    # Prey files stay global — SAINT's prey universe must match across groups
    # so BFDR re-pooling is well-defined.
    shutil.copy(os.path.join(src_dir, "prey.txt"), group_dir)
    if use_imputed_prey:
        shutil.copy(os.path.join(src_dir, "imputed_prey.txt"), group_dir)


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
                        help="0 (none), 1 (prey-specific AFT), 2 (refactored AFT), 3 (one-component AFT)",
                        default="0")

    # Argument for quantification type
    parser.add_argument("--quantType",
                        help="SAINT: quantification type (intensity, spc, LFQ)",
                        default="Intensity")

    # pi estimation options (imputation=2 only)
    parser.add_argument("--pi-method", dest="pi_method",
                        choices=["weighted_average", "single_bait"],
                        default="weighted_average",
                        help="Refactored AFT (imputation=2): how to estimate pi from controls.")
    parser.add_argument("--pi-bait", dest="pi_bait", default=None,
                        help="Required when --pi-method=single_bait: control Bait name to fit.")

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
            if args.pi_method == "single_bait" and not args.pi_bait:
                logger.error("--pi-method=single_bait requires --pi-bait")
                sys.exit(1)
            logger.info("Running refactored AFT imputation (pi_method=%s, pi_bait=%s)...",
                        args.pi_method, args.pi_bait)
            refactored_aft.filter_impute(
                f"{args.scoreInputs}/prey.txt",
                f"{args.scoreInputs}/interaction.txt",
                f"{args.scoreInputs}/",
                args.experimentalDesign,
                impute=True,
                pi_method=args.pi_method,
                pi_bait=args.pi_bait)
        elif args.imputation == "3":
            logger.info("Running one-component AFT imputation...")
            one_component_aft.filter_impute(
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

    # Parse ED once to decide between legacy and grouped SAINT flow
    try:
        ed = ExperimentalDesign(args.experimentalDesign)
    except Exception:
        logger.exception("Failed to parse experimental design for grouping decision")
        sys.exit(1)

    prey_filename = _choose_prey_filename(args.quantType, args.imputation)
    use_imputed_prey = (prey_filename == "imputed_prey.txt")

    if ed.is_grouped():
        groups = ed.get_groups()
        logger.info("Running SAINTexpress in grouped mode (quantType=%s, imputation=%s, %d groups: %s)",
                    args.quantType, args.imputation, len(groups), groups)

        groups_root = os.path.join(args.outputPath, "groups")
        os.makedirs(groups_root, exist_ok=True)

        per_group_dfs = []
        manifest_entries = []
        for group in groups:
            group_dir = os.path.join(groups_root, str(group))
            os.makedirs(group_dir, exist_ok=True)
            _build_group_saint_inputs(args.scoreInputs, group_dir, ed, group, use_imputed_prey)
            _run_saint(group_dir, args.compress_n_rep, args.quantType,
                       args.imputation, prey_filename)

            df = pd.read_csv(os.path.join(group_dir, "list.txt"), sep="\t")
            df["source_group"] = group
            per_group_dfs.append(df)

            tests, ctrls = ed.get_experiments_for_group(group)
            manifest_entries.append({
                "group": group,
                "n_test": len(tests),
                "n_control": len(ctrls),
                "list_txt": f"groups/{group}/list.txt",
            })

        saint_merged = pd.concat(per_group_dfs, ignore_index=True)
        saint_merged = recompute_bfdr(saint_merged)
        saint_merged.to_csv(os.path.join(args.scoreInputs, "list.txt"), sep="\t", index=False)

        with open(os.path.join(groups_root, "manifest.json"), "w") as f:
            json.dump({
                "run_mode": "grouped",
                "groups": manifest_entries,
                "imputation": args.imputation,
                "quant_type": args.quantType,
                "bfdr_repooled": True,
            }, f, indent=2)
        logger.info("Grouped SAINT complete; wrote concatenated list.txt (%d rows) and manifest.json",
                    len(saint_merged))
    else:
        logger.info("Running SAINTexpress (quantType=%s, imputation=%s)...", args.quantType, args.imputation)
        _run_saint(args.scoreInputs, args.compress_n_rep, args.quantType,
                   args.imputation, prey_filename)

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

        # Keep output schema stable: legacy runs get a NA source_group column
        if "source_group" not in merged.columns:
            merged["source_group"] = pd.NA

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
