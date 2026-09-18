import pandas as pd
import numpy as np
import argparse
import functools
import inspect
import os
import sys
from experimental_design import ExperimentalDesign
from protein_groups import ProteinGroups
import shutil
import re
import tempfile
import provenance
from ed_validation import (validate_maxquant_inputs, validate_diann_inputs, validate_fragpipe_inputs,
                           validate_msstats_inputs, validate_pioneer_inputs, PIONEER_METADATA_COLUMNS)
from ed_exceptions import ProxiMateError
from log_config import get_logger, add_file_handler, dataset_log

logger = get_logger(__name__)

# Files a parse entry point produces for the scoring stage to consume.
PARSE_OUTPUTS = ("ED.csv", "prey.txt", "bait.txt", "interaction.txt", "to_CompPASS.csv")

# Parameters of the parse entry points that name a file. Everything else is a
# value (quantType is "LFQ", not a path) and is recorded as a parameter instead.
# Listed explicitly rather than inferred from the value, so that a genuinely
# missing input is still fingerprinted as missing rather than quietly skipped.
PARSE_FILE_PARAMS = frozenset({
    "proteinGroups", "experimentalDesign", "diannMatrix", "pioneerMatrix", "msstatsFile",
    "fp_file", "preyfile", "interactionfile",
})


def _parse_stage(fn):
    """Record a parse entry point in the dataset's run manifest.

    Applied as a decorator so the entry points' bodies stay as they are: the
    output directory and the input files are read off the call signature, which
    every entry point shares in shape.  A parse that fails validation still
    leaves a manifest saying so.
    """
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        bound = inspect.signature(fn).bind(*args, **kwargs)
        bound.apply_defaults()
        output_path = bound.arguments["outputPath"]

        params, inputs = {}, {}
        for name, value in bound.arguments.items():
            if name == "outputPath":
                continue
            if name in PARSE_FILE_PARAMS:
                inputs[name] = value
                params[name] = str(value)
            elif isinstance(value, pd.DataFrame):
                # Callers pass frames as well as paths; a frame's shape is the
                # only part of it worth putting in a manifest.
                params[name] = f"DataFrame{value.shape}"
            else:
                params[name] = value

        # Attach the dataset log around the whole call, not just the part after
        # validation: an entry point that rejects its inputs never reaches its
        # own add_file_handler, and that failure is the one most worth reading.
        with dataset_log(output_path), \
                provenance.stage(output_path, "parse",
                                 entrypoint=f"parse.{fn.__name__}",
                                 params=params) as record:
            for role, path in inputs.items():
                record.add_input(path, role=role)

            result = fn(*args, **kwargs)

            num_expts, num_ctrls = result
            record.metric("n_experiments", num_expts)
            record.metric("n_controls", num_ctrls)
            for produced in PARSE_OUTPUTS:
                path = os.path.join(str(output_path), produced)
                if os.path.exists(path):
                    record.add_output(path)
            return result
    return wrapper

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

@_parse_stage
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

@_parse_stage
def parse_from_saint(bait_df, preyfile, interactionfile, outputPath):

    if not os.path.exists(outputPath):
        logger.info("Creating output directory: %s", outputPath)
        os.makedirs(outputPath)

    add_file_handler(os.path.join(outputPath, "proximate.log"))

    logger.info("Parsing SAINT inputs: prey=%s interaction=%s", preyfile, interactionfile)

    # bait = pd.read_csv(baitfile, sep="\t", header=None, names=["Experiment Name", "Bait", "Type", ])
    prey = pd.read_csv(preyfile, sep="\t", header=None, names=["Prey", "Prey.Name"])
    interaction = pd.read_csv(interactionfile, sep="\t", header=None, names=["Experiment.ID", "Bait", "Prey", "Spectral.Count"])
    logger.info("Read %d preys and %d interaction rows across %d experiments",
                len(prey), len(interaction), interaction["Experiment.ID"].nunique())
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

    # Return information needed for the GUI Datasets table.  Counts exclude
    # controls, matching ExperimentalDesign.num_experiments, which is what the
    # other parse entry points return and what scoring reports for the same
    # dataset.
    n_ctrls = len(new_ed[new_ed["Type"] == "C"]["Experiment Name"].unique())
    n_expts = len(new_ed[new_ed["Type"] != "C"]["Experiment Name"].unique())
    logger.info("Parsed SAINT inputs: %d experiments (%d controls), %d CompPASS rows",
                n_expts, n_ctrls, len(to_compass))

    # bait.txt has no column for the bait's protein ID, so anything keyed on it
    # cannot be derived later.  Left unsaid, the resulting empty annotations
    # read as negative findings rather than as an absent input.
    if not new_ed["Bait ID"].astype(str).str.strip().replace("None", "").any():
        logger.warning(
            "SAINT bait.txt carries no bait protein IDs, so BioGRID and "
            "self-interaction annotation cannot be derived for this dataset.")

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

@_parse_stage
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


def convert_pioneer_to_maxquant_format(pioneer_file, experimental_design):
    """
    Convert Pioneer's protein_groups_wide.tsv to a MaxQuant-like proteinGroups format.

    Decoy groups (``target`` false) and entrapment groups (``entrap_id`` non-zero)
    are not real proteins and are dropped. Run columns are matched to the design
    by exact name.

    :param pioneer_file: path to protein_groups_wide.tsv
    :param experimental_design: ExperimentalDesign object
    :return: DataFrame in MaxQuant-like format
    """
    data = pd.read_csv(pioneer_file, sep="\t")

    real = pd.Series(True, index=data.index)
    if "target" in data.columns:
        real &= data["target"].astype(bool)
    if "entrap_id" in data.columns:
        real &= data["entrap_id"].fillna(0) == 0
    if (~real).any():
        logger.info("Dropping %d decoy/entrapment protein group(s) from the Pioneer table", (~real).sum())
        data = data[real].reset_index(drop=True)

    mq_data = pd.DataFrame()
    mq_data["Majority protein IDs"] = data["protein"]
    mq_data["Protein names"] = data["protein_names"].fillna("") if "protein_names" in data.columns else ""
    mq_data["Gene names"] = data["gene_names"].fillna("") if "gene_names" in data.columns else ""

    # No reverse, site-only or contaminant flags survive the decoy drop above.
    mq_data["Reverse"] = "-"
    mq_data["Only identified by site"] = "-"
    mq_data["Potential contaminant"] = "-"

    # The wide table carries no protein length; nothing reads it for intensity input.
    mq_data["Sequence length"] = 1

    for col in data.columns:
        if col in PIONEER_METADATA_COLUMNS:
            continue
        if col in experimental_design.name2experiment:
            mq_data[f"Intensity {col}"] = data[col].fillna(0)

    return mq_data


@_parse_stage
def parse_pioneer(pioneerMatrix, experimentalDesign, quantType, outputPath):
    """
    Parse Pioneer's protein_groups_wide.tsv and an experimental design file.

    :param pioneerMatrix: path to protein_groups_wide.tsv
    :param experimentalDesign: path to experimental design file
    :param quantType: recorded only; Pioneer input is always parsed as Intensity
    :param outputPath: path for the output directory
    :return: tuple of (num_experiments, num_controls)
    """
    validate_pioneer_inputs(experimentalDesign, pioneerMatrix)

    if not os.path.exists(outputPath):
        logger.info("Creating output directory: %s", outputPath)
        os.makedirs(outputPath)

    add_file_handler(os.path.join(outputPath, "proximate.log"))

    shutil.copy(experimentalDesign, f"{outputPath}/ED.csv")

    experimental_design = ExperimentalDesign(experimentalDesign)
    num_expts = experimental_design.num_experiments
    num_ctrls = experimental_design.num_controls

    mq_format_data = convert_pioneer_to_maxquant_format(pioneerMatrix, experimental_design)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, newline='') as tmp_file:
        tmp_path = tmp_file.name
        mq_format_data.to_csv(tmp_file, sep="\t", index=False)

    shutil.copy(tmp_path, f"{outputPath}/proteinGroups.txt")

    try:
        protein_groups = ProteinGroups(experimental_design, tmp_path,
                                       "Intensity", "Intensity")
        protein_groups.to_SAINT(outputPath)
        protein_groups.to_CompPASS(outputPath)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    return num_expts, num_ctrls


def _read_msstats_csv(msstats_file):
    """Read an MSstats ProteinLevelData.csv with the standard encoding fallback chain."""
    last_err = None
    for encoding in ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252']:
        try:
            df = pd.read_csv(msstats_file, encoding=encoding)
            if not df.empty:
                return df
        except UnicodeDecodeError as e:
            last_err = e
            continue
    raise ValueError(f"Could not read MSstats file {msstats_file}: {last_err}")


def convert_msstats_to_maxquant_format(msstats_file, experimental_design):
    """
    Convert an MSstats::dataProcess() ProteinLevelData.csv to a MaxQuant-like
    proteinGroups format consumable by ProteinGroups.

    LogIntensities are log2(normalized_abundance) — back-transform to linear
    (2 ** LogIntensities) so the rest of the pipeline (which assumes raw
    intensities) sees a normalized intensity matrix. SAINT will re-log
    internally; that round-trip is a no-op but preserves the normalization.
    """
    df = _read_msstats_csv(msstats_file)

    # Label-free filter
    n_heavy = int((df["LABEL"] == "H").sum())
    if n_heavy > 0:
        logger.warning("Dropping %d SILAC heavy (LABEL=='H') rows; only label-free is supported.", n_heavy)
    df = df[df["LABEL"] == "L"].copy()

    # Drop unimputed missing values (only present when MBimpute=FALSE)
    n_na = int(df["LogIntensities"].isna().sum())
    if n_na > 0:
        logger.warning("Dropping %d rows with NA LogIntensities (run dataProcess with MBimpute=TRUE to impute).", n_na)
    df = df.dropna(subset=["LogIntensities"]).copy()

    # Collapse any accidental duplicates on (Protein, originalRUN)
    dup_mask = df.duplicated(subset=["Protein", "originalRUN"], keep=False)
    if dup_mask.any():
        n_dup_groups = df.loc[dup_mask, ["Protein", "originalRUN"]].drop_duplicates().shape[0]
        logger.warning("Found %d (Protein, originalRUN) groups with multiple rows; averaging LogIntensities.",
                       n_dup_groups)
        df = (df.groupby(["Protein", "originalRUN"], as_index=False)
                .agg({"LogIntensities": "mean"}))

    # Back-transform log2 → linear normalized intensity
    df["Intensity"] = np.power(2.0, df["LogIntensities"].astype(float))

    # Wide pivot: Protein × originalRUN, missing → 0
    wide = (df.pivot_table(index="Protein", columns="originalRUN",
                           values="Intensity", aggfunc="mean")
              .fillna(0.0)
              .reset_index())

    mq_data = pd.DataFrame()
    mq_data["Majority protein IDs"] = wide["Protein"].astype(str)
    mq_data["Protein names"] = ""
    mq_data["Gene names"] = wide["Protein"].astype(str)
    mq_data["Reverse"] = "-"
    mq_data["Only identified by site"] = "-"
    mq_data["Potential contaminant"] = "-"
    mq_data["Sequence length"] = 1

    for col in wide.columns:
        if col == "Protein":
            continue
        run_name = str(col)
        if run_name in experimental_design.name2experiment:
            mq_data[f"Intensity {run_name}"] = wide[col].fillna(0.0).astype(float)

    return mq_data


@_parse_stage
def parse_msstats(msstatsFile, experimentalDesign, outputPath):
    """
    Parse an MSstats ProteinLevelData.csv (output of MSstats::dataProcess())
    plus an ED file. The MSstats data is already log2-transformed, normalized,
    and (optionally) imputed; we back-transform to linear intensities and feed
    the existing pipeline. Recommended: run scoring with --imputation 0.
    """
    try:
        ed_df, msstats_df = validate_msstats_inputs(experimentalDesign, msstatsFile)
    except ProxiMateError:
        raise

    if not os.path.exists(outputPath):
        logger.info("Creating output directory: %s", outputPath)
        os.makedirs(outputPath)

    add_file_handler(os.path.join(outputPath, "proximate.log"))

    shutil.copy(experimentalDesign, f"{outputPath}/ED.csv")
    shutil.copy(msstatsFile, f"{outputPath}/ProteinLevelData.csv")

    experimental_design = ExperimentalDesign(experimentalDesign)
    num_expts = experimental_design.num_experiments
    num_ctrls = experimental_design.num_controls

    mq_format_data = convert_msstats_to_maxquant_format(msstatsFile, experimental_design)

    # QC sidecar: preserve MSstats-only columns for downstream inspection
    qc_cols = [c for c in ["Protein", "originalRUN", "GROUP", "SUBJECT",
                           "NumMeasuredFeature", "NumImputedFeature",
                           "MissingPercentage", "more50missing"]
               if c in msstats_df.columns]
    if qc_cols:
        msstats_df[qc_cols].to_csv(os.path.join(outputPath, "msstats_qc.csv"), index=False)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, newline='') as tmp_file:
        tmp_path = tmp_file.name
        mq_format_data.to_csv(tmp_file, sep="\t", index=False)

    shutil.copy(tmp_path, f"{outputPath}/proteinGroups.txt")

    try:
        # MSstats output is always intensity-style; pin quant type accordingly.
        protein_groups = ProteinGroups(experimental_design, tmp_path,
                                       "Intensity", "Intensity")
        protein_groups.to_SAINT(outputPath)
        protein_groups.to_CompPASS(outputPath)
    finally:
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
    matched = set()
    for col in fp_data.columns:
        if col.endswith(fp_suffix):
            # Disambiguate: when in Intensity mode, skip MaxLFQ Intensity columns
            if quant_type == "Intensity" and col.endswith(" MaxLFQ Intensity"):
                continue
            sample_name = col[:-len(fp_suffix)]
            if sample_name in experimental_design.name2experiment:
                mq_data[f"{mq_prefix}{sample_name}"] = fp_data[col].fillna(0)
                matched.add(sample_name)

    # A name that appears in the design but not in the file contributes no
    # quantification at all, which is otherwise indistinguishable from a run of
    # all-zero intensities.
    unmatched = set(experimental_design.name2experiment) - matched
    if unmatched:
        logger.warning("No '%s' column in %s for %d experiment(s) in the design: %s",
                       fp_suffix.strip(), os.path.basename(fp_file),
                       len(unmatched), sorted(unmatched))
    if not matched:
        logger.error("No quantification columns matched the experimental design; "
                     "check that quantType=%r matches the FragPipe export.", quant_type)
    logger.info("Matched %d of %d design experiments to FragPipe columns",
                len(matched), len(experimental_design.name2experiment))

    return mq_data

@_parse_stage
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

    logger.info("Parsing FragPipe file: %s (quantType=%s)", fp_file, quantType)
    logger.info("FragPipe input: %d rows, %d columns", len(fp_df), len(fp_df.columns))

    shutil.copy(experimentalDesign, f"{outputPath}/ED.csv")

    experimental_design = ExperimentalDesign(experimentalDesign)
    num_expts = experimental_design.num_experiments
    num_ctrls = experimental_design.num_controls
    logger.info("Experimental design: %d experiments (%d controls)", num_expts, num_ctrls)

    mq_format_data = convert_fragpipe_to_maxquant_format(fp_file, experimental_design, quantType)
    logger.info("Converted FragPipe to MaxQuant format: %d proteins x %d columns",
                len(mq_format_data), len(mq_format_data.columns))

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

    parser.add_argument("--pioneerMatrix",
                        help="path to Pioneer protein_groups_wide.tsv file",
                        default=None)

    parser.add_argument("--fragpipeFile",
                        help="path to FragPipe combined_protein.tsv file",
                        default=None)

    parser.add_argument("--msstatsFile",
                        help="path to MSstats ProteinLevelData.csv file",
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
    elif args.pioneerMatrix is not None:
        # Pioneer input mode
        if args.experimentalDesign is not None:
            parse_pioneer(args.pioneerMatrix, args.experimentalDesign, args.quantType, args.outputPath)
        else:
            logger.error("No experimental design file provided.")
    elif args.msstatsFile is not None:
        # MSstats input mode
        if args.experimentalDesign is not None:
            parse_msstats(args.msstatsFile, args.experimentalDesign, args.outputPath)
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
            parse_ed_pg(args.proteinGroups, args.experimentalDesign,
                        args.quantType, args.outputPath)
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
        logger.error("No input file provided. Please specify --proteinGroups, --diannMatrix, --pioneerMatrix, --fragpipeFile, --msstatsFile, or --bait/--prey/--interaction.")

if __name__ == "__main__":
    main()