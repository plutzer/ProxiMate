"""
Validation module for ProxiMate experimental design and proteinGroups files.

Performs comprehensive pre-validation before any parsing begins to provide
clear, early feedback about file format and content issues.
"""

import pandas as pd
import os
import re
from log_config import get_logger
from ed_exceptions import (
    EDFileNotFoundError,
    EDFileEmptyError,
    EDFileFormatError,
    EDMissingColumnError,
    EDInvalidTypeError,
    EDInvalidReplicateError,
    EDMissingValueError,
    EDDuplicateExperimentError,
    EDInvalidGroupError,
    PGFileNotFoundError,
    PGFileError,
    PGMissingColumnError,
    FPMissingColumnError,
    EDPGMismatchError
)


GROUP_WILDCARD = "*"
_GROUP_INT_PATTERN = re.compile(r"^[1-9][0-9]*$")


logger = get_logger(__name__)


class EDValidator:
    """
    Comprehensive validator for Experimental Design files.
    Performs all checks BEFORE any parsing begins to provide clear, early feedback.
    """

    REQUIRED_COLUMNS = ["Experiment Name", "Type", "Bait", "Replicate"]
    VALID_TYPES = ["C", "T"]

    @staticmethod
    def validate_file_exists(filepath):
        """Check if file exists and is readable"""
        if not filepath or not os.path.exists(filepath):
            raise EDFileNotFoundError(filepath)

        if os.path.getsize(filepath) == 0:
            raise EDFileEmptyError(filepath)

    @staticmethod
    def validate_file_format(filepath):
        """
        Try to read the file and validate basic format.
        Returns the DataFrame if successful.
        """
        try:
            # Try reading with different encodings
            df = None
            for encoding in ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252']:
                try:
                    df = pd.read_csv(filepath, encoding=encoding)
                    if not df.empty:
                        break
                except UnicodeDecodeError:
                    continue

            if df is None or df.empty:
                raise EDFileEmptyError(filepath)

            return df

        except pd.errors.ParserError as e:
            raise EDFileFormatError(f"CSV parsing error: {str(e)}")
        except Exception as e:
            raise EDFileFormatError(f"Unexpected error reading file: {str(e)}")

    @classmethod
    def validate_required_columns(cls, df):
        """Check all required columns are present"""
        missing = [col for col in cls.REQUIRED_COLUMNS if col not in df.columns]
        if missing:
            raise EDMissingColumnError(missing)

    @classmethod
    def validate_column_values(cls, df):
        """Validate data within each column"""
        # Check for empty values in required columns
        for col in cls.REQUIRED_COLUMNS:
            null_mask = df[col].isnull() | (df[col].astype(str).str.strip() == '')
            if null_mask.any():
                # Get row indices (1-based for user display)
                row_indices = [idx + 2 for idx in df.index[null_mask].tolist()]  # +2 for header row + 0-based indexing
                raise EDMissingValueError(col, row_indices)

        # Validate Type column
        invalid_types = df[~df['Type'].isin(cls.VALID_TYPES)]
        if len(invalid_types) > 0:
            # Get row indices (1-based for user display)
            row_indices = [idx + 2 for idx in invalid_types.index.tolist()]
            raise EDInvalidTypeError(row_indices)

        # Validate Replicate column
        invalid_rows = []
        for idx, val in df['Replicate'].items():
            try:
                int_val = int(val)
                if int_val <= 0:
                    invalid_rows.append(idx + 2)  # 1-based + header row
            except (ValueError, TypeError):
                invalid_rows.append(idx + 2)  # 1-based + header row

        if invalid_rows:
            raise EDInvalidReplicateError(invalid_rows)

        # Check for duplicate Experiment Names
        duplicates = df[df.duplicated(subset=['Experiment Name'], keep=False)]['Experiment Name'].unique().tolist()
        if duplicates:
            raise EDDuplicateExperimentError(duplicates)

    @classmethod
    def validate_group_column(cls, df):
        """Validate the optional Group column for paired-control runs.

        Absent column or all-empty column = legacy run, no validation.
        Grouped run (any non-empty cell): enforce the full rule set.
        """
        if "Group" not in df.columns:
            return

        # Parse each row's cell into (row_number, type, parsed_spec) or collect per-category errors.
        invalid_value_rows = []
        test_multi_rows = []
        test_wildcard_rows = []
        wildcard_mix_rows = []

        parsed = []  # list of (row_number, type, spec) where spec: None | "*" | frozenset[int]

        for idx, row in df.iterrows():
            row_number = idx + 2  # 1-based header + 0-based index
            row_type = row["Type"]
            raw = row.get("Group")

            if pd.isna(raw):
                parsed.append((row_number, row_type, None))
                continue

            s = str(raw).strip()
            if s == "":
                parsed.append((row_number, row_type, None))
                continue

            if s == GROUP_WILDCARD:
                if row_type == "T":
                    test_wildcard_rows.append(row_number)
                parsed.append((row_number, row_type, GROUP_WILDCARD))
                continue

            tokens = [t.strip() for t in s.split(",")]
            has_wildcard = GROUP_WILDCARD in tokens
            explicit = [t for t in tokens if t != GROUP_WILDCARD and t != ""]

            if has_wildcard and explicit:
                wildcard_mix_rows.append(row_number)
                parsed.append((row_number, row_type, None))
                continue

            malformed = False
            values = set()
            for t in explicit:
                if not _GROUP_INT_PATTERN.match(t):
                    invalid_value_rows.append(row_number)
                    malformed = True
                    break
                values.add(int(t))

            if malformed:
                parsed.append((row_number, row_type, None))
                continue

            if not values:
                invalid_value_rows.append(row_number)
                parsed.append((row_number, row_type, None))
                continue

            if row_type == "T" and len(values) > 1:
                test_multi_rows.append(row_number)

            parsed.append((row_number, row_type, frozenset(values)))

        # Raise errors in priority order (most actionable first)
        if invalid_value_rows:
            raise EDInvalidGroupError("invalid_group_value", invalid_value_rows)
        if wildcard_mix_rows:
            raise EDInvalidGroupError("control_wildcard_with_explicit_groups", wildcard_mix_rows)
        if test_wildcard_rows:
            raise EDInvalidGroupError("test_row_wildcard", test_wildcard_rows)
        if test_multi_rows:
            raise EDInvalidGroupError("test_row_multi_group", test_multi_rows)

        # If no row has a non-empty group, this is a legacy run (column present but unused).
        non_empty = [p for p in parsed if p[2] is not None]
        if not non_empty:
            return

        # Grouped run: every test row must declare a group
        missing_on_test = [
            rn for (rn, t, g) in parsed
            if t == "T" and g is None
        ]
        if missing_on_test:
            raise EDInvalidGroupError("test_row_missing_group_in_grouped_run", missing_on_test)

        # Every test group must be covered by at least one control (universal or explicit)
        test_groups = set()
        for rn, t, g in parsed:
            if t == "T" and isinstance(g, frozenset):
                test_groups |= g

        has_universal = any(t == "C" and g == GROUP_WILDCARD for (rn, t, g) in parsed)
        explicit_ctrl_groups = set()
        for rn, t, g in parsed:
            if t == "C" and isinstance(g, frozenset):
                explicit_ctrl_groups |= g

        if not has_universal:
            orphans = sorted(test_groups - explicit_ctrl_groups)
            if orphans:
                raise EDInvalidGroupError("bait_group_has_no_control", orphans)

        # Non-fatal: warn about controls whose declared groups are never referenced by any test
        unreferenced = sorted(explicit_ctrl_groups - test_groups)
        if unreferenced:
            logger.warning(
                "Control row(s) reference groups with no test baits: %s",
                ", ".join(str(g) for g in unreferenced),
            )

    @classmethod
    def validate_ed_file(cls, filepath):
        """
        Complete validation of ED file.
        Returns DataFrame if all validations pass.
        Raises specific exception if any validation fails.
        """
        cls.validate_file_exists(filepath)
        df = cls.validate_file_format(filepath)
        cls.validate_required_columns(df)
        cls.validate_column_values(df)
        cls.validate_group_column(df)

        return df


class PGValidator:
    """Validator for proteinGroups files"""

    REQUIRED_COLUMNS = ["Majority protein IDs", "Gene names", "Reverse",
                       "Only identified by site", "Potential contaminant"]

    @staticmethod
    def validate_file_exists(filepath):
        """Check if file exists and is readable"""
        if not filepath or not os.path.exists(filepath):
            raise PGFileNotFoundError(filepath)

    @staticmethod
    def validate_file_format(filepath):
        """
        Try to read the proteinGroups file.
        Returns the DataFrame if successful (just headers + few rows for speed).
        """
        try:
            # MaxQuant files are tab-separated
            # Read just header + few rows for validation (much faster for large files)
            df = pd.read_csv(filepath, sep="\t", low_memory=False, nrows=10)
            if df.empty:
                raise PGFileError(
                    message="proteinGroups file is empty",
                    user_message="The proteinGroups file is empty",
                    suggestions=["Ensure the file contains data"]
                )
            return df
        except Exception as e:
            raise PGFileError(
                message=f"Error reading proteinGroups: {str(e)}",
                user_message=f"Unable to read proteinGroups file",
                suggestions=[
                    "Ensure this is a valid MaxQuant proteinGroups.txt file",
                    "File should be tab-separated",
                    f"Technical details: {str(e)}"
                ]
            )

    @classmethod
    def validate_required_columns(cls, df):
        """Check required columns are present"""
        missing = [col for col in cls.REQUIRED_COLUMNS if col not in df.columns]
        if missing:
            raise PGMissingColumnError(missing)

    @classmethod
    def validate_pg_file(cls, filepath):
        """
        Complete validation of proteinGroups file.
        Returns DataFrame headers if all validations pass.
        """
        cls.validate_file_exists(filepath)
        df = cls.validate_file_format(filepath)
        cls.validate_required_columns(df)

        return df


class FPValidator:
    """Validator for FragPipe combined_protein.tsv files"""

    REQUIRED_COLUMNS = ["Protein", "Protein ID", "Gene", "Protein Length"]

    @staticmethod
    def validate_file_exists(filepath):
        """Check if file exists and is readable"""
        if not filepath or not os.path.exists(filepath):
            raise PGFileNotFoundError(filepath)

    @staticmethod
    def validate_file_format(filepath):
        """
        Try to read the FragPipe file.
        Returns the DataFrame if successful (just headers + few rows for speed).
        """
        try:
            df = pd.read_csv(filepath, sep="\t", low_memory=False, nrows=10)
            if df.empty:
                raise PGFileError(
                    message="FragPipe file is empty",
                    user_message="The FragPipe combined_protein.tsv file is empty",
                    suggestions=["Ensure the file contains data"]
                )
            return df
        except PGFileError:
            raise
        except Exception as e:
            raise PGFileError(
                message=f"Error reading FragPipe file: {str(e)}",
                user_message="Unable to read FragPipe file",
                suggestions=[
                    "Ensure this is a valid FragPipe combined_protein.tsv file",
                    "File should be tab-separated",
                    f"Technical details: {str(e)}"
                ]
            )

    @classmethod
    def validate_required_columns(cls, df):
        """Check required columns are present"""
        missing = [col for col in cls.REQUIRED_COLUMNS if col not in df.columns]
        if missing:
            raise FPMissingColumnError(missing)

    @classmethod
    def validate_fp_file(cls, filepath):
        """
        Complete validation of FragPipe file.
        Returns DataFrame headers if all validations pass.
        """
        cls.validate_file_exists(filepath)
        df = cls.validate_file_format(filepath)
        cls.validate_required_columns(df)
        return df


class EDPGCrossValidator:
    """Cross-validation between ED and proteinGroups files"""

    @staticmethod
    def validate_experiment_match(ed_df, pg_df, quant_prefix):
        """
        Validate that experiments in ED match columns in proteinGroups.

        Args:
            ed_df: Experimental Design DataFrame
            pg_df: proteinGroups DataFrame (just headers needed)
            quant_prefix: Quantification prefix (e.g., "Intensity ", "LFQ intensity ")
        """
        # Get experiment names from ED
        ed_experiments = set(ed_df['Experiment Name'].unique())

        # Get experiment names from proteinGroups columns
        pg_experiments = set()
        for col in pg_df.columns:
            if col.startswith(quant_prefix):
                exp_name = col.replace(quant_prefix, '', 1)  # Remove prefix
                pg_experiments.add(exp_name)

        # Find mismatches
        ed_only = sorted(list(ed_experiments - pg_experiments))
        pg_only = sorted(list(pg_experiments - ed_experiments))

        # Warn about extra experiments in data file (not an error)
        if pg_only:
            pg_examples = ", ".join(pg_only[:5])
            if len(pg_only) > 5:
                pg_examples += f" (and {len(pg_only) - 5} more)"
            logger.warning("%d experiment(s) in proteinGroups will be ignored (not in ED): %s", len(pg_only), pg_examples)

        # Error only if ED experiments are missing from data file
        if ed_only:
            raise EDPGMismatchError(ed_only, [])

    @staticmethod
    def validate_diann_match(ed_df, diann_df):
        """
        Validate that experiments in ED match columns in DIA-NN matrix.

        DIA-NN columns are raw file paths/names directly.
        """
        # Get experiment names from ED
        ed_experiments = set(ed_df['Experiment Name'].unique())

        # Get all column names from DIA-NN except metadata columns
        metadata_cols = ["Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
                        "N.Sequences", "N.Proteotypic.Sequences"]
        diann_experiments = set(col for col in diann_df.columns if col not in metadata_cols)

        # Find mismatches
        ed_only = sorted(list(ed_experiments - diann_experiments))
        diann_only = sorted(list(diann_experiments - ed_experiments))

        # Warn about extra experiments in data file (not an error)
        if diann_only:
            diann_examples = ", ".join(diann_only[:5])
            if len(diann_only) > 5:
                diann_examples += f" (and {len(diann_only) - 5} more)"
            logger.warning("%d experiment(s) in DIA-NN will be ignored (not in ED): %s", len(diann_only), diann_examples)

        # Error only if ED experiments are missing from data file
        if ed_only:
            raise EDPGMismatchError(ed_only, [])

    @staticmethod
    def validate_msstats_match(ed_df, msstats_df):
        """
        Validate that experiments in ED match originalRUN values in MSstats ProteinLevelData.
        """
        ed_experiments = set(ed_df['Experiment Name'].astype(str).unique())
        msstats_experiments = set(msstats_df['originalRUN'].astype(str).unique())

        ed_only = sorted(list(ed_experiments - msstats_experiments))
        msstats_only = sorted(list(msstats_experiments - ed_experiments))

        if msstats_only:
            examples = ", ".join(msstats_only[:5])
            if len(msstats_only) > 5:
                examples += f" (and {len(msstats_only) - 5} more)"
            logger.warning("%d run(s) in MSstats ProteinLevelData will be ignored (not in ED): %s",
                           len(msstats_only), examples)

        if ed_only:
            raise EDPGMismatchError(ed_only, [])

    @staticmethod
    def validate_fragpipe_match(ed_df, fp_df):
        """
        Validate that experiments in ED match columns in FragPipe combined_protein.tsv.

        FragPipe columns follow the pattern: {SampleName} {QuantSuffix}
        Suffixes are stripped to extract sample names.
        """
        METADATA_COLS = {
            "Protein", "Protein ID", "Entry Name", "Gene", "Protein Length",
            "Organism", "Protein Existence", "Description",
            "Protein Probability", "Top Peptide Probability",
            "Combined Total Peptides", "Combined Spectral Count",
            "Combined Unique Spectral Count", "Combined Total Spectral Count",
            "Indistinguishable Proteins"
        }
        # Ordered longest-first for greedy matching
        QUANT_SUFFIXES = [
            " Unique Spectral Count", " Total Spectral Count",
            " MaxLFQ Intensity", " Spectral Count", " Intensity"
        ]

        fp_experiments = set()
        for col in fp_df.columns:
            if col in METADATA_COLS:
                continue
            for suffix in QUANT_SUFFIXES:
                if col.endswith(suffix):
                    sample_name = col[:-len(suffix)]
                    if sample_name:
                        fp_experiments.add(sample_name)
                    break

        ed_experiments = set(ed_df['Experiment Name'].unique())
        ed_only = sorted(list(ed_experiments - fp_experiments))
        fp_only = sorted(list(fp_experiments - ed_experiments))

        if fp_only:
            fp_examples = ", ".join(fp_only[:5])
            if len(fp_only) > 5:
                fp_examples += f" (and {len(fp_only) - 5} more)"
            logger.warning("%d experiment(s) in FragPipe will be ignored (not in ED): %s", len(fp_only), fp_examples)

        if ed_only:
            raise EDPGMismatchError(ed_only, [])


def validate_maxquant_inputs(ed_file, pg_file, quant_type):
    """
    Validate MaxQuant inputs (ED + proteinGroups).
    Returns (ed_df, pg_df_headers) if successful.
    Raises specific exception if validation fails.

    Args:
        ed_file: Path to experimental design CSV file
        pg_file: Path to proteinGroups.txt file
        quant_type: Quantification type ("Intensity", "LFQ", or "Spectral Counts")

    Returns:
        tuple: (ed_df, pg_df) - Validated DataFrames

    Raises:
        EDFileError: If ED file has validation issues
        PGFileError: If proteinGroups file has validation issues
        EDPGMismatchError: If experiment names don't match
    """
    # Validate ED file
    ed_df = EDValidator.validate_ed_file(ed_file)

    # Validate proteinGroups file
    pg_df = PGValidator.validate_pg_file(pg_file)

    # Determine quantification prefix
    quant_prefix_map = {
        "Intensity": "Intensity ",
        "LFQ": "LFQ intensity ",
        "Spectral Counts": "MS/MS count "
    }
    quant_prefix = quant_prefix_map.get(quant_type, "Intensity ")

    # Cross-validate
    EDPGCrossValidator.validate_experiment_match(ed_df, pg_df, quant_prefix)

    return ed_df, pg_df


def validate_diann_inputs(ed_file, diann_file):
    """
    Validate DIA-NN inputs (ED + matrix).
    Returns (ed_df, diann_df_headers) if successful.

    Args:
        ed_file: Path to experimental design CSV file
        diann_file: Path to DIA-NN report.pg_matrix.tsv file

    Returns:
        tuple: (ed_df, diann_df) - Validated DataFrames

    Raises:
        EDFileError: If ED file has validation issues
        PGFileError: If DIA-NN file has validation issues
        EDPGMismatchError: If experiment names don't match
    """
    # Validate ED file
    ed_df = EDValidator.validate_ed_file(ed_file)

    # Validate DIA-NN file
    try:
        diann_df = pd.read_csv(diann_file, sep="\t", nrows=10)  # Just headers + few rows
        if diann_df.empty:
            raise PGFileError(
                message="DIA-NN matrix is empty",
                user_message="The DIA-NN matrix file is empty",
                suggestions=["Ensure the file contains data"]
            )
    except Exception as e:
        raise PGFileError(
            message=f"Error reading DIA-NN matrix: {str(e)}",
            user_message=f"Unable to read DIA-NN matrix file",
            suggestions=[
                "Ensure this is a valid DIA-NN report.pg_matrix.tsv file",
                "File should be tab-separated",
                f"Technical details: {str(e)}"
            ]
        )

    # Cross-validate
    EDPGCrossValidator.validate_diann_match(ed_df, diann_df)

    return ed_df, diann_df


def validate_msstats_inputs(ed_file, msstats_file):
    """
    Validate MSstats inputs (ED + ProteinLevelData.csv).
    Returns (ed_df, msstats_df) if successful.

    Args:
        ed_file: Path to experimental design CSV file
        msstats_file: Path to MSstats ProteinLevelData.csv file

    Returns:
        tuple: (ed_df, msstats_df) - Validated DataFrames

    Raises:
        EDFileError: If ED file has validation issues
        PGFileError: If MSstats file has validation issues
        EDPGMismatchError: If experiment names don't match
    """
    ed_df = EDValidator.validate_ed_file(ed_file)

    if not msstats_file or not os.path.exists(msstats_file):
        raise PGFileNotFoundError(msstats_file)

    if os.path.getsize(msstats_file) == 0:
        raise PGFileError(
            message="MSstats ProteinLevelData file is empty",
            user_message="The MSstats ProteinLevelData.csv file is empty",
            suggestions=["Ensure the file contains data from MSstats::dataProcess()"]
        )

    msstats_df = None
    last_err = None
    for encoding in ['utf-8', 'utf-8-sig', 'latin-1', 'cp1252']:
        try:
            msstats_df = pd.read_csv(msstats_file, encoding=encoding)
            if not msstats_df.empty:
                break
        except UnicodeDecodeError as e:
            last_err = e
            continue
        except Exception as e:
            last_err = e
            break

    if msstats_df is None or msstats_df.empty:
        raise PGFileError(
            message=f"Error reading MSstats ProteinLevelData: {last_err}",
            user_message="Unable to read MSstats ProteinLevelData.csv file",
            suggestions=[
                "Ensure this is a valid ProteinLevelData CSV produced by MSstats::dataProcess()",
                "File should be comma-separated",
                f"Technical details: {last_err}"
            ]
        )

    required = ["Protein", "originalRUN", "GROUP", "SUBJECT", "LABEL", "LogIntensities"]
    missing = [c for c in required if c not in msstats_df.columns]
    if missing:
        raise PGFileError(
            message=f"MSstats ProteinLevelData missing columns: {missing}",
            user_message=f"The MSstats ProteinLevelData.csv file is missing required column(s): {', '.join(missing)}",
            suggestions=[
                "Ensure the file is the ProteinLevelData output from MSstats::dataProcess()",
                f"Required columns: {', '.join(required)}",
            ]
        )

    EDPGCrossValidator.validate_msstats_match(ed_df, msstats_df)

    return ed_df, msstats_df


def validate_fragpipe_inputs(ed_file, fp_file):
    """
    Validate FragPipe inputs (ED + combined_protein.tsv).
    Returns (ed_df, fp_df_headers) if successful.

    Args:
        ed_file: Path to experimental design CSV file
        fp_file: Path to FragPipe combined_protein.tsv file

    Returns:
        tuple: (ed_df, fp_df) - Validated DataFrames

    Raises:
        EDFileError: If ED file has validation issues
        PGFileError: If FragPipe file has validation issues
        EDPGMismatchError: If experiment names don't match
    """
    ed_df = EDValidator.validate_ed_file(ed_file)
    fp_df = FPValidator.validate_fp_file(fp_file)
    EDPGCrossValidator.validate_fragpipe_match(ed_df, fp_df)
    return ed_df, fp_df
