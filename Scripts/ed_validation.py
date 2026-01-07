"""
Validation module for ProxiMate experimental design and proteinGroups files.

Performs comprehensive pre-validation before any parsing begins to provide
clear, early feedback about file format and content issues.
"""

import pandas as pd
import os
from ed_exceptions import (
    EDFileNotFoundError,
    EDFileEmptyError,
    EDFileFormatError,
    EDMissingColumnError,
    EDInvalidTypeError,
    EDInvalidReplicateError,
    EDMissingValueError,
    EDDuplicateExperimentError,
    PGFileNotFoundError,
    PGFileError,
    PGMissingColumnError,
    EDPGMismatchError
)


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

        if ed_only or pg_only:
            raise EDPGMismatchError(ed_only, pg_only)

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

        if ed_only or diann_only:
            raise EDPGMismatchError(ed_only, diann_only)


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
