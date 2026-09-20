"""
Validation module for ProxiMate experimental design and proteinGroups files.

Performs comprehensive pre-validation before any parsing begins to provide
clear, early feedback about file format and content issues.
"""

import pandas as pd
import os
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
from experimental_design import GROUP_WILDCARD, _parse_group_cell
from protein_groups import apply_column_aliases, get_quant_col_prefix


# Non-run columns of DIA-NN's report.pg_matrix.tsv; every other column is an MS run.
DIANN_METADATA_COLUMNS = frozenset({
    "Protein.Group", "Protein.Names", "Genes", "First.Protein.Description",
    "N.Sequences", "N.Proteotypic.Sequences",
})

# Non-run columns of Pioneer's protein_groups_wide.tsv. Pioneer's output schema
# policy may omit some of them; every other column is an MS run.
PIONEER_METADATA_COLUMNS = frozenset({
    "species", "gene_names", "protein_names", "protein", "target", "entrap_id",
    "global_pg_score", "global_qval",
})

# Non-sample columns of FragPipe's combined_protein.tsv.
FRAGPIPE_METADATA_COLUMNS = frozenset({
    "Protein", "Protein ID", "Entry Name", "Gene", "Protein Length",
    "Organism", "Protein Existence", "Description",
    "Protein Probability", "Top Peptide Probability",
    "Combined Total Peptides", "Combined Spectral Count",
    "Combined Unique Spectral Count", "Combined Total Spectral Count",
    "Indistinguishable Proteins",
})

# FragPipe sample columns are "<sample><suffix>"; ordered longest-first so that
# " Total Spectral Count" is not read as sample "X Total" with " Spectral Count".
FRAGPIPE_QUANT_SUFFIXES = (
    " Unique Spectral Count", " Total Spectral Count",
    " MaxLFQ Intensity", " Spectral Count", " Intensity",
)

CSV_ENCODINGS = ('utf-8', 'utf-8-sig', 'latin-1', 'cp1252')


logger = get_logger(__name__)


def read_csv_any_encoding(filepath):
    """Read a CSV exported from Excel or R, trying each of CSV_ENCODINGS in turn.

    Returns the first non-empty frame; raises ValueError when no encoding yields one.
    """
    last_err = None
    for encoding in CSV_ENCODINGS:
        try:
            df = pd.read_csv(filepath, encoding=encoding)
            if not df.empty:
                return df
        except UnicodeDecodeError as e:
            last_err = e
    raise ValueError(f"Could not read {filepath} as a non-empty CSV: {last_err}")


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
            return read_csv_any_encoding(filepath)
        except pd.errors.ParserError as e:
            raise EDFileFormatError(f"CSV parsing error: {str(e)}")
        except ValueError:
            # No encoding yielded a data row: the file is empty, not malformed.
            raise EDFileEmptyError(filepath)
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
        # Row numbers are reported as the spreadsheet shows them: 1-based, after the header.
        for col in cls.REQUIRED_COLUMNS:
            null_mask = df[col].isnull() | (df[col].astype(str).str.strip() == '')
            if null_mask.any():
                raise EDMissingValueError(col, [idx + 2 for idx in df.index[null_mask]])

        invalid_types = df[~df['Type'].isin(cls.VALID_TYPES)]
        if len(invalid_types) > 0:
            raise EDInvalidTypeError([idx + 2 for idx in invalid_types.index])

        invalid_rows = []
        for idx, val in df['Replicate'].items():
            try:
                if int(val) <= 0:
                    invalid_rows.append(idx + 2)
            except (ValueError, TypeError):
                invalid_rows.append(idx + 2)
        if invalid_rows:
            raise EDInvalidReplicateError(invalid_rows)

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

        # (row_number, type, spec) per row, where spec is None | "*" | frozenset[int];
        # cells the grammar rejects are collected by reason instead.
        parsed = []
        malformed = {}
        for idx, row in df.iterrows():
            row_number = idx + 2
            raw = None if pd.isna(row["Group"]) else row["Group"]
            try:
                spec = _parse_group_cell(raw)
            except EDInvalidGroupError as e:
                malformed.setdefault(e.reason, []).append(row_number)
                continue
            parsed.append((row_number, row["Type"], spec))

        test_wildcard_rows = [rn for rn, t, g in parsed if t == "T" and g == GROUP_WILDCARD]
        test_multi_rows = [rn for rn, t, g in parsed
                           if t == "T" and isinstance(g, frozenset) and len(g) > 1]

        # Raise in priority order (most basic problem first)
        for reason in ("invalid_group_value", "control_wildcard_with_explicit_groups"):
            if malformed.get(reason):
                raise EDInvalidGroupError(reason, malformed[reason])
        if test_wildcard_rows:
            raise EDInvalidGroupError("test_row_wildcard", test_wildcard_rows)
        if test_multi_rows:
            raise EDInvalidGroupError("test_row_multi_group", test_multi_rows)

        # Column present but unused: not a grouped run.
        if all(g is None for _, _, g in parsed):
            return

        missing_on_test = [rn for rn, t, g in parsed if t == "T" and g is None]
        if missing_on_test:
            raise EDInvalidGroupError("test_row_missing_group_in_grouped_run", missing_on_test)

        # Every test group must be covered by a control (universal or explicit)
        test_groups = set()
        explicit_ctrl_groups = set()
        has_universal = False
        for rn, t, g in parsed:
            if t == "T" and isinstance(g, frozenset):
                test_groups |= g
            elif t == "C" and isinstance(g, frozenset):
                explicit_ctrl_groups |= g
            elif t == "C" and g == GROUP_WILDCARD:
                has_universal = True

        if not has_universal:
            orphans = sorted(test_groups - explicit_ctrl_groups)
            if orphans:
                raise EDInvalidGroupError("bait_group_has_no_control", orphans)

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


class _TableValidator:
    """Existence, readability and required-column checks for a tab-separated data table.

    Subclasses set LABEL (for messages), REQUIRED_COLUMNS, MISSING_COLUMN_ERROR and
    HINT (what a valid file is).
    """

    LABEL = "data"
    REQUIRED_COLUMNS = []
    MISSING_COLUMN_ERROR = PGMissingColumnError
    HINT = "Ensure this is a valid tab-separated table"

    @staticmethod
    def validate_file_exists(filepath):
        if not filepath or not os.path.exists(filepath):
            raise PGFileNotFoundError(filepath)

    @classmethod
    def validate_file_format(cls, filepath):
        """Read the header and a few rows (enough to validate, fast for large files)."""
        try:
            df = pd.read_csv(filepath, sep="\t", low_memory=False, nrows=10)
            if df.empty:
                raise PGFileError(
                    message=f"{cls.LABEL} file is empty",
                    user_message=f"The {cls.LABEL} file is empty",
                    suggestions=["Ensure the file contains data"]
                )
            return df
        except PGFileError:
            raise
        except Exception as e:
            raise PGFileError(
                message=f"Error reading {cls.LABEL} file: {str(e)}",
                user_message=f"Unable to read {cls.LABEL} file",
                suggestions=[cls.HINT, "File should be tab-separated",
                             f"Technical details: {str(e)}"]
            )

    @classmethod
    def _columns(cls, df):
        return df.columns

    @classmethod
    def validate_required_columns(cls, df):
        columns = cls._columns(df)
        missing = [col for col in cls.REQUIRED_COLUMNS if col not in columns]
        if missing:
            raise cls.MISSING_COLUMN_ERROR(missing)

    @classmethod
    def validate_file(cls, filepath):
        cls.validate_file_exists(filepath)
        df = cls.validate_file_format(filepath)
        cls.validate_required_columns(df)
        return df


class PGValidator(_TableValidator):
    """Validator for MaxQuant proteinGroups.txt"""

    LABEL = "proteinGroups"
    REQUIRED_COLUMNS = ["Majority protein IDs", "Gene names", "Reverse",
                        "Only identified by site", "Potential contaminant"]
    MISSING_COLUMN_ERROR = PGMissingColumnError
    HINT = "Ensure this is a valid MaxQuant proteinGroups.txt file"

    @classmethod
    def _columns(cls, df):
        """Required columns may appear under their MaxQuant 2.4+ alias names."""
        return apply_column_aliases(df).columns


class FPValidator(_TableValidator):
    """Validator for FragPipe combined_protein.tsv"""

    LABEL = "FragPipe combined_protein.tsv"
    REQUIRED_COLUMNS = ["Protein", "Protein ID", "Gene", "Protein Length"]
    MISSING_COLUMN_ERROR = FPMissingColumnError
    HINT = "Ensure this is a valid FragPipe combined_protein.tsv file"


class _PioneerValidator(_TableValidator):
    LABEL = "Pioneer protein_groups_wide.tsv"
    REQUIRED_COLUMNS = ["protein"]
    HINT = "Ensure this is the protein_groups_wide.tsv written by Pioneer's SearchDIA"

    @classmethod
    def validate_required_columns(cls, df):
        if "protein" not in df.columns:
            raise PGFileError(
                message="Pioneer table lacks a 'protein' column",
                user_message="The Pioneer file has no 'protein' column",
                suggestions=["Upload protein_groups_wide.tsv, not the precursor or long-format table"]
            )


class _DiannValidator(_TableValidator):
    LABEL = "DIA-NN matrix"
    HINT = "Ensure this is a valid DIA-NN report.pg_matrix.tsv file"


class EDPGCrossValidator:
    """Cross-validation between ED and the data file's experiments"""

    @staticmethod
    def _check(ed_df, data_experiments, source):
        """Every design experiment must be in the data; extra data experiments are ignored."""
        ed_experiments = set(ed_df['Experiment Name'].astype(str).unique())
        data_experiments = set(str(e) for e in data_experiments)

        data_only = sorted(data_experiments - ed_experiments)
        if data_only:
            examples = ", ".join(data_only[:5])
            if len(data_only) > 5:
                examples += f" (and {len(data_only) - 5} more)"
            logger.warning("%d experiment(s) in %s will be ignored (not in ED): %s",
                           len(data_only), source, examples)

        ed_only = sorted(ed_experiments - data_experiments)
        if ed_only:
            raise EDPGMismatchError(ed_only, [])

    @classmethod
    def validate_experiment_match(cls, ed_df, pg_df, quant_prefix):
        """proteinGroups experiments are the columns carrying `quant_prefix`."""
        cls._check(ed_df, [col[len(quant_prefix):] for col in pg_df.columns
                           if col.startswith(quant_prefix)], "proteinGroups")

    @classmethod
    def validate_diann_match(cls, ed_df, diann_df):
        """DIA-NN run columns are the raw file paths/names directly."""
        cls._check(ed_df, set(diann_df.columns) - DIANN_METADATA_COLUMNS, "DIA-NN")

    @classmethod
    def validate_pioneer_match(cls, ed_df, pioneer_df):
        """Pioneer run columns are the MS file names without extension."""
        cls._check(ed_df, set(pioneer_df.columns) - PIONEER_METADATA_COLUMNS, "Pioneer")

    @classmethod
    def validate_msstats_match(cls, ed_df, msstats_df):
        """MSstats runs are the originalRUN values."""
        cls._check(ed_df, msstats_df['originalRUN'].unique(), "MSstats ProteinLevelData")

    @classmethod
    def validate_fragpipe_match(cls, ed_df, fp_df):
        """FragPipe sample columns are "<sample><quant suffix>"."""
        samples = set()
        for col in fp_df.columns:
            if col in FRAGPIPE_METADATA_COLUMNS:
                continue
            for suffix in FRAGPIPE_QUANT_SUFFIXES:
                if col.endswith(suffix):
                    if col[:-len(suffix)]:
                        samples.add(col[:-len(suffix)])
                    break
        cls._check(ed_df, samples, "FragPipe")


def validate_maxquant_inputs(ed_file, pg_file, quant_type):
    """
    Validate MaxQuant inputs (ED + proteinGroups).
    Returns (ed_df, pg_df_headers) if successful.

    Raises:
        EDFileError: If ED file has validation issues
        PGFileError: If proteinGroups file has validation issues
        EDPGMismatchError: If experiment names don't match
        ValueError: If quant_type is not a known quantification
    """
    ed_df = EDValidator.validate_ed_file(ed_file)
    pg_df = PGValidator.validate_file(pg_file)
    EDPGCrossValidator.validate_experiment_match(ed_df, pg_df, get_quant_col_prefix(quant_type))
    return ed_df, pg_df


def validate_diann_inputs(ed_file, diann_file):
    """
    Validate DIA-NN inputs (ED + matrix).
    Returns (ed_df, diann_df_headers) if successful.
    """
    ed_df = EDValidator.validate_ed_file(ed_file)
    diann_df = _DiannValidator.validate_file(diann_file)
    EDPGCrossValidator.validate_diann_match(ed_df, diann_df)
    return ed_df, diann_df


def validate_pioneer_inputs(ed_file, pioneer_file):
    """
    Validate Pioneer inputs (ED + protein_groups_wide.tsv).
    Returns (ed_df, pioneer_df_headers) if successful.
    """
    ed_df = EDValidator.validate_ed_file(ed_file)
    pioneer_df = _PioneerValidator.validate_file(pioneer_file)
    EDPGCrossValidator.validate_pioneer_match(ed_df, pioneer_df)
    return ed_df, pioneer_df


MSSTATS_REQUIRED_COLUMNS = ["Protein", "originalRUN", "GROUP", "SUBJECT", "LABEL", "LogIntensities"]


def validate_msstats_inputs(ed_file, msstats_file):
    """
    Validate MSstats inputs (ED + ProteinLevelData.csv).
    Returns (ed_df, msstats_df) if successful; the whole table is read, since its
    runs are rows rather than columns.
    """
    ed_df = EDValidator.validate_ed_file(ed_file)

    if not msstats_file or not os.path.exists(msstats_file):
        raise PGFileNotFoundError(msstats_file)

    try:
        msstats_df = read_csv_any_encoding(msstats_file)
    except Exception as e:
        raise PGFileError(
            message=f"Error reading MSstats ProteinLevelData: {e}",
            user_message="Unable to read MSstats ProteinLevelData.csv file",
            suggestions=[
                "Ensure this is a valid ProteinLevelData CSV produced by MSstats::dataProcess()",
                "File should be comma-separated",
                f"Technical details: {e}"
            ]
        )

    missing = [c for c in MSSTATS_REQUIRED_COLUMNS if c not in msstats_df.columns]
    if missing:
        raise PGFileError(
            message=f"MSstats ProteinLevelData missing columns: {missing}",
            user_message=f"The MSstats ProteinLevelData.csv file is missing required column(s): {', '.join(missing)}",
            suggestions=[
                "Ensure the file is the ProteinLevelData output from MSstats::dataProcess()",
                f"Required columns: {', '.join(MSSTATS_REQUIRED_COLUMNS)}",
            ]
        )

    EDPGCrossValidator.validate_msstats_match(ed_df, msstats_df)

    return ed_df, msstats_df


def validate_fragpipe_inputs(ed_file, fp_file):
    """
    Validate FragPipe inputs (ED + combined_protein.tsv).
    Returns (ed_df, fp_df_headers) if successful.
    """
    ed_df = EDValidator.validate_ed_file(ed_file)
    fp_df = FPValidator.validate_file(fp_file)
    EDPGCrossValidator.validate_fragpipe_match(ed_df, fp_df)
    return ed_df, fp_df
