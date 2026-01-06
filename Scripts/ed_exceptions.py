"""
Custom exceptions for ProxiMate experimental design file parsing.

Each exception provides:
- message: Technical description for logging
- user_message: Non-technical explanation for users
- suggestions: List of actionable fixes
"""


class ProxiMateError(Exception):
    """Base exception for all ProxiMate errors"""

    def __init__(self, message, user_message=None, suggestions=None):
        self.message = message
        self.user_message = user_message or message
        self.suggestions = suggestions or []
        super().__init__(self.message)


class EDFileError(ProxiMateError):
    """Base class for experimental design file errors"""
    pass


class EDFileNotFoundError(EDFileError):
    """ED file does not exist or is not accessible"""

    def __init__(self, filepath):
        message = f"File not found: {filepath}"
        user_message = "The Experimental Design file could not be found"
        suggestions = [
            "Verify the file path is correct",
            "Check file permissions",
            "Ensure the file was uploaded correctly"
        ]
        super().__init__(message, user_message, suggestions)


class EDFileEmptyError(EDFileError):
    """ED file exists but contains no data"""

    def __init__(self, filepath=None):
        message = f"File is empty: {filepath}" if filepath else "File is empty"
        user_message = "The Experimental Design file is empty or has no data rows"
        suggestions = [
            "Ensure the file contains data rows (not just headers)",
            "Check the file wasn't corrupted during upload",
            "Try opening the file to verify it has content"
        ]
        super().__init__(message, user_message, suggestions)


class EDFileFormatError(EDFileError):
    """ED file has incorrect format (encoding, delimiter, etc.)"""

    def __init__(self, details=""):
        message = f"File format error: {details}"
        user_message = "The Experimental Design file format is invalid or corrupted"
        suggestions = [
            "Ensure the file is comma-separated (CSV format)",
            "Check for malformed rows or unmatched quotes",
            "Try opening in Excel and re-saving as CSV (UTF-8)",
            "Ensure the file has a header row with column names"
        ]
        super().__init__(message, user_message, suggestions)


class EDMissingColumnError(EDFileError):
    """ED file is missing required columns"""

    def __init__(self, missing_columns):
        self.missing_columns = missing_columns
        missing_str = ", ".join(missing_columns)
        message = f"Missing required columns: {missing_str}"
        user_message = f"Your Experimental Design file is missing these required columns: {missing_str}"
        suggestions = [
            "Required columns: Experiment Name, Type, Bait, Replicate",
            "Column names are case-sensitive and must match exactly",
            "Check for typos or extra spaces in column headers"
        ]
        super().__init__(message, user_message, suggestions)


class EDInvalidTypeError(EDFileError):
    """Type column contains invalid values (not C or T)"""

    def __init__(self, invalid_rows):
        self.invalid_rows = invalid_rows
        num_invalid = len(invalid_rows)
        message = f"Invalid Type values in {num_invalid} rows"
        user_message = f"Found {num_invalid} row{'s' if num_invalid != 1 else ''} with invalid Type values"

        # Show first 5 problem rows
        row_examples = ", ".join(str(r) for r in invalid_rows[:5])
        if len(invalid_rows) > 5:
            row_examples += f" (and {len(invalid_rows) - 5} more)"

        suggestions = [
            "Type must be either 'C' (control) or 'T' (test/bait)",
            "Type values are case-sensitive",
            f"Problem rows: {row_examples}"
        ]
        super().__init__(message, user_message, suggestions)


class EDInvalidReplicateError(EDFileError):
    """Replicate column contains non-numeric or invalid values"""

    def __init__(self, invalid_rows):
        self.invalid_rows = invalid_rows
        num_invalid = len(invalid_rows)
        message = f"Invalid Replicate values in {num_invalid} rows"
        user_message = f"Found {num_invalid} row{'s' if num_invalid != 1 else ''} with invalid Replicate values"

        # Show first 5 problem rows
        row_examples = ", ".join(str(r) for r in invalid_rows[:5])
        if len(invalid_rows) > 5:
            row_examples += f" (and {len(invalid_rows) - 5} more)"

        suggestions = [
            "Replicate must be a positive integer (1, 2, 3, etc.)",
            "Remove any text or special characters from the Replicate column",
            f"Problem rows: {row_examples}"
        ]
        super().__init__(message, user_message, suggestions)


class EDMissingValueError(EDFileError):
    """Required fields contain empty/null values"""

    def __init__(self, column_name, row_indices):
        self.column_name = column_name
        self.row_indices = row_indices
        num_missing = len(row_indices)
        message = f"Missing values in column '{column_name}' at rows: {row_indices}"
        user_message = f"Column '{column_name}' has empty values in {num_missing} row{'s' if num_missing != 1 else ''}"

        # Show first 5 problem rows
        row_examples = ", ".join(str(r) for r in row_indices[:5])
        if len(row_indices) > 5:
            row_examples += f" (and {len(row_indices) - 5} more)"

        suggestions = [
            f"All rows must have a value for '{column_name}'",
            "Check for blank cells or missing data",
            f"Problem rows: {row_examples}"
        ]
        super().__init__(message, user_message, suggestions)


class EDDuplicateExperimentError(EDFileError):
    """Duplicate Experiment Names found"""

    def __init__(self, duplicates):
        self.duplicates = duplicates
        num_dupes = len(duplicates)
        message = f"Duplicate Experiment Names: {', '.join(duplicates)}"
        user_message = f"Found {num_dupes} duplicate Experiment Name{'s' if num_dupes != 1 else ''}"

        # Show first 3 duplicates
        dupe_examples = ", ".join(duplicates[:3])
        if len(duplicates) > 3:
            dupe_examples += f" (and {len(duplicates) - 3} more)"

        suggestions = [
            "Each Experiment Name must be unique",
            "Use different names for each experiment (e.g., Sample_1, Sample_2)",
            f"Duplicates found: {dupe_examples}"
        ]
        super().__init__(message, user_message, suggestions)


class PGFileError(ProxiMateError):
    """Base class for proteinGroups file errors"""
    pass


class PGFileNotFoundError(PGFileError):
    """proteinGroups file does not exist"""

    def __init__(self, filepath):
        message = f"File not found: {filepath}"
        user_message = "The proteinGroups file could not be found"
        suggestions = [
            "Verify the file path is correct",
            "Check file permissions",
            "Ensure you uploaded the proteinGroups.txt file"
        ]
        super().__init__(message, user_message, suggestions)


class PGMissingColumnError(PGFileError):
    """proteinGroups file missing required columns"""

    def __init__(self, missing_columns):
        self.missing_columns = missing_columns
        missing_str = ", ".join(missing_columns)
        message = f"Missing required columns: {missing_str}"
        user_message = f"Your proteinGroups file is missing these required columns: {missing_str}"
        suggestions = [
            "This should be a standard MaxQuant proteinGroups.txt file",
            "Ensure the file hasn't been modified or filtered",
            "Required columns include: Majority protein IDs, Gene names, Reverse, etc."
        ]
        super().__init__(message, user_message, suggestions)


class EDPGMismatchError(ProxiMateError):
    """Mismatch between ED experiments and proteinGroups columns"""

    def __init__(self, ed_only, pg_only):
        self.ed_only = ed_only
        self.pg_only = pg_only
        message = "Mismatch between Experimental Design and proteinGroups"

        details = []
        if ed_only:
            ed_examples = ", ".join(ed_only[:3])
            if len(ed_only) > 3:
                ed_examples += f" (and {len(ed_only) - 3} more)"
            details.append(f"In ED but not in proteinGroups: {ed_examples}")
        if pg_only:
            pg_examples = ", ".join(pg_only[:3])
            if len(pg_only) > 3:
                pg_examples += f" (and {len(pg_only) - 3} more)"
            details.append(f"In proteinGroups but not in ED: {pg_examples}")

        user_message = "Experiment names don't match between files:\n" + "\n".join(details)
        suggestions = [
            "Experiment Names in ED file must exactly match column names in proteinGroups",
            "Check for typos, extra spaces, or different capitalization",
            "For MaxQuant: column names typically look like 'Intensity [ExperimentName]'",
            "Ensure both files are from the same analysis"
        ]
        super().__init__(message, user_message, suggestions)
