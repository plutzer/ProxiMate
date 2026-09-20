"""
Custom exceptions for ProxiMate experimental design file parsing.

Each exception provides:
- message: Technical description for logging
- user_message: Non-technical explanation for users
- suggestions: List of actionable fixes
"""


def _format_row_examples(items, limit=5):
    """Join up to `limit` items, with an 'and N more' suffix for the rest."""
    examples = ", ".join(str(r) for r in items[:limit])
    if len(items) > limit:
        examples += f" (and {len(items) - limit} more)"
    return examples


def _plural(n):
    return "s" if n != 1 else ""


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
        n = len(invalid_rows)
        message = f"Invalid Type values in {n} rows"
        user_message = f"Found {n} row{_plural(n)} with invalid Type values"
        suggestions = [
            "Type must be either 'C' (control) or 'T' (test/bait)",
            "Type values are case-sensitive",
            f"Problem rows: {_format_row_examples(invalid_rows)}"
        ]
        super().__init__(message, user_message, suggestions)


class EDInvalidReplicateError(EDFileError):
    """Replicate column contains non-numeric or invalid values"""

    def __init__(self, invalid_rows):
        self.invalid_rows = invalid_rows
        n = len(invalid_rows)
        message = f"Invalid Replicate values in {n} rows"
        user_message = f"Found {n} row{_plural(n)} with invalid Replicate values"
        suggestions = [
            "Replicate must be a positive integer (1, 2, 3, etc.)",
            "Remove any text or special characters from the Replicate column",
            f"Problem rows: {_format_row_examples(invalid_rows)}"
        ]
        super().__init__(message, user_message, suggestions)


class EDNoTestExperimentsError(EDFileError):
    """The design declares no Type == "T" rows, so there is nothing to score."""

    def __init__(self):
        super().__init__(
            "Experimental design has no test (Type 'T') rows",
            "The Experimental Design file contains no test experiments",
            ["At least one row must have Type 'T' (a bait experiment)",
             "Check that the Type column is not 'C' on every row"],
        )


# Text for each EDInvalidGroupError reason.  {n} and {s} are the offender count and its
# plural suffix; {examples} is the formatted offender list.
_GROUP_ERROR_TEMPLATES = {
    "invalid_group_value": (
        "Invalid Group values in {n} row{s}",
        "Found {n} row{s} with an invalid Group value",
        ["Group must be a positive integer (1, 2, 3, ...)",
         "For multi-group controls, use a comma-separated list of integers (e.g. '1,2')",
         "Use '*' on a control row to mark it as a universal control",
         "Problem rows: {examples}"],
    ),
    "test_row_multi_group": (
        "Test rows declare multiple groups in {n} row{s}",
        "Found {n} test row{s} that declare more than one group",
        ["Each test bait row must declare exactly one group number",
         "Only control rows may be shared across multiple groups",
         "Problem rows: {examples}"],
    ),
    "test_row_wildcard": (
        "Test rows use '*' in {n} row{s}",
        "Found {n} test row{s} using '*' (wildcard) for Group",
        ["'*' is reserved for control rows (universal controls)",
         "Each test bait must declare a specific group number",
         "Problem rows: {examples}"],
    ),
    "test_row_missing_group_in_grouped_run": (
        "Test rows are missing Group values in a grouped run ({n} row{s})",
        "Found {n} test row{s} with no Group value in a grouped run",
        ["When any row declares a Group, every test row must also declare one",
         "Fill in the Group column for all test rows, or clear it entirely for a standard run",
         "Problem rows: {examples}"],
    ),
    "control_wildcard_with_explicit_groups": (
        "Control rows mix '*' with explicit groups in {n} row{s}",
        "Found {n} control row{s} that mix '*' with specific group numbers",
        ["Use either '*' alone (universal) or a list of group numbers — not both",
         "Example: '*' OR '1,2', but not '*,1'",
         "Problem rows: {examples}"],
    ),
    # Offenders here are group numbers, not row numbers.
    "bait_group_has_no_control": (
        "Test groups with no matching control: {examples}",
        "These test group{s} have no matching control: {examples}",
        ["Every test group must be covered by at least one control row",
         "Assign a control to each missing group, or add a universal control ('*')"],
    ),
}


class EDInvalidGroupError(EDFileError):
    """Group column contains invalid values or inconsistent group assignments."""

    def __init__(self, reason, offenders):
        self.reason = reason
        self.offenders = offenders
        if reason not in _GROUP_ERROR_TEMPLATES:
            raise ValueError(f"Unknown group error reason: {reason!r}")
        message, user_message, suggestions = _GROUP_ERROR_TEMPLATES[reason]
        n = len(offenders)
        fields = {"n": n, "s": _plural(n), "examples": _format_row_examples(offenders)}
        super().__init__(message.format(**fields), user_message.format(**fields),
                         [s.format(**fields) for s in suggestions])


class EDMissingValueError(EDFileError):
    """Required fields contain empty/null values"""

    def __init__(self, column_name, row_indices):
        self.column_name = column_name
        self.row_indices = row_indices
        n = len(row_indices)
        message = f"Missing values in column '{column_name}' at rows: {row_indices}"
        user_message = f"Column '{column_name}' has empty values in {n} row{_plural(n)}"
        suggestions = [
            f"All rows must have a value for '{column_name}'",
            "Check for blank cells or missing data",
            f"Problem rows: {_format_row_examples(row_indices)}"
        ]
        super().__init__(message, user_message, suggestions)


class EDDuplicateExperimentError(EDFileError):
    """Duplicate Experiment Names found"""

    def __init__(self, duplicates):
        self.duplicates = duplicates
        n = len(duplicates)
        message = f"Duplicate Experiment Names: {', '.join(duplicates)}"
        user_message = f"Found {n} duplicate Experiment Name{_plural(n)}"
        suggestions = [
            "Each Experiment Name must be unique",
            "Use different names for each experiment (e.g., Sample_1, Sample_2)",
            f"Duplicates found: {_format_row_examples(duplicates, limit=3)}"
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
            "Required columns include: Majority protein IDs, Gene names, "
            "Reverse (named Decoy by MaxQuant 2.4 and later), etc."
        ]
        super().__init__(message, user_message, suggestions)


class FPMissingColumnError(PGFileError):
    """FragPipe file missing required columns"""

    def __init__(self, missing_columns):
        self.missing_columns = missing_columns
        missing_str = ", ".join(missing_columns)
        message = f"Missing required columns: {missing_str}"
        user_message = f"Your FragPipe file is missing these required columns: {missing_str}"
        suggestions = [
            "This should be a FragPipe combined_protein.tsv file",
            "Ensure the file hasn't been modified or filtered",
            "Required columns include: Protein, Protein ID, Gene, Protein Length"
        ]
        super().__init__(message=message, user_message=user_message, suggestions=suggestions)


class EDPGMismatchError(ProxiMateError):
    """Mismatch between ED experiments and proteinGroups columns"""

    def __init__(self, ed_only, pg_only):
        self.ed_only = ed_only
        self.pg_only = pg_only
        message = "Mismatch between Experimental Design and proteinGroups"

        details = []
        if ed_only:
            details.append(f"In ED but not in proteinGroups: {_format_row_examples(ed_only, limit=3)}")
        if pg_only:
            details.append(f"In proteinGroups but not in ED: {_format_row_examples(pg_only, limit=3)}")

        user_message = "Experiment names don't match between files:\n" + "\n".join(details)
        suggestions = [
            "All Experiment Names in ED must have matching columns in the data file",
            "Extra experiments in the data file are OK (they will be ignored)",
            "Check for typos, extra spaces, or different capitalization in ED file",
            "For MaxQuant: column names look like 'Intensity [ExperimentName]'",
            "For DIA-NN: column names are raw file names",
            "For Pioneer: column names are MS file names without extension, as in protein_groups_wide.tsv"
        ]
        super().__init__(message, user_message, suggestions)
