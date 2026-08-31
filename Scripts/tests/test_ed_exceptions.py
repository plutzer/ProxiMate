"""Tests for the input-error hierarchy.

These exceptions carry the text the GUI shows a user whose file was rejected, and the
class hierarchy decides which ``except`` clause catches what.  ``parse.py`` re-raises on
``ProxiMateError``, so anything outside that tree escapes as an unhandled traceback
instead of reaching the GUI's error panel.
"""

import pytest

from ed_exceptions import (
    EDDuplicateExperimentError,
    EDFileEmptyError,
    EDFileError,
    EDFileFormatError,
    EDFileNotFoundError,
    EDInvalidGroupError,
    EDInvalidReplicateError,
    EDInvalidTypeError,
    EDMissingColumnError,
    EDMissingValueError,
    EDPGMismatchError,
    FPMissingColumnError,
    PGFileError,
    PGFileNotFoundError,
    PGMissingColumnError,
    ProxiMateError,
    _format_row_examples,
)


ALL_ERRORS = [
    (EDFileNotFoundError, ("/tmp/ED.csv",)),
    (EDFileEmptyError, ("/tmp/ED.csv",)),
    (EDFileFormatError, ("bad delimiter",)),
    (EDMissingColumnError, (["Bait"],)),
    (EDInvalidTypeError, ([2, 3],)),
    (EDInvalidReplicateError, ([4],)),
    (EDMissingValueError, ("Bait", [2])),
    (EDDuplicateExperimentError, (["t1"],)),
    (EDInvalidGroupError, ("invalid_group_value", [2])),
    (PGFileNotFoundError, ("/tmp/pg.txt",)),
    (PGMissingColumnError, (["Reverse"],)),
    (FPMissingColumnError, (["Gene"],)),
    (EDPGMismatchError, (["c_1"], [])),
]


# --- the base contract ---------------------------------------------------------

@pytest.mark.parametrize("cls, args", ALL_ERRORS,
                         ids=[c.__name__ for c, _ in ALL_ERRORS])
def test_every_error_is_a_proximate_error(cls, args):
    """parse.py's re-raise and the GUI's handler both key on this base class."""
    assert isinstance(cls(*args), ProxiMateError)


@pytest.mark.parametrize("cls, args", ALL_ERRORS,
                         ids=[c.__name__ for c, _ in ALL_ERRORS])
def test_str_is_the_technical_message(cls, args):
    """Logs get the technical message; the user-facing text is a separate attribute and
    must not leak into a traceback."""
    error = cls(*args)

    assert str(error) == error.message


@pytest.mark.parametrize("cls, args", ALL_ERRORS,
                         ids=[c.__name__ for c, _ in ALL_ERRORS])
def test_every_error_offers_a_user_message_and_suggestions(cls, args):
    error = cls(*args)

    assert error.user_message
    assert error.suggestions


def test_the_user_message_falls_back_to_the_technical_one():
    error = ProxiMateError("something went wrong")

    assert error.user_message == "something went wrong"


def test_suggestions_default_to_empty():
    assert ProxiMateError("m").suggestions == []


# --- hierarchy ------------------------------------------------------------------

@pytest.mark.parametrize("cls, args", [
    (EDFileNotFoundError, ("/tmp/ED.csv",)),
    (EDInvalidTypeError, ([2],)),
    (EDInvalidGroupError, ("invalid_group_value", [2])),
], ids=lambda v: getattr(v, "__name__", ""))
def test_design_file_errors_share_a_base(cls, args):
    assert isinstance(cls(*args), EDFileError)


def test_a_mismatch_is_not_a_design_file_error():
    """It descends straight from ProxiMateError, so a handler catching EDFileError to
    report design problems misses the one raised when the design and the data disagree."""
    error = EDPGMismatchError(["c_1"], [])

    assert isinstance(error, ProxiMateError)
    assert not isinstance(error, EDFileError)


def test_a_fragpipe_column_error_is_a_protein_groups_error():
    """FPMissingColumnError has no FragPipe-specific base, so a handler written for
    proteinGroups problems also catches FragPipe ones."""
    assert isinstance(FPMissingColumnError(["Gene"]), PGFileError)


def test_the_two_file_error_trees_are_disjoint():
    assert not isinstance(PGFileNotFoundError("/tmp/pg.txt"), EDFileError)
    assert not isinstance(EDFileNotFoundError("/tmp/ED.csv"), PGFileError)


# --- messages -------------------------------------------------------------------

def test_a_missing_file_names_the_path():
    assert "/tmp/ED.csv" in EDFileNotFoundError("/tmp/ED.csv").message


def test_an_empty_file_error_tolerates_no_path():
    assert EDFileEmptyError().message == "File is empty"


def test_missing_columns_are_listed_and_kept():
    error = EDMissingColumnError(["Bait", "Replicate"])

    assert "Bait, Replicate" in error.message
    assert error.missing_columns == ["Bait", "Replicate"]


def test_missing_values_name_the_column():
    error = EDMissingValueError("Bait", [2, 5])

    assert "Bait" in error.message
    assert error.row_indices == [2, 5]


def test_duplicates_are_listed():
    error = EDDuplicateExperimentError(["t1", "t2"])

    assert "t1, t2" in error.message
    assert error.duplicates == ["t1", "t2"]


def test_a_mismatch_lists_the_design_only_names():
    error = EDPGMismatchError(["c_1", "c_2"], [])

    assert error.ed_only == ["c_1", "c_2"]
    assert "c_1" in error.user_message


# --- row-example truncation ------------------------------------------------------

def test_short_row_lists_are_shown_in_full():
    assert _format_row_examples([2, 3, 4]) == "2, 3, 4"


def test_long_row_lists_are_truncated_with_a_count():
    assert _format_row_examples([2, 3, 4, 5, 6, 7, 8]) == "2, 3, 4, 5, 6 (and 2 more)"


def test_the_truncation_limit_is_configurable():
    assert _format_row_examples([2, 3, 4], limit=2) == "2, 3 (and 1 more)"


@pytest.mark.parametrize("cls", [EDInvalidTypeError, EDInvalidReplicateError])
def test_row_errors_truncate_at_five(cls):
    error = cls(list(range(2, 12)))

    assert "(and 5 more)" in " ".join(error.suggestions)


def test_duplicate_errors_truncate_at_three():
    """Documented, not fixed: the row-based errors show five examples and this one shows
    three, so the same file yields differently sized lists depending on what was wrong."""
    error = EDDuplicateExperimentError(["a", "b", "c", "d", "e"])

    assert "(and 2 more)" in " ".join(error.suggestions)


# --- EDInvalidGroupError reasons -------------------------------------------------

GROUP_REASONS = [
    "invalid_group_value",
    "test_row_multi_group",
    "test_row_wildcard",
    "test_row_missing_group_in_grouped_run",
    "control_wildcard_with_explicit_groups",
    "bait_group_has_no_control",
]


@pytest.mark.parametrize("reason", GROUP_REASONS)
def test_each_group_reason_produces_its_own_text(reason):
    error = EDInvalidGroupError(reason, [2, 3])

    assert error.reason == reason
    assert error.offenders == [2, 3]
    assert error.suggestions


def test_group_reasons_do_not_share_a_message():
    messages = {EDInvalidGroupError(r, [2]).message for r in GROUP_REASONS}

    assert len(messages) == len(GROUP_REASONS)


def test_an_unknown_reason_falls_back_without_suggestions():
    error = EDInvalidGroupError("something_new", [2])

    assert "something_new" in error.message
    assert error.suggestions == []


def test_the_orphan_group_reason_carries_group_numbers_not_rows():
    """Every other reason lists row numbers; this one lists the groups that have no
    control, so a caller rendering offenders as rows mislabels them."""
    error = EDInvalidGroupError("bait_group_has_no_control", [1, 3])

    assert "1, 3" in error.message


def test_a_single_offender_is_described_in_the_singular():
    assert "1 row " in EDInvalidGroupError("invalid_group_value", [2]).user_message


def test_several_offenders_are_described_in_the_plural():
    assert "2 rows" in EDInvalidGroupError("invalid_group_value", [2, 3]).user_message
