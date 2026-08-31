"""Tests for reading the experimental design file.

``ExperimentalDesign`` is what scoring reads, and it is a second, independent
implementation of rules ``EDValidator`` also enforces -- it uses the csv module where the
validator uses pandas, and the two disagree in places.  The GUI validates first and would
catch most bad input, but the CLI reaches this class directly.

``get_experiments_for_group`` decides which experiments each per-group SAINT run sees, so
a prey's evidence depends on it.
"""

import csv

import pytest

from ed_exceptions import (
    EDFileEmptyError,
    EDFileNotFoundError,
    EDInvalidGroupError,
    EDInvalidReplicateError,
    EDInvalidTypeError,
    EDMissingColumnError,
)
from experimental_design import (
    GROUP_WILDCARD,
    Experiment,
    ExperimentalDesign,
    _parse_group_cell,
)


HEADER = ["Experiment Name", "Type", "Bait", "Replicate", "Bait ID"]


def _write_ed(tmp_path, rows, header=HEADER, name="ED.csv"):
    """Write an ED file from raw cell lists, so a malformed row can be expressed.

    Written through the csv module rather than by joining on commas: a multi-group
    control cell such as "1,2" has to be quoted or it arrives as two columns.
    """
    path = tmp_path / name
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow(["" if c is None else c for c in row])
    return str(path)


# --- _parse_group_cell ---------------------------------------------------------

@pytest.mark.parametrize("raw, expected", [
    (None, None),
    ("", None),
    ("   ", None),
    ("*", GROUP_WILDCARD),
    (" * ", GROUP_WILDCARD),
    ("1", frozenset({1})),
    ("1,2", frozenset({1, 2})),
    (" 1 , 2 ", frozenset({1, 2})),
    ("2,2", frozenset({2})),
    ("10", frozenset({10})),
])
def test_group_cells_parse_to_a_spec(raw, expected):
    assert _parse_group_cell(raw) == expected


@pytest.mark.parametrize("raw", ["0", "01", "-1", "1.0", "a", "1;2", ",", ",,", " , "])
def test_malformed_group_cells_are_rejected(raw):
    with pytest.raises(EDInvalidGroupError) as excinfo:
        _parse_group_cell(raw)

    assert excinfo.value.reason == "invalid_group_value"


def test_a_wildcard_mixed_with_explicit_groups_is_rejected():
    with pytest.raises(EDInvalidGroupError) as excinfo:
        _parse_group_cell("*,1")

    assert excinfo.value.reason == "control_wildcard_with_explicit_groups"


def test_a_repeated_wildcard_reports_the_less_specific_reason():
    """Documented, not fixed: "*,*" leaves no explicit tokens, so the wildcard-mix branch
    does not fire and it falls through to the generic invalid-value message."""
    with pytest.raises(EDInvalidGroupError) as excinfo:
        _parse_group_cell("*,*")

    assert excinfo.value.reason == "invalid_group_value"


def test_the_row_number_is_carried_into_the_error():
    with pytest.raises(EDInvalidGroupError) as excinfo:
        _parse_group_cell("bad", row_context=7)

    assert excinfo.value.offenders == [7]


def test_without_a_row_number_the_error_counts_no_rows():
    """Documented, not fixed: the offenders list is empty, so the message reads
    "in 0 rows" and names none."""
    with pytest.raises(EDInvalidGroupError) as excinfo:
        _parse_group_cell("bad")

    assert excinfo.value.offenders == []
    assert "0 row" in excinfo.value.message


def test_this_parser_does_not_enforce_the_row_type_rules():
    """A test row using "*" and a test row naming two groups are both accepted here;
    only the validator rejects them."""
    assert _parse_group_cell("*") == GROUP_WILDCARD
    assert _parse_group_cell("1,2") == frozenset({1, 2})


# --- Experiment ----------------------------------------------------------------

def test_an_invalid_type_is_rejected():
    with pytest.raises(EDInvalidTypeError):
        Experiment(["Experiment Name", "Type"], ["e1", "X"])


@pytest.mark.parametrize("bad", ["c", "t", "control", ""])
def test_the_type_check_is_case_sensitive_and_exact(bad):
    with pytest.raises(EDInvalidTypeError):
        Experiment(["Experiment Name", "Type"], ["e1", bad])


def test_no_group_column_means_no_group_spec():
    experiment = Experiment(["Experiment Name", "Type"], ["e1", "T"])

    assert experiment.group_spec is None


def test_a_group_column_is_parsed():
    experiment = Experiment(["Experiment Name", "Type", "Group"], ["e1", "T", "2"])

    assert experiment.group_spec == frozenset({2})


def test_a_short_row_silently_drops_its_trailing_columns():
    """dict(zip(...)) stops at the shorter sequence, so a row with fewer cells than the
    header loses the last ones without any error."""
    experiment = Experiment(["Experiment Name", "Type", "Bait"], ["e1", "T"])

    assert "Bait" not in experiment.attributes


def test_a_row_too_short_to_carry_a_type_raises_a_bare_key_error():
    """Documented, not fixed: the truncation happens before the Type check, so this is a
    KeyError rather than one of the file errors the GUI knows how to present."""
    with pytest.raises(KeyError):
        Experiment(["Experiment Name", "Type"], ["e1"])


# --- ExperimentalDesign construction -------------------------------------------

def test_a_missing_file_is_reported(tmp_path):
    with pytest.raises(EDFileNotFoundError):
        ExperimentalDesign(str(tmp_path / "absent.csv"))


def test_an_empty_file_is_reported(tmp_path):
    path = tmp_path / "ED.csv"
    path.write_text("")

    with pytest.raises(EDFileEmptyError):
        ExperimentalDesign(str(path))


def test_a_header_only_file_yields_an_empty_design(tmp_path):
    """Documented, not fixed: only a file with no header at all is reported as empty, so
    a design with no rows parses successfully and reports zero experiments."""
    design = ExperimentalDesign(_write_ed(tmp_path, []))

    assert design.name2experiment == {}
    assert (design.num_experiments, design.num_controls) == (0, 0)


@pytest.mark.parametrize("column", HEADER[:4])
def test_each_required_column_is_checked(tmp_path, column):
    header = [c for c in HEADER if c != column]
    rows = [["t1", "T", "BaitA", 1, "P_A"]]
    trimmed = [[v for h, v in zip(HEADER, rows[0]) if h != column]]

    with pytest.raises(EDMissingColumnError) as excinfo:
        ExperimentalDesign(_write_ed(tmp_path, trimmed, header=header))

    assert excinfo.value.missing_columns == [column]


def test_only_the_first_missing_column_is_reported(tmp_path):
    """The validator collects them all; this parser stops at the first, so a file missing
    several columns takes several corrections when it is read here."""
    with pytest.raises(EDMissingColumnError) as excinfo:
        ExperimentalDesign(_write_ed(tmp_path, [["t1"]], header=["Experiment Name"]))

    assert excinfo.value.missing_columns == ["Type"]


def test_a_bait_id_column_is_not_required(tmp_path):
    """Nothing here needs it, though write_CompPASS and the imputation paths do."""
    design = ExperimentalDesign(_write_ed(tmp_path, [["t1", "T", "BaitA", 1]],
                                          header=HEADER[:4]))

    assert design.num_experiments == 1


# --- counts --------------------------------------------------------------------

@pytest.fixture
def design(tmp_path):
    return ExperimentalDesign(_write_ed(tmp_path, [
        ["t1_1", "T", "BaitA", 1, "P_A"],
        ["t1_2", "T", "BaitA", 2, "P_A"],
        ["t2_1", "T", "BaitB", 1, "P_B"],
        ["c_1", "C", "Ctrl", 1, "P_C"],
    ]))


def test_experiments_are_counted_as_runs_not_baits(design):
    """Two baits over three runs reports three, which is what run.json and the GUI show."""
    assert design.num_experiments == 3


def test_controls_are_counted_separately(design):
    assert design.num_controls == 1


def test_a_non_numeric_replicate_on_a_test_row_is_rejected(tmp_path):
    with pytest.raises(EDInvalidReplicateError) as excinfo:
        ExperimentalDesign(_write_ed(tmp_path, [["t1", "T", "BaitA", "first", "P_A"]]))

    assert excinfo.value.invalid_rows == ["t1"]


def test_a_non_numeric_replicate_on_a_control_row_is_accepted(tmp_path):
    """Replicates are only parsed for test rows, so the same cell passes on a control."""
    design = ExperimentalDesign(_write_ed(tmp_path, [
        ["t1", "T", "BaitA", 1, "P_A"],
        ["c_1", "C", "Ctrl", "first", "P_C"],
    ]))

    assert (design.num_experiments, design.num_controls) == (1, 1)


def test_a_duplicate_experiment_name_silently_overwrites(tmp_path):
    """Documented, not fixed: experiments are keyed by name, so a repeated name keeps
    only the last row and the design reports one fewer experiment than the file lists.
    Only EDValidator rejects duplicates."""
    design = ExperimentalDesign(_write_ed(tmp_path, [
        ["t1", "T", "BaitA", 1, "P_A"],
        ["t1", "T", "BaitB", 1, "P_B"],
    ]))

    assert len(design.name2experiment) == 1
    assert design.name2experiment["t1"].attributes["Bait"] == "BaitB"


# --- grouping ------------------------------------------------------------------

GROUPED_HEADER = HEADER + ["Group"]


@pytest.fixture
def grouped(tmp_path):
    """Two test groups, one control in group 1 only, one universal control."""
    return ExperimentalDesign(_write_ed(tmp_path, [
        ["t1_1", "T", "BaitA", 1, "P_A", "1"],
        ["t1_2", "T", "BaitA", 2, "P_A", "1"],
        ["t2_1", "T", "BaitB", 1, "P_B", "2"],
        ["c_1", "C", "Ctrl1", 1, "P_C", "1"],
        ["c_all", "C", "CtrlAll", 1, "P_D", "*"],
    ], header=GROUPED_HEADER))


def test_an_ungrouped_design_is_not_grouped(design):
    assert design.is_grouped() is False
    assert design.get_groups() == []


def test_a_grouped_design_reports_its_groups(grouped):
    assert grouped.is_grouped() is True
    assert grouped.get_groups() == [1, 2]


def test_wildcard_controls_do_not_create_groups(grouped):
    """Groups come from test rows only, so a universal control adds none."""
    assert 0 not in grouped.get_groups()


def test_controls_alone_can_make_a_design_grouped_with_no_groups(tmp_path):
    """Documented, not fixed: is_grouped looks at every row but get_groups looks only at
    test rows, so a design whose only Group cells are on controls reports grouped mode
    with nothing to iterate."""
    only_controls_grouped = ExperimentalDesign(_write_ed(tmp_path, [
        ["t1_1", "T", "BaitA", 1, "P_A", ""],
        ["c_1", "C", "Ctrl", 1, "P_C", "1"],
    ], header=GROUPED_HEADER))

    assert only_controls_grouped.is_grouped() is True
    assert only_controls_grouped.get_groups() == []


def _names(experiments):
    return {e.attributes["Experiment Name"] for e in experiments}


def test_a_group_sees_its_own_test_experiments(grouped):
    tests, _ = grouped.get_experiments_for_group(1)

    assert _names(tests) == {"t1_1", "t1_2"}


def test_a_group_does_not_see_another_groups_tests(grouped):
    tests, _ = grouped.get_experiments_for_group(2)

    assert _names(tests) == {"t2_1"}


def test_a_universal_control_appears_in_every_group(grouped):
    """A control marked "*" is the background for every bait, so leaving it out of a
    group would remove that group's only shared background."""
    _, group_one = grouped.get_experiments_for_group(1)
    _, group_two = grouped.get_experiments_for_group(2)

    assert "c_all" in _names(group_one)
    assert "c_all" in _names(group_two)


def test_an_explicit_control_appears_only_in_its_group(grouped):
    _, group_one = grouped.get_experiments_for_group(1)
    _, group_two = grouped.get_experiments_for_group(2)

    assert _names(group_one) == {"c_1", "c_all"}
    assert _names(group_two) == {"c_all"}


def test_a_control_naming_several_groups_appears_in_each(tmp_path):
    multi = ExperimentalDesign(_write_ed(tmp_path, [
        ["t1_1", "T", "BaitA", 1, "P_A", "1"],
        ["t2_1", "T", "BaitB", 1, "P_B", "2"],
        ["c_12", "C", "Ctrl", 1, "P_C", "1,2"],
    ], header=GROUPED_HEADER))

    for group in (1, 2):
        _, controls = multi.get_experiments_for_group(group)
        assert _names(controls) == {"c_12"}


def test_an_ungrouped_test_row_belongs_to_no_group(tmp_path):
    """Validation rejects this file, but the CLI can reach the class directly; such a row
    is dropped from every group rather than defaulting into one."""
    partial = ExperimentalDesign(_write_ed(tmp_path, [
        ["t1_1", "T", "BaitA", 1, "P_A", "1"],
        ["t_none", "T", "BaitC", 1, "P_E", ""],
        ["c_all", "C", "Ctrl", 1, "P_C", "*"],
    ], header=GROUPED_HEADER))

    tests, _ = partial.get_experiments_for_group(1)

    assert _names(tests) == {"t1_1"}


def test_an_unknown_group_yields_only_universal_controls(grouped):
    tests, controls = grouped.get_experiments_for_group(99)

    assert tests == []
    assert _names(controls) == {"c_all"}
