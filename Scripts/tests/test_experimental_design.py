"""Tests for reading the experimental design file.

``ExperimentalDesign`` is what scoring reads; the GUI validates first, but the CLI
reaches this class directly.  ``get_experiments_for_group`` decides which experiments
each per-group SAINT run sees, so a prey's evidence depends on it.
"""

import csv

import pytest

from ed_exceptions import (
    EDInvalidGroupError,
    EDInvalidReplicateError,
    EDInvalidTypeError,
)
from experimental_design import (
    GROUP_WILDCARD,
    Experiment,
    ExperimentalDesign,
    _parse_group_cell,
)


HEADER = ["Experiment Name", "Type", "Bait", "Replicate", "Bait ID"]
GROUPED_HEADER = HEADER + ["Group"]


def _write_ed(tmp_path, rows, header=HEADER, name="ED.csv"):
    """Write an ED file through the csv module so a multi-group cell such as "1,2"
    is quoted rather than split into two columns."""
    path = tmp_path / name
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow(["" if c is None else c for c in row])
    return str(path)


def _names(experiments):
    return {e.attributes["Experiment Name"] for e in experiments}


# --- group cells -----------------------------------------------------------------

@pytest.mark.parametrize("raw, expected", [
    (None, None),
    ("   ", None),
    (" * ", GROUP_WILDCARD),
    ("1", frozenset({1})),
    (" 1 , 2 ", frozenset({1, 2})),
    ("2,2", frozenset({2})),
])
def test_group_cells_parse_to_a_spec(raw, expected):
    assert _parse_group_cell(raw) == expected


@pytest.mark.parametrize("raw", ["0", "1.0", "a", ","])
def test_malformed_group_cells_are_rejected(raw):
    with pytest.raises(EDInvalidGroupError) as excinfo:
        _parse_group_cell(raw)

    assert excinfo.value.reason == "invalid_group_value"


# --- Experiment ------------------------------------------------------------------

@pytest.mark.parametrize("bad", ["c", "t", "control", "X", ""])
def test_the_type_check_is_case_sensitive_and_exact(bad):
    with pytest.raises(EDInvalidTypeError):
        Experiment(["Experiment Name", "Type"], ["e1", bad])


# --- counts ----------------------------------------------------------------------

@pytest.fixture
def design(tmp_path):
    return ExperimentalDesign(_write_ed(tmp_path, [
        ["t1_1", "T", "BaitA", 1, "P_A"],
        ["t1_2", "T", "BaitA", 2, "P_A"],
        ["t2_1", "T", "BaitB", 1, "P_B"],
        ["c_1", "C", "Ctrl", 1, "P_C"],
    ]))


def test_experiments_are_counted_as_runs_not_baits(design):
    """Two baits over three runs reports three, which is what run.json and the GUI
    show; controls are counted separately."""
    assert (design.num_experiments, design.num_controls) == (3, 1)


def test_a_non_numeric_replicate_on_a_test_row_is_rejected(tmp_path):
    with pytest.raises(EDInvalidReplicateError) as excinfo:
        ExperimentalDesign(_write_ed(tmp_path, [["t1", "T", "BaitA", "first", "P_A"]]))

    assert excinfo.value.invalid_rows == ["t1"]


# --- grouping --------------------------------------------------------------------

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


def test_a_grouped_design_reports_the_groups_its_test_rows_declare(grouped):
    """Groups come from test rows only, so the universal control adds none."""
    assert grouped.is_grouped() is True
    assert grouped.get_groups() == [1, 2]


def test_a_group_sees_only_its_own_test_experiments(grouped):
    assert _names(grouped.get_experiments_for_group(1)[0]) == {"t1_1", "t1_2"}
    assert _names(grouped.get_experiments_for_group(2)[0]) == {"t2_1"}


def test_a_universal_control_appears_in_every_group(grouped):
    """A control marked "*" is the background for every bait."""
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
    """Validation rejects this file, but the CLI can reach the class directly; such a
    row is dropped from every group rather than defaulting into one."""
    partial = ExperimentalDesign(_write_ed(tmp_path, [
        ["t1_1", "T", "BaitA", 1, "P_A", "1"],
        ["t_none", "T", "BaitC", 1, "P_E", ""],
        ["c_all", "C", "Ctrl", 1, "P_C", "*"],
    ], header=GROUPED_HEADER))

    tests, _ = partial.get_experiments_for_group(1)

    assert _names(tests) == {"t1_1"}
