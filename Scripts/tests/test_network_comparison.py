"""Tests for the bait-vs-bait network comparison calculations.

Covers the calculation half of GUI/network_comparison.py; the plotly and matplotlib
figure builders are presentation and are out of scope.
"""

import numpy as np
import pandas as pd
import pytest

from network_comparison import (
    calculate_volcano_data,
    load_and_filter_bait_data,
    parse_intensity_string,
)


THRESHOLDS = {"SaintScore": 0.7, "BFDR": 0.05, "WD": 1.0, "WDFDR": 0.05}


# --- parse_intensity_string ----------------------------------------------------

def test_parses_pipe_delimited_values():
    assert parse_intensity_string("100.5|200.3|150.2") == [100.5, 200.3, 150.2]


def test_missing_values_are_dropped():
    """A "." marks a missing replicate and is skipped rather than read as zero."""
    assert parse_intensity_string("100.5|.|150.2") == [100.5, 150.2]


def test_whitespace_and_empty_entries_are_ignored():
    assert parse_intensity_string(" 100.5 | | 150.2 ") == [100.5, 150.2]


def test_unparseable_entries_are_skipped():
    assert parse_intensity_string("100.5|not_a_number|150.2") == [100.5, 150.2]


@pytest.mark.parametrize("value", [None, np.nan])
def test_null_input_yields_an_empty_list(value):
    assert parse_intensity_string(value) == []


def test_empty_string_yields_an_empty_list():
    assert parse_intensity_string("") == []


# --- dataset fixture -----------------------------------------------------------

def _write_dataset(tmp_path, score_rows, interaction_rows, name="ds"):
    """Lay out the interaction.txt / ED.csv / annotated_scores.csv trio the
    comparison functions read from <out_dir>/<dataset>/."""
    root = tmp_path / name
    root.mkdir(parents=True)

    with open(root / "interaction.txt", "w", newline="") as handle:
        for experiment, bait, prey, intensity in interaction_rows:
            handle.write(f"{experiment}\t{bait}\t{prey}\t{intensity}\n")

    pd.DataFrame([
        {"Experiment Name": "a_1", "Type": "T", "Bait": "BaitA", "Replicate": 1,
         "Bait ID": "IDA"},
        {"Experiment Name": "a_2", "Type": "T", "Bait": "BaitA", "Replicate": 2,
         "Bait ID": "IDA"},
        {"Experiment Name": "b_1", "Type": "T", "Bait": "BaitB", "Replicate": 1,
         "Bait ID": "IDB"},
        {"Experiment Name": "b_2", "Type": "T", "Bait": "BaitB", "Replicate": 2,
         "Bait ID": "IDB"},
        {"Experiment Name": "c_1", "Type": "C", "Bait": "Ctrl", "Replicate": 1,
         "Bait ID": "IDC"},
    ]).to_csv(root / "ED.csv", index=False)

    pd.DataFrame(score_rows).to_csv(root / "annotated_scores.csv", index=False)
    return str(tmp_path)


def _score(bait, prey, **overrides):
    row = {"Experiment.ID": bait, "Prey.ID": prey, "First_Prey_Gene": prey.lower(),
           "SaintScore": 0.9, "BFDR": 0.01, "WD": 5.0, "WDFDR": 0.01}
    row.update(overrides)
    return row


@pytest.fixture
def dataset(tmp_path):
    """Two baits sharing preys P1 and P2, plus P3 seen only under BaitA."""
    interaction_rows = [
        ("a_1", "BaitA", "P1", 100.0), ("a_2", "BaitA", "P1", 300.0),
        ("a_1", "BaitA", "P2", 48.0),  ("a_2", "BaitA", "P2", 52.0),
        ("a_1", "BaitA", "P3", 10.0),  ("a_2", "BaitA", "P3", 10.0),
        ("b_1", "BaitB", "P1", 25.0),  ("b_2", "BaitB", "P1", 75.0),
        ("b_1", "BaitB", "P2", 49.0),  ("b_2", "BaitB", "P2", 51.0),
        ("c_1", "Ctrl",  "P1", 1.0),
    ]
    score_rows = [_score("BaitA", "P1"), _score("BaitA", "P2"), _score("BaitA", "P3"),
                  _score("BaitB", "P1"), _score("BaitB", "P2", SaintScore=0.1)]
    return _write_dataset(tmp_path, score_rows, interaction_rows)


# --- load_and_filter_bait_data -------------------------------------------------

def test_missing_scores_file_yields_an_empty_frame(tmp_path):
    assert load_and_filter_bait_data("nope", "BaitA", THRESHOLDS,
                                     out_dir=str(tmp_path)).empty


def test_unknown_bait_yields_an_empty_frame(dataset):
    assert load_and_filter_bait_data("ds", "NotABait", THRESHOLDS,
                                     out_dir=dataset).empty


def test_returns_only_rows_passing_all_thresholds(dataset):
    filtered = load_and_filter_bait_data("ds", "BaitB", THRESHOLDS, out_dir=dataset)

    # BaitB's P2 has SaintScore 0.1 and is dropped; P1 remains.
    assert list(filtered["Prey.ID"]) == ["P1"]


# --- calculate_volcano_data ----------------------------------------------------

def test_compares_only_preys_seen_under_both_baits(dataset):
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset)

    # P3 appears only under BaitA.
    assert set(volcano["Prey.ID"]) == {"P1", "P2"}


def test_fold_change_is_log2_of_the_mean_intensity_ratio(dataset):
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset).set_index("Prey.ID")

    # P1: BaitA mean 200, BaitB mean 50.
    assert volcano.loc["P1", "mean_intensity_a"] == pytest.approx(200.0)
    assert volcano.loc["P1", "mean_intensity_b"] == pytest.approx(50.0)
    assert volcano.loc["P1", "log2_fc_ratio"] == pytest.approx(np.log2(4.0))
    # P2 is unchanged between baits.
    assert volcano.loc["P2", "log2_fc_ratio"] == pytest.approx(0.0)


def test_missing_dataset_yields_an_empty_frame(tmp_path):
    assert calculate_volcano_data("nope", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                  out_dir=str(tmp_path)).empty


def test_categories_follow_threshold_passing(dataset):
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset).set_index("Prey.ID")

    assert volcano.loc["P1", "category"] == "Both"
    # BaitB's P2 fails SaintScore, so it passes for A only.
    assert volcano.loc["P2", "category"] == "Network A only"


def test_missing_wdfdr_fails_the_category_threshold(tmp_path):
    """Categories use apply_score_thresholds, which reads a NaN WDFDR as 1.0.

    Scoring with 0 CompPASS iterations leaves WDFDR NaN.  Treating that as passing
    would draw a prey as a hit on the volcano while load_and_filter_bait_data
    excluded it from the very same network.
    """
    interaction_rows = [
        ("a_1", "BaitA", "P1", 90.0), ("a_2", "BaitA", "P1", 110.0),
        ("b_1", "BaitB", "P1", 95.0), ("b_2", "BaitB", "P1", 105.0),
    ]
    score_rows = [_score("BaitA", "P1", WDFDR=np.nan),
                  _score("BaitB", "P1", WDFDR=np.nan)]
    out_dir = _write_dataset(tmp_path, score_rows, interaction_rows)

    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=out_dir).set_index("Prey.ID")
    filtered = load_and_filter_bait_data("ds", "BaitA", THRESHOLDS, out_dir=out_dir)

    assert volcano.loc["P1", "category"] == "Neither"
    assert filtered.empty          # the two agree


def test_categories_agree_with_the_filtered_networks(dataset):
    """Whatever the thresholds, a prey drawn as passing must be one the filtered
    network also keeps."""
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset)
    kept_a = set(load_and_filter_bait_data("ds", "BaitA", THRESHOLDS,
                                           out_dir=dataset)["Prey.ID"])

    drawn_for_a = set(volcano[volcano["category"].isin(["Both", "Network A only"])]["Prey.ID"])
    assert drawn_for_a <= kept_a


# --- documented behavior, not fixed --------------------------------------------

@pytest.mark.filterwarnings("ignore:Precision loss occurred:RuntimeWarning")
def test_prey_absent_from_one_bait_plots_at_the_origin(tmp_path):
    """When either mean intensity is 0 the fold change is forced to 0, so a prey lost
    entirely from one network is drawn at log2 fold change 0 -- the same place as a
    prey with no change at all."""
    interaction_rows = [
        ("a_1", "BaitA", "P1", 400.0), ("a_2", "BaitA", "P1", 600.0),
        ("b_1", "BaitB", "P1", 0.0),   ("b_2", "BaitB", "P1", 0.0),
    ]
    score_rows = [_score("BaitA", "P1"), _score("BaitB", "P1")]
    out_dir = _write_dataset(tmp_path, score_rows, interaction_rows)

    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=out_dir).set_index("Prey.ID")

    assert volcano.loc["P1", "mean_intensity_a"] == pytest.approx(500.0)
    assert volcano.loc["P1", "mean_intensity_b"] == pytest.approx(0.0)
    assert volcano.loc["P1", "log2_fc_ratio"] == 0.0
    assert volcano.loc["P1", "fc_ratio"] == 0


def test_pvalues_carry_no_multiple_testing_correction(dataset):
    """One t-test per shared prey, reported raw; the frame has no adjusted column."""
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset)

    assert "pval" in volcano.columns
    assert not any("adj" in c.lower() or "fdr" in c.lower() for c in volcano.columns)
