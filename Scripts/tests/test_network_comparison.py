"""Tests for the bait-vs-bait network comparison calculations and figures.

Covers the calculation half of GUI/network_comparison.py in depth; the plotly and
matplotlib figure builders get structural smoke tests only.
"""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import ttest_ind

from network_comparison import (
    calculate_volcano_data,
    create_venn_diagram_matplotlib,
    create_volcano_plot,
    create_volcano_plot_matplotlib,
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
        {"Experiment Name": "a_3", "Type": "T", "Bait": "BaitA", "Replicate": 3,
         "Bait ID": "IDA"},
        {"Experiment Name": "b_1", "Type": "T", "Bait": "BaitB", "Replicate": 1,
         "Bait ID": "IDB"},
        {"Experiment Name": "b_2", "Type": "T", "Bait": "BaitB", "Replicate": 2,
         "Bait ID": "IDB"},
        {"Experiment Name": "b_3", "Type": "T", "Bait": "BaitB", "Replicate": 3,
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

def test_every_prey_is_returned_with_a_status(dataset):
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset).set_index("Prey.ID")

    assert set(volcano.index) == {"P1", "P2", "P3"}
    assert volcano.loc["P1", "status"] == "shared"
    assert volcano.loc["P2", "status"] == "shared"
    # P3 appears only under BaitA.
    assert volcano.loc["P3", "status"] == "a_only"


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
    # Side-panel rows are categorized the same way.
    assert volcano.loc["P3", "category"] == "Network A only"


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


# --- t-test scale and multiple-testing correction ------------------------------

def test_ttest_runs_on_log2_intensities(dataset):
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset).set_index("Prey.ID")

    p_log = ttest_ind(np.log2([100.0, 300.0]), np.log2([25.0, 75.0])).pvalue
    p_raw = ttest_ind([100.0, 300.0], [25.0, 75.0]).pvalue

    assert p_log != pytest.approx(p_raw)   # the fixture distinguishes the scales
    assert volcano.loc["P1", "pval"] == pytest.approx(p_log)


def test_pvalues_are_bh_adjusted(dataset):
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset).set_index("Prey.ID")

    p1 = ttest_ind(np.log2([100.0, 300.0]), np.log2([25.0, 75.0])).pvalue
    p2 = ttest_ind(np.log2([48.0, 52.0]), np.log2([49.0, 51.0])).pvalue

    # BH over the two tested preys, written out by hand: the largest p keeps its
    # value; the smaller becomes min(p * n / rank, next adjusted value).
    p_small, p_large = sorted([p1, p2])
    adj_large = min(p_large, 1.0)
    adj_small = min(p_small * 2.0, adj_large)
    expected = {p1: adj_large if p1 == p_large else adj_small,
                p2: adj_large if p2 == p_large else adj_small}

    assert volcano.loc["P1", "pval_adj"] == pytest.approx(expected[p1])
    assert volcano.loc["P2", "pval_adj"] == pytest.approx(expected[p2])
    assert volcano.loc["P1", "neg_log10_pval"] == pytest.approx(-np.log10(expected[p1]))
    assert volcano.loc["P2", "neg_log10_pval"] == pytest.approx(-np.log10(expected[p2]))


def test_zero_replicates_are_excluded_from_the_ttest(tmp_path):
    """A zero intensity is a non-detection: it is dropped before the log2 t-test
    rather than passed through log2(0)."""
    interaction_rows = [
        ("a_1", "BaitA", "P1", 0.0), ("a_2", "BaitA", "P1", 100.0),
        ("a_3", "BaitA", "P1", 400.0),
        ("b_1", "BaitB", "P1", 40.0), ("b_2", "BaitB", "P1", 60.0),
    ]
    score_rows = [_score("BaitA", "P1"), _score("BaitB", "P1")]
    out_dir = _write_dataset(tmp_path, score_rows, interaction_rows)

    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=out_dir).set_index("Prey.ID")

    expected = ttest_ind(np.log2([100.0, 400.0]), np.log2([40.0, 60.0])).pvalue
    assert volcano.loc["P1", "pval"] == pytest.approx(expected)
    assert np.isfinite(volcano.loc["P1", "neg_log10_pval"])


def test_too_few_usable_replicates_yield_nan_pvalues(tmp_path):
    """With <2 nonzero replicates on a side there is no test; the p-values are NaN
    (never a fake 1.0, which would distort the BH ranking) and the prey is drawn
    at the bottom of the volcano."""
    interaction_rows = [
        ("a_1", "BaitA", "P1", 0.0), ("a_2", "BaitA", "P1", 0.0),
        ("a_3", "BaitA", "P1", 100.0),
        ("b_1", "BaitB", "P1", 40.0), ("b_2", "BaitB", "P1", 60.0),
    ]
    score_rows = [_score("BaitA", "P1"), _score("BaitB", "P1")]
    out_dir = _write_dataset(tmp_path, score_rows, interaction_rows)

    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=out_dir).set_index("Prey.ID")

    assert volcano.loc["P1", "status"] == "shared"      # nonzero mean on both sides
    assert np.isnan(volcano.loc["P1", "pval"])
    assert np.isnan(volcano.loc["P1", "pval_adj"])
    assert volcano.loc["P1", "neg_log10_pval"] == 0.0


# --- presence/absence side panels ----------------------------------------------

def test_prey_absent_from_one_bait_goes_to_a_side_panel(tmp_path):
    """A prey with intensity only under bait A carries no fold change or p-value;
    it is routed to the A-only side panel with y = -log10 of its BFDR under A."""
    interaction_rows = [
        ("a_1", "BaitA", "P1", 400.0), ("a_2", "BaitA", "P1", 600.0),
        ("b_1", "BaitB", "P1", 0.0),   ("b_2", "BaitB", "P1", 0.0),
    ]
    score_rows = [_score("BaitA", "P1"), _score("BaitB", "P1")]
    out_dir = _write_dataset(tmp_path, score_rows, interaction_rows)

    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=out_dir).set_index("Prey.ID")

    assert volcano.loc["P1", "status"] == "a_only"
    assert volcano.loc["P1", "mean_intensity_a"] == pytest.approx(500.0)
    assert volcano.loc["P1", "mean_intensity_b"] == pytest.approx(0.0)
    assert np.isnan(volcano.loc["P1", "log2_fc_ratio"])
    assert np.isnan(volcano.loc["P1", "fc_ratio"])
    assert np.isnan(volcano.loc["P1", "pval"])
    assert np.isnan(volcano.loc["P1", "neg_log10_pval"])
    # BFDR 0.01 under BaitA -> y = 2.0
    assert volcano.loc["P1", "neg_log10_bfdr"] == pytest.approx(2.0)


def test_b_only_prey_takes_bfdr_from_bait_b(tmp_path):
    interaction_rows = [
        ("b_1", "BaitB", "P1", 400.0), ("b_2", "BaitB", "P1", 600.0),
        ("a_1", "BaitA", "P2", 90.0),  ("a_2", "BaitA", "P2", 110.0),
        ("b_1", "BaitB", "P2", 95.0),  ("b_2", "BaitB", "P2", 105.0),
    ]
    score_rows = [_score("BaitA", "P1", BFDR=0.5), _score("BaitB", "P1", BFDR=0.001),
                  _score("BaitA", "P2"), _score("BaitB", "P2")]
    out_dir = _write_dataset(tmp_path, score_rows, interaction_rows)

    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=out_dir).set_index("Prey.ID")

    assert volcano.loc["P1", "status"] == "b_only"
    assert volcano.loc["P1", "neg_log10_bfdr"] == pytest.approx(3.0)


def test_side_panel_bfdr_edge_cases(tmp_path):
    """BFDR 0 is floored at 1e-4 (y capped at 4); a missing or NaN BFDR draws at
    the panel floor (y = 0)."""
    interaction_rows = [
        ("a_1", "BaitA", "P1", 400.0), ("a_2", "BaitA", "P1", 600.0),
        ("a_1", "BaitA", "P2", 400.0), ("a_2", "BaitA", "P2", 600.0),
        ("a_1", "BaitA", "P3", 400.0), ("a_2", "BaitA", "P3", 600.0),
        ("b_1", "BaitB", "P4", 90.0),  ("b_2", "BaitB", "P4", 110.0),
        ("a_1", "BaitA", "P4", 95.0),  ("a_2", "BaitA", "P4", 105.0),
    ]
    score_rows = [_score("BaitA", "P1", BFDR=0.0),
                  _score("BaitA", "P2", BFDR=np.nan),
                  # P3 has no annotated_scores row at all
                  _score("BaitA", "P4"), _score("BaitB", "P4")]
    out_dir = _write_dataset(tmp_path, score_rows, interaction_rows)

    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=out_dir).set_index("Prey.ID")

    assert volcano.loc["P1", "neg_log10_bfdr"] == pytest.approx(4.0)
    assert volcano.loc["P2", "neg_log10_bfdr"] == 0.0
    assert volcano.loc["P3", "neg_log10_bfdr"] == 0.0


def test_column_contract_between_shared_and_side_rows(dataset):
    """Shared rows carry the volcano columns and NaN BFDR; side rows carry BFDR
    and NaN volcano columns."""
    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=dataset).set_index("Prey.ID")

    expected_columns = {"First_Prey_Gene", "status", "category",
                        "mean_intensity_a", "mean_intensity_b",
                        "fc_ratio", "log2_fc_ratio",
                        "pval", "pval_adj", "neg_log10_pval", "neg_log10_bfdr"}
    assert expected_columns <= set(volcano.columns)

    shared = volcano[volcano["status"] == "shared"]
    side = volcano[volcano["status"] != "shared"]
    assert shared["neg_log10_bfdr"].isna().all()
    for column in ["fc_ratio", "log2_fc_ratio", "pval", "pval_adj", "neg_log10_pval"]:
        assert side[column].isna().all()


# --- figure smoke tests --------------------------------------------------------

def _volcano_frame():
    """Hand-built calculate_volcano_data output with shared and side-panel rows."""
    return pd.DataFrame([
        {"Prey.ID": "P1", "First_Prey_Gene": "p1", "status": "shared",
         "category": "Both", "mean_intensity_a": 200.0, "mean_intensity_b": 50.0,
         "fc_ratio": 4.0, "log2_fc_ratio": 2.0, "pval": 0.01, "pval_adj": 0.02,
         "neg_log10_pval": -np.log10(0.02), "neg_log10_bfdr": np.nan},
        {"Prey.ID": "P2", "First_Prey_Gene": "p2", "status": "shared",
         "category": "Neither", "mean_intensity_a": 50.0, "mean_intensity_b": 50.0,
         "fc_ratio": 1.0, "log2_fc_ratio": 0.0, "pval": 0.9, "pval_adj": 0.9,
         "neg_log10_pval": -np.log10(0.9), "neg_log10_bfdr": np.nan},
        {"Prey.ID": "P3", "First_Prey_Gene": "p3", "status": "a_only",
         "category": "Network A only", "mean_intensity_a": 500.0,
         "mean_intensity_b": 0.0, "fc_ratio": np.nan, "log2_fc_ratio": np.nan,
         "pval": np.nan, "pval_adj": np.nan, "neg_log10_pval": np.nan,
         "neg_log10_bfdr": 2.0},
        {"Prey.ID": "P4", "First_Prey_Gene": "p4", "status": "b_only",
         "category": "Network B only", "mean_intensity_a": 0.0,
         "mean_intensity_b": 300.0, "fc_ratio": np.nan, "log2_fc_ratio": np.nan,
         "pval": np.nan, "pval_adj": np.nan, "neg_log10_pval": np.nan,
         "neg_log10_bfdr": 3.0},
    ])


def test_plotly_volcano_has_three_panels():
    fig = create_volcano_plot(_volcano_frame(), "BaitA", "BaitB")

    assert fig.layout.xaxis2 is not None and fig.layout.xaxis3 is not None
    # Traces land on all three panels.
    assert {trace.xaxis for trace in fig.data} >= {"x", "x2", "x3"}


def test_plotly_volcano_panel_sides_match_fold_change_direction():
    """Positive log2 FC means higher in bait A, so the A-only strip sits on the
    right and the B-only strip on the left."""
    fig = create_volcano_plot(_volcano_frame(), "BaitA", "BaitB")

    titles = {a.text: a.x for a in fig.layout.annotations if a.text.startswith("Only in")}
    assert titles["Only in BaitB"] < 0.5 < titles["Only in BaitA"]


def test_plotly_volcano_title_is_centered():
    fig = create_volcano_plot(_volcano_frame(), "BaitA", "BaitB")
    assert fig.layout.title.x == 0.5


def test_plotly_volcano_empty_frame_returns_placeholder():
    fig = create_volcano_plot(pd.DataFrame(), "BaitA", "BaitB")

    assert len(fig.data) == 0
    assert len(fig.layout.annotations) == 1


def test_plotly_volcano_jitter_is_reproducible():
    frame = _volcano_frame()
    fig1 = create_volcano_plot(frame, "BaitA", "BaitB")
    fig2 = create_volcano_plot(frame, "BaitA", "BaitB")

    side1 = [tuple(t.x) for t in fig1.data if t.xaxis in ("x", "x3")]
    side2 = [tuple(t.x) for t in fig2.data if t.xaxis in ("x", "x3")]
    assert side1 == side2


def test_matplotlib_volcano_has_three_axes():
    fig = create_volcano_plot_matplotlib(_volcano_frame(), "BaitA", "BaitB")
    assert len(fig.axes) == 3


def _venn_texts(fig):
    return {t.get_text() for ax in fig.axes for t in ax.texts}


def test_venn_region_counts_appear_in_the_diagram():
    set_a = {f"A{i}" for i in range(10)} | {f"S{i}" for i in range(5)}
    set_b = {f"B{i}" for i in range(3)} | {f"S{i}" for i in range(5)}

    fig = create_venn_diagram_matplotlib(set_a, set_b, "BaitA", "BaitB")

    assert {"10", "5", "3"} <= _venn_texts(fig)


def test_venn_circles_are_area_proportional():
    set_a = {f"A{i}" for i in range(100)}
    set_b = {f"B{i}" for i in range(10)}

    fig = create_venn_diagram_matplotlib(set_a, set_b, "BaitA", "BaitB")

    widths = sorted(p.get_extents().width for p in fig.axes[0].patches)
    assert widths[-1] > widths[0] * 1.5


def test_venn_tolerates_an_empty_intersection():
    fig = create_venn_diagram_matplotlib({"A1", "A2"}, {"B1"}, "BaitA", "BaitB")
    assert {"2", "1"} <= _venn_texts(fig)


def test_venn_tolerates_one_empty_set():
    fig = create_venn_diagram_matplotlib({"A1", "A2"}, set(), "BaitA", "BaitB")
    assert "2" in _venn_texts(fig)


def test_venn_with_both_sets_empty_returns_a_placeholder():
    fig = create_venn_diagram_matplotlib(set(), set(), "BaitA", "BaitB")
    assert len(fig.axes) >= 1     # placeholder message, no venn artists
