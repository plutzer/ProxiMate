"""Tests for the bait-vs-bait network comparison calculations and figures."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

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

@pytest.mark.parametrize("value, expected", [
    ("100.5|200.3|150.2", [100.5, 200.3, 150.2]),
    ("100.5|.|150.2", [100.5, 150.2]),          # "." marks a missing replicate
    (" 100.5 | | 150.2 ", [100.5, 150.2]),
    (None, []),
    (np.nan, []),
    ("", []),
])
def test_parse_intensity_string(value, expected):
    assert parse_intensity_string(value) == expected


def test_an_unparseable_entry_is_an_error():
    with pytest.raises(ValueError):
        parse_intensity_string("100.5|not_a_number|150.2")


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


@pytest.fixture
def volcano(dataset):
    return calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                  out_dir=dataset).set_index("Prey.ID")


# --- load_and_filter_bait_data -------------------------------------------------

def test_returns_only_rows_passing_all_thresholds(dataset):
    filtered = load_and_filter_bait_data("ds", "BaitB", THRESHOLDS, out_dir=dataset)

    # BaitB's P2 has SaintScore 0.1 and is dropped; P1 remains.
    assert list(filtered["Prey.ID"]) == ["P1"]


# --- calculate_volcano_data ----------------------------------------------------

def test_every_prey_is_returned_with_a_status(volcano):
    assert set(volcano.index) == {"P1", "P2", "P3"}
    assert volcano.loc["P1", "status"] == "shared"
    assert volcano.loc["P2", "status"] == "shared"
    # P3 appears only under BaitA.
    assert volcano.loc["P3", "status"] == "a_only"


def test_fold_change_is_log2_of_the_mean_intensity_ratio(volcano):
    # P1: BaitA mean 200, BaitB mean 50.
    assert volcano.loc["P1", "mean_intensity_a"] == pytest.approx(200.0)
    assert volcano.loc["P1", "mean_intensity_b"] == pytest.approx(50.0)
    assert volcano.loc["P1", "log2_fc_ratio"] == pytest.approx(np.log2(4.0))
    # P2 is unchanged between baits.
    assert volcano.loc["P2", "log2_fc_ratio"] == pytest.approx(0.0)


def test_categories_follow_threshold_passing(volcano):
    assert volcano.loc["P1", "category"] == "Both"
    # BaitB's P2 fails SaintScore, so it passes for A only.
    assert volcano.loc["P2", "category"] == "Network A only"
    # Side-panel rows are categorized the same way.
    assert volcano.loc["P3", "category"] == "Network A only"


def test_categories_agree_with_the_filtered_networks(dataset, volcano):
    """Whatever the thresholds, a prey drawn as passing must be one the filtered
    network also keeps."""
    kept_a = set(load_and_filter_bait_data("ds", "BaitA", THRESHOLDS,
                                           out_dir=dataset)["Prey.ID"])

    drawn_for_a = set(volcano[volcano["category"].isin(["Both", "Network A only"])].index)
    assert drawn_for_a <= kept_a


# --- t-test scale and multiple-testing correction ------------------------------

def test_ttest_runs_on_log2_intensities(volcano):
    p_log = ttest_ind(np.log2([100.0, 300.0]), np.log2([25.0, 75.0])).pvalue
    p_raw = ttest_ind([100.0, 300.0], [25.0, 75.0]).pvalue

    assert p_log != pytest.approx(p_raw)   # the fixture distinguishes the scales
    assert volcano.loc["P1", "pval"] == pytest.approx(p_log)


def test_pvalues_are_bh_adjusted(volcano):
    p1 = ttest_ind(np.log2([100.0, 300.0]), np.log2([25.0, 75.0])).pvalue
    p2 = ttest_ind(np.log2([48.0, 52.0]), np.log2([49.0, 51.0])).pvalue
    expected = multipletests([p1, p2], method="fdr_bh")[1]

    assert volcano.loc["P1", "pval_adj"] == pytest.approx(expected[0])
    assert volcano.loc["P2", "pval_adj"] == pytest.approx(expected[1])
    assert volcano.loc["P1", "neg_log10_pval"] == pytest.approx(-np.log10(expected[0]))


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

@pytest.mark.parametrize("present, status", [("BaitA", "a_only"), ("BaitB", "b_only")])
def test_prey_absent_from_one_bait_goes_to_a_side_panel(tmp_path, present, status):
    """A prey with intensity under one bait only carries no fold change or p-value;
    it is routed to that bait's side panel with y = -log10 of its BFDR there."""
    absent = "BaitB" if present == "BaitA" else "BaitA"
    prefix = {"BaitA": "a", "BaitB": "b"}
    interaction_rows = [
        (f"{prefix[present]}_1", present, "P1", 400.0),
        (f"{prefix[present]}_2", present, "P1", 600.0),
        (f"{prefix[absent]}_1", absent, "P1", 0.0),
        (f"{prefix[absent]}_2", absent, "P1", 0.0),
    ]
    score_rows = [_score(present, "P1", BFDR=0.001), _score(absent, "P1", BFDR=0.5)]
    out_dir = _write_dataset(tmp_path, score_rows, interaction_rows)

    volcano = calculate_volcano_data("ds", "BaitA", "BaitB", THRESHOLDS, THRESHOLDS,
                                     out_dir=out_dir).set_index("Prey.ID")

    assert volcano.loc["P1", "status"] == status
    for column in ["fc_ratio", "log2_fc_ratio", "pval", "neg_log10_pval"]:
        assert np.isnan(volcano.loc["P1", column])
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


# --- figures -------------------------------------------------------------------

def _volcano_frame():
    """Hand-built calculate_volcano_data output: three shared preys in two categories
    and one prey on each side panel."""
    def shared(prey, category, log2_fc, adj):
        return {"Prey.ID": prey, "First_Prey_Gene": prey.lower(), "status": "shared",
                "category": category, "mean_intensity_a": 200.0, "mean_intensity_b": 50.0,
                "fc_ratio": 2 ** log2_fc, "log2_fc_ratio": log2_fc, "pval": adj / 2,
                "pval_adj": adj, "neg_log10_pval": -np.log10(adj), "neg_log10_bfdr": np.nan}

    def side(prey, status, category, y):
        return {"Prey.ID": prey, "First_Prey_Gene": prey.lower(), "status": status,
                "category": category, "mean_intensity_a": 500.0, "mean_intensity_b": 0.0,
                "fc_ratio": np.nan, "log2_fc_ratio": np.nan, "pval": np.nan,
                "pval_adj": np.nan, "neg_log10_pval": np.nan, "neg_log10_bfdr": y}

    return pd.DataFrame([
        shared("P1", "Both", 2.0, 0.02),
        shared("P2", "Neither", 0.0, 0.9),
        shared("P5", "Both", -1.0, 0.04),
        side("P3", "a_only", "Network A only", 2.0),
        side("P4", "b_only", "Network B only", 3.0),
    ])


def test_plotly_volcano_has_three_panels_with_the_sides_matching_fold_direction():
    """Positive log2 FC means higher in bait A, so the A-only strip sits on the
    right and the B-only strip on the left."""
    fig = create_volcano_plot(_volcano_frame(), "BaitA", "BaitB")

    assert fig.layout.xaxis2 is not None and fig.layout.xaxis3 is not None
    assert {trace.xaxis for trace in fig.data} >= {"x", "x2", "x3"}
    titles = {a.text: a.x for a in fig.layout.annotations if a.text.startswith("Only in")}
    assert titles["Only in BaitB"] < 0.5 < titles["Only in BaitA"]


def test_central_panel_holds_one_trace_per_category_with_its_preys():
    fig = create_volcano_plot(_volcano_frame(), "BaitA", "BaitB")

    central = {trace.name: len(trace.x) for trace in fig.data if trace.xaxis == "x2"}
    assert central == {"Both": 2, "Neither": 1}
    left = {trace.name: len(trace.x) for trace in fig.data if trace.xaxis == "x"}
    right = {trace.name: len(trace.x) for trace in fig.data if trace.xaxis == "x3"}
    assert left == {"Network B only": 1}
    assert right == {"Network A only": 1}


def test_matplotlib_volcano_has_three_axes():
    fig = create_volcano_plot_matplotlib(_volcano_frame(), "BaitA", "BaitB")
    assert len(fig.axes) == 3


def test_venn_region_counts_appear_in_the_diagram():
    set_a = {f"A{i}" for i in range(10)} | {f"S{i}" for i in range(5)}
    set_b = {f"B{i}" for i in range(3)} | {f"S{i}" for i in range(5)}

    fig = create_venn_diagram_matplotlib(set_a, set_b, "BaitA", "BaitB")

    assert {"10", "5", "3"} <= {t.get_text() for ax in fig.axes for t in ax.texts}
