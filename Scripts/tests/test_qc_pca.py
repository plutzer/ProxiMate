"""Tests for the PCA preprocessing pipeline and the experiment/prey PCA figures.

Covers prepare_pca_matrix and its imputation/normalization options, the helpers
detection_counts / reduce_categorical / load_pca_metadata, the plotly figures in
GUI/QC_plots.py, and the matplotlib export twins in GUI/plot_exports.py.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.figure
import numpy as np
import pandas as pd
import pytest

import plot_exports
from QC_plots import (
    detection_counts,
    load_pca_metadata,
    pca_plot,
    prepare_pca_matrix,
    prey_pca_plot,
    reduce_categorical,
)


EXPERIMENTS = [
    ("a_1", "BaitA", "T"),
    ("a_2", "BaitA", "T"),
    ("b_1", "BaitB", "T"),
    ("b_2", "BaitB", "C"),
]

# Prey -> intensity per experiment; None = not observed, 0 = observed as zero.
PREY_VALUES = {
    "P1": [100.0, 110.0, 200.0, 210.0],   # complete
    "P2": [50.0, 55.0, 60.0, 65.0],       # complete
    "P3": [80.0, 90.0, 100.0, 110.0],     # complete
    "P4": [30.0, 35.0, None, None],       # 2 of 4 detections
    "P5": [0.0, 40.0, 45.0, 50.0],        # zero counts as undetected -> 3 of 4
    "P6": [20.0, None, None, None],       # 1 of 4 detections
}


def _write_dataset(tmp_path, prey_values=PREY_VALUES, experiments=EXPERIMENTS):
    rows = []
    for prey, values in prey_values.items():
        for (exp, bait, _), value in zip(experiments, values):
            if value is not None:
                rows.append((exp, bait, prey, value))
    interaction = tmp_path / "interaction.txt"
    with open(interaction, "w", newline="") as fh:
        fh.write("Experiment\tBait\tPrey\tIntensity\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    ed = tmp_path / "ED.csv"
    pd.DataFrame(
        [{"Experiment Name": exp, "Type": typ, "Bait": bait, "Replicate": i + 1}
         for i, (exp, bait, typ) in enumerate(experiments)]
    ).to_csv(ed, index=False)
    return str(interaction), str(ed)


def _old_pipeline(interaction_path):
    """The preprocessing exactly as pca_plot has always done it."""
    df = pd.read_csv(interaction_path, sep="\t", header=0)
    df.columns = ["Experiment", "BaitName", "Prey", "Intensity"]
    data = df.pivot(index="Prey", columns="Experiment", values="Intensity")
    data = data.replace(0, np.nan)
    data = data.dropna(thresh=len(data.columns) * 0.5)
    data = data.apply(lambda row: row.fillna(row.min()), axis=1)
    data = data.apply(lambda row: (row - row.mean()) / row.std(), axis=1)
    return data


# --- prepare_pca_matrix --------------------------------------------------------

def test_defaults_reproduce_the_original_pipeline(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    pd.testing.assert_frame_equal(matrix, _old_pipeline(interaction))


def test_default_detection_filter_drops_sparse_preys(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    assert set(matrix.index) == {"P1", "P2", "P3", "P4", "P5"}


def test_zero_imputation_fills_missing_with_zero(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction, imputation="zero",
                                normalization="none")
    assert matrix.loc["P4", "b_1"] == 0.0
    assert matrix.loc["P5", "a_1"] == 0.0
    assert not matrix.isna().any().any()


def test_drop_imputation_keeps_only_complete_preys(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction, imputation="drop",
                                normalization="none")
    assert set(matrix.index) == {"P1", "P2", "P3"}
    assert not matrix.isna().any().any()


def test_row_min_imputation_fills_with_row_minimum(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction, normalization="none")
    assert matrix.loc["P4", "b_1"] == 30.0
    assert matrix.loc["P5", "a_1"] == 40.0


def test_none_normalization_returns_imputed_values(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction, normalization="none")
    assert matrix.loc["P1", "b_2"] == 210.0


def test_log2_zscore_is_zscore_of_log2_plus_one(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    imputed = prepare_pca_matrix(interaction, normalization="none")
    expected = np.log2(imputed + 1).apply(
        lambda row: (row - row.mean()) / row.std(), axis=1)
    matrix = prepare_pca_matrix(interaction, normalization="log2_zscore")
    pd.testing.assert_frame_equal(matrix, expected)


def test_zero_imputation_with_log2_stays_finite(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction, imputation="zero",
                                normalization="log2_zscore")
    assert np.isfinite(matrix.to_numpy()).all()


def test_min_detection_one_keeps_only_complete_preys(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction, min_detection_frac=1.0,
                                normalization="none")
    assert set(matrix.index) == {"P1", "P2", "P3"}


def test_min_detection_zero_keeps_every_prey(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction, min_detection_frac=0.0,
                                normalization="none")
    assert set(matrix.index) == set(PREY_VALUES)


def test_constant_prey_row_is_dropped_after_zscore(tmp_path):
    prey_values = dict(PREY_VALUES)
    prey_values["PC"] = [7.0, 7.0, 7.0, 7.0]
    interaction, _ = _write_dataset(tmp_path, prey_values=prey_values)
    matrix = prepare_pca_matrix(interaction)
    assert "PC" not in matrix.index
    assert not matrix.isna().any().any()


def test_too_few_preys_raises(tmp_path):
    prey_values = {"P1": PREY_VALUES["P1"], "P2": PREY_VALUES["P2"]}
    interaction, _ = _write_dataset(tmp_path, prey_values=prey_values)
    with pytest.raises(ValueError, match="prey"):
        prepare_pca_matrix(interaction)


def test_too_few_experiments_raises(tmp_path):
    experiments = [("a_1", "BaitA", "T")]
    prey_values = {p: v[:1] for p, v in PREY_VALUES.items()}
    interaction, _ = _write_dataset(tmp_path, prey_values=prey_values,
                                    experiments=experiments)
    with pytest.raises(ValueError, match="experiment"):
        prepare_pca_matrix(interaction)


# --- detection_counts ----------------------------------------------------------

def test_detection_counts_ignore_zeros_and_missing(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    counts = detection_counts(interaction)
    assert counts.to_dict() == {"P1": 4, "P2": 4, "P3": 4,
                                "P4": 2, "P5": 3, "P6": 1}


# --- reduce_categorical --------------------------------------------------------

def test_reduce_categorical_takes_first_semicolon_token():
    s = pd.Series(["Cytosol;Nucleoplasm", "Nucleoplasm"], index=["P1", "P2"])
    reduced = reduce_categorical(s)
    assert reduced.loc["P1"] == "Cytosol"
    assert reduced.loc["P2"] == "Nucleoplasm"


def test_reduce_categorical_collapses_rare_labels_to_other():
    values = ["A"] * 5 + ["B"] * 4 + ["C"] * 3 + ["D"]
    s = pd.Series(values, index=[f"P{i}" for i in range(len(values))])
    reduced = reduce_categorical(s, top_n=2)
    assert set(reduced.unique()) == {"A", "B", "Other"}
    assert (reduced == "Other").sum() == 4


def test_reduce_categorical_maps_missing_to_unknown():
    s = pd.Series(["Cytosol", np.nan], index=["P1", "P2"])
    reduced = reduce_categorical(s)
    assert reduced.loc["P2"] == "Unknown"


# --- load_pca_metadata ---------------------------------------------------------

def test_load_pca_metadata_merges_bait_and_type(tmp_path):
    interaction, ed = _write_dataset(tmp_path)
    metadata = load_pca_metadata(interaction, ed)
    row = metadata[metadata["Experiment"] == "b_2"].iloc[0]
    assert row["BaitName"] == "BaitB"
    assert row["Type"] == "C"
    assert len(metadata) == len(EXPERIMENTS)


# --- pca_plot ------------------------------------------------------------------

def _all_points(fig):
    xs, ys = [], []
    for trace in fig.data:
        xs.extend(trace.x)
        ys.extend(trace.y)
    return xs, ys


def test_pca_plot_has_one_point_per_experiment(tmp_path):
    interaction, ed = _write_dataset(tmp_path)
    fig = pca_plot(interaction, ed)
    xs, _ = _all_points(fig)
    assert len(xs) == len(EXPERIMENTS)


def test_pca_plot_accepts_a_precomputed_matrix(tmp_path):
    interaction, ed = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    from_path = pca_plot(interaction, ed)
    from_matrix = pca_plot(interaction, ed, matrix=matrix)
    assert sorted(_all_points(from_path)[0]) == pytest.approx(
        sorted(_all_points(from_matrix)[0]))


# --- prey_pca_plot -------------------------------------------------------------

def test_prey_pca_uncolored_has_one_point_per_prey(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    fig = prey_pca_plot(matrix)
    assert len(fig.data) == 1
    assert len(fig.data[0].x) == len(matrix)


def test_prey_pca_continuous_coloring(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    counts = detection_counts(interaction).reindex(matrix.index)
    fig = prey_pca_plot(matrix, color_values=counts,
                        color_label="Detections", color_mode="continuous")
    assert len(fig.data) == 1
    assert len(fig.data[0].marker.color) == len(matrix)


def test_prey_pca_categorical_coloring_shows_unknown(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    locations = pd.Series("Cytosol", index=matrix.index)
    locations.iloc[0] = np.nan
    fig = prey_pca_plot(matrix, color_values=reduce_categorical(locations),
                        color_label="Location", color_mode="categorical")
    names = {trace.name for trace in fig.data}
    assert "Unknown" in names
    assert "Cytosol" in names


# --- plot_exports --------------------------------------------------------------

def test_pca_plot_matplotlib_uses_shared_preprocessing(tmp_path, monkeypatch):
    interaction, ed = _write_dataset(tmp_path)
    calls = []
    real = plot_exports.prepare_pca_matrix

    def spy(*args, **kwargs):
        calls.append(args)
        return real(*args, **kwargs)

    monkeypatch.setattr(plot_exports, "prepare_pca_matrix", spy)
    fig = plot_exports.pca_plot_matplotlib(interaction, ed)
    assert isinstance(fig, matplotlib.figure.Figure)
    assert calls


def test_prey_pca_matplotlib_each_color_mode(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    counts = detection_counts(interaction).reindex(matrix.index)
    locations = reduce_categorical(pd.Series("Cytosol", index=matrix.index))
    for values, label, mode in [
        (None, None, "none"),
        (counts, "Detections", "continuous"),
        (locations, "Location", "categorical"),
    ]:
        fig = plot_exports.prey_pca_matplotlib(
            matrix, color_values=values, color_label=label, color_mode=mode)
        assert isinstance(fig, matplotlib.figure.Figure)


def test_prey_pca_continuous_threshold_greys_low_scores(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    scores = pd.Series([0.9, 0.05, 0.0, 0.5, 0.09], index=matrix.index)
    fig = prey_pca_plot(matrix, color_values=scores,
                        color_label="SaintScore", color_mode="continuous",
                        color_threshold=0.1)
    assert len(fig.data) == 2
    by_name = {trace.name: trace for trace in fig.data}
    below = by_name["SaintScore < 0.1"]
    assert len(below.x) == 3
    assert "grey" in below.marker.color or "gray" in below.marker.color
    colored = next(t for t in fig.data if t is not below)
    assert len(colored.x) == 2
    assert list(colored.marker.color) == [0.9, 0.5]


def test_prey_pca_matplotlib_accepts_color_threshold(tmp_path):
    interaction, _ = _write_dataset(tmp_path)
    matrix = prepare_pca_matrix(interaction)
    scores = pd.Series([0.9, 0.05, 0.0, 0.5, 0.09], index=matrix.index)
    fig = plot_exports.prey_pca_matplotlib(
        matrix, color_values=scores, color_label="SaintScore",
        color_mode="continuous", color_threshold=0.1)
    assert isinstance(fig, matplotlib.figure.Figure)
