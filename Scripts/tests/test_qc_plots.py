"""Tests for GUI/QC_plots.py and its matplotlib export twins in GUI/plot_exports.py.

The calculations are covered in depth; figures are only touched where a figure is the
sole way to observe a calculation (the ROC AUC, the scatter's grouping, the PCA point
count).  The data preparation behind each plotly/matplotlib pair is shared, so one data
test covers both renderers and the matplotlib twins get a smoke test.
"""

import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.figure
import numpy as np
import pandas as pd
import pytest

import plot_exports
import provenance
import QC_plots
from QC_plots import (
    apply_score_thresholds,
    calculate_network_degrees,
    calculate_threshold_metrics,
    detection_counts,
    known_status_split,
    load_pca_metadata,
    pca_plot,
    prepare_pca_matrix,
    prey_pca_plot,
    prey_gene_names,
    reduce_categorical,
    roc_plot,
    saint_known_retention,
    saint_scatter_plot,
)


@pytest.fixture(autouse=True)
def isolated_biogrid_cache(monkeypatch):
    """QC_plots caches BioGRID in module-level globals keyed by path.

    Without this, a fake BioGRID written by one test is served to every later test
    that asks for the same path.
    """
    monkeypatch.setattr(QC_plots, "_biogrid_cache", None, raising=False)
    monkeypatch.setattr(QC_plots, "_biogrid_cache_path", None, raising=False)


def _write_biogrid(path, edges):
    pd.DataFrame(
        [{"SWISS-PROT Accessions Interactor A": a,
          "SWISS-PROT Accessions Interactor B": b} for a, b in edges]
    ).to_csv(path, index=False)
    return str(path)


PASSING = {"SaintScore": 0.7, "BFDR": 0.05, "WD": 1.0, "WDFDR": 0.05}


def _score_row(**overrides):
    row = {"SaintScore": 0.9, "BFDR": 0.01, "WD": 5.0, "WDFDR": 0.01}
    row.update(overrides)
    return row


# --- apply_score_thresholds ----------------------------------------------------

def test_all_four_thresholds_must_pass():
    df = pd.DataFrame([
        _score_row(),                     # passes everything
        _score_row(SaintScore=0.5),       # fails SaintScore
        _score_row(BFDR=0.5),             # fails BFDR
        _score_row(WD=0.1),               # fails WD
        _score_row(WDFDR=0.5),            # fails WDFDR
    ])
    assert len(apply_score_thresholds(df, PASSING)) == 1


def test_thresholds_are_inclusive_at_the_boundary():
    df = pd.DataFrame([_score_row(SaintScore=0.7, BFDR=0.05, WD=1.0, WDFDR=0.05)])
    assert len(apply_score_thresholds(df, PASSING)) == 1


def test_missing_wdfdr_fails_the_threshold():
    """Scoring with 0 CompPASS iterations leaves WDFDR NaN; it is read as 1.0, so a
    row with no computed FDR is excluded rather than silently admitted."""
    df = pd.DataFrame([_score_row(WDFDR=np.nan)])

    assert len(apply_score_thresholds(df, PASSING)) == 0
    assert len(apply_score_thresholds(df, {**PASSING, "WDFDR": 1.0})) == 1


# --- calculate_network_degrees -------------------------------------------------

def test_counts_prey_prey_edges_only(tmp_path):
    """Degree counts edges between passing preys; anything touching a bait is dropped."""
    biogrid = _write_biogrid(tmp_path / "bg.csv", [
        ("P1", "P2"),      # prey-prey, counted
        ("P1", "P3"),      # prey-prey, counted
        ("P1", "B1"),      # touches a bait, excluded
        ("P2", "P9"),      # P9 is not a passing prey, excluded
    ])
    passing = pd.DataFrame({"First_ID": ["P1", "P2", "P3"],
                            "Bait.ID": ["B1", "B1", "B1"]})

    degrees = sorted(calculate_network_degrees(passing, biogrid))
    # P1 has two prey-prey partners; P2 and P3 have one each.
    assert degrees == [1, 1, 2]


def test_baits_are_excluded_by_accession_when_the_column_is_present(tmp_path):
    """Symbol-keyed inputs carry the supplied symbol in Bait.ID and the resolved
    accession in Bait_Accession; BioGRID lists the accession."""
    biogrid = _write_biogrid(tmp_path / "biogrid.csv", [("P1", "BACC"), ("P1", "P2")])
    passing = pd.DataFrame({"First_ID": ["P1", "P2"], "Bait.ID": ["BSYM", "BSYM"],
                            "Bait_Accession": ["BACC", "BACC"]})

    assert sorted(calculate_network_degrees(passing, biogrid)) == [1, 1]


def test_an_empty_reference_set_gives_real_zeros_but_a_missing_one_is_absent(tmp_path):
    """A degree of zero is a real answer: a prey with no published partners among the
    other passing preys.  Absent data has to be distinguishable from it, or the metric
    reads as a measurement when no measurement was made."""
    passing = pd.DataFrame({"First_ID": ["P1", "P2"], "Bait.ID": ["B1", "B1"]})

    biogrid = _write_biogrid(tmp_path / "bg.csv", [("X1", "X2")])
    assert calculate_network_degrees(passing, biogrid) == [0, 0]
    assert calculate_network_degrees(passing, str(tmp_path / "absent.csv")) is None


# --- calculate_threshold_metrics -----------------------------------------------

@pytest.fixture
def scores_csv(tmp_path):
    """Ten interactions across two baits; five pass, and known interactions are
    enriched among the passing set."""
    rows = []
    for i in range(1, 6):                      # passing, 4 of 5 known
        rows.append({**_score_row(), "Experiment.ID": "B1" if i <= 3 else "B2",
                     "Prey.ID": f"P{i}", "First_ID": f"P{i}", "Bait.ID": "BID1",
                     "In.BioGRID": i <= 4, "CCO": 0.5})
    for i in range(6, 11):                     # failing, 1 of 5 known
        rows.append({**_score_row(SaintScore=0.1), "Experiment.ID": "B1",
                     "Prey.ID": f"P{i}", "First_ID": f"P{i}", "Bait.ID": "BID1",
                     "In.BioGRID": i == 6, "CCO": 0.5})
    path = tmp_path / "annotated_scores.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def test_counts_before_and_after_filtering(scores_csv):
    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    assert metrics["total_before"] == 10
    assert metrics["total_after"] == 5
    assert metrics["known_before"] == 5
    assert metrics["known_after"] == 4


def test_enrichment_ratio_is_the_ratio_of_known_fractions(scores_csv):
    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    # (4/5) / (5/10)
    assert metrics["enrichment_ratio"] == pytest.approx((4 / 5) / (5 / 10))


def test_median_network_size_is_per_bait(scores_csv):
    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    # B1 keeps three passing preys, B2 keeps two; median of [3, 2].
    assert metrics["median_network_size"] == pytest.approx(2.5)


def test_control_filter_restricts_the_population(scores_csv):
    metrics = calculate_threshold_metrics(scores_csv, PASSING, ctrl_experiments=["B2"])

    assert metrics["total_before"] == 2
    assert metrics["total_after"] == 2


def test_metrics_are_zeroed_when_nothing_passes(scores_csv):
    impossible = {"SaintScore": 1.1, "BFDR": 0.0, "WD": 100.0, "WDFDR": 0.0}
    metrics = calculate_threshold_metrics(scores_csv, impossible)

    assert metrics["total_after"] == 0
    assert metrics["median_network_size"] == 0
    assert metrics["enrichment_ratio"] == 0
    assert metrics["mean_degree"] == 0


def test_mean_degree_is_absent_when_the_reference_set_is_missing(scores_csv):
    """The GUI shows this as "no reference set" rather than as a number."""
    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    assert metrics["mean_degree"] is None


# --- locating the reference set ------------------------------------------------

def _biogrid_for(datasets_dir, organism, edges, filename="biogrid_summary.csv"):
    """Write a BioGRID summary where setup_datasets puts one, for one organism."""
    organism_dir = datasets_dir / organism
    organism_dir.mkdir(parents=True, exist_ok=True)
    return _write_biogrid(organism_dir / filename, edges)


@pytest.fixture
def datasets_dir(tmp_path, monkeypatch):
    """Stand in for the container's /Datasets tree."""
    root = tmp_path / "Datasets"
    root.mkdir()
    monkeypatch.setattr(provenance, "DEFAULT_DATASETS_DIR", str(root))
    return root


def _annotated_with(scores_csv, monkeypatch, **settings):
    """Record an annotate stage in the dataset's manifest with the given settings."""
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260901T120000Z-0badcafe")
    with provenance.stage(os.path.dirname(scores_csv), "annotate") as record:
        record.extra(**settings)


CONNECTED = [("P1", "P2"), ("P1", "P3")]


@pytest.mark.parametrize("case", ["organism", "no_manifest", "explicit", "hcm"])
def test_the_reference_set_follows_how_the_dataset_was_annotated(
        case, datasets_dir, scores_csv, monkeypatch, tmp_path):
    """The summary is built per organism and per Human Cell Map choice.  Reading
    another organism's file would score a mouse experiment against human
    interactions; an explicit path (the command-line annotator takes one) wins over
    the manifest; a dataset with no manifest falls back to human."""
    kwargs = {}
    if case == "organism":
        _annotated_with(scores_csv, monkeypatch, organism="mouse")
        _biogrid_for(datasets_dir, "mouse", CONNECTED)
        _biogrid_for(datasets_dir, "human", [])
    elif case == "no_manifest":
        _biogrid_for(datasets_dir, "human", CONNECTED)
    elif case == "explicit":
        _biogrid_for(datasets_dir, "human", [])
        kwargs["biogrid_path"] = _write_biogrid(tmp_path / "elsewhere.csv", CONNECTED)
    elif case == "hcm":
        _annotated_with(scores_csv, monkeypatch, organism="human", exclude_hcm=True)
        _biogrid_for(datasets_dir, "human", [])
        _biogrid_for(datasets_dir, "human", CONNECTED, "biogrid_summary_no_hcm.csv")

    metrics = calculate_threshold_metrics(scores_csv, PASSING, **kwargs)

    assert metrics["mean_degree"] > 0


# --- roc_plot ------------------------------------------------------------------

def _auc(fig, score):
    trace = next(t for t in fig.data if (t.name or "").startswith(score))
    return float(trace.name.split("= ")[1].rstrip(")"))


def test_fdr_scores_are_negated_so_lower_is_better(tmp_path):
    """SaintScore ranks the known interactions perfectly; BFDR and WDFDR are its
    inverse.  Without the sign flip their AUC would be 0, not 1."""
    rows = []
    for i in range(10):
        known = i < 5
        rows.append({
            "Experiment.ID": "B1",
            "SaintScore": 0.9 if known else 0.1,
            "BFDR": 0.01 if known else 0.9,
            "WD": 5.0 if known else 0.5,
            "WDFDR": 0.01 if known else 0.9,
            "In.BioGRID": known,
            "Multivalidated": known,
        })
    path = tmp_path / "annotated_scores.csv"
    pd.DataFrame(rows).to_csv(path, index=False)

    fig = roc_plot(str(path), "BioGRID")

    for score in ("SaintScore", "BFDR", "WDFDR"):
        assert _auc(fig, score) == pytest.approx(1.0)
    for trace in fig.data:
        if "AUC" in (trace.name or ""):
            assert (trace.x[0], trace.y[0]) == (0.0, 0.0)
            assert (trace.x[-1], trace.y[-1]) == (1.0, 1.0)


def test_nan_truth_values_count_as_not_known(tmp_path):
    """An unannotated prey is a negative rather than dropping out of the curve."""
    rows = []
    for i in range(10):
        known = i < 4
        rows.append({
            "Experiment.ID": "B1",
            "SaintScore": 0.9 if known else 0.1,
            "BFDR": 0.01, "WD": 5.0, "WDFDR": 0.01,
            "In.BioGRID": known if i != 9 else np.nan,
            "Multivalidated": known,
        })
    path = tmp_path / "annotated_scores.csv"
    pd.DataFrame(rows).to_csv(path, index=False)

    fig = roc_plot(str(path), "BioGRID")

    # 4 positives and 6 negatives, the NaN row among the negatives.
    assert _auc(fig, "SaintScore") == pytest.approx(1.0)


# --- saint_known_retention -----------------------------------------------------

def test_known_retention_is_the_known_fraction_above_each_threshold(tmp_path):
    """Known fraction and mean CCO among the rows at or above each SaintScore step."""
    rows = [
        {"Experiment.ID": "B1", "SaintScore": 0.95, "In.BioGRID": True, "CCO": 1.0},
        {"Experiment.ID": "B1", "SaintScore": 0.85, "In.BioGRID": False, "CCO": 0.5},
        {"Experiment.ID": "B1", "SaintScore": 0.30, "In.BioGRID": False, "CCO": 0.0},
        {"Experiment.ID": "B1", "SaintScore": 0.30, "In.BioGRID": True, "CCO": 0.5},
    ]
    path = tmp_path / "annotated_scores.csv"
    pd.DataFrame(rows).to_csv(path, index=False)

    fig = saint_known_retention(str(path))
    retention = dict(zip(np.round(fig.data[0].x, 2), fig.data[0].y))
    cco = dict(zip(np.round(fig.data[1].x, 2), fig.data[1].y))

    assert retention[0.0] == pytest.approx(0.5) and cco[0.0] == pytest.approx(0.5)
    assert retention[0.5] == pytest.approx(0.5) and cco[0.5] == pytest.approx(0.75)
    assert retention[0.9] == pytest.approx(1.0) and cco[0.9] == pytest.approx(1.0)
    # Nothing at or above 1.0: the fraction is reported as 0, not NaN.
    assert retention[1.0] == 0


# --- SAINT score vs fold change -------------------------------------------------

SCATTER_COLUMNS = {
    "Experiment.ID": "B1", "First_Prey_Gene": "GENE", "SaintScore": 0.9, "BFDR": 0.01,
    "FoldChange": 3.0, "In.BioGRID": True, "Multivalidated": False,
}


def _scatter_scores(tmp_path, **quant_columns):
    """Three rows: one in BioGRID only, one not in BioGRID, one multivalidated."""
    rows = [{**SCATTER_COLUMNS, **quant_columns},
            {**SCATTER_COLUMNS, "In.BioGRID": False, **quant_columns},
            {**SCATTER_COLUMNS, "Multivalidated": True, **quant_columns}]
    path = tmp_path / "annotated_scores.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


@pytest.mark.parametrize("columns, avg, ctrl", [
    ({"AvgIntensity": 1.5e6, "ctrlIntensity": "1e5|.|2e5"},
     "Avg Intensity: 1.50e+06", "Avg Ctrl Intensity: 1.50e+05"),
    # SAINTexpress's spectral-count build names these columns AvgSpec and ctrlCounts.
    ({"AvgSpec": 12.5, "ctrlCounts": "2|.|4"}, "Avg Spec: 12.5", "Avg Ctrl Spec: 3.0"),
])
def test_the_scatter_reads_the_quantity_columns_of_the_run_type(tmp_path, columns, avg, ctrl):
    path = _scatter_scores(tmp_path, **columns)

    texts = [t for trace in saint_scatter_plot(path, "B1", 0.7).data for t in trace.text]

    assert len(texts) == 3
    assert all(avg in t for t in texts)
    assert all(ctrl in t for t in texts)


def test_scatter_traces_group_rows_by_known_status(tmp_path):
    """A multivalidated row is drawn once, in the Multivalidated trace, not also in
    the In BioGRID one; the same split feeds the matplotlib export."""
    path = _scatter_scores(tmp_path, AvgIntensity=1.5e6, ctrlIntensity="1e5")
    rows = pd.read_csv(path)
    rows["First_Prey_Gene"] = ["BG", "NONE", "MV"]
    rows.to_csv(path, index=False)

    fig = saint_scatter_plot(path, "B1", 0.7)
    genes = {trace.name: [t.split("</b>")[0].lstrip("<b>") for t in trace.text]
             for trace in fig.data}
    assert genes == {"Not in BioGRID": ["NONE"], "In BioGRID": ["BG"],
                     "Multivalidated": ["MV"]}

    groups = known_status_split(rows)
    assert [list(g["First_Prey_Gene"]) for g in groups] == [["NONE"], ["BG"], ["MV"]]

    without_mv = known_status_split(rows.drop(columns=["Multivalidated"]))
    assert [len(g) for g in without_mv] == [1, 2, 0]
    without_known = known_status_split(rows.drop(columns=["Multivalidated", "In.BioGRID"]))
    assert [len(g) for g in without_known] == [3, 0, 0]


def test_the_matplotlib_scatter_draws_one_series_per_group(tmp_path):
    path = _scatter_scores(tmp_path, AvgIntensity=1.5e6, ctrlIntensity="1e5")

    fig = plot_exports.saint_scatter_matplotlib(path, "B1", 0.7)

    labels = [c.get_label() for c in fig.axes[0].collections]
    assert labels == ["Not in BioGRID", "In BioGRID", "Multivalidated"]


# --- PCA data preparation ------------------------------------------------------

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
    tmp_path.mkdir(exist_ok=True)
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


@pytest.fixture
def dataset(tmp_path):
    return _write_dataset(tmp_path)


def test_default_detection_filter_drops_sparse_preys(dataset):
    interaction, _ = dataset
    assert set(prepare_pca_matrix(interaction).index) == {"P1", "P2", "P3", "P4", "P5"}


def test_zero_imputation_fills_missing_with_zero(dataset):
    interaction, _ = dataset
    matrix = prepare_pca_matrix(interaction, imputation="zero", normalization="none")
    assert matrix.loc["P4", "b_1"] == 0.0
    assert matrix.loc["P5", "a_1"] == 0.0
    assert not matrix.isna().any().any()


def test_drop_imputation_keeps_only_complete_preys(dataset):
    interaction, _ = dataset
    matrix = prepare_pca_matrix(interaction, imputation="drop", normalization="none")
    assert set(matrix.index) == {"P1", "P2", "P3"}
    assert not matrix.isna().any().any()


def test_row_min_imputation_fills_with_row_minimum(dataset):
    interaction, _ = dataset
    matrix = prepare_pca_matrix(interaction, normalization="none")
    assert matrix.loc["P4", "b_1"] == 30.0
    assert matrix.loc["P5", "a_1"] == 40.0


def test_log2_zscore_is_zscore_of_log2_plus_one(dataset):
    interaction, _ = dataset
    imputed = prepare_pca_matrix(interaction, normalization="none")
    expected = np.log2(imputed + 1).apply(
        lambda row: (row - row.mean()) / row.std(), axis=1)
    matrix = prepare_pca_matrix(interaction, normalization="log2_zscore")
    pd.testing.assert_frame_equal(matrix, expected)


def test_zero_imputation_with_log2_stays_finite(dataset):
    interaction, _ = dataset
    matrix = prepare_pca_matrix(interaction, imputation="zero", normalization="log2_zscore")
    assert np.isfinite(matrix.to_numpy()).all()


@pytest.mark.parametrize("fraction, expected", [
    (1.0, {"P1", "P2", "P3"}), (0.0, set(PREY_VALUES))])
def test_min_detection_fraction_selects_preys(dataset, fraction, expected):
    interaction, _ = dataset
    matrix = prepare_pca_matrix(interaction, min_detection_frac=fraction,
                                normalization="none")
    assert set(matrix.index) == expected


def test_constant_prey_row_is_dropped_after_zscore(tmp_path):
    prey_values = dict(PREY_VALUES)
    prey_values["PC"] = [7.0, 7.0, 7.0, 7.0]
    interaction, _ = _write_dataset(tmp_path, prey_values=prey_values)
    matrix = prepare_pca_matrix(interaction)
    assert "PC" not in matrix.index
    assert not matrix.isna().any().any()


def test_too_few_preys_or_experiments_raise(tmp_path):
    interaction, _ = _write_dataset(
        tmp_path / "preys", {"P1": PREY_VALUES["P1"], "P2": PREY_VALUES["P2"]})
    with pytest.raises(ValueError):
        prepare_pca_matrix(interaction)

    interaction, _ = _write_dataset(
        tmp_path / "experiments", {p: v[:1] for p, v in PREY_VALUES.items()},
        experiments=[("a_1", "BaitA", "T")])
    with pytest.raises(ValueError):
        prepare_pca_matrix(interaction)


def test_a_never_detected_prey_is_dropped_before_pca_without_normalization(tmp_path):
    """Zeros are non-detections, so an all-zero prey has nothing to impute from.  It
    must not reach PCA as a row of NaN under the one normalization that does not
    otherwise remove NaN rows."""
    counts = {"P1": [12, 10, 30, 28], "P2": [3, 4, 5, 6], "P3": [8, 9, 11, 10],
              "P4": [0, 0, 0, 0], "P5": [0, 2, 3, 4]}
    interaction, _ = _write_dataset(tmp_path, counts)

    data = prepare_pca_matrix(interaction, min_detection_frac=0.0,
                              imputation="row_min", normalization="none")

    assert set(data.index) == {"P1", "P2", "P3", "P5"}
    assert np.isfinite(data.values).all()


def test_detection_counts_ignore_zeros_and_missing(dataset):
    interaction, _ = dataset
    assert detection_counts(interaction).to_dict() == {
        "P1": 4, "P2": 4, "P3": 4, "P4": 2, "P5": 3, "P6": 1}


def test_reduce_categorical_takes_the_first_token_collapses_rare_labels_and_fills_missing():
    values = ["A;B"] * 5 + ["B"] * 4 + ["C"] * 3 + ["D", np.nan]
    s = pd.Series(values, index=[f"P{i}" for i in range(len(values))])

    reduced = reduce_categorical(s, top_n=2)

    assert set(reduced.unique()) == {"A", "B", "Other", "Unknown"}
    assert reduced.loc["P0"] == "A"
    assert (reduced == "Other").sum() == 4
    assert reduced.loc["P13"] == "Unknown"


def test_load_pca_metadata_merges_bait_and_type(dataset):
    interaction, ed = dataset
    metadata = load_pca_metadata(interaction, ed)
    row = metadata[metadata["Experiment"] == "b_2"].iloc[0]
    assert row["BaitName"] == "BaitB"
    assert row["Type"] == "C"
    assert len(metadata) == len(EXPERIMENTS)


# --- PCA figures ---------------------------------------------------------------

def test_pca_plot_has_one_point_per_experiment(dataset):
    interaction, ed = dataset
    fig = pca_plot(interaction, ed)
    assert sum(len(trace.x) for trace in fig.data) == len(EXPERIMENTS)
    assert isinstance(plot_exports.pca_plot_matplotlib(interaction, ed),
                      matplotlib.figure.Figure)


@pytest.mark.parametrize("mode", ["none", "continuous", "categorical"])
def test_prey_pca_color_modes(dataset, mode):
    interaction, _ = dataset
    matrix = prepare_pca_matrix(interaction)
    values, label = None, None
    if mode == "continuous":
        values, label = detection_counts(interaction).reindex(matrix.index), "Detections"
    elif mode == "categorical":
        locations = pd.Series("Cytosol", index=matrix.index)
        locations.iloc[0] = np.nan
        values, label = reduce_categorical(locations), "Location"

    fig = prey_pca_plot(matrix, color_values=values, color_label=label, color_mode=mode)

    assert sum(len(trace.x) for trace in fig.data) == len(matrix)
    if mode == "continuous":
        assert len(fig.data[0].marker.color) == len(matrix)
    if mode == "categorical":
        assert {trace.name for trace in fig.data} == {"Cytosol", "Unknown"}
    assert isinstance(plot_exports.prey_pca_matplotlib(
        matrix, color_values=values, color_label=label, color_mode=mode),
        matplotlib.figure.Figure)


def test_prey_pca_continuous_threshold_greys_low_scores(dataset):
    interaction, _ = dataset
    matrix = prepare_pca_matrix(interaction)
    scores = pd.Series([0.9, 0.05, 0.0, 0.5, 0.09], index=matrix.index)

    fig = prey_pca_plot(matrix, color_values=scores, color_label="SaintScore",
                        color_mode="continuous", color_threshold=0.1)

    assert len(fig.data) == 2
    below = next(t for t in fig.data if t.name == "SaintScore < 0.1")
    assert len(below.x) == 3
    assert "grey" in below.marker.color or "gray" in below.marker.color
    colored = next(t for t in fig.data if t is not below)
    assert list(colored.marker.color) == [0.9, 0.5]
    assert isinstance(plot_exports.prey_pca_matplotlib(
        matrix, color_values=scores, color_label="SaintScore",
        color_mode="continuous", color_threshold=0.1), matplotlib.figure.Figure)


def test_prey_pca_hover_shows_gene_name_above_accession(dataset, tmp_path):
    interaction, _ = dataset
    matrix = prepare_pca_matrix(interaction)
    prey_file = tmp_path / "prey.txt"
    prey_file.write_text("P1\t100\tGENE1\nP2\t200\tGENE2\n")
    names = prey_gene_names(prey_file)
    assert names.to_dict() == {"P1": "GENE1", "P2": "GENE2"}

    scores = pd.Series([0.9, 0.05, 0.0, 0.5, 0.09], index=matrix.index)
    fig = prey_pca_plot(matrix, color_values=scores, color_label="SaintScore",
                        color_mode="continuous", color_threshold=0.1,
                        gene_names=names)

    for trace in fig.data:
        hover = trace.hovertemplate
        assert hover.startswith("<b>%{")
        assert "</b><br>%{customdata[0]}<br>" in hover
    colored = next(t for t in fig.data if t.name != "SaintScore < 0.1")
    # P1 has a gene name; P4 does not and keeps its accession as the label
    assert list(colored.hovertext) == ["GENE1", "P4"]
    assert [row[0] for row in colored.customdata] == ["P1", "P4"]
