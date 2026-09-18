"""Tests for the QC metric calculations behind the quality-control plots.

Covers the calculation half of GUI/QC_plots.py.  The plotly figures themselves are
presentation and are only touched where a figure is the sole way to observe a
calculation, as with the ROC curves' AUC.
"""

import os

import numpy as np
import pandas as pd
import pytest

import provenance
import QC_plots
from QC_plots import (
    apply_score_thresholds,
    calculate_network_degrees,
    calculate_threshold_metrics,
    roc_plot,
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


def test_empty_interactions_give_no_degrees(tmp_path):
    biogrid = _write_biogrid(tmp_path / "bg.csv", [("P1", "P2")])
    assert calculate_network_degrees(pd.DataFrame(columns=["First_ID", "Bait.ID"]),
                                     biogrid) == []


def test_preys_with_no_edges_get_zero(tmp_path):
    biogrid = _write_biogrid(tmp_path / "bg.csv", [("X1", "X2")])
    passing = pd.DataFrame({"First_ID": ["P1", "P2"], "Bait.ID": ["B1", "B1"]})

    assert calculate_network_degrees(passing, biogrid) == [0, 0]


def test_a_missing_reference_set_is_reported_as_absent_not_as_zeros(tmp_path):
    """A degree of zero is a real answer — a prey with no published partners among the
    other passing preys.  Absent data has to be distinguishable from it, or the metric
    reads as a measurement when no measurement was made."""
    passing = pd.DataFrame({"First_ID": ["P1", "P1", "P2"],
                            "Bait.ID": ["B1", "B2", "B1"]})

    assert calculate_network_degrees(passing, str(tmp_path / "absent.csv")) is None


def test_an_empty_reference_set_still_gives_real_zeros(tmp_path):
    """The distinction above only works if a readable database that happens to share no
    edges still counts as a measurement."""
    biogrid = _write_biogrid(tmp_path / "bg.csv", [("X1", "X2")])
    passing = pd.DataFrame({"First_ID": ["P1", "P2"], "Bait.ID": ["B1", "B1"]})

    assert calculate_network_degrees(passing, biogrid) == [0, 0]


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


def test_enrichment_ratio_is_zero_when_no_known_interactions_exist(tmp_path):
    """With known_before == 0 the ratio is undefined, but 0 is reported -- the same
    value that would mean 'passing interactions are never known'."""
    rows = [{**_score_row(), "Experiment.ID": "B1", "Prey.ID": f"P{i}",
             "First_ID": f"P{i}", "Bait.ID": "BID1", "In.BioGRID": False, "CCO": 0.5}
            for i in range(1, 4)]
    path = tmp_path / "annotated_scores.csv"
    pd.DataFrame(rows).to_csv(path, index=False)

    metrics = calculate_threshold_metrics(str(path), PASSING)

    assert metrics["known_before"] == 0
    assert metrics["enrichment_ratio"] == 0


def test_mean_degree_is_absent_when_the_reference_set_is_missing(scores_csv):
    """The GUI shows this as "no reference set" rather than as a number."""
    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    assert metrics["mean_degree"] is None


# --- locating the reference set ------------------------------------------------

def _biogrid_for(datasets_dir, organism, edges):
    """Write a BioGRID summary where setup_datasets puts one, for one organism."""
    organism_dir = datasets_dir / organism
    organism_dir.mkdir(parents=True, exist_ok=True)
    return _write_biogrid(organism_dir / "biogrid_summary.csv", edges)


@pytest.fixture
def datasets_dir(tmp_path, monkeypatch):
    """Stand in for the container's /Datasets tree."""
    root = tmp_path / "Datasets"
    root.mkdir()
    monkeypatch.setattr(provenance, "DEFAULT_DATASETS_DIR", str(root))
    return root


def test_the_reference_set_comes_from_the_organism_the_dataset_was_annotated_against(
        datasets_dir, scores_csv, monkeypatch):
    """The summary is built per organism.  Reading another organism's file would score
    a mouse experiment against human interactions."""
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260901T120000Z-0badcafe")
    with provenance.stage(os.path.dirname(scores_csv), "annotate") as record:
        record.extra(organism="mouse")
    _biogrid_for(datasets_dir, "mouse", [("P1", "P2"), ("P1", "P3")])
    _biogrid_for(datasets_dir, "human", [])

    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    assert metrics["mean_degree"] > 0


def test_the_top_of_the_datasets_directory_is_not_where_the_summary_is_sought(
        datasets_dir, scores_csv):
    """setup_datasets writes no summary there, so a build looking for one finds nothing
    and every dataset reports a degree of zero."""
    _write_biogrid(datasets_dir / "biogrid_summary.csv",
                   [("P1", "P2"), ("P1", "P3")])

    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    assert metrics["mean_degree"] is None


def test_a_dataset_with_no_manifest_falls_back_to_the_human_reference_set(
        datasets_dir, scores_csv):
    _biogrid_for(datasets_dir, "human", [("P1", "P2"), ("P1", "P3")])

    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    assert metrics["mean_degree"] > 0


def test_an_explicit_reference_set_overrides_the_recorded_organism(
        datasets_dir, scores_csv, tmp_path):
    """The command-line annotator takes an explicit path; the metric honors one too."""
    _biogrid_for(datasets_dir, "human", [])
    explicit = _write_biogrid(tmp_path / "elsewhere.csv", [("P1", "P2"), ("P1", "P3")])

    metrics = calculate_threshold_metrics(scores_csv, PASSING, biogrid_path=explicit)

    assert metrics["mean_degree"] > 0


# --- roc_plot ------------------------------------------------------------------

@pytest.fixture
def roc_scores_csv(tmp_path):
    """SaintScore ranks the known interactions perfectly; BFDR and WDFDR are its
    inverse, so a correct sign flip makes them rank perfectly too."""
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
    return str(path)


def test_fdr_scores_are_negated_so_lower_is_better(roc_scores_csv):
    """BFDR and WDFDR rank inversely to SaintScore.  Without the sign flip their AUC
    would be 0, not 1."""
    fig = roc_plot(roc_scores_csv, "BioGRID")
    aucs = {trace.name.split(" (AUC")[0]: float(trace.name.split("= ")[1].rstrip(")"))
            for trace in fig.data if "AUC" in (trace.name or "")}

    assert aucs["SaintScore"] == pytest.approx(1.0)
    assert aucs["BFDR"] == pytest.approx(1.0)
    assert aucs["WDFDR"] == pytest.approx(1.0)


def test_unknown_known_type_is_rejected(roc_scores_csv):
    with pytest.raises(ValueError, match="Unknown known_type"):
        roc_plot(roc_scores_csv, "NotAThing")


def test_nan_truth_values_count_as_not_known(tmp_path):
    """`x != x` maps NaN to False, so an unannotated prey is a negative rather than
    dropping out of the curve."""
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
    saint = next(t for t in fig.data if (t.name or "").startswith("SaintScore"))

    # 4 positives and 6 negatives, the NaN row among the negatives.
    assert float(saint.name.split("= ")[1].rstrip(")")) == pytest.approx(1.0)


def test_the_reference_set_follows_the_hcm_choice_the_dataset_was_annotated_with(
        datasets_dir, scores_csv, monkeypatch):
    """A run annotated without Human Cell Map evidence is measured against that same
    reduced summary, not the full one beside it."""
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260901T120000Z-0badcafe")
    with provenance.stage(os.path.dirname(scores_csv), "annotate") as record:
        record.extra(organism="human", exclude_hcm=True)
    _biogrid_for(datasets_dir, "human", [("P1", "P2"), ("P1", "P3")])
    _write_biogrid(datasets_dir / "human" / "biogrid_summary_no_hcm.csv", [("P8", "P9")])

    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    assert metrics["mean_degree"] == 0


# --- SAINT score vs fold change -------------------------------------------------

SCATTER_COLUMNS = {
    "Experiment.ID": "B1", "First_Prey_Gene": "GENE", "SaintScore": 0.9, "BFDR": 0.01,
    "FoldChange": 3.0, "In.BioGRID": True, "Multivalidated": False,
}


def _scatter_scores(tmp_path, **quant_columns):
    rows = [{**SCATTER_COLUMNS, **quant_columns},
            {**SCATTER_COLUMNS, "In.BioGRID": False, **quant_columns},
            {**SCATTER_COLUMNS, "Multivalidated": True, **quant_columns}]
    path = tmp_path / "annotated_scores.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _hover_texts(fig):
    return [t for trace in fig.data for t in trace.text]


def test_the_scatter_reads_intensity_columns_from_an_intensity_run(tmp_path):
    path = _scatter_scores(tmp_path, AvgIntensity=1.5e6, ctrlIntensity="1e5|.|2e5")

    texts = _hover_texts(QC_plots.saint_scatter_plot(path, "B1", 0.7))

    assert len(texts) == 3
    assert all("Avg Intensity: 1.50e+06" in t for t in texts)
    assert all("Avg Ctrl Intensity: 1.50e+05" in t for t in texts)


def test_the_scatter_reads_spectral_count_columns_from_a_spc_run(tmp_path):
    """SAINTexpress's spectral-count build names these columns AvgSpec and ctrlCounts."""
    path = _scatter_scores(tmp_path, AvgSpec=12.5, ctrlCounts="2|.|4")

    texts = _hover_texts(QC_plots.saint_scatter_plot(path, "B1", 0.7))

    assert len(texts) == 3
    assert all("Avg Spec: 12.5" in t for t in texts)
    assert all("Avg Ctrl Spec: 3.0" in t for t in texts)


def test_the_scatter_refuses_scores_with_neither_quantity_column(tmp_path):
    path = _scatter_scores(tmp_path)

    with pytest.raises(KeyError, match="AvgIntensity nor AvgSpec"):
        QC_plots.saint_scatter_plot(path, "B1", 0.7)


def test_baits_are_excluded_by_accession_when_the_column_is_present(tmp_path):
    """Symbol-keyed inputs carry the supplied symbol in Bait.ID and the resolved
    accession in Bait_Accession; BioGRID lists the accession."""
    biogrid = _write_biogrid(tmp_path / "biogrid.csv", [("P1", "BACC"), ("P1", "P2")])
    passing = pd.DataFrame({"First_ID": ["P1", "P2"], "Bait.ID": ["BSYM", "BSYM"],
                            "Bait_Accession": ["BACC", "BACC"]})

    assert sorted(calculate_network_degrees(passing, biogrid)) == [1, 1]
