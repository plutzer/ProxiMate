"""Tests for the QC metric calculations behind the quality-control plots.

Covers the calculation half of GUI/QC_plots.py.  The plotly figures themselves are
presentation and are only touched where a figure is the sole way to observe a
calculation, as with the ROC curves' AUC.
"""

import numpy as np
import pandas as pd
import pytest

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


def test_missing_biogrid_file_returns_zeros_per_row_not_per_prey(tmp_path):
    """Documented inconsistency: the missing-file path returns one zero per interaction
    row, while the no-edges path returns one per unique prey.  mean_degree is therefore
    an average over a different denominator depending on which branch ran, and a
    missing database reads as a real degree of zero rather than as absent data."""
    passing = pd.DataFrame({"First_ID": ["P1", "P1", "P2"],
                            "Bait.ID": ["B1", "B2", "B1"]})

    degrees = calculate_network_degrees(passing, str(tmp_path / "absent.csv"))

    assert degrees == [0, 0, 0]
    assert len(degrees) == 3          # rows, though there are only 2 unique preys


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


def test_mean_degree_is_zero_without_the_packaged_biogrid_file(scores_csv):
    """calculate_threshold_metrics calls calculate_network_degrees without a path, so it
    always reads /Datasets/biogrid_summary.csv.  Outside the container that file is
    absent and mean_degree silently reports 0 instead of signalling missing data."""
    metrics = calculate_threshold_metrics(scores_csv, PASSING)

    assert metrics["mean_degree"] == 0


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
