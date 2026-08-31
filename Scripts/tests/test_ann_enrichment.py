"""Tests for the hypergeometric annotation-feature enrichment."""

import warnings

import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from Ann_Enrichment import enrich_foreground, process_refactored, split_and_clean


# --- split_and_clean -----------------------------------------------------------

def test_splits_on_semicolons_and_strips_whitespace():
    assert split_and_clean("Nucleus ; Cytoplasm;  Golgi ") == {
        "Nucleus", "Cytoplasm", "Golgi"}


@pytest.mark.parametrize("value", [np.nan, None, 3.5, ["Nucleus"]])
def test_non_string_input_yields_an_empty_set(value):
    assert split_and_clean(value) == set()


def test_empty_string_yields_an_empty_set():
    assert split_and_clean("") == set()


def test_trailing_digits_are_stripped():
    """Numbered instances of a feature collapse onto one label.

    UniProt numbers repeated features, so the Repeats, Domains and Motifs columns
    carry "WD 1", "WD 2", "WD 3"; enrichment asks whether a protein has WD repeats
    at all, not which one.  GO columns are unaffected either way -- annotator's
    clean_gocc has already stripped the bracketed accession by this point.
    """
    assert split_and_clean("WD 1;WD 2;WD 3") == {"WD"}
    assert split_and_clean("Complex1") == {"Complex"}


def test_annotations_beginning_with_a_digit_are_dropped():
    assert split_and_clean("40S Ribosome;Nucleus") == {"Nucleus"}


def test_all_digit_annotations_collapse_to_nothing():
    assert split_and_clean("12345;Nucleus") == {"Nucleus"}


def test_duplicates_collapse():
    assert split_and_clean("Nucleus;Nucleus;Nucleus") == {"Nucleus"}


# --- enrich_foreground ---------------------------------------------------------

@pytest.fixture
def feature_map():
    """Six of ten proteins carry "Nucleus"; two carry the rarer "Vesicle"."""
    mapping = {f"P{i:02d}": set() for i in range(1, 11)}
    for i in range(1, 7):
        mapping[f"P{i:02d}"] = {"Nucleus"}
    mapping["P07"] = {"Vesicle"}
    mapping["P08"] = {"Vesicle"}
    return mapping


def test_pvalue_and_enrichment_match_the_hypergeometric_definition(feature_map):
    all_ids = set(feature_map)                       # M = 10
    foreground = {"P01", "P02", "P03", "P07"}        # n = 4, Nucleus k = 3

    result = enrich_foreground(foreground, all_ids, feature_map).set_index("Feature")
    row = result.loc["Nucleus"]

    assert row["k"] == 3 and row["n"] == 4 and row["K"] == 6 and row["M"] == 10
    # Survival function at k-1 gives P(X >= k), the one-sided over-representation test.
    assert row["p_value"] == pytest.approx(hypergeom.sf(3 - 1, 10, 6, 4))
    assert row["enrichment"] == pytest.approx((3 / 4) / (6 / 10))


def test_rare_features_are_dropped(feature_map):
    """Vesicle has K = 2 < 5, so it never reaches the test even with k = 2."""
    result = enrich_foreground({"P01", "P02", "P07", "P08"}, set(feature_map), feature_map)

    assert "Nucleus" in set(result["Feature"])
    assert "Vesicle" not in set(result["Feature"])


def test_features_seen_once_in_the_foreground_are_dropped(feature_map):
    """Nucleus clears K >= 5 but a single foreground hit fails the k >= 2 guard."""
    result = enrich_foreground({"P01"}, set(feature_map), feature_map)

    assert "Nucleus" not in set(result["Feature"])


def test_results_are_sorted_by_ascending_pvalue():
    all_ids = {f"P{i:02d}" for i in range(1, 21)}
    mapping = {p: set() for p in all_ids}
    # "Common" is spread across the population; "Focused" concentrates in the foreground.
    for i in range(1, 16):
        mapping[f"P{i:02d}"].add("Common")
    for i in range(1, 7):
        mapping[f"P{i:02d}"].add("Focused")

    result = enrich_foreground({f"P{i:02d}" for i in range(1, 7)}, all_ids, mapping)

    assert len(result) == 2
    assert list(result["p_value"]) == sorted(result["p_value"])


def test_columns_are_stable_when_nothing_passes(feature_map):
    result = enrich_foreground(set(), set(feature_map), feature_map)

    assert len(result) == 0
    assert list(result.columns) == [
        "Feature", "k", "n", "K", "M", "p_value", "enrichment"]


def test_foreground_outside_all_ids_produces_k_greater_than_k_population():
    """feat_population counts only proteins in all_ids while feat_foreground counts every
    foreground protein, so callers must keep the foreground a subset of all_ids.
    process_refactored does; a direct caller need not."""
    all_ids = {f"P{i:02d}" for i in range(1, 6)}
    mapping = {f"P{i:02d}": {"Shared"} for i in range(1, 7)}   # P06 is outside all_ids
    foreground = set(mapping)                                   # includes P06

    row = enrich_foreground(foreground, all_ids, mapping).set_index("Feature").loc["Shared"]

    assert row["k"] == 6
    assert row["K"] == 5
    assert row["k"] > row["K"]


# --- process_refactored --------------------------------------------------------

@pytest.fixture
def annotated_scores():
    """One row per (Experiment.ID, Prey.ID), as annotated_scores.csv carries.

    Preys P01-P06 are nuclear, P07-P08 vesicular, P09-P12 cytoplasmic.  In E1 the
    high-scoring preys are P01-P03 plus P07; in E2 they are P01, P02 and P04-P06.
    """
    localization = {}
    for i in range(1, 7):
        localization[f"P{i:02d}"] = "Nucleus"
    for i in (7, 8):
        localization[f"P{i:02d}"] = "Vesicle"
    for i in range(9, 13):
        localization[f"P{i:02d}"] = "Cytoplasm"

    passing = {"E1": {"P01", "P02", "P03", "P07"},
               "E2": {"P01", "P02", "P04", "P05", "P06"}}

    rows = []
    for experiment in ("E1", "E2"):
        for prey, scl in localization.items():
            rows.append({
                "Experiment.ID": experiment,
                "Prey.ID": prey,
                "SaintScore": 0.9 if prey in passing[experiment] else 0.1,
                "SCL": scl,
            })
    return pd.DataFrame(rows)


def test_threshold_selects_the_foreground(annotated_scores):
    results = process_refactored(annotated_scores, ["SCL"], threshold=0.7)
    e1 = results[results["Bait"] == "E1"].set_index("Feature")

    # E1's foreground is four preys, three of them nuclear, out of twelve overall.
    assert e1.loc["Nucleus", "n"] == 4
    assert e1.loc["Nucleus", "k"] == 3
    assert e1.loc["Nucleus", "K"] == 6
    assert e1.loc["Nucleus", "M"] == 12


def test_a_lower_threshold_widens_the_foreground(annotated_scores):
    narrow = process_refactored(annotated_scores, ["SCL"], threshold=0.7)
    wide = process_refactored(annotated_scores, ["SCL"], threshold=0.0)

    assert narrow[narrow["Bait"] == "E1"]["n"].iloc[0] == 4
    assert wide[wide["Bait"] == "E1"]["n"].iloc[0] == 12


def test_adjusted_pvalues_are_corrected_within_each_bait_and_feature_type(annotated_scores):
    """BH runs per (feature type, experiment) group, not across the whole result."""
    results = process_refactored(annotated_scores, ["SCL"], threshold=0.7)

    for (_, _), group in results.groupby(["Feature_type", "Bait"]):
        expected = multipletests(group["p_value"].tolist(), method="fdr_bh")[1]
        assert np.allclose(group["adj_p"].to_numpy(dtype=float), expected)


def test_result_columns_are_in_a_stable_order(annotated_scores):
    """The GUI and its CSV export read these positionally-familiar columns."""
    results = process_refactored(annotated_scores, ["SCL"], threshold=0.7)

    assert list(results.columns) == [
        "Bait", "Feature", "Feature_type", "k", "n", "K", "M",
        "p_value", "enrichment", "adj_p"]


def test_empty_result_keeps_the_same_columns(annotated_scores):
    empty = process_refactored(annotated_scores, ["SCL"], threshold=1.5)
    populated = process_refactored(annotated_scores, ["SCL"], threshold=0.7)

    assert len(empty) == 0
    assert list(empty.columns) == list(populated.columns)


def test_count_columns_are_integers(annotated_scores):
    """k, n, K and M are counts.  Accumulating results onto an empty seed frame
    would leave them as object dtype."""
    results = process_refactored(annotated_scores, ["SCL"], threshold=0.7)

    for column in ("k", "n", "K", "M"):
        assert results[column].dtype == np.int64, column


def test_produces_no_pandas_warnings(annotated_scores):
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        warnings.simplefilter("error", DeprecationWarning)
        process_refactored(annotated_scores, ["SCL"], threshold=0.7)


def test_carries_bait_and_feature_type_labels(annotated_scores):
    results = process_refactored(annotated_scores, ["SCL"], threshold=0.7)

    assert set(results["Bait"]) == {"E1", "E2"}
    assert set(results["Feature_type"]) == {"SCL"}


def test_experiments_with_no_passing_features_contribute_nothing(annotated_scores):
    """A threshold above every score leaves an empty foreground, so no rows are added."""
    results = process_refactored(annotated_scores, ["SCL"], threshold=1.5)

    assert len(results) == 0


def test_feature_map_keeps_the_last_row_for_a_repeated_prey():
    """feature_map is built with dict(zip(...)) over every row, so when a prey's
    annotation differs between its experiment rows the last one silently wins."""
    rows = []
    for experiment, scl_for_p01 in (("E1", "Nucleus"), ("E2", "Cytoplasm")):
        for i in range(1, 8):
            prey = f"P{i:02d}"
            rows.append({
                "Experiment.ID": experiment,
                "Prey.ID": prey,
                "SaintScore": 0.9,
                "SCL": scl_for_p01 if prey == "P01" else "Nucleus",
            })
    data = pd.DataFrame(rows)

    results = process_refactored(data, ["SCL"], threshold=0.7)
    nucleus = results[results["Feature"] == "Nucleus"].iloc[0]

    # Seven preys, but P01's last row says Cytoplasm, so Nucleus counts only six.
    assert nucleus["M"] == 7
    assert nucleus["K"] == 6
