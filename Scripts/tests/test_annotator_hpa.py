"""Tests for collapsing Human Protein Atlas subcellular locations.

The distributed ``subcellular_location.tsv`` repeats some gene names.  Merging
it onto the scored interactions on gene name without collapsing multiplies every
interaction row whose prey is one of those genes, which silently inflates
network sizes and enrichment counts downstream.
"""

import pandas as pd
import pytest

import annotator


def _hpa(rows):
    return pd.DataFrame(rows, columns=["Gene name", "Main location"])


def test_one_row_per_gene_is_returned():
    collapsed = annotator.collapse_hpa_locations(
        _hpa([("MATR3", "Nucleoplasm"), ("MATR3", "Nucleoplasm"), ("AAA", "Cytosol")]))

    assert collapsed["Gene name"].is_unique
    assert len(collapsed) == 2


def test_identical_duplicate_rows_collapse_to_the_single_value():
    collapsed = annotator.collapse_hpa_locations(
        _hpa([("MATR3", "Nucleoplasm"), ("MATR3", "Nucleoplasm")]))

    assert collapsed.loc[0, "Main location"] == "Nucleoplasm"


def test_conflicting_locations_are_joined_rather_than_dropped():
    """PINX1 ships with two different locations; picking one by row order would
    silently discard a real annotation."""
    collapsed = annotator.collapse_hpa_locations(
        _hpa([("PINX1", "Nucleoli"), ("PINX1", "Nuclear speckles")]))

    assert collapsed.loc[0, "Main location"] == "Nuclear speckles; Nucleoli"


def test_conflicting_locations_are_reported(caplog):
    """A future HPA release introducing a new conflict must announce itself."""
    with caplog.at_level("WARNING", logger="proximate.annotator"):
        annotator.collapse_hpa_locations(
            _hpa([("PINX1", "Nucleoli"), ("PINX1", "Nuclear speckles")]))

    assert "PINX1" in caplog.text


def test_identical_duplicates_are_not_reported_as_conflicts(caplog):
    with caplog.at_level("WARNING", logger="proximate.annotator"):
        annotator.collapse_hpa_locations(
            _hpa([("MATR3", "Nucleoplasm"), ("MATR3", "Nucleoplasm")]))

    assert "MATR3" not in caplog.text


def test_a_single_location_is_unchanged():
    collapsed = annotator.collapse_hpa_locations(_hpa([("AAA", "Cytosol")]))

    assert collapsed.loc[0, "Main location"] == "Cytosol"


def test_a_gene_with_no_location_stays_empty():
    collapsed = annotator.collapse_hpa_locations(_hpa([("AAA", None)]))

    assert pd.isna(collapsed.loc[0, "Main location"])


def test_an_annotation_source_matching_nothing_is_reported(caplog):
    """Zero matches across a whole dataset is far more often an identifier
    mismatch than a real absence of annotation, and reads as a finding."""
    with caplog.at_level("WARNING", logger="proximate.annotator"):
        annotator.check_annotation_coverage(0, 4500, "BioGRID")

    assert "BioGRID" in caplog.text


def test_partial_annotation_coverage_is_not_reported(caplog):
    with caplog.at_level("WARNING", logger="proximate.annotator"):
        annotator.check_annotation_coverage(136, 4500, "BioGRID")

    assert caplog.text == ""


def test_an_empty_dataset_is_not_reported(caplog):
    """Nothing to match is not evidence of a broken join."""
    with caplog.at_level("WARNING", logger="proximate.annotator"):
        annotator.check_annotation_coverage(0, 0, "BioGRID")

    assert caplog.text == ""


def test_merging_collapsed_locations_cannot_duplicate_interactions():
    """The property that actually matters: a left merge on gene name must
    preserve the row count of the scored interactions."""
    scores = pd.DataFrame({
        "Matched_Gene_Name": ["MATR3", "PINX1", "AAA"],
        "SaintScore": [0.1, 0.2, 0.3],
    })
    collapsed = annotator.collapse_hpa_locations(_hpa([
        ("MATR3", "Nucleoplasm"), ("MATR3", "Nucleoplasm"),
        ("PINX1", "Nucleoli"), ("PINX1", "Nuclear speckles"),
        ("AAA", "Cytosol"),
    ]))

    merged = scores.merge(collapsed, left_on="Matched_Gene_Name",
                          right_on="Gene name", how="left")

    assert len(merged) == len(scores)
