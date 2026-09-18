"""Tests for building the BioGRID interaction summary.

``biogrid_summary.csv`` is what tells a user that an interaction their experiment found is
already published.  It is built once when the datasets are prepared, so an error here is
silent at build time and then shows up as a whole column of "not previously reported"
across every dataset scored afterwards.
"""

import numpy as np
import pandas as pd
import pytest

from preprocess_biogrid import preprocess_biogrid


HUMAN, MOUSE = 9606, 10090

ALL_COLUMNS = ["Organism ID Interactor A", "Organism ID Interactor B",
               "Experimental System Type", "Experimental System",
               "SWISS-PROT Accessions Interactor A",
               "SWISS-PROT Accessions Interactor B", "Author", "Publication Source"]


def _interaction(a, b, system="Affinity Capture-MS", system_type="physical",
                 organism_a=HUMAN, organism_b=HUMAN, author="Smith A (2020)",
                 source="PUBMED:1"):
    return {
        "Organism ID Interactor A": organism_a,
        "Organism ID Interactor B": organism_b,
        "Experimental System Type": system_type,
        "Experimental System": system,
        "SWISS-PROT Accessions Interactor A": a,
        "SWISS-PROT Accessions Interactor B": b,
        "Author": author,
        "Publication Source": source,
    }


def _multivalidated(a, b, organism_a=HUMAN, organism_b=HUMAN):
    return {
        "Organism ID Interactor A": organism_a,
        "Organism ID Interactor B": organism_b,
        "SWISS-PROT Accessions Interactor A": a,
        "SWISS-PROT Accessions Interactor B": b,
    }


@pytest.fixture
def build(tmp_path):
    """Write the two BioGRID exports, run the summariser, return the summary frame."""
    def _build(all_rows, mv_rows=(), organism_id=HUMAN):
        all_path = tmp_path / "BIOGRID-ALL.tab3.txt"
        mv_path = tmp_path / "BIOGRID-MV-Physical.tab3.txt"
        pd.DataFrame(list(all_rows), columns=ALL_COLUMNS).to_csv(
            all_path, sep="\t", index=False)
        pd.DataFrame(list(mv_rows), columns=[
            "Organism ID Interactor A", "Organism ID Interactor B",
            "SWISS-PROT Accessions Interactor A",
            "SWISS-PROT Accessions Interactor B"]).to_csv(
            mv_path, sep="\t", index=False)

        preprocess_biogrid(str(all_path), str(mv_path), str(tmp_path), organism_id)
        return pd.read_csv(tmp_path / "biogrid_summary.csv")

    return _build


def _pairs(summary):
    return set(zip(summary["SWISS-PROT Accessions Interactor A"],
                   summary["SWISS-PROT Accessions Interactor B"]))


# --- filtering -----------------------------------------------------------------

def test_a_physical_interaction_within_the_organism_is_kept(build):
    summary = build([_interaction("P1", "P2")])

    assert _pairs(summary) == {("P1", "P2")}


def test_genetic_interactions_are_dropped(build):
    """ProxiMate reports physical proximity; a genetic interaction is not evidence of it."""
    summary = build([_interaction("P1", "P2"),
                     _interaction("P3", "P4", system_type="genetic")])

    assert _pairs(summary) == {("P1", "P2")}


@pytest.mark.parametrize("organism_a, organism_b", [
    (MOUSE, HUMAN), (HUMAN, MOUSE), (MOUSE, MOUSE)])
def test_interactions_outside_the_organism_are_dropped(build, organism_a, organism_b):
    """Both partners must be in the organism; a cross-species record is not a
    within-proteome interaction."""
    summary = build([_interaction("P1", "P2"),
                     _interaction("P3", "P4", organism_a=organism_a,
                                  organism_b=organism_b)])

    assert _pairs(summary) == {("P1", "P2")}


def test_the_organism_is_selectable(build):
    summary = build([_interaction("P1", "P2"),
                     _interaction("M1", "M2", organism_a=MOUSE, organism_b=MOUSE)],
                    organism_id=MOUSE)

    assert _pairs(summary) == {("M1", "M2")}


# --- aggregation ----------------------------------------------------------------

def test_repeated_reports_of_a_pair_collapse_to_one_row(build):
    summary = build([
        _interaction("P1", "P2", author="Smith A (2020)", source="PUBMED:1"),
        _interaction("P1", "P2", author="Jones B (2021)", source="PUBMED:2"),
    ])

    assert len(summary) == 1


def test_the_supporting_evidence_is_joined_rather_than_replaced(build):
    """The count and variety of methods behind a published interaction is what a user
    weighs it by, so keeping only one report would overstate or understate it."""
    summary = build([
        _interaction("P1", "P2", system="Affinity Capture-MS", author="Smith A (2020)"),
        _interaction("P1", "P2", system="Two-hybrid", author="Jones B (2021)"),
    ])

    assert summary.loc[0, "Experimental System"] == "Affinity Capture-MS; Two-hybrid"
    assert summary.loc[0, "Author"] == "Smith A (2020); Jones B (2021)"


def test_repeated_evidence_is_not_deduplicated(build):
    """Two independent reports using the same method are two reports."""
    summary = build([_interaction("P1", "P2", system="Two-hybrid"),
                     _interaction("P1", "P2", system="Two-hybrid")])

    assert summary.loc[0, "Experimental System"] == "Two-hybrid; Two-hybrid"


def test_a_missing_author_is_joined_as_the_text_nan(build):
    """Documented, not fixed: the fields are cast to str before joining, so a record with
    no author contributes the literal "nan" to the list rather than dropping out of it.
    A pair whose only report lacks an author round-trips back to a null instead, because
    the lone "nan" is re-read as one."""
    both = build([_interaction("P1", "P2", author=np.nan),
                  _interaction("P1", "P2", author="Jones B (2021)")])
    alone = build([_interaction("P3", "P4", author=np.nan)])

    assert both.loc[0, "Author"] == "nan; Jones B (2021)"
    assert pd.isna(alone.loc[0, "Author"])


def test_the_two_orientations_of_a_pair_stay_separate(build):
    """Documented, not fixed: pairs are keyed as ordered, and annotation looks up only
    one orientation, so a published interaction recorded the other way round is reported
    as unseen."""
    summary = build([_interaction("P1", "P2"), _interaction("P2", "P1")])

    assert _pairs(summary) == {("P1", "P2"), ("P2", "P1")}


def test_a_pair_missing_an_accession_is_dropped(build):
    """Grouping drops null keys, so an interaction BioGRID could not map to SWISS-PROT
    does not reach the summary."""
    summary = build([_interaction("P1", "P2"), _interaction(np.nan, "P4")])

    assert _pairs(summary) == {("P1", "P2")}


# --- multivalidation --------------------------------------------------------------

def test_a_multivalidated_pair_is_flagged(build):
    summary = build([_interaction("P1", "P2")], [_multivalidated("P1", "P2")])

    assert summary.loc[0, "Multivalidated"] == True  # noqa: E712


def test_a_pair_absent_from_the_multivalidated_set_is_not_flagged(build):
    """The left join leaves a null here, and a plain truth test on a null is True.  Only
    the aggregation's null handling keeps the answer False, so an interaction supported
    once is not presented as independently confirmed."""
    summary = build([_interaction("P1", "P2")], [])

    assert summary.loc[0, "Multivalidated"] == False  # noqa: E712


def test_multivalidation_is_matched_per_pair(build):
    summary = build(
        [_interaction("P1", "P2"), _interaction("P3", "P4")],
        [_multivalidated("P1", "P2")]).set_index(
            "SWISS-PROT Accessions Interactor A")

    assert summary.loc["P1", "Multivalidated"] == True  # noqa: E712
    assert summary.loc["P3", "Multivalidated"] == False  # noqa: E712


def test_a_multivalidated_pair_outside_the_organism_does_not_flag(build):
    summary = build([_interaction("P1", "P2")],
                    [_multivalidated("P1", "P2", organism_a=MOUSE, organism_b=MOUSE)])

    assert summary.loc[0, "Multivalidated"] == False  # noqa: E712


def test_repeated_multivalidated_records_still_flag_once(build):
    summary = build([_interaction("P1", "P2")],
                    [_multivalidated("P1", "P2"), _multivalidated("P1", "P2")])

    assert len(summary) == 1
    assert summary.loc[0, "Multivalidated"] == True  # noqa: E712


# --- output schema -----------------------------------------------------------------

def test_every_summarised_pair_is_marked_as_known(build):
    """The column is the flag annotation joins on; every row here is by definition a
    published interaction."""
    summary = build([_interaction("P1", "P2"), _interaction("P3", "P4")])

    assert summary["In.BioGRID"].all()


def test_the_summary_carries_the_columns_annotation_reads(build):
    summary = build([_interaction("P1", "P2")])

    assert set(summary.columns) == {
        "SWISS-PROT Accessions Interactor A", "SWISS-PROT Accessions Interactor B",
        "Experimental System", "Author", "Publication Source",
        "Multivalidated", "In.BioGRID"}


def test_the_summary_is_written_where_annotation_looks_for_it(tmp_path, build):
    build([_interaction("P1", "P2")])

    assert (tmp_path / "biogrid_summary.csv").exists()


# --- excluding one publication's evidence ---------------------------------------

HCM = "PUBMED:34079125"


def _build_variant(tmp_path, all_rows, mv_rows=(), publication=HCM,
                   output_filename="biogrid_summary_no_hcm.csv"):
    """Run the summariser with one publication's evidence removed.

    The multivalidated export is written with a Publication Source column here, since
    the exclusion reads it from that file too.
    """
    all_path = tmp_path / "BIOGRID-ALL.tab3.txt"
    mv_path = tmp_path / "BIOGRID-MV-Physical.tab3.txt"
    pd.DataFrame(list(all_rows), columns=ALL_COLUMNS).to_csv(all_path, sep="\t", index=False)
    pd.DataFrame(list(mv_rows), columns=[
        "Organism ID Interactor A", "Organism ID Interactor B",
        "SWISS-PROT Accessions Interactor A", "SWISS-PROT Accessions Interactor B",
        "Publication Source"]).to_csv(mv_path, sep="\t", index=False)
    preprocess_biogrid(str(all_path), str(mv_path), str(tmp_path), HUMAN,
                       exclude_publication=publication, output_filename=output_filename)
    return pd.read_csv(tmp_path / output_filename)


def test_a_pair_reported_only_by_the_excluded_publication_is_dropped(tmp_path):
    summary = _build_variant(tmp_path, [_interaction("P1", "P2", source=HCM)])

    assert _pairs(summary) == set()


def test_a_pair_with_other_evidence_keeps_only_that_evidence(tmp_path):
    summary = _build_variant(tmp_path, [
        _interaction("P1", "P2", source=HCM, author="Go CD (2021)",
                     system="Proximity Label-MS"),
        _interaction("P1", "P2", source="PUBMED:1", author="Smith A (2020)",
                     system="Affinity Capture-MS"),
    ])

    assert _pairs(summary) == {("P1", "P2")}
    assert summary.loc[0, "Publication Source"] == "PUBMED:1"
    assert summary.loc[0, "Author"] == "Smith A (2020)"
    assert summary.loc[0, "Experimental System"] == "Affinity Capture-MS"


def test_multivalidation_reported_only_by_the_excluded_publication_is_dropped(tmp_path):
    all_rows = [_interaction("P1", "P2", source="PUBMED:1")]
    mv_rows = [dict(_multivalidated("P1", "P2"), **{"Publication Source": HCM})]

    kept = _build_variant(tmp_path, all_rows, mv_rows, publication=None)
    dropped = _build_variant(tmp_path, all_rows, mv_rows, publication=HCM)

    assert bool(kept.loc[0, "Multivalidated"]) is True
    assert bool(dropped.loc[0, "Multivalidated"]) is False


def test_the_variant_is_written_under_its_own_name(tmp_path, build):
    """Both summaries live in the same directory, so the variant must not overwrite
    the full one."""
    build([_interaction("P1", "P2", source=HCM)])
    variant = _build_variant(tmp_path, [_interaction("P1", "P2", source=HCM)])

    assert _pairs(variant) == set()
    assert _pairs(pd.read_csv(tmp_path / "biogrid_summary.csv")) == {("P1", "P2")}
