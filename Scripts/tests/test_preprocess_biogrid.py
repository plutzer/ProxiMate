"""Tests for building the BioGRID interaction summary.

``biogrid_summary.csv`` is what tells a user that an interaction their experiment found is
already published.  It is built once when the datasets are prepared, so an error here is
silent at build time and then shows up as a whole column of "not previously reported"
across every dataset scored afterwards.
"""

import pandas as pd
import pytest

from preprocess_biogrid import preprocess_biogrid


HUMAN, MOUSE = 9606, 10090

# Publication whose evidence the "no HCM" summary variant leaves out.
HCM = "PUBMED:34079125"

ALL_COLUMNS = ["Organism ID Interactor A", "Organism ID Interactor B",
               "Experimental System Type", "Experimental System",
               "SWISS-PROT Accessions Interactor A",
               "SWISS-PROT Accessions Interactor B", "Author", "Publication Source"]

MV_COLUMNS = ["Organism ID Interactor A", "Organism ID Interactor B",
              "SWISS-PROT Accessions Interactor A",
              "SWISS-PROT Accessions Interactor B", "Publication Source"]


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


def _multivalidated(a, b, organism_a=HUMAN, organism_b=HUMAN, source="PUBMED:1"):
    return {
        "Organism ID Interactor A": organism_a,
        "Organism ID Interactor B": organism_b,
        "SWISS-PROT Accessions Interactor A": a,
        "SWISS-PROT Accessions Interactor B": b,
        "Publication Source": source,
    }


@pytest.fixture
def build(tmp_path):
    """Write the two BioGRID exports, run the summariser, return the summary frame."""
    def _build(all_rows, mv_rows=(), organism_id=HUMAN, exclude_publication=None,
               output_filename="biogrid_summary.csv"):
        all_path = tmp_path / "BIOGRID-ALL.tab3.txt"
        mv_path = tmp_path / "BIOGRID-MV-Physical.tab3.txt"
        pd.DataFrame(list(all_rows), columns=ALL_COLUMNS).to_csv(
            all_path, sep="\t", index=False)
        pd.DataFrame(list(mv_rows), columns=MV_COLUMNS).to_csv(
            mv_path, sep="\t", index=False)

        preprocess_biogrid(str(all_path), str(mv_path), str(tmp_path), organism_id,
                           exclude_publication=exclude_publication,
                           output_filename=output_filename)
        return pd.read_csv(tmp_path / output_filename)

    return _build


def _pairs(summary):
    return set(zip(summary["SWISS-PROT Accessions Interactor A"],
                   summary["SWISS-PROT Accessions Interactor B"]))


# --- filtering -----------------------------------------------------------------

def test_genetic_interactions_are_dropped(build):
    """ProxiMate reports physical proximity; a genetic interaction is not evidence of it."""
    summary = build([_interaction("P1", "P2"),
                     _interaction("P3", "P4", system_type="genetic")])

    assert _pairs(summary) == {("P1", "P2")}


@pytest.mark.parametrize("organism_id, expected", [
    (HUMAN, {("P1", "P2")}), (MOUSE, {("M1", "M2")})])
def test_only_pairs_with_both_partners_in_the_organism_are_kept(build, organism_id,
                                                                expected):
    """A cross-species record is not a within-proteome interaction."""
    summary = build([
        _interaction("P1", "P2"),
        _interaction("M1", "M2", organism_a=MOUSE, organism_b=MOUSE),
        _interaction("X1", "X2", organism_a=MOUSE, organism_b=HUMAN),
        _interaction("X3", "X4", organism_a=HUMAN, organism_b=MOUSE),
    ], organism_id=organism_id)

    assert _pairs(summary) == expected


# --- aggregation ----------------------------------------------------------------

def test_the_supporting_evidence_is_joined_into_one_row_per_pair(build):
    """The count and variety of methods behind a published interaction is what a user
    weighs it by, so keeping only one report would overstate or understate it."""
    summary = build([
        _interaction("P1", "P2", system="Affinity Capture-MS", author="Smith A (2020)",
                     source="PUBMED:1"),
        _interaction("P1", "P2", system="Two-hybrid", author="Jones B (2021)",
                     source="PUBMED:2"),
    ])

    assert len(summary) == 1
    assert summary.loc[0, "Experimental System"] == "Affinity Capture-MS; Two-hybrid"
    assert summary.loc[0, "Author"] == "Smith A (2020); Jones B (2021)"
    assert summary.loc[0, "Publication Source"] == "PUBMED:1; PUBMED:2"


# --- multivalidation --------------------------------------------------------------

def test_multivalidation_is_matched_per_pair(build):
    """The left join leaves a null for a pair absent from the multivalidated set, and a
    plain truth test on a null is True; the aggregation must still report False so an
    interaction supported once is not presented as independently confirmed."""
    summary = build(
        [_interaction("P1", "P2"), _interaction("P3", "P4")],
        [_multivalidated("P1", "P2")]).set_index(
            "SWISS-PROT Accessions Interactor A")

    assert bool(summary.loc["P1", "Multivalidated"]) is True
    assert bool(summary.loc["P3", "Multivalidated"]) is False


def test_a_multivalidated_pair_outside_the_organism_does_not_flag(build):
    summary = build([_interaction("P1", "P2")],
                    [_multivalidated("P1", "P2", organism_a=MOUSE, organism_b=MOUSE)])

    assert bool(summary.loc[0, "Multivalidated"]) is False


def test_repeated_multivalidated_records_still_flag_once(build):
    summary = build([_interaction("P1", "P2")],
                    [_multivalidated("P1", "P2"), _multivalidated("P1", "P2")])

    assert len(summary) == 1
    assert bool(summary.loc[0, "Multivalidated"]) is True


# --- output schema -----------------------------------------------------------------

def test_the_summary_carries_the_columns_annotation_reads(build):
    """In.BioGRID is the flag annotation joins on; every row here is by definition a
    published interaction."""
    summary = build([_interaction("P1", "P2"), _interaction("P3", "P4")])

    assert set(summary.columns) == {
        "SWISS-PROT Accessions Interactor A", "SWISS-PROT Accessions Interactor B",
        "Experimental System", "Author", "Publication Source",
        "Multivalidated", "In.BioGRID"}
    assert summary["In.BioGRID"].all()


# --- excluding one publication's evidence ---------------------------------------

def test_excluding_a_publication_keeps_only_the_other_evidence(build):
    """A pair reported only by the excluded publication drops out; one with other
    evidence keeps exactly that evidence."""
    summary = build([
        _interaction("P1", "P2", source=HCM, author="Go CD (2021)",
                     system="Proximity Label-MS"),
        _interaction("P1", "P2", source="PUBMED:1", author="Smith A (2020)",
                     system="Affinity Capture-MS"),
        _interaction("P3", "P4", source=HCM),
    ], exclude_publication=HCM, output_filename="biogrid_summary_no_hcm.csv")

    assert _pairs(summary) == {("P1", "P2")}
    assert summary.loc[0, "Publication Source"] == "PUBMED:1"
    assert summary.loc[0, "Author"] == "Smith A (2020)"
    assert summary.loc[0, "Experimental System"] == "Affinity Capture-MS"


def test_multivalidation_reported_only_by_the_excluded_publication_is_dropped(build):
    all_rows = [_interaction("P1", "P2", source="PUBMED:1")]
    mv_rows = [_multivalidated("P1", "P2", source=HCM)]

    kept = build(all_rows, mv_rows)
    dropped = build(all_rows, mv_rows, exclude_publication=HCM)

    assert bool(kept.loc[0, "Multivalidated"]) is True
    assert bool(dropped.loc[0, "Multivalidated"]) is False
