"""Tests for the annotation string and lookup helpers.

Every one of these runs inside a ``.apply`` over the whole scored table, so a value one of
them cannot handle stops the annotation of a finished scoring run rather than producing a
blank cell.  Several have no null guard at all, and those cases are pinned here.

The HPA location collapsing and the coverage check are covered by test_annotator_hpa.py.
"""

import numpy as np
import pandas as pd
import pytest

import annotator


# --- get_first_SCL -------------------------------------------------------------

SCL = "SUBCELLULAR LOCATION: "


@pytest.mark.parametrize("value", [None, np.nan, pd.NA])
def test_a_null_location_yields_a_null(value):
    assert pd.isnull(annotator.get_first_SCL(value))


def test_the_first_location_is_taken():
    assert annotator.get_first_SCL(SCL + "Nucleus; Cytoplasm") == "Nucleus"


def test_an_evidence_code_is_stripped():
    assert annotator.get_first_SCL(
        SCL + "Nucleus {ECO:0000255|PROSITE}; Cytoplasm") == "Nucleus"


def test_an_isoform_prefix_is_stripped():
    assert annotator.get_first_SCL(SCL + "[Isoform 2]: Nucleus") == "Nucleus"


def test_a_qualifier_after_a_comma_is_dropped():
    assert annotator.get_first_SCL(
        SCL + "Cell membrane, single-pass type I membrane protein") == "Cell membrane"


def test_a_trailing_period_is_dropped():
    assert annotator.get_first_SCL(SCL + "Cytoplasm.") == "Cytoplasm"


def test_only_the_first_marked_block_is_read():
    """A cell listing several SUBCELLULAR LOCATION blocks keeps the first and discards
    the rest without reporting it."""
    assert annotator.get_first_SCL(
        SCL + "Nucleus. " + SCL + "Mitochondrion") == "Nucleus"


@pytest.mark.parametrize("value", ["", "Nucleus", "no marker here"])
def test_a_cell_without_the_marker_raises(value):
    """Documented, not fixed: an empty string is not null, so the guard misses it and the
    split has no second element.  This runs inside a .apply with no try, so one such cell
    stops the annotation of an already-scored run."""
    with pytest.raises(IndexError):
        annotator.get_first_SCL(value)


# --- clean_gocc and trim_GO_CC -------------------------------------------------

def test_bracketed_accessions_are_removed():
    """Ann_Enrichment depends on this: its split_and_clean strips trailing digits, which
    would otherwise truncate a GO accession into a bare "GO:"."""
    assert annotator.clean_gocc(
        "nucleus [GO:0005634]; cytosol [GO:0005829]") == "nucleus; cytosol"


def test_separators_are_normalized():
    assert annotator.clean_gocc("nucleus ;cytosol") == "nucleus; cytosol"


def test_a_trailing_separator_is_removed():
    assert annotator.clean_gocc("nucleus [GO:0005634];") == "nucleus"


def test_an_unannotated_term_is_left_alone():
    assert annotator.clean_gocc("nucleus") == "nucleus"


def test_duplicate_terms_are_not_collapsed():
    """Deduplication happens downstream, in the set that split_and_clean builds."""
    assert annotator.clean_gocc("nucleus; nucleus") == "nucleus; nucleus"


@pytest.mark.parametrize("value", [None, np.nan])
def test_a_null_go_annotation_yields_a_null(value):
    assert pd.isnull(annotator.trim_GO_CC(value))


def test_the_go_wrapper_delegates_to_the_cleaner():
    assert annotator.trim_GO_CC("nucleus [GO:0005634]") == "nucleus"


def test_the_go_cleaner_has_no_null_guard_of_its_own():
    with pytest.raises(TypeError):
        annotator.clean_gocc(None)


# --- clean_motif and trim_motifs -----------------------------------------------

def test_motif_names_are_extracted_from_their_notes():
    assert annotator.clean_motif('MOTIF 10..20; /note="Nuclear localization signal"') == \
        "Nuclear localization signal"


def test_several_motifs_are_joined():
    assert annotator.clean_motif('/note="First"; /note="Second"') == "First; Second"


def test_a_cell_with_no_note_yields_an_empty_string():
    """Documented, not fixed: the result is "" rather than a null, so a downstream
    isnull check treats an unannotated protein as annotated with nothing."""
    assert annotator.clean_motif("MOTIF 10..20") == ""


@pytest.mark.parametrize("value", [None, np.nan])
def test_a_null_motif_yields_a_null(value):
    assert pd.isnull(annotator.trim_motifs(value))


# --- self_inter ----------------------------------------------------------------

def test_a_prey_matching_the_bait_is_a_self_interaction():
    assert annotator.self_inter("P1", "P1") is True


def test_any_member_of_a_prey_group_can_match():
    assert annotator.self_inter("P1;P2;P3", "P2") is True


def test_an_unmatched_prey_is_not_a_self_interaction():
    assert annotator.self_inter("P1;P2", "P9") is False


def test_group_members_are_not_trimmed():
    """Documented, not fixed: the split is on ";" with no strip, so a group written with
    spaces after the separator does not match."""
    assert annotator.self_inter("P1; P2", "P2") is False


def test_a_null_prey_raises():
    with pytest.raises(AttributeError):
        annotator.self_inter(np.nan, "P1")


# --- prey_is_bait ---------------------------------------------------------------

def test_a_prey_that_is_also_a_bait_is_flagged():
    assert annotator.prey_is_bait("P1", {"P1", "P2"}) is True


def test_any_member_of_a_prey_group_can_be_a_bait():
    assert annotator.prey_is_bait("P9;P2", {"P1", "P2"}) is True


def test_a_prey_that_is_no_bait_is_not_flagged():
    assert annotator.prey_is_bait("P9", {"P1", "P2"}) is False


def test_an_empty_bait_set_flags_nothing():
    assert annotator.prey_is_bait("P1", set()) is False


# --- get_first_pg ---------------------------------------------------------------

def test_the_first_gene_of_a_group_is_taken():
    assert annotator.get_first_pg("G1;G2;G3") == "G1"


def test_a_single_gene_is_returned_unchanged():
    assert annotator.get_first_pg("G1") == "G1"


def test_a_null_gene_raises():
    """Documented, not fixed: there is no guard, so one protein with no gene name stops
    the annotation run."""
    with pytest.raises(AttributeError):
        annotator.get_first_pg(np.nan)


# --- get_match ------------------------------------------------------------------

def test_a_gene_already_known_is_returned_unchanged():
    assert annotator.get_match("MATR3", ["MATR3", "AAA"], []) == "MATR3"


def test_a_synonym_is_resolved_through_uniprot():
    """HPA and the scored table do not always use the same symbol for a gene; the UniProt
    synonym list is what bridges them."""
    resolved = annotator.get_match("OLDNAME", ["NEWNAME"], ["OLDNAME NEWNAME EXTRA"])

    assert resolved == "NEWNAME"


def test_an_unresolvable_gene_yields_nothing():
    """Documented, not fixed: the function falls off the end and returns None, which
    becomes a merge key and yields empty locations rather than an error."""
    assert annotator.get_match("UNKNOWN", ["AAA"], ["BBB CCC"]) is None


def test_a_uniprot_row_with_no_gene_names_is_survivable():
    """The cell is coerced with str(), so a null does not stop the scan."""
    assert annotator.get_match("AAA", ["AAA"], [np.nan]) == "AAA"


# --- get_cco_score ---------------------------------------------------------------

CC_DICT = {"BAIT1": {"PREY1": 0.85, "PREY2": 0.40}}


def test_a_known_pair_returns_its_score():
    assert annotator.get_cco_score("BAIT1", "PREY1", CC_DICT) == 0.85


def test_an_unknown_bait_yields_a_null():
    assert pd.isnull(annotator.get_cco_score("BAIT9", "PREY1", CC_DICT))


def test_an_unknown_prey_yields_a_null():
    assert pd.isnull(annotator.get_cco_score("BAIT1", "PREY9", CC_DICT))


def test_the_first_scored_member_of_a_prey_group_wins():
    assert annotator.get_cco_score("BAIT1", "PREY9;PREY2", CC_DICT) == 0.40


def test_a_prey_group_with_no_scored_member_yields_a_null():
    assert pd.isnull(annotator.get_cco_score("BAIT1", "PREY8;PREY9", CC_DICT))


def test_an_empty_dictionary_yields_a_null():
    assert pd.isnull(annotator.get_cco_score("BAIT1", "PREY1", {}))


def test_identifiers_are_compared_as_text():
    """The GOGO output is parsed into string keys whatever the identifiers looked like
    in the scored table."""
    assert annotator.get_cco_score(1, 2, {"1": {"2": 0.5}}) == 0.5


# --- complex_id -------------------------------------------------------------------

COMPLEXES = {"Complex A": "P1;P2;P3", "Complex B": "P4;P5"}


def test_a_prey_in_a_complex_is_named():
    assert annotator.complex_id("P2", COMPLEXES) == "Complex A"


def test_any_member_of_a_prey_group_can_place_it():
    assert annotator.complex_id("P9;P4", COMPLEXES) == "Complex B"


def test_a_prey_in_no_complex_yields_nothing():
    assert annotator.complex_id("P9", COMPLEXES) is None


def test_a_prey_in_several_complexes_takes_the_first_listed():
    """Documented, not fixed: only one complex is reported, chosen by dictionary order
    rather than by evidence, so a prey shared between complexes loses the others."""
    shared = {"Complex A": "P1;P2", "Complex B": "P2;P3"}

    assert annotator.complex_id("P2", shared) == "Complex A"


def test_an_empty_complex_dictionary_yields_nothing():
    assert annotator.complex_id("P1", {}) is None


# --- organism configuration ---------------------------------------------------------

def test_the_supported_organisms_are_configured():
    assert set(annotator.ORGANISMS) == {"human", "mouse", "yeast"}


@pytest.mark.parametrize("organism, taxonomy_id", [
    ("human", 9606), ("mouse", 10090), ("yeast", 559292)])
def test_each_organism_carries_its_taxonomy_id(organism, taxonomy_id):
    """The same id filters BioGRID and selects the UniProt download, so a wrong one
    yields an empty annotation rather than an error."""
    assert annotator.ORGANISMS[organism]["organism_id"] == taxonomy_id


def test_only_human_has_the_human_specific_databases():
    """HPA and CORUM are human-only, and the annotate step skips both for any organism
    whose flags say so."""
    for organism, config in annotator.ORGANISMS.items():
        expected = organism == "human"
        assert config["has_hpa"] is expected
        assert config["has_corum"] is expected
