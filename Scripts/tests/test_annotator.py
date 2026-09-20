"""Tests for the annotation helpers.

Every string helper runs inside a ``.apply`` over the whole scored table, so its
behavior on the identifier and annotation formats the databases ship is what decides
whether a finished scoring run annotates cleanly.
"""

import numpy as np
import pandas as pd
import pytest

import annotator


SCL = "SUBCELLULAR LOCATION: "


@pytest.mark.parametrize("value, expected", [
    (SCL + "Nucleus; Cytoplasm", "Nucleus"),
    (SCL + "Nucleus {ECO:0000255|PROSITE}; Cytoplasm", "Nucleus"),
    (SCL + "[Isoform 2]: Nucleus", "Nucleus"),
    (SCL + "Cell membrane, single-pass type I membrane protein", "Cell membrane"),
    (SCL + "Cytoplasm.", "Cytoplasm"),
    (SCL + "Nucleus. " + SCL + "Mitochondrion", "Nucleus"),
], ids=["first-of-list", "evidence-code", "isoform-prefix", "qualifier", "period",
        "first-block-only"])
def test_get_first_scl_reads_the_first_location(value, expected):
    assert annotator.get_first_SCL(value) == expected


@pytest.mark.parametrize("value, expected", [
    ("nucleus [GO:0005634]; cytosol [GO:0005829]", "nucleus; cytosol"),
    ("nucleus ;cytosol", "nucleus; cytosol"),
    ("nucleus [GO:0005634];", "nucleus"),
    ("nucleus", "nucleus"),
], ids=["accessions-removed", "separators-normalized", "trailing-separator", "bare-term"])
def test_clean_gocc(value, expected):
    """Ann_Enrichment strips trailing digits from each term, which would truncate a
    GO accession left in place into a bare "GO:"."""
    assert annotator.clean_gocc(value) == expected


@pytest.mark.parametrize("value, expected", [
    ('MOTIF 10..20; /note="Nuclear localization signal"', "Nuclear localization signal"),
    ('/note="First"; /note="Second"', "First; Second"),
], ids=["single", "several"])
def test_clean_motif_extracts_note_names(value, expected):
    assert annotator.clean_motif(value) == expected


@pytest.mark.parametrize("prey_id, bait_id, expected", [
    ("P1", "P1", True),
    ("P1;P2;P3", "P2", True),
    ("P1;P2", "P9", False),
])
def test_self_inter_matches_any_member_of_the_prey_group(prey_id, bait_id, expected):
    assert annotator.self_inter(prey_id, bait_id) is expected


@pytest.mark.parametrize("prey_id, baits, expected", [
    ("P1", {"P1", "P2"}, True),
    ("P9;P2", {"P1", "P2"}, True),
    ("P9", {"P1", "P2"}, False),
    ("P1", set(), False),
])
def test_prey_is_bait_matches_any_member_of_the_prey_group(prey_id, baits, expected):
    assert annotator.prey_is_bait(prey_id, baits) is expected


@pytest.mark.parametrize("item, subcellular, uniprot, expected", [
    ("MATR3", ["MATR3", "AAA"], [], "MATR3"),
    ("OLDNAME", ["NEWNAME"], ["OLDNAME NEWNAME EXTRA"], "NEWNAME"),
    ("AAA", ["AAA"], [np.nan], "AAA"),
], ids=["known", "synonym-via-uniprot", "null-uniprot-row"])
def test_get_match_bridges_symbols_through_uniprot_synonyms(item, subcellular, uniprot,
                                                            expected):
    """HPA and the scored table do not always use the same symbol for a gene; the
    UniProt synonym list is what bridges them."""
    assert annotator.get_match(item, subcellular, uniprot) == expected


CC_DICT = {"BAIT1": {"PREY1": 0.85, "PREY2": 0.40}}


@pytest.mark.parametrize("bait, prey, cc_dict, expected", [
    ("BAIT1", "PREY1", CC_DICT, 0.85),
    ("BAIT9", "PREY1", CC_DICT, np.nan),
    ("BAIT1", "PREY9", CC_DICT, np.nan),
    ("BAIT1", "PREY9;PREY2", CC_DICT, 0.40),
    ("BAIT1", "PREY8;PREY9", CC_DICT, np.nan),
    (1, 2, {"1": {"2": 0.5}}, 0.5),
], ids=["known-pair", "unknown-bait", "unknown-prey", "first-scored-member",
        "no-scored-member", "ids-compared-as-text"])
def test_get_cco_score(bait, prey, cc_dict, expected):
    result = annotator.get_cco_score(bait, prey, cc_dict)
    if np.isnan(expected):
        assert pd.isnull(result)
    else:
        assert result == expected


COMPLEXES = {"Complex A": "P1;P2;P3", "Complex B": "P4;P5"}


@pytest.mark.parametrize("prey_id, expected", [
    ("P2", "Complex A"),
    ("P9;P4", "Complex B"),
    ("P9", None),
])
def test_complex_id_names_the_complex_holding_any_group_member(prey_id, expected):
    assert annotator.complex_id(prey_id, COMPLEXES) == expected


# --- gene symbol resolution ----------------------------------------------------

def _uniprot(rows):
    return pd.DataFrame(rows, columns=["Entry", "Gene Names"])


@pytest.mark.parametrize("value, is_accession", [
    ("P12345", True), ("A0A087X1C5", True), ("Q9Y6K9-2", True),
    ("RUVBL1", False), ("12345", False), ("P1234", False), ("contam_P12345", False),
])
def test_accession_pattern(value, is_accession):
    assert bool(annotator.ACCESSION_RE.match(value)) is is_accession


def test_symbol_map_covers_synonyms_and_leaves_out_shared_symbols():
    """Every symbol in an entry's list maps to it; a symbol listed under two entries
    maps to neither, since choosing one would annotate the wrong protein silently."""
    mapping = annotator.symbol_accession_map(_uniprot([
        ("Q9Y265", "RUVBL1 INO80H NMP238"), ("P1", "SHARED A"), ("P2", "SHARED B"),
        ("P3", None)]))

    assert mapping == {"RUVBL1": "Q9Y265", "INO80H": "Q9Y265", "NMP238": "Q9Y265",
                       "A": "P1", "B": "P2"}


def test_a_protein_group_resolves_element_wise_in_order():
    """Accessions pass through, symbols become accessions, unknown identifiers are kept
    as written."""
    resolved = annotator.resolve_accessions(
        "RUVBL1;P12345; RUVBL2;MYSTERY", {"RUVBL1": "Q9Y265", "RUVBL2": "Q9Y230"})

    assert resolved == "Q9Y265;P12345;Q9Y230;MYSTERY"


def test_unresolved_ids_lists_what_is_still_not_an_accession():
    column = pd.Series(["Q9Y265;MYSTERY", "P12345", "OTHER"])

    assert annotator.unresolved_ids(column) == ["MYSTERY", "OTHER"]


# --- Human Protein Atlas locations ----------------------------------------------

def _hpa(rows):
    return pd.DataFrame(rows, columns=["Gene name", "Main location"])


def test_conflicting_locations_are_joined_rather_than_dropped():
    """PINX1 ships with two different locations; picking one by row order would
    silently discard a real annotation."""
    collapsed = annotator.collapse_hpa_locations(
        _hpa([("PINX1", "Nucleoli"), ("PINX1", "Nuclear speckles")]))

    assert collapsed.loc[0, "Main location"] == "Nuclear speckles; Nucleoli"


def test_merging_collapsed_locations_cannot_duplicate_interactions():
    """``subcellular_location.tsv`` repeats some gene names.  A left merge on gene
    name must preserve the row count of the scored interactions, so the collapsed
    table holds one row per gene with identical repeats reduced to one value."""
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

    assert collapsed["Gene name"].is_unique
    assert len(merged) == len(scores)
    assert merged.set_index("Matched_Gene_Name").loc["MATR3", "Main location"] == "Nucleoplasm"
