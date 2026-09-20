"""Tests for the MaxQuant proteinGroups funnel.

Every input format reaches SAINT and CompPASS through ``ProteinGroups``: the MaxQuant
parser hands it the file directly, and the other parsers convert their vendor export
into a MaxQuant-shaped frame first.  Its filtering, its control-only removal and its
protein-ID truncation therefore decide what every run scores.
"""

import pandas as pd
import pytest

from experimental_design import ExperimentalDesign
from ed_exceptions import EDPGMismatchError, ProxiMateError
from protein_groups import ProteinGroups, get_quant_col_prefix


PG_METADATA = ["Majority protein IDs", "Gene names", "Reverse", "Decoy",
               "Only identified by site", "Potential contaminant", "Sequence length"]

DEFAULT_ED = [
    {"Experiment Name": "t1_1", "Type": "T", "Bait": "BaitA", "Replicate": 1,
     "Bait ID": "P_A"},
    {"Experiment Name": "t1_2", "Type": "T", "Bait": "BaitA", "Replicate": 2,
     "Bait ID": "P_A"},
    {"Experiment Name": "c_1", "Type": "C", "Bait": "Ctrl", "Replicate": 1,
     "Bait ID": "P_C"},
]

INTERACTION_COLUMNS = ["Experiment", "Bait", "Prey", "Value"]


def _protein(ids, genes, intensities, reverse="-", site="-", contaminant="-", length=100):
    """One proteinGroups row; `intensities` maps experiment name to quant value."""
    row = {
        "Majority protein IDs": ids,
        "Gene names": genes,
        "Reverse": reverse,
        "Only identified by site": site,
        "Potential contaminant": contaminant,
        "Sequence length": length,
    }
    row.update(intensities)
    return row


@pytest.fixture
def build(tmp_path):
    """Construct a ProteinGroups from row dicts, writing the files it reads.

    Keys of each row that are not proteinGroups metadata are experiment names; the
    quantification prefix is applied here so tests name experiments rather than columns.
    """
    counter = {"n": 0}

    def _build(protein_rows, ed_rows=None, quant="Intensity"):
        counter["n"] += 1
        suffix = counter["n"]
        ed_rows = DEFAULT_ED if ed_rows is None else ed_rows
        prefix = get_quant_col_prefix(quant)

        ed_path = tmp_path / "ED_{}.csv".format(suffix)
        pd.DataFrame(ed_rows).to_csv(ed_path, index=False)

        renamed = []
        for row in protein_rows:
            out = {k: v for k, v in row.items() if k in PG_METADATA}
            for key, value in row.items():
                if key not in PG_METADATA:
                    out["{}{}".format(prefix, key)] = value
            renamed.append(out)

        pg_path = tmp_path / "proteinGroups_{}.txt".format(suffix)
        pd.DataFrame(renamed).to_csv(pg_path, sep="\t", index=False)

        return ProteinGroups(ExperimentalDesign(str(ed_path)), str(pg_path), quant)

    return _build


@pytest.fixture
def out_dir(tmp_path):
    target = tmp_path / "out"
    target.mkdir()
    return target


# Survivors: P1, P6 (truncated), P7, P8 (gene backfilled).
STANDARD_ROWS = [
    _protein("P1", "G1", {"t1_1": 10, "t1_2": 20, "c_1": 0}),
    _protein("P2rev", "G2", {"t1_1": 10, "t1_2": 10, "c_1": 0}, reverse="+"),
    _protein("P3site", "G3", {"t1_1": 10, "t1_2": 10, "c_1": 0}, site="+"),
    _protein("P4con", "G4", {"t1_1": 10, "t1_2": 10, "c_1": 0}, contaminant="+"),
    _protein("P5ctrl", "G5", {"t1_1": 0, "t1_2": 0, "c_1": 50}),
    _protein("P6a;P6b;P6c;P6d", "G6a;G6b;G6c;G6d", {"t1_1": 5, "t1_2": 5, "c_1": 0}),
    _protein("P7a;P7b;P7c", "G7a;G7b;G7c", {"t1_1": 5, "t1_2": 0, "c_1": 0}),
    _protein("P8", None, {"t1_1": 7, "t1_2": 0, "c_1": 0}),
]

SURVIVING_PROTEINS = 4


@pytest.fixture
def standard(build):
    return build(STANDARD_ROWS)


def _interaction(out_dir):
    return pd.read_csv(out_dir / "interaction.txt", sep="\t", header=None,
                       names=INTERACTION_COLUMNS)


def _comppass(out_dir):
    return pd.read_csv(out_dir / "to_CompPASS.csv")


# --- quantification --------------------------------------------------------------

@pytest.mark.parametrize("quant, prefix", [
    ("LFQ", "LFQ intensity "),
    ("Intensity", "Intensity "),
    ("Spectral Counts", "MS/MS count "),
])
def test_known_quantifications_map_to_a_prefix(quant, prefix):
    assert get_quant_col_prefix(quant) == prefix


@pytest.mark.parametrize("quant", ["spc", "intensity", None])
def test_an_unknown_quantification_is_rejected(quant):
    with pytest.raises(ValueError):
        get_quant_col_prefix(quant)


# --- column discovery ------------------------------------------------------------

def test_quant_columns_are_discovered_for_every_experiment(standard):
    assert standard.quant_cols == ["Intensity t1_1", "Intensity t1_2", "Intensity c_1"]


def test_bait_columns_exclude_controls(standard):
    """The control-only filter sums across these, so a control counted here would keep
    proteins seen in no bait at all."""
    assert standard.bait_cols == ["Intensity t1_1", "Intensity t1_2"]


def test_a_data_column_absent_from_the_design_is_skipped(build):
    rows = [_protein("P1", "G1", {"t1_1": 10, "t1_2": 10, "c_1": 0, "ghost": 5})]

    assert "Intensity ghost" not in build(rows).quant_cols


def test_a_design_experiment_absent_from_the_data_is_rejected(build):
    """Scoring an experiment with no quantification would report zeros as measurements."""
    rows = [_protein("P1", "G1", {"t1_1": 10, "t1_2": 10})]

    with pytest.raises(EDPGMismatchError) as excinfo:
        build(rows)

    assert excinfo.value.ed_only == ["c_1"]


def test_a_design_with_no_test_rows_is_rejected(build):
    """With no bait columns every protein would count as control-only and be dropped."""
    controls_only = [
        {"Experiment Name": "c_1", "Type": "C", "Bait": "Ctrl", "Replicate": 1,
         "Bait ID": "P_C"},
    ]
    rows = [_protein("P1", "G1", {"c_1": 100})]

    with pytest.raises(ProxiMateError):
        build(rows, ed_rows=controls_only)


# --- filtering -------------------------------------------------------------------

def test_flagged_proteins_are_removed(standard):
    kept = set(standard.data["Majority protein IDs"])

    assert {"P2rev", "P3site", "P4con"} & kept == set()


def test_only_an_exact_plus_is_treated_as_a_flag(build):
    """MaxQuant writes "+" or leaves the cell empty; a NaN must not drop the row."""
    rows = [
        _protein("Pkeep", "G", {"t1_1": 5, "t1_2": 5, "c_1": 0}, reverse=None),
        _protein("Pdash", "G", {"t1_1": 5, "t1_2": 5, "c_1": 0}, reverse="-"),
        _protein("Pdrop", "G", {"t1_1": 5, "t1_2": 5, "c_1": 0}, reverse="+"),
    ]

    assert set(build(rows).data["Majority protein IDs"]) == {"Pkeep", "Pdash"}


def test_a_decoy_column_is_read_as_reverse(build):
    """MaxQuant 2.4 and later write the decoy flag under "Decoy"."""
    rows = [
        {**_protein("Pkeep", "G", {"t1_1": 5, "t1_2": 5, "c_1": 0}), "Decoy": "-"},
        {**_protein("Pdrop", "G", {"t1_1": 5, "t1_2": 5, "c_1": 0}), "Decoy": "+"},
    ]
    for row in rows:
        del row["Reverse"]

    pg = build(rows)

    assert set(pg.data["Majority protein IDs"]) == {"Pkeep"}
    assert "Reverse" in pg.data.columns and "Decoy" not in pg.data.columns


def test_proteins_absent_from_every_bait_are_removed(standard):
    """P5ctrl is seen only in the control."""
    assert "P5ctrl" not in set(standard.data["Majority protein IDs"])


def test_a_protein_seen_in_a_single_bait_replicate_is_kept(standard):
    """The criterion is the sum across bait columns, not presence in all of them."""
    assert "P8" in set(standard.data["Short protein IDs"])


# --- identifier truncation ---------------------------------------------------------

@pytest.mark.parametrize("ids, expected", [
    ("A;B;C", "A;B;C"),
    ("A;B;C;D", "A;B;C;+1"),
    ("A;B;C;D;E", "A;B;C;+2"),
])
def test_more_than_three_identifiers_are_truncated_with_a_dropped_count(build, ids,
                                                                        expected):
    rows = [_protein(ids, "G", {"t1_1": 5, "t1_2": 5, "c_1": 0})]

    assert build(rows).data["Short protein IDs"].iloc[0] == expected


def test_genes_are_truncated_independently_of_identifiers(build):
    rows = [_protein("A;B;C;D", "G1;G2", {"t1_1": 5, "t1_2": 5, "c_1": 0})]

    pg = build(rows)

    assert pg.data["Short protein IDs"].iloc[0] == "A;B;C;+1"
    assert pg.data["Short Gene names"].iloc[0] == "G1;G2"


def test_spaces_in_identifiers_become_hyphens(build):
    """SAINT reads whitespace-delimited files, so a space inside an identifier would
    split the row -- whether or not the identifier was also truncated."""
    rows = [
        _protein("A 1;B;C;D", "G", {"t1_1": 5, "t1_2": 5, "c_1": 0}),
        _protein("A 1;B", "G", {"t1_1": 5, "t1_2": 5, "c_1": 0}),
    ]

    assert list(build(rows).data["Short protein IDs"]) == ["A-1;B;C;+1", "A-1;B"]


def test_a_missing_gene_name_is_backfilled_from_the_identifier(standard):
    """Prey rows are labeled by gene; an empty label would read as an unnamed prey."""
    short = dict(zip(standard.data["Majority protein IDs"],
                     standard.data["Short Gene names"]))

    assert short["P8"] == "P8"


def test_the_original_identifier_columns_survive_truncation(standard):
    """Only the Short columns are rewritten; annotation joins on the originals."""
    assert "P6a;P6b;P6c;P6d" in set(standard.data["Majority protein IDs"])


# --- SAINT output files ----------------------------------------------------------

def test_the_prey_file_has_one_row_per_surviving_protein(standard, out_dir):
    standard.to_SAINT(str(out_dir))
    prey = pd.read_csv(out_dir / "prey.txt", sep="\t", header=None, names=["Prey", "Gene"])

    assert len(prey) == SURVIVING_PROTEINS
    assert "P6a;P6b;P6c;+1" in set(prey["Prey"])


def test_a_spectral_count_prey_file_carries_the_sequence_length(build, out_dir):
    """SAINTexpress-spc normalizes counts by protein length and reads it from the
    middle column of prey.txt; intensity input is written without it."""
    rows = [_protein("P1", "G1", {"t1_1": 5, "t1_2": 5, "c_1": 0}, length=321)]

    build(rows, quant="Spectral Counts").to_SAINT(str(out_dir))
    prey = pd.read_csv(out_dir / "prey.txt", sep="\t", header=None)

    assert prey.values.tolist() == [["P1", 321, "G1"]]


def test_an_intensity_prey_file_omits_the_sequence_length(standard, out_dir):
    standard.to_SAINT(str(out_dir))
    prey = pd.read_csv(out_dir / "prey.txt", sep="\t", header=None)

    assert prey.shape[1] == 2


def test_the_bait_file_covers_every_design_row_including_controls(standard, out_dir):
    """bait.txt is built from the design alone, so it is unaffected by protein filtering."""
    standard.to_SAINT(str(out_dir))
    bait = pd.read_csv(out_dir / "bait.txt", sep="\t", header=None,
                       names=["Experiment", "Bait", "Type"])

    assert list(bait["Experiment"]) == ["t1_1", "t1_2", "c_1"]
    assert list(bait["Type"]) == ["T", "T", "C"]


def test_the_interaction_file_is_dense(standard, out_dir):
    """SAINT reads a value for every prey in every experiment; zeros are written here and
    dropped later when filtered_interaction.txt is built."""
    standard.to_SAINT(str(out_dir))
    interaction = _interaction(out_dir)

    assert len(interaction) == SURVIVING_PROTEINS * 3
    assert 0.0 in set(interaction["Value"])


def test_the_interaction_file_names_the_bait_not_the_experiment(standard, out_dir):
    standard.to_SAINT(str(out_dir))
    interaction = _interaction(out_dir)

    assert set(interaction[interaction["Experiment"] == "t1_1"]["Bait"]) == {"BaitA"}
    assert set(interaction[interaction["Experiment"] == "c_1"]["Bait"]) == {"Ctrl"}


# --- CompPASS output -------------------------------------------------------------

def test_the_comppass_file_is_sparse(standard, out_dir):
    """Only strictly positive cells are written -- the complement of interaction.txt's
    dense output."""
    standard.to_SAINT(str(out_dir))
    standard.to_CompPASS(str(out_dir))
    written = _comppass(out_dir)

    assert list(written.columns) == [
        "Experiment.ID", "Replicate", "Bait", "Prey", "Prey.Name", "Spectral.Count"]
    assert len(written) == int((_interaction(out_dir)["Value"] > 0).sum()) == 6


def test_experiment_id_holds_the_bait_name_and_bait_holds_its_protein_id(standard,
                                                                        out_dir):
    """score_compPass flags a self-interaction by comparing Bait against Prey, so both
    must be protein IDs, and it groups on Experiment.ID, which must be the bait name."""
    standard.to_CompPASS(str(out_dir))
    written = _comppass(out_dir)

    assert set(written["Experiment.ID"]) == {"BaitA"}
    assert set(written["Bait"]) == {"P_A"}


def test_the_prey_name_is_the_untruncated_gene_list(standard, out_dir):
    """Prey is the short form that joins against prey.txt; Prey.Name is the full gene
    list, so a group of more than three genes appears abbreviated in one and in full in
    the other."""
    standard.to_CompPASS(str(out_dir))
    written = _comppass(out_dir)

    assert "P6a;P6b;P6c;+1" in set(written["Prey"])
    assert "G6a;G6b;G6c;G6d" in set(written["Prey.Name"].dropna())
