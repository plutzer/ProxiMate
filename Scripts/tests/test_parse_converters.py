"""Tests for the DIA-NN, Pioneer and FragPipe to-MaxQuant converters.

Neither converter raises when it fails to recognize a sample: an unmatched column is
dropped and the run continues with one fewer experiment, which reaches SAINT as an
absence rather than as an error.  These tests pin what each one matches, what it renames,
and what it silently discards.

All three read only ``experimental_design.name2experiment``, so the design here is a
stub rather than a parsed file.
"""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import parse


def _design(*names):
    return SimpleNamespace(name2experiment={name: object() for name in names})


# --- FragPipe ------------------------------------------------------------------

def _fragpipe_frame(quant_columns, proteins=("P1", "P2")):
    """A combined_protein.tsv-shaped frame; `quant_columns` maps column name to values."""
    frame = pd.DataFrame({
        "Protein": ["sp|{}|X_HUMAN".format(p) for p in proteins],
        "Protein ID": list(proteins),
        "Gene": ["G_{}".format(p) for p in proteins],
        "Protein Length": [100 + i for i, _ in enumerate(proteins)],
        "Description": ["desc {}".format(p) for p in proteins],
    })
    for name, values in quant_columns.items():
        frame[name] = values
    return frame


def _write_fragpipe(tmp_path, frame, name="combined_protein.tsv"):
    path = tmp_path / name
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.mark.parametrize("quant_type, fp_suffix, mq_prefix", [
    ("Intensity", " Intensity", "Intensity "),
    ("LFQ", " MaxLFQ Intensity", "LFQ intensity "),
    ("Spectral Counts", " Total Spectral Count", "MS/MS count "),
])
def test_each_quantification_reads_its_own_column(tmp_path, quant_type, fp_suffix,
                                                  mq_prefix):
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({
        "S1{}".format(fp_suffix): [10, 20],
    }))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), quant_type)

    assert "{}S1".format(mq_prefix) in converted.columns
    assert list(converted["{}S1".format(mq_prefix)]) == [10, 20]


def test_intensity_mode_reads_the_plain_column_when_both_are_present(tmp_path):
    """"S1 MaxLFQ Intensity" also ends with " Intensity", but it yields the sample name
    "S1 MaxLFQ", which the design lookup rejects.  That lookup, not the suffix match, is
    what keeps LFQ values out of an Intensity run in the ordinary case."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({
        "S1 Intensity": [10, 20],
        "S1 MaxLFQ Intensity": [11, 21],
    }))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert list(converted["Intensity S1"]) == [10, 20]
    assert "Intensity S1 MaxLFQ" not in converted.columns


def test_a_sample_named_for_maxlfq_is_unreadable_in_intensity_mode(tmp_path, caplog):
    """The explicit skip is reachable only here: a sample whose own name ends in
    " MaxLFQ" owns the column "<name> Intensity", which is spelled exactly like another
    sample's MaxLFQ column.  The skip resolves that collision in favor of the LFQ
    reading, so the sample is reported unmatched rather than quantified."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({
        "S1 MaxLFQ Intensity": [11, 21],
    }))

    with caplog.at_level("WARNING", logger="proximate.parse"):
        converted = parse.convert_fragpipe_to_maxquant_format(
            fp, _design("S1 MaxLFQ"), "Intensity")

    assert "Intensity S1 MaxLFQ" not in converted.columns
    assert "S1 MaxLFQ" in caplog.text


def test_lfq_mode_reads_the_maxlfq_column(tmp_path):
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({
        "S1 Intensity": [10, 20],
        "S1 MaxLFQ Intensity": [11, 21],
    }))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "LFQ")

    assert list(converted["LFQ intensity S1"]) == [11, 21]
    assert "Intensity S1" not in converted.columns


def test_a_sample_absent_from_the_design_is_dropped(tmp_path):
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({
        "S1 Intensity": [10, 20],
        "S2 Intensity": [30, 40],
    }))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert "Intensity S2" not in converted.columns


def test_a_design_experiment_with_no_column_is_reported(tmp_path, caplog):
    """No column at all is otherwise indistinguishable downstream from a run whose
    intensities happened to be zero."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Intensity": [10, 20]}))

    with caplog.at_level("WARNING", logger="proximate.parse"):
        converted = parse.convert_fragpipe_to_maxquant_format(
            fp, _design("S1", "S_missing"), "Intensity")

    assert "S_missing" in caplog.text
    assert "Intensity S_missing" not in converted.columns


def test_matching_nothing_is_reported_as_an_error(tmp_path, caplog):
    """Usually a quantType that does not match the export; the converter still returns a
    frame, so nothing downstream raises."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({
        "S1 Total Spectral Count": [3, 4]}))

    with caplog.at_level("ERROR", logger="proximate.parse"):
        converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"),
                                                              "Intensity")

    assert "No quantification columns matched" in caplog.text
    assert not [c for c in converted.columns if c.startswith("Intensity ")]


def test_an_unknown_quantification_raises(tmp_path):
    """A bare KeyError, not a ProxiMateError -- the value comes from a fixed GUI dropdown
    and the CLI, so an unrecognized one is a programming error rather than user input."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Intensity": [10, 20]}))

    with pytest.raises(KeyError):
        parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "iBAQ")


def test_contaminants_are_flagged_from_the_protein_column(tmp_path):
    """FragPipe marks contaminants with a prefix rather than a dedicated column."""
    frame = _fragpipe_frame({"S1 Intensity": [10, 20]})
    frame.loc[0, "Protein"] = "contam_sp|P1|X_HUMAN"
    fp = _write_fragpipe(tmp_path, frame)

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert list(converted["Potential contaminant"]) == ["+", "-"]


def test_a_missing_protein_value_is_treated_as_a_contaminant(tmp_path):
    """Documented, not fixed: str.startswith yields NaN for a missing value, which
    np.where reads as truthy, so the row is dropped as a contaminant downstream."""
    frame = _fragpipe_frame({"S1 Intensity": [10, 20]})
    frame.loc[0, "Protein"] = np.nan
    fp = _write_fragpipe(tmp_path, frame)

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert converted.loc[0, "Potential contaminant"] == "+"


def test_reverse_and_site_only_are_always_negative(tmp_path):
    """FragPipe reports neither, so the columns exist only to satisfy ProteinGroups."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Intensity": [10, 20]}))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert set(converted["Reverse"]) == {"-"}
    assert set(converted["Only identified by site"]) == {"-"}


def test_missing_quantification_becomes_zero(tmp_path):
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Intensity": [10, np.nan]}))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert list(converted["Intensity S1"]) == [10.0, 0.0]


def test_the_real_protein_length_is_carried_over(tmp_path):
    """FragPipe reports it, unlike DIA-NN."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Intensity": [10, 20]}))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert list(converted["Sequence length"]) == [100, 101]


def test_a_missing_gene_becomes_an_empty_string(tmp_path):
    """ProteinGroups backfills a null gene name from the protein ID, and an empty string
    is not null, so the name stays empty."""
    frame = _fragpipe_frame({"S1 Intensity": [10, 20]})
    frame.loc[0, "Gene"] = np.nan
    fp = _write_fragpipe(tmp_path, frame)

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert converted.loc[0, "Gene names"] == ""


def test_the_description_column_is_optional(tmp_path):
    frame = _fragpipe_frame({"S1 Intensity": [10, 20]}).drop(columns=["Description"])
    fp = _write_fragpipe(tmp_path, frame)

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert set(converted["Protein names"]) == {""}


def test_the_identifier_comes_from_protein_id_not_protein(tmp_path):
    """"Protein" is the full FASTA header; "Protein ID" is the accession the rest of the
    pipeline joins on."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Intensity": [10, 20]}))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert list(converted["Majority protein IDs"]) == ["P1", "P2"]


# --- DIA-NN --------------------------------------------------------------------

def _diann_frame(run_columns, proteins=("P1", "P2")):
    frame = pd.DataFrame({
        "Protein.Group": list(proteins),
        "Protein.Names": ["{}_HUMAN".format(p) for p in proteins],
        "Genes": ["G_{}".format(p) for p in proteins],
        "First.Protein.Description": ["desc"] * len(proteins),
        "N.Sequences": [5] * len(proteins),
        "N.Proteotypic.Sequences": [3] * len(proteins),
    })
    for name, values in run_columns.items():
        frame[name] = values
    return frame


def _write_diann(tmp_path, frame, name="report.pg_matrix.tsv"):
    path = tmp_path / name
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def test_run_columns_gain_an_intensity_prefix(tmp_path):
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10, 20]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert list(converted["Intensity run_a"]) == [10, 20]


def test_metadata_columns_are_never_treated_as_runs(tmp_path):
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10, 20]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert [c for c in converted.columns if c.startswith("Intensity ")] == [
        "Intensity run_a"]


def test_run_matching_is_exact(tmp_path):
    """DIA-NN names its columns after the raw file, commonly a full path.  A design name
    that differs at all yields no column, and nothing is logged."""
    diann = _write_diann(tmp_path, _diann_frame({"D:\\data\\run_a.raw": [10, 20]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert not [c for c in converted.columns if c.startswith("Intensity ")]


def test_a_column_absent_from_the_design_is_dropped(tmp_path):
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10, 20], "run_b": [30, 40]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert "Intensity run_b" not in converted.columns


def test_missing_intensities_become_zero(tmp_path):
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10.0, np.nan]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert list(converted["Intensity run_a"]) == [10.0, 0.0]


def test_the_sequence_length_is_a_placeholder(tmp_path):
    """report.pg_matrix.tsv carries no protein length.  Nothing reads the column today,
    since the prey file that would need it is written without one."""
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10, 20]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert set(converted["Sequence length"]) == {1}


def test_no_protein_is_ever_filtered(tmp_path):
    """DIA-NN reports no reverse, site-only or contaminant flags, so ProteinGroups'
    three filters remove nothing from a DIA-NN run."""
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10, 20]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    for column in ("Reverse", "Only identified by site", "Potential contaminant"):
        assert set(converted[column]) == {"-"}


@pytest.mark.parametrize("column", ["Protein.Names", "Genes"])
def test_the_annotation_columns_are_optional(tmp_path, column):
    frame = _diann_frame({"run_a": [10, 20]}).drop(columns=[column])
    diann = _write_diann(tmp_path, frame)

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    target = "Protein names" if column == "Protein.Names" else "Gene names"
    assert set(converted[target]) == {""}


def test_the_protein_group_column_is_required(tmp_path):
    frame = _diann_frame({"run_a": [10, 20]}).drop(columns=["Protein.Group"])
    diann = _write_diann(tmp_path, frame)

    with pytest.raises(KeyError):
        parse.convert_diann_to_maxquant_format(diann, _design("run_a"))


# --- Pioneer -------------------------------------------------------------------

def _pioneer_frame(run_columns, proteins=("P1", "P2")):
    frame = pd.DataFrame({
        "species": ["HUMAN"] * len(proteins),
        "gene_names": ["G_{}".format(p) for p in proteins],
        "protein_names": ["{} protein".format(p) for p in proteins],
        "protein": list(proteins),
        "target": [True] * len(proteins),
        "entrap_id": [0] * len(proteins),
        "global_pg_score": [0.9] * len(proteins),
        "global_qval": [0.001] * len(proteins),
    })
    for name, values in run_columns.items():
        frame[name] = values
    return frame


def _write_pioneer(tmp_path, frame, name="protein_groups_wide.tsv"):
    path = tmp_path / name
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def test_pioneer_run_columns_gain_an_intensity_prefix(tmp_path):
    pioneer = _write_pioneer(tmp_path, _pioneer_frame({"run_a": [10, 20]}))

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert list(converted["Intensity run_a"]) == [10, 20]
    assert list(converted["Majority protein IDs"]) == ["P1", "P2"]
    assert list(converted["Gene names"]) == ["G_P1", "G_P2"]


def test_pioneer_metadata_columns_are_never_treated_as_runs(tmp_path):
    pioneer = _write_pioneer(tmp_path, _pioneer_frame({"run_a": [10, 20]}))

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert [c for c in converted.columns if c.startswith("Intensity ")] == ["Intensity run_a"]


def test_pioneer_run_matching_is_exact(tmp_path):
    """Pioneer names run columns after the MS file without its extension; a design name
    carrying the extension matches nothing."""
    pioneer = _write_pioneer(tmp_path, _pioneer_frame({"run_a": [10, 20]}))

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a.raw"))

    assert not [c for c in converted.columns if c.startswith("Intensity ")]


def test_a_pioneer_column_absent_from_the_design_is_dropped(tmp_path):
    pioneer = _write_pioneer(tmp_path, _pioneer_frame({"run_a": [10, 20], "run_b": [30, 40]}))

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert "Intensity run_b" not in converted.columns


def test_empty_pioneer_cells_become_zero(tmp_path):
    pioneer = _write_pioneer(tmp_path, _pioneer_frame({"run_a": [10.0, np.nan]}))

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert list(converted["Intensity run_a"]) == [10.0, 0.0]


def test_decoy_and_entrapment_groups_are_dropped(tmp_path):
    frame = _pioneer_frame({"run_a": [10, 20, 30, 40]}, proteins=("P1", "DECOY", "ENTRAP", "P2"))
    frame.loc[1, "target"] = False
    frame.loc[2, "entrap_id"] = 1
    pioneer = _write_pioneer(tmp_path, frame)

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert list(converted["Majority protein IDs"]) == ["P1", "P2"]
    assert list(converted["Intensity run_a"]) == [10, 40]
    for column in ("Reverse", "Only identified by site", "Potential contaminant"):
        assert set(converted[column]) == {"-"}


def test_the_pioneer_flag_columns_are_optional(tmp_path):
    """Pioneer's output schema policy can omit target and entrap_id."""
    frame = _pioneer_frame({"run_a": [10, 20]}).drop(columns=["target", "entrap_id"])
    pioneer = _write_pioneer(tmp_path, frame)

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert len(converted) == 2


def test_empty_pioneer_gene_names_become_blank(tmp_path):
    frame = _pioneer_frame({"run_a": [10, 20]})
    frame.loc[0, "gene_names"] = np.nan
    pioneer = _write_pioneer(tmp_path, frame)

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert list(converted["Gene names"]) == ["", "G_P2"]


def test_the_pioneer_protein_column_is_required(tmp_path):
    frame = _pioneer_frame({"run_a": [10, 20]}).drop(columns=["protein"])
    pioneer = _write_pioneer(tmp_path, frame)

    with pytest.raises(KeyError):
        parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))
