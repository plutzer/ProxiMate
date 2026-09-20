"""Tests for the FragPipe, DIA-NN, Pioneer and MSstats to-MaxQuant converters.

Each reads only ``experimental_design.name2experiment``, so the design here is a stub
rather than a parsed file.
"""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import parse
from ed_exceptions import ProxiMateError


def _design(*names):
    return SimpleNamespace(name2experiment={name: object() for name in names})


def _intensity_columns(frame):
    return [c for c in frame.columns if c.startswith("Intensity ")]


# --- FragPipe --------------------------------------------------------------------

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


BOTH_INTENSITY_COLUMNS = {"S1 Intensity": [10, 20], "S1 MaxLFQ Intensity": [11, 21],
                          "S1 Total Spectral Count": [3, 4]}


@pytest.mark.parametrize("quant_type, mq_column, values", [
    ("Intensity", "Intensity S1", [10, 20]),
    ("LFQ", "LFQ intensity S1", [11, 21]),
    ("Spectral Counts", "MS/MS count S1", [3, 4]),
])
def test_each_quantification_reads_its_own_column(tmp_path, quant_type, mq_column, values):
    """"S1 MaxLFQ Intensity" also ends with " Intensity"; in Intensity mode the plain
    column wins and no "S1 MaxLFQ" experiment is invented."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame(BOTH_INTENSITY_COLUMNS))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), quant_type)

    assert list(converted[mq_column]) == values
    assert len([c for c in converted.columns if c.endswith("S1") or "S1 " in c]) == 1


def test_a_sample_absent_from_the_design_is_dropped(tmp_path):
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({
        "S1 Intensity": [10, 20],
        "S2 Intensity": [30, 40],
    }))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert _intensity_columns(converted) == ["Intensity S1"]


def test_matching_no_design_experiment_is_rejected(tmp_path):
    """Usually a quantType that does not match the export; a frame with no
    quantification would otherwise reach scoring as all-zero runs."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Total Spectral Count": [3, 4]}))

    with pytest.raises(ProxiMateError):
        parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")


def test_contaminants_are_flagged_from_the_protein_column(tmp_path):
    """FragPipe marks contaminants with a prefix rather than a dedicated column."""
    frame = _fragpipe_frame({"S1 Intensity": [10, 20]})
    frame.loc[0, "Protein"] = "contam_sp|P1|X_HUMAN"
    fp = _write_fragpipe(tmp_path, frame)

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert list(converted["Potential contaminant"]) == ["+", "-"]


def test_missing_quantification_becomes_zero(tmp_path):
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Intensity": [10, np.nan]}))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert list(converted["Intensity S1"]) == [10.0, 0.0]


def test_the_identifier_is_the_accession_and_the_length_is_real(tmp_path):
    """"Protein" is the full FASTA header; "Protein ID" is the accession the rest of the
    pipeline joins on.  The protein length feeds the spectral-count prey file."""
    fp = _write_fragpipe(tmp_path, _fragpipe_frame({"S1 Intensity": [10, 20]}))

    converted = parse.convert_fragpipe_to_maxquant_format(fp, _design("S1"), "Intensity")

    assert list(converted["Majority protein IDs"]) == ["P1", "P2"]
    assert list(converted["Sequence length"]) == [100, 101]


# --- DIA-NN ----------------------------------------------------------------------

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


def test_diann_run_columns_gain_an_intensity_prefix(tmp_path):
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10, 20]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert list(converted["Intensity run_a"]) == [10, 20]
    assert list(converted["Majority protein IDs"]) == ["P1", "P2"]


def test_diann_metadata_and_undesigned_columns_are_not_runs(tmp_path):
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10, 20], "run_b": [30, 40]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert _intensity_columns(converted) == ["Intensity run_a"]


def test_diann_missing_intensities_become_zero(tmp_path):
    diann = _write_diann(tmp_path, _diann_frame({"run_a": [10.0, np.nan]}))

    converted = parse.convert_diann_to_maxquant_format(diann, _design("run_a"))

    assert list(converted["Intensity run_a"]) == [10.0, 0.0]


# --- Pioneer ---------------------------------------------------------------------

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


def test_pioneer_metadata_and_undesigned_columns_are_not_runs(tmp_path):
    pioneer = _write_pioneer(tmp_path, _pioneer_frame({"run_a": [10, 20], "run_b": [30, 40]}))

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert _intensity_columns(converted) == ["Intensity run_a"]


def test_empty_pioneer_cells_become_zero(tmp_path):
    pioneer = _write_pioneer(tmp_path, _pioneer_frame({"run_a": [10.0, np.nan]}))

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert list(converted["Intensity run_a"]) == [10.0, 0.0]


def test_decoy_and_entrapment_groups_are_dropped(tmp_path):
    frame = _pioneer_frame({"run_a": [10, 20, 30, 40]},
                           proteins=("P1", "DECOY", "ENTRAP", "P2"))
    frame.loc[1, "target"] = False
    frame.loc[2, "entrap_id"] = 1
    pioneer = _write_pioneer(tmp_path, frame)

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert list(converted["Majority protein IDs"]) == ["P1", "P2"]
    assert list(converted["Intensity run_a"]) == [10, 40]


def test_the_pioneer_flag_columns_are_optional(tmp_path):
    """Pioneer's output schema policy can omit target and entrap_id."""
    frame = _pioneer_frame({"run_a": [10, 20]}).drop(columns=["target", "entrap_id"])
    pioneer = _write_pioneer(tmp_path, frame)

    converted = parse.convert_pioneer_to_maxquant_format(pioneer, _design("run_a"))

    assert len(converted) == 2


# --- MSstats ---------------------------------------------------------------------

def _msstats_rows(*rows):
    """Rows of (Protein, originalRUN, LogIntensities[, LABEL]) as a ProteinLevelData frame."""
    return pd.DataFrame([
        {"Protein": p, "originalRUN": r, "LogIntensities": v,
         "LABEL": row[3] if len(row) > 3 else "L", "GROUP": "g", "SUBJECT": "s"}
        for row in rows for (p, r, v) in [row[:3]]
    ])


def _write_msstats(tmp_path, frame, encoding="utf-8"):
    path = tmp_path / "ProteinLevelData.csv"
    path.write_bytes(frame.to_csv(index=False).encode(encoding))
    return str(path)


def test_log_intensities_are_back_transformed_exactly(tmp_path):
    """MSstats reports log2(normalized abundance); the pipeline expects linear values."""
    path = _write_msstats(tmp_path, _msstats_rows(("P1", "r1", 10.0), ("P1", "r2", 3.5)))

    converted = parse.convert_msstats_to_maxquant_format(path, _design("r1", "r2"))

    assert list(converted["Intensity r1"]) == [2.0 ** 10.0]
    assert list(converted["Intensity r2"]) == [2.0 ** 3.5]


def test_heavy_label_rows_are_dropped(tmp_path):
    path = _write_msstats(tmp_path, _msstats_rows(
        ("P1", "r1", 4.0, "L"), ("P2", "r1", 9.0, "H")))

    converted = parse.convert_msstats_to_maxquant_format(path, _design("r1"))

    assert list(converted["Majority protein IDs"]) == ["P1"]


def test_missing_log_intensities_are_dropped_and_become_zero(tmp_path):
    path = _write_msstats(tmp_path, _msstats_rows(
        ("P1", "r1", 4.0), ("P1", "r2", np.nan), ("P2", "r2", 5.0)))

    converted = parse.convert_msstats_to_maxquant_format(path, _design("r1", "r2"))

    assert converted.set_index("Majority protein IDs").loc["P1", "Intensity r2"] == 0.0


def test_duplicate_protein_run_rows_are_averaged_on_the_log_scale(tmp_path):
    path = _write_msstats(tmp_path, _msstats_rows(("P1", "r1", 2.0), ("P1", "r1", 4.0)))

    converted = parse.convert_msstats_to_maxquant_format(path, _design("r1"))

    assert list(converted["Intensity r1"]) == [2.0 ** 3.0]


def test_runs_absent_from_the_design_are_dropped(tmp_path):
    path = _write_msstats(tmp_path, _msstats_rows(("P1", "r1", 2.0), ("P1", "r2", 4.0)))

    converted = parse.convert_msstats_to_maxquant_format(path, _design("r1"))

    assert _intensity_columns(converted) == ["Intensity r1"]


@pytest.mark.parametrize("encoding", ["utf-8-sig", "latin-1"])
def test_msstats_tables_in_other_encodings_are_read(tmp_path, encoding):
    path = _write_msstats(tmp_path, _msstats_rows(("Café_HUMAN", "r1", 2.0)),
                          encoding=encoding)

    converted = parse.convert_msstats_to_maxquant_format(path, _design("r1"))

    assert list(converted["Majority protein IDs"]) == ["Café_HUMAN"]
