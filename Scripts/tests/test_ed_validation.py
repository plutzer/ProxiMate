"""Tests for pre-parse input validation.

Validation runs before anything is written, and it is the only place a user's mistake
is turned into an actionable message rather than a traceback.  The Group rules are the
substantive ones: they decide which controls back which baits.
"""

import pandas as pd
import pytest

from ed_exceptions import (
    EDDuplicateExperimentError,
    EDFileEmptyError,
    EDFileFormatError,
    EDFileNotFoundError,
    EDInvalidGroupError,
    EDInvalidReplicateError,
    EDInvalidTypeError,
    EDMissingColumnError,
    EDMissingValueError,
    EDPGMismatchError,
    FPMissingColumnError,
    PGFileError,
    PGFileNotFoundError,
    PGMissingColumnError,
    ProxiMateError,
)
from ed_validation import (
    EDPGCrossValidator,
    EDValidator,
    FPValidator,
    PGValidator,
    validate_fragpipe_inputs,
    validate_maxquant_inputs,
    validate_msstats_inputs,
    validate_pioneer_inputs,
)


def _ed_frame(rows=None):
    return pd.DataFrame(rows if rows is not None else [
        {"Experiment Name": "t1_1", "Type": "T", "Bait": "BaitA", "Replicate": 1},
        {"Experiment Name": "t1_2", "Type": "T", "Bait": "BaitA", "Replicate": 2},
        {"Experiment Name": "c_1", "Type": "C", "Bait": "Ctrl", "Replicate": 1},
    ])


def _write(tmp_path, frame, name="ED.csv", sep=","):
    path = tmp_path / name
    frame.to_csv(path, index=False, sep=sep)
    return str(path)


# --- the error contract ------------------------------------------------------------

@pytest.mark.parametrize("cls, args", [
    (EDFileNotFoundError, ("/tmp/ED.csv",)),
    (EDFileEmptyError, ("/tmp/ED.csv",)),
    (EDFileFormatError, ("bad delimiter",)),
    (EDMissingColumnError, (["Bait"],)),
    (EDInvalidTypeError, ([2, 3],)),
    (EDInvalidReplicateError, ([4],)),
    (EDMissingValueError, ("Bait", [2])),
    (EDDuplicateExperimentError, (["t1"],)),
    (EDInvalidGroupError, ("invalid_group_value", [2])),
    (EDInvalidGroupError, ("bait_group_has_no_control", [2])),
    (PGFileNotFoundError, ("/tmp/pg.txt",)),
    (PGMissingColumnError, (["Reverse"],)),
    (FPMissingColumnError, (["Gene"],)),
    (EDPGMismatchError, (["c_1"], [])),
], ids=lambda v: getattr(v, "__name__", None))
def test_every_input_error_is_a_proximate_error_with_a_user_message(cls, args):
    """The GUI's handler keys on the base class and shows user_message."""
    error = cls(*args)

    assert isinstance(error, ProxiMateError)
    assert error.user_message
    assert error.suggestions


# --- file reading ------------------------------------------------------------------

@pytest.mark.parametrize("encoding", ["utf-8-sig", "cp1252"])
def test_designs_exported_from_excel_are_read(tmp_path, encoding):
    path = tmp_path / "ED.csv"
    path.write_bytes(
        "Experiment Name,Type,Bait,Replicate\nrun_é,T,BaitA,1\n".encode(encoding))

    assert len(EDValidator.validate_file_format(str(path))) == 1


def test_a_header_only_file_is_reported_as_empty(tmp_path):
    """A design with no rows is empty, not malformed."""
    path = tmp_path / "ED.csv"
    path.write_text("Experiment Name,Type,Bait,Replicate\n")

    with pytest.raises(EDFileEmptyError):
        EDValidator.validate_file_format(str(path))


# --- columns and values ------------------------------------------------------------

def test_a_missing_required_column_is_reported():
    with pytest.raises(EDMissingColumnError) as excinfo:
        EDValidator.validate_required_columns(_ed_frame().drop(columns=["Bait"]))

    assert excinfo.value.missing_columns == ["Bait"]


@pytest.mark.parametrize("value", [None, "   "])
def test_an_empty_required_value_is_reported(value):
    frame = _ed_frame()
    frame.loc[1, "Bait"] = value

    with pytest.raises(EDMissingValueError) as excinfo:
        EDValidator.validate_column_values(frame)

    assert excinfo.value.column_name == "Bait"


@pytest.mark.parametrize("bad", ["c", "X"])
def test_the_type_column_accepts_only_uppercase_c_and_t(bad):
    frame = _ed_frame()
    frame.loc[1, "Type"] = bad

    with pytest.raises(EDInvalidTypeError):
        EDValidator.validate_column_values(frame)


@pytest.mark.parametrize("bad", [0, "first"])
def test_a_replicate_must_be_a_positive_integer(bad):
    frame = _ed_frame([
        {"Experiment Name": "t1_1", "Type": "T", "Bait": "BaitA", "Replicate": 1},
        {"Experiment Name": "t1_2", "Type": "T", "Bait": "BaitA", "Replicate": bad},
    ])

    with pytest.raises(EDInvalidReplicateError):
        EDValidator.validate_column_values(frame)


def test_a_replicate_written_as_a_float_is_accepted():
    """A blank cell anywhere in the column makes pandas infer float64, and int(2.0)
    succeeds, so an otherwise valid file is not rejected for its dtype."""
    frame = _ed_frame()
    frame["Replicate"] = frame["Replicate"].astype(float)

    EDValidator.validate_column_values(frame)


def test_duplicate_experiment_names_are_rejected():
    """Downstream they would overwrite each other, so one run would vanish."""
    frame = _ed_frame()
    frame.loc[1, "Experiment Name"] = "t1_1"

    with pytest.raises(EDDuplicateExperimentError) as excinfo:
        EDValidator.validate_column_values(frame)

    assert excinfo.value.duplicates == ["t1_1"]


# --- Group column ------------------------------------------------------------------

def _grouped(groups, types=("T", "T", "C")):
    return pd.DataFrame([
        {"Experiment Name": "e{}".format(i), "Type": t, "Bait": "B", "Replicate": 1,
         "Group": g}
        for i, (t, g) in enumerate(zip(types, groups))
    ])


def test_an_entirely_empty_group_column_is_not_checked():
    """A column left over from a template must not force a grouped run."""
    EDValidator.validate_group_column(_grouped([None, None, None]))


@pytest.mark.parametrize("controls", ["*", "1,2"])
def test_controls_covering_every_test_group_pass(controls):
    EDValidator.validate_group_column(_grouped(["1", "2", controls]))


@pytest.mark.parametrize("bad", ["0", "x"])
def test_a_malformed_group_value_is_rejected(bad):
    with pytest.raises(EDInvalidGroupError) as excinfo:
        EDValidator.validate_group_column(_grouped([bad, "2", "*"]))

    assert excinfo.value.reason == "invalid_group_value"


def test_a_test_row_may_not_use_the_wildcard():
    """"*" means "backs every group", which is a statement about a control."""
    with pytest.raises(EDInvalidGroupError) as excinfo:
        EDValidator.validate_group_column(_grouped(["*", "2", "*"]))

    assert excinfo.value.reason == "test_row_wildcard"


def test_a_test_row_may_not_name_several_groups():
    """A bait scored under two sets of controls would appear twice in the output."""
    with pytest.raises(EDInvalidGroupError) as excinfo:
        EDValidator.validate_group_column(_grouped(["1,2", "2", "*"]))

    assert excinfo.value.reason == "test_row_multi_group"


def test_a_control_may_not_mix_the_wildcard_with_group_numbers():
    with pytest.raises(EDInvalidGroupError) as excinfo:
        EDValidator.validate_group_column(_grouped(["1", "2", "*,1"]))

    assert excinfo.value.reason == "control_wildcard_with_explicit_groups"


def test_every_test_row_must_declare_a_group_once_any_row_does():
    """An undeclared bait would be dropped from every group and scored by nothing."""
    with pytest.raises(EDInvalidGroupError) as excinfo:
        EDValidator.validate_group_column(_grouped(["1", None, "*"]))

    assert excinfo.value.reason == "test_row_missing_group_in_grouped_run"


def test_a_test_group_with_no_control_is_rejected():
    """SAINT needs controls; a group without any would score its baits against nothing."""
    with pytest.raises(EDInvalidGroupError) as excinfo:
        EDValidator.validate_group_column(_grouped(["1", "2", "1"]))

    assert excinfo.value.reason == "bait_group_has_no_control"
    assert excinfo.value.offenders == [2]


# --- data-file columns -------------------------------------------------------------

PG_COLUMNS = ["Majority protein IDs", "Gene names", "Reverse",
              "Only identified by site", "Potential contaminant"]
FP_COLUMNS = ["Protein", "Protein ID", "Gene", "Protein Length"]


def _pg_frame():
    return pd.DataFrame([{c: "x" for c in PG_COLUMNS}])


def test_a_protein_groups_column_is_required():
    with pytest.raises(PGMissingColumnError):
        PGValidator.validate_required_columns(_pg_frame().drop(columns=["Gene names"]))


def test_a_decoy_column_satisfies_the_reverse_requirement():
    """MaxQuant 2.4 and later name the column Decoy."""
    PGValidator.validate_required_columns(_pg_frame().rename(columns={"Reverse": "Decoy"}))


def test_a_fragpipe_column_is_required():
    frame = pd.DataFrame([{c: "x" for c in FP_COLUMNS}]).drop(columns=["Gene"])

    with pytest.raises(FPMissingColumnError):
        FPValidator.validate_required_columns(frame)


def test_a_pioneer_table_without_a_protein_column_is_rejected(tmp_path):
    path = tmp_path / "protein_groups_wide.tsv"
    path.write_text("gene_names\tt1_1\tt1_2\tc_1\nG1\t1\t2\t3\n")

    with pytest.raises(PGFileError):
        validate_pioneer_inputs(_write(tmp_path, _ed_frame()), str(path))


def _msstats_frame(runs):
    return pd.DataFrame({
        "Protein": ["P1"] * len(runs), "originalRUN": runs,
        "GROUP": ["g"] * len(runs), "SUBJECT": ["s"] * len(runs),
        "LABEL": ["L"] * len(runs), "LogIntensities": [8.0] * len(runs),
    })


@pytest.mark.parametrize("column", ["Protein", "originalRUN", "LABEL", "LogIntensities"])
def test_an_msstats_table_missing_a_required_column_is_rejected(tmp_path, column):
    frame = _msstats_frame(["t1_1", "t1_2", "c_1"]).drop(columns=[column])
    path = _write(tmp_path, frame, name="ProteinLevelData.csv")

    with pytest.raises(PGFileError) as excinfo:
        validate_msstats_inputs(_write(tmp_path, _ed_frame()), path)

    assert column in excinfo.value.user_message


def test_a_valid_msstats_table_is_returned_with_the_design(tmp_path):
    path = _write(tmp_path, _msstats_frame(["t1_1", "t1_2", "c_1"]),
                  name="ProteinLevelData.csv")

    ed_df, msstats_df = validate_msstats_inputs(_write(tmp_path, _ed_frame()), path)

    assert len(ed_df) == 3
    assert list(msstats_df["originalRUN"]) == ["t1_1", "t1_2", "c_1"]


# --- design vs data-file experiments ----------------------------------------------

def test_a_design_experiment_missing_from_protein_groups_is_rejected():
    pg = pd.DataFrame(columns=["Majority protein IDs", "Intensity t1_1"])

    with pytest.raises(EDPGMismatchError) as excinfo:
        EDPGCrossValidator.validate_experiment_match(_ed_frame(), pg, "Intensity ")

    assert excinfo.value.ed_only == ["c_1", "t1_2"]


def test_an_extra_data_column_does_not_raise():
    """Running a subset of an acquisition is normal; the design decides what is scored."""
    pg = pd.DataFrame(columns=["Intensity t1_1", "Intensity t1_2", "Intensity c_1",
                               "Intensity spare"])

    EDPGCrossValidator.validate_experiment_match(_ed_frame(), pg, "Intensity ")


def test_the_quantification_prefix_selects_the_columns():
    """An LFQ design validated against Intensity columns has no matches at all."""
    pg = pd.DataFrame(columns=["Intensity t1_1", "Intensity t1_2", "Intensity c_1"])

    with pytest.raises(EDPGMismatchError):
        EDPGCrossValidator.validate_experiment_match(_ed_frame(), pg, "LFQ intensity ")


def test_diann_metadata_columns_are_not_treated_as_runs():
    EDPGCrossValidator.validate_diann_match(
        _ed_frame(), pd.DataFrame(columns=["Protein.Group", "Genes", "t1_1", "t1_2", "c_1"]))

    with pytest.raises(EDPGMismatchError) as excinfo:
        EDPGCrossValidator.validate_diann_match(
            _ed_frame(), pd.DataFrame(columns=["Protein.Group", "Genes", "t1_1", "t1_2"]))

    assert excinfo.value.ed_only == ["c_1"]


def test_pioneer_metadata_columns_are_not_treated_as_runs():
    EDPGCrossValidator.validate_pioneer_match(
        _ed_frame(), pd.DataFrame(columns=["species", "gene_names", "protein", "target",
                                           "global_qval", "t1_1", "t1_2", "c_1", "extra"]))

    with pytest.raises(EDPGMismatchError) as excinfo:
        EDPGCrossValidator.validate_pioneer_match(
            _ed_frame(), pd.DataFrame(columns=["protein", "target", "t1_1", "t1_2"]))

    assert excinfo.value.ed_only == ["c_1"]


def test_msstats_run_names_are_compared_as_strings():
    """Numeric-looking run names are read as integers on one side and strings on the
    other, which would otherwise make every name look missing."""
    ed = _ed_frame([
        {"Experiment Name": 1, "Type": "T", "Bait": "B", "Replicate": 1},
        {"Experiment Name": 2, "Type": "C", "Bait": "C", "Replicate": 1},
    ])

    EDPGCrossValidator.validate_msstats_match(ed, pd.DataFrame({"originalRUN": ["1", "2"]}))


@pytest.mark.parametrize("suffix", [
    " Intensity", " MaxLFQ Intensity", " Total Spectral Count",
])
def test_a_fragpipe_quantification_suffix_yields_the_sample_name(suffix):
    """Suffixes are tried longest-first: " Total Spectral Count" also ends with
    " Spectral Count", which would leave "Sample Total" as the name."""
    fp = pd.DataFrame(columns=["Protein", "Protein ID",
                               "t1_1" + suffix, "t1_2" + suffix, "c_1" + suffix])

    EDPGCrossValidator.validate_fragpipe_match(_ed_frame(), fp)


def test_fragpipe_metadata_columns_are_not_treated_as_samples():
    EDPGCrossValidator.validate_fragpipe_match(
        _ed_frame(), pd.DataFrame(columns=["Protein", "Combined Total Spectral Count",
                                           "t1_1 Intensity", "t1_2 Intensity", "c_1 Intensity"]))

    with pytest.raises(EDPGMismatchError) as excinfo:
        EDPGCrossValidator.validate_fragpipe_match(
            _ed_frame(), pd.DataFrame(columns=["Protein", "t1_1 Intensity", "t1_2 Intensity"]))

    assert excinfo.value.ed_only == ["c_1"]


# --- orchestrators -------------------------------------------------------------------

def test_maxquant_validation_returns_both_frames(tmp_path):
    pg = _pg_frame()
    for experiment in ("t1_1", "t1_2", "c_1"):
        pg["Intensity {}".format(experiment)] = 1

    ed_df, pg_df = validate_maxquant_inputs(
        _write(tmp_path, _ed_frame()),
        _write(tmp_path, pg, name="proteinGroups.txt", sep="\t"),
        "Intensity")

    assert list(ed_df["Experiment Name"]) == ["t1_1", "t1_2", "c_1"]
    assert "Intensity t1_1" in pg_df.columns


def test_an_unknown_quantification_is_rejected(tmp_path):
    """Every quantification names a column prefix; one with no prefix cannot be
    validated as anything else."""
    pg = _pg_frame()
    for experiment in ("t1_1", "t1_2", "c_1"):
        pg["Intensity {}".format(experiment)] = 1

    with pytest.raises(ValueError):
        validate_maxquant_inputs(
            _write(tmp_path, _ed_frame()),
            _write(tmp_path, pg, name="proteinGroups.txt", sep="\t"),
            "iBAQ")


def test_fragpipe_validation_returns_both_frames(tmp_path):
    fp = pd.DataFrame([{c: "x" for c in FP_COLUMNS}])
    for experiment in ("t1_1", "t1_2", "c_1"):
        fp["{} Intensity".format(experiment)] = 1

    ed_df, fp_df = validate_fragpipe_inputs(
        _write(tmp_path, _ed_frame()),
        _write(tmp_path, fp, name="combined_protein.tsv", sep="\t"))

    assert len(ed_df) == 3
    assert "Protein ID" in fp_df.columns
