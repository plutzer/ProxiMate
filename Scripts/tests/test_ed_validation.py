"""Tests for pre-parse input validation.

Validation runs before anything is written, and it is the only place a user's mistake is
turned into an actionable message rather than a traceback.  It is also strictly fail-fast:
each rule raises on the first category of problem it finds, so the order the rules run in
decides which of several problems a user is told about.

The Group rules are the substantive ones.  They decide which controls back which baits,
and ``ed_validation`` implements that grammar a second time -- ``experimental_design``
has its own parser, tested alongside it.
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
)
from ed_validation import (
    EDPGCrossValidator,
    EDValidator,
    FPValidator,
    PGValidator,
    validate_diann_inputs,
    validate_fragpipe_inputs,
    validate_maxquant_inputs,
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


# --- validate_file_exists ------------------------------------------------------

def test_a_missing_file_is_reported(tmp_path):
    with pytest.raises(EDFileNotFoundError):
        EDValidator.validate_file_exists(str(tmp_path / "absent.csv"))


@pytest.mark.parametrize("filepath", ["", None])
def test_an_absent_path_is_reported_as_missing(filepath):
    with pytest.raises(EDFileNotFoundError):
        EDValidator.validate_file_exists(filepath)


def test_a_zero_byte_file_is_reported_as_empty(tmp_path):
    path = tmp_path / "ED.csv"
    path.write_text("")

    with pytest.raises(EDFileEmptyError):
        EDValidator.validate_file_exists(str(path))


# --- validate_file_format ------------------------------------------------------

@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig", "latin-1", "cp1252"])
def test_the_supported_encodings_are_read(tmp_path, encoding):
    """Designs are routinely exported from Excel, which does not always write UTF-8."""
    path = tmp_path / "ED.csv"
    path.write_bytes(
        "Experiment Name,Type,Bait,Replicate\nrun_é,T,BaitA,1\n".encode(encoding))

    frame = EDValidator.validate_file_format(str(path))

    assert len(frame) == 1


def test_a_header_only_file_is_reported_as_empty(tmp_path):
    """A design with no rows is empty, not malformed; telling a user their file's format
    is invalid sends them to look for the wrong problem."""
    path = tmp_path / "ED.csv"
    path.write_text("Experiment Name,Type,Bait,Replicate\n")

    with pytest.raises(EDFileEmptyError):
        EDValidator.validate_file_format(str(path))


def test_an_unparseable_file_is_reported_as_a_format_error(tmp_path):
    path = tmp_path / "ED.csv"
    path.write_text('Experiment Name,Type\n"unclosed,T\nmore,T,extra,cells,here\n')

    with pytest.raises(EDFileFormatError):
        EDValidator.validate_file_format(str(path))


# --- validate_required_columns -------------------------------------------------

@pytest.mark.parametrize("column", ["Experiment Name", "Type", "Bait", "Replicate"])
def test_each_required_column_is_checked(column):
    with pytest.raises(EDMissingColumnError):
        EDValidator.validate_required_columns(_ed_frame().drop(columns=[column]))


def test_every_missing_column_is_reported_at_once(_=None):
    """Unlike ExperimentalDesign, which stops at the first, so one pass here fixes the
    whole header."""
    frame = _ed_frame().drop(columns=["Bait", "Replicate"])

    with pytest.raises(EDMissingColumnError) as excinfo:
        EDValidator.validate_required_columns(frame)

    assert excinfo.value.missing_columns == ["Bait", "Replicate"]


def test_optional_columns_are_not_required():
    """Bait ID and Group are both optional here, though write_CompPASS needs Bait ID."""
    EDValidator.validate_required_columns(_ed_frame())


# --- validate_column_values ----------------------------------------------------

def test_a_valid_frame_passes():
    EDValidator.validate_column_values(_ed_frame())


@pytest.mark.parametrize("column", ["Experiment Name", "Type", "Bait", "Replicate"])
def test_an_empty_required_value_is_reported(column):
    frame = _ed_frame()
    frame.loc[1, column] = None

    with pytest.raises(EDMissingValueError) as excinfo:
        EDValidator.validate_column_values(frame)

    assert excinfo.value.column_name == column


def test_a_whitespace_only_value_counts_as_empty():
    frame = _ed_frame()
    frame.loc[1, "Bait"] = "   "

    with pytest.raises(EDMissingValueError):
        EDValidator.validate_column_values(frame)


def test_reported_rows_are_numbered_as_the_spreadsheet_shows_them():
    """One-based and counting the header, so the number matches what a user sees."""
    frame = _ed_frame()
    frame.loc[0, "Bait"] = None

    with pytest.raises(EDMissingValueError) as excinfo:
        EDValidator.validate_column_values(frame)

    assert excinfo.value.row_indices == [2]


@pytest.mark.parametrize("bad", ["X", "c", "t", "Control"])
def test_the_type_column_accepts_only_uppercase_c_and_t(bad):
    frame = _ed_frame()
    frame.loc[1, "Type"] = bad

    with pytest.raises(EDInvalidTypeError):
        EDValidator.validate_column_values(frame)


@pytest.mark.parametrize("bad", [0, -1, "first", "2.0"])
def test_a_replicate_must_be_a_positive_integer(bad):
    """"2.0" is rejected: int() parses it as a string, not as a number."""
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
    """Downstream they would silently overwrite each other, so one run would vanish."""
    frame = _ed_frame()
    frame.loc[1, "Experiment Name"] = "t1_1"

    with pytest.raises(EDDuplicateExperimentError) as excinfo:
        EDValidator.validate_column_values(frame)

    assert excinfo.value.duplicates == ["t1_1"]


# --- validate_group_column -----------------------------------------------------

def _grouped(groups, types=("T", "T", "C")):
    return pd.DataFrame([
        {"Experiment Name": "e{}".format(i), "Type": t, "Bait": "B", "Replicate": 1,
         "Group": g}
        for i, (t, g) in enumerate(zip(types, groups))
    ])


def test_a_design_with_no_group_column_is_not_checked():
    EDValidator.validate_group_column(_ed_frame())


def test_an_entirely_empty_group_column_is_not_checked():
    """A column left over from a template must not force a grouped run."""
    EDValidator.validate_group_column(_grouped([None, None, None]))


def test_a_well_formed_grouped_design_passes():
    EDValidator.validate_group_column(_grouped(["1", "2", "*"]))


def test_explicit_controls_cover_their_groups():
    EDValidator.validate_group_column(
        _grouped(["1", "2", "1,2"]))


@pytest.mark.parametrize("bad", ["0", "01", "-1", "1.0", "x", ","])
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
    """A bait scored under two sets of controls would appear twice in the merged output."""
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


def test_a_universal_control_covers_every_group():
    EDValidator.validate_group_column(_grouped(["1", "7", "*"]))


def test_a_control_for_an_unused_group_is_only_a_warning(caplog):
    """The run is still well-defined -- the control is simply never used -- so this is
    the one Group finding that does not stop the run."""
    with caplog.at_level("WARNING", logger="proximate.ed_validation"):
        EDValidator.validate_group_column(_grouped(["1", "1", "1,9"]))

    assert "9" in caplog.text


# --- rule ordering --------------------------------------------------------------

def test_malformed_values_are_reported_before_row_type_rules():
    """The four Group categories are all collected, then raised in a fixed order, so a
    file with several kinds of problem always reports the most basic one first."""
    with pytest.raises(EDInvalidGroupError) as excinfo:
        EDValidator.validate_group_column(_grouped(["x", "*", "*"]))

    assert excinfo.value.reason == "invalid_group_value"


def test_wildcard_misuse_is_reported_before_multi_group_test_rows():
    with pytest.raises(EDInvalidGroupError) as excinfo:
        EDValidator.validate_group_column(_grouped(["1,2", "1", "*,1"]))

    assert excinfo.value.reason == "control_wildcard_with_explicit_groups"


def test_column_values_are_checked_before_groups(tmp_path):
    """validate_ed_file runs its rules in a fixed order, so a file with both a bad Type
    and a bad Group reports the Type -- the Group cell may well be fine once the row is."""
    frame = _grouped(["1", "2", "*"])
    frame.loc[1, "Type"] = "X"
    frame.loc[1, "Group"] = "bad"

    with pytest.raises(EDInvalidTypeError):
        EDValidator.validate_ed_file(_write(tmp_path, frame))


def test_required_columns_are_checked_before_their_values(tmp_path):
    frame = _ed_frame().drop(columns=["Type"])

    with pytest.raises(EDMissingColumnError):
        EDValidator.validate_ed_file(_write(tmp_path, frame))


def test_a_valid_file_is_returned(tmp_path):
    frame = EDValidator.validate_ed_file(_write(tmp_path, _ed_frame()))

    assert list(frame["Experiment Name"]) == ["t1_1", "t1_2", "c_1"]


# --- PGValidator and FPValidator -------------------------------------------------

PG_COLUMNS = ["Majority protein IDs", "Gene names", "Reverse",
              "Only identified by site", "Potential contaminant"]


def _pg_frame():
    return pd.DataFrame([{c: "x" for c in PG_COLUMNS}])


def test_a_missing_protein_groups_file_is_reported(tmp_path):
    with pytest.raises(PGFileNotFoundError):
        PGValidator.validate_file_exists(str(tmp_path / "absent.txt"))


@pytest.mark.parametrize("column", PG_COLUMNS)
def test_each_protein_groups_column_is_required(column):
    with pytest.raises(PGMissingColumnError):
        PGValidator.validate_required_columns(_pg_frame().drop(columns=[column]))


def test_a_header_only_protein_groups_file_is_reported_as_empty(tmp_path):
    path = tmp_path / "proteinGroups.txt"
    path.write_text("\t".join(PG_COLUMNS) + "\n")

    with pytest.raises(PGFileError) as excinfo:
        PGValidator.validate_file_format(str(path))

    assert "empty" in excinfo.value.user_message.lower()


FP_COLUMNS = ["Protein", "Protein ID", "Gene", "Protein Length"]


@pytest.mark.parametrize("column", FP_COLUMNS)
def test_each_fragpipe_column_is_required(column):
    frame = pd.DataFrame([{c: "x" for c in FP_COLUMNS}]).drop(columns=[column])

    with pytest.raises(FPMissingColumnError):
        FPValidator.validate_required_columns(frame)


def test_a_header_only_fragpipe_file_is_reported_as_empty(tmp_path):
    path = tmp_path / "combined_protein.tsv"
    path.write_text("\t".join(FP_COLUMNS) + "\n")

    with pytest.raises(PGFileError) as excinfo:
        FPValidator.validate_file_format(str(path))

    assert "empty" in excinfo.value.user_message.lower()


# --- cross-validation ------------------------------------------------------------

def test_a_design_experiment_missing_from_protein_groups_is_rejected():
    pg = pd.DataFrame(columns=["Majority protein IDs", "Intensity t1_1"])

    with pytest.raises(EDPGMismatchError) as excinfo:
        EDPGCrossValidator.validate_experiment_match(_ed_frame(), pg, "Intensity ")

    assert excinfo.value.ed_only == ["c_1", "t1_2"]


def test_an_extra_data_column_is_only_a_warning(caplog):
    """Running a subset of an acquisition is normal; the design decides what is scored."""
    pg = pd.DataFrame(columns=["Intensity t1_1", "Intensity t1_2", "Intensity c_1",
                               "Intensity spare"])

    with caplog.at_level("WARNING", logger="proximate.ed_validation"):
        EDPGCrossValidator.validate_experiment_match(_ed_frame(), pg, "Intensity ")

    assert "spare" in caplog.text


def test_the_quantification_prefix_selects_the_columns():
    """An LFQ design validated against Intensity columns has no matches at all."""
    pg = pd.DataFrame(columns=["Intensity t1_1", "Intensity t1_2", "Intensity c_1"])

    with pytest.raises(EDPGMismatchError):
        EDPGCrossValidator.validate_experiment_match(_ed_frame(), pg, "LFQ intensity ")


def test_diann_columns_are_matched_by_name():
    diann = pd.DataFrame(columns=["Protein.Group", "Genes", "t1_1", "t1_2", "c_1"])

    EDPGCrossValidator.validate_diann_match(_ed_frame(), diann)


def test_diann_metadata_columns_are_not_treated_as_runs():
    diann = pd.DataFrame(columns=["Protein.Group", "Genes", "t1_1", "t1_2"])

    with pytest.raises(EDPGMismatchError) as excinfo:
        EDPGCrossValidator.validate_diann_match(_ed_frame(), diann)

    assert excinfo.value.ed_only == ["c_1"]


def test_pioneer_run_columns_are_matched_by_name():
    pioneer = pd.DataFrame(columns=["species", "gene_names", "protein", "target", "global_qval",
                                    "t1_1", "t1_2", "c_1", "extra_run"])

    EDPGCrossValidator.validate_pioneer_match(_ed_frame(), pioneer)


def test_pioneer_metadata_columns_are_not_treated_as_runs():
    pioneer = pd.DataFrame(columns=["protein", "target", "entrap_id", "t1_1", "t1_2"])

    with pytest.raises(EDPGMismatchError) as excinfo:
        EDPGCrossValidator.validate_pioneer_match(_ed_frame(), pioneer)

    assert excinfo.value.ed_only == ["c_1"]


def test_msstats_runs_are_matched_on_the_original_run_column():
    msstats = pd.DataFrame({"originalRUN": ["t1_1", "t1_2", "c_1", "c_1"]})

    EDPGCrossValidator.validate_msstats_match(_ed_frame(), msstats)


def test_msstats_run_names_are_compared_as_strings():
    """Numeric-looking run names are read as integers by pandas on one side and as
    strings on the other, which would make every name look missing."""
    ed = _ed_frame([
        {"Experiment Name": 1, "Type": "T", "Bait": "B", "Replicate": 1},
        {"Experiment Name": 2, "Type": "C", "Bait": "C", "Replicate": 1},
    ])
    msstats = pd.DataFrame({"originalRUN": ["1", "2"]})

    EDPGCrossValidator.validate_msstats_match(ed, msstats)


@pytest.mark.parametrize("suffix", [
    " Intensity", " MaxLFQ Intensity", " Spectral Count",
    " Total Spectral Count", " Unique Spectral Count",
])
def test_every_fragpipe_quantification_suffix_yields_the_sample_name(suffix):
    """The suffixes are tried longest-first: " Total Spectral Count" also ends with
    " Spectral Count", so the shorter one would leave "Sample Total" as the name and no
    design experiment would match."""
    fp = pd.DataFrame(columns=["Protein", "Protein ID",
                               "t1_1" + suffix, "t1_2" + suffix, "c_1" + suffix])

    EDPGCrossValidator.validate_fragpipe_match(_ed_frame(), fp)


def test_fragpipe_metadata_columns_are_not_treated_as_samples():
    fp = pd.DataFrame(columns=["Protein", "Protein ID", "Combined Total Spectral Count",
                               "t1_1 Intensity", "t1_2 Intensity", "c_1 Intensity"])

    EDPGCrossValidator.validate_fragpipe_match(_ed_frame(), fp)


def test_a_design_sample_missing_from_fragpipe_is_rejected():
    fp = pd.DataFrame(columns=["Protein", "t1_1 Intensity", "t1_2 Intensity"])

    with pytest.raises(EDPGMismatchError) as excinfo:
        EDPGCrossValidator.validate_fragpipe_match(_ed_frame(), fp)

    assert excinfo.value.ed_only == ["c_1"]


# --- orchestrators ---------------------------------------------------------------

def test_maxquant_validation_returns_both_frames(tmp_path):
    pg = _pg_frame()
    for experiment in ("t1_1", "t1_2", "c_1"):
        pg["Intensity {}".format(experiment)] = 1

    ed_df, pg_df = validate_maxquant_inputs(
        _write(tmp_path, _ed_frame()),
        _write(tmp_path, pg, name="proteinGroups.txt", sep="\t"),
        "Intensity")

    assert len(ed_df) == 3
    assert "Intensity t1_1" in pg_df.columns


def test_an_unknown_quantification_falls_back_to_intensity(tmp_path):
    """Documented, not fixed: the prefix map has no default beyond "Intensity ", so a
    quant type this function does not know is validated as though it were Intensity."""
    pg = _pg_frame()
    for experiment in ("t1_1", "t1_2", "c_1"):
        pg["Intensity {}".format(experiment)] = 1

    validate_maxquant_inputs(
        _write(tmp_path, _ed_frame()),
        _write(tmp_path, pg, name="proteinGroups.txt", sep="\t"),
        "iBAQ")


def test_the_design_is_validated_before_the_data_file(tmp_path):
    """A design problem is the more actionable of the two, and it is checked with no
    reference to the data file at all."""
    frame = _ed_frame()
    frame.loc[1, "Type"] = "X"

    with pytest.raises(EDInvalidTypeError):
        validate_maxquant_inputs(_write(tmp_path, frame),
                                 str(tmp_path / "absent.txt"), "Intensity")


def test_a_header_only_diann_matrix_is_reported_as_empty(tmp_path):
    path = tmp_path / "report.pg_matrix.tsv"
    path.write_text("Protein.Group\tt1_1\tt1_2\tc_1\n")

    with pytest.raises(PGFileError) as excinfo:
        validate_diann_inputs(_write(tmp_path, _ed_frame()), str(path))

    assert "empty" in excinfo.value.user_message.lower()


def test_a_pioneer_table_without_a_protein_column_is_rejected(tmp_path):
    path = tmp_path / "protein_groups_wide.tsv"
    path.write_text("gene_names\tt1_1\tt1_2\tc_1\nG1\t1\t2\t3\n")

    with pytest.raises(PGFileError) as excinfo:
        validate_pioneer_inputs(_write(tmp_path, _ed_frame()), str(path))

    assert "protein" in excinfo.value.user_message


def test_pioneer_validation_returns_both_frames(tmp_path):
    path = tmp_path / "protein_groups_wide.tsv"
    path.write_text("protein\tt1_1\tt1_2\tc_1\nP1\t1\t2\t3\n")

    ed_df, pioneer_df = validate_pioneer_inputs(_write(tmp_path, _ed_frame()), str(path))

    assert list(pioneer_df["protein"]) == ["P1"]


def test_fragpipe_validation_returns_both_frames(tmp_path):
    fp = pd.DataFrame([{c: "x" for c in FP_COLUMNS}])
    for experiment in ("t1_1", "t1_2", "c_1"):
        fp["{} Intensity".format(experiment)] = 1

    ed_df, fp_df = validate_fragpipe_inputs(
        _write(tmp_path, _ed_frame()),
        _write(tmp_path, fp, name="combined_protein.tsv", sep="\t"))

    assert len(ed_df) == 3
    assert "Protein ID" in fp_df.columns


def test_a_decoy_column_satisfies_the_reverse_requirement():
    """MaxQuant 2.4 and later name the column Decoy."""
    PGValidator.validate_required_columns(_pg_frame().rename(columns={"Reverse": "Decoy"}))


def test_a_file_with_neither_reverse_nor_decoy_names_reverse():
    with pytest.raises(PGMissingColumnError) as excinfo:
        PGValidator.validate_required_columns(_pg_frame().drop(columns=["Reverse"]))

    assert excinfo.value.missing_columns == ["Reverse"]
    assert "Decoy" in " ".join(excinfo.value.suggestions)
