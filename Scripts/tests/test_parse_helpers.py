"""Tests for the standalone helpers in parse.py.

None of these touch proteomics data.  ``validate_name`` guards the dataset names the GUI
turns into directories, ``_read_msstats_csv`` is the encoding fallback every MSstats run
goes through, and ``_parse_stage`` decides what the run manifest records about a parse --
including a parse that failed, which is the manifest most worth reading.
"""

import json

import pandas as pd
import pytest

import parse


# --- validate_name -------------------------------------------------------------

def test_a_valid_name_returns_zero():
    """Success is the int 0, not None and not True.  It is falsy, so callers must
    compare against 0 rather than test truthiness."""
    result = parse.validate_name("dataset_01", [])

    assert result == 0
    assert isinstance(result, int)


@pytest.mark.parametrize("name", ["", None])
def test_an_absent_name_is_rejected(name):
    """The falsy check also keeps None away from re.search, which would raise."""
    assert "empty" in parse.validate_name(name, [])


@pytest.mark.parametrize("name", ["has space", "has\ttab", " leading"])
def test_whitespace_is_rejected(name):
    assert "spaces" in parse.validate_name(name, [])


@pytest.mark.parametrize("name", ["has-hyphen", "has.dot", "has/slash", "café"])
def test_non_alphanumeric_characters_are_rejected(name):
    """The character class is ASCII-only, so an accented letter is rejected even though
    it is a letter."""
    assert "letters, numbers, and underscores" in parse.validate_name(name, [])


def test_the_whitespace_check_precedes_the_character_check():
    """A name failing both reports the space, which is the more actionable of the two."""
    assert "spaces" in parse.validate_name("a -b", [])


@pytest.mark.parametrize("name", ["_", "123", "a_1", "A1"])
def test_leading_digits_and_underscores_are_accepted(name):
    assert parse.validate_name(name, []) == 0


def test_an_existing_name_is_rejected():
    message = parse.validate_name("existing", ["existing", "other"])

    assert "already exists" in message


def test_the_duplicate_check_is_case_sensitive():
    assert parse.validate_name("Existing", ["existing"]) == 0


# --- _read_msstats_csv ---------------------------------------------------------

MSSTATS_TEXT = "Protein,originalRUN,LogIntensities\nCafé_HUMAN,run_01,8.5\n"


@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig", "latin-1", "cp1252"])
def test_every_supported_encoding_round_trips(tmp_path, encoding):
    path = tmp_path / "ProteinLevelData.csv"
    path.write_bytes(MSSTATS_TEXT.encode(encoding))

    frame = parse._read_msstats_csv(str(path))

    assert list(frame.columns) == ["Protein", "originalRUN", "LogIntensities"]
    assert len(frame) == 1


def test_a_latin_1_file_is_recovered_after_utf_8_fails(tmp_path):
    """latin-1 decodes any byte sequence, so it is what actually terminates the chain;
    cp1252 is never reached through the decode path."""
    path = tmp_path / "ProteinLevelData.csv"
    path.write_bytes(MSSTATS_TEXT.encode("latin-1"))

    assert parse._read_msstats_csv(str(path))["Protein"].iloc[0] == "Café_HUMAN"


def test_a_header_only_file_exhausts_the_chain(tmp_path):
    """Readable under every encoding but empty, so no iteration returns and the loop
    falls through to its own error rather than reporting a decoding failure."""
    path = tmp_path / "ProteinLevelData.csv"
    path.write_text("Protein,originalRUN\n")

    with pytest.raises(ValueError, match="Could not read MSstats file"):
        parse._read_msstats_csv(str(path))


def test_a_truly_empty_file_raises_from_pandas(tmp_path):
    """Only UnicodeDecodeError is caught, so an EmptyDataError propagates unwrapped."""
    path = tmp_path / "ProteinLevelData.csv"
    path.write_text("")

    with pytest.raises(pd.errors.EmptyDataError):
        parse._read_msstats_csv(str(path))


# --- _parse_stage --------------------------------------------------------------

@parse._parse_stage
def _records_two_counts(proteinGroups, quantType, outputPath):
    return 3, 2


@parse._parse_stage
def _takes_a_frame(bait_df, preyfile, outputPath):
    return 1, 0


@parse._parse_stage
def _fails(proteinGroups, outputPath):
    raise RuntimeError("parsing blew up")


@parse._parse_stage
def _returns_one_value(proteinGroups, outputPath):
    return 7


def _stage_entry(out_dir):
    """The single stage recorded in the manifest at `out_dir`."""
    document = json.loads((out_dir / "run.json").read_text())
    runs = list(document["runs"].values())
    assert len(runs) == 1
    assert len(runs[0]["stages"]) == 1
    return runs[0]["stages"][0]


@pytest.fixture
def input_file(tmp_path):
    path = tmp_path / "proteinGroups.txt"
    path.write_text("Majority protein IDs\nP1\n")
    return path


def test_the_stage_is_named_for_the_wrapped_function(tmp_path, input_file):
    out = tmp_path / "out"

    _records_two_counts(str(input_file), "LFQ", str(out))

    assert _stage_entry(out)["entrypoint"] == "parse._records_two_counts"


def test_the_returned_counts_become_metrics(tmp_path, input_file):
    out = tmp_path / "out"

    assert _records_two_counts(str(input_file), "LFQ", str(out)) == (3, 2)
    assert _stage_entry(out)["metrics"] == {"n_experiments": 3, "n_controls": 2}


def test_file_parameters_are_recorded_as_inputs(tmp_path, input_file):
    """A checksum of the file is what makes a run reproducible; a value parameter has
    nothing to check."""
    out = tmp_path / "out"

    _records_two_counts(str(input_file), "LFQ", str(out))
    entry = _stage_entry(out)

    assert [i["role"] for i in entry["inputs"]] == ["proteinGroups"]
    assert entry["params"]["quantType"] == "LFQ"


def test_the_output_path_is_not_recorded_as_a_parameter(tmp_path, input_file):
    out = tmp_path / "out"

    _records_two_counts(str(input_file), "LFQ", str(out))

    assert "outputPath" not in _stage_entry(out)["params"]


def test_arguments_bind_by_keyword_as_well_as_by_position(tmp_path, input_file):
    out = tmp_path / "out"

    _records_two_counts(proteinGroups=str(input_file), quantType="LFQ",
                        outputPath=str(out))

    assert _stage_entry(out)["params"]["quantType"] == "LFQ"


def test_a_dataframe_argument_is_recorded_by_shape(tmp_path):
    """parse_from_saint is handed a frame rather than a path.  Its contents do not
    belong in a manifest, but its shape says whether the right thing arrived."""
    out = tmp_path / "out"
    prey = tmp_path / "prey.txt"
    prey.write_text("P1\tG1\n")
    frame = pd.DataFrame({"Experiment Name": ["a", "b"], "Bait": ["X", "Y"]})

    _takes_a_frame(frame, str(prey), str(out))
    entry = _stage_entry(out)

    assert entry["params"]["bait_df"] == "DataFrame(2, 2)"
    assert [i["role"] for i in entry["inputs"]] == ["preyfile"]


def test_an_absent_input_file_is_recorded_as_missing(tmp_path):
    """A path that does not exist is fingerprinted as missing rather than omitted --
    "the input was not there" is itself worth knowing when reading a manifest later."""
    out = tmp_path / "out"
    absent = tmp_path / "never_written.txt"

    _records_two_counts(str(absent), "LFQ", str(out))
    recorded = _stage_entry(out)["inputs"][0]

    assert recorded["missing"] is True
    assert "sha256" not in recorded


def test_a_none_path_fails_before_the_parse_runs(tmp_path):
    """Documented, not fixed: an absent path is tolerated but a None one is not.
    file_record recovers from OSError, which os.path.getsize raises for a path that is
    merely missing, and not from the TypeError it raises for None.  The parse itself
    never starts, and the manifest records the stage as an error."""
    out = tmp_path / "out"

    with pytest.raises(TypeError):
        _records_two_counts(None, "LFQ", str(out))

    assert _stage_entry(out)["status"] == "error"


def test_outputs_are_recorded_only_when_they_exist(tmp_path, input_file):
    out = tmp_path / "out"
    out.mkdir()
    (out / "prey.txt").write_text("P1\tG1\n")

    _records_two_counts(str(input_file), "LFQ", str(out))
    written = {o["path"].replace("\\", "/").rsplit("/", 1)[-1]
               for o in _stage_entry(out)["outputs"]}

    assert written == {"prey.txt"}


def test_a_failing_parse_still_leaves_a_manifest(tmp_path, input_file):
    """The run that failed is the one whose manifest is worth reading."""
    out = tmp_path / "out"

    with pytest.raises(RuntimeError, match="parsing blew up"):
        _fails(str(input_file), str(out))

    entry = _stage_entry(out)
    assert entry["status"] == "error"
    assert entry["error"]["type"] == "RuntimeError"


def test_the_exception_is_not_suppressed(tmp_path, input_file):
    out = tmp_path / "out"

    with pytest.raises(RuntimeError):
        _fails(str(input_file), str(out))


def test_an_entry_point_must_return_two_counts(tmp_path, input_file):
    """Every parse entry point reports (n_experiments, n_controls); the manifest and the
    GUI both read that pair."""
    out = tmp_path / "out"

    with pytest.raises(TypeError):
        _returns_one_value(str(input_file), str(out))
