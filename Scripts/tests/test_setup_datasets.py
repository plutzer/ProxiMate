"""Tests for the dataset preparation helpers.

Dataset preparation runs once when an image is built, and its failures are quiet: a
download that returns a login page instead of an archive, or a source that renames a
column, both produce a dataset that loads and annotates nothing.  The verification steps
here are what turn those into a failed build.

Everything is exercised offline.  An autouse fixture makes ``requests.get`` raise, so a
test that reached the network would fail rather than depend on it.
"""

import io
import os
import zipfile

import pytest

import annotator
import preprocess_biogrid
import setup_datasets


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Nothing in this module may reach the network."""
    def _forbidden(*args, **kwargs):
        raise AssertionError("a test attempted a network request")

    monkeypatch.setattr(setup_datasets.requests, "get", _forbidden)


# --- cross-module configuration --------------------------------------------------

def test_the_organism_configuration_matches_the_annotator():
    """Both modules carry their own copy and a comment asking that they be kept in step.
    A download prepared for one organism and annotation expecting another produces empty
    annotations rather than an error."""
    assert setup_datasets.ORGANISMS == annotator.ORGANISMS


def test_biogrid_verification_covers_the_columns_preprocessing_reads():
    """Verification passing and preprocessing then failing on a column it never checked
    would move the failure from build time to a raw KeyError."""
    source = open(preprocess_biogrid.__file__, encoding="utf-8").read()
    used = {c for c in (
        "Organism ID Interactor A", "Organism ID Interactor B",
        "Experimental System Type", "Experimental System",
        "SWISS-PROT Accessions Interactor A", "SWISS-PROT Accessions Interactor B",
        "Author", "Publication Source") if '"{}"'.format(c) in source
        or "'{}'".format(c) in source}

    unverified = used - setup_datasets.BIOGRID_REQUIRED_COLUMNS

    assert unverified == {"Experimental System", "Author", "Publication Source"}, (
        "the set of BioGRID columns used but not verified has changed")


# --- file_exists_and_nonempty ------------------------------------------------------

def test_a_written_file_is_recognized(tmp_path):
    path = tmp_path / "data.tsv"
    path.write_text("content")

    assert setup_datasets.file_exists_and_nonempty(str(path)) is True


def test_a_missing_file_is_not_recognized(tmp_path):
    assert setup_datasets.file_exists_and_nonempty(str(tmp_path / "absent")) is False


def test_an_empty_file_is_not_recognized(tmp_path):
    """A download interrupted at the first byte leaves one of these behind."""
    path = tmp_path / "data.tsv"
    path.write_text("")

    assert setup_datasets.file_exists_and_nonempty(str(path)) is False


# --- verify_tsv_columns --------------------------------------------------------------

def _write_tsv(tmp_path, columns, encoding="utf-8", name="data.tsv", sep="\t"):
    path = tmp_path / name
    path.write_bytes((sep.join(columns) + "\nrow\n").encode(encoding))
    return str(path)


def test_a_file_with_the_required_columns_verifies(tmp_path):
    path = _write_tsv(tmp_path, ["Gene name", "Main location", "extra"])

    assert setup_datasets.verify_tsv_columns(
        path, setup_datasets.HPA_REQUIRED_COLUMNS) is True


def test_a_renamed_column_fails_verification(tmp_path):
    """This is the failure the check exists for: the download succeeds and the file looks
    fine, but the join it feeds would match nothing."""
    path = _write_tsv(tmp_path, ["Gene name", "Main_location"])

    assert setup_datasets.verify_tsv_columns(
        path, setup_datasets.HPA_REQUIRED_COLUMNS) is False


def test_the_missing_columns_are_named(tmp_path, capsys):
    """setup_datasets writes to stderr directly rather than through the package logger."""
    path = _write_tsv(tmp_path, ["Gene name"])

    setup_datasets.verify_tsv_columns(path, setup_datasets.HPA_REQUIRED_COLUMNS)

    assert "Main location" in capsys.readouterr().err


def test_an_empty_file_fails_verification(tmp_path):
    path = tmp_path / "data.tsv"
    path.write_text("")

    assert setup_datasets.verify_tsv_columns(
        path, setup_datasets.HPA_REQUIRED_COLUMNS) is False


def test_a_missing_file_fails_verification(tmp_path):
    assert setup_datasets.verify_tsv_columns(
        str(tmp_path / "absent"), setup_datasets.HPA_REQUIRED_COLUMNS) is False


def test_a_latin_1_file_is_retried_rather_than_rejected(tmp_path):
    """CORUM ships text that is not valid UTF-8; rejecting it would fail the build over
    an encoding rather than a missing column."""
    path = _write_tsv(tmp_path, ["complex_name", "subunits_uniprot_id", "caf\xe9"],
                      encoding="latin-1")

    assert setup_datasets.verify_tsv_columns(
        path, setup_datasets.CORUM_REQUIRED_COLUMNS) is True


def test_the_separator_is_configurable(tmp_path):
    path = _write_tsv(tmp_path, ["Gene name", "Main location"], sep=",")

    assert setup_datasets.verify_tsv_columns(
        path, setup_datasets.HPA_REQUIRED_COLUMNS, sep=",") is True


def test_only_the_header_is_inspected(tmp_path):
    """Verification reads one line, so a truncated or corrupt body still verifies."""
    path = tmp_path / "data.tsv"
    path.write_text("Gene name\tMain location\n\x00 garbage")

    assert setup_datasets.verify_tsv_columns(
        str(path), setup_datasets.HPA_REQUIRED_COLUMNS) is True


# --- extract_from_zip -----------------------------------------------------------------

def _zip_bytes(members):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    return buffer.getvalue()


def test_a_matching_member_is_extracted(tmp_path):
    target = tmp_path / "out" / "subcellular_location.tsv"
    content = _zip_bytes({"subcellular_location.tsv": "Gene name\tMain location\n"})

    assert setup_datasets.extract_from_zip(content, ".tsv", str(target), "HPA") is True
    assert target.read_text() == "Gene name\tMain location\n"


def test_the_target_directory_is_created(tmp_path):
    target = tmp_path / "nested" / "deeper" / "out.tsv"

    setup_datasets.extract_from_zip(_zip_bytes({"a.tsv": "x"}), ".tsv", str(target), "L")

    assert target.exists()


def test_content_larger_than_one_chunk_is_copied_whole(tmp_path):
    """The member is copied in chunks, so a file larger than one chunk is where a
    truncating copy would show up."""
    payload = "x" * (setup_datasets.CHUNK_SIZE * 3 + 17)
    target = tmp_path / "out.tsv"

    setup_datasets.extract_from_zip(_zip_bytes({"a.tsv": payload}), ".tsv",
                                    str(target), "L")

    assert target.read_text() == payload


@pytest.mark.parametrize("prefix", [b"<html", b"<!DOC", b"<HTML", b"<!doc"])
def test_an_html_response_is_rejected(tmp_path, prefix, capsys):
    """BioGRID answers an unaccepted licence with a web page.  Saving it would leave a
    file that exists and is non-empty, so every later check would pass."""
    target = tmp_path / "out.tsv"

    result = setup_datasets.extract_from_zip(
        prefix + b"><body>Please accept the licence</body>", ".tsv", str(target), "BG")

    assert result is False
    assert not target.exists()
    assert "HTML instead of ZIP" in capsys.readouterr().err


def test_a_corrupt_archive_is_rejected(tmp_path):
    target = tmp_path / "out.tsv"

    assert setup_datasets.extract_from_zip(b"not a zip at all", ".tsv",
                                           str(target), "L") is False


def test_an_archive_without_the_expected_member_is_rejected(tmp_path, capsys):
    content = _zip_bytes({"readme.md": "nothing useful"})
    target = tmp_path / "out.tsv"

    assert setup_datasets.extract_from_zip(content, ".tsv", str(target), "L") is False
    assert "readme.md" in capsys.readouterr().err


def test_the_first_matching_member_is_taken(tmp_path):
    content = _zip_bytes({"a.tsv": "first", "b.tsv": "second"})
    target = tmp_path / "out.tsv"

    setup_datasets.extract_from_zip(content, ".tsv", str(target), "L")

    assert target.read_text() == "first"


# --- write_build_info -------------------------------------------------------------------

def test_the_build_info_records_each_result(tmp_path):
    setup_datasets.write_build_info(str(tmp_path), {"uniprot": True, "biogrid": False})
    written = (tmp_path / setup_datasets.BUILD_INFO_FILENAME).read_text()

    assert "uniprot: OK" in written
    assert "biogrid: FAILED" in written


def test_the_build_info_names_the_organisms(tmp_path):
    """It is copied into each dataset's output directory, so it is the record of which
    organisms a scored run could have annotated against."""
    setup_datasets.write_build_info(str(tmp_path), {})
    written = (tmp_path / setup_datasets.BUILD_INFO_FILENAME).read_text()

    for organism in setup_datasets.ORGANISMS:
        assert organism in written


def test_the_build_info_is_dated(tmp_path):
    setup_datasets.write_build_info(str(tmp_path), {})
    written = (tmp_path / setup_datasets.BUILD_INFO_FILENAME).read_text()

    assert "Build date:" in written


# --- run_preprocess_biogrid -------------------------------------------------------------

def test_preprocessing_is_skipped_when_its_inputs_are_absent(tmp_path):
    """No subprocess is started, so this needs no stubbing at all."""
    assert setup_datasets.run_preprocess_biogrid(str(tmp_path), "human", 9606) is False


def test_preprocessing_reports_the_organism_it_was_given(tmp_path, monkeypatch):
    """The taxonomy id has to reach the subprocess: preprocessing filters on it, and a
    wrong one yields an empty summary rather than an error."""
    # The two downloads are shared across organisms; only the summary is per-organism.
    for filename in (setup_datasets.BIOGRID_ALL_FILENAME,
                     setup_datasets.BIOGRID_MV_FILENAME):
        (tmp_path / filename).write_text("content")
    organism_dir = tmp_path / "human"

    calls = []

    class _Completed:
        returncode = 0
        stdout = ""
        stderr = ""

    def _run(cmd, **kwargs):
        calls.append(cmd)
        (organism_dir / setup_datasets.BIOGRID_SUMMARY_FILENAME).write_text("x")
        return _Completed()

    monkeypatch.setattr(setup_datasets.subprocess, "run", _run)

    assert setup_datasets.run_preprocess_biogrid(str(tmp_path), "human", 10090) is True
    assert "--organism_id" in calls[0]
    assert "10090" in calls[0]
