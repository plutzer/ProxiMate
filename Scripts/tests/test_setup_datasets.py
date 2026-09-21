"""Tests for the dataset preparation helpers.

Dataset preparation runs once when an image is built, and its failures are quiet: a
download that returns a login page instead of an archive, or a source that renames a
column, both produce a dataset that loads and annotates nothing.  The verification steps
here are what turn those into a failed build.

Everything is exercised offline.  An autouse fixture makes ``requests.get`` raise, so a
test that reached the network would fail rather than depend on it.
"""

import io
import zipfile

import pytest

import setup_datasets


@pytest.fixture(autouse=True)
def no_network(monkeypatch):
    """Nothing in this module may reach the network."""
    def _forbidden(*args, **kwargs):
        raise AssertionError("a test attempted a network request")

    monkeypatch.setattr(setup_datasets.requests, "get", _forbidden)


def _zip_bytes(members):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)
    return buffer.getvalue()


# --- verify_tsv_columns --------------------------------------------------------------

@pytest.mark.parametrize("header, encoding, expected", [
    (["Gene name", "Main location", "extra"], "utf-8", True),
    (["Gene name", "Main_location"], "utf-8", False),
    ([], "utf-8", False),
    (["Gene name", "Main location", "caf\xe9"], "latin-1", True),
], ids=["required present", "renamed column", "empty file", "latin-1 retried"])
def test_verify_tsv_columns(tmp_path, header, encoding, expected):
    """A renamed column is the failure the check exists for: the download succeeds and
    the file looks fine, but the join it feeds would match nothing.  CORUM ships text
    that is not valid UTF-8, so a decoding failure is retried rather than rejected."""
    path = tmp_path / "data.tsv"
    path.write_bytes(("\t".join(header) + "\nrow\n").encode(encoding) if header else b"")

    assert setup_datasets.verify_tsv_columns(
        str(path), setup_datasets.HPA_REQUIRED_COLUMNS) is expected


def test_a_missing_file_fails_verification(tmp_path):
    assert setup_datasets.verify_tsv_columns(
        str(tmp_path / "absent"), setup_datasets.HPA_REQUIRED_COLUMNS) is False


# --- extract_from_zip -----------------------------------------------------------------

def test_the_matching_member_is_extracted_whole(tmp_path):
    """The member is copied in chunks into a directory that may not exist yet, so a
    payload larger than one chunk is where a truncating copy would show up."""
    payload = "x" * (setup_datasets.CHUNK_SIZE * 3 + 17)
    target = tmp_path / "out" / "subcellular_location.tsv"
    content = _zip_bytes({"readme.md": "no", "subcellular_location.tsv": payload})

    assert setup_datasets.extract_from_zip(content, ".tsv", str(target), "HPA") is True
    assert target.read_text() == payload


@pytest.mark.parametrize("content", [
    b"<html><body>Please accept the licence</body>",
    b"<!DOCTYPE html><body>Please accept the licence</body>",
    b"not a zip at all",
    _zip_bytes({"readme.md": "nothing useful"}),
], ids=["html", "doctype html", "corrupt", "no matching member"])
def test_unusable_downloads_are_rejected_without_writing(tmp_path, content):
    """BioGRID answers an unaccepted licence with a web page.  Saving it would leave a
    file that exists and is non-empty, so every later check would pass."""
    target = tmp_path / "out.tsv"

    assert setup_datasets.extract_from_zip(content, ".tsv", str(target), "BG") is False
    assert not target.exists()


# --- run_preprocess_biogrid -------------------------------------------------------------

@pytest.mark.parametrize("exclude_hcm, summary, expected_args", [
    (False, setup_datasets.BIOGRID_SUMMARY_FILENAME, []),
    (True, setup_datasets.BIOGRID_NO_HCM_SUMMARY_FILENAME,
     ["--exclude_publication", setup_datasets.HCM_PUBLICATION,
      "--output_filename", setup_datasets.BIOGRID_NO_HCM_SUMMARY_FILENAME]),
], ids=["full", "no HCM"])
def test_preprocessing_passes_the_organism_and_variant_to_the_subprocess(
        tmp_path, monkeypatch, exclude_hcm, summary, expected_args):
    """The taxonomy id has to reach the subprocess: preprocessing filters on it, and a
    wrong one yields an empty summary rather than an error.  The no-HCM variant must
    land beside the full summary without replacing it."""
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
        (organism_dir / summary).write_text("x")
        return _Completed()

    monkeypatch.setattr(setup_datasets.subprocess, "run", _run)

    assert setup_datasets.run_preprocess_biogrid(
        str(tmp_path), "human", 10090, exclude_hcm=exclude_hcm) is True
    cmd = calls[0]
    assert cmd[cmd.index("--organism_id") + 1] == "10090"
    for arg in expected_args:
        assert arg in cmd
