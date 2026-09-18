"""Tests for the session archive: what a downloaded session contains, which uploads are
accepted as one, and how the datasets table survives a restart."""

import os
import zipfile

import pandas as pd
import pytest

import session_archive
from session_archive import (
    DATASETS_TABLE,
    SessionArchiveError,
    extract_session_archive,
    inspect_session_archive,
    load_datasets_table,
    save_datasets_table,
    write_session_archive,
)


def _row(name, scored="Yes"):
    return {"Dataset Name": name, "Input Type": "MaxQuant", "Quant Type": "Intensity",
            "Experiments": 4, "Controls": 2, "Scored": scored, "Imputation": 0,
            "WDFDR iterations": 2}


def _session(out_dir, names):
    """Write a table plus one file per dataset directory, as a session would leave them."""
    for name in names:
        (out_dir / name).mkdir()
        (out_dir / name / "merged.csv").write_text("x")
    save_datasets_table(str(out_dir), pd.DataFrame([_row(n) for n in names]))


# --- the persisted table --------------------------------------------------------------

def test_no_table_means_an_empty_session(tmp_path):
    table = load_datasets_table(str(tmp_path))

    assert table.empty
    assert list(table.columns) == session_archive.DATASET_COLUMNS


def test_the_table_round_trips(tmp_path):
    _session(tmp_path, ["A", "B"])

    assert list(load_datasets_table(str(tmp_path))["Dataset Name"]) == ["A", "B"]


def test_a_row_without_its_directory_is_dropped(tmp_path, caplog):
    """Every tab reads results from the directory; a row without one only breaks them."""
    _session(tmp_path, ["A", "B"])
    (tmp_path / "B" / "merged.csv").unlink()
    (tmp_path / "B").rmdir()

    with caplog.at_level("WARNING"):
        table = load_datasets_table(str(tmp_path))

    assert list(table["Dataset Name"]) == ["A"]
    assert "B" in caplog.text


# --- the download -----------------------------------------------------------------------

def test_the_archive_holds_the_table_and_the_dataset_directories_only(tmp_path):
    out = tmp_path / "Outputs"
    out.mkdir()
    _session(out, ["A"])
    (out / "pca_plot.png").write_text("png")
    (out / "proximate-server.log").write_text("log")
    (out / "A_all_20260101.zip").write_text("zip")
    zip_path = tmp_path / "session.zip"

    write_session_archive(str(out), ["A"], str(zip_path))

    assert set(zipfile.ZipFile(zip_path).namelist()) == {DATASETS_TABLE, "A/merged.csv"}


# --- the upload -------------------------------------------------------------------------

def _archive(tmp_path, members, prefix=""):
    zip_path = tmp_path / "upload.zip"
    with zipfile.ZipFile(zip_path, "w") as zipf:
        for name, content in members.items():
            zipf.writestr(prefix + name, content)
    return str(zip_path)


def _table_csv(names):
    return pd.DataFrame([_row(n) for n in names]).to_csv(index=False)


def test_an_archive_as_downloaded_has_no_prefix(tmp_path):
    path = _archive(tmp_path, {DATASETS_TABLE: _table_csv(["A"]), "A/merged.csv": "x"})

    assert inspect_session_archive(path) == ""


def test_a_rezipped_folder_is_accepted_and_its_folder_stripped(tmp_path):
    """Unpacking a download and zipping the folder back up puts everything one level
    down; that is the layout Windows and macOS produce by default."""
    path = _archive(tmp_path, {DATASETS_TABLE: _table_csv(["A"]), "A/merged.csv": "x"},
                    prefix="ProxiMateSession_20260101/")
    out = tmp_path / "Outputs"
    out.mkdir()

    assert inspect_session_archive(path) == "ProxiMateSession_20260101/"
    table = extract_session_archive(path, str(out))
    assert list(table["Dataset Name"]) == ["A"]
    assert (out / "A" / "merged.csv").read_text() == "x"
    assert not (out / "ProxiMateSession_20260101").exists()


def test_a_per_dataset_results_zip_is_refused(tmp_path):
    path = _archive(tmp_path, {"merged.csv": "x", "annotated_scores.csv": "y", "run.json": "{}"})

    with pytest.raises(SessionArchiveError, match="Not a ProxiMate session archive"):
        inspect_session_archive(path)


def test_a_file_that_is_not_a_zip_is_refused(tmp_path):
    path = tmp_path / "notes.txt"
    path.write_text("hello")

    with pytest.raises(SessionArchiveError, match="not a zip"):
        inspect_session_archive(str(path))


def test_a_member_escaping_the_output_directory_is_refused(tmp_path):
    path = _archive(tmp_path, {DATASETS_TABLE: _table_csv([]), "../etc/passwd": "x"})

    with pytest.raises(SessionArchiveError, match="unsafe path"):
        inspect_session_archive(path)


def test_extraction_drops_rows_whose_directory_is_not_in_the_archive(tmp_path):
    path = _archive(tmp_path, {DATASETS_TABLE: _table_csv(["A", "B"]), "A/merged.csv": "x"})
    out = tmp_path / "Outputs"
    out.mkdir()

    table = extract_session_archive(path, str(out))

    assert list(table["Dataset Name"]) == ["A"]
