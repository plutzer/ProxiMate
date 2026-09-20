"""Tests for the session archive: what a downloaded session contains, which uploads are
accepted as one, and how the datasets table survives a restart."""

import zipfile

import pandas as pd
import pytest

from session_archive import (
    DATASETS_TABLE,
    SessionArchiveError,
    extract_session_archive,
    inspect_session_archive,
    load_datasets_table,
    save_datasets_table,
    write_session_archive,
)


def _row(name):
    return {"Dataset Name": name, "Input Type": "MaxQuant", "Quant Type": "Intensity",
            "Experiments": 4, "Controls": 2, "Scored": "Yes", "Imputation": 0,
            "WDFDR iterations": 2}


def _session(out_dir, names):
    """Write a table plus one file per dataset directory, as a session would leave them."""
    for name in names:
        (out_dir / name).mkdir()
        (out_dir / name / "merged.csv").write_text("x")
    save_datasets_table(str(out_dir), pd.DataFrame([_row(n) for n in names]))


def _archive(tmp_path, members, prefix=""):
    zip_path = tmp_path / "upload.zip"
    with zipfile.ZipFile(zip_path, "w") as zipf:
        for name, content in members.items():
            zipf.writestr(prefix + name, content)
    return str(zip_path)


def _table_csv(names):
    return pd.DataFrame([_row(n) for n in names]).to_csv(index=False)


def test_a_row_without_its_directory_is_dropped_on_load(tmp_path):
    """Every tab reads results from the directory; a row without one only breaks them."""
    _session(tmp_path, ["A", "B"])
    (tmp_path / "B" / "merged.csv").unlink()
    (tmp_path / "B").rmdir()

    assert list(load_datasets_table(str(tmp_path))["Dataset Name"]) == ["A"]


def test_the_archive_holds_the_table_and_the_dataset_directories_only(tmp_path):
    out = tmp_path / "Outputs"
    out.mkdir()
    _session(out, ["A"])
    (out / "A" / "plots").mkdir()
    (out / "A" / "plots" / "pca.png").write_text("png")
    (out / "pca_plot.png").write_text("png")
    (out / "proximate-server.log").write_text("log")
    (out / "A_all_20260101.zip").write_text("zip")
    zip_path = tmp_path / "session.zip"

    write_session_archive(str(out), ["A"], str(zip_path))

    assert set(zipfile.ZipFile(zip_path).namelist()) == {
        DATASETS_TABLE, "A/merged.csv", "A/plots/pca.png"}


@pytest.mark.parametrize("prefix", ["", "ProxiMateSession_20260101/"])
def test_an_archive_is_accepted_with_or_without_a_top_level_folder(tmp_path, prefix):
    """Unpacking a download and zipping the folder back up puts everything one level
    down; that is the layout Windows and macOS produce by default."""
    path = _archive(tmp_path, {DATASETS_TABLE: _table_csv(["A"]), "A/merged.csv": "x"},
                    prefix=prefix)
    out = tmp_path / "Outputs"
    out.mkdir()

    assert inspect_session_archive(path) == prefix
    table = extract_session_archive(path, str(out))
    assert list(table["Dataset Name"]) == ["A"]
    assert (out / "A" / "merged.csv").read_text() == "x"
    if prefix:
        assert not (out / prefix.rstrip("/")).exists()


@pytest.mark.parametrize("members", [
    {"merged.csv": "x", "annotated_scores.csv": "y", "run.json": "{}"},
    None,
    {DATASETS_TABLE: "", "../etc/passwd": "x"},
], ids=["per-dataset results zip", "not a zip", "path traversal"])
def test_unusable_uploads_are_refused(tmp_path, members):
    if members is None:
        path = tmp_path / "notes.txt"
        path.write_text("hello")
        path = str(path)
    else:
        path = _archive(tmp_path, members)

    with pytest.raises(SessionArchiveError):
        inspect_session_archive(path)


def test_extraction_drops_rows_whose_directory_is_not_in_the_archive(tmp_path):
    path = _archive(tmp_path, {DATASETS_TABLE: _table_csv(["A", "B"]), "A/merged.csv": "x"})
    out = tmp_path / "Outputs"
    out.mkdir()

    table = extract_session_archive(path, str(out))

    assert list(table["Dataset Name"]) == ["A"]
