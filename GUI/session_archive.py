"""The session archive: one zip holding every dataset directory plus the datasets table.

The table (``datasets.csv`` at the top of the output directory) is what the GUI shows
as "Datasets in this Session"; it is rewritten whenever the session changes and read
back at server start, so a session outlives the server process.
"""

import os
import zipfile

import pandas as pd

from log_config import get_logger

logger = get_logger(__name__)

DATASETS_TABLE = "datasets.csv"
DATASET_COLUMNS = ['Dataset Name', 'Input Type', 'Quant Type', 'Experiments', 'Controls',
                   'Scored', 'Imputation', 'WDFDR iterations']


class SessionArchiveError(ValueError):
    """The uploaded file is not a usable session archive; the message says why."""


def empty_datasets_table():
    return pd.DataFrame(columns=DATASET_COLUMNS)


def save_datasets_table(out_dir, table):
    table.to_csv(os.path.join(out_dir, DATASETS_TABLE), index=False)


def load_datasets_table(out_dir):
    """The persisted table, less any row whose dataset directory is gone.

    No file means an empty session.  Rows are dropped rather than kept when their
    directory is missing because every tab reads the results from that directory.
    """
    path = os.path.join(out_dir, DATASETS_TABLE)
    if not os.path.isfile(path):
        return empty_datasets_table()
    table = pd.read_csv(path)
    present = table['Dataset Name'].map(lambda name: os.path.isdir(os.path.join(out_dir, str(name))))
    if (~present).any():
        logger.warning("Dropping %d dataset(s) from the session table with no directory: %s",
                       (~present).sum(), ", ".join(table.loc[~present, 'Dataset Name'].astype(str)))
    return table[present].reset_index(drop=True)


def dataset_directories(out_dir):
    """Names of the directories under the output directory, listed or not."""
    return [entry for entry in os.listdir(out_dir) if os.path.isdir(os.path.join(out_dir, entry))]


def write_session_archive(out_dir, names, zip_path):
    """Zip the datasets table and the directory of each named dataset into ``zip_path``.

    Only those go in: the output directory also accumulates exported plots, earlier
    result zips and the server log, none of which belong to a session.
    """
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        zipf.write(os.path.join(out_dir, DATASETS_TABLE), arcname=DATASETS_TABLE)
        for name in names:
            dataset_dir = os.path.join(out_dir, str(name))
            for root, _, files in os.walk(dataset_dir):
                for file in files:
                    abs_file = os.path.join(root, file)
                    zipf.write(abs_file, arcname=os.path.relpath(abs_file, out_dir))


def inspect_session_archive(zip_path):
    """Return the folder prefix under which the archive's contents sit.

    The prefix is "" for an archive as downloaded, or "folder/" for one a user unpacked
    and re-zipped with its top-level folder.  Anything else is refused before the
    current session is touched.
    """
    if not zipfile.is_zipfile(zip_path):
        raise SessionArchiveError("The uploaded file is not a zip archive.")
    with zipfile.ZipFile(zip_path) as zipf:
        members = zipf.namelist()
    for member in members:
        if member.startswith(("/", "\\")) or ".." in member.split("/"):
            raise SessionArchiveError(f"The archive contains an unsafe path: {member}")
    if DATASETS_TABLE in members:
        return ""
    nested = [m for m in members if m.count("/") == 1 and m.endswith("/" + DATASETS_TABLE)]
    if len(nested) == 1:
        return nested[0][:-len(DATASETS_TABLE)]
    raise SessionArchiveError(
        f"Not a ProxiMate session archive: no {DATASETS_TABLE} at its top level. "
        "Per-dataset result zips cannot be loaded as a session.")


def extract_session_archive(zip_path, out_dir):
    """Unpack the archive into ``out_dir``, dropping the folder prefix if it has one.

    Returns the datasets table that arrived, less rows whose directory did not.
    """
    prefix = inspect_session_archive(zip_path)
    with zipfile.ZipFile(zip_path) as zipf:
        for member in zipf.namelist():
            if not member.startswith(prefix) or member.endswith("/"):
                continue
            target = os.path.join(out_dir, member[len(prefix):])
            os.makedirs(os.path.dirname(target), exist_ok=True)
            with zipf.open(member) as src, open(target, "wb") as dst:
                dst.write(src.read())
    return load_datasets_table(out_dir)
