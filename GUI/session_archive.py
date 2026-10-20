"""The session archive: one zip holding every dataset directory plus the datasets table.

The table (``datasets.csv`` at the top of the output directory) is what the GUI shows
as "Datasets in this Session"; it is rewritten whenever the session changes and read
back at server start, so a session outlives the server process.
"""

import json
import os
import zipfile

import pandas as pd

import provenance
from log_config import get_logger

logger = get_logger(__name__)

DATASETS_TABLE = "datasets.csv"
DATASET_COLUMNS = ['Dataset Name', 'Input Type', 'Quant Type', 'Experiments', 'Controls',
                   'Scored', 'Scored at', 'Organism', 'Imputation', 'WDFDR iterations']
# Columns whose values the dataset's run.json also holds; see ``scored_fields``.
MANIFEST_COLUMNS = ('Scored at', 'Organism')


class SessionArchiveError(ValueError):
    """The uploaded file is not a usable session archive; the message says why."""


def empty_datasets_table():
    return pd.DataFrame(columns=DATASET_COLUMNS)


def save_datasets_table(out_dir, table):
    table.to_csv(os.path.join(out_dir, DATASETS_TABLE), index=False)


def scored_fields(dataset_dir):
    """'Scored at' and 'Organism' as the dataset's run.json records them: when its last
    successful annotate stage ended and the organism it annotated against.  Blank for
    a dataset without a finished annotation."""
    path = os.path.join(dataset_dir, provenance.RUN_JSON_FILENAME)
    if not os.path.isfile(path):
        return {'Scored at': '', 'Organism': ''}
    with open(path, encoding='utf-8') as handle:
        document = json.load(handle)
    ended = [stage['ended_utc'] for run in document['runs'].values()
             for stage in run.get('stages', [])
             if stage.get('stage') == 'annotate' and stage.get('status') == 'ok']
    if not ended:
        return {'Scored at': '', 'Organism': ''}
    return {'Scored at': ended[-1][:16].replace('T', ' ') + ' UTC',
            'Organism': provenance.annotation_settings(dataset_dir)['organism']}


def load_datasets_table(out_dir):
    """The persisted table, less any row whose dataset directory is gone.

    No file means an empty session.  Rows are dropped rather than kept when their
    directory is missing because every tab reads the results from that directory.
    A table written without the manifest columns gets them filled from each
    dataset's run.json.
    """
    path = os.path.join(out_dir, DATASETS_TABLE)
    if not os.path.isfile(path):
        return empty_datasets_table()
    table = pd.read_csv(path)
    present = table['Dataset Name'].map(lambda name: os.path.isdir(os.path.join(out_dir, str(name))))
    if (~present).any():
        logger.warning("Dropping %d dataset(s) from the session table with no directory: %s",
                       (~present).sum(), ", ".join(table.loc[~present, 'Dataset Name'].astype(str)))
    table = table[present].reset_index(drop=True)
    missing = [c for c in MANIFEST_COLUMNS if c not in table.columns]
    if missing:
        fields = [scored_fields(os.path.join(out_dir, str(name))) for name in table['Dataset Name']]
        for column in missing:
            table[column] = [f[column] for f in fields]
    return table[DATASET_COLUMNS]


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
