"""The datasets table as one process-wide store.

The table ("Datasets in this Session") is shared by every browser session and by the
MCP server, so it lives at module level rather than in a Shiny session.  Every
mutation is written to ``datasets.csv`` under the output directory and bumps
``version``, which the GUI polls; a row added by any caller appears in every open
session within a second.
"""

import threading

import pandas as pd

import session_archive
from log_config import get_logger

logger = get_logger(__name__)

STATE = {
    'lock': threading.RLock(),
    'version': 0,
    'out_dir': None,
    'table': session_archive.empty_datasets_table(),
}


def configure(out_dir):
    """Point the store at an output directory and read its persisted table."""
    with STATE['lock']:
        STATE['out_dir'] = out_dir
        STATE['table'] = session_archive.load_datasets_table(out_dir)
        STATE['version'] += 1


def version():
    return STATE['version']


def table():
    """A copy of the table; edits to it do not reach the store."""
    with STATE['lock']:
        return STATE['table'].copy()


def names():
    return table()['Dataset Name'].astype(str).tolist()


def scored_names():
    t = table()
    return t.loc[t['Scored'] == 'Yes', 'Dataset Name'].astype(str).tolist()


def row(name):
    t = table()
    hit = t[t['Dataset Name'].astype(str) == str(name)]
    if hit.empty:
        raise KeyError(f"no dataset named {name!r} in the session")
    return hit.iloc[0].to_dict()


def _commit(new_table, what):
    STATE['table'] = new_table.reset_index(drop=True)
    session_archive.save_datasets_table(STATE['out_dir'], STATE['table'])
    STATE['version'] += 1
    logger.info("datasets table: %s (%d rows)", what, len(STATE['table']))


def append(fields):
    """Add one dataset row (a dict keyed by ``session_archive.DATASET_COLUMNS``)."""
    with STATE['lock']:
        if str(fields['Dataset Name']) in names():
            raise ValueError(f"dataset {fields['Dataset Name']!r} is already in the session")
        new_row = pd.DataFrame([[fields.get(c, '') for c in session_archive.DATASET_COLUMNS]],
                               columns=session_archive.DATASET_COLUMNS)
        _commit(pd.concat([STATE['table'], new_row], ignore_index=True),
                f"added {fields['Dataset Name']}")


def update(name, **fields):
    """Set columns on the row named ``name``."""
    with STATE['lock']:
        t = STATE['table'].copy()
        mask = t['Dataset Name'].astype(str) == str(name)
        if not mask.any():
            raise KeyError(f"no dataset named {name!r} in the session")
        for column, value in fields.items():
            t.loc[mask, column] = value
        _commit(t, f"updated {name} ({', '.join(fields)})")


def replace(new_table):
    with STATE['lock']:
        _commit(new_table.copy(), "replaced")


def clear():
    with STATE['lock']:
        _commit(session_archive.empty_datasets_table(), "cleared")
