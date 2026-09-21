"""The datasets table as process-wide state: every mutation lands in datasets.csv and
bumps the version the GUI polls, so a dataset created outside a browser session
still appears in it."""

import pandas as pd
import pytest

import dataset_store as store
import session_archive


@pytest.fixture
def out_dir(tmp_path):
    store.configure(str(tmp_path))
    return tmp_path


def _row(name, scored=''):
    return {'Dataset Name': name, 'Input Type': 'SAINT', 'Quant Type': 'Intensity',
            'Experiments': 2, 'Controls': 1, 'Scored': scored, 'Imputation': '',
            'WDFDR iterations': ''}


def test_configure_loads_the_persisted_table_and_drops_rows_without_a_directory(tmp_path):
    (tmp_path / 'kept').mkdir()
    pd.DataFrame([_row('kept'), _row('gone')]).to_csv(tmp_path / 'datasets.csv', index=False)
    store.configure(str(tmp_path))
    assert store.names() == ['kept']


def test_append_persists_and_bumps_the_version(out_dir):
    before = store.version()
    store.append(_row('ds1'))
    assert store.version() == before + 1
    assert store.names() == ['ds1']
    on_disk = pd.read_csv(out_dir / 'datasets.csv')
    assert on_disk['Dataset Name'].tolist() == ['ds1']


def test_append_refuses_a_duplicate_name(out_dir):
    store.append(_row('ds1'))
    with pytest.raises(ValueError, match='ds1'):
        store.append(_row('ds1'))


def test_update_changes_one_row_and_refuses_an_unknown_name(out_dir):
    store.append(_row('ds1'))
    store.update('ds1', Scored='Yes', Imputation='Default')
    assert store.scored_names() == ['ds1']
    assert store.row('ds1')['Imputation'] == 'Default'
    with pytest.raises(KeyError):
        store.update('nope', Scored='Yes')


def test_replace_and_clear(out_dir):
    store.append(_row('ds1'))
    store.replace(pd.DataFrame([_row('a', 'Yes'), _row('b')]))
    assert store.names() == ['a', 'b']
    store.clear()
    assert store.names() == []
    assert list(pd.read_csv(out_dir / 'datasets.csv').columns) == session_archive.DATASET_COLUMNS


def test_table_returns_a_copy(out_dir):
    store.append(_row('ds1'))
    store.table().loc[0, 'Scored'] = 'Yes'
    assert store.scored_names() == []
