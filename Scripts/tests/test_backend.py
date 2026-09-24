"""The operations layer behind the GUI handlers and the MCP tools.

Dataset operations go through one job lock per dataset and land in the datasets
store; sandbox operations compute from a dataset directory and leave it untouched.
"""

import os
import threading

import pandas as pd
import pytest

import backend
import dataset_store as store


THRESHOLDS = {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0}


@pytest.fixture
def out_dir(tmp_path):
    backend.configure(str(tmp_path))
    return tmp_path


def _saint_files(tmp_path):
    bait = tmp_path / 'bait.txt'
    prey = tmp_path / 'prey.txt'
    inter = tmp_path / 'interaction.txt'
    bait.write_text("t1\tBaitA\tT\nt2\tBaitA\tT\nc1\tCtrl\tC\n")
    prey.write_text("P1\tG1\nP2\tG2\n")
    inter.write_text("t1\tBaitA\tP1\t10\nt2\tBaitA\tP1\t12\nt1\tBaitA\tP2\t3\nc1\tCtrl\tP2\t4\n")
    return {'bait': str(bait), 'prey': str(prey), 'interaction': str(inter)}


def _tree(root):
    """Every file under root with its size, for before/after comparison."""
    return sorted((os.path.relpath(os.path.join(d, f), root), os.path.getsize(os.path.join(d, f)))
                  for d, _, files in os.walk(root) for f in files)


# --- thresholds and paths ------------------------------------------------------

def test_thresholds_must_carry_exactly_the_four_scores():
    assert backend.validate_thresholds({'SaintScore': '0.7', 'BFDR': 0.05, 'WD': 0, 'WDFDR': 1}) == \
        {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0}
    with pytest.raises(ValueError, match='WDFDR'):
        backend.validate_thresholds({'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0})
    with pytest.raises(ValueError, match='Extra'):
        backend.validate_thresholds({**THRESHOLDS, 'Extra': 1})
    with pytest.raises(ValueError, match='SaintScore'):
        backend.validate_thresholds({**THRESHOLDS, 'SaintScore': 1.5})


def test_dataset_paths(out_dir):
    assert backend.dataset_dir('ds') == os.path.join(str(out_dir), 'ds')
    assert backend.results_path('ds').endswith(os.path.join('ds', 'annotated_scores.csv'))
    with pytest.raises(ValueError):
        backend.dataset_dir('../escape')


# --- parse -----------------------------------------------------------------------

def test_parse_saint_creates_the_dataset_and_its_row(out_dir):
    row = backend.run_parse('ds1', 'SAINT', _saint_files(out_dir), 'Spectral Counts', actor='mcp')
    assert row['Dataset Name'] == 'ds1' and row['Experiments'] == 2 and row['Controls'] == 1
    assert store.names() == ['ds1']
    produced = os.listdir(out_dir / 'ds1')
    assert {'ED.csv', 'bait.txt', 'prey.txt', 'interaction.txt', 'to_CompPASS.csv', 'run.json'} <= set(produced)
    ed = pd.read_csv(out_dir / 'ds1' / 'ED.csv', keep_default_na=False)
    assert ed['Bait ID'].tolist() == ['None'] * 3


def test_parse_refuses_a_taken_or_invalid_name(out_dir):
    files = _saint_files(out_dir)
    backend.run_parse('ds1', 'SAINT', files, 'Spectral Counts')
    with pytest.raises(ValueError, match='already'):
        backend.run_parse('ds1', 'SAINT', files, 'Spectral Counts')
    (out_dir / 'leftover').mkdir()
    with pytest.raises(ValueError):
        backend.run_parse('leftover', 'SAINT', files, 'Spectral Counts')
    with pytest.raises(ValueError):
        backend.run_parse('bad name', 'SAINT', files, 'Spectral Counts')


def test_parse_names_the_missing_file(out_dir):
    with pytest.raises(ValueError, match='interaction'):
        backend.run_parse('ds1', 'SAINT', {'bait': 'x', 'prey': 'y'}, 'Spectral Counts')
    with pytest.raises(ValueError, match='format'):
        backend.run_parse('ds1', 'Nope', {}, 'Intensity')


# --- score -----------------------------------------------------------------------

@pytest.fixture
def parsed(out_dir):
    backend.run_parse('ds1', 'SAINT', _saint_files(out_dir), 'Spectral Counts')
    return out_dir / 'ds1'


def _fake_stages(monkeypatch, calls, fail=None):
    """Stand in for the score.py / annotator.py subprocesses."""
    def run(command, dataset_path, run_id):
        calls.append(command)
        script = os.path.basename(command[1])
        if script == fail:
            return 3, 'boom'
        if script == 'score.py':
            pd.DataFrame({'a': [1]}).to_csv(os.path.join(dataset_path, 'merged.csv'), index=False)
        else:
            pd.DataFrame({'a': [1]}).to_csv(os.path.join(dataset_path, 'annotated_scores.csv'), index=False)
        return 0, ''
    monkeypatch.setattr(backend, '_run_stage_subprocess', run)


def test_score_runs_both_stages_with_the_given_settings_and_marks_the_row(out_dir, monkeypatch):
    backend.run_parse('ds1', 'SAINT', _saint_files(out_dir), 'Intensity')
    calls = []
    _fake_stages(monkeypatch, calls)
    result = backend.run_score('ds1', imputation=2, wdfdr_iterations=5, organism='human',
                               exclude_hcm=True, pi_method='single_bait', pi_bait='Ctrl', seed=7)
    score_cmd, ann_cmd = calls
    assert score_cmd[1].endswith('score.py')
    assert score_cmd[score_cmd.index('--quantType') + 1] == 'Intensity'
    assert score_cmd[score_cmd.index('--imputation') + 1] == '2'
    assert score_cmd[score_cmd.index('--n-iterations') + 1] == '5'
    assert score_cmd[score_cmd.index('--pi-bait') + 1] == 'Ctrl'
    assert score_cmd[score_cmd.index('--seed') + 1] == '7'
    assert ann_cmd[1].endswith('annotator.py') and '--excludeHCM' in ann_cmd
    assert store.row('ds1')['Scored'] == 'Yes'
    assert store.row('ds1')['Imputation'] == 'Refactored AFT'
    assert result['run_id']


def test_a_failing_stage_raises_with_its_log_and_leaves_the_row_unscored(parsed, monkeypatch):
    _fake_stages(monkeypatch, [], fail='annotator.py')
    with pytest.raises(backend.StageError, match='boom'):
        backend.run_score('ds1', imputation=0, wdfdr_iterations=0, organism='human', exclude_hcm=False)
    assert store.row('ds1')['Scored'] != 'Yes'


def test_run_info_and_log_tail_report_the_stages_and_the_log(parsed, monkeypatch):
    _fake_stages(monkeypatch, [])
    backend.run_score('ds1', imputation=0, wdfdr_iterations=0, organism='human', exclude_hcm=False)
    info = backend.run_info('ds1')
    # The faked score.py and annotator.py record no stages of their own.
    assert [(s['stage'], s['status']) for run in info['runs'] for s in run['stages']] == [('parse', 'ok')]
    assert info['n_runs'] == len(info['runs']) == 1
    assert 'environment' not in info['runs'][0] and 'traceback' not in str(info)
    assert backend.run_info('ds1', last_n=1) == info
    tail = backend.log_tail('ds1', n_lines=3)
    assert len(tail['lines']) == 3 and tail['n_lines_total'] >= 3
    assert any('Starting scoring' in line for line in backend.log_tail('ds1', n_lines=1000)['lines'])
    with pytest.raises(FileNotFoundError, match='run.json'):
        backend.run_info('ds2')


def test_score_refuses_an_unknown_dataset_or_bad_settings(parsed, monkeypatch):
    _fake_stages(monkeypatch, [])
    with pytest.raises(KeyError):
        backend.run_score('nope', imputation=0, wdfdr_iterations=0, organism='human', exclude_hcm=False)
    with pytest.raises(ValueError, match='imputation'):
        backend.run_score('ds1', imputation=9, wdfdr_iterations=0, organism='human', exclude_hcm=False)
    with pytest.raises(ValueError, match='organism'):
        backend.run_score('ds1', imputation=0, wdfdr_iterations=0, organism='cat', exclude_hcm=False)


@pytest.mark.parametrize("imputation", [1, 2, 3])
def test_score_refuses_imputation_for_spectral_counts(parsed, monkeypatch, imputation):
    """The AFT imputations model missing intensities; the spectral-count SAINT build
    cannot use their output, so the request is refused before any stage runs."""
    calls = []
    _fake_stages(monkeypatch, calls)
    with pytest.raises(ValueError, match='spectral counts'):
        backend.run_score('ds1', imputation=imputation, wdfdr_iterations=0, organism='human',
                          exclude_hcm=False)
    assert calls == []
    assert store.row('ds1')['Scored'] != 'Yes'


def test_a_dataset_being_worked_on_refuses_a_second_job(parsed, monkeypatch):
    started, release = threading.Event(), threading.Event()

    def slow(command, dataset_path, run_id):
        started.set()
        release.wait(5)
        return 3, ''
    monkeypatch.setattr(backend, '_run_stage_subprocess', slow)

    def first():
        with pytest.raises(backend.StageError):
            backend.run_score('ds1', imputation=0, wdfdr_iterations=0, organism='human',
                              exclude_hcm=False, actor='gui')
    worker = threading.Thread(target=first)
    worker.start()
    started.wait(5)
    job = backend.running_jobs()['ds1']
    assert (job['actor'], job['what']) == ('gui', 'score') and job['since']
    with pytest.raises(backend.BusyError, match='gui'):
        backend.run_score('ds1', imputation=0, wdfdr_iterations=0, organism='human', exclude_hcm=False)
    release.set()
    worker.join(5)
    assert backend.running_jobs() == {}


# --- sandbox operations -----------------------------------------------------------

def _score(bait, prey, **overrides):
    row = {"Experiment.ID": bait, "Prey.ID": prey, "First_ID": prey, "Bait.ID": bait,
           "First_Prey_Gene": prey.lower(), "SaintScore": 0.9, "BFDR": 0.01, "WD": 5.0, "WDFDR": 0.01, "In.BioGRID": False,
           "SCL": "Nucleus" if prey != "P3" else "Vesicle"}
    row.update(overrides)
    return row


@pytest.fixture
def scored(out_dir):
    root = out_dir / 'ds'
    root.mkdir()
    with open(root / "interaction.txt", "w", newline="") as handle:
        for experiment, bait, prey, intensity in [
                ("a_1", "BaitA", "P1", 100.0), ("a_2", "BaitA", "P1", 300.0),
                ("a_1", "BaitA", "P2", 48.0), ("a_2", "BaitA", "P2", 52.0),
                ("a_1", "BaitA", "P3", 10.0), ("a_2", "BaitA", "P3", 10.0),
                ("b_1", "BaitB", "P1", 25.0), ("b_2", "BaitB", "P1", 75.0),
                ("b_1", "BaitB", "P2", 49.0), ("b_2", "BaitB", "P2", 51.0)]:
            handle.write(f"{experiment}\t{bait}\t{prey}\t{intensity}\n")
    pd.DataFrame([
        {"Experiment Name": e, "Type": "T", "Bait": b, "Replicate": r, "Bait ID": b}
        for e, b, r in [("a_1", "BaitA", 1), ("a_2", "BaitA", 2), ("b_1", "BaitB", 1), ("b_2", "BaitB", 2)]
    ]).to_csv(root / "ED.csv", index=False)
    pd.DataFrame([_score("BaitA", "P1"), _score("BaitA", "P2"), _score("BaitA", "P3"),
                  _score("BaitB", "P1"), _score("BaitB", "P2", SaintScore=0.1)]
                 ).to_csv(root / "annotated_scores.csv", index=False)
    store.append({'Dataset Name': 'ds', 'Input Type': 'SAINT', 'Quant Type': 'Intensity',
                  'Experiments': 4, 'Controls': 0, 'Scored': 'Yes'})
    return root


def test_threshold_metrics_returns_counts_and_leaves_the_directory_alone(scored, tmp_path):
    biogrid = tmp_path / 'bg.csv'
    biogrid.write_text("SWISS-PROT Accessions Interactor A,SWISS-PROT Accessions Interactor B\nP1,P2\n")
    before = _tree(scored)
    metrics = backend.threshold_metrics('ds', THRESHOLDS, biogrid_path=str(biogrid))
    assert metrics['total_after'] == 4 and metrics['total_before'] == 5
    assert metrics['median_network_size'] == 2
    assert metrics['thresholds'] == THRESHOLDS
    assert _tree(scored) == before


def test_threshold_metrics_can_restrict_to_one_bait(scored, tmp_path):
    biogrid = tmp_path / 'bg.csv'
    biogrid.write_text("SWISS-PROT Accessions Interactor A,SWISS-PROT Accessions Interactor B\n")
    metrics = backend.threshold_metrics('ds', THRESHOLDS, bait='BaitB', biogrid_path=str(biogrid))
    assert metrics['total_before'] == 2 and metrics['total_after'] == 1
    with pytest.raises(ValueError, match='BaitZ'):
        backend.threshold_metrics('ds', THRESHOLDS, bait='BaitZ', biogrid_path=str(biogrid))


def test_feature_analysis_returns_the_table_without_writing_it(scored):
    before = _tree(scored)
    result = backend.feature_analysis('ds', {**THRESHOLDS, 'SaintScore': 0.0}, feature_types=['SCL'])
    assert list(result.columns)[:3] == ['Bait', 'Feature', 'Feature_type']
    assert _tree(scored) == before
    with pytest.raises(ValueError, match='Bogus'):
        backend.feature_analysis('ds', THRESHOLDS, feature_types=['Bogus'])


def test_compare_networks_returns_volcano_rows_and_gene_lists(scored):
    before = _tree(scored)
    out = backend.compare_networks('ds', 'BaitA', 'BaitB', THRESHOLDS, THRESHOLDS)
    assert out['genes'] == {'a_only': ['p2', 'p3'], 'b_only': [], 'both': ['p1']}
    assert set(out['volcano']['Prey.ID']) == {'P1', 'P2', 'P3'}
    assert _tree(scored) == before
    with pytest.raises(ValueError, match='BaitZ'):
        backend.compare_networks('ds', 'BaitA', 'BaitZ', THRESHOLDS, THRESHOLDS)


def test_passing_scores_filters_orders_and_restricts_to_baits(scored):
    rows = backend.passing_scores('ds', THRESHOLDS)
    assert list(zip(rows['Experiment.ID'], rows['Prey.ID'])) == [
        ('BaitA', 'P1'), ('BaitA', 'P2'), ('BaitA', 'P3'), ('BaitB', 'P1')]
    assert list(rows.columns[:5]) == ['Experiment.ID', 'Prey.ID', 'First_ID', 'First_Prey_Gene', 'SaintScore']
    assert list(backend.passing_scores('ds', THRESHOLDS, baits=['BaitB'])['Prey.ID']) == ['P1']
    with pytest.raises(ValueError, match='BaitZ'):
        backend.passing_scores('ds', THRESHOLDS, baits=['BaitZ'])


def test_dataset_info_summarizes_the_row_files_and_baits(scored):
    info = backend.dataset_info('ds')
    assert info['row']['Scored'] == 'Yes'
    assert info['baits'] == ['BaitA', 'BaitB']
    assert info['files']['annotated_scores.csv'] is True
    assert info['files']['Feature_enrichment.csv'] is False
