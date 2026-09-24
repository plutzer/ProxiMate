"""The MCP surface: the operation registry behind call_tool, the five tools, and the
documentation they serve.

Sandbox operations must leave a dataset directory byte-identical; dataset operations
must land in the store the GUI polls; Cytoscape operations must reach the controller
with actor 'mcp'.
"""

import base64
import json
import os

import pandas as pd
import pytest

import backend
import cytoscape_ctl as ctl
import dataset_store as store
import help_text
import mcp_ops
import mcp_registry as registry
import mcp_tools


THRESHOLDS = {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0}
PNG = b'\x89PNG\r\n\x1a\nstub'


@pytest.fixture
def out_dir(tmp_path):
    backend.configure(str(tmp_path))
    return tmp_path


def _tree(root):
    return sorted((os.path.relpath(os.path.join(d, f), root), os.path.getsize(os.path.join(d, f)))
                  for d, _, files in os.walk(root) for f in files)


# --- registry -------------------------------------------------------------------------

def test_every_op_has_a_mode_summary_and_schema():
    assert len(registry.OPS) >= 20
    for name, op in registry.OPS.items():
        assert op.mode in registry.MODES, name
        assert op.summary and len(op.summary) < 200, name
        assert op.schema['type'] == 'object', name
        assert set(op.schema['required']) <= set(op.schema['properties']), name


def test_search_matches_name_summary_and_tags_and_filters_by_mode():
    hits = {h['name'] for h in registry.search('threshold')}
    assert {'threshold_metrics', 'cytoscape_apply_thresholds'} <= hits
    assert all(h['mode'] == 'sandbox' for h in registry.search('', mode='sandbox'))
    assert {h['name'] for h in registry.search('', mode='sandbox')} == \
        {'get_scores', 'threshold_metrics', 'feature_analysis', 'get_prey_annotations', 'compare_networks'}
    assert len(registry.search('')) == len(registry.OPS)
    with pytest.raises(ValueError, match='mode'):
        registry.search('', mode='bogus')


def test_details_carry_the_parameters_with_types_and_defaults():
    d = registry.details('score_dataset')
    props = d['schema']['properties']
    assert props['imputation']['type'] == 'integer'
    assert props['exclude_hcm']['type'] == 'boolean'
    assert props['pi_bait']['default'] is None
    assert 'dataset' in d['schema']['required'] and 'seed' not in d['schema']['required']
    assert d['mode'] == 'dataset' and 'GUI' in d['effects']
    with pytest.raises(KeyError, match='nope'):
        registry.details('nope')


def test_call_validates_arguments_before_running():
    with pytest.raises(ValueError, match='unexpected'):
        registry.call('list_datasets', {'bogus': 1})
    with pytest.raises(ValueError, match='dataset'):
        registry.call('get_dataset_info', {})
    with pytest.raises(ValueError, match='integer'):
        registry.call('score_dataset', {'dataset': 'x', 'imputation': 'two', 'wdfdr_iterations': 0,
                                        'organism': 'human', 'exclude_hcm': False})
    with pytest.raises(KeyError):
        registry.call('nope', {})


# --- ops over a scored dataset -----------------------------------------------------------

def _score(bait, prey, **overrides):
    row = {"Experiment.ID": bait, "Prey.ID": prey, "First_ID": prey, "Bait.ID": bait,
           "First_Prey_Gene": prey.lower(), "SaintScore": 0.9, "BFDR": 0.01, "WD": 5.0,
           "WDFDR": 0.01, "In.BioGRID": False, "GO_CC": "Nucleus" if prey != "P3" else "Vesicle",
           "FoldChange": 3.0, "AvgIntensity": 1e6}
    row.update(overrides)
    return row


@pytest.fixture
def scored(out_dir):
    root = out_dir / 'ds'
    root.mkdir()
    with open(root / "interaction.txt", "w", newline="") as handle:
        for e, b, p, i in [("a_1", "BaitA", "P1", 100.0), ("a_2", "BaitA", "P1", 300.0),
                           ("a_1", "BaitA", "P2", 48.0), ("a_2", "BaitA", "P2", 52.0),
                           ("a_1", "BaitA", "P3", 10.0), ("a_2", "BaitA", "P3", 10.0),
                           ("b_1", "BaitB", "P1", 25.0), ("b_2", "BaitB", "P1", 75.0),
                           ("b_1", "BaitB", "P2", 49.0), ("b_2", "BaitB", "P2", 51.0)]:
            handle.write(f"{e}\t{b}\t{p}\t{i}\n")
    pd.DataFrame([{"Experiment Name": e, "Type": "T", "Bait": b, "Replicate": r, "Bait ID": b}
                  for e, b, r in [("a_1", "BaitA", 1), ("a_2", "BaitA", 2),
                                  ("b_1", "BaitB", 1), ("b_2", "BaitB", 2)]]
                 ).to_csv(root / "ED.csv", index=False)
    pd.DataFrame([_score("BaitA", "P1"), _score("BaitA", "P2"), _score("BaitA", "P3"),
                  _score("BaitB", "P1"), _score("BaitB", "P2", SaintScore=0.1)]
                 ).to_csv(root / "annotated_scores.csv", index=False)
    store.append({'Dataset Name': 'ds', 'Input Type': 'SAINT', 'Quant Type': 'Intensity',
                  'Experiments': 4, 'Controls': 0, 'Scored': 'Yes'})
    return root


def test_read_ops_describe_the_session(scored):
    assert registry.call('list_datasets', {})['datasets'][0]['Dataset Name'] == 'ds'
    info = registry.call('get_dataset_info', {'dataset': 'ds'})
    assert info['baits'] == ['BaitA', 'BaitB']
    status = registry.call('server_status', {})
    assert status['jobs'] == {} and status['network']['dataset'] is None and 'cytoscape' not in status


def test_sandbox_ops_return_json_and_leave_the_dataset_untouched(scored, tmp_path, monkeypatch):
    biogrid = tmp_path / 'bg.csv'
    biogrid.write_text("SWISS-PROT Accessions Interactor A,SWISS-PROT Accessions Interactor B\nP1,P2\n")
    import QC_plots
    monkeypatch.setattr(backend, 'calculate_threshold_metrics',
                        lambda path, t, ctrl_experiments=None, biogrid_path=None:
                        QC_plots.calculate_threshold_metrics(path, t, ctrl_experiments, str(biogrid)))
    before = _tree(scored)
    metrics = registry.call('threshold_metrics', {'dataset': 'ds', 'thresholds': THRESHOLDS})
    features = registry.call('feature_analysis', {'dataset': 'ds', 'feature_types': ['GO_CC'],
                                                  'thresholds': {**THRESHOLDS, 'SaintScore': 0.0}})
    comparison = registry.call('compare_networks', {'dataset': 'ds', 'bait_a': 'BaitA', 'bait_b': 'BaitB',
                                                    'thresholds_a': THRESHOLDS, 'thresholds_b': THRESHOLDS})
    for result in (metrics, features, comparison):
        json.dumps(result)          # every result is plain JSON
    assert metrics['total_after'] == 4
    assert features['columns'][:3] == ['Bait', 'Feature', 'Feature_type']
    assert comparison['genes']['both'] == ['p1'] and comparison['volcano']['n'] == 3
    lists_only = registry.call('compare_networks', {'dataset': 'ds', 'bait_a': 'BaitA', 'bait_b': 'BaitB',
                                                    'thresholds_a': THRESHOLDS, 'thresholds_b': THRESHOLDS,
                                                    'include_volcano': False})
    assert lists_only['genes'] == comparison['genes'] and 'volcano' not in lists_only
    scores = registry.call('get_scores', {'dataset': 'ds', 'thresholds': THRESHOLDS, 'top_n': 2})
    assert scores['n'] == 4 and scores['returned'] == 2 and scores['n_per_bait'] == {'BaitA': 3, 'BaitB': 1}
    assert scores['rows'][0]['First_Prey_Gene'] == 'p1' and scores['rows'][0]['In.BioGRID'] is False
    only_b = registry.call('get_scores', {'dataset': 'ds', 'thresholds': THRESHOLDS, 'baits': ['BaitB']})
    assert [r['Prey.ID'] for r in only_b['rows']] == ['P1']
    assert _tree(scored) == before
    assert store.version() == store.version()


def test_run_info_and_tail_log_read_the_manifest_and_the_log(out_dir):
    files = {'bait': out_dir / 'bait.txt', 'prey': out_dir / 'prey.txt', 'interaction': out_dir / 'interaction.txt'}
    files['bait'].write_text("t1\tBaitA\tT\nc1\tCtrl\tC\n")
    files['prey'].write_text("P1\tG1\n")
    files['interaction'].write_text("t1\tBaitA\tP1\t10\nc1\tCtrl\tP1\t4\n")
    registry.call('parse_dataset', {'dataset': 'new', 'input_format': 'SAINT',
                                    'files': {k: str(v) for k, v in files.items()},
                                    'quant_type': 'Spectral Counts'})
    info = registry.call('get_run_info', {'dataset': 'new'})
    stage = info['runs'][0]['stages'][0]
    assert stage['stage'] == 'parse' and stage['status'] == 'ok' and stage['error'] is None
    assert stage['wall_seconds'] is not None and 'ED.csv' in ''.join(stage['outputs'])
    json.dumps(info)
    tail = registry.call('tail_log', {'dataset': 'new', 'n_lines': 5})
    assert tail['lines'] and tail['n_lines_total'] >= len(tail['lines'])
    with pytest.raises(FileNotFoundError):
        registry.call('tail_log', {'dataset': 'absent'})


def test_score_dataset_defaults_are_the_scoring_cards(monkeypatch):
    params = registry.details('score_dataset')['schema']['properties']
    assert {k: params[k]['default'] for k in ('imputation', 'wdfdr_iterations', 'organism', 'exclude_hcm')} == \
        {'imputation': 0, 'wdfdr_iterations': 1000, 'organism': 'human', 'exclude_hcm': False}
    seen = {}
    monkeypatch.setattr(backend, 'run_score', lambda name, *a, **k: seen.update(args=a, **k) or {'run_id': 'r'})
    registry.call('score_dataset', {'dataset': 'ds'})
    assert seen['args'] == (0, 1000, 'human', False)


def test_prey_annotations_summarize_each_prey_across_baits(scored, tmp_path):
    before = sorted(os.listdir(scored))
    result = registry.call('get_prey_annotations', {'dataset': 'ds', 'thresholds': THRESHOLDS})
    assert sorted(os.listdir(scored)) == before
    rows = {r['First_ID']: r for r in result['rows']}
    assert result['n'] == 3 and result['unmatched'] == []
    # P1 passes under both baits; BaitB's P2 row scores 0.1 and fails
    assert rows['P1']['passing_baits'] == ['BaitA', 'BaitB']
    assert rows['P2']['passing_baits'] == ['BaitA'] and rows['P2']['n_baits_seen'] == 2
    assert rows['P3']['GO_CC'] == 'Vesicle' and rows['P3']['known_baits'] == []
    assert rows['P1']['max_saint'] == 0.9 and rows['P1']['max_fold_change'] == 3.0
    assert result['rows'][0]['First_ID'] == 'P1'  # most passing baits first
    # first_SCL, Main location and Human_Complex are absent from this annotation
    assert 'Human_Complex' not in result['columns']

    # ids match accessions or gene symbols, case-insensitively; absent ones are listed
    result = registry.call('get_prey_annotations', {'dataset': 'ds', 'thresholds': THRESHOLDS,
                                                    'ids': ['p1', 'P3', 'NOPE']})
    assert [r['First_ID'] for r in result['rows']] == ['P1', 'P3']
    assert result['unmatched'] == ['NOPE']


def test_feature_analysis_truncates_to_top_n(scored):
    result = registry.call('feature_analysis', {'dataset': 'ds', 'feature_types': ['GO_CC'],
                                                'thresholds': {**THRESHOLDS, 'SaintScore': 0.0},
                                                'top_n': 1})
    assert result['n'] >= result['returned'] == min(1, result['n'])


def test_parse_dataset_lands_in_the_store_and_records_the_actor(out_dir):
    files = {'bait': out_dir / 'bait.txt', 'prey': out_dir / 'prey.txt', 'interaction': out_dir / 'interaction.txt'}
    files['bait'].write_text("t1\tBaitA\tT\nc1\tCtrl\tC\n")
    files['prey'].write_text("P1\tG1\n")
    files['interaction'].write_text("t1\tBaitA\tP1\t10\nc1\tCtrl\tP1\t4\n")
    before = store.version()
    result = registry.call('parse_dataset', {'dataset': 'new', 'input_format': 'SAINT',
                                             'files': {k: str(v) for k, v in files.items()},
                                             'quant_type': 'Spectral Counts'})
    assert result['Dataset Name'] == 'new' and store.version() > before
    manifest = json.load(open(out_dir / 'new' / 'run.json'))
    stages = [s for run in manifest['runs'].values() for s in run['stages']]
    assert stages[0]['stage'] == 'parse' and stages[0]['status'] == 'ok'
    log = registry.activity()
    assert log[-1]['op'] == 'parse_dataset' and log[-1]['actor'] == 'mcp' and log[-1]['ok']


def test_load_session_needs_confirmation(scored, tmp_path):
    with pytest.raises(ValueError, match='confirm'):
        registry.call('load_session', {'zip_path': str(tmp_path / 'x.zip')})


# --- cytoscape ops reach the controller with actor mcp ---------------------------------------

@pytest.fixture
def drawn(monkeypatch):
    calls = []
    stub = {'selected': ['P1'], 'positions': {'P1': (0.0, 0.0), 'P2': (1.0, 1.0), 'BA': (2.0, 2.0)}}
    monkeypatch.setattr(ctl.cy, 'selected_nodes', lambda net: list(stub['selected']))
    monkeypatch.setattr(ctl.cy, 'current_positions', lambda net: dict(stub['positions']))
    monkeypatch.setattr(ctl.cy, 'select_nodes', lambda net, names, add=False: calls.append(('select', list(names), add)))
    monkeypatch.setattr(ctl.cy, 'clear_selection', lambda net: calls.append('clear'))
    monkeypatch.setattr(ctl.cy, 'render_png', lambda net, height: calls.append(('render', height)) or PNG)
    monkeypatch.setattr(ctl.cy, 'set_positions', lambda net, pos: calls.append(('positions', dict(pos))))
    monkeypatch.setattr(ctl.cy, 'update_edge_columns', lambda net, frame: calls.append(('edges', frame)))
    monkeypatch.setattr(ctl.cy, 'update_node_columns', lambda net, frame: calls.append(('nodes', frame)))
    nodes = pd.DataFrame({'id': ['BA', 'P1', 'P2'], 'symbol': ['GA', 'G1', 'G2'], 'accession': ['QA', 'P1', 'P2'],
                          'role': ['bait', 'prey', 'prey'], 'x': [0.0] * 3, 'y': [0.0] * 3})
    edges = pd.DataFrame({'name': ['e1', 'e2'], 'source': ['BA', 'BA'], 'target': ['P1', 'P2'],
                          'interaction': ['proximity'] * 2, 'visible': [True, True],
                          'SaintScore': [0.9, 0.8], 'BFDR': [0.01, 0.02], 'WD': [2.0, 1.0],
                          'WDFDR': [0.01, 0.01], 'width': [2.0, 2.0], 'kind': ['proximity'] * 2})
    ctl.STATE.update(dataset='ds', title='ProxiMate: ds', net_suid=1, nodes=nodes, edges=edges,
                     thresholds=dict(THRESHOLDS),
                     style={'width_source': 'abundance', 'literature_weighted': False, 'biogrid_scope': 'all'})
    stub['calls'] = calls
    yield stub
    ctl.STATE.update(dataset=None, title=None, net_suid=None, nodes=None, edges=None,
                     thresholds=None, style=None)


def test_cytoscape_ops_act_on_the_drawn_network_as_mcp(drawn):
    status = registry.call('server_status', {})
    assert status['network']['dataset'] == 'ds' and status['network']['thresholds'] == THRESHOLDS
    sel = registry.call('cytoscape_read_selection', {})
    assert sel['nodes'][0]['id'] == 'P1' and sel['edges'][0]['SaintScore'] == 0.9
    assert registry.call('cytoscape_select_nodes', {'ids': ['G2'], 'add': True}) == {'selected': ['P2'], 'added': True}
    assert drawn['calls'][-1] == ('select', ['P2'], True)
    assert registry.call('cytoscape_get_positions', {'ids': ['P1']}) == {'positions': {'P1': [0.0, 0.0]}}
    nodes = registry.call('cytoscape_list_nodes', {})['nodes']
    assert nodes[0] == {'id': 'BA', 'accession': 'QA', 'symbol': 'GA', 'role': 'bait'}
    assert [n['id'] for n in registry.call('cytoscape_list_nodes', {'role': 'prey'})['nodes']] == ['P1', 'P2']
    with pytest.raises(ValueError, match='role'):
        registry.call('cytoscape_list_nodes', {'role': 'edge'})
    assert registry.call('cytoscape_move_nodes', {'positions': {'P1': [5, 5]}}) == {'moved': 1}
    assert registry.call('cytoscape_clear_selection', {}) == {'selected': []}
    assert drawn['calls'][-1] == 'clear'
    drawn['positions'].update(P2=(1.0, 1.0), BA=(0.0, 0.0))
    related = registry.call('cytoscape_select_related', {'seed': 'GA', 'relation': 'satellites', 'include_seed': True})
    assert related == {'selected': ['BA', 'P1', 'P2'], 'added': False}
    assert registry.call('cytoscape_select_related', {'seed': 'BA', 'relation': 'singletons'})['selected'] == ['P1', 'P2']
    with pytest.raises(ValueError, match='needs a bait'):
        registry.call('cytoscape_select_related', {'seed': 'P1', 'relation': 'satellites'})
    with pytest.raises(ValueError, match='clear_selection'):
        registry.call('cytoscape_select_nodes', {'ids': []})
    assert ctl.snapshot()['log'][-1]['actor'] == 'mcp'


def test_view_network_renders_the_fitted_network_without_writing(drawn, out_dir):
    image = mcp_tools.view_network(height=300)
    assert image.data == PNG and image._mime_type == 'image/png'
    assert drawn['calls'][-1] == ('render', 300)
    assert not (out_dir / 'ds' / 'cytoscape').exists()
    assert ctl.snapshot()['log'][-1] == {**ctl.snapshot()['log'][-1], 'actor': 'mcp', 'op': 'view_image'}


def test_cytoscape_send_records_its_settings_in_the_manifest(scored, drawn, monkeypatch):
    seen = {}

    def fake_draw(dataset, scores_path, thresholds, **kwargs):
        seen.update(dataset=dataset, thresholds=thresholds, **kwargs)
        return ctl.snapshot()
    monkeypatch.setattr(ctl, 'draw', fake_draw)
    monkeypatch.setattr(mcp_ops.provenance, 'biogrid_summary_path', lambda organism, exclude_hcm=False: '/no/biogrid.csv')
    with pytest.raises(ValueError, match="'ds' is drawn; pass replace=true"):
        registry.call('cytoscape_send', {'dataset': 'ds', 'thresholds': THRESHOLDS})
    assert not seen
    result = registry.call('cytoscape_send', {'dataset': 'ds', 'thresholds': THRESHOLDS, 'baits': ['BaitA'],
                                              'layout': 'grid', 'width_source': 'SaintScore', 'replace': True})
    assert seen['actor'] == 'mcp' and seen['baits'] == ['BaitA'] and seen['layout'] == 'grid'
    assert result['network']['dataset'] == 'ds'
    manifest = json.load(open(scored / 'run.json'))
    stage = [s for run in manifest['runs'].values() for s in run['stages']][-1]
    assert stage['stage'] == 'cytoscape' and stage['entrypoint'] == 'mcp.cytoscape_send'
    assert stage['params']['thresholds'] == THRESHOLDS and stage['params']['width_source'] == 'SaintScore'
    with pytest.raises(ValueError, match='width_source'):
        registry.call('cytoscape_send', {'dataset': 'ds', 'thresholds': THRESHOLDS, 'width_source': 'nope',
                                         'replace': True})


def test_status_names_the_network_in_cytoscapes_window(drawn, monkeypatch):
    monkeypatch.setattr(ctl, 'health', lambda: {'ok': True, 'url': 'u', 'version': '3.10', 'error': None})
    monkeypatch.setattr(ctl.cy, 'current_network_title', lambda: 'ProxiMate: other')
    status = registry.call('server_status', {'probe': True})
    assert status['cytoscape']['ok'] is True
    assert status['current_network'] == {'title': 'ProxiMate: other', 'is_drawn': False}
    monkeypatch.setattr(ctl.cy, 'current_network_title', lambda: 'ProxiMate: ds')
    assert registry.call('server_status', {'probe': True})['current_network']['is_drawn'] is True
    monkeypatch.setattr(ctl, 'health', lambda: {'ok': False, 'url': 'u', 'version': None, 'error': 'down'})
    assert 'current_network' not in registry.call('server_status', {'probe': True})


# --- the five tools --------------------------------------------------------------------------

def test_tool_wrappers_return_json_envelopes(scored):
    assert mcp_tools.search_tools('dataset', mode='read')[0]['mode'] == 'read'
    assert mcp_tools.get_tool_details('list_datasets')['name'] == 'list_datasets'
    ok = mcp_tools.call_tool('list_datasets', {})
    assert ok['ok'] is True and ok['result']['datasets'][0]['Dataset Name'] == 'ds' and ok['run_id']
    bad = mcp_tools.call_tool('get_dataset_info', {'dataset': 'missing'})
    assert bad['ok'] is False and 'missing' in bad['error'] and bad['error_type'] == 'KeyError'


def test_gui_documentation_covers_every_tab_and_every_tooltip():
    toc = mcp_tools.get_gui_documentation()
    assert set(toc['sections']) == set(help_text.SECTIONS)
    assert all(toc['overview'][tab] for tab in toc['sections'])
    covered = set()
    for tab in help_text.SECTIONS:
        page = mcp_tools.get_gui_documentation(tab)
        assert page['section'] == tab and page['text'].startswith('## ')
        covered |= set(page['controls'])
    assert covered == set(help_text.TOOLTIPS)
    cy = mcp_tools.get_gui_documentation('Cytoscape')
    assert 'Select Loners' in cy['text']            # the README section rides along
    with pytest.raises(KeyError):
        mcp_tools.get_gui_documentation('Nope')


# --- the server over a real streamable-http connection ------------------------------------------

@pytest.fixture
def mcp_url(scored):
    import socket
    import threading
    import uvicorn
    import server

    with socket.socket() as probe:
        probe.bind(('127.0.0.1', 0))
        port = probe.getsockname()[1]
    mcp = server.build_mcp(host='127.0.0.1', port=port)
    config = uvicorn.Config(mcp.streamable_http_app(), host='127.0.0.1', port=port, log_level='warning')
    uv = uvicorn.Server(config)
    thread = threading.Thread(target=uv.run, daemon=True)
    thread.start()
    import requests
    for _ in range(100):
        try:
            requests.get(f'http://127.0.0.1:{port}/api/health', timeout=1).raise_for_status()
            break
        except requests.RequestException:
            import time
            time.sleep(0.05)
    yield f'http://127.0.0.1:{port}'
    uv.should_exit = True
    thread.join(5)


def test_the_five_tools_answer_over_http(mcp_url, drawn):
    import anyio
    import requests
    from mcp import ClientSession
    from mcp.client.streamable_http import streamablehttp_client

    health = requests.get(mcp_url + '/api/health', timeout=5).json()
    assert health['datasets'] == ['ds'] and health['operations'] == len(registry.OPS)

    async def drive():
        async with streamablehttp_client(mcp_url + '/mcp') as (read, write, _):
            async with ClientSession(read, write) as session:
                await session.initialize()
                tools = await session.list_tools()
                names = sorted(t.name for t in tools.tools)
                found = await session.call_tool('search_tools', {'query': 'threshold_metrics'})
                details = await session.call_tool('get_tool_details', {'name': 'threshold_metrics'})
                called = await session.call_tool('call_tool', {'name': 'list_datasets', 'arguments': {}})
                docs = await session.call_tool('get_gui_documentation', {'section': 'Downloads'})
                picture = await session.call_tool('view_network', {'height': 300})
                return names, found, details, called, docs, picture

    def payload(result):
        assert not result.isError, result.content
        if result.structuredContent is not None:
            return result.structuredContent.get('result', result.structuredContent)
        return json.loads(result.content[0].text)

    names, found, details, called, docs, picture = anyio.run(drive)
    assert names == ['call_tool', 'get_gui_documentation', 'get_tool_details', 'search_tools', 'view_network']
    assert 'threshold_metrics' in json.dumps(payload(found))
    assert payload(details)['mode'] == 'sandbox'
    assert payload(called)['ok'] is True
    assert payload(called)['result']['datasets'][0]['Dataset Name'] == 'ds'
    assert payload(docs)['section'] == 'Downloads'
    assert not picture.isError, picture.content
    block = picture.content[0]
    assert block.type == 'image' and block.mimeType == 'image/png'
    assert base64.b64decode(block.data) == PNG


def test_activity_entries_are_sequenced_and_name_their_target(scored):
    start = registry.last_seq()
    registry.call('get_dataset_info', {'dataset': 'ds'})
    with pytest.raises(KeyError):
        registry.call('get_dataset_info', {'dataset': 'nope'})
    new = registry.activity_since(start)
    assert [e['seq'] for e in new] == [start + 1, start + 2]
    assert new[0]['detail'] == 'ds' and new[0]['ok']
    assert new[1]['detail'].startswith('nope: KeyError') and not new[1]['ok']
    assert registry.last_seq() == start + 2
    assert registry.activity_since(registry.last_seq()) == []
    registry.call('get_scores', {'dataset': 'ds', 'thresholds': THRESHOLDS, 'baits': ['BaitA'], 'top_n': 2})
    assert registry.activity(1)[0]['detail'] == f"ds thresholds={THRESHOLDS}, baits=['BaitA'], top_n=2 -> 2 rows"


def test_cytoscape_activity_entries_name_the_relation_and_the_count(drawn):
    registry.call('cytoscape_select_related', {'seed': 'BA', 'relation': 'singletons'})
    assert registry.activity(1)[0]['detail'] == 'seed=BA, relation=singletons -> 2 selected'
    registry.call('cytoscape_list_nodes', {'role': 'prey'})
    assert registry.activity(1)[0]['detail'] == 'role=prey -> 2 nodes'
