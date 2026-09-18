"""Tests for the Cytoscape controller with CyREST stubbed out.

The state is seeded from ``cytoscape_net.build`` as ``draw`` would leave it, and every
``cytoscape_p4c`` call the operations make is recorded instead of sent.
"""

import pandas as pd
import pytest

import cytoscape_ctl as ctl
import cytoscape_net as cn


PASSING = {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0}


def _row(bait, bait_acc, prey, gene, saint=0.9, intensity=1e6):
    return {'Experiment.ID': bait, 'Bait_Accession': bait_acc, 'First_ID': prey,
            'First_Prey_Gene': gene, 'SaintScore': saint, 'BFDR': 0.01, 'WD': 2.0,
            'WDFDR': 0.01, 'FoldChange': 3.0, 'AvgIntensity': intensity}


@pytest.fixture
def scores():
    return pd.DataFrame([
        _row('BaitA', 'PA', 'P1', 'G1', intensity=1e5),
        _row('BaitA', 'PA', 'P2', 'G2', intensity=1e7),
        _row('BaitA', 'PA', 'P4', 'G4'),
        _row('BaitB', 'PB', 'P2', 'G2'),
        _row('BaitB', 'PB', 'P5', 'G5'),
    ])


@pytest.fixture
def cytoscape(monkeypatch):
    """Stub every CyREST touch; ``calls`` records them in order."""
    calls = []
    stub = {'selected': [], 'positions': {}}
    monkeypatch.setattr(ctl.cy, 'selected_nodes', lambda net: list(stub['selected']))
    monkeypatch.setattr(ctl.cy, 'current_positions', lambda net: dict(stub['positions']))
    monkeypatch.setattr(ctl.cy, 'select_nodes', lambda net, names, add=False: calls.append(('select', list(names), add)))
    monkeypatch.setattr(ctl.cy, 'set_positions', lambda net, pos: calls.append(('positions', dict(pos))))
    monkeypatch.setattr(ctl.cy, 'update_edge_columns', lambda net, frame: calls.append(('edges', frame)))
    monkeypatch.setattr(ctl.cy, 'update_node_columns', lambda net, frame: calls.append(('nodes', frame)))
    stub['calls'] = calls
    return stub


@pytest.fixture
def drawn(scores, tmp_path, cytoscape):
    """The state after a draw with a literature edge P1-P2, without touching Cytoscape."""
    biogrid = tmp_path / 'biogrid_summary.csv'
    pd.DataFrame([{'SWISS-PROT Accessions Interactor A': 'P1', 'SWISS-PROT Accessions Interactor B': 'P2',
                   'Publication Source': 'PUBMED:1; PUBMED:2', 'Multivalidated': False}]).to_csv(biogrid, index=False)
    nodes, edges = cn.build(scores, PASSING, biogrid_path=str(biogrid))
    ctl.STATE.update(dataset='d', title='ProxiMate: d', net_suid=1, nodes=nodes, edges=edges,
                     thresholds=dict(PASSING),
                     style={'width_source': 'abundance', 'literature_weighted': False, 'biogrid_scope': 'all'})
    yield nodes, edges
    ctl.STATE.update(dataset=None, title=None, net_suid=None, nodes=None, edges=None,
                     thresholds=None, style=None)


def _last(calls, kind):
    return [payload for name, *payload in calls if name == kind][-1]


# --- edge style -------------------------------------------------------------------------------

def test_restyling_pushes_width_and_visibility_without_moving_anything(drawn, cytoscape):
    changed = ctl.restyle_edges('uniform', False, 'multivalidated')

    frame = _last(cytoscape['calls'], 'edges')[0]
    assert changed == len(frame) > 0
    assert set(frame.columns) == {'name', 'visible', 'width'}
    assert not frame.loc[frame['name'] == 'P1 (literature) P2', 'visible'].item()
    assert ctl.STATE['edges']['width'][ctl.STATE['edges']['interaction'] == 'proximity'].nunique() == 1
    assert ctl.snapshot()['style'] == {'width_source': 'uniform', 'literature_weighted': False,
                                       'biogrid_scope': 'multivalidated'}
    assert not any(name == 'positions' for name, *_ in cytoscape['calls'])


def test_rethresholding_keeps_the_drawn_biogrid_scope(drawn, cytoscape):
    ctl.restyle_edges('abundance', False, 'multivalidated')

    ctl.apply_thresholds(PASSING)

    edges = ctl.STATE['edges']
    assert not edges.loc[edges['interaction'] == 'literature', 'visible'].any()


# --- selection tools --------------------------------------------------------------------------

def test_loners_need_exactly_one_selected_bait(drawn, cytoscape):
    cytoscape['selected'] = ['P1']
    with pytest.raises(ValueError, match='exactly one bait'):
        ctl.select_loners()

    cytoscape['selected'] = ['PA', 'PB']
    with pytest.raises(ValueError, match='BaitA, BaitB'):
        ctl.select_loners()


def test_loners_are_selected_with_their_bait(drawn, cytoscape):
    cytoscape['selected'] = ['PA', 'P2']

    chosen = ctl.select_loners()

    assert chosen == ['PA', 'P4']
    assert _last(cytoscape['calls'], 'select') == [['PA', 'P4'], False]


def test_satellites_use_the_positions_cytoscape_holds(drawn, cytoscape):
    cytoscape['selected'] = ['PA']
    cytoscape['positions'] = {'PA': (0, 0), 'PB': (100, 0), 'P1': (0, 5), 'P2': (10, 0),
                              'P4': (5, 5), 'P5': (100, 5)}

    chosen = ctl.select_satellites()

    assert chosen == ['PA', 'P1', 'P4', 'P2']


def test_related_selection_can_add_to_the_current_one(drawn, cytoscape):
    chosen = ctl.select_related('BaitB', 'interactors', add=True, min_saint=0.5)

    assert chosen == ['P2', 'P5']
    assert _last(cytoscape['calls'], 'select') == [['P2', 'P5'], True]


def test_an_empty_relation_is_an_error_not_a_silent_deselect(drawn, cytoscape):
    with pytest.raises(ValueError, match='no interactors'):
        ctl.select_related('BaitB', 'interactors', min_saint=0.99)
    assert not any(name == 'select' for name, *_ in cytoscape['calls'])


# --- clustering ---------------------------------------------------------------------------------

def test_clustering_recolours_numbers_and_moves_only_the_selection(drawn, cytoscape):
    cytoscape['selected'] = ['PA', 'P1', 'P2', 'P4']
    cytoscape['positions'] = {i: (10.0 * k, 50.0) for k, i in enumerate(['PA', 'P1', 'P2', 'P4', 'PB', 'P5'])}

    result = ctl.cluster_selection(resolution=1.0, seed=1)

    assert result['n'] == 4 and sum(result['sizes']) == 4
    frame = _last(cytoscape['calls'], 'nodes')[0]
    assert set(frame['name']) == {'PA', 'P1', 'P2', 'P4'}
    assert set(frame.columns) == {'name', 'community', 'fill'}
    moved = _last(cytoscape['calls'], 'positions')[0]
    assert set(moved) == {'PA', 'P1', 'P2', 'P4'}
    assert all(x >= 0 and y >= 50 for x, y in moved.values())
    nodes = ctl.STATE['nodes'].set_index('id')
    assert nodes.loc['PB', 'community'] != nodes.loc['PB', 'community']   # NaN: untouched
    assert nodes.loc['PA', 'fill'] in cn.COMMUNITY_FILL


def test_a_second_clustering_numbers_above_the_first(drawn, cytoscape):
    cytoscape['positions'] = {i: (10.0 * k, 0.0) for k, i in enumerate(['PA', 'P1', 'P2', 'P4', 'PB', 'P5'])}
    cytoscape['selected'] = ['PA', 'P1', 'P2', 'P4']
    ctl.cluster_selection()
    first = ctl.STATE['nodes'].set_index('id').loc[['PA', 'P1', 'P2', 'P4'], 'community'].max()

    cytoscape['selected'] = ['PB', 'P5', 'P2', 'P4']
    ctl.cluster_selection()

    second = ctl.STATE['nodes'].set_index('id').loc[['PB', 'P5'], 'community'].min()
    assert second > first


def test_clustering_nothing_is_an_error(drawn, cytoscape):
    with pytest.raises(ValueError, match='nothing is selected'):
        ctl.cluster_selection()
