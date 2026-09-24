"""Tests for the Cytoscape tab: the node/edge tables (cytoscape_net), the controller
with CyREST stubbed out (cytoscape_ctl), and the CyREST helpers (cytoscape_p4c).

Nothing here may set a bypass: a locked view shows nothing on screen and stops the
mouse, so every bypass setter is made to raise for every test.
"""

import itertools

import numpy as np
import pandas as pd
import pytest

import cytoscape_ctl as ctl
import cytoscape_net as cn
import cytoscape_p4c as cy


PASSING = {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0}

BYPASS_SETTERS = ('set_node_position_bypass', 'set_network_zoom_bypass',
                  'set_network_center_bypass', 'set_node_property_bypass',
                  'set_edge_property_bypass', 'set_network_property_bypass')


@pytest.fixture(autouse=True)
def no_bypass(monkeypatch):
    """Every bypass setter raises, so any code path that reaches one fails the test."""
    def boom(*a, **k):
        raise AssertionError("a bypass setter was called")
    for name in BYPASS_SETTERS:
        if hasattr(cy.p4c, name):
            monkeypatch.setattr(cy.p4c, name, boom)


def _row(bait, bait_acc, prey, gene, saint=0.9, bfdr=0.01, wd=2.0, wdfdr=0.01,
         intensity=1e6, known=False):
    return {'Experiment.ID': bait, 'Bait_Accession': bait_acc, 'First_ID': prey,
            'First_Prey_Gene': gene, 'SaintScore': saint, 'BFDR': bfdr, 'WD': wd,
            'WDFDR': wdfdr, 'FoldChange': 3.0, 'AvgIntensity': intensity, 'In.BioGRID': known}


@pytest.fixture
def scores():
    return pd.DataFrame([
        _row('BaitA', 'PA', 'P1', 'G1', intensity=1e5, known=True),
        _row('BaitA', 'PA', 'P2', 'G2', intensity=1e7),
        _row('BaitA', 'PA', 'P3', 'G3', saint=0.2),          # fails SaintScore
        _row('BaitB', 'PB', 'P2', 'G2', intensity=1e6),
        _row('BaitB', 'PB', 'PA', 'GA'),                     # the other bait as a prey
    ])


def _biogrid(tmp_path, pairs):
    """A summary of ``(a, b)`` or ``(a, b, publications, multivalidated)`` rows."""
    path = tmp_path / 'biogrid_summary.csv'
    rows = []
    for a, b, *rest in pairs:
        pubs, mv = (rest + [1, False])[:2]
        rows.append({'SWISS-PROT Accessions Interactor A': a, 'SWISS-PROT Accessions Interactor B': b,
                     'Publication Source': '; '.join(f'PUBMED:{k}' for k in range(pubs)),
                     'Multivalidated': mv})
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _corum(tmp_path, complexes):
    """A CORUM table of ``(name, [accessions])`` rows."""
    path = tmp_path / 'corum_humanComplexes.txt'
    pd.DataFrame([{'complex_name': n, 'subunits_uniprot_id': ';'.join(s), 'organism': 'Human'}
                  for n, s in complexes]).to_csv(path, sep='\t', index=False)
    return str(path)


# =============================================================================================
# cytoscape_net: nodes
# =============================================================================================

def test_preys_are_keyed_on_accession_and_baits_on_name(scores):
    nodes, _ = cn.build(scores, PASSING, prey_prey=False)

    by_id = nodes.set_index('id')
    assert by_id.loc['BaitA', 'role'] == 'bait' and by_id.loc['BaitA', 'accession'] == 'PA'
    assert by_id.loc['P1', 'role'] == 'prey' and by_id.loc['P1', 'symbol'] == 'G1'
    assert 'P3' not in by_id.index


def test_a_prey_that_is_also_a_bait_keeps_one_bait_node(scores):
    nodes, _ = cn.build(scores, PASSING, prey_prey=False)

    assert (nodes['accession'] == 'PA').sum() == 1
    assert nodes.set_index('id').loc['BaitA', 'symbol'] == 'BaitA'


def test_two_constructs_of_one_protein_stay_separate_baits(scores):
    scores['Bait_Accession'] = 'PA'                     # BaitB is a second tag on the same protein
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'PA', 'GA')   # each bait detects itself

    nodes, edges = cn.build(scores, PASSING, prey_prey=False)

    baits = nodes[nodes['role'] == 'bait'].set_index('id')
    assert list(baits.index) == ['BaitA', 'BaitB'] and set(baits['accession']) == {'PA'}
    assert 'PA' not in set(nodes['id'])
    # neither self-detection draws an edge, and P2 gets one edge per bait
    assert not (edges['target'].isin(['BaitA', 'BaitB'])).any()
    assert sorted(edges.loc[edges['target'] == 'P2', 'source']) == ['BaitA', 'BaitB']


def test_the_bait_filter_restricts_both_tables(scores):
    nodes, edges = cn.build(scores, PASSING, baits=['BaitB'], prey_prey=False)

    assert set(nodes['id']) == {'BaitB', 'P2', 'PA'}
    assert set(edges['source']) == {'BaitB'}


def test_bait_id_is_used_when_no_accession_column_exists(scores):
    df = scores.rename(columns={'Bait_Accession': 'Bait.ID'})

    nodes, _ = cn.build(df, PASSING, prey_prey=False)

    assert 'PA' in set(nodes['accession'])


@pytest.mark.parametrize('policy, expected', [
    ('all', {'BaitA', 'G1', 'G2'}), ('baits', {'BaitA', ''}), ('none', {''})])
def test_label_policies(scores, policy, expected):
    nodes, _ = cn.build(scores, PASSING, baits=['BaitA'], prey_prey=False, label_policy=policy)

    assert set(nodes['display_label']) == expected


@pytest.mark.parametrize('kwargs', [
    {'label_policy': 'some'}, {'width_source': 'Entropy'}])
def test_an_unknown_option_raises(scores, kwargs):
    with pytest.raises(ValueError):
        cn.build(scores, PASSING, prey_prey=False, **kwargs)


# =============================================================================================
# cytoscape_net: edges
# =============================================================================================

def test_only_passing_interactions_become_edges_and_carry_their_scores(scores):
    _, edges = cn.build(scores, PASSING, prey_prey=False)

    assert len(edges) == 4
    assert set(edges['interaction']) == {'proximity'}
    assert {'SaintScore', 'BFDR', 'FoldChange', 'WD', 'WDFDR', 'AvgIntensity', 'In.BioGRID'} <= set(edges.columns)
    assert edges['visible'].all()
    assert edges['name'].iloc[0] == 'BaitA (proximity) P1'


def test_widths_follow_log_abundance(scores):
    _, edges = cn.build(scores, PASSING, prey_prey=False)

    by_target = edges.set_index(['source', 'target'])['width']
    assert by_target[('BaitA', 'P1')] < by_target[('BaitB', 'P2')] < by_target[('BaitA', 'P2')]


def test_spectral_counts_drive_widths_when_intensity_is_absent(scores):
    df = scores.rename(columns={'AvgIntensity': 'AvgSpec'})

    _, edges = cn.build(df, PASSING, prey_prey=False)

    assert edges['width'].nunique() > 1


def test_a_uniform_width_or_no_abundance_column_gives_one_width(scores):
    _, edges = cn.build(scores, PASSING, prey_prey=False, width_source='uniform')
    assert edges['width'].nunique() == 1

    _, edges = cn.build(scores.drop(columns=['AvgIntensity']), PASSING, prey_prey=False)
    assert edges['width'].nunique() == 1


@pytest.mark.parametrize('source, thinner, thicker', [
    ('SaintScore', ('BaitB', 'BaitA'), ('BaitA', 'P1')),
    ('WD', ('BaitA', 'P1'), ('BaitB', 'P2')),
    ('FoldChange', ('BaitA', 'P1'), ('BaitB', 'P2')),
])
def test_other_width_sources_follow_their_column(scores, source, thinner, thicker):
    scores.loc[scores['First_ID'] == 'P2', 'WD'] = 9.0
    scores.loc[scores['First_ID'] == 'P2', 'FoldChange'] = 30.0
    scores.loc[scores['First_ID'] == 'P1', 'SaintScore'] = 1.0
    scores.loc[scores['First_ID'] == 'BaitA', 'SaintScore'] = 0.7

    _, edges = cn.build(scores, PASSING, prey_prey=False, width_source=source)

    by_pair = edges.set_index(['source', 'target'])['width']
    assert by_pair[thinner] < by_pair[thicker]


def test_nothing_passing_raises(scores):
    with pytest.raises(ValueError):
        cn.build(scores, {**PASSING, 'SaintScore': 0.99}, prey_prey=False)


# --- prey-prey edges ------------------------------------------------------------------------

def test_literature_edges_join_passing_preys_but_never_a_bait(scores, tmp_path):
    biogrid = _biogrid(tmp_path, [('P1', 'P2'), ('P2', 'P1'), ('P1', 'PA'), ('P1', 'P9')])

    _, edges = cn.build(scores, PASSING, biogrid_path=biogrid)

    lit = edges[edges['interaction'] == 'literature']
    assert len(lit) == 1
    assert set(lit[['source', 'target']].iloc[0]) == {'P1', 'P2'}
    assert lit['name'].iloc[0] == 'P1 (literature) P2'


def test_an_unreadable_biogrid_summary_yields_no_literature_edges(scores, tmp_path):
    _, edges = cn.build(scores, PASSING, biogrid_path=str(tmp_path / 'absent.csv'))

    assert set(edges['interaction']) == {'proximity'}


def test_pairs_carry_the_larger_publication_count_and_any_multivalidation(scores, tmp_path):
    biogrid = _biogrid(tmp_path, [('P1', 'P2', 3, False), ('P2', 'P1', 1, True)])

    pairs = cn.prey_prey_pairs({'P1', 'P2'}, {'BaitA', 'BaitB'}, biogrid)

    assert len(pairs) == 1
    assert pairs['n_publications'].iloc[0] == 3 and bool(pairs['multivalidated'].iloc[0])


def test_literature_edges_thicken_with_publications_only_when_weighted(scores, tmp_path):
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'P4', 'G4')
    biogrid = _biogrid(tmp_path, [('P1', 'P2', 40, False), ('P2', 'P4', 1, False)])

    _, edges = cn.build(scores, PASSING, biogrid_path=biogrid)
    lit = edges[edges['interaction'] == 'literature']
    assert set(lit['width']) == {cn.LITERATURE_WIDTH[0]}

    _, edges = cn.build(scores, PASSING, biogrid_path=biogrid, literature_weighted=True)
    by_pair = edges[edges['interaction'] == 'literature'].set_index(['source', 'target'])['width']
    assert by_pair[('P2', 'P4')] < by_pair[('P1', 'P2')] <= cn.LITERATURE_WIDTH[1]


def test_the_multivalidated_scope_hides_single_evidence_pairs(scores, tmp_path):
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'P4', 'G4')
    biogrid = _biogrid(tmp_path, [('P1', 'P2', 3, True), ('P2', 'P4', 1, False)])
    _, edges = cn.build(scores, PASSING, biogrid_path=biogrid)

    visible = cn.edge_visibility(edges, PASSING, biogrid_scope='multivalidated')

    assert visible[edges['name'] == 'P1 (literature) P2'].all()
    assert not visible[edges['name'] == 'P2 (literature) P4'].any()
    assert visible[edges['interaction'] == 'proximity'].all()


# --- complex edges ----------------------------------------------------------------------

@pytest.fixture
def complexes(tmp_path):
    return _corum(tmp_path, [
        ('Trio', ['P1', 'P2', 'P9-2']),          # P9 is not drawn: 2 of 3 drawn, BaitA has 2/3
        ('Pair', ['P1', 'P2']),                  # only two subunits
        ('Wide', ['P1', 'P2', 'P4', 'X1', 'X2', 'X3']),   # 3 drawn but no bait covers half
        ('Mine', ['PA', 'P1', 'P2', 'P4']),      # bait subunit: BaitA covers 4/4
    ])


def test_only_recovered_complexes_draw_and_link_every_drawn_pair(scores, complexes):
    scores.loc[len(scores)] = _row('BaitB', 'PB', 'P4', 'G4')

    nodes, edges = cn.build(scores, PASSING, prey_prey=False, corum_path=complexes,
                            corum_min_members=3, corum_min_fraction=0.5)

    complex_edges = edges[edges['interaction'] == 'complex']
    # Mine qualifies; PA-P1 and PA-P2 are proximity edges and are left out; P4 is BaitB's
    # prey but a subunit all the same.
    assert set(complex_edges['name']) == {'P1 (complex) P2', 'P1 (complex) P4', 'P2 (complex) P4', 'P4 (complex) BaitA'}
    assert set(complex_edges['complex_names']) == {'Mine'}
    assert nodes.set_index('id').loc['P4', 'complexes'] == 'Mine'


def test_lower_criteria_admit_more_complexes(scores, complexes):
    _, edges = cn.build(scores, PASSING, prey_prey=False, corum_path=complexes,
                        corum_min_members=2, corum_min_fraction=0.6)

    names = set(';'.join(edges.loc[edges['interaction'] == 'complex', 'complex_names']).split(';'))
    assert names == {'Trio', 'Pair', 'Mine'}


def test_a_complex_edge_follows_its_ends_when_thresholds_tighten(scores, complexes):
    _, edges = cn.build(scores, PASSING, prey_prey=False, corum_path=complexes)

    visible = cn.edge_visibility(edges, {**PASSING, 'SaintScore': 0.95})

    assert not visible[edges['interaction'] == 'complex'].any()


# --- re-thresholding in place -------------------------------------------------------------

def test_tightening_hides_edges_and_literature_edges_follow_their_preys(scores, tmp_path):
    biogrid = _biogrid(tmp_path, [('P1', 'P2')])
    _, edges = cn.build(scores, PASSING, biogrid_path=biogrid)

    visible = cn.edge_visibility(edges, {**PASSING, 'SaintScore': 0.95})

    assert not visible.any()
    edges['visible'] = visible
    faded = cn.node_alpha(pd.DataFrame({'id': ['BaitA', 'P1']}), edges)
    assert list(faded) == [cn.NODE_ALPHA['faded']] * 2


@pytest.mark.parametrize('new, tighter', [
    ({**PASSING, 'SaintScore': 0.9}, True),
    ({**PASSING, 'BFDR': 0.1}, False),
    (PASSING, True),
])
def test_only_tighter_or_equal_thresholds_can_be_applied_in_place(new, tighter):
    assert cn.is_tighter_or_equal(new, PASSING) is tighter


def test_zwidth_maps_a_constant_to_the_midpoint_and_missing_or_zero_to_the_floor():
    assert list(cn.zwidth(pd.Series([10.0, 10.0, 10.0]), (1.0, 3.0))) == [2.0, 2.0, 2.0]

    widths = cn.zwidth(pd.Series([np.nan, 0.0, 10.0, 1000.0]), (1.0, 3.0))
    assert list(widths[:2]) == [1.0, 1.0]
    assert 1.0 <= widths[2] < widths[3] <= 3.0


# --- questions about the drawn tables ------------------------------------------------------

@pytest.fixture
def drawn(scores, tmp_path):
    """BaitA -> P1, P2; BaitB -> P2, PA(=BaitA); literature P1-P2."""
    return cn.build(scores, PASSING, biogrid_path=_biogrid(tmp_path, [('P1', 'P2', 5, True)]))


def test_degrees_count_visible_edges_of_every_kind(drawn):
    nodes, edges = drawn

    assert cn.degrees(edges)['P2'] == 3        # PA, PB and the literature pair
    edges.loc[edges['interaction'] == 'literature', 'visible'] = False
    assert cn.degrees(edges)['P2'] == 2


def test_a_loner_has_the_bait_as_its_only_neighbour(scores, tmp_path):
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'P4', 'G4')
    nodes, edges = cn.build(scores, PASSING, biogrid_path=_biogrid(tmp_path, [('P1', 'P2')]))

    # P1 has a literature partner, P2 a second bait, P4 nothing else.
    assert cn.loners(nodes, edges, 'BaitA') == ['P4']
    edges.loc[edges['interaction'] == 'literature', 'visible'] = False
    assert cn.loners(nodes, edges, 'BaitA') == ['P1', 'P4']


def test_satellites_are_only_baits_preys_plus_nearer_shared_ones(drawn):
    nodes, edges = drawn
    near_a = {'BaitA': (0, 0), 'BaitB': (100, 0), 'P1': (0, 10), 'P2': (10, 0)}
    near_b = {**near_a, 'P2': (90, 0)}

    assert cn.satellites(nodes, edges, 'BaitA', near_a) == (['P1'], ['P2'])
    assert cn.satellites(nodes, edges, 'BaitA', near_b) == (['P1'], [])


def test_a_seed_resolves_by_symbol_or_accession_case_insensitively(drawn):
    nodes, _ = drawn

    assert cn.resolve_node(nodes, 'baita')['id'] == 'BaitA'
    assert cn.resolve_node(nodes, 'p1')['symbol'] == 'G1'
    with pytest.raises(ValueError):
        cn.resolve_node(nodes, 'nobody')


def test_interactors_honour_the_cuts(drawn):
    nodes, edges = drawn

    assert cn.related(nodes, edges, 'BaitA', 'interactors') == ['P1', 'P2']
    assert cn.related(nodes, edges, 'BaitA', 'interactors', min_abundance=1e6) == ['P2']
    assert cn.related(nodes, edges, 'BaitA', 'interactors', min_saint=0.95) == []


def test_singletons_are_the_preys_with_no_other_bait(drawn):
    nodes, edges = drawn

    assert cn.related(nodes, edges, 'BaitA', 'singletons') == ['P1']
    assert cn.related(nodes, edges, 'BaitB', 'singletons') == ['BaitA']


def test_partners_are_reference_neighbours_above_the_publication_floor(drawn):
    nodes, edges = drawn

    assert cn.related(nodes, edges, 'G1', 'partners') == ['P2']
    assert cn.related(nodes, edges, 'G1', 'partners', min_publications=10) == []


def test_cocomplex_follows_the_complex_layer(scores, complexes):
    nodes, edges = cn.build(scores, PASSING, prey_prey=False, corum_path=complexes)

    assert cn.related(nodes, edges, 'G1', 'cocomplex') == ['BaitA', 'P2']


# --- clustering ------------------------------------------------------------------------------

def _two_cliques():
    """Two four-cliques of literature edges joined by one weak proximity edge."""
    rows = []
    for group in ('A', 'B'):
        members = [f'{group}{k}' for k in range(4)]
        for s, t in itertools.combinations(members, 2):
            rows.append({'source': s, 'target': t, 'interaction': 'literature', 'n_publications': 20})
    rows.append({'source': 'A0', 'target': 'B0', 'interaction': 'proximity', 'AvgIntensity': 1.0})
    edges = pd.DataFrame(rows)
    edges['visible'] = True
    return edges


def test_two_cliques_joined_by_a_weak_edge_split_in_two():
    edges = _two_cliques()

    membership = cn.cluster(edges, list(edges['source']) + list(edges['target']), seed=3)

    assert membership.nunique() == 2
    assert membership['A0'] == membership['A3'] and membership['B0'] == membership['B3']
    assert membership['A0'] != membership['B0']


def test_zero_literature_weight_leaves_only_the_proximity_edge():
    edges = _two_cliques()
    ids = sorted(set(edges['source']) | set(edges['target']))

    membership = cn.cluster(edges, ids, literature_weight=0.0)

    # With no reference edges, the six untouched nodes are singletons and A0-B0 pair up.
    assert membership['A0'] == membership['B0']
    assert membership.nunique() == 7


def test_packing_keeps_communities_apart_and_members_on_one_circle():
    membership = pd.Series([0, 0, 0, 1, 1, 2], index=['a', 'b', 'c', 'd', 'e', 'f'])
    positions = {i: (10.0 * k + 100, 200.0) for k, i in enumerate(membership.index)}

    placed = cn.pack_communities(membership, positions)

    assert set(placed) == set(membership.index)
    assert all(x >= 100 and y >= 200 for x, y in placed.values())
    centres = {}
    for c in (0, 1, 2):
        members = [i for i in membership.index if membership[i] == c]
        pts = np.array([placed[i] for i in members])
        centre = pts.mean(axis=0)
        radii = np.linalg.norm(pts - centre, axis=1)
        assert np.allclose(radii, radii[0])
        centres[c] = (centre, radii[0])
    for a, b in itertools.combinations(centres, 2):
        (ca, ra), (cb, rb) = centres[a], centres[b]
        assert np.linalg.norm(ca - cb) >= ra + rb


# =============================================================================================
# cytoscape_ctl: the controller over a stubbed CyREST
# =============================================================================================

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
def controlled(scores, tmp_path, cytoscape):
    """The state after a draw of BaitA -> P1, P2, P4 and BaitB -> P2, PA with a literature
    edge P1-P2, without touching Cytoscape."""
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'P4', 'G4')
    biogrid = _biogrid(tmp_path, [('P1', 'P2', 2, False)])
    nodes, edges = cn.build(scores, PASSING, biogrid_path=biogrid)
    ctl.STATE.update(dataset='d', title='ProxiMate: d', net_suid=1, nodes=nodes, edges=edges,
                     thresholds=dict(PASSING),
                     style={'width_source': 'abundance', 'literature_weighted': False, 'biogrid_scope': 'all'})
    yield nodes, edges
    ctl.STATE.update(dataset=None, title=None, net_suid=None, nodes=None, edges=None,
                     thresholds=None, style=None)


def _last(calls, kind):
    return [payload for name, *payload in calls if name == kind][-1]


def test_restyling_pushes_width_and_visibility_without_moving_anything(controlled, cytoscape):
    changed = ctl.restyle_edges('uniform', False, 'multivalidated')

    frame = _last(cytoscape['calls'], 'edges')[0]
    assert changed == len(frame) > 0
    assert set(frame.columns) == {'name', 'visible', 'width'}
    assert not frame.loc[frame['name'] == 'P1 (literature) P2', 'visible'].item()
    assert ctl.STATE['edges']['width'][ctl.STATE['edges']['interaction'] == 'proximity'].nunique() == 1
    assert ctl.snapshot()['style'] == {'width_source': 'uniform', 'literature_weighted': False,
                                       'biogrid_scope': 'multivalidated'}
    assert not any(name == 'positions' for name, *_ in cytoscape['calls'])


def test_rethresholding_keeps_the_drawn_biogrid_scope(controlled, cytoscape):
    ctl.restyle_edges('abundance', False, 'multivalidated')

    ctl.apply_thresholds(PASSING)

    edges = ctl.STATE['edges']
    assert not edges.loc[edges['interaction'] == 'literature', 'visible'].any()


def test_loners_are_selected_with_their_bait(controlled, cytoscape):
    cytoscape['selected'] = ['BaitA', 'P2']

    chosen = ctl.select_loners()

    assert chosen == ['BaitA', 'P4']
    assert _last(cytoscape['calls'], 'select') == [['BaitA', 'P4'], False]


def test_related_selection_can_add_to_the_current_one(controlled, cytoscape):
    chosen = ctl.select_related('BaitB', 'interactors', add=True, min_saint=0.5)

    assert chosen == ['BaitA', 'P2']
    assert _last(cytoscape['calls'], 'select') == [['BaitA', 'P2'], True]


def test_satellites_are_a_relation_with_an_explicit_seed(controlled, cytoscape):
    cytoscape['positions'] = {'BaitA': (0, 0), 'BaitB': (100, 0), 'P1': (0, 10), 'P2': (10, 0),
                              'P4': (0, 20), 'PA': (100, 10)}

    assert ctl.select_related('BaitA', 'satellites') == ['P1', 'P4', 'P2']
    assert ctl.select_related('baita', 'satellites', include_seed=True) == ['BaitA', 'P1', 'P4', 'P2']
    assert _last(cytoscape['calls'], 'select') == [['BaitA', 'P1', 'P4', 'P2'], False]
    with pytest.raises(ValueError, match='needs a bait'):
        ctl.select_related('P1', 'satellites')


def test_an_empty_relation_is_an_error_not_a_silent_deselect(controlled, cytoscape):
    with pytest.raises(ValueError):
        ctl.select_related('BaitB', 'interactors', min_saint=0.99)
    assert not any(name == 'select' for name, *_ in cytoscape['calls'])


@pytest.fixture
def clusterable(scores, tmp_path, cytoscape):
    """``controlled`` with a fourth and fifth prey of BaitA, enough to cluster without
    the baits."""
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'P4', 'G4')
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'P5', 'G5')
    biogrid = _biogrid(tmp_path, [('P1', 'P2', 2, False)])
    nodes, edges = cn.build(scores, PASSING, biogrid_path=biogrid)
    ctl.STATE.update(dataset='d', title='ProxiMate: d', net_suid=1, nodes=nodes, edges=edges,
                     thresholds=dict(PASSING),
                     style={'width_source': 'abundance', 'literature_weighted': False, 'biogrid_scope': 'all'})
    cytoscape['positions'] = {i: (10.0 * k, 50.0) for k, i in enumerate(['BaitA', 'P1', 'P2', 'P4', 'P5', 'BaitB'])}
    yield nodes, edges
    ctl.STATE.update(dataset=None, title=None, net_suid=None, nodes=None, edges=None,
                     thresholds=None, style=None)


def test_clustering_recolours_numbers_and_moves_only_the_selected_preys(clusterable, cytoscape):
    cytoscape['selected'] = ['BaitA', 'P1', 'P2', 'P4', 'P5']

    result = ctl.cluster_selection(resolution=1.0, seed=1)

    assert result['n'] == 4 and sum(result['sizes']) == 4 and result['baits_left'] == ['BaitA']
    frame = _last(cytoscape['calls'], 'nodes')[0]
    assert set(frame['name']) == {'P1', 'P2', 'P4', 'P5'}
    assert set(frame.columns) == {'name', 'community', 'fill'}
    moved = _last(cytoscape['calls'], 'positions')[0]
    assert set(moved) == {'P1', 'P2', 'P4', 'P5'}
    assert all(x >= 0 and y >= 50 for x, y in moved.values())
    nodes = ctl.STATE['nodes'].set_index('id')
    assert np.isnan(nodes.loc['BaitA', 'community']) and np.isnan(nodes.loc['BaitB', 'community'])
    assert nodes.loc['P1', 'fill'] in cn.COMMUNITY_FILL
    cytoscape['selected'] = ['BaitA', 'BaitB']
    with pytest.raises(ValueError, match='only baits'):
        ctl.cluster_selection()


def test_a_second_clustering_numbers_above_the_first(clusterable, cytoscape):
    cytoscape['selected'] = ['P1', 'P2', 'P4', 'P5']
    ctl.cluster_selection()
    first = ctl.STATE['nodes'].set_index('id').loc[['P1', 'P2', 'P4', 'P5'], 'community'].max()

    cytoscape['selected'] = ['BaitB', 'P2', 'P4', 'P5', 'P1']
    ctl.cluster_selection()

    second = ctl.STATE['nodes'].set_index('id').loc[['P2', 'P4'], 'community'].min()
    assert second > first


# =============================================================================================
# cytoscape_p4c: the CyREST helpers
# =============================================================================================

class _Response:
    def __init__(self, code=200, content=b'', payload=None):
        self.status_code = code
        self.content = content
        self._payload = payload or {}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise cy.requests.HTTPError(str(self.status_code))

    def json(self):
        return self._payload


def test_the_probe_reports_a_down_desktop_as_an_error_not_an_exception(monkeypatch):
    def refused(*a, **k):
        raise cy.requests.ConnectionError('refused')
    monkeypatch.setattr(cy.requests, 'get', refused)

    result = cy.probe()

    assert result['ok'] is False and 'refused' in result['error']


def test_selecting_clears_first_unless_adding(monkeypatch):
    calls = []
    monkeypatch.setattr(cy.p4c, 'clear_selection', lambda **k: calls.append('clear'))
    monkeypatch.setattr(cy.p4c, 'select_nodes', lambda names, **k: calls.append(('select', names, k['preserve_current_selection'])))

    cy.select_nodes(2, ['a', 'b'])
    cy.select_nodes(2, ['c'], add=True)

    assert calls == ['clear', ('select', ['a', 'b'], False), ('select', ['c'], True)]


def test_positions_are_written_as_plain_view_values(monkeypatch):
    seen = {}
    monkeypatch.setattr(cy, 'node_suids', lambda net: {'a': 11, 'b': 12})
    monkeypatch.setattr(cy, 'view_suid', lambda net: 7)
    monkeypatch.setattr(cy.requests, 'put', lambda url, **k: seen.update(url=url, **k) or _Response())

    n = cy.set_positions(3, {'a': (1.5, 2), 'b': (-3, 4)})

    assert n == 2
    assert seen['url'].endswith('/networks/3/views/7/nodes')
    assert 'bypass' not in seen['url']
    assert {p['visualProperty'] for p in seen['json'][0]['view']} == {'NODE_X_LOCATION', 'NODE_Y_LOCATION'}


def test_unlock_clears_only_the_locks_that_are_held(monkeypatch):
    cleared = []
    monkeypatch.setattr(cy, 'locked_view_properties', lambda net: ['NETWORK_SCALE_FACTOR'])
    monkeypatch.setattr(cy.p4c, 'clear_network_property_bypass',
                        lambda vp, **k: cleared.append(vp))
    assert cy.unlock(3) == ['NETWORK_SCALE_FACTOR']
    assert cleared == ['NETWORK_SCALE_FACTOR']

    monkeypatch.setattr(cy, 'locked_view_properties', lambda net: [])

    def boom(*a, **k):
        raise AssertionError("cleared a lock that was not held")
    monkeypatch.setattr(cy.p4c, 'clear_network_property_bypass', boom)
    assert cy.unlock(3) == []


def test_proximate_networks_are_told_apart_by_title(monkeypatch):
    monkeypatch.setattr(cy.p4c, 'get_network_list', lambda: ['ProxiMate: a', 'mine', 'ProxiMate: b'])
    assert cy.proximate_networks() == ['ProxiMate: a', 'ProxiMate: b']
    monkeypatch.setattr(cy.p4c, 'get_network_count', lambda: 0)
    assert cy.current_network_title() is None
    monkeypatch.setattr(cy.p4c, 'get_network_count', lambda: 2)
    monkeypatch.setattr(cy.p4c, 'get_network_name', lambda: 'mine')
    assert cy.current_network_title() == 'mine'


def test_draw_removes_every_earlier_proximate_network(scores, tmp_path, cytoscape, monkeypatch):
    deleted = []
    monkeypatch.setattr(ctl.cy, 'proximate_networks', lambda: ['ProxiMate: old', 'ProxiMate: d'])
    monkeypatch.setattr(ctl.p4c, 'get_network_suid', lambda title: {'ProxiMate: old': 7, 'ProxiMate: d': 8}[title])
    monkeypatch.setattr(ctl.p4c, 'delete_network', lambda suid: deleted.append(suid))
    monkeypatch.setattr(ctl.p4c, 'create_network_from_data_frames', lambda *a, **k: 9)
    monkeypatch.setattr(ctl.cy, 'apply_passthrough_style', lambda *a: None)
    monkeypatch.setattr(ctl.p4c, 'layout_network', lambda *a, **k: None)
    monkeypatch.setattr(ctl.cy, 'unlock', lambda net: [])
    monkeypatch.setattr(ctl.p4c, 'fit_content', lambda **k: None)
    scores.to_csv(tmp_path / 's.csv', index=False)
    try:
        snap = ctl.draw('d', str(tmp_path / 's.csv'), PASSING, biogrid_path=_biogrid(tmp_path, [('P1', 'P2')]))
        assert deleted == [7, 8] and snap['net_suid'] == 9 and snap['title'] == 'ProxiMate: d'
    finally:
        ctl.STATE.update(dataset=None, title=None, net_suid=None, nodes=None, edges=None,
                         thresholds=None, style=None)


def test_render_unlocks_and_fits_before_capturing(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(cy, 'unlock', lambda net: calls.append('unlock') or ['NETWORK_SCALE_FACTOR'])
    monkeypatch.setattr(cy.p4c, 'fit_content', lambda **k: calls.append('fit'))
    monkeypatch.setattr(cy.requests, 'get',
                        lambda url, **k: calls.append(('png', k['params'])) or _Response(content=b'PNG'))

    assert cy.render_png(3, height=600) == b'PNG'
    assert calls == ['unlock', 'fit', ('png', {'h': 600})]
    path = cy.export_png(3, str(tmp_path / 'net.png'), height=600)
    assert open(path, 'rb').read() == b'PNG'


# --- actor, explicit selection and positions (the MCP surface) --------------------------

def test_mutations_record_their_actor(controlled, cytoscape):
    ctl.select_nodes(['P1', 'P2'], actor='mcp')
    ctl.set_edge_visibility('show_all')
    log = ctl.snapshot()['log']
    assert [e['actor'] for e in log[-2:]] == ['mcp', 'gui']
    assert log[-2]['op'] == 'select_nodes'


def test_select_nodes_resolves_symbols_and_refuses_unknown_ones(controlled, cytoscape):
    chosen = ctl.select_nodes(['G1', 'P2'], add=True)
    assert chosen == ['P1', 'P2']
    assert cytoscape['calls'][-1] == ('select', ['P1', 'P2'], True)
    with pytest.raises(ValueError, match='ZZ'):
        ctl.select_nodes(['P1', 'ZZ'])
    with pytest.raises(ValueError, match='clear_selection'):
        ctl.select_nodes([])


def test_get_positions_reads_cytoscape_for_the_drawn_nodes(controlled, cytoscape):
    nodes, _ = controlled
    cytoscape['positions'] = {i: (1.0, 2.0) for i in nodes['id']}
    cytoscape['positions'].update({'P2': (3.0, 4.0), 'stray': (9.0, 9.0)})
    assert set(ctl.get_positions()) == set(nodes['id'])
    assert ctl.get_positions(['G2']) == {'P2': [3.0, 4.0]}
    del cytoscape['positions']['BaitA']
    with pytest.raises(ValueError, match='BaitA'):
        ctl.get_positions(['BaitA'])       # drawn, but Cytoscape reported no position


def test_move_nodes_writes_positions_and_keeps_them_in_the_state(controlled, cytoscape):
    nodes, _ = controlled
    n = ctl.move_nodes({'G1': [10, 20], 'P2': (30, 40)}, actor='mcp')
    assert n == 2
    assert cytoscape['calls'][-1] == ('positions', {'P1': (10.0, 20.0), 'P2': (30.0, 40.0)})
    assert nodes.set_index('id').loc['P1', ['x', 'y']].tolist() == [10.0, 20.0]
    with pytest.raises(ValueError):
        ctl.move_nodes({'P1': [1]})
