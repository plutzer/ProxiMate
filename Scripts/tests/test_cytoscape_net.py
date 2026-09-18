"""Tests for building the node and edge tables the Cytoscape tab draws."""

import itertools

import numpy as np
import pandas as pd
import pytest

import cytoscape_net as cn


PASSING = {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0}


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


# --- nodes --------------------------------------------------------------------------------

def test_nodes_are_keyed_on_accession_with_baits_as_diamonds(scores):
    nodes, _ = cn.build(scores, PASSING, prey_prey=False)

    by_id = nodes.set_index('id')
    assert by_id.loc['PA', 'role'] == 'bait' and by_id.loc['PA', 'shape'] == 'DIAMOND'
    assert by_id.loc['P1', 'role'] == 'prey' and by_id.loc['P1', 'symbol'] == 'G1'
    assert 'P3' not in by_id.index


def test_a_prey_that_is_also_a_bait_keeps_one_bait_node(scores):
    nodes, _ = cn.build(scores, PASSING, prey_prey=False)

    assert (nodes['id'] == 'PA').sum() == 1
    assert nodes.set_index('id').loc['PA', 'symbol'] == 'BaitA'


def test_the_bait_filter_restricts_both_tables(scores):
    nodes, edges = cn.build(scores, PASSING, baits=['BaitB'], prey_prey=False)

    assert set(nodes['id']) == {'PB', 'P2', 'PA'}
    assert set(edges['source']) == {'PB'}


def test_bait_id_is_used_when_no_accession_column_exists(scores):
    df = scores.rename(columns={'Bait_Accession': 'Bait.ID'})

    nodes, _ = cn.build(df, PASSING, prey_prey=False)

    assert 'PA' in set(nodes['id'])


@pytest.mark.parametrize('policy, expected', [
    ('all', {'BaitA', 'G1', 'G2'}), ('baits', {'BaitA', ''}), ('none', {''})])
def test_label_policies(scores, policy, expected):
    nodes, _ = cn.build(scores, PASSING, baits=['BaitA'], prey_prey=False, label_policy=policy)

    assert set(nodes['display_label']) == expected


def test_an_unknown_label_policy_raises(scores):
    with pytest.raises(ValueError, match='label policy'):
        cn.build(scores, PASSING, prey_prey=False, label_policy='some')


# --- edges --------------------------------------------------------------------------------

def test_only_passing_interactions_become_edges_and_carry_their_scores(scores):
    _, edges = cn.build(scores, PASSING, prey_prey=False)

    assert len(edges) == 4
    assert set(edges['interaction']) == {'proximity'}
    assert {'SaintScore', 'BFDR', 'FoldChange', 'WD', 'WDFDR', 'AvgIntensity', 'In.BioGRID'} <= set(edges.columns)
    assert edges['visible'].all()
    assert edges['name'].iloc[0] == 'PA (proximity) P1'


def test_a_missing_wdfdr_fails_the_thresholds(scores):
    scores.loc[0, 'WDFDR'] = np.nan

    _, edges = cn.build(scores, {**PASSING, "WDFDR": 0.5}, prey_prey=False)

    assert 'P1' not in set(edges['target'])


def test_widths_follow_log_abundance(scores):
    _, edges = cn.build(scores, PASSING, prey_prey=False)

    by_target = edges.set_index(['source', 'target'])['width']
    assert by_target[('PA', 'P1')] < by_target[('PB', 'P2')] < by_target[('PA', 'P2')]


def test_spectral_counts_drive_widths_when_intensity_is_absent(scores):
    df = scores.rename(columns={'AvgIntensity': 'AvgSpec'})

    _, edges = cn.build(df, PASSING, prey_prey=False)

    assert edges['width'].nunique() > 1


def test_no_abundance_column_gives_one_width(scores):
    _, edges = cn.build(scores.drop(columns=['AvgIntensity']), PASSING, prey_prey=False)

    assert edges['width'].nunique() == 1


@pytest.mark.parametrize('source, thinner, thicker', [
    ('SaintScore', ('PB', 'PA'), ('PA', 'P1')),      # 0.9 for both PA rows, PB->PA row is 0.9 too
    ('WD', ('PA', 'P1'), ('PB', 'P2')),
    ('FoldChange', ('PA', 'P1'), ('PB', 'P2')),
])
def test_other_width_sources_follow_their_column(scores, source, thinner, thicker):
    scores.loc[scores['First_ID'] == 'P2', 'WD'] = 9.0
    scores.loc[scores['First_ID'] == 'P2', 'FoldChange'] = 30.0
    scores.loc[scores['First_ID'] == 'P1', 'SaintScore'] = 1.0
    scores.loc[scores['First_ID'] == 'PA', 'SaintScore'] = 0.7

    _, edges = cn.build(scores, PASSING, prey_prey=False, width_source=source)

    by_pair = edges.set_index(['source', 'target'])['width']
    assert by_pair[thinner] < by_pair[thicker]


def test_saint_widths_are_linear_over_the_band(scores):
    scores['SaintScore'] = [1.0, 0.75, 0.2, 0.75, 1.0]

    _, edges = cn.build(scores, PASSING, prey_prey=False, width_source='SaintScore')

    lo, hi = cn.PROXIMITY_WIDTH
    assert set(edges['width'].round(3)) == {hi, round(lo + 0.75 * (hi - lo), 3)}


def test_a_uniform_width_source_gives_one_width(scores):
    _, edges = cn.build(scores, PASSING, prey_prey=False, width_source='uniform')

    assert edges['width'].nunique() == 1


def test_an_unknown_width_source_raises(scores):
    with pytest.raises(ValueError, match='width source'):
        cn.build(scores, PASSING, prey_prey=False, width_source='Entropy')


def test_a_width_source_without_its_column_names_it(scores):
    with pytest.raises(KeyError, match='FoldChange'):
        cn.build(scores.drop(columns=['FoldChange']), PASSING, prey_prey=False, width_source='FoldChange')


def test_a_missing_required_column_names_itself(scores):
    with pytest.raises(KeyError, match='BFDR'):
        cn.build(scores.drop(columns=['BFDR']), PASSING)


def test_nothing_passing_raises(scores):
    with pytest.raises(ValueError, match='no interactions pass'):
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

    pairs = cn.prey_prey_pairs({'P1', 'P2'}, {'PA', 'PB'}, biogrid)

    assert len(pairs) == 1
    assert pairs['n_publications'].iloc[0] == 3 and bool(pairs['multivalidated'].iloc[0])


def test_unweighted_literature_edges_sit_at_the_band_floor(scores, tmp_path):
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'P4', 'G4')
    biogrid = _biogrid(tmp_path, [('P1', 'P2', 40, False), ('P2', 'P4', 1, False)])

    _, edges = cn.build(scores, PASSING, biogrid_path=biogrid)

    lit = edges[edges['interaction'] == 'literature']
    assert set(lit['width']) == {cn.LITERATURE_WIDTH[0]}


def test_weighted_literature_edges_thicken_with_publications(scores, tmp_path):
    scores.loc[len(scores)] = _row('BaitA', 'PA', 'P4', 'G4')
    biogrid = _biogrid(tmp_path, [('P1', 'P2', 40, False), ('P2', 'P4', 1, False)])

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


def test_an_unknown_biogrid_scope_raises(scores):
    _, edges = cn.build(scores, PASSING, prey_prey=False)

    with pytest.raises(ValueError, match='BioGRID scope'):
        cn.edge_visibility(edges, PASSING, biogrid_scope='some')


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
    assert set(complex_edges['name']) == {'P1 (complex) P2', 'P1 (complex) P4', 'P2 (complex) P4', 'P4 (complex) PA'}
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
    faded = cn.node_alpha(pd.DataFrame({'id': ['PA', 'P1']}), edges)
    assert list(faded) == [cn.NODE_ALPHA['faded']] * 2


def test_the_literature_edge_survives_while_both_preys_keep_an_edge(scores, tmp_path):
    biogrid = _biogrid(tmp_path, [('P1', 'P2')])
    _, edges = cn.build(scores, PASSING, biogrid_path=biogrid)

    visible = cn.edge_visibility(edges, PASSING)

    assert visible.all()


@pytest.mark.parametrize('new, tighter', [
    ({**PASSING, 'SaintScore': 0.9}, True),
    ({**PASSING, 'BFDR': 0.1}, False),
    (PASSING, True),
])
def test_only_tighter_or_equal_thresholds_can_be_applied_in_place(new, tighter):
    assert cn.is_tighter_or_equal(new, PASSING) is tighter


# --- zwidth ------------------------------------------------------------------------------------

def test_a_constant_column_maps_to_the_band_midpoint():
    widths = cn.zwidth(pd.Series([10.0, 10.0, 10.0]), (1.0, 3.0))

    assert list(widths) == [2.0, 2.0, 2.0]


def test_missing_and_zero_values_map_to_the_band_floor():
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
    assert cn.loners(nodes, edges, 'PA') == ['P4']
    edges.loc[edges['interaction'] == 'literature', 'visible'] = False
    assert cn.loners(nodes, edges, 'PA') == ['P1', 'P4']


def test_satellites_are_only_baits_preys_plus_nearer_shared_ones(drawn):
    nodes, edges = drawn
    near_a = {'PA': (0, 0), 'PB': (100, 0), 'P1': (0, 10), 'P2': (10, 0)}
    near_b = {**near_a, 'P2': (90, 0)}

    assert cn.satellites(nodes, edges, 'PA', near_a) == (['P1'], ['P2'])
    assert cn.satellites(nodes, edges, 'PA', near_b) == (['P1'], [])


def test_a_seed_resolves_by_symbol_or_accession_case_insensitively(drawn):
    nodes, _ = drawn

    assert cn.resolve_node(nodes, 'baita')['id'] == 'PA'
    assert cn.resolve_node(nodes, 'p1')['symbol'] == 'G1'
    with pytest.raises(ValueError, match='not in the drawn network'):
        cn.resolve_node(nodes, 'nobody')


def test_interactors_honour_the_cuts(drawn):
    nodes, edges = drawn

    assert cn.related(nodes, edges, 'BaitA', 'interactors') == ['P1', 'P2']
    assert cn.related(nodes, edges, 'BaitA', 'interactors', min_abundance=1e6) == ['P2']
    assert cn.related(nodes, edges, 'BaitA', 'interactors', min_saint=0.95) == []


def test_singletons_are_the_preys_with_no_other_bait(drawn):
    nodes, edges = drawn

    assert cn.related(nodes, edges, 'BaitA', 'singletons') == ['P1']
    assert cn.related(nodes, edges, 'BaitB', 'singletons') == ['PA']


def test_a_prey_seed_cannot_have_interactors(drawn):
    nodes, edges = drawn

    with pytest.raises(ValueError, match='needs a bait'):
        cn.related(nodes, edges, 'G1', 'interactors')


def test_partners_are_reference_neighbours_above_the_publication_floor(drawn):
    nodes, edges = drawn

    assert cn.related(nodes, edges, 'G1', 'partners') == ['P2']
    assert cn.related(nodes, edges, 'G1', 'partners', min_publications=10) == []


def test_cocomplex_needs_the_complex_layer(scores, drawn, complexes):
    nodes, edges = drawn
    with pytest.raises(ValueError, match='without CORUM'):
        cn.related(nodes, edges, 'G1', 'cocomplex')

    nodes, edges = cn.build(scores, PASSING, prey_prey=False, corum_path=complexes)

    assert cn.related(nodes, edges, 'G1', 'cocomplex') == ['P2', 'PA']


def test_an_unknown_relation_raises(drawn):
    nodes, edges = drawn

    with pytest.raises(ValueError, match='relation'):
        cn.related(nodes, edges, 'BaitA', 'friends')


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


def test_clustering_is_reproducible_for_a_seed():
    edges = _two_cliques()
    ids = sorted(set(edges['source']) | set(edges['target']))

    assert cn.cluster(edges, ids, seed=5).equals(cn.cluster(edges, ids, seed=5))


def test_zero_literature_weight_leaves_only_the_proximity_edge():
    edges = _two_cliques()
    ids = sorted(set(edges['source']) | set(edges['target']))

    membership = cn.cluster(edges, ids, literature_weight=0.0)

    # With no reference edges, the six untouched nodes are singletons and A0-B0 pair up.
    assert membership['A0'] == membership['B0']
    assert membership.nunique() == 7


def test_too_small_a_selection_cannot_be_clustered():
    with pytest.raises(ValueError, match='at least 4'):
        cn.cluster(_two_cliques(), ['A0', 'A1', 'A2'])


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


def test_community_fills_cycle_the_palette():
    membership = pd.Series([0, 1, len(cn.COMMUNITY_FILL)], index=['a', 'b', 'c'])

    fills = cn.community_fill(membership)

    assert fills['a'] == fills['c'] == cn.COMMUNITY_FILL[0] and fills['b'] == cn.COMMUNITY_FILL[1]
