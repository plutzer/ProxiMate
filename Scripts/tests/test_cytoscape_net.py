"""Tests for building the node and edge tables the Cytoscape tab draws."""

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
    path = tmp_path / 'biogrid_summary.csv'
    pd.DataFrame([{'SWISS-PROT Accessions Interactor A': a,
                   'SWISS-PROT Accessions Interactor B': b} for a, b in pairs]).to_csv(path, index=False)
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
