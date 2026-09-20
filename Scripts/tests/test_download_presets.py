"""Tests for GUI/download_presets.py: the Downloads tab preset registry and
the pure export-builder functions behind it."""

import os

import numpy as np
import pandas as pd
import pytest

from download_presets import (DEFAULT_CUSTOM_COLUMNS, PRESETS, ColumnGroup,
                              Preset, available_presets, build_annotated_table,
                              build_custom_table, build_cytoscape_edges,
                              build_enrichment_table, build_gene_list,
                              build_prohits_table, effective_selection,
                              resolve_columns, saint_input_files,
                              write_gene_list)

# Thresholds where every score family participates; row 3 has WDFDR NaN.
THRESHOLDS = {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0}


@pytest.fixture
def human_scores():
    """Small annotated_scores frame: 2 baits x 4 preys, one WDFDR NaN row,
    human-only columns (Main location) present."""
    return pd.DataFrame({
        'Prey.ID': ['P1;P1b', 'P2', 'P3', 'P4', 'P1;P1b', 'P2', 'P5', 'P6'],
        'PreyGene': ['G1;G1b', 'G2', 'G3', 'G4', 'G1;G1b', 'G2', 'G5', 'G6'],
        'Experiment.ID': ['BaitA'] * 4 + ['BaitB'] * 4,
        'Bait.ID': ['BA'] * 4 + ['BB'] * 4,
        'First_ID': ['P1', 'P2', 'P3', 'P4', 'P1', 'P2', 'P5', 'P6'],
        'First_Prey_Gene': ['G1', 'G2', 'G3', 'G4', 'G1', 'G2', 'G5', 'G6'],
        'AvgIntensity': [100.0, 200, 300, 400, 150, 250, 350, 450],
        'AvePSM': [10.0, 20, 30, 40, 15, 25, 35, 45],
        'SaintScore': [0.95, 0.8, 0.3, 0.99, 0.9, 0.1, 0.85, 0.7],
        'BFDR': [0.0, 0.01, 0.5, 0.0, 0.02, 0.8, 0.01, 0.04],
        'FoldChange': [5.0, 4, 1, 6, 5, 0.5, 4, 3],
        'WD': [2.0, 1.5, 0.1, 3.0, 2.5, 0.05, 1.8, 1.2],
        'WDFDR': [0.01, 0.02, 0.9, np.nan, 0.01, 0.95, 0.03, 0.02],
        'first_SCL': ['Nucleus'] * 8,
        'Main location': ['Nucleoplasm'] * 8,
        'In.BioGRID': [True, False] * 4,
        'Prey_Is_Bait': [False] * 8,
        'GO_CC': ['nucleus'] * 8,
    })


@pytest.fixture
def mouse_scores(human_scores):
    """Same frame without the human-only and BioGRID columns."""
    return human_scores.drop(columns=['Main location', 'In.BioGRID'])


@pytest.fixture
def enrichment_frame():
    return pd.DataFrame({
        'Bait': ['BaitA', 'BaitA', 'BaitB'],
        'Feature': ['nucleus', 'kinase', 'membrane'],
        'Feature_type': ['GO_CC', 'Domains', 'GO_CC'],
        'k': [5, 3, 2], 'n': [50, 50, 40], 'K': [10, 6, 8], 'M': [500, 500, 500],
        'p_value': [0.001, 0.01, 0.2],
        'enrichment': [5.0, 4.0, 1.2],
        'adj_p': [0.01, 0.05, 0.4],
    })


def _touch(dirpath, *names):
    for name in names:
        with open(os.path.join(dirpath, name), 'w', newline='') as fh:
            fh.write("x\n")


# -------------------------------------------------------- available_presets

SAINT_FILES = ('bait.txt', 'prey.txt', 'interaction.txt')


@pytest.mark.parametrize('files, expected', [
    ((), []),
    (('ED.csv',), ['ed']),
    (('ED.csv', *SAINT_FILES), ['ed', 'saint_inputs']),
    (('ED.csv', *SAINT_FILES, 'annotated_scores.csv'),
     ['ed', 'saint_inputs', 'annotated', 'cytoscape', 'genelist', 'prohits', 'custom']),
    (('ED.csv', *SAINT_FILES, 'annotated_scores.csv', 'Feature_enrichment.csv'),
     ['ed', 'saint_inputs', 'enrichment', 'annotated', 'cytoscape', 'genelist',
      'prohits', 'custom']),
])
def test_available_presets_follow_the_files_present(tmp_path, files, expected):
    _touch(tmp_path, *files)
    assert available_presets(str(tmp_path)) == expected


# ---------------------------------------------------------- resolve_columns

def test_resolve_columns_df_order(human_scores):
    cols = resolve_columns(human_scores.columns, PRESETS['annotated'],
                           ['identifiers', 'saint'])
    assert cols == ['Prey.ID', 'PreyGene', 'Experiment.ID', 'Bait.ID',
                    'First_ID', 'First_Prey_Gene', 'AvgIntensity',
                    'SaintScore', 'BFDR', 'FoldChange']


def test_resolve_columns_dedupes_overlap():
    preset = Preset(key='x', label='x', kind='table', required_files=(),
                    uses_thresholds=False,
                    groups=(ColumnGroup('g1', 'G1', ('A', 'B')),
                            ColumnGroup('g2', 'G2', ('B', 'C'))))
    assert resolve_columns(['A', 'B', 'C'], preset, ['g1', 'g2']) == ['A', 'B', 'C']


# ----------------------------------------------------- build_annotated_table

def test_annotated_table_thresholds(human_scores):
    out = build_annotated_table(human_scores, ['identifiers', 'saint'], THRESHOLDS)
    assert len(out) == 6            # rows 2 and 5 fail; NaN WDFDR passes at 1.0
    assert 'G3' not in out['First_Prey_Gene'].values
    assert list(out.columns)[0] == 'Prey.ID'


def test_annotated_table_no_thresholds(human_scores):
    out = build_annotated_table(human_scores, ['identifiers'], None)
    assert len(out) == 8


def test_annotated_table_no_columns_raises(mouse_scores):
    """A group whose columns are all absent from the dataset resolves to nothing."""
    with pytest.raises(ValueError):
        build_annotated_table(mouse_scores, ['biogrid'], None)


# ---------------------------------------------------- build_enrichment_table

def test_enrichment_table_projection_keeps_rows(enrichment_frame):
    out = build_enrichment_table(enrichment_frame, ['identifiers', 'statistics'])
    assert list(out.columns) == ['Bait', 'Feature', 'Feature_type',
                                 'p_value', 'enrichment', 'adj_p']
    assert len(out) == len(enrichment_frame)


# ---------------------------------------------------- build_cytoscape_edges

def test_cytoscape_edges_shape(human_scores, mouse_scores):
    out = build_cytoscape_edges(human_scores, THRESHOLDS)
    assert list(out.columns)[:2] == ['source', 'target']
    assert set(out['source']) == {'BaitA', 'BaitB'}
    assert 'G1' in out['target'].values
    assert len(out) == 6
    assert 'SaintScore' in out.columns and 'In.BioGRID' in out.columns

    # Attributes are limited to the columns the dataset carries.
    assert 'In.BioGRID' not in build_cytoscape_edges(mouse_scores, THRESHOLDS).columns


def test_cytoscape_edges_preygene_fallback(human_scores):
    df = human_scores.drop(columns=['First_Prey_Gene'])
    out = build_cytoscape_edges(df, THRESHOLDS)
    assert 'G1;G1b' in out['target'].values


# --------------------------------------------------------- build_gene_list

def test_gene_list_pooled(human_scores):
    out = build_gene_list(human_scores, THRESHOLDS, mode='pooled')
    assert list(out.columns) == ['Gene']
    # G1 passes under both baits but appears once; G3 fails thresholds
    assert list(out['Gene']) == ['G1', 'G2', 'G4', 'G5', 'G6']


def test_gene_list_per_bait(human_scores):
    out = build_gene_list(human_scores, THRESHOLDS, mode='per_bait')
    assert list(out.columns) == ['Bait', 'Gene']
    assert list(out[out['Bait'] == 'BaitA']['Gene']) == ['G1', 'G2', 'G4']
    assert list(out[out['Bait'] == 'BaitB']['Gene']) == ['G1', 'G5', 'G6']


def test_gene_list_drops_missing_genes(human_scores):
    df = human_scores.copy()
    df.loc[0, 'First_Prey_Gene'] = np.nan
    df.loc[4, 'First_Prey_Gene'] = np.nan
    out = build_gene_list(df, THRESHOLDS, mode='pooled')
    assert 'G1' not in out['Gene'].values


# ------------------------------------------------------ build_prohits_table

def test_prohits_headers_and_abundance_column(human_scores):
    out = build_prohits_table(human_scores, THRESHOLDS)
    assert list(out.columns) == ['Bait', 'Prey', 'PreyGene', 'Abundance',
                                 'SaintScore', 'BFDR']
    assert out['Abundance'].tolist() == human_scores.loc[
        [0, 1, 3, 4, 6, 7], 'AvePSM'].tolist()

    out = build_prohits_table(human_scores, THRESHOLDS, abundance_col='AvgIntensity')
    assert out['Abundance'].iloc[0] == 100.0


# ------------------------------------------------------- build_custom_table

def test_custom_table_reports_missing_columns(human_scores):
    out, missing = build_custom_table(human_scores, ['SaintScore', 'NotACol'],
                                      THRESHOLDS)
    assert list(out.columns) == ['SaintScore']
    assert missing == ['NotACol']


def test_custom_table_empty_selection_uses_defaults(human_scores):
    out, missing = build_custom_table(human_scores, [], THRESHOLDS)
    assert list(out.columns) == DEFAULT_CUSTOM_COLUMNS
    assert missing == []


# --------------------------------------- saint_input_files / writers

def test_saint_input_files_reports_missing(tmp_path):
    _touch(tmp_path, 'bait.txt', 'prey.txt', 'interaction.txt')
    paths, missing = saint_input_files(
        str(tmp_path), ['bait', 'prey', 'interaction', 'filtered_interaction'])
    assert [os.path.basename(p) for p in paths] == ['bait.txt', 'prey.txt',
                                                    'interaction.txt']
    assert all(os.path.isabs(p) for p in paths)
    assert missing == ['filtered_interaction']


@pytest.mark.parametrize('frame, expected', [
    (pd.DataFrame({'Gene': ['G1', 'G2']}), "G1\nG2\n"),
    (pd.DataFrame({'Bait': ['A', 'B'], 'Gene': ['G1', 'G2']}), "A\tG1\nB\tG2\n"),
])
def test_write_gene_list(tmp_path, frame, expected):
    path = str(tmp_path / 'genes.txt')
    write_gene_list(frame, path)
    with open(path) as fh:
        assert fh.read() == expected


# ------------------------------------------------------ effective_selection

def test_effective_selection_keeps_valid_keys():
    assert effective_selection(PRESETS['annotated'],
                               ['identifiers', 'bait']) == ['identifiers']


@pytest.mark.parametrize('preset, selected', [
    ('saint_inputs', ['identifiers', 'saint']),   # keys of another preset
    ('annotated', []),
])
def test_effective_selection_falls_back_to_defaults(preset, selected):
    defaults = {'saint_inputs': ['bait', 'prey', 'interaction'],
                'annotated': ['identifiers', 'saint', 'comppass']}
    assert effective_selection(PRESETS[preset], selected) == defaults[preset]


def test_effective_selection_preset_without_groups():
    assert effective_selection(PRESETS['cytoscape'], ['identifiers']) == []
