"""Tests for GUI/download_presets.py: the Downloads tab preset registry and
the pure export-builder functions behind it."""

import os
import zipfile

import numpy as np
import pandas as pd
import pytest

from download_presets import (DEFAULT_CUSTOM_COLUMNS, PRESETS, ColumnGroup,
                              Preset, available_presets, build_annotated_table,
                              build_custom_table, build_cytoscape_edges,
                              build_enrichment_table, build_gene_list,
                              build_prohits_table, resolve_columns,
                              saint_input_files, write_gene_list, zip_files)

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


# ---------------------------------------------------------------- registry

def test_registry_has_eight_presets_in_order():
    assert list(PRESETS) == ['ed', 'saint_inputs', 'enrichment', 'annotated',
                             'cytoscape', 'genelist', 'prohits', 'custom']


def test_registry_kinds_and_extensions():
    assert PRESETS['saint_inputs'].kind == 'files'
    assert PRESETS['saint_inputs'].extension == '.zip'
    assert PRESETS['genelist'].kind == 'genelist'
    assert PRESETS['genelist'].extension == '.txt'
    for key in ('ed', 'enrichment', 'annotated', 'cytoscape', 'prohits', 'custom'):
        assert PRESETS[key].kind == 'table'
        assert PRESETS[key].extension == '.csv'
    for key in ('ed', 'saint_inputs', 'enrichment'):
        assert not PRESETS[key].uses_thresholds
    for key in ('annotated', 'cytoscape', 'genelist', 'prohits', 'custom'):
        assert PRESETS[key].uses_thresholds


def test_annotated_group_defaults():
    groups = {g.key: g for g in PRESETS['annotated'].groups}
    assert set(groups) == {'identifiers', 'saint', 'comppass', 'uniprot', 'go',
                           'localization', 'biogrid', 'topology'}
    for key in ('identifiers', 'saint', 'comppass'):
        assert groups[key].default
    for key in ('uniprot', 'go', 'localization', 'biogrid', 'topology'):
        assert not groups[key].default


def test_saint_inputs_file_choices():
    files = {f.key: f for f in PRESETS['saint_inputs'].files}
    assert set(files) == {'bait', 'prey', 'interaction', 'filtered_interaction',
                          'imputed_prey', 'imputed_params'}
    for key in ('bait', 'prey', 'interaction'):
        assert files[key].default
    for key in ('filtered_interaction', 'imputed_prey', 'imputed_params'):
        assert not files[key].default


# -------------------------------------------------------- available_presets

def test_available_presets_empty_dir(tmp_path):
    assert available_presets(str(tmp_path)) == []


def test_available_presets_ed_only(tmp_path):
    _touch(tmp_path, 'ED.csv')
    assert available_presets(str(tmp_path)) == ['ed']


def test_available_presets_unscored(tmp_path):
    _touch(tmp_path, 'ED.csv', 'bait.txt', 'prey.txt', 'interaction.txt')
    assert available_presets(str(tmp_path)) == ['ed', 'saint_inputs']


def test_available_presets_scored(tmp_path):
    _touch(tmp_path, 'ED.csv', 'bait.txt', 'prey.txt', 'interaction.txt',
           'annotated_scores.csv')
    assert available_presets(str(tmp_path)) == [
        'ed', 'saint_inputs', 'annotated', 'cytoscape', 'genelist', 'prohits',
        'custom']


def test_available_presets_with_enrichment(tmp_path):
    _touch(tmp_path, 'ED.csv', 'bait.txt', 'prey.txt', 'interaction.txt',
           'annotated_scores.csv', 'Feature_enrichment.csv')
    assert available_presets(str(tmp_path)) == [
        'ed', 'saint_inputs', 'enrichment', 'annotated', 'cytoscape',
        'genelist', 'prohits', 'custom']


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


def test_resolve_columns_absent_group_is_empty(mouse_scores):
    assert resolve_columns(mouse_scores.columns, PRESETS['annotated'],
                           ['biogrid']) == []


def test_resolve_columns_empty_selection(human_scores):
    assert resolve_columns(human_scores.columns, PRESETS['annotated'], []) == []


# ----------------------------------------------------- build_annotated_table

def test_annotated_table_thresholds(human_scores):
    out = build_annotated_table(human_scores, ['identifiers', 'saint'], THRESHOLDS)
    assert len(out) == 6            # rows 2 and 5 fail; NaN WDFDR passes at 1.0
    assert 'G3' not in out['First_Prey_Gene'].values
    assert list(out.columns)[0] == 'Prey.ID'


def test_annotated_table_wdfdr_nan_fails_tight_threshold(human_scores):
    thresholds = {'SaintScore': 0.0, 'BFDR': 1.0, 'WD': 0.0, 'WDFDR': 0.05}
    out = build_annotated_table(human_scores, ['identifiers'], thresholds)
    assert 'G4' not in out['First_Prey_Gene'].values
    assert len(out) == 5


def test_annotated_table_no_thresholds(human_scores):
    out = build_annotated_table(human_scores, ['identifiers'], None)
    assert len(out) == 8


def test_annotated_table_no_columns_raises(mouse_scores):
    with pytest.raises(ValueError, match="column"):
        build_annotated_table(mouse_scores, ['biogrid'], None)


# ---------------------------------------------------- build_enrichment_table

def test_enrichment_table_projection_keeps_rows(enrichment_frame):
    out = build_enrichment_table(enrichment_frame, ['identifiers', 'statistics'])
    assert list(out.columns) == ['Bait', 'Feature', 'Feature_type',
                                 'p_value', 'enrichment', 'adj_p']
    assert len(out) == len(enrichment_frame)


# ---------------------------------------------------- build_cytoscape_edges

def test_cytoscape_edges_shape(human_scores):
    out = build_cytoscape_edges(human_scores, THRESHOLDS)
    assert list(out.columns)[:2] == ['source', 'target']
    assert set(out['source']) == {'BaitA', 'BaitB'}
    assert 'G1' in out['target'].values
    assert len(out) == 6
    assert 'SaintScore' in out.columns and 'In.BioGRID' in out.columns


def test_cytoscape_edges_attrs_limited_to_present(mouse_scores):
    out = build_cytoscape_edges(mouse_scores, THRESHOLDS)
    assert 'In.BioGRID' not in out.columns


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

def test_prohits_headers_and_default_abundance(human_scores):
    out = build_prohits_table(human_scores, THRESHOLDS)
    assert list(out.columns) == ['Bait', 'Prey', 'PreyGene', 'Abundance',
                                 'SaintScore', 'BFDR']
    assert out['Abundance'].tolist() == human_scores.loc[
        [0, 1, 3, 4, 6, 7], 'AvePSM'].tolist()


def test_prohits_abundance_switch(human_scores):
    out = build_prohits_table(human_scores, THRESHOLDS,
                              abundance_col='AvgIntensity')
    assert out['Abundance'].iloc[0] == 100.0


def test_prohits_unknown_abundance_raises(human_scores):
    with pytest.raises(ValueError, match="[Aa]bundance"):
        build_prohits_table(human_scores, THRESHOLDS, abundance_col='Nope')


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


def test_custom_table_applies_thresholds(human_scores):
    out, _ = build_custom_table(human_scores, ['SaintScore'], THRESHOLDS)
    assert len(out) == 6


# --------------------------------------- saint_input_files / writers / zip

def test_saint_input_files_reports_missing(tmp_path):
    _touch(tmp_path, 'bait.txt', 'prey.txt', 'interaction.txt')
    paths, missing = saint_input_files(
        str(tmp_path), ['bait', 'prey', 'interaction', 'filtered_interaction'])
    assert [os.path.basename(p) for p in paths] == ['bait.txt', 'prey.txt',
                                                    'interaction.txt']
    assert all(os.path.isabs(p) for p in paths)
    assert missing == ['filtered_interaction']


def test_zip_files_round_trip(tmp_path):
    _touch(tmp_path, 'a.txt', 'b.txt')
    paths = [str(tmp_path / 'a.txt'), str(tmp_path / 'b.txt')]
    zip_path = str(tmp_path / 'out.zip')
    assert zip_files(paths, zip_path) == zip_path
    with zipfile.ZipFile(zip_path) as zf:
        assert sorted(zf.namelist()) == ['a.txt', 'b.txt']
        assert zf.read('a.txt') == b"x\n"


def test_zip_files_empty_raises(tmp_path):
    with pytest.raises(ValueError):
        zip_files([], str(tmp_path / 'out.zip'))


def test_write_gene_list_pooled(tmp_path):
    path = str(tmp_path / 'genes.txt')
    write_gene_list(pd.DataFrame({'Gene': ['G1', 'G2']}), path)
    with open(path) as fh:
        assert fh.read() == "G1\nG2\n"


def test_write_gene_list_per_bait(tmp_path):
    path = str(tmp_path / 'genes.txt')
    write_gene_list(pd.DataFrame({'Bait': ['A', 'B'], 'Gene': ['G1', 'G2']}),
                    path)
    with open(path) as fh:
        assert fh.read() == "A\tG1\nB\tG2\n"


# ------------------------------------------------------ effective_selection

def test_effective_selection_keeps_valid_keys():
    from download_presets import effective_selection
    assert effective_selection(PRESETS['annotated'],
                               ['identifiers', 'bait']) == ['identifiers']


def test_effective_selection_stale_keys_fall_back_to_defaults():
    from download_presets import effective_selection
    assert effective_selection(PRESETS['saint_inputs'],
                               ['identifiers', 'saint']) == [
        'bait', 'prey', 'interaction']


def test_effective_selection_empty_falls_back_to_defaults():
    from download_presets import effective_selection
    assert effective_selection(PRESETS['annotated'], []) == [
        'identifiers', 'saint', 'comppass']


def test_effective_selection_preset_without_groups():
    from download_presets import effective_selection
    assert effective_selection(PRESETS['cytoscape'], ['identifiers']) == []
