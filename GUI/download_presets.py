"""
Download preset registry and export builders for the Downloads tab.

Pure logic only: functions take DataFrames/paths and return DataFrames, file
lists, or paths, so everything is testable without a Shiny session. app.py owns
the reactive wiring.

Column groups list *candidate* columns; every projection intersects them with
the columns actually present in the dataset, because annotated_scores.csv
varies by organism (HPA/CORUM are human-only) and by pipeline options.
"""

import os
import zipfile
from dataclasses import dataclass
from typing import Sequence

import pandas as pd

from QC_plots import apply_score_thresholds
from log_config import get_logger

logger = get_logger(__name__)

DEFAULT_CUSTOM_COLUMNS = ["Experiment.ID", "Prey.ID", "SaintScore", "BFDR"]


@dataclass(frozen=True)
class ColumnGroup:
    key: str                    # stable id used as checkbox value
    label: str                  # checkbox label
    columns: tuple              # candidate columns, intersected with the df
    default: bool = True        # checked by default


@dataclass(frozen=True)
class FileChoice:
    key: str
    label: str
    filename: str               # relative to the dataset directory
    default: bool = True


@dataclass(frozen=True)
class Preset:
    key: str
    label: str
    kind: str                   # 'table' | 'files' | 'genelist'
    required_files: tuple       # all must exist in the dataset dir
    uses_thresholds: bool
    groups: tuple = ()          # ColumnGroups for table presets with checkboxes
    files: tuple = ()           # FileChoices for 'files' presets
    extension: str = ".csv"


_ANNOTATED_GROUPS = (
    ColumnGroup('identifiers', 'Identifiers',
                ('Prey.ID', 'PreyGene', 'Experiment.ID', 'Bait.ID', 'First_ID',
                 'First_Prey_Gene', 'Matched_Gene_Name', 'Entry', 'Entry Name',
                 'source_group')),
    ColumnGroup('saint', 'SAINT scores',
                ('Intensity', 'IntensitySum', 'AvgIntensity', 'NumReplicates',
                 'ctrlIntensity', 'AvgP', 'MaxP', 'TopoAvgP', 'TopoMaxP',
                 'SaintScore', 'OddsScore', 'FoldChange', 'BFDR')),
    ColumnGroup('comppass', 'CompPASS scores',
                ('AvePSM', 'N_Saw', 'Entropy', 'Mean', 'SD', 'Z', 'WD',
                 'N_Exp_With_Prey', 'WD_pval', 'WDFDR', 'Self.Interaction',
                 'Self.Only')),
    ColumnGroup('uniprot', 'UniProt annotations',
                ('Protein names', 'Gene Names', 'Subcellular location [CC]',
                 'Sequence similarities', 'Zinc finger', 'Motif',
                 'Protein families', 'Region', 'Repeat', 'Coiled coil',
                 'Compositional bias', 'Domain [CC]', 'Domain [FT]',
                 'Involvement in disease', 'Post-translational modification',
                 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains'),
                default=False),
    ColumnGroup('go', 'GO terms',
                ('GO_CC', 'GO_BP', 'GO_MF',
                 'Gene Ontology (cellular component)',
                 'Gene Ontology (biological process)',
                 'Gene Ontology (molecular function)', 'Gene Ontology IDs'),
                default=False),
    ColumnGroup('localization', 'Localization (HPA/CORUM/UniProt)',
                ('first_SCL', 'Main location', 'Gene name', 'Human_Complex'),
                default=False),
    ColumnGroup('biogrid', 'BioGRID evidence',
                ('SWISS-PROT Accessions Interactor A',
                 'SWISS-PROT Accessions Interactor B', 'Experimental System',
                 'Author', 'Publication Source', 'Multivalidated',
                 'In.BioGRID'),
                default=False),
    ColumnGroup('topology', 'Topology / similarity',
                ('Prey_Is_Bait', 'Self-Interaction', 'CCO'),
                default=False),
)

_ENRICHMENT_GROUPS = (
    ColumnGroup('identifiers', 'Identifiers', ('Bait', 'Feature', 'Feature_type')),
    ColumnGroup('counts', 'Counts', ('k', 'n', 'K', 'M')),
    ColumnGroup('statistics', 'Statistics', ('p_value', 'enrichment', 'adj_p')),
)

_SAINT_INPUT_FILES = (
    FileChoice('bait', 'Bait file', 'bait.txt'),
    FileChoice('prey', 'Prey file', 'prey.txt'),
    FileChoice('interaction', 'Interaction file', 'interaction.txt'),
    FileChoice('filtered_interaction', 'Filtered interactions',
               'filtered_interaction.txt', default=False),
    FileChoice('imputed_prey', 'Imputed prey file', 'imputed_prey.txt',
               default=False),
    FileChoice('imputed_params', 'Imputation parameters', 'imputed_params.csv',
               default=False),
)

# Edge attributes offered to Cytoscape, in output order; limited to those
# present in the dataset.
_CYTOSCAPE_ATTRS = ('SaintScore', 'BFDR', 'FoldChange', 'WD', 'WDFDR',
                    'AvgIntensity', 'In.BioGRID')

PRESETS = {p.key: p for p in (
    Preset('ed', 'Experimental Design', 'table', ('ED.csv',),
           uses_thresholds=False),
    Preset('saint_inputs', 'SAINT Inputs', 'files',
           ('bait.txt', 'prey.txt', 'interaction.txt'),
           uses_thresholds=False, files=_SAINT_INPUT_FILES, extension='.zip'),
    Preset('enrichment', 'Enriched Features', 'table',
           ('Feature_enrichment.csv',),
           uses_thresholds=False, groups=_ENRICHMENT_GROUPS),
    Preset('annotated', 'Annotated Scores', 'table', ('annotated_scores.csv',),
           uses_thresholds=True, groups=_ANNOTATED_GROUPS),
    Preset('cytoscape', 'Cytoscape Edge Table', 'table',
           ('annotated_scores.csv',), uses_thresholds=True),
    Preset('genelist', 'Gene List (STRING/g:Profiler)', 'genelist',
           ('annotated_scores.csv',), uses_thresholds=True, extension='.txt'),
    Preset('prohits', 'ProHits-viz Table', 'table', ('annotated_scores.csv',),
           uses_thresholds=True),
    Preset('custom', 'Custom Columns', 'table', ('annotated_scores.csv',),
           uses_thresholds=True),
)}


def available_presets(dataset_dir):
    """Preset keys whose required files all exist in dataset_dir, in registry
    order."""
    return [key for key, preset in PRESETS.items()
            if all(os.path.exists(os.path.join(dataset_dir, f))
                   for f in preset.required_files)]


def effective_selection(preset, selected_keys):
    """Selected keys restricted to the preset's own group/file keys, falling
    back to the preset's defaults when nothing valid remains. The UI's checkbox
    state can briefly belong to a previously shown preset (the update round-trips
    through the browser), so builders must never see foreign keys."""
    items = preset.files if preset.kind == "files" else preset.groups
    valid = {item.key for item in items}
    kept = [k for k in selected_keys if k in valid]
    if kept or not items:
        return kept
    return [item.key for item in items if item.default]


def resolve_columns(df_columns, preset, selected_group_keys):
    """Union of the selected groups' candidate columns, intersected with
    df_columns, deduplicated, in df_columns order."""
    selected = set(selected_group_keys)
    candidates = set()
    for group in preset.groups:
        if group.key in selected:
            candidates.update(group.columns)
    return [c for c in df_columns if c in candidates]


def _apply_thresholds(df, thresholds):
    if thresholds is None:
        return df
    return apply_score_thresholds(df, thresholds)


def build_annotated_table(df, selected_group_keys, thresholds):
    """Threshold-filtered projection of annotated scores onto the selected
    column groups. Raises ValueError when nothing resolves."""
    columns = resolve_columns(df.columns, PRESETS['annotated'],
                              selected_group_keys)
    if not columns:
        raise ValueError(
            "No columns to export: the selected groups match no column in this "
            "dataset. Select at least one applicable group.")
    return _apply_thresholds(df, thresholds)[columns]


def build_enrichment_table(df, selected_group_keys):
    """Projection of Feature_enrichment results onto the selected groups.
    Feature-level data: score thresholds never apply."""
    columns = resolve_columns(df.columns, PRESETS['enrichment'],
                              selected_group_keys)
    if not columns:
        raise ValueError(
            "No columns to export: select at least one column group.")
    return df[columns]


def _target_gene_column(df):
    if 'First_Prey_Gene' in df.columns:
        return 'First_Prey_Gene'
    return 'PreyGene'


def build_cytoscape_edges(df, thresholds):
    """Edge table for Cytoscape's Import Network from Table: source (bait),
    target (prey gene), then the score columns present as edge attributes."""
    filtered = _apply_thresholds(df, thresholds)
    edges = pd.DataFrame({
        'source': filtered['Experiment.ID'],
        'target': filtered[_target_gene_column(filtered)],
    })
    for col in _CYTOSCAPE_ATTRS:
        if col in filtered.columns:
            edges[col] = filtered[col].values
    return edges


def build_gene_list(df, thresholds, mode="pooled"):
    """Prey genes passing the thresholds. 'pooled': single unique sorted Gene
    column; 'per_bait': Bait/Gene rows, deduplicated within bait."""
    filtered = _apply_thresholds(df, thresholds)
    genes = filtered[_target_gene_column(filtered)]
    genes = genes[genes.notna() & (genes.astype(str).str.strip() != "")]
    if mode == "pooled":
        return pd.DataFrame({'Gene': sorted(genes.unique())})
    if mode == "per_bait":
        out = pd.DataFrame({'Bait': filtered.loc[genes.index, 'Experiment.ID'],
                            'Gene': genes})
        out = out.drop_duplicates().sort_values(['Bait', 'Gene'])
        return out.reset_index(drop=True)
    raise ValueError(f"Unknown gene list mode: {mode!r}")


def build_prohits_table(df, thresholds, abundance_col="AvePSM"):
    """Minimal ProHits-viz-style table: Bait, Prey, PreyGene, Abundance,
    SaintScore, BFDR."""
    if abundance_col not in df.columns:
        raise ValueError(
            f"Abundance column {abundance_col!r} not present in this dataset.")
    filtered = _apply_thresholds(df, thresholds)
    return pd.DataFrame({
        'Bait': filtered['Experiment.ID'],
        'Prey': filtered['Prey.ID'],
        'PreyGene': filtered[_target_gene_column(filtered)],
        'Abundance': filtered[abundance_col],
        'SaintScore': filtered['SaintScore'],
        'BFDR': filtered['BFDR'],
    }).reset_index(drop=True)


def build_custom_table(df, columns, thresholds):
    """Threshold-filtered projection onto the requested columns. Columns absent
    from the dataset are reported back, not raised, so a selection made on one
    dataset degrades gracefully on another. Empty selection falls back to
    DEFAULT_CUSTOM_COLUMNS."""
    requested = list(columns) or list(DEFAULT_CUSTOM_COLUMNS)
    present = [c for c in requested if c in df.columns]
    missing = [c for c in requested if c not in df.columns]
    if missing:
        logger.warning("Custom download: columns not in dataset: %s", missing)
    return _apply_thresholds(df, thresholds)[present], missing


def saint_input_files(dataset_dir, selected_file_keys):
    """(existing absolute paths for the selected file keys, keys of selected
    files that are absent). Missing files are reported, not raised: optional
    inputs (imputed/filtered) legitimately don't exist for every run."""
    choices = {f.key: f for f in PRESETS['saint_inputs'].files}
    paths, missing = [], []
    for key in selected_file_keys:
        path = os.path.abspath(os.path.join(dataset_dir, choices[key].filename))
        if os.path.exists(path):
            paths.append(path)
        else:
            missing.append(key)
    return paths, missing


def zip_files(file_paths, zip_path):
    """Zip the given files under their basenames. Raises ValueError on an
    empty list."""
    if not file_paths:
        raise ValueError("No files to zip.")
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for path in file_paths:
            zf.write(path, os.path.basename(path))
    logger.info("Wrote %d files to %s", len(file_paths), zip_path)
    return zip_path


def write_gene_list(genes_df, path):
    """Write a gene-list frame as plain text: one gene per line, or
    Bait<TAB>Gene lines when a Bait column is present."""
    with open(path, 'w', newline='') as fh:
        if 'Bait' in genes_df.columns:
            for _, row in genes_df.iterrows():
                fh.write(f"{row['Bait']}\t{row['Gene']}\n")
        else:
            for gene in genes_df['Gene']:
                fh.write(f"{gene}\n")
    return path
