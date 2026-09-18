"""Build the node and edge tables a thresholded ProxiMate network is drawn from.

Pure functions over annotated_scores: no Cytoscape, no I/O beyond the cached BioGRID
reader.  Every visual encoding is computed here into a column and mapped by a
passthrough in Cytoscape, so changing a look later is a column update.
"""

import numpy as np
import pandas as pd

from QC_plots import _load_biogrid_cached, apply_score_thresholds
from download_presets import _CYTOSCAPE_ATTRS

THRESHOLD_COLUMNS = ('SaintScore', 'BFDR', 'WD', 'WDFDR')
LABEL_POLICIES = ('all', 'baits', 'none')

# Node encodings by role.
NODE_SHAPE = {'bait': 'DIAMOND', 'prey': 'ELLIPSE'}
NODE_SIZE = {'bait': 60.0, 'prey': 35.0}
NODE_FILL = {'bait': '#E69F00', 'prey': '#B8C4D6'}
NODE_BORDER = {'bait': 3.0, 'prey': 1.0}
NODE_BORDER_COLOR = '#2F3B4A'
LABEL_SIZE = {'bait': 16, 'prey': 11}
NODE_ALPHA = {'visible': 255, 'faded': 60}

# Edge encodings by kind.  Widths are z-scores of log10 abundance clipped to ±Z_CLIP
# and mapped onto a band, so orders of magnitude read as a gradation, not a spike.
EDGE_COLOR = {'proximity': '#0072B2', 'literature': '#7D838C'}
EDGE_ALPHA = {'proximity': 220, 'literature': 120}
PROXIMITY_WIDTH = (1.0, 8.0)
LITERATURE_WIDTH = 1.5
Z_CLIP = 2.5


def zwidth(values, band, clip=Z_CLIP):
    """Map values through log10 then z-score onto ``band = (lo, hi)``.

    NaNs and non-positive values map to the band floor; a constant column maps to the
    band midpoint.
    """
    lo, hi = band
    x = np.log10(values.astype(float).where(values > 0))
    sigma = float(np.nanstd(x))
    if np.isnan(sigma) or sigma == 0:
        z = pd.Series(np.zeros(len(x)), index=values.index)
    else:
        z = ((x - float(np.nanmean(x))) / sigma).clip(-clip, clip)
    return (lo + (z + clip) / (2 * clip) * (hi - lo)).fillna(lo)


def _require(df, columns):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"annotated scores lack the column(s) {', '.join(missing)}")


def _abundance_column(df):
    """The per-interaction abundance SAINT wrote, by input type."""
    for column in ('AvgIntensity', 'AvgSpec'):
        if column in df.columns:
            return column
    return None


def bait_id_column(df):
    return 'Bait_Accession' if 'Bait_Accession' in df.columns else 'Bait.ID'


def prey_symbol_column(df):
    return 'First_Prey_Gene' if 'First_Prey_Gene' in df.columns else 'PreyGene'


def prey_prey_pairs(prey_ids, bait_ids, biogrid_path):
    """BioGRID pairs with both interactors among the preys and neither a bait.

    The same selection ``QC_plots.calculate_network_degrees`` counts degrees over.
    Returns an empty frame when the summary cannot be read.
    """
    biogrid = _load_biogrid_cached(biogrid_path) if biogrid_path else None
    if biogrid is None:
        return pd.DataFrame(columns=['source', 'target'])
    a = biogrid['SWISS-PROT Accessions Interactor A']
    b = biogrid['SWISS-PROT Accessions Interactor B']
    keep = a.isin(prey_ids) & b.isin(prey_ids) & ~a.isin(bait_ids) & ~b.isin(bait_ids) & (a != b)
    pairs = pd.DataFrame({'source': a[keep].to_numpy(), 'target': b[keep].to_numpy()})
    # BioGRID keys pairs in one orientation; the summary may still carry both.
    ordered = pd.DataFrame(np.sort(pairs.to_numpy(), axis=1), columns=['source', 'target'])
    return ordered.drop_duplicates().reset_index(drop=True)


def display_labels(nodes, policy):
    if policy not in LABEL_POLICIES:
        raise ValueError(f"unknown label policy {policy!r}; expected one of {', '.join(LABEL_POLICIES)}")
    if policy == 'all':
        return nodes['symbol']
    if policy == 'none':
        return pd.Series([''] * len(nodes), index=nodes.index)
    return nodes['symbol'].where(nodes['role'] == 'bait', '')


def build(df, thresholds, baits=None, prey_prey=True, biogrid_path=None, label_policy='all'):
    """(nodes, edges) for the interactions passing ``thresholds``.

    Nodes are keyed on UniProt accession (``id``); a prey that is also a bait keeps
    one node with the bait role.  Bait-prey edges carry every score column present
    and their thresholds; prey-prey edges are BioGRID pairs among the passing preys.
    """
    _require(df, ('Experiment.ID', 'First_ID', *THRESHOLD_COLUMNS))
    bait_col = bait_id_column(df)
    symbol_col = prey_symbol_column(df)
    _require(df, (bait_col, symbol_col))

    passing = apply_score_thresholds(df, thresholds)
    if baits:
        passing = passing[passing['Experiment.ID'].isin(baits)]
    passing = passing[passing['First_ID'].notna() & passing[bait_col].notna()]
    if passing.empty:
        raise ValueError("no interactions pass the thresholds for the selected bait(s)")

    baits_tbl = (passing[['Experiment.ID', bait_col]].drop_duplicates()
                 .rename(columns={bait_col: 'id', 'Experiment.ID': 'symbol'}))
    baits_tbl['id'] = baits_tbl['id'].astype(str)
    baits_tbl['role'] = 'bait'
    bait_ids = set(baits_tbl['id'])

    prey_cols = {'First_ID': 'id', symbol_col: 'symbol'}
    for extra in ('first_SCL', 'Main location', 'Human_Complex', 'GO_CC', 'In.BioGRID'):
        if extra in passing.columns:
            prey_cols[extra] = extra
    preys_tbl = passing[list(prey_cols)].rename(columns=prey_cols).drop_duplicates('id')
    preys_tbl['id'] = preys_tbl['id'].astype(str)
    preys_tbl = preys_tbl[~preys_tbl['id'].isin(bait_ids)]
    preys_tbl['role'] = 'prey'

    nodes = pd.concat([baits_tbl, preys_tbl], ignore_index=True)
    nodes['symbol'] = nodes['symbol'].fillna(nodes['id']).astype(str)
    nodes['shape'] = nodes['role'].map(NODE_SHAPE)
    nodes['size'] = nodes['role'].map(NODE_SIZE)
    nodes['fill'] = nodes['role'].map(NODE_FILL)
    nodes['border'] = nodes['role'].map(NODE_BORDER)
    nodes['border_color'] = NODE_BORDER_COLOR
    nodes['label_size'] = nodes['role'].map(LABEL_SIZE)
    nodes['display_label'] = display_labels(nodes, label_policy)
    nodes['node_alpha'] = NODE_ALPHA['visible']
    if 'In.BioGRID' in nodes.columns:
        nodes['In.BioGRID'] = nodes['In.BioGRID'].eq(True)

    edges = pd.DataFrame({'source': passing[bait_col].astype(str).to_numpy(),
                          'target': passing['First_ID'].astype(str).to_numpy()})
    edges['interaction'] = 'proximity'
    for col in _CYTOSCAPE_ATTRS:
        if col in passing.columns:
            edges[col] = passing[col].to_numpy()
    abundance = _abundance_column(passing)
    edges['width'] = (zwidth(passing[abundance], PROXIMITY_WIDTH).to_numpy()
                      if abundance else np.mean(PROXIMITY_WIDTH))
    edges['color'] = EDGE_COLOR['proximity']
    edges['alpha'] = EDGE_ALPHA['proximity']

    if prey_prey:
        pairs = prey_prey_pairs(set(preys_tbl['id']), bait_ids, biogrid_path)
        if len(pairs):
            pairs['interaction'] = 'literature'
            pairs['width'] = LITERATURE_WIDTH
            pairs['color'] = EDGE_COLOR['literature']
            pairs['alpha'] = EDGE_ALPHA['literature']
            edges = pd.concat([edges, pairs], ignore_index=True)

    edges['visible'] = True
    edges['name'] = [f'{s} ({k}) {t}' for s, k, t in
                     zip(edges['source'], edges['interaction'], edges['target'])]
    return nodes.reset_index(drop=True), edges


def is_tighter_or_equal(new, drawn):
    """Whether ``new`` thresholds keep a subset of what ``drawn`` kept."""
    return (new['SaintScore'] >= drawn['SaintScore'] and new['BFDR'] <= drawn['BFDR']
            and new['WD'] >= drawn['WD'] and new['WDFDR'] <= drawn['WDFDR'])


def edge_visibility(edges, thresholds):
    """The ``visible`` column under new thresholds, for edges already drawn.

    A proximity edge is visible when its scores pass; a literature edge when both
    of its preys keep at least one visible proximity edge.
    """
    proximity = edges['interaction'] == 'proximity'
    passing_index = apply_score_thresholds(edges[proximity], thresholds).index
    visible = pd.Series(False, index=edges.index)
    visible[passing_index] = True
    lit = edges[~proximity]
    if len(lit):
        live = set(edges.loc[visible & proximity, 'target'])
        visible[lit.index] = lit['source'].isin(live) & lit['target'].isin(live)
    return visible


def node_alpha(nodes, edges):
    """Nodes touching no visible edge fade; the others stay opaque."""
    live_edges = edges[edges['visible'].astype(bool)]
    live = set(live_edges['source']) | set(live_edges['target'])
    return pd.Series(np.where(nodes['id'].isin(live), NODE_ALPHA['visible'], NODE_ALPHA['faded']),
                     index=nodes.index)
