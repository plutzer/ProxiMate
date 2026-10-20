"""Build the node and edge tables a thresholded ProxiMate network is drawn from.

Pure functions over annotated_scores: no Cytoscape, no I/O beyond the cached BioGRID
and CORUM readers.  Every visual encoding is computed here into a column and mapped by
a passthrough in Cytoscape, so changing a look later is a column update.  The helpers
after ``node_alpha`` answer questions about the drawn tables (degrees, loners,
relations, communities) and never touch Cytoscape either.
"""

import itertools
import random
import re

import numpy as np
import pandas as pd

from QC_plots import _load_biogrid_cached, apply_score_thresholds
from cytoscape_p4c import edge_name
from download_presets import CYTOSCAPE_ATTRS

THRESHOLD_COLUMNS = ('SaintScore', 'BFDR', 'WD', 'WDFDR')
LABEL_POLICIES = ('all', 'baits', 'none')
EDGE_WIDTH_SOURCES = ('abundance', 'SaintScore', 'WD', 'FoldChange', 'uniform')
BIOGRID_SCOPES = ('all', 'multivalidated')
RELATIONS = ('interactors', 'singletons', 'partners', 'cocomplex')
# A send that would draw more visible edges than this asks for confirmation first.
SEND_WARN_EDGES = 5000
# Decimals the selection table shows per score column.
SELECTION_DECIMALS = {'SaintScore': 3, 'BFDR': 3, 'WD': 3, 'WDFDR': 3, 'FoldChange': 2}

# A style preset is every look the tables encode: node encodings by role (sizes are
# (width, height)), edge encodings by kind, and the send options that go with it in the
# GUI.  ``known_border`` is the border width of a prey whose interaction with a drawn
# bait is BioGRID-multivalidated; None leaves known preys unmarked.  Edge widths are
# z-scores of a log quantity clipped to ±Z_CLIP and mapped onto the preset's band, so
# orders of magnitude read as a gradation, not a spike.
STYLE_PRESETS = {
    'abundance_weighted': {
        'label': 'Abundance/literature weighted',
        'description': 'Orange diamond baits and blue-grey preys; bait-prey width follows '
                       'abundance and BioGRID edges thicken with publications.',
        'node_shape': {'bait': 'DIAMOND', 'prey': 'ELLIPSE'},
        'node_size': {'bait': (60.0, 60.0), 'prey': (35.0, 35.0)},
        'node_fill': {'bait': '#E69F00', 'prey': '#B8C4D6'},
        'node_border': {'bait': 3.0, 'prey': 1.0},
        'node_border_color': '#2F3B4A',
        'label_size': {'bait': 16, 'prey': 11},
        'label_color': {'bait': '#1F2933', 'prey': '#1F2933'},
        'known_border': None,
        'edge_color': {'proximity': '#0072B2', 'literature': '#7D838C', 'complex': '#000000'},
        'edge_alpha': {'proximity': 220, 'literature': 120, 'complex': 200},
        'proximity_width': (1.0, 8.0), 'literature_width': (1.0, 4.0), 'complex_width': 2.0,
        'send_options': {'width_source': 'abundance', 'prey_prey': True, 'biogrid_scope': 'all',
                         'literature_weighted': True, 'corum': False, 'layout': 'force-directed'},
    },
    'highlight_knowns': {
        'label': 'Highlight knowns',
        'description': 'Black ellipse baits and grey rounded preys with a red border on preys '
                       'whose bait interaction BioGRID marks multivalidated; grey edges, black '
                       'complex edges, one edge width.',
        'node_shape': {'bait': 'ELLIPSE', 'prey': 'ROUND_RECTANGLE'},
        'node_size': {'bait': (100.0, 60.0), 'prey': (60.0, 20.0)},
        'node_fill': {'bait': '#000000', 'prey': '#A6A6A6'},
        'node_border': {'bait': 0.0, 'prey': 0.0},
        'node_border_color': '#B30000',
        'label_size': {'bait': 12, 'prey': 8},
        'label_color': {'bait': '#FFFFFF', 'prey': '#000000'},
        'known_border': 4.0,
        'edge_color': {'proximity': '#888888', 'literature': '#CCCCCC', 'complex': '#000000'},
        'edge_alpha': {'proximity': 255, 'literature': 255, 'complex': 255},
        # Under the preset's own send options (uniform width, unweighted literature) every
        # edge is 2.5: the band midpoint, the literature floor and the complex width.
        'proximity_width': (1.0, 4.0), 'literature_width': (2.5, 4.0), 'complex_width': 2.5,
        'send_options': {'width_source': 'uniform', 'prey_prey': True, 'biogrid_scope': 'multivalidated',
                         'literature_weighted': False, 'corum': True, 'layout': 'kamada-kawai'},
    },
}
# The table columns a preset writes; each is mapped by a passthrough in Cytoscape.
NODE_LOOK = ['shape', 'node_width', 'node_height', 'fill', 'border', 'border_color',
             'label_size', 'label_color']
EDGE_LOOK = ['color', 'alpha', 'width']

NODE_ALPHA = {'visible': 255, 'faded': 60}
# Fills for clustered nodes, one per community, cycling past twelve.
COMMUNITY_FILL = ('#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F', '#EDC948',
                  '#B07AA1', '#FF9DA7', '#9C755F', '#BAB0AC', '#86BCB6', '#D37295')
Z_CLIP = 2.5

# Clustering: proximity edges pull by banded abundance; a selection smaller than
# MIN_CLUSTER_NODES carries no structure to partition.
CLUSTER_WEIGHT = (0.5, 1.5)
MIN_CLUSTER_NODES = 4
# Centre-to-centre spacing of packed nodes, above the largest node size.
NODE_SPACING = 70.0
# How clustered communities are arranged: packed circles, or an igraph layout of each
# community inside its own disc with the discs placed by the edges between communities.
CLUSTER_LAYOUTS = {'circles': None, 'force-directed': 'layout_fruchterman_reingold',
                   'kamada-kawai': 'layout_kamada_kawai'}

ISOFORM = re.compile(r'-\d+$')


def zwidth(values, band, clip=Z_CLIP, log_fn=np.log10):
    """Map values through ``log_fn`` then z-score onto ``band = (lo, hi)``.

    NaNs and non-positive values map to the band floor; a constant column maps to the
    band midpoint.
    """
    lo, hi = band
    x = log_fn(values.astype(float).where(values > 0))
    sigma = float(np.nanstd(x)) if x.notna().sum() > 1 else 0.0
    if sigma == 0:
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


# --- edge widths ------------------------------------------------------------------------

def proximity_widths(edges, source, band):
    """Width per bait-prey edge from the quantity ``source`` names, over ``band``.

    ``abundance`` bands log10 AvgIntensity/AvgSpec (one width when neither column
    exists); ``WD`` and ``FoldChange`` band their own log10; ``SaintScore`` is linear
    over 0-1; ``uniform`` is the band midpoint.
    """
    if source not in EDGE_WIDTH_SOURCES:
        raise ValueError(f"unknown edge width source {source!r}; "
                         f"expected one of {', '.join(EDGE_WIDTH_SOURCES)}")
    lo, hi = band
    if source == 'abundance':
        column = _abundance_column(edges)
        if column is None:
            source = 'uniform'
    if source == 'uniform':
        return pd.Series(np.mean(band), index=edges.index)
    if source != 'abundance':
        column = source
        _require(edges, (column,))
    if source == 'SaintScore':
        return lo + edges[column].astype(float).fillna(0).clip(0, 1) * (hi - lo)
    return zwidth(edges[column], band)


def literature_widths(edges, weighted, band):
    """Width per BioGRID edge: the band floor, or banded log1p publication count."""
    if not weighted:
        return pd.Series(band[0], index=edges.index, dtype=float)
    return zwidth(edges['n_publications'], band, log_fn=np.log1p)


def edge_widths(edges, preset, width_source='abundance', literature_weighted=False):
    """The ``width`` column for a mixed edge table under a preset's bands; complex
    edges are one width."""
    width = pd.Series(preset['complex_width'], index=edges.index, dtype=float)
    proximity = edges['interaction'] == 'proximity'
    width[proximity] = proximity_widths(edges[proximity], width_source, preset['proximity_width'])
    literature = edges['interaction'] == 'literature'
    if literature.any():
        width[literature] = literature_widths(edges[literature], literature_weighted,
                                              preset['literature_width'])
    return width


def style_tables(nodes, edges, preset, width_source='abundance', literature_weighted=False):
    """Write the ``NODE_LOOK`` and ``EDGE_LOOK`` columns of ``preset`` onto drawn tables,
    in place.  A preset with ``known_border`` needs the ``Multivalidated`` edge column and
    gives that border to every prey with a multivalidated proximity edge."""
    if preset not in STYLE_PRESETS:
        raise ValueError(f"unknown style preset {preset!r}; expected one of {', '.join(STYLE_PRESETS)}")
    look = STYLE_PRESETS[preset]
    role = nodes['role']
    nodes['shape'] = role.map(look['node_shape'])
    nodes['node_width'] = role.map(lambda r: look['node_size'][r][0])
    nodes['node_height'] = role.map(lambda r: look['node_size'][r][1])
    nodes['fill'] = role.map(look['node_fill'])
    nodes['border'] = role.map(look['node_border'])
    nodes['border_color'] = look['node_border_color']
    nodes['label_size'] = role.map(look['label_size'])
    nodes['label_color'] = role.map(look['label_color'])
    if look['known_border'] is not None:
        _require(edges, ('Multivalidated',))
        proximity = edges[edges['interaction'] == 'proximity']
        known = set(proximity.loc[proximity['Multivalidated'].eq(True), 'target'])
        nodes.loc[(role == 'prey') & nodes['id'].isin(known), 'border'] = look['known_border']
    kind = edges['interaction']
    edges['color'] = kind.map(look['edge_color'])
    edges['alpha'] = kind.map(look['edge_alpha'])
    edges['width'] = edge_widths(edges, look, width_source, literature_weighted)
    return nodes, edges


# --- reference layers --------------------------------------------------------------------

def _publication_counts(sources):
    """Distinct publications per summary row, from the semicolon-joined column."""
    return sources.fillna('').astype(str).map(
        lambda s: len({p.strip() for p in s.split(';') if p.strip()}))


def prey_prey_pairs(prey_ids, bait_ids, biogrid_path):
    """BioGRID pairs with both interactors among the preys and neither a bait.

    The same selection ``QC_plots.calculate_network_degrees`` counts degrees over.
    Each pair carries ``n_publications`` and ``multivalidated``; a pair the summary
    lists in both orientations keeps the larger count.  Returns an empty frame when
    the summary cannot be read.
    """
    columns = ['source', 'target', 'n_publications', 'multivalidated']
    biogrid = _load_biogrid_cached(biogrid_path) if biogrid_path else None
    if biogrid is None:
        return pd.DataFrame(columns=columns)
    a = biogrid['SWISS-PROT Accessions Interactor A']
    b = biogrid['SWISS-PROT Accessions Interactor B']
    keep = a.isin(prey_ids) & b.isin(prey_ids) & ~a.isin(bait_ids) & ~b.isin(bait_ids) & (a != b)
    if not keep.any():
        return pd.DataFrame(columns=columns)
    ends = np.sort(np.column_stack([a[keep].to_numpy(), b[keep].to_numpy()]), axis=1)
    pairs = pd.DataFrame({
        'source': ends[:, 0], 'target': ends[:, 1],
        'n_publications': _publication_counts(biogrid.loc[keep, 'Publication Source']).to_numpy(),
        'multivalidated': biogrid.loc[keep, 'Multivalidated'].fillna(False).astype(bool).to_numpy()})
    return (pairs.groupby(['source', 'target'], as_index=False)
                 .agg(n_publications=('n_publications', 'max'), multivalidated=('multivalidated', 'any')))


_corum_cache = {}


def _load_corum_cached(path):
    """One row per CORUM complex: ``complex_name`` and its ``subunits`` accessions.

    Isoform suffixes (``P05067-4``) are dropped so subunits match the node
    ``accession`` column.  The file is read once per path.
    """
    if path not in _corum_cache:
        raw = pd.read_table(path, encoding='latin-1', usecols=['complex_name', 'subunits_uniprot_id'])
        subunits = [sorted({ISOFORM.sub('', s.strip()) for s in str(cell).split(';') if s.strip()})
                    for cell in raw['subunits_uniprot_id']]
        _corum_cache[path] = pd.DataFrame({'complex_name': raw['complex_name'], 'subunits': subunits})
    return _corum_cache[path]


def complex_pairs(node_ids, bait_preys, corum_path, min_members=3, min_fraction=0.5):
    """CORUM pairs among the drawn proteins, from the complexes the screen recovered.

    A complex qualifies when at least ``min_members`` of its subunits are drawn and
    one bait, counted with its preys, covers at least ``min_fraction`` of the full
    membership.  Every drawn pair of a qualifying complex is returned, including
    subunits recovered only under other baits.  Returns ``(pairs, membership)``:
    pairs with ``source < target``, ``complex_names`` and ``n_complexes``; membership
    maps each drawn accession to the semicolon-joined names of its qualifying complexes.
    ``node_ids`` and ``bait_preys`` are accessions.
    """
    ids = set(node_ids)
    complexes = _load_corum_cached(corum_path)
    pairs, membership = {}, {}
    for name, subunits in zip(complexes['complex_name'], complexes['subunits']):
        members = set(subunits)
        drawn = sorted(members & ids)
        if len(drawn) < min_members:
            continue
        share = max((len(members & preys) / len(members) for preys in bait_preys.values()), default=0.0)
        if share < min_fraction:
            continue
        for u in drawn:
            membership.setdefault(u, []).append(name)
        for a, b in itertools.combinations(drawn, 2):
            pairs.setdefault((a, b), []).append(name)
    frame = pd.DataFrame([{'source': a, 'target': b, 'complex_names': ';'.join(names),
                           'n_complexes': len(names)} for (a, b), names in sorted(pairs.items())],
                         columns=['source', 'target', 'complex_names', 'n_complexes'])
    return frame, {u: ';'.join(names) for u, names in membership.items()}


def display_labels(nodes, policy):
    if policy not in LABEL_POLICIES:
        raise ValueError(f"unknown label policy {policy!r}; expected one of {', '.join(LABEL_POLICIES)}")
    if policy == 'all':
        return nodes['symbol']
    if policy == 'none':
        return pd.Series([''] * len(nodes), index=nodes.index)
    return nodes['symbol'].where(nodes['role'] == 'bait', '')


# --- build -----------------------------------------------------------------------------------

def build(df, thresholds, baits=None, prey_prey=True, biogrid_path=None, label_policy='all',
          width_source='abundance', literature_weighted=False, corum_path=None,
          corum_min_members=3, corum_min_fraction=0.5, preset='abundance_weighted'):
    """(nodes, edges) for the interactions passing ``thresholds``, styled by ``preset``.

    Prey nodes are keyed on UniProt accession and bait nodes on bait name (``id``), so
    two constructs of one protein stay separate baits; ``accession`` holds the protein
    of every node.  A prey that is also a bait keeps only its bait node(s) and a bait
    never draws its own protein.  Bait-prey edges carry every score column present
    and their thresholds; literature edges are BioGRID pairs among the passing preys;
    complex edges are CORUM pairs among the drawn proteins when ``corum_path`` is given,
    and the node table then carries each protein's qualifying ``complexes``.
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
                 .rename(columns={'Experiment.ID': 'id', bait_col: 'accession'}).astype(str))
    baits_tbl['symbol'] = baits_tbl['id']
    baits_tbl['role'] = 'bait'
    bait_accessions = set(baits_tbl['accession'])
    # accession -> the node ids carrying it: a prey is its own accession, a bait
    # accession names every bait construct of that protein
    bait_nodes = baits_tbl.groupby('accession')['id'].agg(list)

    prey_cols = {'First_ID': 'id', symbol_col: 'symbol'}
    for extra in ('first_SCL', 'Main location', 'Human_Complex', 'GO_CC', 'In.BioGRID'):
        if extra in passing.columns:
            prey_cols[extra] = extra
    preys_tbl = passing[list(prey_cols)].rename(columns=prey_cols).drop_duplicates('id')
    preys_tbl['id'] = preys_tbl['id'].astype(str)
    preys_tbl = preys_tbl[~preys_tbl['id'].isin(bait_accessions)]
    preys_tbl['accession'] = preys_tbl['id']
    preys_tbl['role'] = 'prey'

    nodes = pd.concat([baits_tbl, preys_tbl], ignore_index=True)
    nodes['symbol'] = nodes['symbol'].fillna(nodes['id']).astype(str)
    nodes['display_label'] = display_labels(nodes, label_policy)
    nodes['node_alpha'] = NODE_ALPHA['visible']
    if 'In.BioGRID' in nodes.columns:
        nodes['In.BioGRID'] = nodes['In.BioGRID'].eq(True)

    # A bait detecting its own accession is a self-interaction and draws nothing; a bait
    # detecting another bait's protein connects to every construct of that bait.
    edges = pd.DataFrame({'source': passing['Experiment.ID'].astype(str).to_numpy(),
                          'target': passing['First_ID'].astype(str).to_numpy(),
                          'accession': passing[bait_col].astype(str).to_numpy()})
    edges['interaction'] = 'proximity'
    for col in (*CYTOSCAPE_ATTRS, 'AvgSpec', 'Multivalidated'):
        if col in passing.columns:
            edges[col] = passing[col].to_numpy()
    edges = edges[edges['target'] != edges['accession']].drop(columns='accession')
    edges['target'] = edges['target'].map(lambda t: bait_nodes.get(t, [t]))
    edges = edges.explode('target', ignore_index=True)

    if prey_prey:
        pairs = prey_prey_pairs(set(preys_tbl['id']), bait_accessions, biogrid_path)
        if len(pairs):
            pairs['interaction'] = 'literature'
            edges = pd.concat([edges, pairs], ignore_index=True)

    if corum_path:
        preys_of = {str(b): set(t.astype(str)) | {str(b)}
                    for b, t in passing.groupby(bait_col)['First_ID']}
        pairs, membership = complex_pairs(set(nodes['accession']), preys_of, corum_path,
                                          corum_min_members, corum_min_fraction)
        for end in ('source', 'target'):
            pairs[end] = pairs[end].map(lambda a: bait_nodes.get(a, [a]))
            pairs = pairs.explode(end, ignore_index=True)
        proximity = edges[edges['interaction'] == 'proximity']
        drawn = {frozenset(p) for p in zip(proximity['source'], proximity['target'])}
        pairs = pairs[[frozenset((s, t)) not in drawn for s, t in zip(pairs['source'], pairs['target'])]]
        nodes['complexes'] = nodes['accession'].map(membership).fillna('')
        if len(pairs):
            pairs['interaction'] = 'complex'
            edges = pd.concat([edges, pairs], ignore_index=True)

    edges['visible'] = True
    edges['name'] = [edge_name(s, k, t) for s, k, t in
                     zip(edges['source'], edges['interaction'], edges['target'])]
    return style_tables(nodes.reset_index(drop=True), edges, preset, width_source, literature_weighted)


def is_tighter_or_equal(new, drawn):
    """Whether ``new`` thresholds keep a subset of what ``drawn`` kept."""
    return (new['SaintScore'] >= drawn['SaintScore'] and new['BFDR'] <= drawn['BFDR']
            and new['WD'] >= drawn['WD'] and new['WDFDR'] <= drawn['WDFDR'])


def edge_visibility(edges, thresholds, biogrid_scope='all'):
    """The ``visible`` column under new thresholds, for edges already drawn.

    A proximity edge is visible when its scores pass; a literature or complex edge
    when both of its ends keep at least one visible proximity edge.  Under the
    ``multivalidated`` scope a literature edge also needs BioGRID's multivalidation.
    """
    if biogrid_scope not in BIOGRID_SCOPES:
        raise ValueError(f"unknown BioGRID scope {biogrid_scope!r}; "
                         f"expected one of {', '.join(BIOGRID_SCOPES)}")
    proximity = edges['interaction'] == 'proximity'
    passing_index = apply_score_thresholds(edges[proximity], thresholds).index
    visible = pd.Series(False, index=edges.index)
    visible[passing_index] = True
    live_proximity = edges[visible & proximity]
    live = set(live_proximity['source']) | set(live_proximity['target'])
    visible = visible.where(proximity, edges['source'].isin(live) & edges['target'].isin(live))
    literature = edges['interaction'] == 'literature'
    if biogrid_scope == 'multivalidated' and literature.any():
        visible &= ~literature | edges['multivalidated'].eq(True)
    return visible


def network_counts(nodes, edges, thresholds, biogrid_scope='all'):
    """What drawing ``build``'s tables shows: bait, prey and node counts, the visible
    edges of each kind under ``biogrid_scope`` and their total (``edges``), and the
    BioGRID edges the scope sends hidden (``hidden_edges``)."""
    visible = edge_visibility(edges, thresholds, biogrid_scope)
    kinds = edges.loc[visible, 'interaction'].value_counts()
    roles = nodes['role'].value_counts()
    return {'baits': int(roles.get('bait', 0)), 'preys': int(roles.get('prey', 0)),
            'nodes': len(nodes), 'proximity_edges': int(kinds.get('proximity', 0)),
            'literature_edges': int(kinds.get('literature', 0)),
            'complex_edges': int(kinds.get('complex', 0)),
            'edges': int(visible.sum()), 'hidden_edges': int((~visible).sum())}


def selection_table(detail, nodes):
    """The edges around a selection as the Cytoscape tab lists them: a ``Target gene``
    column after ``target`` and the scores rounded by ``SELECTION_DECIMALS``."""
    table = detail.round(SELECTION_DECIMALS)
    symbols = dict(zip(nodes['id'], nodes['symbol']))
    table.insert(table.columns.get_loc('target') + 1, 'Target gene', table['target'].map(symbols))
    return table


def node_alpha(nodes, edges):
    """Nodes touching no visible edge fade; the others stay opaque."""
    live_edges = edges[edges['visible'].astype(bool)]
    live = set(live_edges['source']) | set(live_edges['target'])
    return pd.Series(np.where(nodes['id'].isin(live), NODE_ALPHA['visible'], NODE_ALPHA['faded']),
                     index=nodes.index)


# --- questions about the drawn tables ------------------------------------------------------
# Hidden edges do not count anywhere below: a hidden edge failed the thresholds or was
# hidden on purpose, and either way it is not part of the network on screen.

def _visible(edges):
    return edges[edges['visible'].astype(bool)]


def degrees(edges):
    """Visible-edge degree per node id, any kind of edge."""
    live = _visible(edges)
    return pd.concat([live['source'], live['target']]).value_counts()


def _baits_of(edges, prey_ids):
    """prey id -> the set of baits with a visible proximity edge to it."""
    live = _visible(edges)
    proximity = live[(live['interaction'] == 'proximity') & live['target'].isin(prey_ids)]
    return proximity.groupby('target')['source'].agg(set)


def loners(nodes, edges, bait):
    """Preys whose only visible neighbour, of any kind, is ``bait``.

    A subunit of a drawn complex or a prey with a literature partner is never a
    loner, so a loner set can be moved without pulling anything else along.
    """
    live = _visible(edges)
    touching = live[(live['source'] == bait) | (live['target'] == bait)]
    prey = set(nodes.loc[nodes['role'] == 'prey', 'id'])
    mine = (set(touching['source']) | set(touching['target'])) & prey
    degree = degrees(edges)
    return sorted(p for p in mine if degree.get(p, 0) == 1)


def satellites(nodes, edges, bait, positions):
    """The preys that belong beside ``bait``, over proximity edges alone.

    Returns ``(only, nearer)``: the preys whose only bait it is, and the two-bait preys
    that lie closer to it than to their other bait by ``positions`` ({id: (x, y)}).
    """
    prey = set(nodes.loc[nodes['role'] == 'prey', 'id'])
    baits_of = _baits_of(edges, prey)
    mine = [p for p, b in baits_of.items() if bait in b]
    only = sorted(p for p in mine if baits_of[p] == {bait})

    def d2(a, b):
        (x1, y1), (x2, y2) = positions[a], positions[b]
        return (x1 - x2) ** 2 + (y1 - y2) ** 2

    nearer = sorted(p for p in mine if len(baits_of[p]) == 2
                    and d2(p, bait) < d2(p, next(iter(baits_of[p] - {bait}))))
    return only, nearer


def resolve_node(nodes, query):
    """The drawn node ``query`` names, by accession or symbol, case-insensitively."""
    q = str(query).strip().lower()
    hit = nodes[(nodes['id'].str.lower() == q) | (nodes['symbol'].str.lower() == q)]
    if len(hit) == 0:
        raise ValueError(f"{query!r} is not in the drawn network")
    if len(hit) > 1:
        raise ValueError(f"{query!r} names {len(hit)} nodes: {', '.join(hit['id'])}")
    return hit.iloc[0]


def related(nodes, edges, seed, relation, min_saint=None, max_bfdr=None,
            min_abundance=None, min_publications=None):
    """The ids ``relation`` names for the node ``seed``.

    interactors  preys of a bait passing the SAINT / BFDR / abundance cuts
    singletons   preys whose only bait is the seed
    partners     literature or complex neighbours of any node; ``min_publications``
                 applies to the literature ones
    cocomplex    other subunits of the drawn complexes the seed belongs to
    """
    if relation not in RELATIONS:
        raise ValueError(f"unknown relation {relation!r}; expected one of {', '.join(RELATIONS)}")
    node = resolve_node(nodes, seed)
    sid = node['id']
    live = _visible(edges)
    if relation in ('interactors', 'singletons'):
        if node['role'] != 'bait':
            raise ValueError(f"{node['symbol']} is a prey; {relation} needs a bait")
        proximity = live[live['interaction'] == 'proximity']
        mine = proximity[proximity['source'] == sid]
        if relation == 'interactors':
            keep = pd.Series(True, index=mine.index)
            if min_saint is not None:
                keep &= mine['SaintScore'] >= min_saint
            if max_bfdr is not None:
                keep &= mine['BFDR'] <= max_bfdr
            if min_abundance is not None:
                keep &= mine[_abundance_column(mine)] >= min_abundance
            mine = mine[keep]
        else:
            baits_of = proximity.groupby('target')['source'].agg(set)
            only = {t for t, b in baits_of.items() if b == {sid}}
            mine = mine[mine['target'].isin(only)]
        return sorted(set(mine['target']))
    if relation == 'partners':
        reference = live[live['interaction'].isin(('literature', 'complex'))]
        if min_publications is not None:
            reference = reference[(reference['interaction'] != 'literature')
                                  | (reference['n_publications'] >= min_publications)]
        touching = reference[(reference['source'] == sid) | (reference['target'] == sid)]
        return sorted((set(touching['source']) | set(touching['target'])) - {sid})
    if 'complexes' not in nodes.columns:
        raise ValueError("the network was drawn without CORUM complex edges")
    names = nodes['complexes'].fillna('').astype(str)
    mine = set(str(node['complexes']).split(';')) - {''}
    return sorted(i for i, c in zip(nodes['id'], names)
                  if i != sid and c and set(c.split(';')) & mine)


def _weighted_graph(edges, ids, literature_weight):
    """An igraph over the visible edges among ``ids`` with edge attribute ``weight``:
    proximity edges weigh their banded log abundance, literature edges
    ``literature_weight`` times log1p publications, complex edges ``literature_weight``.
    A weight of zero drops the edge."""
    import igraph

    row = {i: k for k, i in enumerate(ids)}
    live = _visible(edges)
    among = live[live['source'].isin(row) & live['target'].isin(row)]
    weight = pd.Series(float(literature_weight), index=among.index)
    proximity = among['interaction'] == 'proximity'
    abundance = _abundance_column(among)
    weight[proximity] = (zwidth(among.loc[proximity, abundance], CLUSTER_WEIGHT) if abundance
                         else np.mean(CLUSTER_WEIGHT))
    literature = among['interaction'] == 'literature'
    weight[literature] = literature_weight * np.log1p(among.loc[literature, 'n_publications'].astype(float))
    among, weight = among[weight > 0], weight[weight > 0]
    graph = igraph.Graph(n=len(ids), edges=[(row[s], row[t]) for s, t in zip(among['source'], among['target'])])
    graph.es['weight'] = weight.tolist()
    return graph


def cluster(edges, ids, resolution=1.0, seed=17, literature_weight=1.0):
    """Leiden communities (modularity) over the visible edges among ``ids``, weighted
    as in ``_weighted_graph``.  Returns a Series of community numbers indexed by id,
    numbered by size, largest first."""
    import igraph

    ids = list(dict.fromkeys(str(i) for i in ids))
    if len(ids) < MIN_CLUSTER_NODES:
        raise ValueError(f"a selection needs at least {MIN_CLUSTER_NODES} nodes to cluster; got {len(ids)}")
    graph = _weighted_graph(edges, ids, literature_weight)
    igraph.set_random_number_generator(random.Random(seed))
    parts = graph.community_leiden(objective_function='modularity', resolution=resolution,
                                   weights='weight' if graph.ecount() else None, n_iterations=-1)
    membership = pd.Series(parts.membership, index=ids)
    by_size = sorted(membership.unique(), key=lambda c: (-int((membership == c).sum()), c))
    return membership.map({c: k for k, c in enumerate(by_size)})


def pack_communities(membership, positions):
    """Positions that put each community on its own circle, the circles tiled in rows
    over the box the nodes occupy now.

    A community's circle grows with its size so members sit ``NODE_SPACING`` apart; a
    lone member sits at its circle's centre.  Rows wrap at the box's current width,
    or at the widest circle when that is wider.  Returns ``{id: (x, y)}``.
    """
    ids = list(membership.index)
    x0 = min(positions[i][0] for i in ids)
    y0 = min(positions[i][1] for i in ids)
    circles = []
    for c in sorted(membership.unique()):
        members = [i for i in ids if membership[i] == c]
        radius = max(NODE_SPACING / 2, len(members) * NODE_SPACING / (2 * np.pi))
        circles.append((members, radius))
    width = max(max(positions[i][0] for i in ids) - x0, *(2 * r for _, r in circles))
    gap = NODE_SPACING / 2

    placed = {}
    x, y, row_height = x0, y0, 0.0
    for members, radius in circles:
        if x > x0 and x + 2 * radius > x0 + width:
            x, y, row_height = x0, y + row_height + gap, 0.0
        cx, cy = x + radius, y + radius
        if len(members) == 1:
            placed[members[0]] = (cx, cy)
        else:
            for k, i in enumerate(members):
                angle = 2 * np.pi * k / len(members)
                placed[i] = (cx + radius * np.cos(angle), cy + radius * np.sin(angle))
        x += 2 * radius + gap
        row_height = max(row_height, 2 * radius)
    return placed


def _spread(points, radii, gap, rounds=500):
    """Push points apart until each pair sits at least ``radii[a] + radii[b] + gap``
    apart.  Each round a point moves by the mean of its half-overlaps, so a point
    overlapping many others does not overshoot."""
    points = np.array(points, float)
    need = radii[:, None] + radii[None, :] + gap
    np.fill_diagonal(need, 0.0)
    for _ in range(rounds):
        v = points[:, None, :] - points[None, :, :]
        d = np.maximum(np.linalg.norm(v, axis=2), 1e-6)
        over = np.clip(need - d, 0.0, None)
        if not over.any():
            break
        hits = np.maximum((over > 0).sum(axis=1), 1)
        points += ((over / d)[:, :, None] * v).sum(axis=1) / (2 * hits[:, None])
    return points


def _shape(graph, method, seed, spacing=NODE_SPACING):
    """Centered coordinates of ``graph`` by the igraph layout ``CLUSTER_LAYOUTS[method]``,
    scaled so the median nearest-neighbor distance is ``spacing`` and spread so no two
    nodes sit closer than ``NODE_SPACING``.  Kamada-Kawai reads weights as lengths, so
    it takes their inverse."""
    import igraph
    from scipy.spatial import cKDTree

    if graph.vcount() == 1:
        return np.zeros((1, 2))
    weights = np.array(graph.es['weight'], float) if graph.ecount() else None
    if weights is not None and method == 'kamada-kawai':
        weights = 1 / weights
    igraph.set_random_number_generator(random.Random(seed))
    xy = np.array(getattr(graph, CLUSTER_LAYOUTS[method])(weights=weights).coords, float)
    xy -= xy.mean(axis=0)
    nearest = cKDTree(xy).query(xy, k=2)[0][:, 1]
    xy *= spacing / max(np.median(nearest), 1e-6)
    return _spread(xy, np.full(len(xy), NODE_SPACING / 2), 0.0)


def layout_communities(membership, edges, positions, method, seed=17, literature_weight=1.0):
    """Positions that lay each community out by ``method`` over its own edges inside a
    disc, the discs placed by the same layout over the summed weight between communities
    and pushed apart until none overlap.  The result keeps the top-left corner of the
    box the nodes occupy now.  Returns ``{id: (x, y)}``."""
    if method not in CLUSTER_LAYOUTS or CLUSTER_LAYOUTS[method] is None:
        raise ValueError(f"cluster layout must be one of {tuple(CLUSTER_LAYOUTS)}; got {method!r}")
    ids = list(membership.index)
    graph = _weighted_graph(edges, ids, literature_weight)
    communities = sorted(membership.unique())
    members = [[k for k, i in enumerate(ids) if membership[i] == c] for c in communities]
    shapes = [_shape(graph.induced_subgraph(m), method, seed) for m in members]
    radii = np.array([np.linalg.norm(shape, axis=1).max() + NODE_SPACING / 2 for shape in shapes])

    quotient = graph.copy()
    quotient.contract_vertices([communities.index(membership[i]) for i in ids])
    quotient.simplify(combine_edges='sum')
    centres = _spread(_shape(quotient, method, seed, 2 * np.median(radii)), radii, NODE_SPACING / 2)

    xy = np.zeros((len(ids), 2))
    for m, centre, shape in zip(members, centres, shapes):
        xy[m] = centre + shape
    xy += np.array([min(positions[i][0] for i in ids), min(positions[i][1] for i in ids)]) - xy.min(axis=0)
    return {i: (float(x), float(y)) for i, (x, y) in zip(ids, xy)}


def community_fill(membership):
    """A fill colour per clustered node, cycling through ``COMMUNITY_FILL``."""
    return membership.map(lambda c: COMMUNITY_FILL[int(c) % len(COMMUNITY_FILL)])
