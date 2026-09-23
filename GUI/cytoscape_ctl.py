"""The Cytoscape controller: one state, one lock, every operation on the drawn network.

A Cytoscape desktop is a single shared thing, so the state lives at module level
rather than in a Shiny session: two browser tabs and the MCP server see the same
network, and every mutation bumps ``version`` for the GUI to poll.  Nothing outside
this module touches Cytoscape or the state.

Every mutating operation takes an ``actor`` ('gui' or 'mcp') that lands in the activity
log, so a person at the GUI can see what an agent did to the drawn network.

Lock order: ``cy_lock`` serializes every touch of Cytoscape; ``lock`` guards the
dict.  Take ``cy_lock`` first when both are needed.
"""

import collections
import contextlib
import datetime
import os
import threading

import numpy as np
import pandas as pd

import cytoscape_net as net
import cytoscape_p4c as cy
from cytoscape_p4c import p4c
from log_config import get_logger

logger = get_logger(__name__)

STYLE = 'ProxiMate'
COLLECTION = 'ProxiMate'

STATE = {
    'lock': threading.RLock(),
    'cy_lock': threading.Lock(),
    'version': 0,
    'log': collections.deque(maxlen=200),   # {'ts', 'actor', 'op', 'detail'}
    'dataset': None,
    'title': None,
    'net_suid': None,
    'nodes': None,        # DataFrame as drawn
    'edges': None,        # DataFrame as drawn, 'visible' kept current
    'thresholds': None,   # the thresholds the edges were drawn at
    'style': None,        # {'width_source', 'literature_weighted', 'biogrid_scope'} as drawn
    'busy': '',
}

PASSTHROUGH = [
    ('NODE_LABEL', 'display_label'), ('NODE_SHAPE', 'shape'), ('NODE_SIZE', 'size'),
    ('NODE_FILL_COLOR', 'fill'), ('NODE_BORDER_WIDTH', 'border'),
    ('NODE_BORDER_PAINT', 'border_color'), ('NODE_LABEL_FONT_SIZE', 'label_size'),
    ('NODE_TRANSPARENCY', 'node_alpha'),
    ('EDGE_WIDTH', 'width'), ('EDGE_STROKE_UNSELECTED_PAINT', 'color'),
    ('EDGE_TRANSPARENCY', 'alpha'), ('EDGE_VISIBLE', 'visible'),
]
STYLE_DEFAULTS = {
    'NODE_SHAPE': 'ELLIPSE', 'NODE_SIZE': net.NODE_SIZE['prey'],
    'NODE_FILL_COLOR': net.NODE_FILL['prey'], 'NODE_BORDER_WIDTH': net.NODE_BORDER['prey'],
    'NODE_BORDER_PAINT': net.NODE_BORDER_COLOR, 'NODE_LABEL_FONT_SIZE': net.LABEL_SIZE['prey'],
    'NODE_LABEL_COLOR': '#1F2933', 'EDGE_WIDTH': 2.0, 'NETWORK_BACKGROUND_PAINT': '#FFFFFF',
}

VISIBILITY_ACTIONS = ('hide_selected', 'show_selected', 'hide_unselected', 'show_all')


@contextlib.contextmanager
def _mutate(op, detail='', actor='gui'):
    with STATE['lock']:
        yield STATE
        STATE['version'] += 1
        STATE['log'].append({'ts': datetime.datetime.now().isoformat(timespec='seconds'),
                             'actor': actor, 'op': op, 'detail': detail})
    logger.info("cytoscape %s (%s): %s", op, actor, detail)


@contextlib.contextmanager
def _busy(text):
    STATE['busy'] = text
    try:
        yield
    finally:
        STATE['busy'] = ''


def _net():
    n = STATE['net_suid']
    if n is None:
        raise RuntimeError("no ProxiMate network is in Cytoscape; send one first")
    return n


def snapshot():
    """The scalar state, copied under the lock, for the tab to render."""
    with STATE['lock']:
        edges = STATE['edges']
        return {
            'version': STATE['version'], 'dataset': STATE['dataset'], 'title': STATE['title'],
            'net_suid': STATE['net_suid'], 'busy': STATE['busy'],
            'n_nodes': 0 if STATE['nodes'] is None else len(STATE['nodes']),
            'n_edges': 0 if edges is None else len(edges),
            'n_hidden': 0 if edges is None else int((~edges['visible'].astype(bool)).sum()),
            'thresholds': dict(STATE['thresholds']) if STATE['thresholds'] else None,
            'style': dict(STATE['style']) if STATE['style'] else None,
            'log': list(STATE['log'])[-20:],
        }


def health():
    """Whether Cytoscape answers, without touching the state."""
    return cy.probe()


def layout_names():
    with STATE['cy_lock']:
        return sorted(p4c.get_layout_names())


# --- operations ------------------------------------------------------------------------

def draw(dataset, scores_path, thresholds, baits=None, prey_prey=True, biogrid_path=None,
         label_policy='all', layout='force-directed', width_source='abundance',
         literature_weighted=False, biogrid_scope='all', corum_path=None,
         corum_min_members=3, corum_min_fraction=0.5, actor='gui'):
    """Build the thresholded network and draw it, replacing the previous one."""
    df = pd.read_csv(scores_path)
    nodes, edges = net.build(df, thresholds, baits=baits, prey_prey=prey_prey,
                             biogrid_path=biogrid_path, label_policy=label_policy,
                             width_source=width_source, literature_weighted=literature_weighted,
                             corum_path=corum_path, corum_min_members=corum_min_members,
                             corum_min_fraction=corum_min_fraction)
    edges['visible'] = net.edge_visibility(edges, thresholds, biogrid_scope).to_numpy()
    nodes['node_alpha'] = net.node_alpha(nodes, edges).to_numpy()
    style = {'width_source': width_source, 'literature_weighted': bool(literature_weighted),
             'biogrid_scope': biogrid_scope}
    title = f'ProxiMate: {dataset}'
    with _busy(f'drawing {title}'), STATE['cy_lock']:
        if title in p4c.get_network_list():
            p4c.delete_network(p4c.get_network_suid(title))
        suid = p4c.create_network_from_data_frames(
            nodes, edges.drop(columns=['name']), title=title, collection=COLLECTION)
        cy.apply_passthrough_style(STYLE, suid, STYLE_DEFAULTS, PASSTHROUGH)
        p4c.layout_network(layout, network=suid)
        cy.unlock(suid)
        p4c.fit_content(network=suid)
    with _mutate('draw', f'{title}: {len(nodes)} nodes, {len(edges)} edges, layout {layout}',
                 actor=actor):
        STATE.update(dataset=dataset, title=title, net_suid=suid, nodes=nodes, edges=edges,
                     thresholds=dict(thresholds), style=style)
    return snapshot()


def restyle_edges(width_source, literature_weighted, biogrid_scope, actor='gui'):
    """Recompute edge widths and the BioGRID scope on the drawn network, in place."""
    suid = _net()
    edges = STATE['edges']
    width = net.edge_widths(edges, width_source, literature_weighted)
    want = net.edge_visibility(edges, STATE['thresholds'], biogrid_scope)
    changed = edges.index[(want != edges['visible'].astype(bool)) | (width != edges['width'])]
    with STATE['lock']:
        edges['width'] = width.to_numpy()
    _push_visibility(suid, want, changed, extra_columns=['width'])
    with _mutate('restyle_edges', f'width by {width_source}, literature '
                 f'{"weighted" if literature_weighted else "unweighted"}, BioGRID {biogrid_scope}: '
                 f'{len(changed)} edge(s) changed', actor=actor):
        STATE['style'] = {'width_source': width_source, 'literature_weighted': bool(literature_weighted),
                          'biogrid_scope': biogrid_scope}
    return len(changed)


def read_selection():
    """The selected nodes with their symbols, roles and the scores of their edges."""
    with STATE['cy_lock']:
        selected = cy.selected_nodes(_net())
    if not selected:
        raise ValueError("nothing is selected in Cytoscape")
    nodes, edges = STATE['nodes'], STATE['edges']
    chosen = nodes[nodes['id'].isin(selected)][['id', 'symbol', 'accession', 'role']]
    touching = edges[edges['source'].isin(selected) | edges['target'].isin(selected)]
    score_cols = [c for c in ('SaintScore', 'BFDR', 'FoldChange', 'WD', 'WDFDR') if c in edges.columns]
    detail = touching[['source', 'target', 'interaction', 'visible', *score_cols]]
    return chosen.reset_index(drop=True), detail.reset_index(drop=True)


def set_edge_visibility(action, actor='gui'):
    """Hide or show edges against the Cytoscape selection, as a column update.

    ``hide_selected``/``show_selected`` act on edges touching a selected node;
    ``hide_unselected`` hides every other edge; ``show_all`` restores everything.
    """
    if action not in VISIBILITY_ACTIONS:
        raise ValueError(f"action must be one of {', '.join(VISIBILITY_ACTIONS)}")
    suid = _net()
    edges = STATE['edges']
    if action == 'show_all':
        want = pd.Series(True, index=edges.index)
    else:
        with STATE['cy_lock']:
            selected = cy.selected_nodes(suid)
        if not selected:
            raise ValueError("nothing is selected in Cytoscape")
        touching = edges['source'].isin(selected) | edges['target'].isin(selected)
        current = edges['visible'].astype(bool)
        want = {'hide_selected': current & ~touching,
                'show_selected': current | touching,
                'hide_unselected': current & touching}[action]
    changed = edges.index[want != edges['visible'].astype(bool)]
    _push_visibility(suid, want, changed)
    with _mutate('set_edge_visibility', f'{action}: {len(changed)} edge(s) changed', actor=actor):
        pass
    return len(changed)


def apply_thresholds(thresholds, actor='gui'):
    """Re-threshold the drawn network in place, keeping the user's layout.

    Only thresholds at least as tight as the ones drawn can be applied this way; a
    looser set needs edges that were never drawn, so it asks for a new send.
    """
    suid = _net()
    if not net.is_tighter_or_equal(thresholds, STATE['thresholds']):
        raise ValueError("these thresholds are looser than the drawn network's; "
                         "send the network again to add interactions")
    edges = STATE['edges']
    want = net.edge_visibility(edges, thresholds, STATE['style']['biogrid_scope'])
    changed = edges.index[want != edges['visible'].astype(bool)]
    _push_visibility(suid, want, changed)
    with _mutate('apply_thresholds', f'{int(want.sum())} of {len(edges)} edges visible', actor=actor):
        STATE['thresholds'] = dict(thresholds)
    return int((~want).sum())


def _push_visibility(suid, want, changed, extra_columns=()):
    """Set ``visible`` to ``want`` on the ``changed`` edges, plus any ``extra_columns``
    already updated in the state, and refade the nodes."""
    edges = STATE['edges']
    frame = pd.DataFrame({'name': edges.loc[changed, 'name'].to_numpy(),
                          'visible': want[changed].astype(bool).to_numpy()})
    for column in extra_columns:
        frame[column] = edges.loc[changed, column].to_numpy()
    with STATE['lock']:
        edges['visible'] = want.astype(bool).to_numpy()
        alpha = net.node_alpha(STATE['nodes'], edges)
        STATE['nodes']['node_alpha'] = alpha.to_numpy()
    with STATE['cy_lock']:
        cy.update_edge_columns(suid, frame)
        cy.update_node_columns(suid, pd.DataFrame({'name': STATE['nodes']['id'].to_numpy(),
                                                   'node_alpha': alpha.astype(int).to_numpy()}))


def sync_positions(actor='gui'):
    """Record the positions the user dragged nodes to."""
    suid = _net()
    with STATE['cy_lock']:
        positions = cy.current_positions(suid)
    with _mutate('sync_positions', f'{len(positions)} node(s)', actor=actor):
        nodes = STATE['nodes']
        nodes['x'] = nodes['id'].map(lambda i: positions.get(i, (None, None))[0])
        nodes['y'] = nodes['id'].map(lambda i: positions.get(i, (None, None))[1])
    return len(positions)


# --- selection tools ---------------------------------------------------------------------

def _selected_bait():
    """The one bait in the Cytoscape selection; refuses none or several."""
    with STATE['cy_lock']:
        selected = cy.selected_nodes(_net())
    nodes = STATE['nodes']
    chosen = nodes[nodes['id'].isin(selected)]
    baits = chosen.loc[chosen['role'] == 'bait', 'id'].tolist()
    if len(baits) != 1:
        named = ', '.join(chosen.loc[chosen['role'] == 'bait', 'symbol']) or 'no bait'
        raise ValueError(f"select exactly one bait in Cytoscape (selected: {named})")
    return baits[0]


def _select(op, ids, detail, add=False, actor='gui'):
    with STATE['cy_lock']:
        cy.select_nodes(_net(), ids, add=add)
    with _mutate(op, detail, actor=actor):
        pass
    return ids


def _resolve_ids(queries):
    """Drawn node ids for a list of accessions or symbols; unknown names raise."""
    if not queries:
        raise ValueError("no nodes named")
    nodes = STATE['nodes']
    return [str(net.resolve_node(nodes, q)['id']) for q in queries]


def select_nodes(ids, add=False, actor='gui'):
    """Make these nodes (accessions or symbols) the Cytoscape selection; ``add`` keeps
    what is already selected.  Returns the resolved ids."""
    _net()
    chosen = _resolve_ids(ids)
    return _select('select_nodes', chosen, f'{len(chosen)} node(s)' + (' added' if add else ''),
                   add=add, actor=actor)


def list_nodes(role=None):
    """The drawn nodes' id, accession, symbol and role, all or one role only."""
    _net()
    nodes = STATE['nodes'][['id', 'accession', 'symbol', 'role']]
    if role is not None:
        if role not in ('bait', 'prey'):
            raise ValueError(f"role must be bait or prey, got {role!r}")
        nodes = nodes[nodes['role'] == role]
    return nodes.astype(str).to_dict('records')


def get_positions(ids=None):
    """``{id: [x, y]}`` as Cytoscape has the drawn nodes now, all or the named ones."""
    suid = _net()
    with STATE['cy_lock']:
        positions = cy.current_positions(suid)
    drawn = [str(i) for i in STATE['nodes']['id']]
    wanted = _resolve_ids(ids) if ids else drawn
    missing = [i for i in wanted if i not in positions]
    if missing:
        raise ValueError(f"Cytoscape reported no position for {', '.join(missing[:8])}")
    return {i: [positions[i][0], positions[i][1]] for i in wanted}


def move_nodes(positions, actor='gui'):
    """Move nodes to absolute ``{id or symbol: [x, y]}`` positions as plain values, so
    they stay draggable; the state keeps the new coordinates.  Returns the count."""
    suid = _net()
    targets = {}
    for query, xy in dict(positions).items():
        if not isinstance(xy, (list, tuple)) or len(xy) != 2:
            raise ValueError(f"position for {query!r} must be [x, y], got {xy!r}")
        node_id = _resolve_ids([query])[0]
        targets[node_id] = (float(xy[0]), float(xy[1]))
    with STATE['cy_lock']:
        cy.set_positions(suid, targets)
    with _mutate('move_nodes', f'{len(targets)} node(s)', actor=actor):
        nodes = STATE['nodes']
        chosen = nodes['id'].isin(targets)
        nodes.loc[chosen, 'x'] = nodes.loc[chosen, 'id'].map(lambda i: targets[i][0]).to_numpy()
        nodes.loc[chosen, 'y'] = nodes.loc[chosen, 'id'].map(lambda i: targets[i][1]).to_numpy()
    return len(targets)


def select_loners(actor='gui'):
    """Select the one selected bait together with its loners, so they move as one."""
    bait = _selected_bait()
    found = net.loners(STATE['nodes'], STATE['edges'], bait)
    symbol = _symbols()[bait]
    if not found:
        raise ValueError(f"{symbol} has no loners")
    return _select('select_loners', [bait] + found, f'{symbol}: {len(found)} loner(s)', actor=actor)


def select_satellites(actor='gui'):
    """Select the one selected bait with its satellites: the preys whose only bait it
    is, and the two-bait preys that currently sit nearer to it than to the other."""
    bait = _selected_bait()
    with STATE['cy_lock']:
        positions = cy.current_positions(_net())
    only, nearer = net.satellites(STATE['nodes'], STATE['edges'], bait, positions)
    symbol = _symbols()[bait]
    if not only and not nearer:
        raise ValueError(f"{symbol} has no satellites")
    return _select('select_satellites', [bait] + only + nearer,
                   f'{symbol}: {len(only)} only-bait prey(s), {len(nearer)} nearer of two', actor=actor)


def select_related(seed, relation, add=False, actor='gui', **cuts):
    """Select what ``relation`` names for ``seed`` (``net.related``); ``add`` keeps the
    current selection."""
    _net()
    ids = net.related(STATE['nodes'], STATE['edges'], seed, relation, **cuts)
    if not ids:
        raise ValueError(f"no {relation} for {seed} under these cuts")
    cut = ', '.join(f'{k} {v}' for k, v in cuts.items() if v is not None) or 'no cuts'
    return _select('select_related', ids, f'{relation} of {seed} ({cut}): {len(ids)}'
                   + (' added' if add else ''), add=add, actor=actor)


def _symbols():
    return dict(zip(STATE['nodes']['id'], STATE['nodes']['symbol']))


def cluster_selection(resolution=1.0, seed=17, literature_weight=1.0, actor='gui'):
    """Leiden over the selected nodes, then re-pack them by community inside the box
    they occupy.  Only the selected nodes move; they recolour by community and carry
    it in the node table."""
    suid = _net()
    with STATE['cy_lock']:
        selected = cy.selected_nodes(suid)
        positions = cy.current_positions(suid)
    if not selected:
        raise ValueError("nothing is selected in Cytoscape")
    with _busy(f'clustering {len(selected)} selected nodes'):
        membership = net.cluster(STATE['edges'], selected, resolution, seed, literature_weight)
        placed = net.pack_communities(membership, positions)
    nodes = STATE['nodes']
    # Numbers continue above any community already assigned, so two clustered regions
    # never share one.
    if 'community' in nodes.columns:
        membership = membership + int(nodes['community'].fillna(-1).max()) + 1
    fill = net.community_fill(membership)
    with STATE['lock']:
        if 'community' not in nodes.columns:
            nodes['community'] = np.nan
        chosen = nodes['id'].isin(membership.index)
        nodes.loc[chosen, 'community'] = nodes.loc[chosen, 'id'].map(membership).to_numpy()
        nodes.loc[chosen, 'fill'] = nodes.loc[chosen, 'id'].map(fill).to_numpy()
        nodes.loc[chosen, 'x'] = nodes.loc[chosen, 'id'].map(lambda i: placed[i][0]).to_numpy()
        nodes.loc[chosen, 'y'] = nodes.loc[chosen, 'id'].map(lambda i: placed[i][1]).to_numpy()
    with STATE['cy_lock']:
        cy.update_node_columns(suid, pd.DataFrame({'name': membership.index,
                                                   'community': membership.astype(int).to_numpy(),
                                                   'fill': fill.to_numpy()}))
        cy.set_positions(suid, placed)
    sizes = membership.value_counts().sort_index().tolist()
    with _mutate('cluster_selection', f'{len(membership)} nodes -> {len(sizes)} communities '
                 f'{sizes} (resolution {resolution}, seed {seed})', actor=actor):
        pass
    return {'n': len(membership), 'n_communities': len(sizes), 'sizes': sizes}


def export_image(out_dir, height=2000, actor='gui'):
    """Write the drawn network as a PNG under ``<out_dir>/cytoscape/``."""
    suid = _net()
    target = os.path.join(out_dir, 'cytoscape')
    os.makedirs(target, exist_ok=True)
    stamp = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
    path = os.path.join(target, f'network_{stamp}.png')
    with _busy('exporting image'), STATE['cy_lock']:
        cy.export_png(suid, path, height=height)
    with _mutate('export_image', path, actor=actor):
        pass
    return path


def unlock(actor='gui'):
    """Release any lock a view picked up, restoring pan and zoom."""
    with STATE['cy_lock']:
        held = cy.unlock(_net())
    with _mutate('unlock', ', '.join(held) or 'nothing was locked', actor=actor):
        pass
    return held
