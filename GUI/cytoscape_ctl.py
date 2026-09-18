"""The Cytoscape controller: one state, one lock, every operation on the drawn network.

A Cytoscape desktop is a single shared thing, so the state lives at module level
rather than in a Shiny session: two browser tabs see the same network, and every
mutation bumps ``version`` for the GUI to poll.  Nothing outside this module touches
Cytoscape or the state.

Lock order: ``cy_lock`` serializes every touch of Cytoscape; ``lock`` guards the
dict.  Take ``cy_lock`` first when both are needed.
"""

import collections
import contextlib
import datetime
import os
import threading

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
    'log': collections.deque(maxlen=200),   # {'ts', 'op', 'detail'}
    'dataset': None,
    'title': None,
    'net_suid': None,
    'nodes': None,        # DataFrame as drawn
    'edges': None,        # DataFrame as drawn, 'visible' kept current
    'thresholds': None,   # the thresholds the edges were drawn at
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
def _mutate(op, detail=''):
    with STATE['lock']:
        yield STATE
        STATE['version'] += 1
        STATE['log'].append({'ts': datetime.datetime.now().isoformat(timespec='seconds'),
                             'op': op, 'detail': detail})
    logger.info("cytoscape %s: %s", op, detail)


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
         label_policy='all', layout='force-directed'):
    """Build the thresholded network and draw it, replacing the previous one."""
    df = pd.read_csv(scores_path)
    nodes, edges = net.build(df, thresholds, baits=baits, prey_prey=prey_prey,
                             biogrid_path=biogrid_path, label_policy=label_policy)
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
    with _mutate('draw', f'{title}: {len(nodes)} nodes, {len(edges)} edges, layout {layout}'):
        STATE.update(dataset=dataset, title=title, net_suid=suid, nodes=nodes, edges=edges,
                     thresholds=dict(thresholds))
    return snapshot()


def read_selection():
    """The selected nodes with their symbols, roles and the scores of their edges."""
    with STATE['cy_lock']:
        selected = cy.selected_nodes(_net())
    if not selected:
        raise ValueError("nothing is selected in Cytoscape")
    nodes, edges = STATE['nodes'], STATE['edges']
    chosen = nodes[nodes['id'].isin(selected)][['id', 'symbol', 'role']]
    touching = edges[edges['source'].isin(selected) | edges['target'].isin(selected)]
    score_cols = [c for c in ('SaintScore', 'BFDR', 'FoldChange', 'WD', 'WDFDR') if c in edges.columns]
    detail = touching[['source', 'target', 'interaction', 'visible', *score_cols]]
    return chosen.reset_index(drop=True), detail.reset_index(drop=True)


def set_edge_visibility(action):
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
    with _mutate('set_edge_visibility', f'{action}: {len(changed)} edge(s) changed'):
        pass
    return len(changed)


def apply_thresholds(thresholds):
    """Re-threshold the drawn network in place, keeping the user's layout.

    Only thresholds at least as tight as the ones drawn can be applied this way; a
    looser set needs edges that were never drawn, so it asks for a new send.
    """
    suid = _net()
    if not net.is_tighter_or_equal(thresholds, STATE['thresholds']):
        raise ValueError("these thresholds are looser than the drawn network's; "
                         "send the network again to add interactions")
    edges = STATE['edges']
    want = net.edge_visibility(edges, thresholds)
    changed = edges.index[want != edges['visible'].astype(bool)]
    _push_visibility(suid, want, changed)
    with _mutate('apply_thresholds', f'{int(want.sum())} of {len(edges)} edges visible'):
        STATE['thresholds'] = dict(thresholds)
    return int((~want).sum())


def _push_visibility(suid, want, changed):
    edges = STATE['edges']
    frame = pd.DataFrame({'name': edges.loc[changed, 'name'].to_numpy(),
                          'visible': want[changed].astype(bool).to_numpy()})
    with STATE['lock']:
        edges['visible'] = want.astype(bool).to_numpy()
        alpha = net.node_alpha(STATE['nodes'], edges)
        STATE['nodes']['node_alpha'] = alpha.to_numpy()
    with STATE['cy_lock']:
        cy.update_edge_columns(suid, frame)
        cy.update_node_columns(suid, pd.DataFrame({'name': STATE['nodes']['id'].to_numpy(),
                                                   'node_alpha': alpha.astype(int).to_numpy()}))


def sync_positions():
    """Record the positions the user dragged nodes to."""
    suid = _net()
    with STATE['cy_lock']:
        positions = cy.current_positions(suid)
    with _mutate('sync_positions', f'{len(positions)} node(s)'):
        nodes = STATE['nodes']
        nodes['x'] = nodes['id'].map(lambda i: positions.get(i, (None, None))[0])
        nodes['y'] = nodes['id'].map(lambda i: positions.get(i, (None, None))[1])
    return len(positions)


def export_image(out_dir, height=2000):
    """Write the drawn network as a PNG under ``<out_dir>/cytoscape/``."""
    suid = _net()
    target = os.path.join(out_dir, 'cytoscape')
    os.makedirs(target, exist_ok=True)
    stamp = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')
    path = os.path.join(target, f'network_{stamp}.png')
    with _busy('exporting image'), STATE['cy_lock']:
        cy.export_png(suid, path, height=height)
    with _mutate('export_image', path):
        pass
    return path


def unlock():
    """Release any lock a view picked up, restoring pan and zoom."""
    with STATE['cy_lock']:
        held = cy.unlock(_net())
    with _mutate('unlock', ', '.join(held) or 'nothing was locked'):
        pass
    return held
