"""Thin CyREST helpers under the Cytoscape tab.

Import this before py4cytoscape anywhere in the GUI: it fixes the CyREST address and
the detail-log directory py4cytoscape reads from the environment at import, and it
replaces py4cytoscape's request function with one that has a deadline.

Nothing here locks a visual property.  Positions and the viewport are written onto the
view as plain values so the mouse keeps working; ``unlock`` is the repair for a view
that arrived locked another way.
"""

import logging
import os
import tempfile

import requests

from log_config import get_logger

logger = get_logger(__name__)


def cytoscape_base_url():
    """Where CyREST answers: an explicit override, else the host from a container."""
    url = os.environ.get('PROXIMATE_CYTOSCAPE_URL')
    if url:
        return url.rstrip('/')
    in_docker = os.path.exists('/.dockerenv')
    return 'http://host.docker.internal:1234/v1' if in_docker else 'http://127.0.0.1:1234/v1'


BASE_URL = cytoscape_base_url()
CY_TIMEOUT = float(os.environ.get('PROXIMATE_CYTOSCAPE_TIMEOUT', '120'))
CY_RETRIES = 3

# py4cytoscape reads both of these once, at import.
os.environ.setdefault('DEFAULT_BASE_URL', BASE_URL)
os.environ.setdefault('PY4CYTOSCAPE_DETAIL_LOGGER_DIR', os.path.join(
    os.environ.get('PROXIMATE_LOG_DIR', tempfile.gettempdir()), 'py4cytoscape'))

import py4cytoscape as p4c  # noqa: E402
from py4cytoscape import commands as _commands  # noqa: E402


def harden_cyrest(timeout=CY_TIMEOUT, max_tries=CY_RETRIES):
    """Give every CyREST request a deadline and a short retry budget.

    py4cytoscape sends requests with no timeout and retries a refused connection ten
    times under exponential backoff, which is minutes of a frozen handler when the
    desktop is down.  ``commands._do_request`` resolves this name at call time, so
    replacing the attribute routes every call through it.
    """
    import backoff

    @backoff.on_exception(backoff.expo, requests.exceptions.ConnectionError,
                          max_tries=max_tries)
    def _do_request_local(method, url, **kwargs):
        kwargs.setdefault('timeout', timeout)
        _commands.log_http_request(method, url, **kwargs)
        r = requests.request(method, url, **kwargs)
        _commands.log_http_result(r)
        return r

    _commands._do_request_local = _do_request_local


def quiet_p4c():
    """Keep py4cytoscape's DEBUG detail on its own file handler, off ours."""
    for name in ('py4...', 'py4...S'):
        logging.getLogger(name).propagate = False


harden_cyrest()
quiet_p4c()


def probe(timeout=5):
    """Whether CyREST answers at BASE_URL, without going through py4cytoscape."""
    try:
        r = requests.get(BASE_URL + '/version', timeout=timeout)
        r.raise_for_status()
        version = r.json().get('cytoscapeVersion')
        return {'ok': True, 'url': BASE_URL, 'version': version, 'error': None}
    except (requests.RequestException, ValueError) as e:
        return {'ok': False, 'url': BASE_URL, 'version': None, 'error': str(e)}


# --- networks and views ----------------------------------------------------------------

def node_suids(network):
    names = p4c.get_table_columns('node', ['name'], network=network)['name']
    return {name: suid for suid, name in names.items()}


def view_suid(network):
    return p4c.get_network_views(network=network)[0]


def current_positions(network):
    """{node name: (x, y)} as Cytoscape has them now."""
    pos = p4c.get_node_position(network=network).astype(float)
    return {str(name): (float(row.x), float(row.y)) for name, row in pos.iterrows()}


def selected_nodes(network):
    """Names of the selected nodes; an empty list when nothing is selected.

    py4cytoscape answers None rather than [] for an empty selection.
    """
    return [str(n) for n in (p4c.get_selected_nodes(node_suids=False, network=network) or [])]


def edge_name(source, interaction, target):
    """The name Cytoscape gives an edge created from a data frame."""
    return f'{source} ({interaction}) {target}'


def update_edge_columns(network, frame):
    """Push columns onto existing edges in place, keyed on the edge name.

    ``frame`` carries ``name`` plus the columns to set.  Nothing moves and nothing is
    rebuilt; a passthrough mapping on the column does the rest.
    """
    if not len(frame):
        return 0
    p4c.load_table_data(frame, data_key_column='name', table='edge',
                        table_key_column='name', network=network)
    return len(frame)


def update_node_columns(network, frame):
    """Push columns onto existing nodes in place, keyed on the node name."""
    if not len(frame):
        return 0
    p4c.load_table_data(frame, data_key_column='name', table='node',
                        table_key_column='name', network=network)
    return len(frame)


def apply_passthrough_style(style, network, defaults, passthrough):
    """Recreate ``style`` with one passthrough mapping per (property, column).

    Every look is computed into a table column in Python; Cytoscape only maps.  A
    style cannot be deleted while applied, so the network steps aside first.
    """
    if style in p4c.get_visual_style_names():
        p4c.set_visual_style('default', network=network)
        p4c.delete_visual_style(style)
    p4c.create_visual_style(style, defaults=defaults)
    p4c.set_visual_style(style, network=network)
    for prop, column in passthrough:
        p4c.update_style_mapping(style, p4c.map_visual_property(prop, column, 'p'))


# --- locks -------------------------------------------------------------------------------

# The view properties Cytoscape will hold a locked value for.  A locked
# NETWORK_SCALE_FACTOR is a view the mouse wheel cannot zoom.
VIEW_LOCKS = ('NETWORK_SCALE_FACTOR', 'NETWORK_CENTER_X_LOCATION', 'NETWORK_CENTER_Y_LOCATION')


def locked_view_properties(network):
    """The view properties currently holding a locked value.

    A lock shows nothing on screen; the bypass endpoint answers 200 when one is set
    and 404 when none is.
    """
    view = view_suid(network)
    return [vp for vp in VIEW_LOCKS
            if requests.get(f'{BASE_URL}/networks/{network}/views/{view}/network/{vp}/bypass',
                            timeout=CY_TIMEOUT).status_code == 200]


def unlock(network):
    """Release any view lock.  Only the properties holding one are cleared."""
    held = locked_view_properties(network)
    for vp in held:
        p4c.clear_network_property_bypass(vp, network=network)
    return held


# --- export --------------------------------------------------------------------------------

def export_png(network, path, height=2000):
    """Write the whole network as a PNG of the given pixel height.

    CyREST's view-PNG endpoint fits the network to the height and returns the bytes,
    so the image is written from this process and no path is handed to Cytoscape.
    ``export_image`` would rasterize the window instead, capped by the screen.
    """
    p4c.fit_content(network=network)
    r = requests.get(f'{BASE_URL}/networks/{network}/views/first.png',
                     params={'h': int(height)}, timeout=CY_TIMEOUT)
    r.raise_for_status()
    with open(path, 'wb') as fh:
        fh.write(r.content)
    return path
