"""Tests for the CyREST helpers, with py4cytoscape and requests stubbed out.

What matters is the shape of what reaches Cytoscape and, above all, that nothing here
ever sets a bypass: a locked view shows nothing on screen and stops the mouse.
"""

import pytest

import cytoscape_p4c as cy


BYPASS_SETTERS = ('set_node_position_bypass', 'set_network_zoom_bypass',
                  'set_network_center_bypass', 'set_node_property_bypass',
                  'set_edge_property_bypass', 'set_network_property_bypass')


@pytest.fixture(autouse=True)
def no_bypass(monkeypatch):
    """Every bypass setter raises, so any code path that reaches one fails the test."""
    def boom(*a, **k):
        raise AssertionError("a bypass setter was called")
    for name in BYPASS_SETTERS:
        if hasattr(cy.p4c, name):
            monkeypatch.setattr(cy.p4c, name, boom)


class _Response:
    def __init__(self, code=200, content=b'', payload=None):
        self.status_code = code
        self.content = content
        self._payload = payload or {}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise cy.requests.HTTPError(str(self.status_code))

    def json(self):
        return self._payload


# --- address ------------------------------------------------------------------------------

def test_an_explicit_url_wins(monkeypatch):
    monkeypatch.setenv('PROXIMATE_CYTOSCAPE_URL', 'http://cy.example:1234/v1/')

    assert cy.cytoscape_base_url() == 'http://cy.example:1234/v1'


def test_a_container_reaches_the_host(monkeypatch):
    monkeypatch.delenv('PROXIMATE_CYTOSCAPE_URL', raising=False)
    monkeypatch.setattr(cy.os.path, 'exists', lambda p: p == '/.dockerenv')

    assert cy.cytoscape_base_url() == 'http://host.docker.internal:1234/v1'


def test_the_probe_reports_a_down_desktop_as_an_error_not_an_exception(monkeypatch):
    def refused(*a, **k):
        raise cy.requests.ConnectionError('refused')
    monkeypatch.setattr(cy.requests, 'get', refused)

    result = cy.probe()

    assert result['ok'] is False and 'refused' in result['error']


def test_the_probe_reads_the_version(monkeypatch):
    monkeypatch.setattr(cy.requests, 'get',
                        lambda *a, **k: _Response(payload={'cytoscapeVersion': '3.10.2'}))

    assert cy.probe()['version'] == '3.10.2'


# --- selection and edges ---------------------------------------------------------------------

def test_an_empty_selection_is_an_empty_list(monkeypatch):
    monkeypatch.setattr(cy.p4c, 'get_selected_nodes', lambda **k: None)

    assert cy.selected_nodes(1) == []


def test_edge_names_follow_cytoscapes_convention():
    assert cy.edge_name('PA', 'proximity', 'P1') == 'PA (proximity) P1'


def test_selecting_clears_first_unless_adding(monkeypatch):
    calls = []
    monkeypatch.setattr(cy.p4c, 'clear_selection', lambda **k: calls.append('clear'))
    monkeypatch.setattr(cy.p4c, 'select_nodes', lambda names, **k: calls.append(('select', names, k['preserve_current_selection'])))

    cy.select_nodes(2, ['a', 'b'])
    cy.select_nodes(2, ['c'], add=True)

    assert calls == ['clear', ('select', ['a', 'b'], False), ('select', ['c'], True)]


def test_positions_are_written_as_plain_view_values(monkeypatch):
    seen = {}
    monkeypatch.setattr(cy, 'node_suids', lambda net: {'a': 11, 'b': 12})
    monkeypatch.setattr(cy, 'view_suid', lambda net: 7)
    monkeypatch.setattr(cy.requests, 'put', lambda url, **k: seen.update(url=url, **k) or _Response())

    n = cy.set_positions(3, {'a': (1.5, 2), 'b': (-3, 4)})

    assert n == 2
    assert seen['url'].endswith('/networks/3/views/7/nodes')
    assert seen['json'][0] == {'SUID': 11, 'view': [{'visualProperty': 'NODE_X_LOCATION', 'value': 1.5},
                                                     {'visualProperty': 'NODE_Y_LOCATION', 'value': 2.0}]}
    assert 'bypass' not in seen['url']


def test_edge_columns_are_pushed_keyed_on_name(monkeypatch):
    import pandas as pd
    seen = {}
    monkeypatch.setattr(cy.p4c, 'load_table_data', lambda frame, **k: seen.update(k, n=len(frame)))

    n = cy.update_edge_columns(5, pd.DataFrame({'name': ['a (x) b'], 'visible': [False]}))

    assert n == 1
    assert seen['table'] == 'edge' and seen['table_key_column'] == 'name' and seen['network'] == 5


# --- style ---------------------------------------------------------------------------------

def test_an_applied_style_steps_aside_before_it_is_replaced(monkeypatch):
    calls = []
    monkeypatch.setattr(cy.p4c, 'get_visual_style_names', lambda: ['ProxiMate'])
    for name in ('set_visual_style', 'delete_visual_style', 'create_visual_style'):
        monkeypatch.setattr(cy.p4c, name, lambda *a, _n=name, **k: calls.append((_n, a)))
    monkeypatch.setattr(cy.p4c, 'map_visual_property', lambda prop, col, kind: (prop, col, kind))
    monkeypatch.setattr(cy.p4c, 'update_style_mapping', lambda style, m: calls.append(('map', m)))

    cy.apply_passthrough_style('ProxiMate', 9, {'NODE_SIZE': 1}, [('NODE_LABEL', 'display_label')])

    assert calls[0] == ('set_visual_style', ('default',))
    assert calls[1] == ('delete_visual_style', ('ProxiMate',))
    assert ('map', ('NODE_LABEL', 'display_label', 'p')) in calls


# --- locks ------------------------------------------------------------------------------------

def test_a_lock_is_read_from_the_bypass_endpoints_status(monkeypatch):
    monkeypatch.setattr(cy, 'view_suid', lambda net: 7)
    monkeypatch.setattr(cy.requests, 'get', lambda url, **k: _Response(
        200 if url.endswith('NETWORK_SCALE_FACTOR/bypass') else 404))

    assert cy.locked_view_properties(3) == ['NETWORK_SCALE_FACTOR']


def test_unlock_clears_only_the_locks_that_are_held(monkeypatch):
    cleared = []
    monkeypatch.setattr(cy, 'locked_view_properties', lambda net: ['NETWORK_SCALE_FACTOR'])
    monkeypatch.setattr(cy.p4c, 'clear_network_property_bypass',
                        lambda vp, **k: cleared.append(vp))

    assert cy.unlock(3) == ['NETWORK_SCALE_FACTOR']
    assert cleared == ['NETWORK_SCALE_FACTOR']


def test_unlock_leaves_a_free_view_alone(monkeypatch):
    monkeypatch.setattr(cy, 'locked_view_properties', lambda net: [])

    def boom(*a, **k):
        raise AssertionError("cleared a lock that was not held")
    monkeypatch.setattr(cy.p4c, 'clear_network_property_bypass', boom)

    assert cy.unlock(3) == []


# --- export -----------------------------------------------------------------------------------

def test_the_export_fetches_the_bytes_at_the_requested_height(monkeypatch, tmp_path):
    seen = {}
    monkeypatch.setattr(cy.p4c, 'fit_content', lambda **k: None)
    monkeypatch.setattr(cy.requests, 'get',
                        lambda url, **k: seen.update(url=url, **k) or _Response(content=b'PNG'))
    path = tmp_path / 'net.png'

    cy.export_png(4, str(path), height=1500)

    assert path.read_bytes() == b'PNG'
    assert seen['url'].endswith('/networks/4/views/first.png') and seen['params'] == {'h': 1500}
