"""The four MCP tools: search_tools, get_tool_details, call_tool, get_gui_documentation.

The agent's context holds only these four; every operation lives behind ``call_tool``
in ``mcp_registry`` (populated by importing ``mcp_ops``).  ``server.py`` registers
these functions with FastMCP; they are plain functions here so tests drive them
without a transport.
"""

import os
import re

import help_text
import log_config
import mcp_ops  # noqa: F401  (registers the operations)
import mcp_registry as registry
from log_config import get_logger

logger = get_logger(__name__)

GUI_DIR = os.path.dirname(os.path.abspath(__file__))
DOCS_PATH = os.path.join(GUI_DIR, 'gui_docs.md')
README_PATH = os.path.join(os.path.dirname(GUI_DIR), 'README.md')

# README sections that complement a tab's page.
README_SECTIONS = {'Cytoscape': ('Cytoscape',), 'Network Scoring': ('Logs and provenance',
                                                                     'Excluding Human Cell Map evidence')}


def search_tools(query: str = '', mode: str = None) -> list:
    """Find operations for call_tool by a word in their name, summary or tags.

    An empty query lists every operation.  mode narrows to read, sandbox (computes and
    returns without touching the GUI or the dataset), dataset (creates or scores a
    dataset; the GUI updates) or cytoscape (acts on the drawn network).
    """
    return registry.search(query, mode)


def get_tool_details(name: str) -> dict:
    """Full description of one operation: what it does, its effect on the GUI and on
    disk, and its parameter schema with types and defaults."""
    return registry.details(name)


def call_tool(name: str, arguments: dict = None) -> dict:
    """Run one operation with a JSON object of arguments.

    Returns {ok: true, result, run_id} or {ok: false, error, error_type}.  Arguments
    are validated first: unknown or missing keys fail before anything runs.  Every
    threshold argument is an explicit {SaintScore, BFDR, WD, WDFDR} object.
    """
    run_id = log_config.new_run_id()
    try:
        with log_config.run_context(run_id):
            result = registry.call(name, arguments or {}, actor='mcp')
    except Exception as e:
        logger.warning("mcp %s failed: %s: %s", name, type(e).__name__, e)
        return {'ok': False, 'error': str(e), 'error_type': type(e).__name__, 'run_id': run_id}
    return {'ok': True, 'result': result, 'run_id': run_id}


def _split_docs(text):
    """{section title: markdown} for each ``## `` heading, plus the preamble under ''."""
    parts = re.split(r'^## ', text, flags=re.M)
    sections = {'': parts[0].strip()}
    for chunk in parts[1:]:
        title, _, body = chunk.partition('\n')
        sections[title.strip()] = '## ' + title.strip() + '\n' + body.strip()
    return sections


def _readme_section(title):
    with open(README_PATH, encoding='utf-8') as handle:
        text = handle.read()
    match = re.search(r'^#{2,3} ' + re.escape(title) + r'\s*$\n(.*?)(?=^#{2,3} |\Z)', text, re.M | re.S)
    if not match:
        raise KeyError(f"README has no section {title!r}")
    return match.group(1).strip()


def get_gui_documentation(section: str = None) -> dict:
    """What the user sees in the ProxiMate GUI, tab by tab, to help them use it.

    With no section: the table of contents with a paragraph per tab.  With a tab name
    (Network Scoring, Data Thresholding, Protein Feature Analysis, Network Comparison,
    Cytoscape, Downloads): that tab's walkthrough, the tooltip of every control on it
    keyed by control, and any README section that explains it further.
    """
    with open(DOCS_PATH, encoding='utf-8') as handle:
        docs = _split_docs(handle.read())
    tabs = list(help_text.SECTIONS)
    missing = [t for t in tabs if t not in docs]
    if missing:
        raise RuntimeError(f"gui_docs.md lacks sections for {missing}")
    if section is None:
        overview = {tab: docs[tab].split('\n', 1)[1].strip().split('\n\n')[0] for tab in tabs}
        return {'preamble': docs[''], 'sections': tabs, 'overview': overview}
    if section not in help_text.SECTIONS:
        raise KeyError(f"no tab named {section!r}; tabs are {tabs}")
    text = docs[section]
    for title in README_SECTIONS.get(section, ()):
        text += '\n\n### From the README: ' + title + '\n\n' + _readme_section(title)
    controls = {key: help_text.TOOLTIPS[key] for key in help_text.SECTIONS[section]}
    return {'section': section, 'text': text, 'controls': controls}
