"""The ProxiMate server: the Shiny GUI and the MCP tools in one process.

    python3 GUI/server.py

- GUI at                     http://localhost:3838        (PROXIMATE_GUI_PORT)
- MCP (streamable-http) at   http://localhost:3839/mcp   (PROXIMATE_MCP_PORT)
- health at                  http://localhost:3839/api/health

Both surfaces call the same ``backend`` over the same in-process state (the datasets
store and the Cytoscape controller), so a person at the GUI and an agent over MCP see
one session.  The MCP port carries no authentication: publish it on the loopback
interface only.
"""

import functools
import os
import sys
import threading
import time

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'Scripts'))

import log_config  # noqa: E402
from log_config import get_logger  # noqa: E402

logger = get_logger(__name__)

GUI_HOST = os.environ.get('PROXIMATE_GUI_HOST', '0.0.0.0')
GUI_PORT = int(os.environ.get('PROXIMATE_GUI_PORT', '3838'))
MCP_HOST = os.environ.get('PROXIMATE_MCP_HOST', '0.0.0.0')
MCP_PORT = int(os.environ.get('PROXIMATE_MCP_PORT', '3839'))
STARTED = time.time()


def threaded(fn):
    """Run a tool body in a worker thread.

    FastMCP calls a plain function inline on the event loop, so a scoring run or a
    Cytoscape round trip inside it would stall every other request on this port for
    as long as it runs.  The wrapper keeps the signature FastMCP reads for the schema.
    """
    from starlette.concurrency import run_in_threadpool

    @functools.wraps(fn)
    async def run(*args, **kwargs):
        return await run_in_threadpool(functools.partial(fn, *args, **kwargs))
    return run


def build_mcp(host=MCP_HOST, port=MCP_PORT):
    """The FastMCP server with the five tools and the health route registered."""
    from mcp.server.fastmcp import FastMCP
    from starlette.requests import Request
    from starlette.responses import JSONResponse

    import backend
    import cytoscape_ctl
    import dataset_store
    import mcp_registry
    import mcp_tools

    mcp = FastMCP('proximate', host=host, port=port)
    for tool in (mcp_tools.search_tools, mcp_tools.get_tool_details, mcp_tools.call_tool,
                 mcp_tools.get_gui_documentation, mcp_tools.view_network):
        mcp.tool()(threaded(tool))

    @mcp.custom_route('/api/health', methods=['GET'])
    async def health(_request: Request) -> JSONResponse:
        """The process's state without touching Cytoscape."""
        snap = cytoscape_ctl.snapshot()
        return JSONResponse({
            'ok': True, 'uptime_s': int(time.time() - STARTED), 'out_dir': backend.OUT_DIR,
            'datasets': dataset_store.names(), 'datasets_version': dataset_store.version(),
            'jobs': backend.running_jobs(),
            'cytoscape': {'dataset': snap['dataset'], 'version': snap['version'], 'busy': snap['busy']},
            'operations': len(mcp_registry.OPS), 'activity': mcp_registry.activity(5),
        })
    return mcp


def _run_gui():
    import uvicorn

    from app import app
    uvicorn.run(app, host=GUI_HOST, port=GUI_PORT, log_level='warning')


def main():
    log_config.setup_logging()
    import app  # noqa: F401  (configures the backend from PROXIMATE_OUTPUT_DIR and logs startup)
    mcp = build_mcp()
    threading.Thread(target=_run_gui, daemon=True).start()
    logger.info("GUI on http://%s:%d  MCP on http://%s:%d/mcp", GUI_HOST, GUI_PORT, MCP_HOST, MCP_PORT)
    mcp.run(transport='streamable-http')


if __name__ == '__main__':
    main()
