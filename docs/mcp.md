# Agent access (MCP)

An agent such as Claude Code or Codex can drive ProxiMate alongside you. The
container serves an [MCP](https://modelcontextprotocol.io) endpoint from the same
process as the GUI, so the agent and your browser share one session.

## Connect an agent

1. Start the container with port 3839 published:

   ```
   docker run -p 3838:3838 -p 127.0.0.1:3839:3839 \
     --mount type=bind,source=<output_dir>,target=/Outputs \
     plutzer/proximate
   ```

   **Keep the `127.0.0.1:` prefix.** The endpoint has no authentication.

2. Register the endpoint once. The GUI sidebar shows the same commands.

   ```
   claude mcp add --transport http proximate http://localhost:3839/mcp   # Claude Code
   codex mcp add proximate --url http://localhost:3839/mcp               # Codex CLI
   ```

   Inside this repository, `.mcp.json` already registers it for Claude Code.
   Codex also accepts an entry in `~/.codex/config.toml`:

   ```toml
   [mcp_servers.proximate]
   url = "http://localhost:3839/mcp"
   ```

   Any other MCP client that speaks streamable HTTP connects to the same URL.

3. Check the endpoint, then ask the agent for help with a dataset:

   ```
   curl localhost:3839/api/health
   ```

   The reply lists the datasets, running jobs and drawn network without touching
   Cytoscape.

## What the agent sees

Five tools:

| Tool | Does |
| --- | --- |
| `search_tools` | Lists the operations |
| `get_tool_details` | Describes one operation and its arguments |
| `call_tool` | Runs one operation |
| `get_gui_documentation` | Explains the GUI tab by tab, so the agent can guide you |
| `view_network` | Returns a picture of the drawn Cytoscape network, fitted to show all of it; no file is written |

Every operation takes its thresholds and settings as explicit arguments. None
reads the GUI's sliders.

## Operations

**read** — no effect on the GUI.

`list_datasets`, `get_dataset_info`, `get_run_info`, `tail_log`, `server_status`,
`cytoscape_read_selection`, `cytoscape_list_nodes`, `cytoscape_get_positions`,
`cytoscape_list_layouts`, `cytoscape_list_style_presets`

**sandbox** — computes and returns a result. Nothing is written under the dataset
and the tabs' settings stay as you left them.

`get_scores`, `threshold_metrics`, `feature_analysis`, `get_prey_annotations`,
`compare_networks`, `cytoscape_estimate`

**dataset** — the Datasets table and the sidebar's active-dataset select update
within a second.

`parse_dataset`, `score_dataset`, `load_session`, `delete_dataset`, `rename_dataset`

- A dataset being scored refuses a second job, from either side.
- `parse_dataset` reads input files by their paths inside the container, so
  mount the folder that holds them.
- `load_session` and `delete_dataset` are destructive and need `confirm`.

**cytoscape** — the drawn network changes while you watch. The Cytoscape tab's
activity table lists each operation with actor `mcp`.

`cytoscape_send`, `cytoscape_apply_thresholds`, `cytoscape_restyle`,
`cytoscape_select_nodes`, `cytoscape_select_related`, `cytoscape_clear_selection`,
`cytoscape_set_edge_visibility`, `cytoscape_move_nodes`,
`cytoscape_cluster_selection`, `cytoscape_save_layout`, `cytoscape_load_layout`,
`cytoscape_export_image`

Sends, exports and layout saves from MCP are recorded in the dataset's `run.json`
with their full arguments.

## Ports

| Variable | Sets |
| --- | --- |
| `PROXIMATE_MCP_PORT`, `PROXIMATE_MCP_HOST` | where the MCP endpoint listens (default port 3839) |
| `PROXIMATE_GUI_PORT`, `PROXIMATE_GUI_HOST` | where the GUI listens (default port 3838) |
