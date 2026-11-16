# Cytoscape setup

The Cytoscape tab draws networks in a Cytoscape desktop running on the same
computer as ProxiMate. Start Cytoscape first; the tab's *Status* panel says
whether it answers.

## Connect the container to Cytoscape

Cytoscape stays on your computer. The container reaches it through CyREST on
port 1234.

| Docker setup | What to do |
| --- | --- |
| Docker Desktop (Windows, macOS) | Nothing. |
| Native Linux engine | Add `--add-host host.docker.internal:host-gateway` to `docker run`. |
| Linux with `--network host` | Set `PROXIMATE_CYTOSCAPE_URL=http://127.0.0.1:1234/v1`. |

| Variable | Default | Effect |
| --- | --- | --- |
| `PROXIMATE_CYTOSCAPE_URL` | `http://host.docker.internal:1234/v1` in a container, `http://127.0.0.1:1234/v1` otherwise | Where the tab reaches CyREST |
| `PROXIMATE_CYTOSCAPE_TIMEOUT` | `120` | Seconds one CyREST request may take before the tab reports an error |

## Where files go

| File | Location |
| --- | --- |
| Exported PNGs | `<output_dir>/<dataset>/cytoscape/` |
| Saved layouts | `<output_dir>/<dataset>/cytoscape/layouts/<name>.csv` |
| py4cytoscape's request log | `$PROXIMATE_LOG_DIR/py4cytoscape/` |

A layout file has three columns, `id,x,y`. `id` is the node name Cytoscape
shows: the UniProt accession for a prey, the bait name for a bait. Layouts
travel with the session zip.

## Using the tab

- [Tutorial 5: Networks and export](tutorials/5-networks-and-export.md) walks
  through sending, styling, selecting and clustering a network.
- The [GUI reference](../GUI/gui_docs.md#cytoscape) describes every control:
  style presets, CORUM complex edges, selection tools, clustering and layouts.
