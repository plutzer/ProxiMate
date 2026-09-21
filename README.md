# ProxiMate
All-in-one GUI and scripts for analyzing proximity labelling data.

## How to run the tool
### Running the GUI locally (easiest)
1. Download and install Docker: https://www.docker.com/products/docker-desktop/
2. Run the pre-built container from dockerhub, exposing the 3838 port and mounting
   a directory for the results:
    ```
    docker run -p 3838:3838 \
      --mount type=bind,source=<native_path_to_output_directory>,target=/Outputs \
      plutzer/proximate
    ```
   Without the mount, everything the tool writes — datasets, logs and run manifests —
   lives inside the container and is lost when it is removed.
3. Access the GUI through a web browswer at localhost:3838

### Running the docker container interactively (experienced users)
The backend of the application can be accessed interactively by overriding the command to start the shiny app:
1. Start the docker container interactively. In order to analyze files, you'll need to mount a directory to the container:
    `docker run -it --mount type=bind,source=<native_path_to_data_directory>,target=<working_directory_within_container> plutzer/proximate /bin/bash`
2. Within the docker container, python, R, perl, or SAINT scripts can be run manually.

### Running scripts on a high-performance computing cluster (experienced users)
The docker container contains a shell script to run the entire pipeline. Use the `--format` flag to specify your input type:

**MaxQuant:**
```
docker run --mount type=bind,source=<data_dir>,target=<container_dir> plutzer/proximate \
  /bin/bash /run_pipeline.sh --format maxquant \
  <ED_file> <PG_file> <quant_type> <output_dir> <n_iterations> <imputation>
```

**DIA-NN:**
```
docker run --mount type=bind,source=<data_dir>,target=<container_dir> plutzer/proximate \
  /bin/bash /run_pipeline.sh --format diann \
  <ED_file> <matrix_file> <output_dir> <n_iterations> <imputation>
```

**FragPipe:**
```
docker run --mount type=bind,source=<data_dir>,target=<container_dir> plutzer/proximate \
  /bin/bash /run_pipeline.sh --format fragpipe \
  <ED_file> <FP_file> <quant_type> <output_dir> <n_iterations> <imputation>
```

**SAINT:**
```
docker run --mount type=bind,source=<data_dir>,target=<container_dir> plutzer/proximate \
  /bin/bash /run_pipeline.sh --format saint \
  <bait_file> <prey_file> <interaction_file> <quant_type> <output_dir> <n_iterations> <imputation>
```

**Options:**
- `--organism`: `human` (default), `mouse`, or `yeast` — add before `--format` if needed
- `quant_type`: `Intensity`, `LFQ`, or `Spectral Counts` (DIA-NN and Pioneer always use Intensity; FragPipe supports all three)
- `imputation`: `0` (none), `1` (prey-specific AFT), `2` (refactored AFT), or `3` (one-component AFT)
- `--seed`: CompPASS permutation seed, so WD p-values reproduce between runs

The pipeline stops at the first stage that fails, rather than carrying on with
missing inputs.

## Logs and provenance

Three artifacts record what the tool did.

| File | Scope | Contents |
| --- | --- | --- |
| `<output_dir>/run.json` | one dataset | parameters, input checksums, row counts, which SAINTexpress binary ran, package versions, pass/fail per stage |
| `<output_dir>/proximate.log` | one dataset | timestamped log from parsing, scoring and annotation |
| `<output_dir>/log.txt` | one dataset, CLI only | raw console transcript, including output the logger never sees |
| `$PROXIMATE_LOG_DIR/proximate-server.log` | all datasets | the server's operational log, rotated at 10 MB |

Every line and every manifest entry carries a **run ID** shared by the GUI and the
parse, score and annotate processes it launches, so one run can be traced across
all four. Attach `run.json` to a bug report.

**Environment variables**

| Variable | Default | Effect |
| --- | --- | --- |
| `LOG_LEVEL` | `INFO` | `DEBUG` adds per-prey imputation detail. Affects ProxiMate's own loggers only. |
| `PROXIMATE_LOG_DIR` | `/Outputs` (container), the output directory (CLI) | where the operational log is written |
| `PROXIMATE_OUTPUT_DIR` | `/Outputs` | where the GUI stores datasets |
| `PROXIMATE_RUN_ID` | minted per run | set it to correlate an external job with a ProxiMate run |
| `PROXIMATE_CYTOSCAPE_URL` | `http://host.docker.internal:1234/v1` in a container, `http://127.0.0.1:1234/v1` otherwise | where the Cytoscape tab reaches CyREST |
| `PROXIMATE_CYTOSCAPE_TIMEOUT` | `120` | seconds any one CyREST request may take before the tab reports an error |

## Cytoscape

The Cytoscape tab sends the interactions passing a set of thresholds into a Cytoscape
desktop running on the same machine and keeps a link to it: read the selection back
with its scores, hide or show edges around it, tighten the thresholds without
disturbing a hand layout, and save a PNG under `<output_dir>/<dataset>/cytoscape/`.
Each send replaces the previous ProxiMate network. Cytoscape stays on the host; the
container only talks to CyREST on port 1234, so on Docker Desktop no extra flag is
needed. On a native Linux engine add `--add-host host.docker.internal:host-gateway`
to `docker run`, or run with `--network host` and set `PROXIMATE_CYTOSCAPE_URL` to
`http://127.0.0.1:1234/v1`. py4cytoscape's own request log lands in
`$PROXIMATE_LOG_DIR/py4cytoscape/`.

**Edge style.** Bait–prey edge width encodes abundance by default (log-banded
intensity or spectral counts); the *Edge width* select swaps in SAINT score, WD score,
fold change, or a uniform width. Prey–prey BioGRID edges can be limited to pairs BioGRID
marks multivalidated, and thickened by the number of publications behind them. *Apply
Edge Style* pushes these onto the drawn network as a style update, so a hand layout
survives.

**Complexes.** For human datasets, *CORUM complex edges* draws black edges between drawn
proteins that are subunits of one curated complex the screen recovered: at least *Min
subunits drawn* of its members are in the network and one bait, counted with its preys,
covers at least *Min share by one bait* of the full membership. The criteria apply when
the network is sent.

**Selection tools.** With one bait selected in Cytoscape, *Select Loners* picks it with
the preys whose only visible neighbor it is and *Select Satellites* adds the two-bait
preys that currently sit nearer to it. *Select by relation* names a seed and selects its
interactors (with SAINT, BFDR and abundance cuts), singletons, BioGRID or complex
partners, or CORUM co-complex members, replacing or adding to the selection. *Cluster
and Repack* runs Leiden over the selected nodes, colors them by community and re-packs
each community on its own circle inside the box the selection occupies; only the
selected nodes move, and the community number lands in the node table.

**Building with a version stamp.** The image contains no git repository, so the
version recorded in `run.json` comes from a build argument:

```
docker build --build-arg PROXIMATE_VERSION=$(git rev-parse --short HEAD) \
  -t plutzer/proximate:latest .
```

## Agent access (MCP)

The container also serves an [MCP](https://modelcontextprotocol.io) endpoint at
`http://localhost:3839/mcp` (streamable HTTP) from the same process as the GUI, so an
agent such as Claude Code and a person at the browser share one session. Publish the
port on the loopback interface only; it carries no authentication:

```
docker run -p 3838:3838 -p 127.0.0.1:3839:3839   --mount type=bind,source=<native_path_to_output_directory>,target=/Outputs   plutzer/proximate
```

The repository's `.mcp.json` points Claude Code at that address; from any other
directory, register it once with

```
claude mcp add --transport http proximate http://localhost:3839/mcp
```

(the GUI's sidebar shows the same command). The agent sees four
tools: `search_tools` and `get_tool_details` describe the operations, `call_tool` runs
one, and `get_gui_documentation` explains the GUI tab by tab so the agent can help a
person use it. Every operation takes its thresholds and settings as explicit arguments;
none reads what the GUI's sliders show. Operations are classed by what they touch:

| Mode | Operations | Effect on a person at the GUI |
| --- | --- | --- |
| read | `list_datasets`, `get_dataset_info`, `server_status`, `cytoscape_status`, `cytoscape_read_selection`, `cytoscape_get_positions` | none |
| sandbox | `threshold_metrics`, `feature_analysis`, `compare_networks` | none: results are returned, nothing is written under the dataset, and the tab's own settings stay as the user left them |
| dataset | `parse_dataset`, `score_dataset`, `load_session` | the dataset table and dropdowns update within a second; a dataset being scored by either side refuses a second job; `load_session` is destructive and needs `confirm` |
| cytoscape | `cytoscape_send`, `cytoscape_apply_thresholds`, `cytoscape_restyle`, `cytoscape_select_*`, `cytoscape_set_edge_visibility`, `cytoscape_move_nodes`, `cytoscape_cluster_selection`, `cytoscape_sync_positions`, `cytoscape_export_image`, `cytoscape_unlock` | the drawn network changes under the mouse; the Cytoscape tab's status shows the thresholds it was drawn at and the activity panel lists each operation as `[mcp]` |

Sends and exports from MCP are recorded in the dataset's `run.json` with their full
argument set. `curl localhost:3839/api/health` reports the datasets, running jobs and
drawn network without touching Cytoscape. `PROXIMATE_MCP_PORT`, `PROXIMATE_MCP_HOST`,
`PROXIMATE_GUI_PORT` and `PROXIMATE_GUI_HOST` move the ports.

## Annotations and Databases:
Every release ships the BioGRID, UniProt and Human Protein Atlas snapshot downloaded
on the day it was built; CORUM is tracked in the repository. The download dates are
in `/Datasets/build_info.txt` inside the image, in each release's notes on GitHub, and
in every `run.json`.

### Excluding Human Cell Map evidence

The Human Cell Map (Go et al. 2021, PubMed 34079125) is a large BioID screen deposited in
BioGRID. Scoring a proximity-labeling experiment against it counts hits as "known" on the
strength of the same kind of assay, so the container also builds
`/Datasets/human/biogrid_summary_no_hcm.csv`: the human summary with every evidence line
from that paper removed. Interactions with other evidence are kept, with the HCM entries
dropped from their evidence lists. Select it with the "Exclude Human Cell Map evidence"
checkbox in the scoring panel (human only) or `--exclude-hcm` on `run_pipeline.sh`. The
choice is recorded in `run.json`, and the QC tab's known-interaction and network-degree
metrics use the same file the run was annotated against.

### Updating databases manually (experienced users):
Databases can be downloaded automatically by running `python3 Scripts/setup_datasets.py --output-dir Datasets`. Use `--skip` to exclude specific databases (e.g., `--skip corum`). See `python3 Scripts/setup_datasets.py --help` for all options.

Alternatively, you can manually assemble the following files in a `/Datasets` subdirectory inside the ProxiMate parent directory and re-build the docker container. File names will need to match or be changed in the Dockerfile before building.


`BIOGRID-ALL.tab3.txt` - downloaded from [BioGRID](https://downloads.thebiogrid.org/BioGRID)

`BIOGRID-MV-Physical.tab3.txt` - downloaded from [BioGRID](https://downloads.thebiogrid.org/BioGRID)

`uniprot_anns.tsv` - tsv-formatted annotations for human proteins from [UniProt](https://www.uniprot.org/uniprotkb?query=%28proteome%3AUP000005640%29&facets=reviewed%3Atrue)

`subcellular_location.tsv` - downloaded from [Human Protein Atlas](https://www.proteinatlas.org/humanproteome/subcellular/data#locations)

`corum_humanComplexes.txt` - downloaded from [CORUM](https://mips.helmholtz-muenchen.de/corum/download): Human Complexes.




## Releases

Images are published for `linux/amd64` and `linux/arm64` (Apple Silicon runs it
natively) to two registries; either name pulls the same image:

```
docker pull plutzer/proximate:v0.2.0          # Docker Hub
docker pull ghcr.io/plutzer/proximate:v0.2.0  # GitHub Container Registry
```

`latest` follows the newest release. Versions are `vMAJOR.MINOR.PATCH`; what changed
in each is in `CHANGELOG.md` and on the
[Releases page](https://github.com/plutzer/ProxiMate/releases), where the notes
also list the image digest and the annotation database dates. Cite a version with
the DOI Zenodo mints for each release, or with `CITATION.cff`.

Publishing a release on GitHub runs `.github/workflows/release.yml`: it downloads
the databases once, builds both architectures natively, joins them under the
release tag on both registries, and appends the digest and `build_info.txt` to the
release notes. `release.sh` performs the same build from a local checkout of the
tag when GitHub Actions is not an option:

```
./release.sh v0.2.0                  # both architectures, pushed to Docker Hub
./release.sh v0.2.0 --no-push        # amd64 only, loaded into local Docker
```

## Common Errors

### SAINT
Permission Denied:
    Can occur when python tries to run a subprocess on a directory rather than a file. Is the SAINT directory correct?

Invalid delimiter:
    SAINT throws this error when the file path is incorrect
