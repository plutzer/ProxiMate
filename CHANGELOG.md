# Changelog

All notable, user-visible changes to ProxiMate are listed here, newest first. The
format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and versions
follow [Semantic Versioning](https://semver.org/): a patch release fixes bugs, a
minor release adds features, a major release changes results or file formats in a
way that is not backwards compatible.

Each pull request adds one line under **Unreleased**. Cutting a release renames that
section to the version and date and starts a fresh Unreleased section above it.

## [Unreleased]

## [0.2.1] - 2026-09-24

### Changed
- Example datasets live in `examples/` (TIP49 spectral counts from the SAINTexpress
  distribution and the Human Cell Map LFQ set), each with a README; the SAINTexpress
  sources live under `saint/` as `upstream/` and `patches/`.
- `THIRD_PARTY_LICENSES.md` lists the licenses of the bundled SAINTexpress, Boost,
  NLopt, GOGO and CORUM.
- Hovering a point in the prey PCA shows the gene name in bold above the accession.
- The Protein Feature Analysis parameters card takes a quarter of the row, leaving the
  heatmap three quarters.
- Scoring a spectral-count dataset with AFT imputation is refused up front with a
  message naming the reason; the GUI offers only Default for such datasets.
- Two MCP operations: `get_prey_annotations` summarizes each prey of a scored dataset
  (passing baits at given thresholds, BioGRID partners, best scores, localization, GO
  CC, complex) for all preys or a list of accessions or symbols, and
  `cytoscape_list_nodes` lists the drawn nodes with accession, symbol and role.
- A fifth MCP tool, `view_network`, returns a picture of the drawn Cytoscape network
  that the agent sees directly, and `cytoscape_clear_selection` deselects everything.
- MCP operations `get_scores` (the interactions passing a threshold set), `get_run_info` and
  `tail_log` (a dataset's run manifest and log); `score_dataset` defaults to the
  Scoring card's settings; `compare_networks` can omit the volcano table;
  `cytoscape_send` is refused while a network is drawn unless `replace` is true, and
  `server_status` with `probe` names the network in Cytoscape's window.
- The MCP operation list is pruned to what an agent can use: `cytoscape_status` folds
  into `server_status`; loners and satellites become `singletons` (with
  `include_seed`) and `satellites` relations of `cytoscape_select_related` with an
  explicit seed; `cytoscape_sync_positions`, `cytoscape_unlock` and the base64
  `cytoscape_view_image` are gone (the GUI keeps its buttons, `view_network` its picture).
- Cluster and Repack leaves selected baits where they are and moves only the preys.
- The MCP activity entries name the arguments an operation was called with and what
  it returned, such as the relation and the number of nodes selected.
- Annotation builds its HPA, CORUM and GO lookups once and evaluates them per distinct
  prey or bait-prey pair, and GOGO scores only the cellular-component ontology; a
  212k-row dataset annotates in under three minutes instead of fifteen.

### Fixed
- The `Self-Interaction` column in `annotated_scores.csv` compares each row's prey
  against its own bait; it used to reflect a single arbitrary bait.
- The feature enrichment panel no longer shows the last exported PNG figure when it has
  nothing to draw.
- Cytoscape PNG exports release any view lock before fitting, so a view zoomed by a
  selection no longer crops the picture.
- Sending a second dataset to Cytoscape removes the earlier ProxiMate network instead
  of leaving it behind stripped of its style.
- Release builds no longer fail while preprocessing BioGRID; the download step now installs pandas.

### Removed
- The SAINTexpress precompiled binaries, its demo input files and an unused copy of its
  spectral-count sources are no longer tracked; the image builds both binaries from source.

## [0.2.0] - 2026-09-21

### Added
- Cytoscape bait nodes are keyed on bait name, so two constructs of one protein draw as
  separate baits; every MCP dataset operation names its dataset with `dataset`.
- README, CLAUDE.md and the GUI sidebar give the Claude Code and Codex CLI commands
  that connect an agent to the MCP endpoint.
- Cytoscape tab with live link to a Cytoscape desktop: send a thresholded network,
  read the selection back with scores, hide and show edges, re-threshold without
  disturbing a hand layout, export a PNG.
- Cytoscape edge style controls: choose what edge width encodes (abundance, SAINT
  score, WD score, fold change, uniform), weight BioGRID edges by publication count,
  limit BioGRID edges to multivalidated pairs.
- CORUM complex edges in the Cytoscape network, with minimum-subunit and
  minimum-coverage criteria (human datasets).
- Cytoscape selection tools: loners and satellites of one bait, select by relation
  (interactors, singletons, partners, co-complex members), Leiden clustering of the
  selection with re-packing by community.
- MCP server on port 3839, served from the same process as the GUI: four tools
  (`search_tools`, `get_tool_details`, `call_tool`, `get_gui_documentation`) give an
  agent parsing, scoring, threshold metrics, feature analysis, network comparison and
  the Cytoscape controls, with every threshold passed explicitly and every action
  labeled `mcp` in the Cytoscape activity panel. The datasets table is shared across
  browser sessions and updates within a second of a change from either side.
- Pioneer input format.
- Tooltips on every control.
- Tests and logging for the scoring and annotation functions.
- Automated container builds for amd64 and arm64 on every GitHub Release, with the
  annotation database snapshot recorded in the release notes.

### Changed
- PCA and download tabs reworked.
- Network comparison and thresholding improved.

### Removed
- The prebuilt x86-64 SAINTexpress binary that the image copied but never ran.

## [0.1.11] - 2026-06-23

Last version published to Docker Hub by hand, before this changelog existed. Older
versions are listed at https://hub.docker.com/r/plutzer/proximate/tags.

[Unreleased]: https://github.com/plutzer/ProxiMate/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/plutzer/ProxiMate/releases/tag/v0.2.1
[0.2.0]: https://github.com/plutzer/ProxiMate/releases/tag/v0.2.0
[0.1.11]: https://hub.docker.com/r/plutzer/proximate/tags
