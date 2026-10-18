# ProxiMate GUI user-experience overhaul — specification

This document is the hand-off for a series of implementation sessions. Each numbered
chunk below is self-contained, leaves the app working, and is committed by the
repository owner after review. A session implements exactly one chunk, verifies it
as described under *Verification*, and stops.

Decisions recorded here were made with the owner during a review of the running GUI
(Docker image built from this tree, `examples/TIP49_spc` parsed as SAINT input with
quant type "Spectral Counts", scored with 50 WDFDR iterations, then every tab
exercised, including a live Cytoscape 3.10.2). Items marked **Assumption** were not
asked and may be overturned by the owner at the start of a session.

## Ground rules for every chunk

- Shiny for Python (core API, not Express), Bootstrap default theme. No theme
  library. One small stylesheet `GUI/www/proximate.css` (created in chunk 1, extended
  in chunk 2) for spacing, section headings and button groups; no inline `style=`
  strings for anything the stylesheet already covers.
- Fail fast. No silent fallbacks: an error is logged through `notify()` in
  `GUI/app.py` and shown to the user with the detail needed to act.
- Fewer lines beats cleverness. Prefer deleting duplicated handlers over
  abstracting them.
- Comments and docstrings describe the code as it now is, never the change or the
  plan (see `~/.claude/rules/code-comments.md`). Do not mention this spec or chunk
  numbers in the source.
- The GUI and the MCP tools share `GUI/backend.py` and `GUI/dataset_store.py`.
  Any new dataset action (delete, rename) or new sandbox computation (network size
  estimate) is implemented in `backend.py` and registered in `GUI/mcp_ops.py` with
  the right mode. MCP operations keep taking thresholds as explicit arguments; they
  never read the shared GUI state.
- `GUI/help_text.py` (`TOOLTIPS`, `SECTIONS`) and `GUI/gui_docs.md` feed the MCP
  `get_gui_documentation` tool. `Scripts/tests/test_mcp.py` asserts every tooltip
  key is listed in `SECTIONS` and tooltips are at most 250 characters. Every chunk
  that adds, removes or renames a control updates both files in the same change.
- Tests: `pytest` (config in `pytest.ini`, tests in `Scripts/tests/`). New pure
  logic gets unit tests there. Acceptance for each chunk is a Chrome walkthrough
  (below). TDD is not required.
- American spelling in all user-facing text.

## Verification (every chunk)

```bash
docker build -t plutzer/proximate .
docker rm -f proximate_ux; docker run -d --name proximate_ux -p 3840:3838 -p 127.0.0.1:3841:3839 plutzer/proximate
curl -s http://127.0.0.1:3841/api/health
```

Port 3838 may be taken by the owner's `proximate_dev` container; use 3840. Open
`http://localhost:3840` in Chrome (Claude in Chrome extension), parse
`examples/TIP49_spc/{bait,prey,interaction}.txt` as format SAINT, quant type
"Spectral Counts", dataset name `TIP49`, score with WDFDR iterations 50, then walk
every tab the chunk touches and look at screenshots of the rendered result. The
`examples/HCM_LFQ` set lacks its interaction file and cannot be parsed. Cytoscape
chunks need Cytoscape desktop running on the host (the container reaches it at
`host.docker.internal:1234`). Run `pytest` before finishing.

## Current layout (for orientation)

`GUI/app.py` (about 2740 lines) holds `app_ui` (six `ui.nav_panel`s inside
`ui.page_navbar` with a left `ui.sidebar`) and one `server()` function. Tabs:
Network Scoring, Data Thresholding, Protein Feature Analysis, Network Comparison,
Cytoscape, Downloads. Each tab has its own "Select Dataset" `input_select` and
(except Network Scoring) its own four threshold sliders with three preset buttons.
Datasets are rows in `dataset_store` (columns from
`session_archive.DATASET_COLUMNS`: Dataset Name, Input Type, Quant Type,
Experiments, Controls, Scored, Imputation, WDFDR iterations); the GUI polls the
store version once a second. Long jobs run inline under `ui.Progress`, blocking the
session.

---

## Chunk 1 — Defects

**Status: completed.**

All observed on the current build. Fix each; no redesign in this chunk.

1. **Stale imputation choices.** `update_imputation_choices` runs on
   `input.score_dataset` change, before the one-second `datasets()` poll has the new
   row, finds `quant.empty` and returns, so a freshly parsed spectral-count dataset
   still offers the two AFT options. Fix: read the row from `dataset_store.row()`
   directly (not the polled frame) and also re-run when `datasets()` changes.
   Acceptance: parse TIP49 and the radio shows only "Default" without reloading.
2. **Preset highlight.** Every tab's "Stringent" button carries `btn-outline-primary`
   while the others carry `btn-outline-secondary`, which reads as "active" though the
   sliders hold other values. Fix now by giving all preset buttons the same class;
   chunk 2 replaces them with a computed active state.
3. **Scroll offset across tabs.** Switching nav tabs keeps the window scroll, so a
   tab opens mid-page. Scroll to top on tab change (a few lines of JS on
   `shown.bs.tab`, included via `ui.tags.script` in the page header).
4. **Stuck tooltip.** The "Run Feature Analysis" tooltip stays open after the click.
   Same for other `ui.tooltip`-wrapped action buttons. Hide tooltips on click (JS,
   same script as item 3).
5. **Raw exception in Cytoscape status.** `cy_status` prints the `requests`
   exception text (HTTPConnectionPool ... Errno 111). Show one sentence: "Cytoscape is
   not reachable at {url}. Start Cytoscape on this machine or set
   PROXIMATE_CYTOSCAPE_URL." and put the exception text in a collapsed
   `ui.accordion_panel("Details")`. The health dict already carries `error`.
6. **Bait B defaults to Bait A.** In Network Comparison, when the bait list loads,
   select the second bait for B when there are two or more.
7. **Scatter axis label.** `saint_scatter_plot` in `GUI/QC_plots.py` labels the x
   axis "Fold Change (log2)" but plots `FoldChange` as stored (values to 500). Label
   it "Fold change" and switch the axis to log scale (`xaxis_type="log"`) so the
   spread near 1 is readable; same in `saint_scatter_matplotlib`
   (`GUI/plot_exports.py`).
8. **"Upload complete" clipped.** The Shiny file-input progress bar text is cut off
   under each input. CSS: give `.progress` under `.shiny-input-container` enough
   height, or hide the text and keep the bar.
   *Done by creating `GUI/www/proximate.css` (loaded with `ui.include_css` from the
   page header) rather than an inline style; chunk 2 extends that file.*
9. **Dataset name not cleared.** After a successful parse clear `dataset_name` and
   the ED/bait grids (Shiny cannot reset a file input). Acceptance: a second parse
   without renaming does not collide with the first.
10. **Empty "Include" label.** In Downloads, the `dl_groups` checkbox group shows its
    label with no boxes for presets without groups. Hide the whole control when the
    preset has no groups (render it as `output_ui`).
11. **Volcano render freeze.** `create_volcano_plot` in `GUI/network_comparison.py`
    froze the renderer for about 30 s on TIP49 (three panels, many traces). Use
    `go.Scattergl` for the point traces and one trace per category rather than per
    point; keep hover text. Acceptance: ACTR5 vs RUVBL1 renders in under 3 s.
12. **Feature analysis feedback.** No success notification was visible after the
    run. Make `notify()` fire on completion with the row count written to
    `Feature_enrichment.csv`, and show an error notification on failure.

Tests: unit test for the imputation-choices helper if it is extracted; unit test
that the volcano builder emits at most one trace per category per panel.

---

## Chunk 2 — Sidebar: active dataset and shared thresholds

### Active dataset

- The sidebar gains, at the top, `input_select("active_dataset")` listing every
  dataset in the store (scored or not), with a small badge "scored" / "parsed" next
  to it (rendered as `output_ui`). All tabs read this one value. Remove the per-tab
  dataset selects (`qc_dataset`, `feature_dataset`, `cy_dataset`,
  `download_dataset`, `score_dataset`, `comp_dataset_a`, `comp_dataset_b`).
- Network Comparison compares two baits of the active dataset. Cross-dataset
  comparison is dropped (the volcano already required one dataset). Say so in
  `gui_docs.md`.
- Selecting a row in the session table on the first tab (DataGrid
  `selection_mode="row"`) sets the active dataset; parsing a dataset sets it too.
- When the active dataset is not scored, tabs that need scores show their empty-state
  message ("'{name}' is parsed but not scored. Score it on the Data & Scoring tab.")
  instead of controls. When no dataset exists the existing empty states remain.

### Shared thresholds

- Below the active dataset, four `input_slider`s (SAINT ≥, BFDR ≤, WD ≥, WDFDR ≤;
  same ranges and steps as today) and three preset buttons. One set for the whole
  app; remove the sliders and preset handlers from Data Thresholding, Protein Feature
  Analysis, Cytoscape, Downloads, and the duplicated preset effects
  (`apply_qc_*`, `apply_fa_*`, `apply_dl_*`, `_make_preset`).
- A single reactive `thresholds()` calc returns the dict every consumer uses.
- **Presets** (replaces Stringent / Moderate / Relaxed everywhere, including
  `help_text.TOOLTIPS` and `gui_docs.md`):

  | Preset       | SAINT ≥ | BFDR ≤ | WD ≥ | WDFDR ≤ |
  | ------------ | ------- | ------ | ---- | ------- |
  | SAINT only   | 0.9     | 0.05   | 0.0  | 1.0     |
  | WDFDR only   | 0.0     | 1.0    | 0.0  | 0.05    |
  | Both         | 0.9     | 0.05   | 0.0  | 0.05    |

  **Assumption:** WD is left at 0 (disabled) in all three, since the paper's sets
  are defined on SAINT/BFDR and WDFDR alone. Define the table once as a module
  constant `PRESETS` in `GUI/thresholds.py` (new, tiny) so the tooltips and the
  handlers read the same numbers.
- **Startup default:** the sliders start at the "SAINT only" values.
  **Assumption:** chosen because a dataset scored with 0 WDFDR iterations has no
  WDFDR values and a WDFDR cut would empty it. The existing `wdfdr_warning` moves to
  the sidebar under the sliders and shows whenever the active dataset was scored
  with 0 iterations and WDFDR < 1.
- **Active state is computed:** the preset whose four values equal the current
  sliders gets `btn-primary`; otherwise none does. Render the button row as
  `output_ui` so the class follows the values.
- Network Comparison keeps its A and B threshold sets but they initialize from the
  shared set when the tab's bait lists load, plus a "Use shared thresholds" button
  on each card that copies the current shared values.
- Export tab: an `input_switch("dl_apply_thresholds")`, **off by default**, whose
  label reads "Apply thresholds (SAINT ≥ 0.9, BFDR ≤ 0.05, WD ≥ 0, WDFDR ≤ 1.0)"
  with live values. Off exports unfiltered rows, as today's "No Filter". Remove the
  Export tab's own sliders and the four preset buttons there.
- Feature Enrichment uses the shared set; its 0.9 / 1.0 / 0 / 1.0 defaults go away.

### Sidebar structure

Order: active dataset (with badge), thresholds + presets, WDFDR warning, then two
collapsed `ui.accordion` panels: "Agent access" (the two MCP connect commands) and
"About and feedback" (bug report and repository links). The version label becomes
a one-line muted footer. Sidebar width about 300 px; sliders must not wrap their
value bubbles.

### Tests

`Scripts/tests/test_thresholds.py`: `PRESETS` values; `active_preset(values)` returns
the matching key or `None`; each preset dict passes `backend.validate_thresholds`.

### Acceptance walkthrough

Change a slider in the sidebar and see the Quality Control metrics, the Feature
Enrichment foreground, the Cytoscape send and the Export preview (with the switch
on) all follow it. Pick a dataset by clicking its row; every tab follows.

---

## Chunk 3 — Data & Scoring tab, dataset actions, progress modal

### Layout

Top to bottom:

1. **Session table card** (full width). Columns: Dataset, Input, Quant, Experiments,
   Controls, State, Scored at, Organism, Imputation, WDFDR iterations. *State* is one
   of `parsed`, `scored`, `enriched` (enriched = `Feature_enrichment.csv` exists);
   *Scored at* and *Organism* come from `run.json` via `backend.run_info` /
   `provenance.annotation_settings`. Store the new fields in the dataset row
   (`session_archive.DATASET_COLUMNS` grows; `load_session` must accept archives
   without them). Row selection sets the active dataset (chunk 2).
   Buttons beside the table: *Delete selected*, *Rename selected*, *Clear All
   Datasets*, *Download Session Zip*; the session-zip upload + *Load Session* stay
   in a smaller card to the right.
2. **Step 1 — Parse** card and **Step 2 — Score** card side by side, each headed
   with its step number and a one-line purpose. The score card's dataset is the
   active dataset (no select); its header shows "Scoring: {name}".

### Dataset actions (backend + MCP)

- `backend.delete_dataset(name, actor)`: takes the dataset's job lock
  (`BusyError` if running), removes the directory and the store row, logs it. If the
  deleted dataset is drawn in Cytoscape, leave the drawing alone but clear the
  controller's dataset reference so status reads "drawn dataset no longer exists".
- `backend.rename_dataset(old, new, actor)`: validates `new` with
  `parse.validate_name`, refuses an existing name, takes the lock, renames the
  directory, rewrites the Dataset Name in the store row and in `run.json`, and
  renames the Cytoscape network title if that dataset is drawn.
- Register both in `mcp_ops.py` with mode `dataset`; delete requires
  `confirm=True` like `load_session`. GUI: Delete opens a `ui.modal` confirm naming
  the dataset; Rename opens a modal with a text input.

### Parse step improvements

- **Clearer file requirements:** under the format select, a short muted paragraph
  per format stating the expected file and which columns are read (text lives in
  `help_text.py` as `FORMAT_NOTES[format]` so `get_gui_documentation` can serve
  it). A *Download ED template* `download_button` writes a CSV with the header
  `Experiment Name,Type,Bait,Replicate,Bait ID,Group` and two example rows.
- **Post-parse summary:** after a successful parse, show a summary strip at the
  bottom of the parse card: experiments, controls, preys, interactions written, and
  a *Score now* button that focuses the Score card. Values come from `run_parse`'s
  return (extend it to return counts; `provenance` already records row counts).

### Progress modal for long jobs

Applies to parse, score, feature enrichment and Cytoscape send.

- Jobs run in a worker thread (`threading.Thread`) started from the handler; the
  handler returns immediately. The job lock in `backend._job` already prevents
  double runs. The thread writes its stage, result or exception into a lock-guarded
  dict that the session reads with `reactive.poll` once a second (the same pattern
  that surfaces MCP jobs today).
- A `ui.modal` (not dismissible while running) shows: the job title ("Scoring
  TIP49"), a vertical stage list with the current stage highlighted (taken from the
  `progress` callback messages already emitted by `backend.run_parse` /
  `run_score`), elapsed time, and the last 15 lines of the dataset's
  `proximate.log` (`backend.log_tail`, monospace, auto-scrolled).
- On success the modal closes and `notify()` reports it. On failure the modal stays
  open with the error text and the log tail, and a Close button.
- The MCP banner for agent jobs stays as is.

### Tests

`test_backend.py`: delete and rename (happy path, busy refusal, name collision,
invalid name). `test_session_archive.py`: an archive without the new columns loads.
`test_mcp.py`: the two new operations appear in search.

### Acceptance walkthrough

Parse TIP49 and watch the modal stages; see the summary strip; click *Score now*;
watch the scoring modal with the log tail; rename to `TIP49b`; delete it; confirm
Cytoscape status copes if it was drawn.

---

## Chunk 4 — Quality Control tab: results table and plots

### Results table

Below the metrics row, a card "Passing interactions" with a `render.DataGrid`
(filters on, sortable, height about 400 px) of `backend.passing_scores(active,
thresholds, baits)` restricted to the QC bait selector (All or one bait). Columns,
in order: Bait, Prey gene, Prey ID, SAINT, BFDR, WD, WDFDR, Fold change, Abundance
(AvgSpec or AvgIntensity, whichever exists), In BioGRID, Multivalidated. Numbers
rounded to 3 decimals (format in the frame, not with CSS). A *Download CSV* button
writes the same rows. Row count shown above the grid.

### Experiment PCA legend

In `QC_plots.pca_plot`: controls become one grey series "Control" regardless of
bait name; test baits keep a color each; symbol by Type stays. Legend stays
horizontal below the plot but is limited to one entry per bait (not per
experiment). On datasets with more than 20 baits, hide the legend by default and add
a "Show legend" checkbox in the card header.

### Empty states

Each plot card renders a short message and the action that fills it instead of a
blank area: scatter ("Select one bait above to see its SAINT score against fold
change"), feature heatmap ("Run Feature Enrichment to see the top features"),
volcano and Venn ("Choose two baits and click Compare"). Implement as `output_ui`
that returns the message or the plot widget; no fixed 650 px height when empty.

### Heatmap readability

`Ann_Enrichment.plot_results`: figure height = 0.28 in × number of rows + 2 in
(minimum 5 in); row label font 9; labels truncated at 60 characters; the
downloadable enrichment table carries the full term and its GO id. The Shiny
`output_plot` height follows the figure height (`render.plot` with `height`
computed from `num_features`).

### Scatter

Done in chunk 1 (log axis). Here add a note under the plot: "Dashed line: SAINT
threshold from the sidebar."

### Tests

`test_qc_plots.py`: pca legend has one entry per test bait plus one "Control";
`test_backend.py`: passing_scores columns and the rounding helper.

---

## Chunk 5 — Cytoscape tab

### Activity panel

Replace the `output_text_verbatim` with a `render.DataGrid` (or a simple HTML table
via `output_ui`) with columns Time, Actor, Operation, Detail, newest first, text
wrapping, height about 300 px. Source stays `cytoscape_ctl.snapshot()['log']`.

### Pre-send estimate and warning

- `backend.estimate_network(name, thresholds, baits, prey_prey, biogrid_scope,
  corum, corum_min_members, corum_min_fraction)` returns `{'baits', 'preys',
  'nodes', 'proximity_edges', 'literature_edges', 'complex_edges'}` by calling
  `cytoscape_net.build` and counting. Register as a `sandbox` MCP operation
  `cytoscape_estimate`.
- The Network card shows "Will draw: N baits, X nodes, Y edges" under the send
  options, recomputed with a 500 ms debounce when thresholds, baits or edge options
  change. On error (no passing interactions) it shows the message in muted red.
- Sends above `SEND_WARN_EDGES = 5000` edges open a confirm modal stating the
  counts and suggesting "Multivalidated only" BioGRID scope or fewer baits.
  **Assumption:** 5000 is the limit; a module constant in `cytoscape_net.py`.

### Controls grouped into collapsibles

The middle card becomes a `ui.accordion` with panels: Status (open), Selection
(open), Select by relation (closed), Cluster the selection (closed), Saved layouts
(closed), Export (closed). The left Network card groups its options under headings:
Network (baits, thresholds summary read from the sidebar), Style, Edges, Complexes,
then the three action buttons. The thresholds themselves are in the sidebar
(chunk 2); the Network card shows them read-only in one line.

### Selection table polish

Round SaintScore, BFDR, WD, WDFDR to 3 decimals, FoldChange to 2; add a "Target
gene" column from the node table's symbol; grid height 400 px; a *Download
selection CSV* button; the row count in the panel header.

### Tests

`test_cytoscape.py`: estimate counts on a small synthetic frame; rounding helper.

### Acceptance walkthrough

With Cytoscape running: estimate updates as thresholds change; sending all baits
with "All pairs" triggers the warning; the send draws; select by relation with seed
RUVBL1 then Read Selection fills the polished table; the activity table lists the
operations.

---

## Chunk 6 — Names, stylesheet, documentation pass

- Rename tabs: Data & Scoring, Quality Control, Feature Enrichment, Compare Baits,
  Cytoscape, Export. Update every reference: `app.py`, `help_text.py` SECTIONS keys
  and tooltip texts that name a tab, `gui_docs.md` headings and prose, `mcp_ops.py`
  / `mcp_tools.py` summaries, `network_comparison.py` strings, `README.md`,
  `CLAUDE.md` GUI section, `CHANGELOG.md` entry. `test_mcp.py` must still pass.
- Step numbers appear on the first two cards only, not in tab titles.
- Stylesheet pass: consistent card padding, section heading style, button-group
  spacing, muted help text class; remove inline `style=` strings that the
  stylesheet now covers.
- Documentation pass: re-read `gui_docs.md` end to end against the running GUI and
  make it describe the current app (shared thresholds, active dataset, presets by
  name, progress modal, dataset actions, estimate, renamed tabs).
- Final full walkthrough of all six tabs with TIP49.

---

## Out of scope (decided)

Live validation of the experimental design table while editing; auto-generating the
design from the quant file; a theme library or brand styling; Playwright tests;
cross-dataset bait comparison.
