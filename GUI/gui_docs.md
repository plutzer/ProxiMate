# The ProxiMate GUI

ProxiMate scores proximity-labeling proteomics screens. A run moves left to right
through the navigation bar: parse the raw quantification into a **dataset**, score it,
then explore the scored network in the remaining tabs. Every dataset is a directory
under the output directory holding the SAINT inputs, `merged.csv` (SAINT and CompPASS
scores), `annotated_scores.csv` (with BioGRID, UniProt, HPA and CORUM annotation),
`proximate.log` and `run.json` (the run manifest). The "Datasets in this Session"
table on the first tab lists every dataset and whether it is scored. The sidebar names
the **active dataset**, which every tab works on, and holds the four score thresholds
every tab reads.

## Sidebar

The *Active dataset* select lists every dataset in the session with a badge saying
whether it is *parsed* or *scored*. Clicking a row of the session table on the Network
Scoring tab picks that dataset, and parsing a dataset makes it active. Tabs that need
scores show a message instead of their controls while the active dataset is only parsed.

Below it, four sliders set the score thresholds every tab reads: SAINT score (higher is
better), BFDR (lower is stricter), CompPASS WD score (higher is better; 0 disables), and
WDFDR (lower is stricter; 1 disables). An interaction passes only when it clears all
four. Three preset buttons set all four at once: *SAINT only* (SAINT ≥ 0.9, BFDR ≤ 0.05,
WDFDR off), *WDFDR only* (WDFDR ≤ 0.05, SAINT off) and *Both* (SAINT ≥ 0.9, BFDR ≤ 0.05,
WDFDR ≤ 0.05); WD stays at 0 in all three. The button whose values the sliders hold is
highlighted. The sliders start at *SAINT only*. A warning appears when the active
dataset was scored with 0 WDFDR iterations and the WDFDR slider is below 1, because such
a dataset has no WDFDR values and nothing would pass.

Two collapsed panels hold the commands that connect an agent (Claude Code or Codex CLI)
to this server's MCP endpoint and the bug-report and repository links. The version
label sits at the bottom.

## Network Scoring

Three cards. The **Datasets in this Session** table lists every dataset with its input
and quant type, experiment and control counts, its *State* (`parsed`, `scored`, or
`enriched` once Feature Enrichment has been run), when it was scored and against which
organism, and its imputation and WDFDR settings. Clicking a row makes that dataset
active. *Delete selected* removes the active dataset and its directory after a
confirmation; *Rename selected* gives it a new name everywhere (directory, table, run
manifest and the drawn Cytoscape network); *Clear All Datasets* empties the session; and
*Download Session Zip* saves every dataset with its results as one archive, which the
smaller card beside the table restores with *Load Session* (replacing the current
session).

**Step 1 — Parse** takes a dataset name, an input format (MaxQuant, DIA-NN, Pioneer,
FragPipe, MSstats or SAINT files) and the matching uploads. A note under the format
select says which file the format needs and which of its columns are read, and
*Download ED template* writes an experimental design CSV with the expected header
(Experiment Name, Type, Bait, Replicate, Bait ID, Group) and two example rows: one row
per experiment, with T for a test bait or C for a control, the bait and replicate, and
optionally the bait's UniProt ID and a control group. The design table is editable
before parsing and the edited table is what gets parsed. *Parse Data* writes the SAINT
and CompPASS inputs, adds the dataset to the table and makes it active. When it
finishes, a summary strip under the button gives the experiment, control, prey and
interaction counts and a *Score now* button that jumps to the Score card.

**Step 2 — Score** scores the active dataset (its name is in the card header) and picks
the organism, the imputation method for missing control values (AFT imputation applies
to intensity data; a spectral-count dataset offers Default only), how π is estimated
(imputation 2 only), the sigma floor's minimum observation count (imputation 2 and 3;
default 4, 0 disables it), the number of WDFDR permutations, and for human data whether
Human Cell Map evidence is excluded from BioGRID. *Score Data* runs SAINTexpress and
CompPASS, then annotates the result; the row's State becomes `scored`.

Parsing, scoring, feature enrichment and sending to Cytoscape each run under a progress
modal that lists the job's stages with the current one marked, the elapsed time and the
last lines of the dataset's `proximate.log`. On success the modal closes and a
notification reports the outcome; on failure it stays open with the error and the log
so the cause can be read before closing it. The page stays usable for other tabs, but a
second job cannot start until the first finishes.

## Data Thresholding

Quality control for the active dataset. Two PCA plots, available as soon as the dataset
is parsed, show experiments positioned by their protein profile (replicates should
cluster) and preys by their pattern across experiments, with controls for imputation,
normalization and a minimum detection fraction. Once the dataset is scored, the sidebar
thresholds drive three metrics: the median number of preys per bait that pass, how much
more often passing interactions are already known in BioGRID than before filtering, and
the mean number of known BioGRID partners each passing prey has among the other passing
preys. A bait selector restricts the metrics to one bait, and a per-bait scatter shows
each prey's SAINT score against fold change (log scale) with the sidebar's SAINT cutoff
as a dashed line. In the experiment PCA, controls form one grey series and each test
bait has one color and one legend entry; with more than 20 baits the legend is hidden
until *Show legend* is ticked. Below the metrics, *Passing interactions* lists every
interaction passing the sidebar thresholds for the selected bait (or all baits): bait,
prey gene and UniProt ID, SAINT, BFDR, WD, WDFDR, fold change, abundance (AvgSpec or
AvgIntensity), and the BioGRID and multivalidated flags, rounded to 3 decimals. The
table filters and sorts, and *Download CSV* saves the same rows. Plots can be exported
as PNG.

## Protein Feature Analysis

For each bait, the preys passing the sidebar thresholds (the foreground) are tested for
over-representation of protein features against every prey detected in the run: Gene
Ontology terms (cellular component, biological process, molecular function) and UniProt
sequence features (motifs, regions, repeats, compositions, domains). *Run Analysis*
computes the table and writes `Feature_enrichment.csv` into the dataset directory; the
heatmap shows the top features of one feature type (labels are cut at 60 characters;
the figure grows with the number of rows), and the table can be downloaded with p-value
and enrichment cutoffs. The download carries each feature's full name and, for Gene
Ontology terms, its GO id.

## Network Comparison

Compares two baits of the active dataset, each with its own threshold set. Both sets
start from the sidebar thresholds when the bait lists load, and *Use shared thresholds*
copies the sidebar's current values into one set. Baits from different datasets cannot
be compared. A three-panel
volcano plot shows preys present under both baits (log2 fold change of mean intensity
against the adjusted t-test p-value) flanked by preys detected under only one bait,
plotted by their SAINT BFDR. A Venn diagram counts preys passing thresholds in only A,
only B, or both, and the gene lists behind each region are shown as text to copy.

## Cytoscape

Sends the active dataset's interactions passing the sidebar thresholds into a Cytoscape
desktop on the same machine and keeps a link to it. The *Style preset* select picks the look: *Abundance/
literature weighted* draws orange diamond baits and blue-grey preys with bait-prey width
following abundance and BioGRID edges thickened by publications; *Highlight knowns* draws
black ellipse baits and grey rounded preys, puts a red border on preys whose interaction
with a drawn bait is BioGRID-multivalidated, and uses grey edges of one width with black
complex edges. Picking a preset also sets the send options it is designed for, which the
user can still change: which prey-prey BioGRID edges to add (optionally only multivalidated
pairs, thickened by publication count), CORUM complex edges for human data, what edge width
encodes and the layout algorithm. The other send options choose which baits to draw and the
node label policy. Each send replaces the previous ProxiMate network. With a network drawn:
*Apply Thresholds* hides edges that fail tighter thresholds without redrawing, so a hand
layout survives (loosening needs a new send); *Apply Style* applies the preset and edge
options in place; *Read Selection* lists the selected nodes with the
scores of their edges; the hide/show buttons act on edges around the selection; the
selection tools pick a bait's loners or satellites or the nodes a relation names for a
seed; *Cluster and Repack* runs Leiden over the selected preys and arranges them by
community with the chosen cluster layout (circles, or a force-directed or Kamada-Kawai
layout of each community inside its own disc), coloring communities of two or more nodes
unless *Recolor communities* is off, and leaving selected baits in place; *Save Layout*
writes every drawn node's
position to a named `id,x,y` file under the dataset's `cytoscape/layouts/` folder, merging
into an earlier save of that name so one layout can cover several bait subsets; *Load
Layout* moves the drawn nodes a saved layout names back to those positions and leaves
the rest where they are; *Export PNG* writes an image under the dataset's `cytoscape/`
folder; *Unlock* frees a view that stopped panning. The status line names the drawn
dataset and the thresholds it was drawn at, and the activity panel lists every
operation with who did it (`gui` or `mcp`).

## Downloads

Exports for the active dataset. Ready-made presets (annotated scores, Cytoscape edge
table, gene lists pooled or per bait, a ProHits-viz table, SAINT inputs, enriched
features) and a custom column picker. The *Apply thresholds* switch, off by default,
filters score-based exports by the sidebar thresholds its label names; off, every row
is exported. A batch export zips everything for the dataset.
