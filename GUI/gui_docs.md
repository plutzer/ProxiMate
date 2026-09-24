# The ProxiMate GUI

ProxiMate scores proximity-labeling proteomics screens. A run moves left to right
through the navigation bar: parse the raw quantification into a **dataset**, score it,
then explore the scored network in the remaining tabs. Every dataset is a directory
under the output directory holding the SAINT inputs, `merged.csv` (SAINT and CompPASS
scores), `annotated_scores.csv` (with BioGRID, UniProt, HPA and CORUM annotation),
`proximate.log` and `run.json` (the run manifest). The "Datasets in this Session"
table at the top of the first tab lists every dataset and whether it is scored.

The same four score thresholds appear on most tabs and always mean the same thing:
SAINT score (higher is better), BFDR (lower is stricter), CompPASS WD score (higher is
better; 0 disables), and WDFDR (lower is stricter; 1.0 disables). An interaction passes
only when it clears all four. The preset buttons set all four at once.

## Network Scoring

Two cards. **Data Parsing** takes a dataset name, an input format (MaxQuant, DIA-NN,
Pioneer, FragPipe, MSstats or SAINT files) and the matching uploads, plus an
experimental design CSV that names each experiment, its type (T for a test bait, C for a
control), its bait and replicate, and optionally the bait's UniProt ID. The design table
is editable before parsing and the edited table is what gets parsed. *Parse Data* writes
the SAINT and CompPASS inputs and adds a row to the session table.

**Scoring** picks a parsed dataset, the organism, the imputation method for missing
control values (AFT imputation applies to intensity data; a spectral-count dataset
offers Default only), how π is estimated (imputation 2 only), the number of WDFDR
permutations, and for human data whether Human Cell Map evidence is excluded from
BioGRID. *Score Data* runs SAINTexpress and CompPASS, then annotates the result; the
row's Scored column flips to Yes. Scoring blocks the page until it finishes.

The session controls save every dataset as one zip (*Download Session*), restore such a
zip (*Upload Session Zip*, which replaces the current session), or remove everything
(*Clear Datasets*).

## Data Thresholding

Quality control for one scored dataset. The four threshold sliders drive three metrics:
the median number of preys per bait that pass, how much more often passing interactions
are already known in BioGRID than before filtering, and the mean number of known BioGRID
partners each passing prey has among the other passing preys. A bait selector restricts
the metrics to one bait. Two PCA plots show experiments positioned by their protein
profile (replicates should cluster) and preys by their pattern across experiments, with
controls for imputation, normalization and a minimum detection fraction. A per-bait
scatter shows each prey's scores against the SAINT cutoff. Plots can be exported as PNG.

## Protein Feature Analysis

For each bait, the preys passing the four thresholds (the foreground) are tested for
over-representation of protein features against every prey detected in the run: Gene
Ontology terms (cellular component, biological process, molecular function) and UniProt
sequence features (motifs, regions, repeats, compositions, domains). *Run Analysis*
computes the table and writes `Feature_enrichment.csv` into the dataset directory; the
heatmap shows the top features of one feature type, and the table can be downloaded
with p-value and enrichment cutoffs.

## Network Comparison

Compares two baits of one dataset, each with its own threshold set. A three-panel
volcano plot shows preys present under both baits (log2 fold change of mean intensity
against the adjusted t-test p-value) flanked by preys detected under only one bait,
plotted by their SAINT BFDR. A Venn diagram counts preys passing thresholds in only A,
only B, or both, and the gene lists behind each region are shown as text to copy.

## Cytoscape

Sends the interactions passing a threshold set into a Cytoscape desktop on the same
machine and keeps a link to it. The send options choose which baits to draw, whether to
add grey prey-prey edges from BioGRID (optionally only multivalidated pairs, thickened by
publication count), black CORUM complex edges for human data, node label policy, the
layout algorithm and what edge width encodes. Each send replaces the previous ProxiMate
network. With a network drawn: *Apply Thresholds* hides edges that fail tighter
thresholds without redrawing, so a hand layout survives (loosening needs a new send);
*Apply Edge Style* restyles in place; *Read Selection* lists the selected nodes with the
scores of their edges; the hide/show buttons act on edges around the selection; the
selection tools pick a bait's loners or satellites or the nodes a relation names for a
seed; *Cluster and Repack* runs Leiden over the selected preys and re-packs them by
community, leaving selected baits in place; *Record Positions* stores
where nodes were dragged; *Export PNG* writes an image under the dataset's `cytoscape/`
folder; *Unlock* frees a view that stopped panning. The status line names the drawn
dataset and the thresholds it was drawn at, and the activity panel lists every
operation with who did it (`gui` or `mcp`).

## Downloads

Exports for one dataset. Ready-made presets (annotated scores, Cytoscape edge table,
gene lists pooled or per bait, a ProHits-viz table, SAINT inputs, enriched features) and
a custom column picker. The four thresholds filter score-based exports only; a batch
export zips everything for the dataset.
