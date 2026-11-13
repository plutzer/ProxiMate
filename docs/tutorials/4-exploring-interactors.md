# 4. Exploring interactors

Two tabs interpret a scored dataset. **Feature Enrichment** asks what each bait's
passing preys have in common, and **Compare Baits** sets two baits side by side. This
tutorial uses `HCM_LFQ_Default` at the *SAINT only* thresholds, with the cell-junction
baits CTNNB1 (β-catenin) and CTNNA1 (α-catenin).

## Feature Enrichment

For each bait, the preys passing the sidebar thresholds (the *foreground*) are tested
for over-represented features against every prey detected in the run. Features are:

- Gene Ontology terms: cellular component *GO_CC*, biological process *GO_BP*, molecular
  function *GO_MF*.
- UniProt sequence features: motifs, regions, repeats, compositions, domains.

1. Open **Feature Enrichment** and click **Run Feature Enrichment** in the *Foreground*
   card. The job takes about half a minute on this dataset and writes
   `Feature_enrichment.csv` into the dataset's directory. The dataset's State in the
   *Datasets* table becomes `enriched`.
2. In the *Top features* card, set *Select Feature Type* to `GO_CC` and *Number of
   Features to Display* to 15 (the default is 30).

![Heatmap of the top 15 GO cellular-component features across baits](img/4-enrichment-heatmap.png)

Rows are features, columns are baits, and color is log2 enrichment. Both axes are
clustered, so baits with similar localization sit together. With all 192 HCM baits only
some column labels are shown. Read the heatmap for its blocks (a group of baits sharing
*nuclear envelope*, another sharing *adherens junction*) and use the download below for
single baits. *Export Heatmap PNG* saves the figure as
`HCM_LFQ_Default_GO_CC_top15_heatmap.png`.

### Download the results for one bait

![Download filters set to CTNNB1](img/4-enrichment-download.png)

Under *Download enrichment results*, set *Bait* to **CTNNB1 (BirA-FLAG)**, leave
*Feature Type* at *All*, and keep the defaults *Max Adjusted p-value* 0.05 and *Min
Enrichment* 2. **Download Filtered Enrichment Results** saves
`enrichment_HCM_LFQ_Default_<date>_<time>.csv` with one row per feature:

| Column | Meaning |
| --- | --- |
| feature name, GO id | the feature |
| `k` of `n` | foreground preys carrying it, out of the foreground |
| `K` of `M` | detected preys carrying it, out of all detected preys |
| p-value, enrichment, adjusted p-value | the test result |

For CTNNB1 the top rows are *cell junction* (16 of its 31 passing preys, about 31-fold
enriched), *adherens junction*, *anchoring junction* and *catenin complex*: the passing
preys are where β-catenin is expected to be.

## Compare Baits

Open **Compare Baits**. The *Network A* and *Network B* cards each pick a bait and carry
their own four thresholds, with the same number boxes and sliders as the sidebar. They
start from the sidebar values; *Use shared thresholds* copies the sidebar's current
values into that card.

![Network A set to CTNNB1 and Network B to CTNNA1](img/4-compare-setup.png)

Set **Bait A** to *CTNNB1 (BirA-FLAG)* and **Bait B** to *CTNNA1 (C-BirA-FLAG)*, then
click **Compare**.

### Volcano plot

![Three-panel volcano plot of CTNNB1 against CTNNA1](img/4-compare-volcano.png)

- The **central panel** holds preys detected with both baits: log2 of the ratio of mean
  intensities (positive means more with CTNNB1) against the BH-adjusted t-test p-value.
  The dotted red line marks adjusted p = 0.05.
- The **flanking strips** hold preys detected with only one bait (left: only CTNNA1;
  right: only CTNNB1), plotted by their SAINT BFDR because there is no ratio to take.
- Color says where a prey passes the thresholds: *Network A only* (orange), *Network B
  only*, *Both* (purple) or *Neither* (grey). No prey passes for CTNNA1 alone here, so
  that category is missing from the legend.

*Export PNG* under the plot names the file after the dataset, both baits and both
threshold sets:
`HCM_LFQ_Default_CTNNB1_BirA-FLAG_vs_CTNNA1_C-BirA-FLAG_volcano_A-SAINT0.9_BFDR0.05_WD0_WDFDR1_B-SAINT0.9_BFDR0.05_WD0_WDFDR1.png`.
The Venn diagram's *Export PNG* uses the same pattern with `venn` in place of `volcano`.

### Overlap and gene lists

![Venn diagram of passing preys and the shared gene list](img/4-compare-venn.png)

The Venn diagram counts preys passing in A only, in both, and in B only: 18, 13 and 0.
CTNNB1 has 31 passing preys and CTNNA1 has 13, all shared with CTNNB1. The *Gene Lists*
tabs hold the genes behind each region as plain text to copy. *Both Networks* lists the
13 shared ones, among them CTNND1, JUP, APC, PARD3 and TJP1.

## What to look for

- **Enrichment should match the bait's biology.** CTNNB1's top terms are cell-junction
  terms. If a bait's top terms disagree with its known localization, recheck its QC.
- **Counts behind fold changes.** In the enrichment CSV, a large enrichment from `k` = 2
  or 3 preys is weaker evidence than a moderate one from 15.
- **Shared preys mean a shared neighborhood.** α- and β-catenin bind each other at
  adherens junctions, so a large shared set is expected. The CTNNB1-only preys
  (for example ERBIN, MAGI1, PLEKHA7, AFDN) point to what β-catenin reaches that
  α-catenin's tag does not.
- **Compare like with like.** The two cards' thresholds are independent; check that they
  match before reading differences in network size.

Next: [5. Networks and export](5-networks-and-export.md). Reference:
[`GUI/gui_docs.md`](../../GUI/gui_docs.md#feature-enrichment).
