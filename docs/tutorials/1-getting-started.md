# 1. Getting started

Learn the scores ProxiMate reports and the sidebar that every tab reads. The Human Cell
Map datasets ship with the image, so there is nothing to upload.

Start the GUI as shown in the [tutorial index](README.md) and open
http://localhost:3838.

## Concepts in plain language

**Baits and controls.** In a proximity-labeling screen (BioID, TurboID, APEX) each
*bait* protein is fused to a labeling enzyme. Mass spectrometry identifies the proteins
labeled near it: the *preys*. *Controls* are runs with the enzyme alone or fused to an
unrelated protein. A prey counts as a bait's interactor only if it is clearly more
abundant with that bait than in the controls. In the experimental design file, test runs
are type `T` and controls are type `C`.

| Score | What it measures | Direction |
| --- | --- | --- |
| **SAINT score** | Probability (0 to 1) that a prey is a true interactor, from its abundance with the bait against the controls. Computed by SAINTexpress. | Higher is better; the ProxiMate default cutoff is 0.9 |
| **BFDR** | Bayesian false discovery rate: the expected fraction of false interactions among those kept at a SAINT cutoff. | Lower is stricter; 0.05 is common |
| **WD** | CompPASS score for preys that are abundant, reproducible across replicates, and seen with few other baits. It catches "sticky" preys labeled by almost every bait, which SAINT can miss when controls are few. | Higher is better |
| **WDFDR** | False discovery rate for WD, estimated by reshuffling the data (the *WDFDR iterations* set at scoring). A dataset scored with 0 iterations has no WDFDR values. | Lower is stricter |

**The Human Cell Map (HCM).** A BioID screen of 192 baits covering the major compartments
of HEK293 cells (Go et al. 2021, PubMed 34079125). The image ships three scored versions
of the same input data, which differ only in imputation:

| Dataset | What it is |
| --- | --- |
| `HCM_LFQ_Default` | MaxQuant LFQ, 382 experiments (86 controls), scored without imputation. Used in these tutorials. |
| `HCM_LFQ_1comp` | Missing control values filled by one-component AFT imputation. |
| `HCM_LFQ_2comp` | Two-component AFT imputation. |

## The sidebar

![Sidebar with HCM_LFQ_Default active and the default thresholds](img/1-sidebar.png)

**Active dataset.** Every tab works on the dataset named here. Pick `HCM_LFQ_Default`
from the select, or click its row in the *Datasets* table on the Data & Scoring tab. The
badge reads *parsed* or *scored*; tabs that need scores show a message instead of their
controls until the dataset is scored.

**Thresholds.** An interaction *passes* only when it clears all four cutoffs:

| Threshold | Default | Meaning |
| --- | --- | --- |
| SAINT ≥ | 0.9 | minimum SAINT score |
| BFDR ≤ | 0.05 | maximum BFDR |
| WD ≥ | 0 | minimum WD score (0 switches the filter off) |
| WDFDR ≤ | 1 | maximum WDFDR (1 switches the filter off) |

Each cutoff is a number box with a slider. A typed value applies when you press Enter or
leave the box, a dragged slider when you release it. An out-of-range value shows a
warning under the boxes, and every output that reads the thresholds stays empty until
you fix it.

Quality Control metrics, Feature Enrichment, the Cytoscape network and thresholded
exports all follow these thresholds, so one change reaches every tab. Type `0.95` into
the SAINT box and press Enter. The slider follows, and no preset button is highlighted
because the values match none of them.

![SAINT threshold box set to 0.95 with its slider following](img/1-threshold-box.png)

**Presets.** Three buttons set all four thresholds together. The button matching the
current values is highlighted:

- *SAINT only*: SAINT ≥ 0.9, BFDR ≤ 0.05, WDFDR off. This is the starting state.
- *WDFDR only*: WDFDR ≤ 0.05, SAINT off.
- *Both*: SAINT ≥ 0.9, BFDR ≤ 0.05 and WDFDR ≤ 0.05.

Click **Both**. SAINT returns to 0.9, WDFDR moves to 0.05 and *Both* is highlighted:

![Thresholds after clicking the Both preset](img/1-preset-both.png)

Fewer interactions pass, because each must satisfy SAINT and CompPASS. Click **SAINT
only** to return to the defaults the remaining tutorials use.

Below the presets are two collapsed panels. *Agent access* holds the commands for
connecting an AI agent to ProxiMate's MCP server. *About and feedback* holds the *Report
a bug or request a new feature* link and the repository link.

## What to look for

- **Blank tab?** Check the badge: a *parsed* dataset has no scores yet.
- **Fewer interactions than expected?** Thresholds are shared, so look at the sidebar
  first. A stray WD or WDFDR cutoff is the usual cause.
- **Which preset?** *SAINT only* suits a screen with many controls, like HCM. *Both* is
  stricter and trims promiscuous preys that SAINT lets through.
- **WDFDR without iterations.** With WDFDR below 1 on a dataset scored with 0 WDFDR
  iterations, a warning appears and nothing passes.

Next: [2. Parse and score](2-parse-and-score.md). The full reference for every tab is in
[`GUI/gui_docs.md`](../../GUI/gui_docs.md).
