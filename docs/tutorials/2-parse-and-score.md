# 2. Parse and score

Turn a MaxQuant search into a scored dataset on the **Data & Scoring** tab:

- *Step 1 — Parse* reads the quantification and the experimental design.
- *Step 2 — Score* runs SAINT and CompPASS and annotates the result.
- **Duplicate selected** then makes a copy with one control left out, to score again.

The small example dataset is bundled with the GUI. To re-create the full Human Cell Map
datasets instead, [`examples/HCM_LFQ/`](../../examples/HCM_LFQ) holds the experimental
design and SAINT bait and prey files, and [`examples/README.md`](../../examples/README.md)
says where to get the MaxQuant output.

## Get the example dataset

![Step 1 header with the Download ED template and Download example dataset buttons](img/2-example-download.png)

In the *Step 1 — Parse* card, click **Download example dataset** and unzip
`proximate_example_HCM.zip`. It holds two files:

- `proteinGroups.txt`: a MaxQuant protein groups table (1.8 MB) with Intensity, LFQ and
  MS/MS count columns.
- `ED.csv`: its experimental design.

The data are a subset of the Human Cell Map (Go et al. 2021), 28 runs in all: eight
baits with two replicates each, plus twelve controls (eight of BirA-FLAG alone and four
of GFP-BirA-FLAG).

| Bait | Compartment |
| --- | --- |
| COIL | Cajal body |
| LMNA | nuclear lamina |
| NUP50 | nuclear pore |
| SEC61B (named `SEC61B_Nterm` in the design) | endoplasmic reticulum membrane |
| LAMP1 | lysosome |
| PEX14 | peroxisome |
| COX4I1 | mitochondrial inner membrane |
| RAB5A | early endosome |

## The experimental design

The experimental design (ED) is a CSV with one row per MS run:

| Column | Meaning |
| --- | --- |
| Experiment Name | The run's name, exactly as in the quantification file's column headers |
| Type | `T` for a bait (test) run, `C` for a control |
| Bait | The bait's name; replicates of one bait share it |
| Replicate | Replicate number |
| Bait ID | The bait's UniProt accession (optional; used by CompPASS and annotation) |
| Group | Optional control group |

Two rows of the example's `ED.csv`, a control and a bait run:

```
580_GFP__BirAFLAG_Go_RIPA_10102013.raw,C,GFP (BirA-FLAG),1,,
458_LMNA_BirAFLAG_Go_RIPA_06092013_BR1.raw,T,LMNA (BirA-FLAG),1,P02545,
```

Bait names carry the tag as a suffix, `LMNA (BirA-FLAG)`, and appear that way in every
select of the GUI. **Download ED template** writes an empty design with this header and
two example rows.

## Step 1 — Parse

![Step 1 with the example files uploaded and LFQ selected](img/2-parse-step1.png)

1. **Dataset Name**: `HCM_example`. Use letters, digits and underscores only.
2. **Input Format**: *MaxQuant*. The note under it says which file the format needs and
   which columns are read. The other formats are DIA-NN, Pioneer, FragPipe, MSstats and
   SAINT files.
3. **MaxQuant proteinGroups.txt file**: upload `proteinGroups.txt`.
4. **Experimental Design File**: upload `ED.csv`. Its 28 rows appear in the table on the
   right. Edit cells there before parsing, for example to fix a bait name or fill in a
   missing Bait ID; the edited table is what gets parsed.
5. **Quantification Type**: *LFQ*. MaxQuant offers *Intensity*, *LFQ* and *Spectral
   Counts*, each reading a different column per experiment.
6. Click **Parse Data**.

A progress modal lists the stages, the elapsed time and the end of the dataset's log:

![Progress modal while parsing](img/2-progress-modal.png)

The log shows what the parser removed: reverse hits, proteins only identified by site,
potential contaminants, and proteins found only in controls. Parsing the example takes a
second or two. The same modal appears for scoring, duplicating, Feature Enrichment and
Cytoscape sends, and stays open with the error and the log if a job fails.

A summary strip then appears under the button and the new dataset becomes active:

![Parse summary: experiments, controls, preys and interactions](img/2-parse-summary.png)

Check the counts: 16 experiments (the bait runs) and 12 controls, matching the `T` and
`C` rows of the design, 2477 preys kept, and 69356 experiment-prey quantities written.
In *Step 2 — Score*, *Dataset to score* now names `HCM_example`.

## Step 2 — Score

![Step 2 with Dataset to score set to HCM_example and the default settings](img/2-score-step2.png)

| Setting | What it does |
| --- | --- |
| **Dataset to score** | The dataset **Score Data** scores. It follows the active dataset and a newly parsed one. Picking another here changes only what gets scored, not the active dataset. |
| **Organism** | *Human* here. Picks the annotation databases (BioGRID, UniProt, Human Protein Atlas, CORUM). |
| **Exclude Human Cell Map evidence** | Annotates against a copy of BioGRID without the HCM paper's evidence. HCM is itself a BioID screen, so its evidence makes proximity hits look more "known" than they are. Tick it when benchmarking against BioGRID; this tutorial leaves it off. |
| **Imputation Method** | How missing control values are filled before SAINT. *Default* fills nothing. The two AFT options model values that fell below detection and apply to intensity data only. The preloaded `HCM_LFQ_1comp` and `HCM_LFQ_2comp` are the full HCM data scored with them. |
| **WDFDR Iterations** | How many data shuffles estimate the WD false discovery rate. 1000 is the default; 0 skips WDFDR. |

Click **Score Data**. The modal walks through *Running SAINT and CompPASS*, *Adding
protein annotation* and *Updating the datasets table*. The example takes about 20
seconds, and the sidebar badge changes from *parsed* to *scored*.

## Duplicate a dataset

**Duplicate selected** parses a copy of a dataset from its stored input, under a new
name and with an edited design. The original and its scores stay as they are. Use it
here to test how much the scores depend on one control, as you might for an outlier in
the experiment PCA (see [3. Quality control](3-quality-control.md)).

Click the `HCM_example` row in the *Datasets* table to make it active, then click
**Duplicate selected**:

![Duplicate dialog with one control left out](img/2-duplicate-modal.png)

- **Name of the copy**: `HCM_example_11controls`.
- **Leave out experiments**: pick `580_GFP__BirAFLAG_Go_RIPA_10102013.raw`, the first GFP
  control. The list offers every experiment of the design.
- The table below is the design. Double-click a cell to edit *Type*, *Bait*,
  *Replicate*, *Bait ID* or *Group*. Experiment names are fixed, because they must match
  the columns of the stored `proteinGroups.txt`.

Click **Duplicate**. The copy is parsed, becomes the active dataset and starts unscored;
*Dataset to score* switches to it:

![Step 2 with Dataset to score set to HCM_example_11controls](img/2-duplicate-score.png)

A notification reports the copy's size: 16 experiments and 11 controls. Click **Score
Data** again with the same settings. Then open **Quality Control** and switch the active
dataset between the two to compare their metrics at the same thresholds.

## The datasets table

![Datasets with HCM_example and HCM_example_11controls scored](img/2-datasets-table.png)

*Datasets* has one row per dataset: input format, quant type, experiment and control
counts, *State* (`parsed`, `scored`, or `enriched` after Feature Enrichment), when and
against which organism it was scored, and, scrolling right, its imputation and WDFDR
settings. Here it holds the three preloaded HCM datasets and the two example datasets.
Clicking a row makes that dataset active.

The buttons under the table act on the active (highlighted) dataset or the whole
session:

- **Delete selected** removes the dataset and its directory after a confirmation. A
  network drawn from it stays in Cytoscape.
- **Rename selected** renames the dataset everywhere: directory, table, run manifest and
  the drawn Cytoscape network.
- **Duplicate selected** opens the dialog described above.
- **Clear All Datasets** empties the session.
- **Download Session Zip** saves every dataset and its results as one archive, over a
  gigabyte with the three HCM datasets. The *Load a saved session* card next to the
  table restores an archive with **Load Session**, replacing the current session.

To clean up, select the `HCM_example_11controls` row, click **Delete selected** and
confirm, then do the same for `HCM_example`.

## What to look for

- **Counts match the design.** The summary strip's experiment and control counts should
  equal the `T` and `C` rows in your ED. A mismatch usually means an Experiment Name
  does not match a column header in the quantification file.
- **The filters are plausible.** In the parse log, "found only in controls" removes
  background proteins. If it removes most proteins, check that the Type column is not
  swapped.
- **The quant type is the one you meant.** LFQ, Intensity and Spectral Counts give
  different scores from the same file. The table's *Quant* column records which was used.
- **Scoring settings are recorded.** The *Imputation* and *WDFDR* columns say how each
  dataset was scored.
- **Duplicate rather than re-upload.** To drop a run or try another scoring setting,
  change a copy; the original keeps its scores for comparison.

Next: [3. Quality control](3-quality-control.md). Reference:
[`GUI/gui_docs.md`](../../GUI/gui_docs.md#data--scoring).
