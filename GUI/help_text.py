"""Central registry of in-app help text.

Every user-facing tooltip string lives in TOOLTIPS, one entry per line, each
tagged with a "# tooltip" comment — search for that marker to edit the text.
Write entries as short phrases, about 10-20 words: what the control is, what
each option means, and one caveat when a user would otherwise get a wrong
result.  gui_docs.md carries the longer explanation.
"""

from shiny import ui

from thresholds import PRESETS, describe

TOOLTIPS = {
    # ---- Sidebar: the active dataset and the thresholds every tab reads ----
    "active_dataset": "The dataset every tab works on. Pick here or click a row in the Datasets table.",  # tooltip
    "saintscore": "SAINT confidence that a prey is a true interactor, 0 to 1. Higher is stricter; 0.7-0.9 is common.",  # tooltip
    "bfdr": "Bayesian false discovery rate: estimated share of false interactions kept. Lower is stricter; 0.05 is common.",  # tooltip
    "wd": "CompPASS WD score: high for abundant preys seen with few baits. Higher is stricter; 0 turns it off.",  # tooltip
    "wdfdr": "False discovery rate of the WD score. Lower is stricter; 1 turns it off. Needs WDFDR iterations > 0 at scoring.",  # tooltip
    "preset_saint_only": f"Sets {describe(PRESETS['SAINT only'])}: SAINT-confident interactions; CompPASS ignored.",  # tooltip
    "preset_wdfdr_only": f"Sets {describe(PRESETS['WDFDR only'])}: WD-specific interactions; SAINT ignored.",  # tooltip
    "preset_both": f"Sets {describe(PRESETS['Both'])}: interactions passing both the SAINT and WDFDR cuts.",  # tooltip
    # ---- Data & Scoring tab ----
    "dataset_name": "Short name for this analysis. Labels the dataset in the app and in downloaded files.",  # tooltip
    "input_format": "Software that produced the quantification file. Sets which upload fields appear.",  # tooltip
    "quant_type": "Measurement to score. Intensity: summed raw signal. LFQ: normalized across runs. Spectral Counts: spectra per protein.",  # tooltip
    "pg_file": "proteinGroups.txt from MaxQuant's combined/txt folder.",  # tooltip
    "diann_matrix_file": "report.pg_matrix.tsv, DIA-NN's protein-group matrix.",  # tooltip
    "pioneer_matrix_file": "protein_groups_wide.tsv from Pioneer's SearchDIA. Run columns: MS file names without extension.",  # tooltip
    "fragpipe_file": "combined_protein.tsv from the FragPipe output folder.",  # tooltip
    "msstats_file": "ProteinLevelData.csv exported by MSstats after summarization.",  # tooltip
    "ed_template": "Experimental design CSV with the header and two example rows.",  # tooltip
    "example_data": "Zip of a small Human Cell Map subset (Go et al. 2021): MaxQuant proteinGroups.txt and ED.csv, 8 baits, 12 controls. Parse as MaxQuant, any quant type.",  # tooltip
    "ed_file": "CSV, one row per run: Experiment Name, Type (T = test, C = control), Bait, Replicate, Bait ID. Optional Group pairs tests with controls; a control may list several groups or * for all.",  # tooltip
    "ed_file_msstats": "Experimental design CSV. Each Experiment Name must match an originalRUN value in ProteinLevelData.csv.",  # tooltip
    "saint_bait": "SAINT bait.txt, one row per experiment: experiment name, bait name, T (test) or C (control).",  # tooltip
    "saint_prey": "SAINT prey.txt, one row per protein: protein ID, sequence length, gene name.",  # tooltip
    "saint_interaction": "SAINT interaction.txt, one row per protein per experiment: experiment name, bait, prey, quantity.",  # tooltip
    "score_dataset": "Dataset that Score Data scores. Follows the active dataset; changing it here does not change the active dataset.",  # tooltip
    "organism": "Species of the samples. Selects the annotation databases (BioGRID, UniProt, HPA, CORUM).",  # tooltip
    "exclude_hcm": "Annotate against BioGRID without Human Cell Map evidence (Go et al. 2021). HCM is a BioID screen, so it inflates 'known' hits. Other evidence is kept.",  # tooltip
    "imputation_method": "Fills missing control values before SAINT. Default: no fill. AFT: models below-detection values; intensity data only. MSstats input: use Default.",  # tooltip
    "pi_method": "π: share of missing values that are true absences, not below-detection signal. Estimate from all control baits or from one.",  # tooltip
    "pi_bait": "Control bait used to fit π. Pick one with several replicates.",  # tooltip
    "aft_min_obs": "Sigma floor for the AFT fit. Preys seen in fewer runs than this use the dataset's median per-prey spread. 4 recommended; 0 turns it off.",  # tooltip
    "wdfdr_iterations": "Data shuffles used to estimate WDFDR. More: slower, more stable. 0 skips WDFDR.",  # tooltip
    "delete_dataset": "Deletes the active dataset and all its results, after confirmation. Its Cytoscape network stays drawn.",  # tooltip
    "rename_dataset": "Renames the active dataset everywhere: directory, Datasets table, run manifest, drawn Cytoscape network.",  # tooltip
    "duplicate_dataset": "Parses a copy of the active dataset under a new name. Leave out experiments or edit the design in the dialog. Copy starts unscored; original unchanged.",  # tooltip
    "clear_datasets": "Deletes every dataset and result in this session. Download a session zip first to keep them.",  # tooltip
    "download_session": "One zip of every dataset and its results. Restore it later with Upload Session Zip.",  # tooltip
    "session_file": "A session zip downloaded earlier. Loading it replaces the current session.",  # tooltip
    # ---- Quality Control tab ----
    "pca_imputation": "Fill for proteins missing from some experiments. Row minimum: the protein's smallest value. Zero: 0. Drop: keep only fully observed preys.",  # tooltip
    "pca_normalization": "Scaling before PCA. Z-score: centers each protein. log2 + Z-score: compresses large ranges first; use for raw intensities. None: values as-is.",  # tooltip
    "pca_min_detection": "Keep proteins detected in at least this fraction of experiments. Higher removes sparse proteins that distort the PCA.",  # tooltip
    "experiment_pca": "One point per experiment, placed by protein profile. Replicates should cluster; a stray point may be a failed run. Controls are grey.",  # tooltip
    "prey_pca": "One point per prey, placed by its pattern across experiments. Similar preys sit together.",  # tooltip
    "prey_pca_color": "Color for the prey points, such as number of experiments detected in.",  # tooltip
    "qc_bait": "Limit plots and metrics to one bait. All pools every bait.",  # tooltip
    "metric_network_size": "Median preys per bait passing the thresholds: the typical network size.",  # tooltip
    "metric_enrichment": "How much more often kept interactions are known in BioGRID than before filtering. Higher: thresholds favor real interactions.",  # tooltip
    "qc_results": "Interactions passing the sidebar thresholds, for the selected bait or all. Filter and sort here; Download CSV saves the same rows.",  # tooltip
    "metric_degree": "Mean known BioGRID partners per kept prey among the other kept preys. Dense networks often reflect real complexes.",  # tooltip
    # ---- Feature Enrichment tab ----
    "feature_analysis": "Tests each bait for over-represented protein features: preys passing all four thresholds against all detected preys.",  # tooltip
    "feature_type": "Feature category. GO_BP / GO_CC / GO_MF: Gene Ontology Biological Process, Cellular Component, Molecular Function. Others: UniProt sequence features.",  # tooltip
    "num_features": "Number of top-scoring features in the heatmap.",  # tooltip
    "download_pvalue_threshold": "Keep features with adjusted p-value at or below this. 0.05 is common.",  # tooltip
    "download_enrichment_threshold": "Keep features at least this many times more frequent than expected by chance.",  # tooltip
    # ---- Compare Baits tab ----
    "comp_bait": "Bait whose network forms this side of the comparison.",  # tooltip
    "comp_use_shared": "Copy the sidebar's four thresholds into this network's boxes.",  # tooltip
    "volcano_plot": "Prey abundance between the two baits. x: fold change. y: statistical confidence. Points far up and to a side differ most reliably.",  # tooltip
    "venn": "Preys passing thresholds in network A only, B only, or both.",  # tooltip
    "gene_lists": "Gene names behind each Venn region, ready to copy.",  # tooltip
    # ---- Cytoscape tab ----
    "cy_baits": "Baits to draw. Empty draws every bait.",  # tooltip
    "cy_prey_prey": "Add grey edges between preys that BioGRID reports as interacting, so complexes show as clusters.",  # tooltip
    "cy_labels": "Node names to draw. Baits only keeps large networks readable; click a node in Cytoscape for its name.",  # tooltip
    "cy_bait_labels": "Bait label text. Experiment name: keeps two constructs of one protein apart. Gene symbol: shorter. Apply Style relabels in place.",  # tooltip
    "cy_layout": "Cytoscape layout applied on send. Re-layout in Cytoscape afterward if needed.",  # tooltip
    "cy_send": "Draw the network passing the thresholds in Cytoscape. Replaces the previous ProxiMate network.",  # tooltip
    "cy_rethreshold": "Hide edges failing the current thresholds without redrawing; hand layout survives. Loosening needs a new send.",  # tooltip
    "cy_read_selection": "List the nodes selected in Cytoscape with their edges' scores.",  # tooltip
    "cy_download_selection": "Save the selection table as shown, rounded scores included, as CSV.",  # tooltip
    "cy_hide_selected": "Hide every edge touching a selected node.",  # tooltip
    "cy_show_selected": "Show every edge touching a selected node.",  # tooltip
    "cy_hide_unselected": "Hide every edge not touching a selected node, leaving the selection's neighborhood.",  # tooltip
    "cy_show_all": "Show every edge again.",  # tooltip
    "cy_layout_name": "Saved layout to load, or a new name to save under. Letters, digits and underscores.",  # tooltip
    "cy_save_layout": "Save every drawn node's position to the named layout. Saving again updates drawn nodes and keeps the rest, so one layout can cover several bait subsets.",  # tooltip
    "cy_load_layout": "Move drawn nodes to their saved positions and fit the view. Nodes not in the layout stay put.",  # tooltip
    "cy_export": "Save the drawn network as a PNG in the dataset's cytoscape folder.",  # tooltip
    "cy_unlock": "Release a view that no longer pans or zooms.",  # tooltip
    "cy_edge_width": "What bait-prey edge width encodes. Abundance, WD, fold change: log-scale bands. SAINT: linear, 0 to 1. Uniform: one width.",  # tooltip
    "cy_biogrid_scope": "BioGRID pairs drawn as prey-prey edges: all reported pairs, or only multivalidated ones (more than one study or system).",  # tooltip
    "cy_lit_weighted": "Thicken BioGRID edges by publication count. Off: all thin. Most pairs have one paper, so only well-studied pairs stand out.",  # tooltip
    "cy_corum": "Black edges between drawn subunits of a CORUM complex the screen recovered. Human only; the two criteria below decide which complexes count.",  # tooltip
    "cy_corum_min_members": "Draw a complex only when at least this many subunits are in the network.",  # tooltip
    "cy_corum_min_fraction": "Draw a complex only when one bait, with its preys, covers this share of all subunits. 0.5 = half.",  # tooltip
    "cy_style_preset": "Base look for the network. Also sets the send options it is designed for (edges, width, layout); change them before sending if needed.",  # tooltip
    "cy_restyle": "Apply preset, bait labels, edge width, BioGRID scope and publication weighting in place. Nothing moves; cluster colors stay. New prey-prey or complex edges need a new send.",  # tooltip
    "cy_select_loners": "With one bait selected: select it plus the preys whose only visible neighbor it is, so the group drags as one.",  # tooltip
    "cy_select_satellites": "With one bait selected: select it, its own preys, and two-bait preys sitting nearer to it than to their other bait.",  # tooltip
    "cy_rel_seed": "Starting bait or protein: experiment name, gene symbol or UniProt accession. If two drawn baits share a gene or accession, type the experiment name.",  # tooltip
    "cy_rel_kind": "Interactors, singletons: need a bait seed. Partners: BioGRID or complex neighbors of any node. Co-complex members: need CORUM edges drawn.",  # tooltip
    "cy_rel_cuts": "Optional cuts. SAINT, BFDR, abundance: interactors. Publications: BioGRID partners. Blank: no cut.",  # tooltip
    "cy_rel_replace": "Make the related nodes the Cytoscape selection.",  # tooltip
    "cy_rel_add": "Add the related nodes to the current Cytoscape selection.",  # tooltip
    "cy_cl_resolution": "Leiden resolution. Above 1: more, smaller communities. Below 1: fewer, larger.",  # tooltip
    "cy_cl_seed": "Random seed for Leiden, so a selection clusters the same way again.",  # tooltip
    "cy_cl_lit_weight": "Pull of a BioGRID or complex edge relative to a bait-prey edge. 0: bait-prey edges only.",  # tooltip
    "cy_cl_layout": "Arrangement of each community. Circles: one ring each, tiled in rows. Force-directed, Kamada-Kawai: laid out by its edges inside a disc; better for large selections.",  # tooltip
    "cy_cl_recolor": "Give each community of two or more nodes its own color. Single nodes keep theirs.",  # tooltip
    "cy_cluster": "Run Leiden on the selected preys, color communities, and arrange them inside the selection's box. Only selected preys move; baits stay put.",  # tooltip
    # ---- Export tab ----
    "dl_apply_thresholds": "On: score-based exports keep only rows passing the sidebar thresholds. Off: every row. Experimental Design, SAINT Inputs and Enriched Features are never filtered.",  # tooltip
    "dl_preset": "Export formats for common next steps. A preset appears only when the dataset has the files it needs.",  # tooltip
    "dl_groups": "Include or drop whole families of columns (files, for SAINT Inputs).",  # tooltip
    "custom_columns": "Columns to export. Type to search; click a selected column and press Delete to remove it.",  # tooltip
    "dl_genelist_mode": "Pooled: one deduplicated list of genes passing thresholds. Per bait: two columns, bait and gene.",  # tooltip
    "dl_prohits_abundance": "Abundance column ProHits-viz uses for dot size: AvePSM (average spectral counts) or AvgIntensity.",  # tooltip
    "batch_export": "One zip for this dataset: merged.csv, annotated_scores.csv, Feature_enrichment.csv (if run), and the SAINT input files.",  # tooltip
}


# What each input format needs, shown under the format select and served with the
# documentation.  Every format but SAINT also takes the experimental design CSV.
FORMAT_NOTES = {
    "MaxQuant": "Needs proteinGroups.txt (combined/txt folder) and the experimental design CSV.",
    "DIA-NN": "Needs report.pg_matrix.tsv and the experimental design CSV.",
    "Pioneer": "Needs protein_groups_wide.tsv (SearchDIA) and the experimental design CSV.",
    "FragPipe": "Needs combined_protein.tsv and the experimental design CSV.",
    "MSstats": "Needs ProteinLevelData.csv and the experimental design CSV. Already log2-transformed, normalized and imputed: score with Imputation Method = Default.",
    "SAINT": "Needs bait.txt, prey.txt and interaction.txt in SAINTexpress format.",
}


def tip(label, key):
    """Label followed by an info icon; hovering the icon shows TOOLTIPS[key].

    Missing keys raise KeyError so a bad reference fails at page build, not
    silently in the browser.
    """
    return ui.span(
        label,
        " ",
        ui.tooltip(ui.span("ⓘ", class_="info-icon"), TOOLTIPS[key]),
    )


# Which tab each tooltip belongs to, for documentation served outside the page.
SECTIONS = {
    "Sidebar": ["active_dataset", "saintscore", "bfdr", "wd", "wdfdr", "preset_saint_only",
                "preset_wdfdr_only", "preset_both"],
    "Data & Scoring": [
        "dataset_name", "input_format", "quant_type", "pg_file", "diann_matrix_file",
        "pioneer_matrix_file", "fragpipe_file", "msstats_file", "ed_file", "ed_file_msstats",
        "saint_bait", "saint_prey", "saint_interaction", "score_dataset", "organism", "exclude_hcm",
        "imputation_method", "pi_method", "pi_bait", "aft_min_obs", "wdfdr_iterations",
        "delete_dataset", "rename_dataset", "duplicate_dataset", "clear_datasets", "download_session", "session_file",
        "ed_template", "example_data",
    ],
    "Quality Control": [
        "pca_imputation", "pca_normalization",
        "pca_min_detection", "experiment_pca", "prey_pca", "prey_pca_color", "qc_bait",
        "metric_network_size", "metric_enrichment", "metric_degree", "qc_results",
    ],
    "Feature Enrichment": [
        "feature_analysis", "feature_type", "num_features", "download_pvalue_threshold",
        "download_enrichment_threshold",
    ],
    "Compare Baits": ["comp_bait", "comp_use_shared", "volcano_plot", "venn", "gene_lists"],
    "Cytoscape": [key for key in TOOLTIPS if key.startswith("cy_")],
    "Export": [
        "dl_apply_thresholds", "dl_preset", "dl_groups", "custom_columns", "dl_genelist_mode",
        "dl_prohits_abundance", "batch_export",
    ],
}
