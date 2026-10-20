"""Central registry of in-app help text.

Every user-facing tooltip string lives in TOOLTIPS, one entry per line, each
tagged with a "# tooltip" comment — search for that marker to edit the text.
Keep entries to one or two plain-language sentences (tests cap them at 250
characters).
"""

from shiny import ui

from thresholds import PRESETS, describe

TOOLTIPS = {
    # ---- Sidebar: the active dataset and the thresholds every tab reads ----
    "active_dataset": "The dataset every tab works on. Pick one here or click its row in the session table; parsing a dataset makes it active.",  # tooltip
    "saintscore": "SAINT's confidence that a prey is a true interactor of the bait, from 0 to 1. Higher is better; 0.7-0.9 are common cutoffs.",  # tooltip
    "bfdr": "Bayesian False Discovery Rate: the estimated fraction of false interactions among those kept at this cutoff. Lower is stricter; 0.05 is a common cutoff.",  # tooltip
    "wd": "CompPASS WD score: high for preys that are abundant and seen with few baits (specific). Higher is better; 0 disables this filter.",  # tooltip
    "wdfdr": "False discovery rate for the WD score, estimated by shuffling the data. Lower is stricter; 1.0 disables this filter. Only meaningful when the dataset was scored with WDFDR iterations > 0.",  # tooltip
    "preset_saint_only": f"Sets {describe(PRESETS['SAINT only'])}: keeps the interactions SAINT calls confident and ignores the CompPASS scores.",  # tooltip
    "preset_wdfdr_only": f"Sets {describe(PRESETS['WDFDR only'])}: keeps the interactions the CompPASS WD score calls specific and ignores SAINT.",  # tooltip
    "preset_both": f"Sets {describe(PRESETS['Both'])}: keeps only interactions that pass both the SAINT and the WDFDR cut.",  # tooltip
    # ---- Data & Scoring tab ----
    "dataset_name": "A short name for this analysis. It labels the dataset everywhere in the app and in downloaded files.",  # tooltip
    "input_format": "The software that produced your quantification file. This sets which upload fields appear below.",  # tooltip
    "quant_type": "Which measurement to score. Intensity = summed raw signal; LFQ = label-free quantification, normalized across runs; Spectral Counts = number of spectra identifying each protein.",  # tooltip
    "pg_file": "The proteinGroups.txt table from MaxQuant's combined/txt output folder.",  # tooltip
    "diann_matrix_file": "The report.pg_matrix.tsv protein-group matrix written by DIA-NN.",  # tooltip
    "pioneer_matrix_file": "The protein_groups_wide.tsv table written by Pioneer's SearchDIA. Run columns are the MS file names without extension.",  # tooltip
    "fragpipe_file": "The combined_protein.tsv table from a FragPipe output folder.",  # tooltip
    "msstats_file": "The ProteinLevelData.csv table exported by MSstats after summarization.",  # tooltip
    "ed_template": "A CSV with the experimental design header and two example rows to fill in.",  # tooltip
    "ed_file": "CSV describing your experiments. Columns: Experiment Name, Type (T = test, C = control), Bait, Replicate, Bait ID. Optional Group column pairs tests with controls; controls may list several groups or * for all.",  # tooltip
    "ed_file_msstats": "Same experimental design CSV, but each Experiment Name must match an originalRUN value in ProteinLevelData.csv.",  # tooltip
    "saint_bait": "SAINT bait.txt: one row per experiment - experiment name, bait name, and T (test) or C (control).",  # tooltip
    "saint_prey": "SAINT prey.txt: one row per protein - protein ID, sequence length, and gene name.",  # tooltip
    "saint_interaction": "SAINT interaction.txt: one row per protein per experiment - experiment name, bait, prey, and quantity.",  # tooltip
    "organism": "The species your samples came from. Chooses which annotation databases (BioGRID, UniProt, HPA, CORUM) are used for annotation and enrichment.",  # tooltip
    "exclude_hcm": "Annotate against a BioGRID copy with all Human Cell Map (Go et al. 2021) evidence removed. HCM is itself a BioID screen, so it makes proximity-labeling hits look more 'known' than they are. Interactions with other evidence are kept.",  # tooltip
    "imputation_method": "How missing control values are filled in before SAINT scoring. Default = no fill-in; the AFT options model the chance an intensity is missing because it fell below detection, so they apply to intensity data only. Use Default for MSstats input.",  # tooltip
    "pi_method": "π is the estimated share of missing values that are true absences rather than below-detection signals. Choose to estimate it from all control baits or fit it from a single one.",  # tooltip
    "pi_bait": "The control bait whose replicates are used to fit π. Pick one with several replicates.",  # tooltip
    "aft_min_obs": "Sigma floor for the AFT fit. A prey seen in fewer than this many runs is fitted with the dataset's median per-prey spread instead of its own, which keeps its imputed control mean from being pinned to the bait intensity in small designs. 4 is recommended; 0 turns the floor off.",  # tooltip
    "wdfdr_iterations": "Number of data shuffles used to estimate the WD score's false discovery rate (WDFDR). More iterations are slower but more stable; 0 skips WDFDR entirely.",  # tooltip
    "delete_dataset": "Removes the active dataset and every result in its directory, after a confirmation. A network drawn from it stays in Cytoscape.",  # tooltip
    "rename_dataset": "Gives the active dataset a new name everywhere: its directory, the session table, its run manifest and the drawn Cytoscape network.",  # tooltip
    "clear_datasets": "Removes every dataset in this session, including results stored on the server. Download a session zip first if you want to keep them.",  # tooltip
    "download_session": "Saves every dataset in this session, with its results, as one zip you can restore later with Upload Session Zip.",  # tooltip
    "session_file": "A session zip downloaded earlier. Loading it replaces the current session with all of its datasets and results.",  # tooltip
    # ---- Quality Control tab ----
    "pca_imputation": "How to handle proteins not detected in every experiment. Row minimum = fill with that protein's smallest observed value; Zero = fill with 0; Drop = keep only fully observed preys.",  # tooltip
    "pca_normalization": "Scaling applied before PCA. Z-score centers each protein; log2 + Z-score first compresses large intensity ranges (good for raw intensities); None uses values as-is.",  # tooltip
    "pca_min_detection": "Keep only proteins detected in at least this fraction of experiments. Raising it removes sparsely observed proteins that can distort the PCA.",  # tooltip
    "experiment_pca": "Each point is one experiment, positioned by its overall protein profile. Replicates of the same bait should cluster; a stray point may be a failed run. Controls are one grey series.",  # tooltip
    "prey_pca": "Each point is one prey protein, positioned by its pattern across experiments. Preys that behave alike sit close together.",  # tooltip
    "prey_pca_color": "Color the prey points, e.g. by how many experiments each prey was detected in.",  # tooltip
    "qc_bait": "Restrict the plots and metrics to one bait, or choose All to pool every bait.",  # tooltip
    "metric_network_size": "Median number of preys per bait passing the current thresholds - the size of a typical network at these settings.",  # tooltip
    "metric_enrichment": "How much more often the kept interactions are already known in BioGRID, compared to before filtering. Higher suggests the thresholds favor real interactions.",  # tooltip
    "qc_results": "Every interaction passing the sidebar thresholds, for the bait selected above or all baits. Filter and sort in the table; Download CSV saves the same rows.",  # tooltip
    "metric_degree": "Average number of known BioGRID partners each kept prey has among the other kept preys. Densely connected networks often reflect real complexes.",  # tooltip
    # ---- Feature Enrichment tab ----
    "feature_analysis": "Tests whether protein features are over-represented: preys passing all four thresholds (the foreground) are compared against all detected preys.",  # tooltip
    "feature_type": "Category of features to display. GO_BP / GO_CC / GO_MF = Gene Ontology Biological Process, Cellular Component, Molecular Function; the rest are UniProt sequence features.",  # tooltip
    "num_features": "How many of the top-scoring features to show in the heatmap.",  # tooltip
    "download_pvalue_threshold": "Keep features whose p-value, adjusted for testing many features at once, is at or below this. 0.05 is a common cutoff.",  # tooltip
    "download_enrichment_threshold": "Keep features at least this many times more frequent in the network than expected by chance.",  # tooltip
    # ---- Compare Baits tab ----
    "comp_bait": "The bait whose interaction network forms this side of the comparison.",  # tooltip
    "comp_use_shared": "Copy the sidebar's four thresholds into this network's sliders.",  # tooltip
    "volcano_plot": "Compares prey abundance between the two baits: x = fold change, y = statistical confidence. Points far up and to either side differ most reliably.",  # tooltip
    "venn": "Counts of preys passing thresholds in only network A, only network B, or both.",  # tooltip
    "gene_lists": "The gene names behind each region of the Venn diagram, ready to copy into other tools.",  # tooltip
    # ---- Cytoscape tab ----
    "cy_baits": "Baits whose networks to draw. Leave empty to draw every bait in the dataset.",  # tooltip
    "cy_prey_prey": "Also draw grey edges between preys that BioGRID reports as interacting, so complexes show as clusters.",  # tooltip
    "cy_labels": "Which node names to draw. Baits only keeps large networks readable; click a node in Cytoscape to see its name.",  # tooltip
    "cy_layout": "The Cytoscape layout algorithm applied when the network is sent. You can re-layout in Cytoscape afterward.",  # tooltip
    "cy_send": "Build the network passing the thresholds and draw it in Cytoscape, replacing the previous ProxiMate network.",  # tooltip
    "cy_rethreshold": "Hide edges that fail the current thresholds without redrawing, so your hand layout survives. Loosening needs a new send.",  # tooltip
    "cy_read_selection": "List the nodes selected in Cytoscape with the scores of their edges.",  # tooltip
    "cy_download_selection": "Save the selection table as shown, with its rounded scores, as a CSV file.",  # tooltip
    "cy_hide_selected": "Hide every edge touching a node selected in Cytoscape.",  # tooltip
    "cy_show_selected": "Show every edge touching a node selected in Cytoscape.",  # tooltip
    "cy_hide_unselected": "Hide every edge that does not touch a selected node, leaving the selection's neighborhood.",  # tooltip
    "cy_show_all": "Show every edge again.",  # tooltip
    "cy_layout_name": "A layout saved for the drawn dataset to load, or a new name to save under. Letters, digits and underscores.",  # tooltip
    "cy_save_layout": "Write where every drawn node sits to the named layout in the dataset's cytoscape/layouts folder. Saving again updates the drawn nodes and keeps rows for nodes not drawn, so one layout can cover several bait subsets.",  # tooltip
    "cy_load_layout": "Move the drawn nodes the layout names back to their saved positions and fit the view. Nodes the layout does not name stay where they are.",  # tooltip
    "cy_export": "Save the drawn network as a PNG in the dataset's cytoscape folder under the output directory.",  # tooltip
    "cy_unlock": "Release a view that no longer pans or zooms with the mouse.",  # tooltip
    "cy_edge_width": "What a bait-prey edge's width encodes. Abundance, WD and fold change are banded on a log scale; SAINT is linear from 0 to 1; Uniform draws every edge alike.",  # tooltip
    "cy_biogrid_scope": "Which BioGRID pairs count as prey-prey edges: every reported pair, or only those BioGRID marks multivalidated (seen in more than one study or system).",  # tooltip
    "cy_lit_weighted": "Thicken a BioGRID edge with the number of publications behind it. Off, every BioGRID edge is thin. Most pairs have one paper, so only well-studied pairs stand out.",  # tooltip
    "cy_corum": "Draw black edges between drawn proteins that are subunits of one CORUM complex the screen recovered. Human datasets only; the two criteria below decide which complexes count.",  # tooltip
    "cy_corum_min_members": "A complex draws only when at least this many of its subunits are in the network.",  # tooltip
    "cy_corum_min_fraction": "A complex draws only when one bait, counted with its preys, covers at least this share of the complex's full membership. 0.5 means half the subunits.",  # tooltip
    "cy_style_preset": "A base look for the network. Picking one also sets the send options it is designed for (edges, width, layout); Send draws with both, Apply Style recolors the drawn network in place.",  # tooltip
    "cy_restyle": "Apply the style preset, edge width, BioGRID scope and publication weighting to the drawn network as a style update. Nothing moves and clustered nodes keep their colors; adding prey-prey or complex edges needs a new send.",  # tooltip
    "cy_select_loners": "With exactly one bait selected in Cytoscape, select it with the preys whose only visible neighbor of any kind it is, so the group drags as one.",  # tooltip
    "cy_select_satellites": "With exactly one bait selected, select it with its own preys and the two-bait preys that currently sit nearer to it than to their other bait.",  # tooltip
    "cy_rel_seed": "The bait or protein the relation starts from, by name or UniProt accession, as drawn in Cytoscape.",  # tooltip
    "cy_rel_kind": "Interactors and singletons need a bait seed. Partners are BioGRID or complex neighbors of any node. Co-complex members need the CORUM layer drawn.",  # tooltip
    "cy_rel_cuts": "Optional cuts on the relation: SAINT, BFDR and abundance apply to interactors, publications to BioGRID partners. Leave blank for no cut.",  # tooltip
    "cy_rel_replace": "Make the related nodes the Cytoscape selection.",  # tooltip
    "cy_rel_add": "Add the related nodes to whatever is already selected in Cytoscape.",  # tooltip
    "cy_cl_resolution": "Leiden resolution. Above 1 splits the selection into more, smaller communities; below 1 merges them.",  # tooltip
    "cy_cl_seed": "Random seed for the Leiden run, so the same selection clusters the same way again.",  # tooltip
    "cy_cl_lit_weight": "How much a BioGRID or complex edge pulls relative to a bait-prey edge. 0 clusters on bait-prey edges alone.",  # tooltip
    "cy_cl_layout": "How each community is arranged after clustering. Circles puts each on its own ring, tiled in rows. Force-directed and Kamada-Kawai lay each community out by its own edges inside a disc sized by its membership and place the discs by the edges between communities; better for large selections.",  # tooltip
    "cy_cl_recolor": "Color each community of two or more nodes with its own fill. A community of one node always keeps its color.",  # tooltip
    "cy_cluster": "Run Leiden over the selected preys, color communities of two or more nodes, and arrange each community by the cluster layout inside the box the selection occupies. Only selected preys move; selected baits stay put.",  # tooltip
    # ---- Export tab ----
    "dl_apply_thresholds": "On, score-based exports (Annotated Scores, Cytoscape, Gene List, ProHits-viz, Custom) keep only rows passing the sidebar thresholds; off exports every row. Experimental Design, SAINT Inputs and Enriched Features are never filtered.",  # tooltip
    "dl_preset": "Ready-made export formats for common next steps. A preset appears only when the dataset has the files it needs.",  # tooltip
    "dl_groups": "Toggle whole families of related columns (or files, for SAINT Inputs) in and out of the export.",  # tooltip
    "custom_columns": "Pick exactly the columns you want. Type to search; click a selected column and press Delete to remove it.",  # tooltip
    "dl_genelist_mode": "Pooled = one deduplicated list of every gene passing thresholds; Per bait = a two-column table of bait and gene.",  # tooltip
    "dl_prohits_abundance": "The measurement written to the Abundance column that ProHits-viz uses for dot size: AvePSM (average spectral counts) or AvgIntensity.",  # tooltip
    "batch_export": "One zip with everything for this dataset: merged data, annotated scores, enrichment results (if run), and the SAINT input files.",  # tooltip
}


# What each input format needs, shown under the format select and served with the
# documentation.  Every format but SAINT also takes the experimental design CSV.
FORMAT_NOTES = {
    "MaxQuant": "Upload proteinGroups.txt from the combined/txt folder plus the experimental design CSV. Read: Majority protein IDs, Gene names, Reverse, Potential contaminant, and the Intensity, LFQ intensity or MS/MS count column of each experiment named in the design.",
    "DIA-NN": "Upload report.pg_matrix.tsv plus the experimental design CSV. Read: Protein.Group, Protein.Names, Genes and one quantity column per run; the design's Experiment Name values must match the run column headers exactly. Intensity only.",
    "Pioneer": "Upload protein_groups_wide.tsv from SearchDIA plus the experimental design CSV. Read: the protein group, gene and one quantity column per run, named by the MS file without its extension. Intensity only.",
    "FragPipe": "Upload combined_protein.tsv plus the experimental design CSV. Read: Protein ID, Gene and each experiment's Intensity, MaxLFQ Intensity or Spectral Count column.",
    "MSstats": "Upload ProteinLevelData.csv as exported by MSstats plus the experimental design CSV. Read: Protein, originalRUN and LogIntensities; each Experiment Name must match an originalRUN value. Already log2, normalized and imputed, so score with Imputation = Default.",
    "SAINT": "Upload bait.txt, prey.txt and interaction.txt as SAINTexpress takes them (tab-separated, no header). The design is derived from bait.txt; bait protein IDs can be typed into the table before parsing.",
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
        "saint_bait", "saint_prey", "saint_interaction", "organism", "exclude_hcm",
        "imputation_method", "pi_method", "pi_bait", "aft_min_obs", "wdfdr_iterations",
        "delete_dataset", "rename_dataset", "clear_datasets", "download_session", "session_file",
        "ed_template",
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
