from functools import partial
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui, run_app
import plotly.graph_objects as go
import plotly.express as px
from shinywidgets import output_widget, render_widget, render_plotly
import logging
import platform
import pandas as pd
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'Scripts'))
from ed_exceptions import ProxiMateError
import parse
import subprocess
import zipfile
import tempfile
import datetime
import shutil
import threading
import time
from QC_plots import pca_plot, prepare_pca_matrix, prey_pca_plot as plot_prey_pca, prey_gene_names, detection_counts, reduce_categorical, saint_known_retention, roc_plot, saint_scatter_plot as plot_saint_scatter, calculate_threshold_metrics
from Ann_Enrichment import plot_results, go_term_ids, GO_SOURCE_COLUMNS
from network_comparison import (
    load_and_filter_bait_data,
    calculate_volcano_data,
    create_volcano_plot,
    create_venn_diagram_matplotlib,
    create_volcano_plot_matplotlib
)
from plot_exports import pca_plot_matplotlib, prey_pca_matplotlib, saint_scatter_matplotlib
import download_presets as dp
from download_presets import DEFAULT_CUSTOM_COLUMNS
from help_text import tip, TOOLTIPS, FORMAT_NOTES
from thresholds import PRESETS, SLIDERS, DEFAULT, NO_FILTER, preset_id, active_preset, describe
import session_archive
import backend
import dataset_store
import mcp_registry
import cytoscape_ctl
from cytoscape_net import SEND_WARN_EDGES, STYLE_PRESETS, selection_table
import log_config
import provenance
from log_config import get_logger

logger = get_logger(__name__)

out_dir = os.environ.get("PROXIMATE_OUTPUT_DIR", "/Outputs")
backend.configure(out_dir)
# Shown in the sidebar.  The browser reaches the MCP port on the same host as the GUI,
# so "localhost" is right whenever the port is published alongside 3838.
MCP_URL = f"http://localhost:{os.environ.get('PROXIMATE_MCP_PORT', '3839')}/mcp"
MCP_ADD_COMMANDS = {
    'Claude Code': f"claude mcp add --transport http proximate {MCP_URL}",
    'Codex CLI': f"codex mcp add proximate --url {MCP_URL}",
}

STYLESHEET = os.path.join(os.path.dirname(os.path.abspath(__file__)), "www", "proximate.css")

# The sidebar layout's main pane is the scroll container, so a tab opened while scrolled
# down would start mid-page.  A clicked button keeps focus, which keeps its bslib tooltip
# open; the click hides the tooltip and mutes it until the pointer or focus leaves.
PAGE_SCRIPT = """
document.addEventListener('shown.bs.tab', () => {
  document.querySelectorAll('.bslib-sidebar-layout > .main').forEach(el => el.scrollTop = 0);
});
document.addEventListener('click', e => {
  const trigger = e.target.closest('bslib-tooltip > *');
  if (!trigger) return;
  trigger.dataset.tooltipMuted = '1';
  const tip = bootstrap.Tooltip.getInstance(trigger);
  if (tip) tip.hide();
});
document.addEventListener('show.bs.tooltip', e => {
  if (e.target.dataset.tooltipMuted) e.preventDefault();
});
for (const event of ['mouseleave', 'focusout']) {
  document.addEventListener(event, e => {
    if (e.target.dataset) delete e.target.dataset.tooltipMuted;
  }, true);
}
"""

# More test baits than this and the experiment PCA hides its legend until asked.
PCA_LEGEND_MAX_BAITS = 20
EMPTY_BAITS = pd.DataFrame(columns=["Experiment Name", "Bait", "Type", "Bait ID"])
EMPTY_ED = pd.DataFrame(columns=["Experiment Name", "Type", "Bait", "Replicate", "Bait ID", "Group"])


def threshold_sliders(prefix, values):
    """The four cutoff sliders, with ids ``{prefix}_{suffix}``, starting at ``values``."""
    return [ui.input_slider(f"{prefix}_{suffix}", tip(label, suffix), min=0.0, max=high,
                            value=values[key], step=step)
            for suffix, key, label, high, step in SLIDERS]


app_ui = ui.page_navbar(
    ui.nav_spacer(),
    ui.nav_panel("Data & Scoring",
                 ui.layout_columns(
                    ui.card(
                        ui.card_header("Datasets in this Session"),
                        ui.output_data_frame("render_datasets"),
                        ui.div(
                            ui.tooltip(ui.input_action_button("delete_dataset", "Delete selected",
                                                              class_="btn-sm btn-outline-danger"),
                                       TOOLTIPS["delete_dataset"]),
                            ui.tooltip(ui.input_action_button("rename_dataset", "Rename selected",
                                                              class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["rename_dataset"]),
                            ui.tooltip(ui.input_action_button("clear_datasets", "Clear All Datasets",
                                                              class_="btn-sm btn-outline-danger"),
                                       TOOLTIPS["clear_datasets"]),
                            ui.tooltip(ui.download_button("download_session", "Download Session Zip",
                                                          class_="btn-sm"),
                                       TOOLTIPS["download_session"]),
                            class_="button-row"),
                    ),
                    ui.card(
                        ui.card_header("Load a saved session"),
                        ui.input_file("session_file", tip("Upload Session Zip", "session_file")),
                        ui.input_action_button("upload_session", "Load Session", class_="btn-sm"),
                    ),
                    col_widths=[9, 3]
                 ),
                 ui.layout_columns(
                    ui.card(
                        ui.card_header(
                            ui.span("Step 1 — Parse", class_="step-title"),
                            ui.span("Turn quantification files and an experimental design into a dataset.",
                                    class_="step-purpose")),
                        ui.layout_columns(
                            ui.input_text("dataset_name",
                                        tip("Dataset Name", "dataset_name"),
                                        placeholder="No spaces or special characters (/ \\ : * ? \" < > |)"),
                            ui.input_select("input_format", tip("Input Format", "input_format"),
                                           choices=["MaxQuant", "DIA-NN", "Pioneer", "FragPipe", "MSstats", "SAINT"],
                                           selected="MaxQuant"),
                            col_widths=[4, 8]
                        ),
                        ui.output_ui("format_note"),
                        ui.tooltip(ui.download_button("ed_template", "Download ED template",
                                                      class_="btn-sm btn-outline-secondary mb-3"),
                                   TOOLTIPS["ed_template"]),
                        ui.panel_conditional(
                            "input.input_format === 'MaxQuant'",
                            ui.layout_columns(
                                ui.div(
                                    ui.input_file("pg_file", tip("MaxQuant proteinGroups.txt file", "pg_file")),
                                    ui.input_file("ed_file", tip("Experimental Design File", "ed_file")),
                                    ui.input_select("quant_type", tip("Quantification Type", "quant_type"),
                                                  choices=["Intensity", "LFQ", "Spectral Counts"],
                                                  selected="Intensity")
                                ),
                                ui.output_data_frame("ed_table_mq"),
                                col_widths=[4, 8]
                            )
                        ),
                        ui.panel_conditional(
                            "input.input_format === 'DIA-NN'",
                            ui.layout_columns(
                                ui.div(
                                    ui.input_file("diann_matrix_file", tip("DIA-NN report.pg_matrix.tsv file", "diann_matrix_file")),
                                    ui.input_file("ed_file", tip("Experimental Design File", "ed_file"))
                                ),
                                ui.output_data_frame("ed_table_diann"),
                                col_widths=[4, 8]
                            )
                        ),
                        ui.panel_conditional(
                            "input.input_format === 'Pioneer'",
                            ui.layout_columns(
                                ui.div(
                                    ui.input_file("pioneer_matrix_file", tip("Pioneer protein_groups_wide.tsv file", "pioneer_matrix_file")),
                                    ui.input_file("ed_file", tip("Experimental Design File", "ed_file"))
                                ),
                                ui.output_data_frame("ed_table_pioneer"),
                                col_widths=[4, 8]
                            )
                        ),
                        ui.panel_conditional(
                            "input.input_format === 'FragPipe'",
                            ui.layout_columns(
                                ui.div(
                                    ui.input_file("fragpipe_file", tip("FragPipe combined_protein.tsv file", "fragpipe_file")),
                                    ui.input_file("ed_file", tip("Experimental Design File", "ed_file")),
                                    ui.input_select("quant_type", tip("Quantification Type", "quant_type"),
                                                  choices=["Intensity", "LFQ", "Spectral Counts"],
                                                  selected="Intensity")
                                ),
                                ui.output_data_frame("ed_table_fragpipe"),
                                col_widths=[4, 8]
                            )
                        ),
                        ui.panel_conditional(
                            "input.input_format === 'MSstats'",
                            ui.layout_columns(
                                ui.div(
                                    ui.input_file("msstats_file", tip("MSstats ProteinLevelData.csv file", "msstats_file")),
                                    ui.input_file("ed_file", tip("Experimental Design File", "ed_file_msstats")),
                                    ui.tags.small(
                                        "Note: MSstats data is already log2-transformed, normalized, and imputed. Set Imputation Method to 'Default' (0).",
                                        class_="help-text d-block mt-2"
                                    )
                                ),
                                ui.output_data_frame("ed_table_msstats"),
                                col_widths=[4, 8]
                            )
                        ),
                        ui.panel_conditional(
                            "input.input_format === 'SAINT'",
                            ui.layout_columns(
                                ui.div(
                                    ui.input_file("bait", tip("SAINT bait.txt file", "saint_bait")),
                                    ui.input_file("prey", tip("SAINT prey.txt file", "saint_prey")),
                                    ui.input_file("interaction", tip("SAINT interaction.txt file", "saint_interaction")),
                                    ui.input_select("quant_type", tip("Quantification Type", "quant_type"),
                                                  choices=["Intensity", "LFQ", "Spectral Counts"],
                                                  selected="Intensity")
                                ),
                                ui.output_data_frame("bait_table"),
                                col_widths=[4, 8]
                            )
                        ),
                        ui.input_action_button("parse_data", "Parse Data"),
                        ui.output_ui("parse_summary"),
                    ),
                    ui.card(
                        ui.card_header(
                            ui.span("Step 2 — Score", class_="step-title"),
                            ui.span(ui.output_text("score_header", inline=True), class_="step-purpose")),
                        # Organism selection for annotation databases.
                        # To add a new organism, add a choice here and sync with
                        # ORGANISMS config in Scripts/annotator.py and setup_datasets.py
                        ui.input_select("organism", tip("Organism", "organism"),
                            choices={"human": "Human (H. sapiens)",
                                     "mouse": "Mouse (M. musculus)",
                                     "yeast": "Yeast (S. cerevisiae)"},
                            selected="human"),
                        ui.panel_conditional(
                            "input.organism === 'human'",
                            ui.input_checkbox("exclude_hcm",
                                tip("Exclude Human Cell Map evidence", "exclude_hcm"),
                                value=False),
                        ),
                        ui.input_radio_buttons("imputation_method", tip("Imputation Method", "imputation_method"),
                                              choices={0: "Default", 
                                                    #    1: "Prey-specific",
                                                       2: "Two-component AFT",
                                                       3: "One-component AFT"}),
                        ui.panel_conditional(
                            "String(input.imputation_method) === '2'",
                            ui.input_radio_buttons(
                                "pi_method", tip("π Estimation Method", "pi_method"),
                                choices={"weighted_average": "Weighted avg of controls (≥3 replicates)",
                                         "single_bait": "Single control bait"},
                                selected="weighted_average"),
                            ui.panel_conditional(
                                "input.pi_method === 'single_bait'",
                                ui.input_select("pi_bait", tip("Control Bait for π Fit", "pi_bait"), choices=[])),
                        ),
                        ui.panel_conditional(
                            "String(input.imputation_method) === '2' || String(input.imputation_method) === '3'",
                            ui.input_numeric("aft_min_obs", tip("Sigma Floor: Min. Observations", "aft_min_obs"),
                                             value=4, min=0, step=1),
                        ),
                        ui.input_numeric("wdfdr_iterations", tip("WDFDR Iterations", "wdfdr_iterations"), value=1000),
                        ui.input_action_button("score_data", "Score Data"),
                        id="score_card",
                    ),
                    col_widths=[8, 4]
                 ),
            ),
    ui.nav_panel("Quality Control",
                    ui.output_ui("empty_state_thresholding"),
                    # The PCA plots read the parsed data, so they show for any dataset.
                    ui.panel_conditional(
                    "output.active_state !== 'none'",
                    # Preprocessing controls shared by both PCA plots below
                    ui.card(
                        ui.card_header("PCA Preprocessing"),
                        ui.layout_columns(
                            ui.input_select("pca_imputation", tip("Imputation", "pca_imputation"),
                                            choices={"row_min": "Row minimum",
                                                     "zero": "Zero",
                                                     "drop": "Drop incomplete preys"},
                                            selected="row_min"),
                            ui.input_select("pca_normalization", tip("Normalization", "pca_normalization"),
                                            choices={"zscore": "Z-score",
                                                     "log2_zscore": "log2 + Z-score",
                                                     "none": "None"},
                                            selected="zscore"),
                            ui.input_slider("pca_min_detection", tip("Min detection fraction", "pca_min_detection"),
                                            min=0.0, max=1.0, value=0.5, step=0.05),
                            col_widths=(4, 4, 4),
                        ),
                    ),

                    # Row 1: experiment-level and prey-level PCA side by side
                    ui.layout_columns(
                        ui.card(
                            ui.card_header(tip("Experiment PCA", "experiment_pca"),
                                           ui.output_ui("pca_legend_toggle"),
                                           class_="d-flex justify-content-between align-items-center"),
                            output_widget("raw_pca_plot"),
                            ui.download_button("download_pca_plot", "Export PNG", class_="btn-sm"),
                        ),
                        ui.card(
                            ui.card_header(tip("Prey PCA", "prey_pca")),
                            ui.layout_columns(
                                ui.input_select("prey_pca_color", tip("Color by", "prey_pca_color"),
                                                choices={"none": "None",
                                                         "detection": "Detection count"},
                                                selected="none"),
                                ui.panel_conditional(
                                    "input.prey_pca_color == 'saint'",
                                    ui.input_select("prey_pca_bait", "Bait", choices=[]),
                                ),
                                col_widths=(6, 6),
                            ),
                            output_widget("prey_pca_plot"),
                            ui.download_button("download_prey_pca_plot", "Export PNG", class_="btn-sm"),
                        ),
                        col_widths=(6, 6),
                    ),
                    ),

                    # Row 2: SAINT scatter | metrics, at the sidebar thresholds
                    ui.panel_conditional(
                    "output.active_state === 'scored'",
                    ui.input_select("qc_bait", tip("Select QC Bait", "qc_bait"), choices=["All"]),
                    ui.layout_columns(
                        ui.card(
                            ui.card_header("SAINT Score vs Fold Change"),
                            ui.output_ui("scatter_body"),
                            ui.download_button("download_scatter_plot", "Export PNG", class_="btn-sm"),
                        ),
                        ui.div(
                            ui.output_ui("metric_bait_label"),
                            ui.output_ui("metric_network_size"),
                            ui.output_ui("metric_enrichment"),
                            ui.output_ui("metric_degree"),
                        ),
                        col_widths=(8, 4),
                    ),
                    ui.card(
                        ui.card_header(tip("Passing interactions", "qc_results")),
                        ui.div(ui.output_text("qc_results_count", inline=True),
                               ui.download_button("download_qc_results", "Download CSV", class_="btn-sm"),
                               class_="d-flex justify-content-between align-items-center"),
                        ui.output_data_frame("qc_results"),
                    ),

                    # Keep these plots in code but hide them (for potential future use)
                    ui.panel_conditional(
                        "false",  # Never show
                        ui.card(
                            ui.card_header("Known Physical Interactions"),
                            output_widget("known_retention_plot")
                        ),
                        ui.card(
                            ui.card_header("ROC Curve for BioGRID Interactions"),
                            ui.input_radio_buttons("roc_known_type", "True Positive Type",
                                                   choices={"BioGRID":"All Physical Interactions",
                                                           "Multivalidated":"Multivalidated Physical Interactions"}),
                            output_widget("roc_curve_plot")
                        ),
                    ),
                    ),
    ),
    ui.nav_panel("Feature Enrichment",
                ui.output_ui("empty_state_feature_analysis"),
                ui.panel_conditional(
                "output.active_state === 'scored'",
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Foreground"),
                        ui.p("Preys passing the sidebar thresholds form the foreground tested for enrichment.",
                             class_="help-text"),
                        ui.tooltip(ui.input_action_button("feature_analysis", "Run Feature Enrichment"),
                                   TOOLTIPS["feature_analysis"]),
                    ),
                    ui.card(
                        ui.card_header("Top features"),
                        ui.input_select("feature_type", tip("Select Feature Type", "feature_type"), choices=["GO_CC", "GO_BP", "GO_MF", "Motifs", "Regions", "Repeats", "Compositions", "Domains"]),
                        ui.input_numeric("num_features", tip("Number of Features to Display", "num_features"), value=30, min=1, max=100),
                        ui.output_ui("feature_heatmap_body"),
                        ui.download_button("download_heatmap", "Export Heatmap PNG", class_="btn-sm"),
                        ui.p("Download Enrichment Results", class_="section-heading"),
                        ui.layout_columns(
                            ui.input_select("download_feature_type", tip("Feature Type", "feature_type"),
                                          choices=["All", "GO_CC", "GO_BP", "GO_MF", "Motifs", "Regions", "Repeats", "Compositions", "Domains"]),
                            ui.input_select("download_bait_filter", "Bait", choices=["All"]),
                            col_widths=(6, 6)
                        ),
                        ui.layout_columns(
                            ui.input_slider("download_pvalue_threshold", tip("Max Adjusted p-value", "download_pvalue_threshold"),
                                          min=0.0, max=1.0, value=0.05, step=0.01),
                            ui.input_slider("download_enrichment_threshold", tip("Min Enrichment", "download_enrichment_threshold"),
                                          min=0.0, max=10.0, value=2.0, step=0.1),
                            col_widths=(6, 6)
                        ),
                        ui.download_button("download_enrichment", "Download Filtered Enrichment Results"),
                    ),
                    col_widths=(3, 9),
                ),
                ),
    ),
    ui.nav_panel("Compare Baits",
        ui.output_ui("empty_state_network_comparison"),
        ui.panel_conditional(
        "output.active_state === 'scored'",
        # Two baits of the active dataset, each with its own threshold set.  The sets
        # start from the sidebar's whenever the bait lists load.
        ui.layout_columns(
            ui.card(
                ui.card_header("Network A"),
                ui.input_select("comp_bait_a", tip("Bait A", "comp_bait"), choices=[]),
                *threshold_sliders("comp_a", DEFAULT),
                ui.tooltip(ui.input_action_button("comp_use_shared_a", "Use shared thresholds",
                                                  class_="btn-sm btn-outline-secondary"),
                           TOOLTIPS["comp_use_shared"]),
            ),
            ui.card(
                ui.card_header("Network B"),
                ui.input_select("comp_bait_b", tip("Bait B", "comp_bait"), choices=[]),
                *threshold_sliders("comp_b", DEFAULT),
                ui.tooltip(ui.input_action_button("comp_use_shared_b", "Use shared thresholds",
                                                  class_="btn-sm btn-outline-secondary"),
                           TOOLTIPS["comp_use_shared"]),
            ),
            col_widths=(6, 6),
        ),

        # Compare button
        ui.div(
            ui.input_action_button("compare_networks", "Compare", class_="btn-primary"),
            class_="centered-action"
        ),

        # Middle section: Volcano plot
        ui.card(
            ui.card_header(tip("Differential Abundance Volcano Plot", "volcano_plot")),
            ui.output_ui("volcano_body"),
            ui.download_button("download_volcano_plot", "Export PNG", class_="btn-sm"),
            ui.p("Flanking strips hold preys detected under only one bait "
                 "(y = -log10 BFDR); the central panel holds shared preys "
                 "(y = -log10 BH-adjusted p).",
                 class_="help-text mt-2"),
        ),

        # Bottom section: Venn diagram and gene lists
        ui.layout_columns(
            ui.card(
                ui.card_header(tip("Network Overlap", "venn")),
                ui.output_ui("venn_body"),
                ui.download_button("download_venn_diagram", "Export PNG", class_="btn-sm"),
            ),
            ui.card(
                ui.card_header(tip("Gene Lists", "gene_lists")),
                ui.navset_tab(
                    ui.nav_panel("Network A Only",
                        ui.output_text_verbatim("genes_a_only"),
                    ),
                    ui.nav_panel("Network B Only",
                        ui.output_text_verbatim("genes_b_only"),
                    ),
                    ui.nav_panel("Both Networks",
                        ui.output_text_verbatim("genes_both"),
                    ),
                ),
            ),
            col_widths=(6, 6),
        ),
        ),
    ),
    ui.nav_panel("Cytoscape",
        ui.output_ui("empty_state_cytoscape"),
        ui.panel_conditional(
        "output.active_state === 'scored'",
        ui.layout_columns(
            ui.card(
                ui.card_header("Network"),
                ui.p("Network", class_="section-heading"),
                ui.input_selectize("cy_baits", tip("Baits", "cy_baits"), choices=[], multiple=True),
                ui.output_ui("cy_thresholds_line"),
                ui.p("Style", class_="section-heading"),
                ui.input_select("cy_style_preset", tip("Style preset", "cy_style_preset"),
                                choices={key: preset["label"] for key, preset in STYLE_PRESETS.items()}),
                ui.input_select("cy_labels", tip("Labels", "cy_labels"),
                                choices={"all": "All nodes", "baits": "Baits only", "none": "None"}),
                ui.input_select("cy_layout", tip("Layout", "cy_layout"),
                                choices=["force-directed", "kamada-kawai", "circular", "grid",
                                         "hierarchical", "degree-circle"]),
                ui.input_select("cy_edge_width", tip("Edge width", "cy_edge_width"),
                                choices={"abundance": "Abundance (intensity / spectral counts)",
                                         "SaintScore": "SAINT score", "WD": "WD score",
                                         "FoldChange": "Fold change", "uniform": "Uniform"}),
                ui.p("Edges", class_="section-heading"),
                ui.input_checkbox("cy_prey_prey", tip("Prey-prey BioGRID edges", "cy_prey_prey"), value=True),
                ui.input_select("cy_biogrid_scope", tip("BioGRID edges", "cy_biogrid_scope"),
                                choices={"all": "All pairs", "multivalidated": "Multivalidated only"}),
                ui.input_checkbox("cy_lit_weighted", tip("Weight BioGRID edges by publications", "cy_lit_weighted"),
                                  value=True),
                ui.p("Complexes", class_="section-heading"),
                ui.input_checkbox("cy_corum", tip("CORUM complex edges", "cy_corum"), value=False),
                ui.layout_columns(
                    ui.input_numeric("cy_corum_min_members", tip("Min subunits drawn", "cy_corum_min_members"),
                                     value=3, min=2, step=1),
                    ui.input_numeric("cy_corum_min_fraction", tip("Min share by one bait", "cy_corum_min_fraction"),
                                     value=0.5, min=0.0, max=1.0, step=0.05),
                    col_widths=(6, 6),
                ),
                ui.output_ui("cy_estimate"),
                ui.div(
                    ui.tooltip(ui.input_action_button("cy_send", "Send to Cytoscape", class_="btn-primary"),
                               TOOLTIPS["cy_send"]),
                    ui.tooltip(ui.input_action_button("cy_rethreshold", "Re-apply Thresholds", class_="btn-outline-secondary"),
                               TOOLTIPS["cy_rethreshold"]),
                    ui.tooltip(ui.input_action_button("cy_restyle", "Apply Style", class_="btn-outline-secondary"),
                               TOOLTIPS["cy_restyle"]),
                    class_="button-row",
                ),
            ),
            ui.card(
                ui.card_header("Cytoscape"),
                ui.accordion(
                    ui.accordion_panel(
                        "Status",
                        ui.output_ui("cy_status"),
                        ui.div(
                            ui.input_action_button("cy_probe", "Check Connection", class_="btn-sm btn-outline-secondary"),
                            ui.tooltip(ui.input_action_button("cy_unlock", "Unlock View", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_unlock"]),
                            class_="button-row",
                        ),
                        value="status",
                    ),
                    ui.accordion_panel(
                        ui.span("Selection ", ui.output_text("cy_selection_count", inline=True)),
                        ui.div(
                            ui.tooltip(ui.input_action_button("cy_read_selection", "Read Selection", class_="btn-sm btn-outline-primary"),
                                       TOOLTIPS["cy_read_selection"]),
                            ui.tooltip(ui.input_action_button("cy_hide_selected", "Hide Edges", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_hide_selected"]),
                            ui.tooltip(ui.input_action_button("cy_show_selected", "Show Edges", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_show_selected"]),
                            ui.tooltip(ui.input_action_button("cy_hide_unselected", "Hide Others", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_hide_unselected"]),
                            ui.tooltip(ui.input_action_button("cy_show_all", "Show All", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_show_all"]),
                            ui.tooltip(ui.input_action_button("cy_select_loners", "Select Loners", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_select_loners"]),
                            ui.tooltip(ui.input_action_button("cy_select_satellites", "Select Satellites", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_select_satellites"]),
                            class_="button-row mb-3",
                        ),
                        ui.output_data_frame("cy_selection_table"),
                        ui.div(
                            ui.tooltip(ui.download_button("cy_download_selection", "Download selection CSV", class_="btn-sm"),
                                       TOOLTIPS["cy_download_selection"]),
                            class_="button-row",
                        ),
                        value="selection",
                    ),
                    ui.accordion_panel(
                        "Select by relation",
                        ui.layout_columns(
                            ui.input_text("cy_rel_seed", tip("Seed", "cy_rel_seed"), placeholder="bait or protein"),
                            ui.input_select("cy_rel_kind", tip("Relation", "cy_rel_kind"),
                                            choices={"interactors": "Interactors of a bait",
                                                     "singletons": "Singletons of a bait (its only preys)",
                                                     "partners": "BioGRID or complex partners",
                                                     "cocomplex": "CORUM co-complex members"}),
                            col_widths=(5, 7),
                        ),
                        ui.layout_columns(
                            ui.input_numeric("cy_rel_saint", tip("Min SAINT", "cy_rel_cuts"), value=None, min=0, max=1, step=0.05),
                            ui.input_numeric("cy_rel_bfdr", "Max BFDR", value=None, min=0, max=1, step=0.01),
                            ui.input_numeric("cy_rel_abundance", "Min abundance", value=None, min=0),
                            ui.input_numeric("cy_rel_pubs", "Min publications", value=None, min=0, step=1),
                            col_widths=(3, 3, 3, 3),
                        ),
                        ui.div(
                            ui.tooltip(ui.input_action_button("cy_rel_replace", "Replace Selection", class_="btn-sm btn-outline-primary"),
                                       TOOLTIPS["cy_rel_replace"]),
                            ui.tooltip(ui.input_action_button("cy_rel_add", "Add to Selection", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_rel_add"]),
                            class_="button-row",
                        ),
                        value="relation",
                    ),
                    ui.accordion_panel(
                        "Cluster the selection",
                        ui.layout_columns(
                            ui.input_numeric("cy_cl_resolution", tip("Resolution", "cy_cl_resolution"), value=1.0, min=0.1, max=5, step=0.1),
                            ui.input_numeric("cy_cl_seed", tip("Seed", "cy_cl_seed"), value=17, min=0, step=1),
                            ui.input_numeric("cy_cl_lit_weight", tip("Reference weight", "cy_cl_lit_weight"), value=1.0, min=0, step=0.25),
                            col_widths=(4, 4, 4),
                        ),
                        ui.layout_columns(
                            ui.input_select("cy_cl_layout", tip("Cluster layout", "cy_cl_layout"),
                                            choices={"circles": "Circles", "force-directed": "Force-directed",
                                                     "kamada-kawai": "Kamada-Kawai"}),
                            ui.input_checkbox("cy_cl_recolor", tip("Recolor communities", "cy_cl_recolor"), value=True),
                            col_widths=(6, 6),
                        ),
                        ui.div(
                            ui.tooltip(ui.input_action_button("cy_cluster", "Cluster and Repack", class_="btn-sm btn-outline-primary"),
                                       TOOLTIPS["cy_cluster"]),
                            class_="button-row",
                        ),
                        value="cluster",
                    ),
                    ui.accordion_panel(
                        "Saved layouts",
                        ui.input_selectize("cy_layout_name", tip("Saved layout", "cy_layout_name"), choices=[],
                                           options={"create": True, "placeholder": "pick a saved layout or type a new name"}),
                        ui.div(
                            ui.tooltip(ui.input_action_button("cy_save_layout", "Save Layout", class_="btn-sm btn-outline-primary"),
                                       TOOLTIPS["cy_save_layout"]),
                            ui.tooltip(ui.input_action_button("cy_load_layout", "Load Layout", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_load_layout"]),
                            class_="button-row",
                        ),
                        value="layouts",
                    ),
                    ui.accordion_panel(
                        "Export",
                        ui.div(
                            ui.tooltip(ui.input_action_button("cy_export", "Export PNG", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["cy_export"]),
                            class_="button-row",
                        ),
                        value="export",
                    ),
                    open=["status", "selection"],
                ),
            ),
            col_widths=(4, 8),
        ),
        ui.card(
            ui.card_header("Activity"),
            ui.output_ui("cy_activity"),
        ),
        ),
    ),
    ui.nav_panel("Export",
                    ui.output_ui("empty_state_downloads"),
                    ui.panel_conditional(
                    "output.active_state !== 'none'",
                    ui.layout_columns(
                        ui.card(
                            ui.card_header("Download Builder"),
                            # The label carries the live sidebar values.
                            ui.tooltip(ui.input_switch("dl_apply_thresholds",
                                                       f"Apply thresholds ({describe(DEFAULT)})", value=False),
                                       TOOLTIPS["dl_apply_thresholds"]),
                            ui.input_radio_buttons("dl_preset", tip("Preset", "dl_preset"),
                                choices={key: preset.label for key, preset in dp.PRESETS.items()}),
                            ui.output_ui("dl_groups_ui"),
                            ui.panel_conditional("input.dl_preset === 'custom'",
                                ui.input_selectize("custom_columns", tip("Select Columns", "custom_columns"),
                                    choices=DEFAULT_CUSTOM_COLUMNS, multiple=True,
                                    selected=DEFAULT_CUSTOM_COLUMNS),
                            ),
                            ui.panel_conditional("input.dl_preset === 'genelist'",
                                ui.input_radio_buttons("dl_genelist_mode", tip("Gene list mode", "dl_genelist_mode"),
                                    choices={"pooled": "Pooled unique genes", "per_bait": "Per bait"}),
                            ),
                            ui.panel_conditional("input.dl_preset === 'prohits'",
                                ui.input_select("dl_prohits_abundance", tip("Abundance column", "dl_prohits_abundance"),
                                    choices=["AvePSM", "AvgIntensity"]),
                            ),
                        ),
                        ui.card(
                            ui.card_header("Preview"),
                            ui.output_text("download_row_count"),
                            ui.output_data_frame("dl_preview_table"),
                            ui.download_button("download_preset", "Download"),
                        ),
                        col_widths=(4,8)
                    ),
                    ui.card(
                        ui.card_header(tip("Batch Export", "batch_export")),
                        ui.p("Download all results for a dataset as a ZIP file. Includes merged.csv, annotated_scores.csv, Feature_enrichment.csv (if available), and SAINT input files."),
                        ui.download_button("download_batch", "Download All Results (ZIP)"),
                    ),
                    ),
    ),
    # The sidebar holds what every tab reads: the active dataset and the thresholds.
    sidebar=ui.sidebar(
        ui.input_select("active_dataset",
                        ui.span(tip("Active dataset", "active_dataset"), " ",
                                ui.output_ui("active_dataset_badge", inline=True)),
                        choices=[]),
        *threshold_sliders("threshold", DEFAULT),
        ui.output_ui("preset_buttons"),
        ui.output_ui("wdfdr_warning"),
        ui.accordion(
            ui.accordion_panel(
                "Agent access",
                *[ui.div(ui.p(f"Connect {agent} to this server with:", class_="small mb-1"),
                         ui.tags.pre(command, class_="mcp-command"))
                  for agent, command in MCP_ADD_COMMANDS.items()],
            ),
            ui.accordion_panel(
                "About and feedback",
                ui.p("Your feedback is invaluable in helping us improve the tool.", class_="small"),
                ui.tags.ul(
                    ui.tags.li(ui.a("Report a bug", href="https://github.com/plutzer/ProxiMate/issues",
                                    target="_blank")),
                    ui.tags.li(ui.a("GitHub repository", href="https://github.com/plutzer/ProxiMate",
                                    target="_blank")),
                    class_="small",
                ),
            ),
            open=False,
        ),
        # Resolved once: the build cannot change while the server is running, and this
        # is the identifier a bug report has to quote to be reproducible.
        ui.div(provenance.version_label(), class_="help-text"),
        width=300,
    ),
    title="ProxiMate",
    header=ui.TagList(ui.include_css(STYLESHEET), ui.tags.script(PAGE_SCRIPT),
                      ui.output_ui("agent_banner"),
                      # Read by the tabs' conditional panels; never shown.
                      ui.div(ui.output_text("active_state"), class_="d-none")),
)



def notify(message, type="message", duration=5, exc_info=False):
    """Show a Shiny notification and log the same text.

    Notifications are raised only through here.  A message shown to the user
    that leaves no trace on the server gives a later support request nothing to
    work from, and the toast itself is gone as soon as it is dismissed.
    """
    level = {"error": logging.ERROR, "warning": logging.WARNING}.get(type, logging.INFO)
    logger.log(level, "notification [%s]: %s", type,
               " | ".join(str(message).splitlines()), exc_info=exc_info)
    ui.notification_show(message, type=type, duration=duration)


def format_error_notification(error):
    """
    Format a ProxiMateError into a user-friendly notification message.
    Returns formatted string with message and suggestions.
    """
    if isinstance(error, ProxiMateError):
        msg_parts = [error.user_message]
        if error.suggestions:
            msg_parts.append("\n\nHow to fix:")
            for i, suggestion in enumerate(error.suggestions, 1):
                msg_parts.append(f"{i}. {suggestion}")
        return "\n".join(msg_parts)
    return str(error)


def server(input: Inputs, output: Outputs, session: Session):
    # The table lives in dataset_store, shared by every browser session and the MCP
    # server; the session learns about changes by polling its version counter.
    @reactive.poll(lambda: dataset_store.version(), 1.0)
    def datasets():
        return dataset_store.table()

    # Agent activity.  The MCP server shares this process, so its running jobs and
    # finished operations are polled and surfaced here: a banner while a job runs and a
    # notification when a dataset or Cytoscape operation starts, finishes or fails.
    @reactive.poll(lambda: backend.running_jobs(), 1.0)
    def agent_jobs():
        return {name: job for name, job in backend.running_jobs().items() if job['actor'] == 'mcp'}

    @reactive.poll(lambda: mcp_registry.last_seq(), 1.0)
    def agent_seq():
        return mcp_registry.last_seq()

    announced = {'jobs': set(), 'seq': mcp_registry.last_seq()}

    @reactive.effect
    def announce_agent_jobs():
        jobs = agent_jobs()
        for name in set(jobs) - announced['jobs']:
            notify(f"Agent started {jobs[name]['what']} of '{name}'.", type="message", duration=8)
        announced['jobs'] = set(jobs)

    @reactive.effect
    def announce_agent_operations():
        seq = agent_seq()
        for entry in mcp_registry.activity_since(announced['seq']):
            if entry['mode'] in ('dataset', 'cytoscape'):
                if entry['ok']:
                    notify(f"Agent finished {entry['op']} {entry['detail']}.".replace(' .', '.'),
                           type="message", duration=8)
                else:
                    notify(f"Agent {entry['op']} failed: {entry['detail']}", type="error", duration=None)
        announced['seq'] = seq

    @render.ui
    def agent_banner():
        jobs = agent_jobs()
        if not jobs:
            return None
        lines = [f"{job['what']} of '{name}' (since {job['since'][11:]})" for name, job in jobs.items()]
        return ui.div("Agent working: " + "; ".join(lines) + ". The dataset is locked until it finishes.",
                      class_="agent-banner")

    feature_enrichment = reactive.Value(pd.DataFrame())
    last_parse = reactive.Value(None)

    @render.data_frame
    def render_datasets():
        table = datasets()
        feature_enrichment()  # the State column follows a finished enrichment run
        enriched = [os.path.isfile(os.path.join(out_dir, str(name), 'Feature_enrichment.csv'))
                    for name in table['Dataset Name']]
        state = ['enriched' if e else 'scored' if scored == 'Yes' else 'parsed'
                 for e, scored in zip(enriched, table['Scored'])]
        shown = table.rename(columns={'Dataset Name': 'Dataset', 'Input Type': 'Input',
                                      'Quant Type': 'Quant'}).drop(columns=['Scored'])
        shown.insert(5, 'State', state)
        return render.DataGrid(shown, selection_mode="row")

    @render.ui
    def format_note():
        return ui.p(FORMAT_NOTES[input.input_format()], class_="help-text")

    @render.download(filename="experimental_design_template.csv")
    def ed_template():
        yield "Experiment Name,Type,Bait,Replicate,Bait ID,Group\n"
        yield "RUVBL1_1,T,RUVBL1,1,Q9Y265,\n"
        yield "GFP_1,C,GFP,1,None,\n"

    @render.ui
    def parse_summary():
        row = last_parse.get()
        if row is None or row['Dataset Name'] not in available_datasets():
            return None
        counts = (('Experiments', 'experiments'), ('Controls', 'controls'),
                  ('n_preys', 'preys'), ('n_interactions', 'interactions written'))
        return ui.div(
            ui.span(f"Parsed '{row['Dataset Name']}':", class_="fw-semibold"),
            *[ui.span(f"{row[key]} {label}", class_="badge text-bg-light") for key, label in counts],
            ui.tags.button("Score now", type="button", class_="btn btn-sm btn-primary ms-auto",
                           onclick="document.getElementById('score_card').scrollIntoView("
                                   "{behavior: 'smooth', block: 'start'}); "
                                   "document.getElementById('score_data').focus();"),
            class_="parse-summary")

    @reactive.calc
    def available_datasets():
        datasets()
        return dataset_store.names()

    @reactive.calc
    def scored_datasets():
        datasets()
        return dataset_store.scored_names()

    # The active dataset: one sidebar select that every tab reads.
    @reactive.calc
    def active():
        """The active dataset's name, or None while the select names nothing in the store."""
        name = input.active_dataset()
        return name if name and name in available_datasets() else None

    @reactive.calc
    def active_scored():
        """The active dataset's name when it is scored, else None."""
        name = active()
        return name if name in scored_datasets() else None

    @output(suspend_when_hidden=False)
    @render.text
    def active_state():
        """'none', 'parsed' or 'scored'; the tabs' conditional panels key on it."""
        if not active():
            return "none"
        return "scored" if active_scored() else "parsed"

    @render.ui
    def active_dataset_badge():
        if not active():
            return None
        if active_scored():
            return ui.span("scored", class_="badge bg-success")
        return ui.span("parsed", class_="badge bg-secondary")

    @reactive.effect
    @reactive.event(datasets)
    def update_active_dataset():
        names = available_datasets()
        current = input.active_dataset()
        if current not in names:
            current = names[0] if names else None
        ui.update_select("active_dataset", choices=names, selected=current)

    @reactive.effect
    @reactive.event(input.render_datasets_cell_selection)
    def select_active_from_table():
        rows = render_datasets.cell_selection()["rows"]
        if rows:
            ui.update_select("active_dataset", selected=str(datasets().iloc[rows[0]]['Dataset Name']))

    # The shared thresholds.  Sliders are addressed by prefix so the Compare Baits
    # tab's two private sets reuse the same helpers.
    def read_thresholds(prefix):
        return {key: getattr(input, f"{prefix}_{suffix}")() for suffix, key, *_ in SLIDERS}

    def set_thresholds(prefix, values):
        for suffix, key, *_ in SLIDERS:
            ui.update_slider(f"{prefix}_{suffix}", value=values[key])

    @reactive.calc
    def thresholds():
        return read_thresholds("threshold")

    # The row re-renders whenever the sliders move, so the highlighted button follows
    # them.  A click sets ``input.preset`` as an event rather than counting on an
    # action button: a re-rendered action button restarts its count at 0, so its next
    # click would repeat a value the server has already seen and be dropped.
    @render.ui
    def preset_buttons():
        current = active_preset(thresholds())
        return ui.div(
            *[ui.tooltip(ui.tags.button(
                  name, id=preset_id(name), type="button",
                  class_="btn btn-sm " + ("btn-primary" if name == current else "btn-outline-secondary"),
                  onclick=f"Shiny.setInputValue('preset', {name!r}, {{priority: 'event'}})"),
                  TOOLTIPS[preset_id(name)], placement="bottom")  # not over the other buttons
              for name in PRESETS],
            class_="preset-row mb-2")

    @reactive.effect
    @reactive.event(input.preset)
    def apply_preset():
        set_thresholds("threshold", PRESETS[input.preset()])

    @render.ui
    def wdfdr_warning():
        name = active_scored()
        if not name or thresholds()['WDFDR'] >= 1.0:
            return None
        if float(dataset_store.row(name)['WDFDR iterations']) > 0:
            return None
        return ui.div(f"⚠ '{name}' was scored with 0 WDFDR iterations, so it has no WDFDR values "
                      "and nothing passes a WDFDR cut below 1.", class_="wdfdr-warning")

    def empty_state(purpose, needs_scores=True):
        """What a tab shows in place of its controls until a dataset that fits is active."""
        if not available_datasets():
            return ui.div(ui.h4("No datasets in this session"),
                          ui.p(f"Parse and score a dataset on the Data & Scoring tab to {purpose}."),
                          class_="empty-state")
        if needs_scores and active() and not active_scored():
            return ui.div(ui.p(f"'{active()}' is parsed but not scored. "
                               "Score it on the Data & Scoring tab."), class_="empty-state")
        return None

    @render.ui
    def empty_state_thresholding():
        return empty_state("view quality metrics")

    @render.ui
    def empty_state_feature_analysis():
        return empty_state("test its networks for enriched protein features")

    @render.ui
    def empty_state_network_comparison():
        return empty_state("compare the networks of two baits")

    @render.ui
    def empty_state_cytoscape():
        return empty_state("send its network to Cytoscape")

    @render.ui
    def empty_state_downloads():
        return empty_state("download its results", needs_scores=False)

    def plot_placeholder(text):
        """What a plot card shows until it has something to draw."""
        return ui.p(text, class_="plot-placeholder")

    def dataset_baits(name):
        """The baits of a scored dataset, in file order."""
        scores = pd.read_csv(backend.results_path(name), usecols=['Experiment.ID'])
        return scores['Experiment.ID'].astype(str).unique().tolist()

    saint_baits = reactive.Value(EMPTY_BAITS)
    ed_dataframe = reactive.Value(EMPTY_ED)

    @render.data_frame
    def bait_table():
        # Return the bait table for SAINT input
        return render.DataGrid(saint_baits.get(), editable=True)

    @render.data_frame
    def ed_table_mq():
        # Return the ED table for MaxQuant input
        return render.DataGrid(ed_dataframe.get(), editable=True)

    @render.data_frame
    def ed_table_diann():
        # Return the ED table for DIA-NN input
        return render.DataGrid(ed_dataframe.get(), editable=True)

    @render.data_frame
    def ed_table_pioneer():
        return render.DataGrid(ed_dataframe.get(), editable=True)

    @render.data_frame
    def ed_table_fragpipe():
        # Return the ED table for FragPipe input
        return render.DataGrid(ed_dataframe.get(), editable=True)

    @render.data_frame
    def ed_table_msstats():
        # Return the ED table for MSstats input
        return render.DataGrid(ed_dataframe.get(), editable=True)

    @reactive.effect
    @reactive.event(input.bait)
    def update_bait_table():
        # Check to see if the bait file has been uploaded
        if input.bait.get():
            # Read the bait file and set it to the saint_baits reactive value
            # Read the bait file into a DataFrame

            # Catch bad input files here

            # Check to make sure the columns are correct and don't have any missing values
            bait_path = input.bait.get()[0]['datapath']
            try:
                baits = pd.read_csv(bait_path, sep="\t", header=None, index_col=None,
                                    names=["Experiment Name", "Bait", "Type"])
                # Set the column before publishing: mutating the frame afterwards
                # through .get() changes it without notifying dependants.
                baits['Bait ID'] = 'None'  # Default value for Bait ID
                saint_baits.set(baits)
                logger.info("Read SAINT bait file: %d experiments", len(baits))
            except Exception as e:
                notify(f"Could not read the bait file: {e}\n\n"
                       "Expected a tab-separated file with columns: "
                       "Experiment Name, Bait, Type.",
                       type="error", duration=10, exc_info=True)

    @reactive.effect
    @reactive.event(input.ed_file)
    def update_ed_table():
        # Check to see if the ED file has been uploaded
        if input.ed_file.get():
            # Read the ED file and set it to the ed_dataframe reactive value
            try:
                # Read the ED file into a DataFrame
                ed_df = pd.read_csv(input.ed_file.get()[0]['datapath'])

                # Validate required columns
                required_cols = ["Experiment Name", "Type", "Bait", "Replicate"]
                missing_cols = [col for col in required_cols if col not in ed_df.columns]
                if missing_cols:
                    notify(
                        f"ED file is missing required columns: {', '.join(missing_cols)}",
                        type="error"
                    )
                    return

                # Add Bait ID column if not present
                if "Bait ID" not in ed_df.columns:
                    ed_df['Bait ID'] = 'None'

                # Add Group column if not present (optional, for paired controls)
                if "Group" not in ed_df.columns:
                    ed_df['Group'] = ''

                # Validate Type values
                invalid_types = ed_df[~ed_df['Type'].isin(['C', 'T'])]
                if len(invalid_types) > 0:
                    notify(
                        f"ED file contains invalid Type values. Must be 'C' or 'T'.",
                        type="error"
                    )
                    return

                ed_dataframe.set(ed_df)
            except Exception as e:
                notify(
                    f"Error reading ED file: {str(e)}",
                    type="error",
                    exc_info=True
                )

    @reactive.effect
    @reactive.event(input.input_format)
    def clear_tables_on_format_change():
        # Clear tables when format changes to avoid showing stale data
        if input.input_format.get() == "SAINT":
            ed_dataframe.set(EMPTY_ED)
        else:
            saint_baits.set(EMPTY_BAITS)

    # Long jobs (parse, score, feature enrichment, Cytoscape send) run in a worker
    # thread so the session keeps answering.  The thread writes its stage and outcome
    # into ``job``; the session polls it to drive the progress modal and, on success,
    # runs the completion callback on its own side.
    job = {'lock': threading.Lock(), 'version': 0, 'running': False, 'title': '', 'dataset': None,
           'stages': (), 'stage': None, 'started': 0.0, 'result': None, 'error': None,
           'on_success': None}

    def start_job(title, dataset, stages, work, on_success):
        """Run ``work(progress)`` in a thread under the progress modal; ``on_success``
        receives its return value.  ``stages`` are the messages ``work`` reports."""
        with job['lock']:
            busy = job['title'] if job['running'] else None
            if not busy:
                job.update(running=True, title=title, dataset=dataset, stages=tuple(stages),
                           stage=None, started=time.time(), result=None, error=None,
                           on_success=on_success, version=job['version'] + 1)
        if busy:
            notify(f"Wait for '{busy}' to finish first.", type="error")
            return

        def progress(message, value):
            with job['lock']:
                job.update(stage=message, version=job['version'] + 1)

        def run():
            try:
                outcome = {'result': work(progress)}
            except Exception as e:
                logger.exception("%s failed", title)
                outcome = {'error': e}
            with job['lock']:
                job.update(running=False, version=job['version'] + 1, **outcome)

        threading.Thread(target=run, daemon=True, name=title).start()
        ui.modal_show(ui.modal(ui.output_ui("job_modal"), title=title, footer=None,
                               easy_close=False, size="l"))

    @reactive.poll(lambda: (job['version'], int(time.time()) if job['running'] else 0), 1.0)
    def job_state():
        with job['lock']:
            return {key: value for key, value in job.items() if key != 'lock'}

    @render.ui
    def job_modal():
        state = job_state()
        reached = state['stages'].index(state['stage']) if state['stage'] in state['stages'] else -1
        marks = {'done': '✓', 'current': '●', 'failed': '✗', 'pending': '○'}
        items = []
        for i, stage in enumerate(state['stages']):
            status = ('done' if i < reached else 'pending' if i > reached
                      else 'failed' if state['error'] is not None else 'current')
            items.append(ui.tags.li(f"{marks[status]} {stage}", class_=status))
        elapsed = int(time.time() - state['started'])
        parts = [ui.tags.ul(*items, class_="job-stages"),
                 ui.p(f"Elapsed {elapsed // 60}:{elapsed % 60:02d}", class_="help-text")]
        if state['error'] is not None:
            parts.append(ui.div(format_error_notification(state['error']), class_="job-error"))
        if state['dataset'] and os.path.isfile(os.path.join(out_dir, state['dataset'], 'proximate.log')):
            lines = backend.log_tail(state['dataset'], 15)['lines']
            parts.append(ui.tags.pre("\n".join(lines), class_="job-log"))
            parts.append(ui.tags.script(
                "document.querySelectorAll('.job-log').forEach(el => el.scrollTop = el.scrollHeight);"))
        if not state['running']:
            parts.append(ui.div(ui.modal_button("Close"), class_="text-end"))
        return ui.div(*parts)

    @reactive.effect
    def finish_job():
        state = job_state()
        if state['running'] or state['on_success'] is None:
            return
        with job['lock']:
            job['on_success'] = None
        if state['error'] is None:
            ui.modal_remove()
            state['on_success'](state['result'])

    @reactive.effect
    @reactive.event(input.parse_data)
    def parse_data():
        dataset_name = input.dataset_name.get()
        input_format = input.input_format.get()

        def uploaded(field):
            files = getattr(input, field).get()
            return files[0]['datapath'] if files else None

        grids = {'MaxQuant': ed_table_mq, 'DIA-NN': ed_table_diann, 'Pioneer': ed_table_pioneer,
                 'FragPipe': ed_table_fragpipe, 'MSstats': ed_table_msstats}
        if input_format == 'SAINT':
            files = {'bait': uploaded('bait'), 'prey': uploaded('prey'),
                     'interaction': uploaded('interaction')}
            bait_frame, ed_frame = bait_table.data_view(), None
        else:
            main = {'MaxQuant': 'pg_file', 'DIA-NN': 'diann_matrix_file', 'Pioneer': 'pioneer_matrix_file',
                    'FragPipe': 'fragpipe_file', 'MSstats': 'msstats_file'}[input_format]
            key = backend.INPUT_FORMATS[input_format][0][0]
            files = {key: uploaded(main), 'ed': uploaded('ed_file')}
            bait_frame, ed_frame = None, grids[input_format].data_view()

        quant_type = input.quant_type.get()

        def work(progress):
            return backend.run_parse(dataset_name, input_format, files, quant_type, actor='gui',
                                     ed_frame=ed_frame, bait_frame=bait_frame, progress=progress)

        def done(row):
            notify(f"Parsed '{dataset_name}': {row['Experiments']} experiments, {row['Controls']} "
                   f"controls, {row['n_preys']} preys, {row['n_interactions']} interactions.")
            last_parse.set(row)
            # Make the new dataset active. The choices are passed here too because the
            # polled refresh of the select has not run yet.
            ui.update_select("active_dataset", choices=dataset_store.names(), selected=dataset_name)
            # A file input cannot be reset from the server; the name and grids can.
            ui.update_text("dataset_name", value="")
            saint_baits.set(EMPTY_BAITS)
            ed_dataframe.set(EMPTY_ED)

        start_job(f"Parsing {dataset_name or 'dataset'}", dataset_name, backend.PARSE_STAGES, work, done)

    # Delete and rename confirm in a modal.  The confirm buttons set their input as an
    # event (see ``preset_buttons``), since each modal re-creates them.
    @reactive.effect
    @reactive.event(input.delete_dataset)
    def ask_delete():
        name = active()
        if not name:
            notify("No dataset is active.", type="error")
            return
        ui.modal_show(ui.modal(
            ui.p(f"Delete '{name}' and every result in its directory? This cannot be undone."),
            title="Delete dataset", easy_close=True,
            footer=ui.div(
                ui.modal_button("Cancel"),
                ui.tags.button("Delete", type="button", class_="btn btn-danger ms-2",
                               onclick="Shiny.setInputValue('confirm_delete', 1, {priority: 'event'})"))))

    @reactive.effect
    @reactive.event(input.confirm_delete)
    def confirm_delete():
        name = active()
        ui.modal_remove()
        try:
            backend.delete_dataset(name, actor='gui')
            notify(f"Deleted '{name}'.")
        except Exception as e:
            notify(f"Could not delete '{name}': {e}", type="error", duration=None, exc_info=True)

    @reactive.effect
    @reactive.event(input.rename_dataset)
    def ask_rename():
        name = active()
        if not name:
            notify("No dataset is active.", type="error")
            return
        ui.modal_show(ui.modal(
            ui.input_text("rename_to", f"New name for '{name}'", value=name),
            title="Rename dataset", easy_close=True,
            footer=ui.div(
                ui.modal_button("Cancel"),
                ui.tags.button("Rename", type="button", class_="btn btn-primary ms-2",
                               onclick="Shiny.setInputValue('confirm_rename', 1, {priority: 'event'})"))))

    @reactive.effect
    @reactive.event(input.confirm_rename)
    def confirm_rename():
        old, new = active(), input.rename_to.get().strip()
        ui.modal_remove()
        try:
            backend.rename_dataset(old, new, actor='gui')
            ui.update_select("active_dataset", choices=dataset_store.names(), selected=new)
            notify(f"Renamed '{old}' to '{new}'.")
        except Exception as e:
            notify(f"Could not rename '{old}': {e}", type="error", duration=None, exc_info=True)

    @reactive.effect
    @reactive.event(input.clear_datasets)
    def clear_datasets():
        try:
            backend.clear_datasets()
        except backend.BusyError as e:
            notify(f"Cannot clear the session: {e}", type="error", duration=None)

    @render.download_button(
        filename=lambda: f"ProxiMateSession_{datetime.datetime.now().strftime('%Y%m%d')}.zip")
    def download_session():
        try:
            # One fixed path, overwritten per download, so archives do not accumulate.
            zip_path = os.path.join(tempfile.gettempdir(), "ProxiMateSession.zip")
            session_archive.write_session_archive(
                out_dir, dataset_store.names(), zip_path)
            logger.info("Session archive written to %s", zip_path)
            return zip_path
        except Exception as e:
            notify(f"Could not build the session archive: {e}",
                   type="error", duration=None, exc_info=True)
            return None

    @reactive.effect
    @reactive.event(input.upload_session)
    def upload_session():
        uploaded = input.session_file.get()
        if not uploaded:
            notify("No session archive selected.", type="error")
            return
        zip_path = uploaded[0]['datapath']

        try:
            table = backend.load_session(zip_path, actor='gui')
            notify(f"Session restored: {len(table)} dataset(s).", type="message")
        except (session_archive.SessionArchiveError, backend.BusyError) as e:
            notify(str(e), type="error", duration=None)
        except Exception as e:
            logger.exception("Failed to restore session from %s", zip_path)
            notify(f"Could not restore the session archive: {e}", type="error",
                   duration=None)


    @render.text
    def score_header():
        if not active():
            return "Scoring: no dataset in the session yet."
        return f"Scoring: {active()} with SAINT and CompPASS, then annotating it."

    @reactive.effect
    @reactive.event(input.active_dataset, datasets)
    def update_imputation_choices():
        dataset = active()
        if not dataset:
            return
        choices = backend.imputation_choices(dataset_store.row(dataset)['Quant Type'])
        selected = input.imputation_method.get()
        if int(selected) not in choices:
            selected = "0"
        ui.update_radio_buttons("imputation_method", choices=choices, selected=selected)

    @reactive.effect
    @reactive.event(input.active_dataset, input.imputation_method, input.pi_method)
    def update_pi_bait_choices():
        """Populate pi_bait dropdown from the active dataset's saved ED,
        filtered to control baits with >= 3 replicates."""
        if str(input.imputation_method.get()) != "2":
            return
        if input.pi_method.get() != "single_bait":
            return
        dataset = active()
        if not dataset:
            ui.update_select("pi_bait", choices=[])
            return
        ed_path = os.path.join(out_dir, dataset, "ED.csv")
        if not os.path.exists(ed_path):
            ui.update_select("pi_bait", choices=[])
            return
        try:
            ed = pd.read_csv(ed_path)
            controls = ed[ed['Type'] == 'C']
            counts = controls.groupby('Bait').size()
            eligible = sorted(counts[counts >= 3].index.tolist())
            ui.update_select("pi_bait", choices=eligible,
                             selected=eligible[0] if eligible else None)
        except Exception:
            # An empty dropdown here is indistinguishable from "no eligible
            # controls", so the reason has to be recorded somewhere.
            logger.exception("Could not read control baits from %s", ed_path)
            ui.update_select("pi_bait", choices=[])

    # Scoring data
    @reactive.effect
    @reactive.event(input.score_data)
    def score_data():
        dataset_name = active()
        if not dataset_name:
            notify("Parse a dataset first.", type="error")
            return
        imputation = int(input.imputation_method.get())
        pi_method = input.pi_method.get() if imputation == 2 else None
        settings = dict(imputation=imputation, wdfdr_iterations=input.wdfdr_iterations.get(),
                        organism=input.organism.get(), exclude_hcm=input.exclude_hcm.get(),
                        pi_method=pi_method,
                        pi_bait=input.pi_bait.get() if pi_method == 'single_bait' else None,
                        aft_min_obs=input.aft_min_obs.get() if imputation in (2, 3) else None)

        def work(progress):
            return backend.run_score(dataset_name, actor='gui', progress=progress, **settings)

        def done(result):
            notify(f"Scored and annotated '{dataset_name}' ({result['organism']}, "
                   f"{result['imputation']} imputation, {result['wdfdr_iterations']} WDFDR iterations).")

        start_job(f"Scoring {dataset_name}", dataset_name, backend.SCORE_STAGES, work, done)

    # Quality Control tab
    @reactive.Calc
    def pca_matrix_cached():
        # Shared by the experiment and prey PCA plots and their PNG exports, so
        # the preprocessing runs once per settings change and all four agree.
        dataset_name = active()
        if not dataset_name:
            return None
        interaction_path = os.path.join(out_dir, dataset_name, "interaction.txt")
        if not os.path.exists(interaction_path):
            return None
        return prepare_pca_matrix(
            interaction_path,
            min_detection_frac=input.pca_min_detection(),
            imputation=input.pca_imputation(),
            normalization=input.pca_normalization(),
        )

    def _prey_pca_color_data(dataset_name, mode, bait, prey_index):
        """Color data for the prey PCA: (values, label, mode, threshold) as
        prey_pca_plot expects. Annotation-based modes fall back to uncolored
        when the dataset has no annotated_scores.csv yet. SaintScore coloring
        greys out preys below 0.1 so the color scale is spent on candidate
        interactors."""
        if mode == "detection":
            interaction_path = os.path.join(out_dir, dataset_name, "interaction.txt")
            counts = detection_counts(interaction_path).reindex(prey_index)
            return counts, "Detections", "continuous", None
        if mode in ("saint", "hpa", "scl"):
            scores_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")
            if not os.path.exists(scores_path):
                return None, None, "none", None
            scores = pd.read_csv(scores_path)
            if mode == "saint":
                if not bait:
                    return None, None, "none", None
                s = (scores[scores['Experiment.ID'] == bait]
                     .set_index('Prey.ID')['SaintScore'])
                # A prey never scored for this bait is a non-interactor: score 0
                return s.reindex(prey_index).fillna(0.0), f"SaintScore ({bait})", "continuous", 0.1
            column = "Main location" if mode == "hpa" else "first_SCL"
            if column not in scores.columns:
                return None, None, "none", None
            s = (scores.drop_duplicates('Prey.ID')
                 .set_index('Prey.ID')[column])
            label = "HPA Main location" if mode == "hpa" else "UniProt localization"
            return reduce_categorical(s.reindex(prey_index)), label, "categorical", None
        return None, None, "none", None

    @render_plotly
    def raw_pca_plot():
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Generating Plots...", value=25)
            dataset_name = active()
            if not dataset_name:
                return None

            interaction_path = os.path.join(out_dir, dataset_name, "interaction.txt")
            ed_path = os.path.join(out_dir, dataset_name, "ED.csv")
            if not (os.path.exists(interaction_path) and os.path.exists(ed_path)):
                return None

            matrix = pca_matrix_cached()
            if matrix is None:
                return None
            show_legend = input.pca_show_legend() if many_pca_baits() else True
            return pca_plot(interaction_path, ed_path, matrix=matrix, show_legend=show_legend)

    @reactive.calc
    def many_pca_baits():
        """Whether the active dataset has more test baits than the PCA legend shows by default."""
        ed = pd.read_csv(os.path.join(out_dir, active(), "ED.csv"))
        return ed.loc[ed['Type'] != 'C', 'Bait'].nunique() > PCA_LEGEND_MAX_BAITS

    @render.ui
    def pca_legend_toggle():
        if not active() or not many_pca_baits():
            return None
        return ui.input_checkbox("pca_show_legend", "Show legend", value=False)

    @render_plotly
    def prey_pca_plot():
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Generating Plots...", value=25)
            dataset_name = active()
            if not dataset_name:
                return None
            matrix = pca_matrix_cached()
            if matrix is None:
                return None
            values, label, mode, threshold = _prey_pca_color_data(
                dataset_name, input.prey_pca_color(), input.prey_pca_bait(),
                matrix.index)
            gene_names = prey_gene_names(os.path.join(out_dir, dataset_name, "prey.txt"))
            return plot_prey_pca(matrix, color_values=values,
                                 color_label=label, color_mode=mode,
                                 color_threshold=threshold, gene_names=gene_names)


    def qc_experiments():
        """The baits the QC metrics and hidden plots cover: None for all, else one."""
        return None if input.qc_bait.get() == "All" else [input.qc_bait.get()]

    @render_widget
    def known_retention_plot():
        if not active_scored():
            return None
        return saint_known_retention(backend.results_path(active_scored()), ctrl_experiments=qc_experiments())

    @render_widget
    def roc_curve_plot():
        if not active_scored():
            return None
        return roc_plot(backend.results_path(active_scored()), known_type=input.roc_known_type.get(),
                        ctrl_experiments=qc_experiments())

    @render.ui
    def metric_bait_label():
        bait_selection = input.qc_bait.get()
        if bait_selection == "All":
            label_text = "Metrics for All Baits"
        else:
            label_text = f"Metrics for {bait_selection}"
        return ui.p(label_text, class_="section-heading")

    @reactive.calc
    def qc_metrics():
        """The three threshold metrics for the active scored dataset, or None."""
        if not active_scored():
            return None
        return calculate_threshold_metrics(backend.results_path(active_scored()), thresholds(),
                                           ctrl_experiments=qc_experiments())

    @render.ui
    def metric_network_size():
        metrics = qc_metrics()
        if metrics is None:
            return ui.value_box(tip("Median Network Size", "metric_network_size"), "No data", showcase=None)
        return ui.value_box("Median Network Size", f"{metrics['median_network_size']:.1f}",
                            showcase=None, theme="primary")

    @render.ui
    def metric_enrichment():
        metrics = qc_metrics()
        if metrics is None:
            return ui.value_box(tip("Known Enrichment", "metric_enrichment"), "No data", showcase=None)
        return ui.value_box("Known Enrichment", f"{metrics['enrichment_ratio']:.2f}x",
                            showcase=None, theme="success")

    @render.ui
    def metric_degree():
        metrics = qc_metrics()
        if metrics is None:
            return ui.value_box(tip("Mean Prey-Prey Degree", "metric_degree"), "No data", showcase=None)
        # A degree of zero is a real result, so an unreadable BioGRID summary is named
        # rather than averaged into one.
        if metrics['mean_degree'] is None:
            return ui.value_box(
                "Mean Prey-Prey Degree",
                "No reference set",
                ui.p("No BioGRID summary for this dataset's organism.",
                     class_="small mb-0"),
                showcase=None
            )
        return ui.value_box("Mean Prey-Prey Degree", f"{metrics['mean_degree']:.1f}",
                            showcase=None, theme="info")

    @reactive.calc
    def qc_results_frame():
        """The QC tab's passing interactions for the QC bait selection, as displayed."""
        rows = backend.passing_scores(active_scored(), thresholds(), qc_experiments())
        return backend.results_table(rows)

    @render.text
    def qc_results_count():
        if not active_scored():
            return ""
        return f"{len(qc_results_frame())} interactions pass the sidebar thresholds"

    @render.data_frame
    def qc_results():
        if not active_scored():
            return None
        return render.DataGrid(qc_results_frame(), filters=True, width="100%", height="400px")

    @render.download_button(filename=lambda: f"{active_scored()}_passing_interactions.csv")
    def download_qc_results():
        if not active_scored():
            notify("Score a dataset first.", type="error")
            return None
        path = os.path.join(out_dir, "passing_interactions.csv")
        qc_results_frame().to_csv(path, index=False)
        return path

    # Plot export download handlers (using matplotlib for PNG export)
    @render.download_button(filename="pca_plot.png")
    def download_pca_plot():
        dataset_name = active()
        if not dataset_name:
            notify("No dataset selected.", type="error")
            return None
        interaction_path = os.path.join(out_dir, dataset_name, "interaction.txt")
        ed_path = os.path.join(out_dir, dataset_name, "ED.csv")
        if not (os.path.exists(interaction_path) and os.path.exists(ed_path)):
            notify("Required files not found.", type="error")
            return None
        fig = pca_plot_matplotlib(interaction_path, ed_path, matrix=pca_matrix_cached())
        filepath = os.path.join(out_dir, "pca_plot.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.download_button(filename="prey_pca_plot.png")
    def download_prey_pca_plot():
        dataset_name = active()
        if not dataset_name:
            notify("No dataset selected.", type="error")
            return None
        matrix = pca_matrix_cached()
        if matrix is None:
            notify("Required files not found.", type="error")
            return None
        values, label, mode, threshold = _prey_pca_color_data(
            dataset_name, input.prey_pca_color(), input.prey_pca_bait(),
            matrix.index)
        fig = prey_pca_matplotlib(matrix, color_values=values,
                                  color_label=label, color_mode=mode,
                                  color_threshold=threshold)
        filepath = os.path.join(out_dir, "prey_pca_plot.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.download_button(filename="scatter_plot.png")
    def download_scatter_plot():
        dataset_name = active_scored()
        bait_selection = input.qc_bait.get()
        if not dataset_name or bait_selection == "All":
            notify("Select a specific bait to export the scatter plot.", type="error")
            return None
        fig = saint_scatter_matplotlib(backend.results_path(dataset_name), bait_selection,
                                       thresholds()['SaintScore'])
        filepath = os.path.join(out_dir, "scatter_plot.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.download_button(filename="heatmap.png")
    def download_heatmap():
        if isinstance(feature_heatmap(), str):
            notify(feature_heatmap(), type="error")
            return None
        # A fresh figure: the on-screen one has been resized to its card.
        feature_data = pd.read_csv(os.path.join(out_dir, active_scored(), "Feature_enrichment.csv"))
        heatmap = plot_results(feature_data, input.feature_type(), input.num_features() or 30)
        filepath = os.path.join(out_dir, "heatmap.png")
        heatmap.savefig(filepath, dpi=150, bbox_inches='tight')
        return filepath

    @render.download_button(filename="volcano_plot.png")
    def download_volcano_plot():
        bait_a = input.comp_bait_a.get()
        bait_b = input.comp_bait_b.get()
        if not all([active_scored(), bait_a, bait_b]):
            notify("Select baits first.", type="error")
            return None
        # Use cached volcano data
        volcano_data = comp_volcano_data_cached()
        if volcano_data.empty:
            notify("No data available for volcano plot.", type="error")
            return None
        # Use matplotlib version for export
        fig = create_volcano_plot_matplotlib(volcano_data, bait_a, bait_b)
        filepath = os.path.join(out_dir, "volcano_plot.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.download_button(filename="venn_diagram.png")
    def download_venn_diagram():
        bait_a = input.comp_bait_a.get()
        bait_b = input.comp_bait_b.get()
        if not all([bait_a, bait_b]):
            notify("Select baits first.", type="error")
            return None
        # Use cached filtered data to create sets
        data_a = comp_filtered_data_a_cached()
        data_b = comp_filtered_data_b_cached()
        set_a = set(data_a['Prey.ID'].unique()) if len(data_a) > 0 else set()
        set_b = set(data_b['Prey.ID'].unique()) if len(data_b) > 0 else set()
        if len(set_a) == 0 and len(set_b) == 0:
            notify("No data available for Venn diagram.", type="error")
            return None
        # Use matplotlib version for clean export
        fig = create_venn_diagram_matplotlib(set_a, set_b, bait_a, bait_b)
        filepath = os.path.join(out_dir, "venn_diagram.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.ui
    def scatter_body():
        if input.qc_bait() == "All":
            return plot_placeholder("Select one bait above to see its SAINT score against fold change.")
        return ui.div(output_widget("saint_scatter_plot"),
                      ui.p("Dashed line: SAINT threshold from the sidebar.", class_="help-text"))

    @render_widget
    def saint_scatter_plot():
        dataset_name = active_scored()
        bait_selection = input.qc_bait()
        if not dataset_name or bait_selection == "All":
            return None
        return plot_saint_scatter(backend.results_path(dataset_name), bait_selection,
                                  thresholds()['SaintScore'])

    # Feature Enrichment tab
    @reactive.effect
    @reactive.event(input.feature_analysis)
    def feature_analysis():
        dataset_name = active_scored()
        if not dataset_name:
            notify("Score a dataset first.", type="error")
            return
        cutoffs = thresholds()

        def work(progress):
            return backend.write_feature_enrichment(dataset_name, cutoffs, actor='gui', progress=progress)

        def done(result):
            feature_enrichment.set(result)
            notify(f"Feature enrichment for '{dataset_name}': {len(result)} rows written to "
                   "Feature_enrichment.csv")

        start_job(f"Feature enrichment for {dataset_name}", dataset_name, backend.FEATURE_STAGES,
                  work, done)

    @reactive.calc
    def feature_heatmap():
        """The feature heatmap for the active dataset, or the message shown in its place."""
        _ = feature_enrichment.get()
        dataset = active_scored()
        if not dataset or not os.path.exists(os.path.join(out_dir, dataset, "Feature_enrichment.csv")):
            return "Run Feature Enrichment to see the top features."
        feature_data = pd.read_csv(os.path.join(out_dir, dataset, "Feature_enrichment.csv"))
        feature_type = input.feature_type()
        try:
            return plot_results(feature_data, feature_type, input.num_features() or 30)
        except ValueError as e:
            # Raised when too few features pass to draw a map
            logger.warning("Cannot plot %s enrichment for '%s': %s", feature_type, dataset, e)
            return (f"No {feature_type} features pass adjusted p ≤ 0.05 and enrichment ≥ 2. "
                    "Try another feature type or looser thresholds.")

    @render.ui
    def feature_heatmap_body():
        heatmap = feature_heatmap()
        if isinstance(heatmap, str):
            return plot_placeholder(heatmap)
        # 96 px per inch, so the card is as tall as the figure was laid out for.
        return ui.output_plot("feature_enrichment_plot", height=f"{heatmap.figure.get_figheight() * 96:.0f}px")

    @render.plot
    def feature_enrichment_plot():
        heatmap = feature_heatmap()
        return None if isinstance(heatmap, str) else heatmap

    @reactive.effect
    @reactive.event(input.active_dataset, datasets, feature_enrichment)
    def update_download_bait_filter():
        """Update the bait dropdown for enrichment download based on available data."""
        dataset = active_scored()
        if not dataset:
            ui.update_select("download_bait_filter", choices=["All"])
            return

        feature_file = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(feature_file):
            ui.update_select("download_bait_filter", choices=["All"])
            return

        try:
            feature_data = pd.read_csv(feature_file)
            baits = feature_data['Bait'].unique().tolist()
            baits.insert(0, "All")
            ui.update_select("download_bait_filter", choices=baits)
        except Exception:
            logger.exception("Could not read baits from %s", feature_file)
            ui.update_select("download_bait_filter", choices=["All"])

    @render.download_button()
    def download_enrichment():
        """Download filtered enrichment results."""
        dataset = active_scored()
        if not dataset:
            notify("No dataset selected.", type="error")
            return None

        feature_file = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(feature_file):
            notify("No enrichment results available. Please run feature analysis first.", type="error")
            return None

        # Load enrichment data
        feature_data = pd.read_csv(feature_file)

        # Apply filters
        filtered_data = feature_data.copy()

        # Filter by feature type
        feature_type_filter = input.download_feature_type.get()
        if feature_type_filter != "All":
            filtered_data = filtered_data[filtered_data['Feature_type'] == feature_type_filter]

        # Filter by bait
        bait_filter = input.download_bait_filter.get()
        if bait_filter != "All":
            filtered_data = filtered_data[filtered_data['Bait'] == bait_filter]

        # Filter by adjusted p-value
        pvalue_threshold = input.download_pvalue_threshold.get()
        filtered_data = filtered_data[filtered_data['adj_p'] <= pvalue_threshold]

        # Filter by enrichment
        enrichment_threshold = input.download_enrichment_threshold.get()
        filtered_data = filtered_data[filtered_data['enrichment'] >= enrichment_threshold]

        if filtered_data.empty:
            notify("No results match the current filters. Try adjusting the thresholds.", type="warning")
            return None

        # The feature names are UniProt's GO terms with their ids trimmed; put the ids back.
        go_ids = go_term_ids(pd.read_csv(backend.results_path(dataset), usecols=list(GO_SOURCE_COLUMNS.values())))
        filtered_data.insert(2, 'GO_id', [go_ids.get(key) for key in
                                          zip(filtered_data['Feature_type'], filtered_data['Feature'])])

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"enrichment_{dataset}_{timestamp}.csv"
        savepath = os.path.join(out_dir, filename)
        filtered_data.to_csv(savepath, index=False)
        return savepath

    @reactive.effect
    @reactive.event(input.active_dataset, datasets)
    def update_qc_bait():
        baits = dataset_baits(active_scored()) if active_scored() else []
        ui.update_select("qc_bait", choices=["All"] + baits)

    @reactive.effect
    @reactive.event(input.active_dataset, datasets)
    def update_prey_pca_color_choices():
        # Annotation-based color options exist only once the dataset is scored;
        # HPA localization only when the organism has HPA data (human).
        choices = {"none": "None", "detection": "Detection count"}
        baits = []
        if active_scored():
            scores_path = backend.results_path(active_scored())
            cols = pd.read_csv(scores_path, nrows=0).columns
            choices["saint"] = "SaintScore (per bait)"
            if "Main location" in cols:
                choices["hpa"] = "HPA Main location"
            if "first_SCL" in cols:
                choices["scl"] = "UniProt localization"
            baits = dataset_baits(active_scored())
        ui.update_select("prey_pca_color", choices=choices)
        ui.update_select("prey_pca_bait", choices=baits)

    @reactive.effect
    @reactive.event(input.active_dataset, datasets)
    def update_pca_normalization_default():
        # Intensity-scale data defaults to log2; spectral counts stay linear.
        # Fires only on dataset change, so a manual override sticks after it.
        if not active():
            return
        quant = dataset_store.row(active())['Quant Type']
        default = "log2_zscore" if quant in ("Intensity", "LFQ") else "zscore"
        ui.update_select("pca_normalization", selected=default)

    # Compare Baits tab: two baits of the active dataset, each with its own
    # threshold set.

    @reactive.Calc
    def comp_volcano_data_cached():
        dataset = active_scored()
        bait_a = input.comp_bait_a()
        bait_b = input.comp_bait_b()
        if not all([dataset, bait_a, bait_b]):
            return pd.DataFrame()
        return calculate_volcano_data(dataset, bait_a, bait_b, read_thresholds("comp_a"),
                                      read_thresholds("comp_b"), out_dir)

    @reactive.Calc
    def comp_filtered_data_a_cached():
        if not active_scored() or not input.comp_bait_a():
            return pd.DataFrame()
        return load_and_filter_bait_data(active_scored(), input.comp_bait_a(),
                                         read_thresholds("comp_a"), out_dir)

    @reactive.Calc
    def comp_filtered_data_b_cached():
        if not active_scored() or not input.comp_bait_b():
            return pd.DataFrame()
        return load_and_filter_bait_data(active_scored(), input.comp_bait_b(),
                                         read_thresholds("comp_b"), out_dir)

    @reactive.effect
    @reactive.event(input.active_dataset, datasets)
    def update_comp_baits():
        """Bait lists for the active dataset, B starting on the second bait, and both
        threshold sets starting from the shared one."""
        baits = dataset_baits(active_scored()) if active_scored() else []
        ui.update_select("comp_bait_a", choices=baits)
        ui.update_select("comp_bait_b", choices=baits, selected=baits[1] if len(baits) > 1 else None)
        for side in ("comp_a", "comp_b"):
            set_thresholds(side, thresholds())

    for _side in ("a", "b"):
        def _make_use_shared(side):
            @reactive.effect
            @reactive.event(getattr(input, f"comp_use_shared_{side}"))
            def _apply():
                set_thresholds(f"comp_{side}", thresholds())
        _make_use_shared(_side)

    @render.ui
    def volcano_body():
        if not input.compare_networks():
            return plot_placeholder("Choose two baits and click Compare.")
        return output_widget("volcano_plot")

    @render.ui
    def venn_body():
        if not input.compare_networks():
            return plot_placeholder("Choose two baits and click Compare.")
        return ui.output_plot("venn_diagram")

    @render_widget
    @reactive.event(input.compare_networks)
    def volcano_plot():
        """Render volcano plot comparing two baits."""
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Comparing networks...", detail="Initializing", value=0)

            bait_a = input.comp_bait_a()
            bait_b = input.comp_bait_b()
            if not all([active_scored(), bait_a, bait_b]):
                fig = go.Figure()
                fig.add_annotation(
                    text="Please select two baits to compare",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font=dict(size=16)
                )
                fig.update_layout(
                    xaxis=dict(visible=False),
                    yaxis=dict(visible=False),
                    height=500
                )
                return fig

            progress.set(message="Comparing networks...", detail="Loading cached data", value=33)

            # Get cached volcano data
            volcano_data = comp_volcano_data_cached()

            progress.set(message="Comparing networks...", detail="Generating visualization", value=66)

            # Create plot
            fig = create_volcano_plot(volcano_data, bait_a, bait_b)

            progress.set(value=100)
            return fig

    @render.plot
    @reactive.event(input.compare_networks)
    def venn_diagram():
        """Render Venn diagram showing network overlap."""
        import matplotlib.pyplot as plt

        bait_a = input.comp_bait_a()
        bait_b = input.comp_bait_b()
        if not all([active_scored(), bait_a, bait_b]):
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, "Select two baits to view overlap",
                   ha='center', va='center', fontsize=14, color='#666')
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
            return fig

        # Get cached filtered data
        data_a = comp_filtered_data_a_cached()
        data_b = comp_filtered_data_b_cached()

        # Get sets of Prey.IDs
        set_a = set(data_a['Prey.ID'].unique()) if len(data_a) > 0 else set()
        set_b = set(data_b['Prey.ID'].unique()) if len(data_b) > 0 else set()

        # Create Venn diagram using matplotlib version
        fig = create_venn_diagram_matplotlib(set_a, set_b, bait_a, bait_b)

        return fig

    def _gene_list(region, label, empty):
        if not all([active_scored(), input.comp_bait_a(), input.comp_bait_b()]):
            return "Select two baits to view gene lists"
        genes = backend.compare_sets(comp_filtered_data_a_cached(), comp_filtered_data_b_cached())[region]
        if not genes:
            return empty
        return f"{label} ({len(genes)} genes):\n\n" + "\n".join(genes)

    @render.text
    @reactive.event(input.compare_networks)
    def genes_a_only():
        return _gene_list('a_only', "Network A only", "No unique genes in Network A")

    @render.text
    @reactive.event(input.compare_networks)
    def genes_b_only():
        return _gene_list('b_only', "Network B only", "No unique genes in Network B")

    @render.text
    @reactive.event(input.compare_networks)
    def genes_both():
        return _gene_list('both', "Both networks", "No shared genes between networks")

    @reactive.effect
    def update_dl_threshold_label():
        ui.update_switch("dl_apply_thresholds", label=f"Apply thresholds ({describe(thresholds())})")

    @reactive.Effect
    @reactive.event(input.active_dataset, datasets, feature_enrichment)
    def update_dl_presets():
        """Offer only the presets whose required files exist for this dataset."""
        dataset_name = active()
        avail = (dp.available_presets(os.path.join(out_dir, dataset_name))
                 if dataset_name else [])
        # A radio group cannot render zero choices (Shiny force-selects the
        # first), so with no dataset or no usable preset fall back to the full
        # registry; previews stay empty until a dataset provides the files.
        choices = ({key: dp.PRESETS[key].label for key in avail}
                   or {key: preset.label for key, preset in dp.PRESETS.items()})
        selected = input.dl_preset.get()
        if selected not in choices:
            selected = next(iter(choices))
        ui.update_radio_buttons("dl_preset", choices=choices, selected=selected)

    @render.ui
    @reactive.event(input.dl_preset, input.active_dataset)
    def dl_groups_ui():
        """The Include checkboxes for the preset's column groups or file choices,
        preserving still-valid picks across switches; nothing for a preset
        without items."""
        preset = dp.PRESETS.get(input.dl_preset.get())
        items = () if preset is None else (
            preset.files if preset.kind == "files" else preset.groups)
        if not items:
            return None
        choices = {item.key: item.label for item in items}
        with reactive.isolate():
            picked = input.dl_groups.get() if "dl_groups" in input else ()
        current = [k for k in picked if k in choices]
        selected = current or [item.key for item in items if item.default]
        return ui.input_checkbox_group("dl_groups", tip("Include", "dl_groups"),
                                       choices=choices, selected=selected)

    @reactive.Effect
    @reactive.event(input.active_dataset, datasets)
    def update_selectize_custom_columns():
        """Offer the active dataset's full column list in the Custom preset,
        keeping the user's still-valid picks. Client-side options only: the
        column list is small and server-side selectize never delivers options
        to the browser in this app."""
        dataset_name = active()
        scores_path = (os.path.join(out_dir, dataset_name, "annotated_scores.csv")
                       if dataset_name else "")
        if not scores_path or not os.path.exists(scores_path):
            ui.update_selectize("custom_columns", choices=DEFAULT_CUSTOM_COLUMNS,
                                selected=DEFAULT_CUSTOM_COLUMNS)
            return
        cols = sorted(pd.read_csv(scores_path, nrows=0).columns.tolist())
        keep = ([c for c in input.custom_columns.get() if c in cols]
                or DEFAULT_CUSTOM_COLUMNS)
        ui.update_selectize("custom_columns", choices=cols, selected=keep)

    dl_result = reactive.Value(pd.DataFrame())
    dl_result_total = reactive.Value(0)
    dl_result_files = reactive.Value([])

    @reactive.Calc
    def cached_download_data():
        """Cache the active dataset's annotated_scores.csv for the Export tab.
        Only re-reads when the active dataset changes, not on slider changes."""
        dataset = active()
        if not dataset:
            return pd.DataFrame()
        dataset_path = os.path.join(out_dir, dataset, "annotated_scores.csv")
        if not os.path.exists(dataset_path):
            return pd.DataFrame()
        return pd.read_csv(dataset_path)

    @reactive.Calc
    def cached_enrichment_data():
        """Cache Feature_enrichment.csv for the active dataset; empty frame
        when the Feature Enrichment tab has not produced one."""
        dataset = active()
        if not dataset:
            return pd.DataFrame()
        path = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(path):
            return pd.DataFrame()
        return pd.read_csv(path)

    def _dl_thresholds():
        return thresholds() if input.dl_apply_thresholds() else NO_FILTER

    @render.data_frame
    def dl_preview_table():
        """Build the selected preset's export and preview it. The result is
        stored in dl_result / dl_result_files so the download handler writes
        exactly what is previewed."""
        dataset = active()
        preset_key = input.dl_preset.get()
        preset = dp.PRESETS.get(preset_key)
        dl_result.set(pd.DataFrame())
        dl_result_total.set(0)
        dl_result_files.set([])
        if not dataset or preset is None:
            return pd.DataFrame()
        dataset_dir = os.path.join(out_dir, dataset)
        # Checkbox state can still belong to the previously shown preset while
        # the update_dl_groups round-trip is in flight; sanitize it.
        groups = dp.effective_selection(preset, input.dl_groups.get())
        try:
            if preset_key == 'ed':
                ed_path = os.path.join(dataset_dir, "ED.csv")
                if not os.path.exists(ed_path):
                    return pd.DataFrame()
                result = pd.read_csv(ed_path)
                dl_result_total.set(len(result))
            elif preset_key == 'saint_inputs':
                paths, missing = dp.saint_input_files(dataset_dir, groups)
                if missing:
                    notify("Not produced by this run: " + ", ".join(missing),
                           type="warning")
                dl_result_files.set(paths)
                dl_result_total.set(len(paths))
                return pd.DataFrame({
                    'File': [os.path.basename(p) for p in paths],
                    'Size (KB)': [round(os.path.getsize(p) / 1024, 1)
                                  for p in paths],
                })
            elif preset_key == 'enrichment':
                df = cached_enrichment_data()
                if df.empty:
                    return pd.DataFrame()
                result = dp.build_enrichment_table(df, groups)
                dl_result_total.set(len(df))
            else:
                df = cached_download_data()
                if df.empty:
                    return pd.DataFrame()
                dl_result_total.set(len(df))
                thresholds = _dl_thresholds()
                if preset_key == 'annotated':
                    result = dp.build_annotated_table(df, groups, thresholds)
                elif preset_key == 'cytoscape':
                    result = dp.build_cytoscape_edges(df, thresholds)
                elif preset_key == 'genelist':
                    result = dp.build_gene_list(df, thresholds,
                                                mode=input.dl_genelist_mode())
                elif preset_key == 'prohits':
                    result = dp.build_prohits_table(
                        df, thresholds,
                        abundance_col=input.dl_prohits_abundance())
                else:  # custom
                    result, missing = dp.build_custom_table(
                        df, list(input.custom_columns.get()), thresholds)
                    if missing:
                        notify("Columns not in this dataset: " + ", ".join(missing),
                               type="warning")
        except ValueError as err:
            notify(str(err), type="error")
            return pd.DataFrame()
        dl_result.set(result)
        return result

    @render.text
    def download_row_count():
        """Preset-aware summary line above the preview."""
        preset = dp.PRESETS.get(input.dl_preset.get())
        total = dl_result_total.get()
        if preset is None or total == 0:
            return ""
        if preset.kind == 'files':
            return f"{total} files selected"
        if preset.kind == 'genelist':
            return f"{len(dl_result.get())} genes from {total} interactions"
        if preset.uses_thresholds:
            return f"Showing {len(dl_result.get())} of {total} interactions"
        return f"{len(dl_result.get())} rows"

    @render.download_button()
    def download_preset():
        dataset = active()
        preset = dp.PRESETS.get(input.dl_preset.get())
        if not dataset or preset is None:
            notify("No dataset selected.", type="error")
            return None
        if preset.kind == 'files':
            if not dl_result_files.get():
                notify("No files selected to download.", type="error")
                return None
        elif dl_result.get().empty:
            notify("Nothing to download with the current selection.", type="error")
            return None
        notify("Preparing download...", type="message", duration=2)
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        savepath = os.path.join(
            out_dir, f"{dataset}_{preset.key}_{timestamp}{preset.extension}")
        if preset.extension == '.zip':
            dp.zip_files(dl_result_files.get(), savepath)
        elif preset.extension == '.txt':
            dp.write_gene_list(dl_result.get(), savepath)
        else:
            dl_result.get().to_csv(savepath, index=False)
        return savepath

    @render.download_button()
    def download_batch():
        """Download all results for a dataset as a ZIP file."""
        dataset = active()
        if not dataset:
            notify("No dataset selected.", type="error")
            return None

        notify("Preparing ZIP file...", type="message", duration=2)
        dataset_dir = os.path.join(out_dir, dataset)
        if not os.path.exists(dataset_dir):
            notify("Dataset directory not found.", type="error")
            return None

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        zip_filename = f"{dataset}_all_{timestamp}.zip"
        zip_path = os.path.join(out_dir, zip_filename)

        # List of files to include in the ZIP. The log, the run manifest and the
        # database build stamp travel with the results so a recipient can see how
        # they were produced.
        files_to_include = [
            'merged.csv',
            'annotated_scores.csv',
            'Feature_enrichment.csv',
            'bait.txt',
            'prey.txt',
            'interaction.txt',
            'ED.csv',
            'run.json',
            'proximate.log',
            'build_info.txt'
        ]

        included_files = []
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
            for filename in files_to_include:
                filepath = os.path.join(dataset_dir, filename)
                if os.path.exists(filepath):
                    zf.write(filepath, filename)
                    included_files.append(filename)

        if not included_files:
            notify("No files found to include in ZIP.", type="error")
            os.remove(zip_path)
            return None

        notify(f"ZIP created with {len(included_files)} files.", type="message", duration=3)
        return zip_path

    # Cytoscape tab.  The controller's state is shared by every browser session, so
    # the tab learns about changes by polling its version counter.
    @reactive.poll(lambda: cytoscape_ctl.STATE['version'], 1.0)
    def cy_version():
        return cytoscape_ctl.STATE['version']

    cy_health = reactive.Value(None)
    cy_selection = reactive.Value(pd.DataFrame())

    def cy_call(doing, fn, *args, **kwargs):
        """Run a controller operation, turning its failure into a notification."""
        try:
            return fn(*args, **kwargs)
        except Exception as e:
            logger.exception("Cytoscape: could not %s", doing)
            notify(f"Could not {doing}: {e}", type="error", duration=None)
            return None

    @reactive.effect
    @reactive.event(input.active_dataset, datasets)
    def update_cy_baits():
        baits = sorted(dataset_baits(active_scored())) if active_scored() else []
        ui.update_selectize("cy_baits", choices=baits, selected=[])

    @reactive.effect
    @reactive.event(input.cy_style_preset, ignore_init=True)
    def apply_style_preset_options():
        """A preset brings the send options it is designed for; the user can still change them."""
        opts = STYLE_PRESETS[input.cy_style_preset.get()]["send_options"]
        ui.update_select("cy_edge_width", selected=opts["width_source"])
        ui.update_checkbox("cy_prey_prey", value=opts["prey_prey"])
        ui.update_select("cy_biogrid_scope", selected=opts["biogrid_scope"])
        ui.update_checkbox("cy_lit_weighted", value=opts["literature_weighted"])
        ui.update_checkbox("cy_corum", value=opts["corum"])
        ui.update_select("cy_layout", selected=opts["layout"])

    @reactive.effect
    @reactive.event(input.cy_probe)
    def cy_probe():
        cy_health.set(cytoscape_ctl.health())

    @render.ui
    def cy_status():
        cy_version()
        health = cy_health.get()
        if health is None:
            health = cytoscape_ctl.health()
            cy_health.set(health)
        snap = cytoscape_ctl.snapshot()
        if health['ok']:
            line = ui.p(ui.span("● ", class_="status-ok"),
                        f"Cytoscape {health['version']} at {health['url']}")
        else:
            line = ui.div(
                ui.p(ui.span("● ", class_="status-bad"),
                     f"Cytoscape is not reachable at {health['url']}. Start Cytoscape on this "
                     "machine or set PROXIMATE_CYTOSCAPE_URL.", class_="status-error"),
                ui.accordion(ui.accordion_panel("Details", ui.tags.pre(health['error'])),
                             open=False, class_="mb-3"))
        if snap['net_suid'] is None:
            drawn = ui.p("No ProxiMate network drawn yet.", class_="help-text")
        elif snap['dataset'] is None:
            drawn = ui.p(f"{snap['title']} is still drawn, but its dataset no longer exists.",
                         class_="status-error")
        else:
            drawn = ui.p(f"{snap['title']}: {snap['n_nodes']} nodes, {snap['n_edges']} edges "
                         f"({snap['n_hidden']} hidden)" + (f" — {snap['busy']}" if snap['busy'] else ""))
        return ui.div(line, drawn)

    @render.ui
    def cy_thresholds_line():
        return ui.p(f"Thresholds (sidebar): {describe(thresholds())}", class_="help-text")

    @reactive.calc
    def cy_estimate_args():
        """The ``backend.estimate_network`` arguments the send options currently give."""
        return dict(name=active_scored(), thresholds=thresholds(),
                    baits=list(input.cy_baits.get() or []), prey_prey=input.cy_prey_prey.get(),
                    biogrid_scope=input.cy_biogrid_scope.get(), corum=input.cy_corum.get(),
                    corum_min_members=int(input.cy_corum_min_members.get() or 3),
                    corum_min_fraction=float(input.cy_corum_min_fraction.get() or 0.0))

    def estimate_text(counts):
        text = (f"{counts['baits']} baits, {counts['nodes']} nodes, {counts['edges']} edges "
                f"({counts['proximity_edges']} bait-prey, {counts['literature_edges']} BioGRID, "
                f"{counts['complex_edges']} complex)")
        if counts['hidden_edges']:
            text += f", plus {counts['hidden_edges']} BioGRID edges hidden by the scope"
        return text

    # The estimate waits until the send options have been still for half a second, so
    # dragging a slider builds the network once rather than at every step.
    cy_estimate_request = reactive.Value((None, 0.0))

    @reactive.effect
    def cy_queue_estimate():
        cy_estimate_request.set((cy_estimate_args(), time.monotonic() + 0.5))

    @render.ui
    def cy_estimate():
        args, due = cy_estimate_request.get()
        wait = due - time.monotonic()
        if wait > 0:
            reactive.invalidate_later(wait)
            req(False, cancel_output=True)
        if not args or not args['name']:
            return None
        try:
            counts = backend.estimate_network(**args)
        except ValueError as e:
            return ui.p(str(e), class_="cy-estimate error")
        return ui.p("Will draw: " + estimate_text(counts), class_="cy-estimate")

    @reactive.effect
    @reactive.event(input.cy_send)
    def cy_send():
        args = cy_estimate_args()
        if not args['name']:
            notify("Score a dataset first.", type="error")
            return
        try:
            counts = backend.estimate_network(**args)
        except ValueError as e:
            notify(f"Could not send to Cytoscape: {e}", type="error")
            return
        if counts['edges'] <= SEND_WARN_EDGES:
            cy_start_send()
            return
        ui.modal_show(ui.modal(
            ui.p(f"This network has {estimate_text(counts)}: more than {SEND_WARN_EDGES} edges, "
                 "which Cytoscape draws slowly and which is hard to read."),
            ui.p("Consider the \"Multivalidated only\" BioGRID scope or fewer baits."),
            title="Large network", easy_close=True,
            footer=ui.div(
                ui.modal_button("Cancel"),
                ui.tags.button("Send anyway", type="button", class_="btn btn-primary ms-2",
                               onclick="Shiny.setInputValue('confirm_cy_send', 1, {priority: 'event'})"))))

    @reactive.effect
    @reactive.event(input.confirm_cy_send)
    def confirm_cy_send():
        ui.modal_remove()
        cy_start_send()

    def cy_start_send():
        dataset_name = active_scored()
        dataset_path = os.path.join(out_dir, dataset_name)
        try:
            biogrid_path, corum_path = backend.network_sources(dataset_name, input.cy_corum.get())
        except ValueError as e:
            notify(f"Could not send to Cytoscape: {e}", type="error")
            return
        cutoffs = thresholds()
        options = dict(baits=list(input.cy_baits.get() or []), prey_prey=input.cy_prey_prey.get(),
                       biogrid_path=biogrid_path, label_policy=input.cy_labels.get(),
                       layout=input.cy_layout.get(), width_source=input.cy_edge_width.get(),
                       literature_weighted=input.cy_lit_weighted.get(),
                       biogrid_scope=input.cy_biogrid_scope.get(), corum_path=corum_path,
                       corum_min_members=int(input.cy_corum_min_members.get() or 3),
                       corum_min_fraction=float(input.cy_corum_min_fraction.get() or 0.0),
                       preset=input.cy_style_preset.get())

        def work(progress):
            return cytoscape_ctl.draw(dataset_name, os.path.join(dataset_path, "annotated_scores.csv"),
                                      cutoffs, progress=progress, **options)

        def done(snap):
            notify(f"Drew {snap['n_nodes']} nodes and {snap['n_edges']} edges in Cytoscape.")

        start_job(f"Sending {dataset_name} to Cytoscape", dataset_name, cytoscape_ctl.SEND_STAGES,
                  work, done)

    @reactive.effect
    @reactive.event(input.cy_rethreshold)
    def cy_rethreshold():
        hidden = cy_call("re-apply the thresholds", cytoscape_ctl.apply_thresholds, thresholds())
        if hidden is not None:
            notify(f"Thresholds applied: {hidden} edge(s) hidden.")

    @reactive.effect
    @reactive.event(input.cy_restyle)
    def cy_restyle():
        changed = cy_call("apply the style", cytoscape_ctl.restyle, input.cy_style_preset.get(),
                          input.cy_edge_width.get(), input.cy_lit_weighted.get(),
                          input.cy_biogrid_scope.get())
        if changed is not None:
            notify(f"Style applied: {changed} edge(s) changed.")

    def cy_show_chosen(ids):
        """Put a fresh selection in the table so the tab shows what Cytoscape now has."""
        nodes = cytoscape_ctl.STATE['nodes']
        cy_selection.set(nodes[nodes['id'].isin(ids)][['id', 'symbol', 'role']].reset_index(drop=True))

    @reactive.effect
    @reactive.event(input.cy_select_loners)
    def cy_select_loners():
        chosen = cy_call("select the loners", cytoscape_ctl.select_loners)
        if chosen:
            cy_show_chosen(chosen)
            notify(f"Selected the bait and its {len(chosen) - 1} loner(s).")

    @reactive.effect
    @reactive.event(input.cy_select_satellites)
    def cy_select_satellites():
        chosen = cy_call("select the satellites", cytoscape_ctl.select_satellites)
        if chosen:
            cy_show_chosen(chosen)
            notify(f"Selected the bait and its {len(chosen) - 1} satellite(s).")

    def cy_select_related(add):
        seed = (input.cy_rel_seed.get() or "").strip()
        if not seed:
            notify("Name a seed bait or protein first.", type="error")
            return
        cuts = {'min_saint': input.cy_rel_saint.get(), 'max_bfdr': input.cy_rel_bfdr.get(),
                'min_abundance': input.cy_rel_abundance.get(), 'min_publications': input.cy_rel_pubs.get()}
        chosen = cy_call("select by relation", cytoscape_ctl.select_related, seed,
                         input.cy_rel_kind.get(), add=add, **cuts)
        if chosen:
            cy_show_chosen(chosen)
            notify(f"Selected {len(chosen)} {input.cy_rel_kind.get()} of {seed}"
                   + (" (added to the selection)." if add else "."))

    @reactive.effect
    @reactive.event(input.cy_rel_replace)
    def cy_rel_replace():
        cy_select_related(add=False)

    @reactive.effect
    @reactive.event(input.cy_rel_add)
    def cy_rel_add():
        cy_select_related(add=True)

    @reactive.effect
    @reactive.event(input.cy_cluster)
    def cy_cluster():
        resolution, seed, weight = (input.cy_cl_resolution.get(), input.cy_cl_seed.get(),
                                    input.cy_cl_lit_weight.get())
        if resolution is None or seed is None or weight is None:
            notify("Fill in resolution, seed and reference weight first.", type="error")
            return
        result = cy_call("cluster the selection", cytoscape_ctl.cluster_selection,
                         resolution=float(resolution), seed=int(seed), literature_weight=float(weight),
                         layout=input.cy_cl_layout.get(), recolor=input.cy_cl_recolor.get())
        if result:
            notify(f"{result['n']} preys clustered into {result['n_communities']} communities "
                   f"(sizes {result['sizes']}), {result['recolored']} recolored; "
                   f"{len(result['baits_left'])} bait(s) left in place.")

    @reactive.effect
    @reactive.event(input.cy_read_selection)
    def cy_read_selection():
        result = cy_call("read the selection", cytoscape_ctl.read_selection)
        if result:
            chosen, detail = result
            cy_selection.set(selection_table(detail, cytoscape_ctl.STATE['nodes']) if len(detail) else chosen)
            notify(f"{len(chosen)} node(s) selected, touching {len(detail)} edge(s).")

    for _button, _action in (("cy_hide_selected", "hide_selected"),
                             ("cy_show_selected", "show_selected"),
                             ("cy_hide_unselected", "hide_unselected"),
                             ("cy_show_all", "show_all")):
        def _make_visibility(button, action):
            @reactive.effect
            @reactive.event(getattr(input, button))
            def _apply():
                changed = cy_call(f"{action.replace('_', ' ')} edges",
                                  cytoscape_ctl.set_edge_visibility, action)
                if changed is not None:
                    notify(f"{changed} edge(s) changed.")
        _make_visibility(_button, _action)

    def cy_layout_listing():
        """The drawn dataset and its saved layouts; save and load act on the drawn network."""
        dataset = cytoscape_ctl.STATE['dataset']
        names = cytoscape_ctl.saved_layouts(os.path.join(out_dir, dataset)) if dataset else []
        return dataset, tuple(names)

    # Polled, so the list follows draws and saves from the GUI or MCP, and only a change
    # in the listing re-renders the control: a name being typed is never wiped.
    @reactive.poll(cy_layout_listing, 2.0)
    def cy_saved_layouts():
        return list(cy_layout_listing()[1])

    @reactive.effect
    def update_cy_layout_names():
        names = cy_saved_layouts()
        with reactive.isolate():
            current = input.cy_layout_name.get()
        # The empty choice is the placeholder, so no layout is picked until the user picks one.
        ui.update_selectize("cy_layout_name", choices=[''] + names,
                            selected=current if current in names else '')

    def cy_layout_target():
        """``(snapshot, dataset_path, name)`` for the drawn network, or None after a notification."""
        snap = cytoscape_ctl.snapshot()
        name = (input.cy_layout_name.get() or "").strip()
        if snap['dataset'] is None:
            notify("Send a network to Cytoscape first.", type="error")
            return None
        if not name:
            notify("Pick or type a layout name first.", type="error")
            return None
        return snap, os.path.join(out_dir, snap['dataset']), name

    @reactive.effect
    @reactive.event(input.cy_save_layout)
    def cy_save_layout():
        target = cy_layout_target()
        if not target:
            return
        snap, dataset_path, name = target

        def save_with_record():
            with provenance.stage(dataset_path, "cytoscape", entrypoint="cytoscape_ctl.save_layout",
                                  params={'name': name}) as record:
                result = cytoscape_ctl.save_layout(dataset_path, name)
                record.add_output(result['path'], role="layout", rows=result['n_rows'])
                record.extra(thresholds=snap['thresholds'])
                return result

        result = cy_call("save the layout", save_with_record)
        if result:
            notify(f"Layout {name!r} saved: {result['written']} drawn node(s), {result['n_rows']} in the file.")

    @reactive.effect
    @reactive.event(input.cy_load_layout)
    def cy_load_layout():
        target = cy_layout_target()
        if not target:
            return
        _, dataset_path, name = target
        result = cy_call("load the layout", cytoscape_ctl.load_layout, dataset_path, name)
        if result:
            notify(f"Layout {name!r}: {result['moved']} node(s) placed, "
                   f"{result['unplaced']} left where they were.")

    @reactive.effect
    @reactive.event(input.cy_export)
    def cy_export():
        snap = cytoscape_ctl.snapshot()
        if snap['dataset'] is None:
            notify("Send a network to Cytoscape first.", type="error")
            return
        dataset_path = os.path.join(out_dir, snap['dataset'])

        def export_with_record():
            with provenance.stage(dataset_path, "cytoscape",
                                  entrypoint="cytoscape_ctl.export_image") as record:
                path = cytoscape_ctl.export_image(dataset_path)
                record.add_output(path, role="image")
                record.extra(thresholds=snap['thresholds'])
                return path

        path = cy_call("export the image", export_with_record)
        if path:
            notify(f"Image written to {path}")

    @reactive.effect
    @reactive.event(input.cy_unlock)
    def cy_unlock():
        held = cy_call("unlock the view", cytoscape_ctl.unlock)
        if held is not None:
            notify("Released " + (", ".join(held) if held else "nothing; the view was not locked."))

    @render.data_frame
    def cy_selection_table():
        if cy_selection.get().empty:
            return None
        return render.DataGrid(cy_selection.get(), width="100%", height="400px")

    @render.text
    def cy_selection_count():
        n = len(cy_selection.get())
        return f"({n} rows)" if n else ""

    @render.download_button(filename=lambda: f"{cytoscape_ctl.STATE['dataset'] or 'cytoscape'}_selection.csv")
    def cy_download_selection():
        if cy_selection.get().empty:
            notify("Read a selection first.", type="error")
            return None
        path = os.path.join(out_dir, "cytoscape_selection.csv")
        cy_selection.get().to_csv(path, index=False)
        return path

    @render.ui
    def cy_activity():
        cy_version()
        entries = cytoscape_ctl.snapshot()['log']
        if not entries:
            return ui.p("No Cytoscape activity yet.", class_="help-text")
        header = ui.tags.tr(*[ui.tags.th(h) for h in ("Time", "Actor", "Operation", "Detail")])
        rows = [ui.tags.tr(ui.tags.td(e['ts'][11:]), ui.tags.td(e['actor']), ui.tags.td(e['op']),
                           ui.tags.td(e['detail']))
                for e in reversed(entries)]
        return ui.div(ui.tags.table(ui.tags.thead(header), ui.tags.tbody(*rows), class_="table table-sm"),
                      class_="cy-activity")

def log_startup():
    """Record the configuration the server came up with.

    Written at import so it also appears when an ASGI server loads this module
    rather than running it as a script.  Without it the operational log stays
    empty until someone acts, and there is no way to confirm from the logs which
    version is serving or where it is writing.
    """
    version = provenance.proximate_version()
    logger.info("ProxiMate starting: version=%s (%s), python=%s",
                version["version"], version["source"], platform.python_version())
    logger.info("Datasets in %s; operational log in %s; LOG_LEVEL=%s",
                out_dir, os.environ.get("PROXIMATE_LOG_DIR", log_config.DEFAULT_LOG_DIR),
                logging.getLevelName(logging.getLogger(log_config.PACKAGE).level))
    build_info = provenance.read_build_info()
    if build_info:
        first_line = next((l for l in build_info.splitlines() if l.startswith("Build date")), None)
        logger.info("Annotation databases: %s", first_line or "build_info.txt present")
    else:
        logger.warning("No /Datasets/build_info.txt; annotation database versions are unknown")


log_startup()

app = App(app_ui, server)

if __name__ == "__main__":
    run_app(app, host=os.environ.get("PROXIMATE_GUI_HOST", "0.0.0.0"),
            port=int(os.environ.get("PROXIMATE_GUI_PORT", "3838")))
