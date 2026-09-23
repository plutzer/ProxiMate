from functools import partial
from shiny import App, Inputs, Outputs, Session, reactive, render, ui, run_app
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
from QC_plots import pca_plot, prepare_pca_matrix, prey_pca_plot as plot_prey_pca, prey_gene_names, detection_counts, reduce_categorical, saint_known_retention, roc_plot, saint_scatter_plot as plot_saint_scatter, calculate_threshold_metrics
from Ann_Enrichment import process_refactored, plot_results
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
from help_text import tip, TOOLTIPS
import session_archive
import backend
import dataset_store
import mcp_registry
import cytoscape_ctl
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

app_ui = ui.page_navbar(
    ui.nav_spacer(),
    ui.nav_panel("Network Scoring",
                 ui.layout_columns(
                    ui.card(
                        ui.card_header("Data Parsing"),
                        ui.layout_columns(
                            ui.input_text("dataset_name",
                                        tip("Dataset Name", "dataset_name"),
                                        placeholder="No spaces or special characters (/ \\ : * ? \" < > |)"),
                            ui.input_select("input_format", tip("Input Format", "input_format"),
                                           choices=["MaxQuant", "DIA-NN", "Pioneer", "FragPipe", "MSstats", "SAINT"],
                                           selected="MaxQuant"),
                            col_widths=[4, 8]
                        ),
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
                                        class_="text-muted small d-block mt-2"
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
                        ui.input_action_button("parse_data", "Parse Data")
                    ),
                    ui.card(
                        ui.card_header('Scoring Parameters'),
                        ui.input_select("score_dataset", "Select Dataset", choices=[]),
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
                        ui.input_numeric("wdfdr_iterations", tip("WDFDR Iterations", "wdfdr_iterations"), value=1000),
                        ui.input_action_button("score_data", "Score Data")
                    ),
                    col_widths=[8,4]
                 ),
                 # Render text if parse button has executed successfully
                 ui.card(
                     ui.card_header(
                        ui.h1('Datasets in this Session'),
                        ui.layout_columns(
                            ui.card(
                                ui.tooltip(ui.input_action_button("clear_datasets", "Clear All Datasets"),
                                           TOOLTIPS["clear_datasets"]),
                                ui.tooltip(ui.download_button("download_session", "Download Session Zip"),
                                           TOOLTIPS["download_session"]),
                            ),
                            ui.card(
                                ui.input_file("session_file", tip("Upload Session Zip", "session_file")),
                                ui.input_action_button("upload_session", "Load Session"),
                            ),
                        )
                    ),
                    ui.output_data_frame("render_datasets"),
                ),
            ),
    ui.nav_panel("Data Thresholding",
                    # Top row: Dataset selector only
                    ui.input_select("qc_dataset", "Select Dataset", choices=[]),
                    ui.output_ui("empty_state_thresholding"),

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
                            ui.card_header(tip("Experiment PCA", "experiment_pca")),
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

                    # Row 2: threshold controls | SAINT scatter | metrics
                    ui.layout_columns(
                        ui.card(
                            ui.card_header("Threshold Settings"),
                            ui.input_select("qc_bait", tip("Select QC Bait", "qc_bait"), choices=["All"]),
                            ui.input_slider("threshold_saintscore", tip("SAINT Score Threshold", "saintscore"),
                                          min=0.0, max=1.0, value=0.7, step=0.01),
                            ui.input_slider("threshold_bfdr", tip("BFDR Threshold", "bfdr"),
                                          min=0.0, max=1.0, value=0.05, step=0.01),
                            ui.input_slider("threshold_wd", tip("WD Score Threshold", "wd"),
                                          min=0.0, max=10.0, value=0.0, step=0.1),
                            ui.input_slider("threshold_wdfdr", tip("WDFDR Threshold", "wdfdr"),
                                          min=0.0, max=1.0, value=1.0, step=0.01),
                            ui.p("Presets:", style="margin-top: 15px; margin-bottom: 5px; font-weight: 500;"),
                            # d-flex, not layout_columns: layout_columns collapses to
                            # stacked full-width rows inside a narrow card
                            ui.div(
                                ui.tooltip(ui.input_action_button("qc_preset_stringent", "Stringent", class_="btn-sm btn-outline-primary"),
                                           TOOLTIPS["preset_stringent"]),
                                ui.tooltip(ui.input_action_button("qc_preset_moderate", "Moderate", class_="btn-sm btn-outline-secondary"),
                                           TOOLTIPS["preset_moderate"]),
                                ui.tooltip(ui.input_action_button("qc_preset_relaxed", "Relaxed", class_="btn-sm btn-outline-secondary"),
                                           TOOLTIPS["preset_relaxed"]),
                                class_="d-flex gap-2 mb-3",
                            ),
                            # ui.p("Note: Thresholds are shown as reference lines on plots. Data is not filtered.",
                            #      style="font-style: italic; color: #666; margin-top: 10px;"),
                            ui.output_ui("wdfdr_warning"),
                        ),
                        ui.card(
                            ui.card_header("SAINT Score vs Fold Change"),
                            output_widget("saint_scatter_plot"),
                            ui.download_button("download_scatter_plot", "Export PNG", class_="btn-sm"),
                        ),
                        ui.div(
                            ui.output_ui("metric_bait_label"),
                            ui.output_ui("metric_network_size"),
                            ui.output_ui("metric_enrichment"),
                            ui.output_ui("metric_degree"),
                        ),
                        col_widths=(4, 5, 3),
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
    ui.nav_panel("Protein Feature Analysis",
                ui.output_ui("empty_state_feature_analysis"),
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Parameters for Feature Analysis"),
                        ui.input_select("feature_dataset", "Select Dataset", choices=[]), # Need this to be dynamic
                        ui.p("Preys passing all four scores form the foreground tested for enrichment.",
                             class_="text-muted small"),
                        ui.input_slider("fa_threshold_saintscore", tip("SAINT Score (≥)", "saintscore"),
                                      min=0.0, max=1.0, value=0.9, step=0.01),
                        ui.input_slider("fa_threshold_bfdr", tip("BFDR (≤)", "bfdr"),
                                      min=0.0, max=1.0, value=1.0, step=0.01),
                        ui.input_slider("fa_threshold_wd", tip("WD Score (≥)", "wd"),
                                      min=0.0, max=10.0, value=0.0, step=0.1),
                        ui.input_slider("fa_threshold_wdfdr", tip("WDFDR (≤)", "wdfdr"),
                                      min=0.0, max=1.0, value=1.0, step=0.01),
                        ui.p("Presets:", style="margin-top: 15px; margin-bottom: 5px; font-weight: 500;"),
                        # Flex row rather than layout_columns, whose columns collapse to
                        # full-width stacked rows at this card's width.
                        ui.div(
                            ui.tooltip(ui.input_action_button("fa_preset_stringent", "Stringent",
                                                              class_="btn-sm btn-outline-primary flex-fill"),
                                       TOOLTIPS["preset_stringent"]),
                            ui.tooltip(ui.input_action_button("fa_preset_moderate", "Moderate",
                                                              class_="btn-sm btn-outline-secondary flex-fill"),
                                       TOOLTIPS["preset_moderate"]),
                            ui.tooltip(ui.input_action_button("fa_preset_relaxed", "Relaxed",
                                                              class_="btn-sm btn-outline-secondary flex-fill"),
                                       TOOLTIPS["preset_relaxed"]),
                            class_="d-flex gap-2 mb-3",
                        ),
                        ui.tooltip(ui.input_action_button("feature_analysis", "Run Feature Analysis"),
                                   TOOLTIPS["feature_analysis"]),
                    ),
                    ui.card(
                        ui.card_header("Feature Enrichment Analysis"),
                        ui.input_select("feature_type", tip("Select Feature Type", "feature_type"), choices=["GO_CC", "GO_BP", "GO_MF", "Motifs", "Regions", "Repeats", "Compositions", "Domains"]),
                        ui.input_numeric("num_features", tip("Number of Features to Display", "num_features"), value=30, min=1, max=100),
                        # Thirty feature rows need the height to stay readable.
                        ui.output_plot("feature_enrichment_plot", height="650px"),
                        ui.download_button("download_heatmap", "Export Heatmap PNG", class_="btn-sm"),
                        ui.hr(),
                        ui.h5("Download Enrichment Results"),
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
                ),
    ),
    ui.nav_panel("Network Comparison",
        ui.output_ui("empty_state_network_comparison"),
        # Top section: Bait selectors for A and B
        ui.layout_columns(
            # Bait A selector card
            ui.card(
                ui.card_header("Network A"),
                ui.input_select("comp_dataset_a", "Dataset A", choices=[]),
                ui.input_select("comp_bait_a", tip("Bait A", "comp_bait"), choices=[]),
                ui.input_slider("comp_saintscore_a", tip("SAINT Score (≥)", "saintscore"),
                              min=0.0, max=1.0, value=0.7, step=0.01),
                ui.input_slider("comp_bfdr_a", tip("BFDR (≤)", "bfdr"),
                              min=0.0, max=1.0, value=0.05, step=0.01),
                ui.input_slider("comp_wd_a", tip("WD Score (≥)", "wd"),
                              min=0.0, max=10.0, value=0.0, step=0.1),
                # WDFDR compares a prey's WD across experiments, so among replicate
                # baits each prey passes in only its max-WD experiment; filtering on
                # it by default would empty the venn overlap.  Default = no filter.
                ui.input_slider("comp_wdfdr_a", tip("WDFDR (≤)", "wdfdr"),
                              min=0.0, max=1.0, value=1.0, step=0.01),
            ),
            # Bait B selector card
            ui.card(
                ui.card_header("Network B"),
                ui.input_select("comp_dataset_b", "Dataset B", choices=[]),
                ui.input_select("comp_bait_b", tip("Bait B", "comp_bait"), choices=[]),
                ui.input_slider("comp_saintscore_b", tip("SAINT Score (≥)", "saintscore"),
                              min=0.0, max=1.0, value=0.7, step=0.01),
                ui.input_slider("comp_bfdr_b", tip("BFDR (≤)", "bfdr"),
                              min=0.0, max=1.0, value=0.05, step=0.01),
                ui.input_slider("comp_wd_b", tip("WD Score (≥)", "wd"),
                              min=0.0, max=10.0, value=0.0, step=0.1),
                ui.input_slider("comp_wdfdr_b", tip("WDFDR (≤)", "wdfdr"),
                              min=0.0, max=1.0, value=1.0, step=0.01),
            ),
            col_widths=(6, 6),
        ),

        # Compare button
        ui.div(
            ui.input_action_button("compare_networks", "Compare Networks", class_="btn-primary"),
            style="text-align: center; margin: 20px 0;"
        ),

        # Middle section: Volcano plot
        ui.card(
            ui.card_header(tip("Differential Abundance Volcano Plot", "volcano_plot")),
            output_widget("volcano_plot"),
            ui.download_button("download_volcano_plot", "Export PNG", class_="btn-sm"),
            ui.p("Volcano plot only shown when baits are from the same dataset. "
                 "Flanking strips hold preys detected under only one bait "
                 "(y = -log10 BFDR); the central panel holds shared preys "
                 "(y = -log10 BH-adjusted p).",
                 class_="text-muted small mt-2"),
        ),

        # Bottom section: Venn diagram and gene lists
        ui.layout_columns(
            ui.card(
                ui.card_header(tip("Network Overlap", "venn")),
                ui.output_plot("venn_diagram"),
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
    ui.nav_panel("Cytoscape",
        ui.output_ui("empty_state_cytoscape"),
        ui.layout_columns(
            ui.card(
                ui.card_header("Network"),
                ui.input_select("cy_dataset", "Select Dataset", choices=[]),
                ui.input_selectize("cy_baits", tip("Baits", "cy_baits"), choices=[], multiple=True),
                ui.input_slider("cy_threshold_saintscore", tip("SAINT Score Threshold", "saintscore"),
                                min=0.0, max=1.0, value=0.7, step=0.01),
                ui.input_slider("cy_threshold_bfdr", tip("BFDR Threshold", "bfdr"),
                                min=0.0, max=1.0, value=0.05, step=0.01),
                ui.input_slider("cy_threshold_wd", tip("WD Score Threshold", "wd"),
                                min=0.0, max=10.0, value=0.0, step=0.1),
                ui.input_slider("cy_threshold_wdfdr", tip("WDFDR Threshold", "wdfdr"),
                                min=0.0, max=1.0, value=1.0, step=0.01),
                ui.div(
                    ui.tooltip(ui.input_action_button("cy_preset_stringent", "Stringent", class_="btn-sm btn-outline-primary"),
                               TOOLTIPS["preset_stringent"]),
                    ui.tooltip(ui.input_action_button("cy_preset_moderate", "Moderate", class_="btn-sm btn-outline-secondary"),
                               TOOLTIPS["preset_moderate"]),
                    ui.tooltip(ui.input_action_button("cy_preset_relaxed", "Relaxed", class_="btn-sm btn-outline-secondary"),
                               TOOLTIPS["preset_relaxed"]),
                    class_="d-flex gap-2 mb-3",
                ),
                ui.input_select("cy_labels", tip("Labels", "cy_labels"),
                                choices={"all": "All nodes", "baits": "Baits only", "none": "None"}),
                ui.input_select("cy_layout", tip("Layout", "cy_layout"),
                                choices=["force-directed", "kamada-kawai", "circular", "grid",
                                         "hierarchical", "degree-circle"]),
                ui.p("Edge style", style="font-weight: 500; margin-bottom: 5px;"),
                ui.input_select("cy_edge_width", tip("Edge width", "cy_edge_width"),
                                choices={"abundance": "Abundance (intensity / spectral counts)",
                                         "SaintScore": "SAINT score", "WD": "WD score",
                                         "FoldChange": "Fold change", "uniform": "Uniform"}),
                ui.input_checkbox("cy_prey_prey", tip("Prey-prey BioGRID edges", "cy_prey_prey"), value=True),
                ui.input_select("cy_biogrid_scope", tip("BioGRID edges", "cy_biogrid_scope"),
                                choices={"all": "All pairs", "multivalidated": "Multivalidated only"}),
                ui.input_checkbox("cy_lit_weighted", tip("Weight BioGRID edges by publications", "cy_lit_weighted"),
                                  value=False),
                ui.p("Complexes", style="font-weight: 500; margin-bottom: 5px;"),
                ui.input_checkbox("cy_corum", tip("CORUM complex edges", "cy_corum"), value=False),
                ui.layout_columns(
                    ui.input_numeric("cy_corum_min_members", tip("Min subunits drawn", "cy_corum_min_members"),
                                     value=3, min=2, step=1),
                    ui.input_numeric("cy_corum_min_fraction", tip("Min share by one bait", "cy_corum_min_fraction"),
                                     value=0.5, min=0.0, max=1.0, step=0.05),
                    col_widths=(6, 6),
                ),
                ui.div(
                    ui.tooltip(ui.input_action_button("cy_send", "Send to Cytoscape", class_="btn-primary"),
                               TOOLTIPS["cy_send"]),
                    ui.tooltip(ui.input_action_button("cy_rethreshold", "Re-apply Thresholds", class_="btn-outline-secondary"),
                               TOOLTIPS["cy_rethreshold"]),
                    ui.tooltip(ui.input_action_button("cy_restyle", "Apply Edge Style", class_="btn-outline-secondary"),
                               TOOLTIPS["cy_restyle"]),
                    class_="d-flex flex-wrap gap-2 mt-2",
                ),
            ),
            ui.card(
                ui.card_header("Cytoscape"),
                ui.output_ui("cy_status"),
                ui.input_action_button("cy_probe", "Check Connection", class_="btn-sm btn-outline-secondary"),
                ui.hr(),
                ui.p("Selection", style="font-weight: 500; margin-bottom: 5px;"),
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
                    class_="d-flex flex-wrap gap-2 mb-3",
                ),
                ui.p("Select by relation", style="font-weight: 500; margin-bottom: 5px;"),
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
                    class_="d-flex flex-wrap gap-2 mb-3",
                ),
                ui.p("Cluster the selection", style="font-weight: 500; margin-bottom: 5px;"),
                ui.layout_columns(
                    ui.input_numeric("cy_cl_resolution", tip("Resolution", "cy_cl_resolution"), value=1.0, min=0.1, max=5, step=0.1),
                    ui.input_numeric("cy_cl_seed", tip("Seed", "cy_cl_seed"), value=17, min=0, step=1),
                    ui.input_numeric("cy_cl_lit_weight", tip("Reference weight", "cy_cl_lit_weight"), value=1.0, min=0, step=0.25),
                    col_widths=(4, 4, 4),
                ),
                ui.div(
                    ui.tooltip(ui.input_action_button("cy_cluster", "Cluster and Repack", class_="btn-sm btn-outline-primary"),
                               TOOLTIPS["cy_cluster"]),
                    class_="d-flex flex-wrap gap-2 mb-3",
                ),
                ui.p("Network", style="font-weight: 500; margin-bottom: 5px;"),
                ui.div(
                    ui.tooltip(ui.input_action_button("cy_sync", "Sync Positions", class_="btn-sm btn-outline-secondary"),
                               TOOLTIPS["cy_sync"]),
                    ui.tooltip(ui.input_action_button("cy_export", "Export PNG", class_="btn-sm btn-outline-secondary"),
                               TOOLTIPS["cy_export"]),
                    ui.tooltip(ui.input_action_button("cy_unlock", "Unlock View", class_="btn-sm btn-outline-secondary"),
                               TOOLTIPS["cy_unlock"]),
                    class_="d-flex flex-wrap gap-2 mb-3",
                ),
                ui.output_data_frame("cy_selection_table"),
            ),
            ui.card(
                ui.card_header("Activity"),
                ui.output_text_verbatim("cy_activity"),
            ),
            col_widths=(4, 5, 3),
        ),
    ),
    ui.nav_panel("Downloads",
                    ui.input_select("download_dataset", "Select Dataset", choices=[]),
                    ui.output_ui("empty_state_downloads"),
                    ui.card(
                        ui.card_header(tip("Filter Data Before Download", "dl_filter_card")),
                        ui.layout_columns(
                            ui.input_slider("dl_threshold_saintscore", tip("SAINT Score (≥)", "saintscore"),
                                          min=0.0, max=1.0, value=0.0, step=0.01),
                            ui.input_slider("dl_threshold_bfdr", tip("BFDR (≤)", "bfdr"),
                                          min=0.0, max=1.0, value=1.0, step=0.01),
                            ui.input_slider("dl_threshold_wd", tip("WD Score (≥)", "wd"),
                                          min=0.0, max=10.0, value=0.0, step=0.1),
                            ui.input_slider("dl_threshold_wdfdr", tip("WDFDR (≤)", "wdfdr"),
                                          min=0.0, max=1.0, value=1.0, step=0.01),
                            col_widths=(3, 3, 3, 3)
                        ),
                        ui.layout_columns(
                            ui.tooltip(ui.input_action_button("dl_preset_stringent", "Stringent", class_="btn-sm btn-outline-primary"),
                                       TOOLTIPS["preset_stringent"]),
                            ui.tooltip(ui.input_action_button("dl_preset_moderate", "Moderate", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["preset_moderate"]),
                            ui.tooltip(ui.input_action_button("dl_preset_relaxed", "Relaxed", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["preset_relaxed"]),
                            ui.tooltip(ui.input_action_button("dl_preset_none", "No Filter", class_="btn-sm btn-outline-secondary"),
                                       TOOLTIPS["preset_none"]),
                            col_widths=(3, 3, 3, 3)
                        ),
                    ui.p("Thresholds apply only to score-based presets (Annotated Scores, "
                             "Cytoscape, Gene List, ProHits-viz, Custom). Set all thresholds to "
                             "their default values (0.0/1.0) to download unfiltered data.",
                             class_="text-muted small"),
                    ),
                    ui.layout_columns(
                        ui.card(
                            ui.card_header("Download Builder"),
                            ui.input_radio_buttons("dl_preset", tip("Preset", "dl_preset"),
                                choices={key: preset.label for key, preset in dp.PRESETS.items()}),
                            ui.input_checkbox_group("dl_groups", tip("Include", "dl_groups"), choices=[]),
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
                    )
    ),
    sidebar=ui.sidebar(
        ui.h4("ProxiMate"),
        # Resolved once: the build cannot change while the server is running, and this
        # is the identifier a bug report has to quote to be reproducible.
        ui.div(provenance.version_label(), class_="text-muted small"),
        ui.hr(),
        ui.p(
            "Your feedback is invaluable in helping us improve the tool."
        ),
        ui.p(ui.strong("Get in touch:"), style="margin-bottom: 5px;"),
        ui.tags.ul(
            ui.tags.li(
                ui.a("Report a bug", href="https://github.com/plutzer/ProxiMate/issues", target="_blank"),
            ),
            ui.tags.li(
                ui.a("GitHub Repository", href="https://github.com/plutzer/ProxiMate", target="_blank"),
            ),
        ),
        ui.hr(),
        ui.p(ui.strong("Agent access:"), style="margin-bottom: 5px;"),
        *[ui.div(ui.p(f"Connect {agent} to this server with:", class_="small", style="margin-bottom: 5px;"),
                 ui.tags.pre(command, style="white-space: pre-wrap; word-break: break-all; font-size: 0.75em;"))
          for agent, command in MCP_ADD_COMMANDS.items()],
        ui.hr(),
    ),
    title="ProxiMate",
    header=ui.output_ui("agent_banner"),
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
                      style="background: #fff3cd; color: #664d03; padding: 8px 16px; "
                            "border-bottom: 1px solid #ffe69c;")

    # Function to render the datasets table
    @render.data_frame
    def render_datasets():
        return render.DataGrid(datasets())
    
    saint_baits = reactive.Value(pd.DataFrame(
        columns=["Experiment Name", "Bait", "Type", "Bait ID"]
    ))

    ed_dataframe = reactive.Value(pd.DataFrame(
        columns=["Experiment Name", "Type", "Bait", "Replicate", "Bait ID", "Group"]
    ))

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
        if input.input_format.get() in ["MaxQuant", "DIA-NN", "Pioneer", "FragPipe", "MSstats"]:
            # Clear SAINT bait table
            saint_baits.set(pd.DataFrame(columns=["Experiment Name", "Bait", "Type", "Bait ID"]))
        elif input.input_format.get() == "SAINT":
            # Clear ED table
            ed_dataframe.set(pd.DataFrame(columns=["Experiment Name", "Type", "Bait", "Replicate", "Bait ID"]))

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

        run_id = log_config.new_run_id()
        try:
            with ui.Progress(min=0, max=1) as progress:
                backend.run_parse(dataset_name, input_format, files, input.quant_type.get(),
                                  actor='gui', ed_frame=ed_frame, bait_frame=bait_frame,
                                  run_id=run_id,
                                  progress=lambda message, value: progress.set(value, message=message))
            notify(f"Successfully parsed dataset '{dataset_name}'", type="message", duration=5)
        except ProxiMateError as e:
            notify(format_error_notification(e), type="error", duration=None)
        except (ValueError, backend.BusyError) as e:
            notify(f"Parser: {e}", type="error", duration=None)
        except FileNotFoundError as e:
            notify(f"File not found: {e}", type="error", duration=10, exc_info=True)
        except PermissionError as e:
            notify(f"Permission denied accessing file: {e}", type="error", duration=10, exc_info=True)
        except pd.errors.ParserError as e:
            notify(f"Error parsing file: {e}\n\nEnsure files are in correct format.",
                   type="error", duration=10, exc_info=True)
        except Exception as e:
            logger.exception("Unexpected error parsing dataset '%s'", dataset_name)
            notify(
                f"An unexpected error occurred while parsing '{dataset_name}':\n{e}\n\n"
                f"The full error was written to {dataset_name}/proximate.log (run {run_id}).",
                type="error", duration=None)



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


    @reactive.effect
    @reactive.event(input.score_dataset, input.imputation_method, input.pi_method)
    def update_pi_bait_choices():
        """Populate pi_bait dropdown from the selected dataset's saved ED,
        filtered to control baits with >= 3 replicates."""
        if str(input.imputation_method.get()) != "2":
            return
        if input.pi_method.get() != "single_bait":
            return
        dataset = input.score_dataset.get()
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
        dataset_name = input.score_dataset.get()
        imputation = int(input.imputation_method.get())
        pi_method = input.pi_method.get() if imputation == 2 else None
        pi_bait = input.pi_bait.get() if pi_method == 'single_bait' else None
        try:
            with ui.Progress(min=0, max=1) as progress:
                backend.run_score(dataset_name, imputation, input.wdfdr_iterations.get(),
                                  input.organism.get(), input.exclude_hcm.get(),
                                  pi_method=pi_method, pi_bait=pi_bait, actor='gui',
                                  progress=lambda message, value: progress.set(
                                      value, message="Scoring data", detail=message))
            notify(f"Successfully scored and annotated dataset '{dataset_name}'",
                   type="message", duration=5)
        except (backend.StageError, backend.BusyError, ValueError, KeyError) as e:
            notify(str(e), type="error", duration=None)
        except Exception as e:
            logger.exception("Unexpected error during scoring of '%s'", dataset_name)
            notify(f"Unexpected error during scoring: {e}", type="error", duration=None)

    # Quality controls tab
    @reactive.Calc
    def pca_matrix_cached():
        # Shared by the experiment and prey PCA plots and their PNG exports, so
        # the preprocessing runs once per settings change and all four agree.
        dataset_name = input.qc_dataset()
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
            dataset_name = input.qc_dataset.get()
            if not dataset_name:
                return None

            interaction_path = os.path.join(out_dir, dataset_name, "interaction.txt")
            ed_path = os.path.join(out_dir, dataset_name, "ED.csv")
            if not (os.path.exists(interaction_path) and os.path.exists(ed_path)):
                return None

            matrix = pca_matrix_cached()
            if matrix is None:
                return None
            fig = pca_plot(interaction_path, ed_path, matrix=matrix)

            return fig

    @render_plotly
    def prey_pca_plot():
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Generating Plots...", value=25)
            dataset_name = input.qc_dataset.get()
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


    @render_widget
    def known_retention_plot():
        # Get the selected dataset
        dataset_name = input.qc_dataset.get()
        if not dataset_name:
            return None
        
        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        ctrls = None
        if input.qc_bait.get() != "All":
            ctrls = [input.qc_bait.get()]

        if not os.path.exists(results_path):
            ctrls = None
            return None
        else:
            # Call the saint_known_retention function to generate the plot
            fig = saint_known_retention(results_path, ctrl_experiments=ctrls)

            # Return the figure widget
            return fig
        
    @render_widget
    def roc_curve_plot():
        # Get the selected dataset
        dataset_name = input.qc_dataset.get()
        if not dataset_name:
            return None
        
        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        if not os.path.exists(results_path):
            return None
        
        ctrls = None
        if input.qc_bait.get() != "All":
            ctrls = [input.qc_bait.get()]

        fig = roc_plot(results_path, known_type=input.roc_known_type.get(), ctrl_experiments=ctrls) # Add selected ctrls

        # Return the figure widget
        return fig

    @render.ui
    def metric_bait_label():
        bait_selection = input.qc_bait.get()
        if bait_selection == "All":
            label_text = "Metrics for All Baits"
        else:
            label_text = f"Metrics for {bait_selection}"
        return ui.h5(label_text, style="margin-bottom: 15px; margin-top: 0px; font-weight: 600;")

    @render.ui
    def metric_network_size():
        # Get the selected dataset
        dataset_name = input.qc_dataset.get()
        if not dataset_name:
            return ui.value_box(tip("Median Network Size", "metric_network_size"), "No data", showcase=None)

        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        if not os.path.exists(results_path):
            return ui.value_box(tip("Median Network Size", "metric_network_size"), "No data", showcase=None)

        # Get threshold values
        thresholds = {
            'SaintScore': input.threshold_saintscore.get(),
            'BFDR': input.threshold_bfdr.get(),
            'WD': input.threshold_wd.get(),
            'WDFDR': input.threshold_wdfdr.get()
        }

        # Filter by bait if not "All"
        ctrls = None
        if input.qc_bait.get() != "All":
            ctrls = [input.qc_bait.get()]

        # Call the metrics calculation function
        metrics = calculate_threshold_metrics(results_path, thresholds, ctrl_experiments=ctrls)

        return ui.value_box(
            "Median Network Size",
            f"{metrics['median_network_size']:.1f}",
            showcase=None,
            theme="primary"
        )

    @render.ui
    def metric_enrichment():
        # Get the selected dataset
        dataset_name = input.qc_dataset.get()
        if not dataset_name:
            return ui.value_box(tip("Known Enrichment", "metric_enrichment"), "No data", showcase=None)

        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        if not os.path.exists(results_path):
            return ui.value_box(tip("Known Enrichment", "metric_enrichment"), "No data", showcase=None)

        # Get threshold values
        thresholds = {
            'SaintScore': input.threshold_saintscore.get(),
            'BFDR': input.threshold_bfdr.get(),
            'WD': input.threshold_wd.get(),
            'WDFDR': input.threshold_wdfdr.get()
        }

        # Filter by bait if not "All"
        ctrls = None
        if input.qc_bait.get() != "All":
            ctrls = [input.qc_bait.get()]

        # Call the metrics calculation function
        metrics = calculate_threshold_metrics(results_path, thresholds, ctrl_experiments=ctrls)

        return ui.value_box(
            "Known Enrichment",
            f"{metrics['enrichment_ratio']:.2f}x",
            showcase=None,
            theme="success"
        )

    @render.ui
    def metric_degree():
        # Get the selected dataset
        dataset_name = input.qc_dataset.get()
        if not dataset_name:
            return ui.value_box(tip("Mean Prey-Prey Degree", "metric_degree"), "No data", showcase=None)

        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        if not os.path.exists(results_path):
            return ui.value_box(tip("Mean Prey-Prey Degree", "metric_degree"), "No data", showcase=None)

        # Get threshold values
        thresholds = {
            'SaintScore': input.threshold_saintscore.get(),
            'BFDR': input.threshold_bfdr.get(),
            'WD': input.threshold_wd.get(),
            'WDFDR': input.threshold_wdfdr.get()
        }

        # Filter by bait if not "All"
        ctrls = None
        if input.qc_bait.get() != "All":
            ctrls = [input.qc_bait.get()]

        # Call the metrics calculation function
        metrics = calculate_threshold_metrics(results_path, thresholds, ctrl_experiments=ctrls)

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

        return ui.value_box(
            "Mean Prey-Prey Degree",
            f"{metrics['mean_degree']:.1f}",
            showcase=None,
            theme="info"
        )

    # Empty state messaging for tabs
    @render.ui
    def empty_state_thresholding():
        if not scored_datasets():
            return ui.div(
                ui.h4("No Scored Datasets Available"),
                ui.p("Upload and score a dataset in the Network Scoring tab to view quality metrics and threshold controls."),
                style="text-align: center; padding: 40px; color: #666; background-color: #f8f9fa; border-radius: 8px; margin-bottom: 20px;"
            )
        return None

    @render.ui
    def empty_state_feature_analysis():
        if not scored_datasets():
            return ui.div(
                ui.h4("No Scored Datasets Available"),
                ui.p("Upload and score a dataset in the Network Scoring tab, then run feature analysis to view protein feature enrichment."),
                style="text-align: center; padding: 40px; color: #666; background-color: #f8f9fa; border-radius: 8px; margin-bottom: 20px;"
            )
        return None

    @render.ui
    def empty_state_network_comparison():
        if not scored_datasets():
            return ui.div(
                ui.h4("No Scored Datasets Available"),
                ui.p("Upload and score a dataset in the Network Scoring tab to compare protein interaction networks between baits."),
                style="text-align: center; padding: 40px; color: #666; background-color: #f8f9fa; border-radius: 8px; margin-bottom: 20px;"
            )
        return None

    @render.ui
    def empty_state_downloads():
        if not scored_datasets():
            return ui.div(
                ui.h4("No Scored Datasets Available"),
                ui.p("Upload and score a dataset in the Network Scoring tab to download filtered results."),
                style="text-align: center; padding: 40px; color: #666; background-color: #f8f9fa; border-radius: 8px; margin-bottom: 20px;"
            )
        return None

    # WDFDR warning when dataset scored with 0 iterations
    @render.ui
    def wdfdr_warning():
        dataset = input.qc_dataset.get()
        if not dataset:
            return None
        results_path = os.path.join(out_dir, dataset, "annotated_scores.csv")
        if not os.path.exists(results_path):
            return None
        try:
            df = pd.read_csv(results_path)
            if 'WDFDR' in df.columns and df['WDFDR'].isna().any():
                return ui.div(
                    ui.span("⚠ ", style="color: orange;"),
                    "WDFDR values are missing (dataset scored with 0 iterations). WDFDR threshold will not filter data.",
                    style="color: orange; font-size: 0.9em; margin-top: 5px;"
                )
        except Exception:
            # This check exists to warn about unusable results; if it cannot run,
            # silence would read as "nothing to warn about".
            logger.exception("Could not check WDFDR completeness in %s", results_path)
            return ui.div(
                ui.span("⚠ ", style="color: orange;"),
                "Could not read the results file to check WDFDR values.",
                style="color: orange; font-size: 0.9em; margin-top: 5px;"
            )
        return None

    # Threshold presets for Data Thresholding tab
    # Preset values (SaintScore, BFDR, WD, WDFDR): Stringent (0.9, 0.01, 2.0, 0.05), Moderate (0.7, 0.05, 1.0, 0.1), Relaxed (0.5, 0.1, 0.0, 1.0)
    @reactive.effect
    @reactive.event(input.qc_preset_stringent)
    def apply_qc_stringent_preset():
        ui.update_slider("threshold_saintscore", value=0.9)
        ui.update_slider("threshold_bfdr", value=0.01)
        ui.update_slider("threshold_wd", value=2.0)
        ui.update_slider("threshold_wdfdr", value=0.05)

    @reactive.effect
    @reactive.event(input.qc_preset_moderate)
    def apply_qc_moderate_preset():
        ui.update_slider("threshold_saintscore", value=0.7)
        ui.update_slider("threshold_bfdr", value=0.05)
        ui.update_slider("threshold_wd", value=1.0)
        ui.update_slider("threshold_wdfdr", value=0.1)

    @reactive.effect
    @reactive.event(input.qc_preset_relaxed)
    def apply_qc_relaxed_preset():
        ui.update_slider("threshold_saintscore", value=0.5)
        ui.update_slider("threshold_bfdr", value=0.1)
        ui.update_slider("threshold_wd", value=0.0)
        ui.update_slider("threshold_wdfdr", value=1.0)

    # Threshold presets for Protein Feature Analysis tab
    @reactive.effect
    @reactive.event(input.fa_preset_stringent)
    def apply_fa_stringent_preset():
        ui.update_slider("fa_threshold_saintscore", value=0.9)
        ui.update_slider("fa_threshold_bfdr", value=0.01)
        ui.update_slider("fa_threshold_wd", value=2.0)
        ui.update_slider("fa_threshold_wdfdr", value=0.05)

    @reactive.effect
    @reactive.event(input.fa_preset_moderate)
    def apply_fa_moderate_preset():
        ui.update_slider("fa_threshold_saintscore", value=0.7)
        ui.update_slider("fa_threshold_bfdr", value=0.05)
        ui.update_slider("fa_threshold_wd", value=1.0)
        ui.update_slider("fa_threshold_wdfdr", value=0.1)

    @reactive.effect
    @reactive.event(input.fa_preset_relaxed)
    def apply_fa_relaxed_preset():
        ui.update_slider("fa_threshold_saintscore", value=0.5)
        ui.update_slider("fa_threshold_bfdr", value=0.1)
        ui.update_slider("fa_threshold_wd", value=0.0)
        ui.update_slider("fa_threshold_wdfdr", value=1.0)

    # Threshold presets for Downloads tab
    @reactive.effect
    @reactive.event(input.dl_preset_stringent)
    def apply_dl_stringent_preset():
        ui.update_slider("dl_threshold_saintscore", value=0.9)
        ui.update_slider("dl_threshold_bfdr", value=0.01)
        ui.update_slider("dl_threshold_wd", value=2.0)
        ui.update_slider("dl_threshold_wdfdr", value=0.05)

    @reactive.effect
    @reactive.event(input.dl_preset_moderate)
    def apply_dl_moderate_preset():
        ui.update_slider("dl_threshold_saintscore", value=0.7)
        ui.update_slider("dl_threshold_bfdr", value=0.05)
        ui.update_slider("dl_threshold_wd", value=1.0)
        ui.update_slider("dl_threshold_wdfdr", value=0.1)

    @reactive.effect
    @reactive.event(input.dl_preset_relaxed)
    def apply_dl_relaxed_preset():
        ui.update_slider("dl_threshold_saintscore", value=0.5)
        ui.update_slider("dl_threshold_bfdr", value=0.1)
        ui.update_slider("dl_threshold_wd", value=0.0)
        ui.update_slider("dl_threshold_wdfdr", value=1.0)

    @reactive.effect
    @reactive.event(input.dl_preset_none)
    def apply_dl_no_filter_preset():
        ui.update_slider("dl_threshold_saintscore", value=0.0)
        ui.update_slider("dl_threshold_bfdr", value=1.0)
        ui.update_slider("dl_threshold_wd", value=0.0)
        ui.update_slider("dl_threshold_wdfdr", value=1.0)

    # Plot export download handlers (using matplotlib for PNG export)
    @render.download_button(filename="pca_plot.png")
    def download_pca_plot():
        dataset_name = input.qc_dataset.get()
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
        dataset_name = input.qc_dataset.get()
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
        dataset_name = input.qc_dataset.get()
        bait_selection = input.qc_bait.get()
        if not dataset_name or bait_selection == "All":
            notify("Select a specific bait to export the scatter plot.", type="error")
            return None
        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")
        saintscore_threshold = input.threshold_saintscore.get()
        fig = saint_scatter_matplotlib(results_path, bait_selection, saintscore_threshold)
        filepath = os.path.join(out_dir, "scatter_plot.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.download_button(filename="heatmap.png")
    def download_heatmap():
        dataset = input.feature_dataset.get()
        if not dataset:
            notify("No dataset selected.", type="error")
            return None
        feature_file = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(feature_file):
            notify("Run feature analysis first.", type="error")
            return None
        feature_data = pd.read_csv(feature_file)
        feature_type = input.feature_type.get() or 'GO_CC'
        num_features = input.num_features.get() or 30
        try:
            fig = plot_results(feature_data, feature_type, num_features)
            filepath = os.path.join(out_dir, "heatmap.png")
            fig.savefig(filepath, dpi=150, bbox_inches='tight')
            return filepath
        except ValueError as e:
            notify(f"Insufficient data to generate heatmap: {e}", type="error")
            return None

    @render.download_button(filename="volcano_plot.png")
    def download_volcano_plot():
        dataset_a = input.comp_dataset_a.get()
        dataset_b = input.comp_dataset_b.get()
        bait_a = input.comp_bait_a.get()
        bait_b = input.comp_bait_b.get()
        if not all([dataset_a, dataset_b, bait_a, bait_b]):
            notify("Select datasets and baits first.", type="error")
            return None
        if dataset_a != dataset_b:
            notify("Volcano plot requires baits from the same dataset.", type="error")
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

    @render_widget
    def saint_scatter_plot():
        # Get the selected dataset
        dataset_name = input.qc_dataset.get()
        if not dataset_name:
            return None

        # Only show for individual baits, not "All"
        bait_selection = input.qc_bait.get()
        if bait_selection == "All":
            # Return a message figure
            fig = go.Figure()
            fig.add_annotation(
                text="Please select a specific bait to view this plot",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            fig.update_layout(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False)
            )
            return fig

        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        if not os.path.exists(results_path):
            return None

        # Get threshold value for reference line
        saintscore_threshold = input.threshold_saintscore.get()

        # Call the scatter plot function
        fig = plot_saint_scatter(results_path, bait_selection, saintscore_threshold)

        return fig

    feature_enrichment = reactive.Value(pd.DataFrame())


    # Feature analysis tab
    @reactive.effect
    @reactive.event(input.feature_analysis)
    def feature_analysis():
        dataset_name = input.feature_dataset.get()
        try:
            with ui.Progress(min=0, max=100) as progress, \
                    log_config.dataset_log(os.path.join(out_dir, dataset_name)):
                progress.set(message="Running protein feature analysis", value=5)
                thresholds = {
                    'SaintScore': input.fa_threshold_saintscore.get(),
                    'BFDR': input.fa_threshold_bfdr.get(),
                    'WD': input.fa_threshold_wd.get(),
                    'WDFDR': input.fa_threshold_wdfdr.get()
                }
                logger.info("Starting feature enrichment for '%s' (thresholds=%s)",
                            dataset_name, thresholds)

                # Get the selected dataset
                progress.set(message="Running protein feature analysis", detail="Loading dataset...", value=20)
                dataset = pd.read_csv(os.path.join(out_dir, dataset_name, "annotated_scores.csv"))

                progress.set(message="Running protein feature analysis", detail="Processing data...", value=50)
                result = process_refactored(
                    dataset,
                    columns_for_analysis = ['GO_CC', 'GO_BP', 'GO_MF', 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains'],
                    thresholds = thresholds
                )

                progress.set(message="Running protein feature analysis", detail="Generating plots...", value=80)
                # Store the results in the reactive value
                feature_enrichment.set(result)

                progress.set(message="Running protein feature analysis", detail="Saving results...", value=90)
                # Save the results to the dataset directory
                output_path = os.path.join(out_dir, dataset_name, "Feature_enrichment.csv")
                result.to_csv(output_path, index=False)
                logger.info("Feature enrichment written to %s (%d rows)", output_path, len(result))

                progress.set(message="Running protein feature analysis", detail="Done!", value=100)
        except Exception as e:
            notify(f"Feature analysis failed for '{dataset_name}': {e}",
                   type="error", duration=None, exc_info=True)

    @render.plot
    def feature_enrichment_plot():
        # Trigger re-render when feature analysis completes
        _ = feature_enrichment.get()

        # Set the feature enrichment to whatever dataset is selected
        dataset = input.feature_dataset.get()
        if not dataset:
            return None

        feature_file = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(feature_file):
            return None

        # Load the feature enrichment data from file
        feature_data = pd.read_csv(feature_file)

        feature_type = input.feature_type.get()
        if not feature_type:
            feature_type = 'GO_CC'  # Default feature type if none is selected
        # Plot the results for a specific feature type, e.g., 'Domains'

        num_features = input.num_features.get() if input.num_features.get() else 30  # Default to 30 if not set
        if num_features < 1:
            num_features = 30

        # Call the plot_results function to generate the heatmap
        try:
            heatmap = plot_results(feature_data, feature_type, num_features)
            return heatmap
        except ValueError as e:
            # Handle case where there aren't enough features to cluster
            logger.warning("Cannot plot %s enrichment for '%s': %s",
                           feature_type, dataset, e)
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, f'Insufficient data to generate plot for {feature_type}\n\nTry selecting a different feature type or lowering the SAINT threshold.',
                    ha='center', va='center', fontsize=14, wrap=True)
            ax.axis('off')
            return fig

    @reactive.effect
    @reactive.event(input.feature_dataset, feature_enrichment)
    def update_download_bait_filter():
        """Update the bait dropdown for enrichment download based on available data."""
        dataset = input.feature_dataset.get()
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
        dataset = input.feature_dataset.get()
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

        # Save to temp file and return
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"enrichment_{dataset}_{timestamp}.csv"
        savepath = os.path.join(out_dir, filename)
        filtered_data.to_csv(savepath, index=False)
        return savepath

    @reactive.Calc
    def available_datasets():
        datasets()
        return dataset_store.names()

    @reactive.Calc
    def scored_datasets():
        datasets()
        return dataset_store.scored_names()

    @reactive.Effect
    @reactive.event(datasets)
    def update_score_dataset():
        # Update the dropdown choices dynamically
        available_choices = available_datasets()
        scored_choices = scored_datasets()
        ui.update_select("score_dataset", choices=available_choices)
        ui.update_select("qc_dataset", choices=available_choices)
        ui.update_select("download_dataset", choices=available_choices)
        ui.update_select("feature_dataset", choices=scored_choices)
        ui.update_select("comp_dataset_a", choices=scored_choices)
        ui.update_select("comp_dataset_b", choices=scored_choices)
        ui.update_select("cy_dataset", choices=scored_choices)

    @reactive.effect
    @reactive.event(input.qc_dataset, datasets)
    def update_qc_bait():
        # Update the bait dropdown choices dynamically based on the selected dataset
        dataset_name = input.qc_dataset()
        if dataset_name:
            # Check to see if the dataset has been scored
            try:
                scores = pd.read_csv(os.path.join(out_dir, dataset_name, "annotated_scores.csv"))
                baits = scores['Experiment.ID'].unique().tolist()
                baits.insert(0, "All")  # Add "All" option
                ui.update_select("qc_bait", choices=baits)
            except FileNotFoundError:
                # Expected before a dataset has been scored.
                logger.debug("No annotated_scores.csv for %s yet", dataset_name)
                ui.update_select("qc_bait", choices=["All"])  # Reset to default if file not found
            except Exception:
                logger.exception("Could not read baits for dataset %s", dataset_name)
                ui.update_select("qc_bait", choices=["All"])
        else:
            ui.update_select("qc_bait", choices=["All"])  # Reset to default if no dataset is selected

    @reactive.effect
    @reactive.event(input.qc_dataset, datasets)
    def update_prey_pca_color_choices():
        # Annotation-based color options exist only once the dataset is scored;
        # HPA localization only when the organism has HPA data (human).
        choices = {"none": "None", "detection": "Detection count"}
        dataset_name = input.qc_dataset()
        baits = []
        if dataset_name:
            scores_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")
            if os.path.exists(scores_path):
                try:
                    cols = pd.read_csv(scores_path, nrows=0).columns
                    choices["saint"] = "SaintScore (per bait)"
                    if "Main location" in cols:
                        choices["hpa"] = "HPA Main location"
                    if "first_SCL" in cols:
                        choices["scl"] = "UniProt localization"
                    baits = pd.read_csv(scores_path, usecols=['Experiment.ID'])['Experiment.ID'].unique().tolist()
                except Exception:
                    logger.exception("Could not read %s for prey PCA color options", scores_path)
        ui.update_select("prey_pca_color", choices=choices)
        ui.update_select("prey_pca_bait", choices=baits)

    @reactive.effect
    @reactive.event(input.qc_dataset, datasets)
    def update_pca_normalization_default():
        # Intensity-scale data defaults to log2; spectral counts stay linear.
        # Fires only on dataset change, so a manual override sticks after it.
        dataset_name = input.qc_dataset()
        if not dataset_name:
            return
        df = datasets()
        quant = df.loc[df['Dataset Name'] == dataset_name, 'Quant Type']
        if quant.empty:
            return
        default = "log2_zscore" if quant.values[0] in ("Intensity", "LFQ") else "zscore"
        ui.update_select("pca_normalization", selected=default)

    # Network Comparison Tab - Reactive Effects and Renderers

    # Cached data loaders for Network Comparison
    @reactive.Calc
    def comp_volcano_data_cached():
        """Cache volcano data calculation to avoid recomputing."""
        # Get all inputs
        dataset_a = input.comp_dataset_a()
        dataset_b = input.comp_dataset_b()
        bait_a = input.comp_bait_a()
        bait_b = input.comp_bait_b()

        # Return empty if not all selected or datasets don't match
        if not all([dataset_a, dataset_b, bait_a, bait_b]) or dataset_a != dataset_b:
            return pd.DataFrame()

        # Get thresholds
        thresholds_a = {
            'SaintScore': input.comp_saintscore_a(),
            'BFDR': input.comp_bfdr_a(),
            'WD': input.comp_wd_a(),
            'WDFDR': input.comp_wdfdr_a()
        }
        thresholds_b = {
            'SaintScore': input.comp_saintscore_b(),
            'BFDR': input.comp_bfdr_b(),
            'WD': input.comp_wd_b(),
            'WDFDR': input.comp_wdfdr_b()
        }

        # Calculate and return volcano data
        return calculate_volcano_data(
            dataset_a, bait_a, bait_b,
            thresholds_a, thresholds_b,
            out_dir
        )

    @reactive.Calc
    def comp_filtered_data_a_cached():
        """Cache filtered data for bait A."""
        dataset_a = input.comp_dataset_a()
        bait_a = input.comp_bait_a()

        if not dataset_a or not bait_a:
            return pd.DataFrame()

        thresholds_a = {
            'SaintScore': input.comp_saintscore_a(),
            'BFDR': input.comp_bfdr_a(),
            'WD': input.comp_wd_a(),
            'WDFDR': input.comp_wdfdr_a()
        }

        return load_and_filter_bait_data(dataset_a, bait_a, thresholds_a, out_dir)

    @reactive.Calc
    def comp_filtered_data_b_cached():
        """Cache filtered data for bait B."""
        dataset_b = input.comp_dataset_b()
        bait_b = input.comp_bait_b()

        if not dataset_b or not bait_b:
            return pd.DataFrame()

        thresholds_b = {
            'SaintScore': input.comp_saintscore_b(),
            'BFDR': input.comp_bfdr_b(),
            'WD': input.comp_wd_b(),
            'WDFDR': input.comp_wdfdr_b()
        }

        return load_and_filter_bait_data(dataset_b, bait_b, thresholds_b, out_dir)

    @reactive.effect
    @reactive.event(input.comp_dataset_a, datasets)
    def update_comp_bait_a():
        """Update bait A dropdown based on selected dataset."""
        dataset_name = input.comp_dataset_a()
        if dataset_name:
            try:
                scores = pd.read_csv(os.path.join(out_dir, dataset_name, "annotated_scores.csv"))
                baits = scores['Experiment.ID'].unique().tolist()
                ui.update_select("comp_bait_a", choices=baits)
            except FileNotFoundError:
                # Expected before a dataset has been scored.
                logger.debug("No annotated_scores.csv for %s yet", dataset_name)
                ui.update_select("comp_bait_a", choices=[])
            except Exception:
                logger.exception("Could not read baits for dataset %s", dataset_name)
                ui.update_select("comp_bait_a", choices=[])
        else:
            ui.update_select("comp_bait_a", choices=[])

    @reactive.effect
    @reactive.event(input.comp_dataset_b, datasets)
    def update_comp_bait_b():
        """Update bait B dropdown based on selected dataset."""
        dataset_name = input.comp_dataset_b()
        if dataset_name:
            try:
                scores = pd.read_csv(os.path.join(out_dir, dataset_name, "annotated_scores.csv"))
                baits = scores['Experiment.ID'].unique().tolist()
                ui.update_select("comp_bait_b", choices=baits)
            except FileNotFoundError:
                # Expected before a dataset has been scored.
                logger.debug("No annotated_scores.csv for %s yet", dataset_name)
                ui.update_select("comp_bait_b", choices=[])
            except Exception:
                logger.exception("Could not read baits for dataset %s", dataset_name)
                ui.update_select("comp_bait_b", choices=[])
        else:
            ui.update_select("comp_bait_b", choices=[])

    @render_widget
    @reactive.event(input.compare_networks)
    def volcano_plot():
        """Render volcano plot comparing two baits."""
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Comparing networks...", detail="Initializing", value=0)

            # Get selections
            dataset_a = input.comp_dataset_a()
            dataset_b = input.comp_dataset_b()
            bait_a = input.comp_bait_a()
            bait_b = input.comp_bait_b()

            # Validate inputs
            if not all([dataset_a, dataset_b, bait_a, bait_b]):
                fig = go.Figure()
                fig.add_annotation(
                    text="Please select datasets and baits for comparison",
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

            # Check if datasets match
            if dataset_a != dataset_b:
                fig = go.Figure()
                fig.add_annotation(
                    text="Volcano plot only available when comparing baits from the same dataset.<br>Please select the same dataset for both networks.",
                    xref="paper", yref="paper",
                    x=0.5, y=0.5, showarrow=False,
                    font=dict(size=14, color="orange")
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

        # Get selections
        dataset_a = input.comp_dataset_a()
        dataset_b = input.comp_dataset_b()
        bait_a = input.comp_bait_a()
        bait_b = input.comp_bait_b()

        # Validate inputs
        if not all([dataset_a, dataset_b, bait_a, bait_b]):
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.5, 0.5, "Select datasets and baits to view overlap",
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
        if not all([input.comp_dataset_a(), input.comp_dataset_b(), input.comp_bait_a(), input.comp_bait_b()]):
            return "Select datasets and baits to view gene lists"
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

    @reactive.Effect
    @reactive.event(input.download_dataset, datasets)
    def update_dl_presets():
        """Offer only the presets whose required files exist for this dataset."""
        dataset_name = input.download_dataset.get()
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

    @reactive.Effect
    @reactive.event(input.dl_preset, input.download_dataset)
    def update_dl_groups():
        """Populate the Include checkboxes with the preset's column groups or
        file choices, preserving still-valid picks across switches."""
        preset = dp.PRESETS.get(input.dl_preset.get())
        items = () if preset is None else (
            preset.files if preset.kind == "files" else preset.groups)
        if not items:
            ui.update_checkbox_group("dl_groups", choices=[], selected=[])
            return
        choices = {item.key: item.label for item in items}
        current = [k for k in input.dl_groups.get() if k in choices]
        selected = current or [item.key for item in items if item.default]
        ui.update_checkbox_group("dl_groups", choices=choices, selected=selected)

    @reactive.Effect
    @reactive.event(input.download_dataset, datasets)
    def update_selectize_custom_columns():
        """Offer the selected dataset's full column list in the Custom preset,
        keeping the user's still-valid picks. Client-side options only: the
        column list is small and server-side selectize never delivers options
        to the browser in this app."""
        dataset_name = input.download_dataset.get()
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
        """Cache the annotated_scores.csv file for the selected dataset in Downloads tab.
        Only re-reads when the dataset selection changes, not on slider changes."""
        dataset = input.download_dataset.get()
        if not dataset:
            return pd.DataFrame()
        dataset_path = os.path.join(out_dir, dataset, "annotated_scores.csv")
        if not os.path.exists(dataset_path):
            return pd.DataFrame()
        return pd.read_csv(dataset_path)

    @reactive.Calc
    def cached_enrichment_data():
        """Cache Feature_enrichment.csv for the selected dataset; empty frame
        when the Feature Analysis tab has not produced one."""
        dataset = input.download_dataset.get()
        if not dataset:
            return pd.DataFrame()
        path = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(path):
            return pd.DataFrame()
        return pd.read_csv(path)

    def _dl_thresholds():
        return {
            'SaintScore': input.dl_threshold_saintscore(),
            'BFDR': input.dl_threshold_bfdr(),
            'WD': input.dl_threshold_wd(),
            'WDFDR': input.dl_threshold_wdfdr()
        }

    @render.data_frame
    def dl_preview_table():
        """Build the selected preset's export and preview it. The result is
        stored in dl_result / dl_result_files so the download handler writes
        exactly what is previewed."""
        dataset = input.download_dataset.get()
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
        dataset = input.download_dataset.get()
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
        dataset = input.download_dataset.get()
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

    def cy_thresholds():
        return {'SaintScore': input.cy_threshold_saintscore.get(),
                'BFDR': input.cy_threshold_bfdr.get(),
                'WD': input.cy_threshold_wd.get(),
                'WDFDR': input.cy_threshold_wdfdr.get()}

    @render.ui
    def empty_state_cytoscape():
        if not scored_datasets():
            return ui.div(
                ui.h4("No Scored Datasets Available"),
                ui.p("Score a dataset in the Network Scoring tab to send its network to Cytoscape."),
                style="text-align: center; padding: 40px; color: #666; background-color: #f8f9fa; border-radius: 8px; margin-bottom: 20px;"
            )
        return None

    @reactive.effect
    @reactive.event(input.cy_dataset, datasets)
    def update_cy_baits():
        dataset_name = input.cy_dataset.get()
        baits = []
        if dataset_name:
            try:
                scores = pd.read_csv(os.path.join(out_dir, dataset_name, "annotated_scores.csv"),
                                     usecols=['Experiment.ID'])
                baits = sorted(scores['Experiment.ID'].astype(str).unique())
            except Exception:
                logger.exception("Could not list baits for %s", dataset_name)
        ui.update_selectize("cy_baits", choices=baits, selected=[])

    for _key, _values in (("stringent", (0.9, 0.01, 2.0, 0.05)),
                          ("moderate", (0.7, 0.05, 1.0, 0.1)),
                          ("relaxed", (0.5, 0.1, 0.0, 1.0))):
        def _make_preset(values, key):
            @reactive.effect
            @reactive.event(getattr(input, f"cy_preset_{key}"))
            def _apply():
                for name, value in zip(("saintscore", "bfdr", "wd", "wdfdr"), values):
                    ui.update_slider(f"cy_threshold_{name}", value=value)
        _make_preset(_values, _key)

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
            line = ui.p(ui.span("● ", style="color: green;"),
                        f"Cytoscape {health['version']} at {health['url']}")
        else:
            line = ui.p(ui.span("● ", style="color: red;"),
                        f"No Cytoscape at {health['url']}: {health['error']}",
                        ui.br(), "Start Cytoscape on this machine, or set PROXIMATE_CYTOSCAPE_URL.",
                        style="color: #a33;")
        if snap['net_suid'] is None:
            drawn = ui.p("No ProxiMate network drawn yet.", style="color: #666;")
        else:
            drawn = ui.p(f"{snap['title']}: {snap['n_nodes']} nodes, {snap['n_edges']} edges "
                         f"({snap['n_hidden']} hidden)" + (f" — {snap['busy']}" if snap['busy'] else ""))
        return ui.div(line, drawn)

    @reactive.effect
    @reactive.event(input.cy_send)
    def cy_send():
        dataset_name = input.cy_dataset.get()
        if not dataset_name:
            notify("Select a scored dataset first.", type="error")
            return
        dataset_path = os.path.join(out_dir, dataset_name)
        settings = provenance.annotation_settings(dataset_path)
        organism = settings["organism"]
        biogrid_path = provenance.biogrid_summary_path(
            organism, exclude_hcm=settings["exclude_hcm"])
        corum_path = None
        if input.cy_corum.get():
            from setup_datasets import CORUM_FILENAME, ORGANISMS
            if ORGANISMS[organism]["has_corum"]:
                corum_path = os.path.join(provenance.DEFAULT_DATASETS_DIR, CORUM_FILENAME)
            else:
                notify(f"CORUM covers human complexes only; drawing the {organism} network without them.",
                       type="warning")
        snap = cy_call("send the network to Cytoscape", cytoscape_ctl.draw,
                       dataset_name, os.path.join(dataset_path, "annotated_scores.csv"),
                       cy_thresholds(), baits=list(input.cy_baits.get() or []),
                       prey_prey=input.cy_prey_prey.get(), biogrid_path=biogrid_path,
                       label_policy=input.cy_labels.get(), layout=input.cy_layout.get(),
                       width_source=input.cy_edge_width.get(),
                       literature_weighted=input.cy_lit_weighted.get(),
                       biogrid_scope=input.cy_biogrid_scope.get(), corum_path=corum_path,
                       corum_min_members=int(input.cy_corum_min_members.get() or 3),
                       corum_min_fraction=float(input.cy_corum_min_fraction.get() or 0.0))
        if snap:
            notify(f"Drew {snap['n_nodes']} nodes and {snap['n_edges']} edges in Cytoscape.")

    @reactive.effect
    @reactive.event(input.cy_rethreshold)
    def cy_rethreshold():
        hidden = cy_call("re-apply the thresholds", cytoscape_ctl.apply_thresholds, cy_thresholds())
        if hidden is not None:
            notify(f"Thresholds applied: {hidden} edge(s) hidden.")

    @reactive.effect
    @reactive.event(input.cy_restyle)
    def cy_restyle():
        changed = cy_call("apply the edge style", cytoscape_ctl.restyle_edges,
                          input.cy_edge_width.get(), input.cy_lit_weighted.get(),
                          input.cy_biogrid_scope.get())
        if changed is not None:
            notify(f"Edge style applied: {changed} edge(s) changed.")

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
                         resolution=float(resolution), seed=int(seed), literature_weight=float(weight))
        if result:
            notify(f"{result['n']} nodes clustered into {result['n_communities']} communities "
                   f"(sizes {result['sizes']}).")

    @reactive.effect
    @reactive.event(input.cy_read_selection)
    def cy_read_selection():
        result = cy_call("read the selection", cytoscape_ctl.read_selection)
        if result:
            chosen, detail = result
            cy_selection.set(detail if len(detail) else chosen)
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

    @reactive.effect
    @reactive.event(input.cy_sync)
    def cy_sync():
        n = cy_call("read node positions", cytoscape_ctl.sync_positions)
        if n is not None:
            notify(f"Recorded positions of {n} node(s).")

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
        return render.DataGrid(cy_selection.get(), height="250px")

    @render.text
    def cy_activity():
        cy_version()
        entries = cytoscape_ctl.snapshot()['log']
        if not entries:
            return "No Cytoscape activity yet."
        return "\n".join(f"{e['ts'][11:]}  [{e['actor']}] {e['op']}: {e['detail']}"
                         for e in reversed(entries))

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
