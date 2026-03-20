from functools import partial
from shiny import App, Inputs, Outputs, Session, reactive, render, ui, run_app
import plotly.graph_objects as go
import plotly.express as px
from shinywidgets import output_widget, render_widget, render_plotly
import pandas as pd
import os
import sys
sys.path.append('/Scripts')
from ed_exceptions import ProxiMateError
import parse
import subprocess
import zipfile
import tempfile
import datetime
import shutil
from QC_plots import pca_plot, saint_known_retention, roc_plot, saint_scatter_plot as plot_saint_scatter, calculate_threshold_metrics, apply_score_thresholds
from Ann_Enrichment import process_refactored, plot_results
from network_comparison import (
    load_and_filter_bait_data,
    calculate_volcano_data,
    create_volcano_plot,
    create_venn_diagram,
    create_venn_diagram_matplotlib,
    create_volcano_plot_matplotlib
)
from plot_exports import pca_plot_matplotlib, saint_scatter_matplotlib
import py4cytoscape as p4c


out_dir = "/Outputs"

app_ui = ui.page_navbar(
    ui.nav_spacer(),
    ui.nav_panel("Network Scoring",
                 ui.layout_columns(
                    ui.card(
                        ui.card_header("Data Parsing"),
                        ui.layout_columns(
                            ui.input_text("dataset_name",
                                        "Dataset Name",
                                        placeholder="No spaces or special characters (/ \\ : * ? \" < > |)"),
                            ui.input_select("input_format", "Input Format",
                                           choices=["MaxQuant", "DIA-NN", "SAINT"],
                                           selected="MaxQuant"),
                            col_widths=[4, 8]
                        ),
                        ui.panel_conditional(
                            "input.input_format === 'MaxQuant'",
                            ui.layout_columns(
                                ui.div(
                                    ui.input_file("pg_file", "MaxQuant proteinGroups.txt file"),
                                    ui.tooltip(
                                        ui.input_file("ed_file", "Experimental Design File"),
                                        "ED Format: Experiment Name, Type, Bait, Replicate, Bait ID"
                                    ),
                                    ui.input_select("quant_type", "Quantification Type",
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
                                    ui.input_file("diann_matrix_file", "DIA-NN report.pg_matrix.tsv file"),
                                    ui.tooltip(
                                        ui.input_file("ed_file", "Experimental Design File"),
                                        "ED Format: Experiment Name, Type, Bait, Replicate, Bait ID"
                                    )
                                ),
                                ui.output_data_frame("ed_table_diann"),
                                col_widths=[4, 8]
                            )
                        ),
                        ui.panel_conditional(
                            "input.input_format === 'SAINT'",
                            ui.layout_columns(
                                ui.div(
                                    ui.input_file("bait", "SAINT bait.txt file"),
                                    ui.input_file("prey", "SAINT prey.txt file"),
                                    ui.input_file("interaction", "SAINT interaction.txt file"),
                                    ui.input_select("quant_type", "Quantification Type",
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
                        ui.input_select("organism", "Organism",
                            choices={"human": "Human (H. sapiens)",
                                     "mouse": "Mouse (M. musculus)",
                                     "yeast": "Yeast (S. cerevisiae)"},
                            selected="human"),
                        ui.input_radio_buttons("imputation_method", "Imputation Method",
                                              choices={0: "Default", 1: "Prey-specific"}),
                        ui.input_numeric("wdfdr_iterations", "WDFDR Iterations", value=1000),
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
                                ui.input_action_button("clear_datasets", "Clear All Datasets"),
                                ui.download_button("download_session", "Download Session Zip"),
                            ),
                            ui.card(
                                ui.input_file("session_file", "Upload Session Zip"),
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

                    # Row 1: PCA and Threshold controls (50% width each)
                    ui.layout_columns(
                        ui.card(
                            ui.card_header("Raw Data PCA"),
                            output_widget("raw_pca_plot"),
                            ui.download_button("download_pca_plot", "Export PNG", class_="btn-sm"),
                        ),
                        ui.card(
                            ui.card_header("Threshold Settings"),
                            ui.input_select("qc_bait", "Select Control Bait", choices=["All"]),
                            ui.input_slider("threshold_saintscore", "SAINT Score Threshold",
                                          min=0.0, max=1.0, value=0.7, step=0.01),
                            ui.input_slider("threshold_bfdr", "BFDR Threshold",
                                          min=0.0, max=1.0, value=0.05, step=0.01),
                            ui.input_slider("threshold_wd", "WD Score Threshold",
                                          min=0.0, max=10.0, value=0.0, step=0.1),
                            ui.input_slider("threshold_wdfdr", "WDFDR Threshold",
                                          min=0.0, max=1.0, value=0.05, step=0.01),
                            ui.p("Presets:", style="margin-top: 15px; margin-bottom: 5px; font-weight: 500;"),
                            ui.layout_columns(
                                ui.input_action_button("qc_preset_stringent", "Stringent", class_="btn-sm btn-outline-primary"),
                                ui.input_action_button("qc_preset_moderate", "Moderate", class_="btn-sm btn-outline-secondary"),
                                ui.input_action_button("qc_preset_relaxed", "Relaxed", class_="btn-sm btn-outline-secondary"),
                                col_widths=(4, 4, 4)
                            ),
                            ui.p("Note: Thresholds are shown as reference lines on plots. Data is not filtered.",
                                 style="font-style: italic; color: #666; margin-top: 10px;"),
                            ui.output_ui("wdfdr_warning"),
                        ),
                        col_widths=(6, 6),
                    ),

                    # Row 2: SAINT scatter plot and metrics side-by-side
                    ui.layout_columns(
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
                        col_widths=(6, 6),
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
                        ui.input_slider("saint_threshold", "Saint Threshold", min=0.0, max=1.0, value=0.9),
                        ui.input_action_button("feature_analysis", "Run Feature Analysis"),
                    ),
                    ui.card(
                        ui.card_header("Feature Enrichment Analysis"),
                        ui.input_select("feature_type", "Select Feature Type", choices=["GO_CC", "GO_BP", "GO_MF", "Motifs", "Regions", "Repeats", "Compositions", "Domains"]),
                        ui.input_numeric("num_features", "Number of Features to Display", value=30, min=1, max=100),
                        ui.output_plot("feature_enrichment_plot"),
                        ui.download_button("download_heatmap", "Export Heatmap PNG", class_="btn-sm"),
                        ui.hr(),
                        ui.h5("Download Enrichment Results"),
                        ui.layout_columns(
                            ui.input_select("download_feature_type", "Feature Type",
                                          choices=["All", "GO_CC", "GO_BP", "GO_MF", "Motifs", "Regions", "Repeats", "Compositions", "Domains"]),
                            ui.input_select("download_bait_filter", "Bait", choices=["All"]),
                            col_widths=(6, 6)
                        ),
                        ui.layout_columns(
                            ui.input_slider("download_pvalue_threshold", "Max Adjusted p-value",
                                          min=0.0, max=1.0, value=0.05, step=0.01),
                            ui.input_slider("download_enrichment_threshold", "Min Enrichment",
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
                ui.input_select("comp_bait_a", "Bait A", choices=[]),
                ui.input_slider("comp_saintscore_a", "SAINT Score (≥)",
                              min=0.0, max=1.0, value=0.7, step=0.01),
                ui.input_slider("comp_bfdr_a", "BFDR (≤)",
                              min=0.0, max=1.0, value=0.05, step=0.01),
                ui.input_slider("comp_wd_a", "WD Score (≥)",
                              min=0.0, max=10.0, value=0.0, step=0.1),
                ui.input_slider("comp_wdfdr_a", "WDFDR (≤)",
                              min=0.0, max=1.0, value=0.05, step=0.01),
            ),
            # Bait B selector card
            ui.card(
                ui.card_header("Network B"),
                ui.input_select("comp_dataset_b", "Dataset B", choices=[]),
                ui.input_select("comp_bait_b", "Bait B", choices=[]),
                ui.input_slider("comp_saintscore_b", "SAINT Score (≥)",
                              min=0.0, max=1.0, value=0.7, step=0.01),
                ui.input_slider("comp_bfdr_b", "BFDR (≤)",
                              min=0.0, max=1.0, value=0.05, step=0.01),
                ui.input_slider("comp_wd_b", "WD Score (≥)",
                              min=0.0, max=10.0, value=0.0, step=0.1),
                ui.input_slider("comp_wdfdr_b", "WDFDR (≤)",
                              min=0.0, max=1.0, value=0.05, step=0.01),
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
            ui.card_header("Differential Abundance Volcano Plot"),
            output_widget("volcano_plot"),
            ui.download_button("download_volcano_plot", "Export PNG", class_="btn-sm"),
            ui.p("Volcano plot only shown when baits are from the same dataset.",
                 style="font-style: italic; color: #666; margin-top: 10px;"),
        ),

        # Bottom section: Venn diagram and gene lists
        ui.layout_columns(
            ui.card(
                ui.card_header("Network Overlap"),
                ui.output_plot("venn_diagram"),
                ui.download_button("download_venn_diagram", "Export PNG", class_="btn-sm"),
            ),
            ui.card(
                ui.card_header("Gene Lists"),
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
        ui.layout_columns(
            ui.card(
                ui.card_header("Cytoscape Connection Test"),
                ui.p("This tab tests communication between ProxiMate and Cytoscape via py4cytoscape."),
                ui.p("Requirements:", style="font-weight: bold; margin-top: 15px;"),
                ui.tags.ul(
                    ui.tags.li("Cytoscape must be running on your host machine"),
                    ui.tags.li("Windows/Mac Docker Desktop: Use standard docker run -p 3838:3838"),
                    ui.tags.li("Native Linux: Use docker run --network host"),
                ),
                ui.hr(),
                ui.input_action_button("test_create_node", "Create Test Node", class_="btn-primary"),
                ui.input_action_button("test_delete_node", "Delete Test Node", class_="btn-danger",
                                      style="margin-left: 10px;"),
                ui.hr(),
                ui.output_text_verbatim("cytoscape_status"),
            ),
            col_widths=(12,),
        ),
    ),
    ui.nav_panel("Downloads",
                    ui.input_select("download_dataset", "Select Dataset", choices=[]),
                    ui.output_ui("empty_state_downloads"),
                    ui.card(
                        ui.card_header("Filter Data Before Download"),
                        ui.layout_columns(
                            ui.input_slider("dl_threshold_saintscore", "SAINT Score (≥)",
                                          min=0.0, max=1.0, value=0.0, step=0.01),
                            ui.input_slider("dl_threshold_bfdr", "BFDR (≤)",
                                          min=0.0, max=1.0, value=1.0, step=0.01),
                            ui.input_slider("dl_threshold_wd", "WD Score (≥)",
                                          min=0.0, max=10.0, value=0.0, step=0.1),
                            ui.input_slider("dl_threshold_wdfdr", "WDFDR (≤)",
                                          min=0.0, max=1.0, value=1.0, step=0.01),
                            col_widths=(3, 3, 3, 3)
                        ),
                        ui.layout_columns(
                            ui.input_action_button("dl_preset_stringent", "Stringent", class_="btn-sm btn-outline-primary"),
                            ui.input_action_button("dl_preset_moderate", "Moderate", class_="btn-sm btn-outline-secondary"),
                            ui.input_action_button("dl_preset_relaxed", "Relaxed", class_="btn-sm btn-outline-secondary"),
                            ui.input_action_button("dl_preset_none", "No Filter", class_="btn-sm btn-outline-secondary"),
                            col_widths=(3, 3, 3, 3)
                        ),
                        ui.p("Set all thresholds to their default values (0.0/1.0) to download unfiltered data.",
                             style="font-style: italic; color: #666;"),
                    ),
                    ui.layout_columns(
                        ui.card(
                            ui.card_header("Create a Custom Dataset"),
                            "Select Columns for the Custom Dataset",
                            ui.input_selectize("custom_columns", "Select Columns", choices=["Experiment.ID", "Prey.ID", "SaintScore", "BFDR"], multiple=True, selected=["Experiment.ID", "Prey.ID", "SaintScore", "BFDR"]),
                        ),
                        ui.card(
                            ui.card_header("Custom Dataset"),
                            ui.output_text("download_row_count"),
                            ui.output_data_frame("custom_table"),
                            ui.download_button("download_custom_dataset", "Download Custom Dataset"),
                        ),
                        col_widths=(4,8)
                    ),
                    ui.card(
                        ui.card_header("Batch Export"),
                        ui.p("Download all results for a dataset as a ZIP file. Includes merged.csv, annotated_scores.csv, Feature_enrichment.csv (if available), and SAINT input files."),
                        ui.download_button("download_batch", "Download All Results (ZIP)"),
                    )
    ),
    sidebar=ui.sidebar(
        ui.h4("ProxiMate Beta"),
        ui.hr(),
        ui.p(
            "Welcome! This is a ",
            ui.strong("pre-release beta version"),
            " of ProxiMate.",
        ),
        ui.p(
            "Features may change, and you may encounter bugs. "
            "Your feedback is invaluable in helping us improve the tool."
        ),
        ui.hr(),
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
        ui.p(
            "Thank you for testing ProxiMate!",
            style="font-style: italic; color: #666;"
        ),
    ),
    title="ProxiMate",
)


def get_cytoscape_base_url():
    """
    Get the Cytoscape base URL based on platform.

    For Docker containers:
    - Docker Desktop (Windows/Mac): host.docker.internal:1234
    - Native Linux: Use --network host and localhost:1234
    """
    # Check if running in Docker by looking for /.dockerenv
    in_docker = os.path.exists('/.dockerenv')

    if in_docker:
        # Running in Docker container
        # Try host.docker.internal first (works on Docker Desktop for Windows/Mac)
        # On native Linux, this won't resolve unless using --add-host or --network host
        return "http://host.docker.internal:1234/v1"
    else:
        # Running natively (for local development)
        return "http://127.0.0.1:1234/v1"


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
    datasets = reactive.Value(pd.DataFrame(
        columns=['Dataset Name', 'Input Type', 'Quant Type', 'Experiments', 'Controls', 'Scored', 'Imputation', 'WDFDR iterations']
    ))

    # Function to render the datasets table
    @render.data_frame
    def render_datasets():
        return render.DataGrid(datasets.get())
    
    saint_baits = reactive.Value(pd.DataFrame(
        columns=["Experiment Name", "Bait", "Type", "Bait ID"]
    ))

    ed_dataframe = reactive.Value(pd.DataFrame(
        columns=["Experiment Name", "Type", "Bait", "Replicate", "Bait ID"]
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

    @reactive.effect
    @reactive.event(input.bait)
    def update_bait_table():
        # Check to see if the bait file has been uploaded
        if input.bait.get():
            # Read the bait file and set it to the saint_baits reactive value
            # Read the bait file into a DataFrame

            # Catch bad input files here

            # Check to make sure the columns are correct and don't have any missing values
            saint_baits.set(pd.read_csv(input.bait.get()[0]['datapath'], sep="\t", header=None, index_col=None, names=["Experiment Name", "Bait", "Type"]))
            # Add a new column for Bait ID
            saint_baits.get()['Bait ID'] = 'None'  # Default value for Bait ID

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
                    ui.notification_show(
                        f"ED file is missing required columns: {', '.join(missing_cols)}",
                        type="error"
                    )
                    return

                # Add Bait ID column if not present
                if "Bait ID" not in ed_df.columns:
                    ed_df['Bait ID'] = 'None'

                # Validate Type values
                invalid_types = ed_df[~ed_df['Type'].isin(['C', 'T'])]
                if len(invalid_types) > 0:
                    ui.notification_show(
                        f"ED file contains invalid Type values. Must be 'C' or 'T'.",
                        type="error"
                    )
                    return

                ed_dataframe.set(ed_df)
            except Exception as e:
                ui.notification_show(
                    f"Error reading ED file: {str(e)}",
                    type="error"
                )

    @reactive.effect
    @reactive.event(input.input_format)
    def clear_tables_on_format_change():
        # Clear tables when format changes to avoid showing stale data
        if input.input_format.get() in ["MaxQuant", "DIA-NN"]:
            # Clear SAINT bait table
            saint_baits.set(pd.DataFrame(columns=["Experiment Name", "Bait", "Type", "Bait ID"]))
        elif input.input_format.get() == "SAINT":
            # Clear ED table
            ed_dataframe.set(pd.DataFrame(columns=["Experiment Name", "Type", "Bait", "Replicate", "Bait ID"]))

    @reactive.effect
    @reactive.event(input.parse_data)
    def parse_data():
        # Check if the dataset name is valid
        dataset_name = input.dataset_name.get()
        check_result = parse.validate_name(dataset_name, datasets.get()['Dataset Name'].tolist())
        if check_result != 0:
            ui.notification_show(
                    f"Parser: {check_result}",
                    type="error",
                )
            return "Error: " + check_result

        # Get the input format (MaxQuant, DIA-NN, or SAINT)
        input_format = input.input_format.get()
        output_path = out_dir + '/' + dataset_name

        try:
            with ui.Progress(min=0, max=1) as progress:
                if input_format == "MaxQuant":
                    progress.set(message="Parsing MaxQuant inputs", value=0.25)

                    # Check if files are uploaded
                    if not input.pg_file.get() or not input.ed_file.get():
                        ui.notification_show(
                            "Please upload both proteinGroups.txt and Experimental Design files",
                            type="error"
                        )
                        return "Error: Missing files"

                    progress.set(0.45)

                    n_exp, n_ctrl = parse.parse_ed_pg(
                        input.pg_file.get()[0]['datapath'],
                        input.ed_file.get()[0]['datapath'],
                        input.quant_type.get(),
                        output_path
                    )

                    # Overwrite ED.csv with edited table data
                    ed_df = ed_table_mq.data_view()
                    ed_df.to_csv(f"{output_path}/ED.csv", index=False)

                    progress.set(0.85)

                    # Update the datasets dataframe
                    new_row = pd.DataFrame(
                        [[dataset_name, 'MaxQuant', input.quant_type.get(), n_exp, n_ctrl, '', '', '']],
                        columns=datasets.get().columns
                    )
                    updated_datasets = pd.concat([datasets.get(), new_row], ignore_index=True)
                    datasets.set(updated_datasets)
                    progress.set(1.0)

                elif input_format == "DIA-NN":
                    progress.set(message="Parsing DIA-NN inputs", value=0.25)

                    # Check if files are uploaded
                    if not input.diann_matrix_file.get() or not input.ed_file.get():
                        ui.notification_show(
                            "Please upload both DIA-NN matrix and Experimental Design files",
                            type="error"
                        )
                        return "Error: Missing files"

                    progress.set(0.45)

                    n_exp, n_ctrl = parse.parse_diann(
                        input.diann_matrix_file.get()[0]['datapath'],
                        input.ed_file.get()[0]['datapath'],
                        "Intensity",  # DIA-NN always uses intensity
                        output_path
                    )

                    # Overwrite ED.csv with edited table data
                    ed_df = ed_table_diann.data_view()
                    ed_df.to_csv(f"{output_path}/ED.csv", index=False)

                    progress.set(0.85)

                    # Update the datasets dataframe
                    new_row = pd.DataFrame(
                        [[dataset_name, 'DIA-NN', 'Intensity', n_exp, n_ctrl, '', '', '']],
                        columns=datasets.get().columns
                    )
                    updated_datasets = pd.concat([datasets.get(), new_row], ignore_index=True)
                    datasets.set(updated_datasets)
                    progress.set(1.0)

                elif input_format == "SAINT":
                    progress.set(message="Parsing SAINT inputs", value=0.25)

                    # Check if files are uploaded
                    if not input.bait.get() or not input.prey.get() or not input.interaction.get():
                        ui.notification_show(
                            "Please upload all three SAINT files (bait, prey, interaction)",
                            type="error"
                        )
                        return "Error: Missing files"

                    progress.set(0.45)

                    # Create the output directory if it doesn't exist
                    if not os.path.exists(output_path):
                        os.makedirs(output_path)

                    n_expts, n_ctrls = parse.parse_from_saint(
                        bait_table.data_view(),
                        input.prey.get()[0]['datapath'],
                        input.interaction.get()[0]['datapath'],
                        output_path
                    )

                    # Copy the bait, prey, and interaction files to the output directory
                    shutil.copy(input.bait.get()[0]['datapath'], output_path + '/bait.txt')
                    shutil.copy(input.prey.get()[0]['datapath'], output_path + '/prey.txt')
                    shutil.copy(input.interaction.get()[0]['datapath'], output_path + '/interaction.txt')

                    progress.set(0.85)

                    # Update the datasets dataframe
                    new_row = pd.DataFrame([[dataset_name, 'SAINT', input.quant_type.get(), n_expts, n_ctrls, '', '', '']], columns=datasets.get().columns)

                    updated_datasets = pd.concat([datasets.get(), new_row], ignore_index=True)
                    datasets.set(updated_datasets)
                    progress.set(1.0)

            # Success notification
            ui.notification_show(
                f"Successfully parsed dataset '{dataset_name}'",
                type="message",
                duration=5
            )

        except ProxiMateError as e:
            # Handle our custom exceptions with user-friendly messages
            error_msg = format_error_notification(e)
            ui.notification_show(
                error_msg,
                type="error",
                duration=None  # Keep error visible until dismissed
            )
            return f"Error: {e.user_message}"

        except FileNotFoundError as e:
            ui.notification_show(
                f"File not found: {str(e)}",
                type="error",
                duration=10
            )
            return "Error: File not found"

        except PermissionError as e:
            ui.notification_show(
                f"Permission denied accessing file: {str(e)}",
                type="error",
                duration=10
            )
            return "Error: Permission denied"

        except pd.errors.ParserError as e:
            ui.notification_show(
                f"Error parsing file: {str(e)}\n\nEnsure files are in correct format.",
                type="error",
                duration=10
            )
            return "Error: File parsing failed"

        except Exception as e:
            # Catch-all for unexpected errors
            import traceback
            traceback.print_exc()
            ui.notification_show(
                f"An unexpected error occurred:\n{str(e)}\n\nPlease check the console for details.",
                type="error",
                duration=None
            )
            return f"Error: {str(e)}"

        return "Parsed!"



    def do_clear_datasets():
        datasets.set(pd.DataFrame(columns=['Dataset Name', 'Input Type', 'Quant Type', 'Experiments', 'Controls', 'Scored', 'Imputation', 'WDFDR iterations']))
        # Clear the output folder
        for root, dirs, files in os.walk(out_dir):
            for file in files:
                abs_file = os.path.join(root, file)
                os.remove(abs_file)
            for dir in dirs:
                abs_dir = os.path.join(root, dir)
                shutil.rmtree(abs_dir)

    @reactive.effect
    @reactive.event(input.clear_datasets)
    def clear_datasets():
        do_clear_datasets()

    @render.download()
    def download_session():
        # Save the current state of the datasets dataframe to a CSV file
        datasets.get().to_csv(out_dir + "/datasets.csv", index=False)

        file_prefix = f"ProxiMateSession_{datetime.datetime.now().strftime('%Y%m%d')}"
        tmp_zip = tempfile.NamedTemporaryFile(prefix=file_prefix, suffix=".zip", delete=False)
        with zipfile.ZipFile(tmp_zip, "w", zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(out_dir):
                for file in files:
                    abs_file = os.path.join(root, file)
                    # Write the file using a relative path
                    zipf.write(abs_file, arcname=os.path.relpath(abs_file, out_dir))
        tmp_zip.close()
        # Return the path of the zip file for download
        return tmp_zip.name
    
    @reactive.effect
    @reactive.event(input.upload_session)
    def upload_session():
        # Start by clearing the current datasets
        do_clear_datasets()

        # Check if the upload has occurred
        uploaded = input.session_file.get()
        if uploaded:
            # uploaded is a list of dicts; use the first file
            zip_path = uploaded[0]['datapath']
            # Choose a destination directory (for example, your session directory)
            session_dest = "/Outputs"
            # Extract the zip file to the destination
            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                zip_ref.extractall(session_dest)

            # Load the datasets.csv file into the datasets reactive value
            datasets.set(pd.read_csv(os.path.join(session_dest, "datasets.csv")))


    # Scoring data
    @reactive.effect
    @reactive.event(input.score_data)
    def score_data():
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Scoring data", detail="Gathering inputs...", value=0)

            # Get the quant type from the datasets dataframe
            quant_type = datasets.get().loc[datasets.get()['Dataset Name'] == input.score_dataset.get(), 'Quant Type'].values[0]

            progress.set(message="Scoring data", detail="Running SAINT...", value=25)
            # Code for scoring goes here
            result = subprocess.run([
                "python3",
                "/Scripts/score.py",
                "--experimentalDesign",
                out_dir + '/' + input.score_dataset.get() + "/ED.csv",
                "--scoreInputs",
                out_dir + '/' + input.score_dataset.get(),
                "--outputPath",
                out_dir + '/' + input.score_dataset.get(),
                "--n-iterations",
                str(input.wdfdr_iterations.get()),
                "--imputation",
                input.imputation_method.get(),
                "--quantType",
                quant_type,
            ], capture_output=True, text=True)

            progress.set(message="Scoring data", detail="Adding protein annotation...", value=65)

            # Add script for annotation here
            ann_result = subprocess.run([
                "python3",
                "/Scripts/annotator.py",
                "--organism",
                input.organism.get(),
                "--scoreFile",
                out_dir + '/' + input.score_dataset.get() + "/merged.csv",
                "--outputDir",
                out_dir + '/' + input.score_dataset.get(),
            ], capture_output=True, text=True)

            progress.set(message="Scoring data", detail="Updating datasets...", value=90)

            # Update the datasets dataframe
            curr_dataset = datasets.get().copy()
            # Edit the row for the current dataset

            imp_mapping = {0: 'Default', 1: 'Prey-specific'}

            curr_dataset.loc[curr_dataset['Dataset Name'] == input.score_dataset.get(), 'Scored'] = 'Yes'
            curr_dataset.loc[curr_dataset['Dataset Name'] == input.score_dataset.get(), 'Imputation'] = imp_mapping[int(input.imputation_method.get())]
            curr_dataset.loc[curr_dataset['Dataset Name'] == input.score_dataset.get(), 'WDFDR iterations'] = input.wdfdr_iterations.get()
            datasets.set(curr_dataset)
            progress.set(message="Scoring data", detail="Done!", value=100)

    # Quality controls tab
    @render_plotly
    def raw_pca_plot():
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Generating Plots...", value=25)
            # Get the selected dataset
            dataset_name = input.qc_dataset.get()
            if not dataset_name:
                return None
            
            # Build the file paths
            interaction_path = os.path.join(out_dir, dataset_name, "interaction.txt")
            ed_path = os.path.join(out_dir, dataset_name, "ED.csv")

            # Check if the files exist
            if not (os.path.exists(interaction_path) and os.path.exists(ed_path)):
                return None
            
            # Generate the PCA plot
            fig = pca_plot(interaction_path, ed_path)

            return fig
    
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
            return ui.value_box("Median Network Size", "No data", showcase=None)

        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        if not os.path.exists(results_path):
            return ui.value_box("Median Network Size", "No data", showcase=None)

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
            return ui.value_box("Known Enrichment", "No data", showcase=None)

        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        if not os.path.exists(results_path):
            return ui.value_box("Known Enrichment", "No data", showcase=None)

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
            return ui.value_box("Mean Prey-Prey Degree", "No data", showcase=None)

        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

        if not os.path.exists(results_path):
            return ui.value_box("Mean Prey-Prey Degree", "No data", showcase=None)

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
            pass
        return None

    # Threshold presets for Data Thresholding tab
    # Preset values: Stringent (0.9, 0.01, 2.0, 0.01), Moderate (0.7, 0.05, 1.0, 0.05), Relaxed (0.5, 0.1, 0.0, 0.1)
    @reactive.effect
    @reactive.event(input.qc_preset_stringent)
    def apply_qc_stringent_preset():
        ui.update_slider("threshold_saintscore", value=0.9)
        ui.update_slider("threshold_bfdr", value=0.01)
        ui.update_slider("threshold_wd", value=2.0)
        ui.update_slider("threshold_wdfdr", value=0.01)

    @reactive.effect
    @reactive.event(input.qc_preset_moderate)
    def apply_qc_moderate_preset():
        ui.update_slider("threshold_saintscore", value=0.7)
        ui.update_slider("threshold_bfdr", value=0.05)
        ui.update_slider("threshold_wd", value=1.0)
        ui.update_slider("threshold_wdfdr", value=0.05)

    @reactive.effect
    @reactive.event(input.qc_preset_relaxed)
    def apply_qc_relaxed_preset():
        ui.update_slider("threshold_saintscore", value=0.5)
        ui.update_slider("threshold_bfdr", value=0.1)
        ui.update_slider("threshold_wd", value=0.0)
        ui.update_slider("threshold_wdfdr", value=0.1)

    # Threshold presets for Downloads tab
    @reactive.effect
    @reactive.event(input.dl_preset_stringent)
    def apply_dl_stringent_preset():
        ui.update_slider("dl_threshold_saintscore", value=0.9)
        ui.update_slider("dl_threshold_bfdr", value=0.01)
        ui.update_slider("dl_threshold_wd", value=2.0)
        ui.update_slider("dl_threshold_wdfdr", value=0.01)

    @reactive.effect
    @reactive.event(input.dl_preset_moderate)
    def apply_dl_moderate_preset():
        ui.update_slider("dl_threshold_saintscore", value=0.7)
        ui.update_slider("dl_threshold_bfdr", value=0.05)
        ui.update_slider("dl_threshold_wd", value=1.0)
        ui.update_slider("dl_threshold_wdfdr", value=0.05)

    @reactive.effect
    @reactive.event(input.dl_preset_relaxed)
    def apply_dl_relaxed_preset():
        ui.update_slider("dl_threshold_saintscore", value=0.5)
        ui.update_slider("dl_threshold_bfdr", value=0.1)
        ui.update_slider("dl_threshold_wd", value=0.0)
        ui.update_slider("dl_threshold_wdfdr", value=0.1)

    @reactive.effect
    @reactive.event(input.dl_preset_none)
    def apply_dl_no_filter_preset():
        ui.update_slider("dl_threshold_saintscore", value=0.0)
        ui.update_slider("dl_threshold_bfdr", value=1.0)
        ui.update_slider("dl_threshold_wd", value=0.0)
        ui.update_slider("dl_threshold_wdfdr", value=1.0)

    # Plot export download handlers (using matplotlib for PNG export)
    @render.download(filename="pca_plot.png")
    def download_pca_plot():
        dataset_name = input.qc_dataset.get()
        if not dataset_name:
            ui.notification_show("No dataset selected.", type="error")
            return None
        interaction_path = os.path.join(out_dir, dataset_name, "interaction.txt")
        ed_path = os.path.join(out_dir, dataset_name, "ED.csv")
        if not (os.path.exists(interaction_path) and os.path.exists(ed_path)):
            ui.notification_show("Required files not found.", type="error")
            return None
        fig = pca_plot_matplotlib(interaction_path, ed_path)
        filepath = os.path.join(out_dir, "pca_plot.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.download(filename="scatter_plot.png")
    def download_scatter_plot():
        dataset_name = input.qc_dataset.get()
        bait_selection = input.qc_bait.get()
        if not dataset_name or bait_selection == "All":
            ui.notification_show("Select a specific bait to export the scatter plot.", type="error")
            return None
        results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")
        saintscore_threshold = input.threshold_saintscore.get()
        fig = saint_scatter_matplotlib(results_path, bait_selection, saintscore_threshold)
        filepath = os.path.join(out_dir, "scatter_plot.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.download(filename="heatmap.png")
    def download_heatmap():
        dataset = input.feature_dataset.get()
        if not dataset:
            ui.notification_show("No dataset selected.", type="error")
            return None
        feature_file = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(feature_file):
            ui.notification_show("Run feature analysis first.", type="error")
            return None
        feature_data = pd.read_csv(feature_file)
        feature_type = input.feature_type.get() or 'GO_CC'
        num_features = input.num_features.get() or 30
        try:
            fig = plot_results(feature_data, feature_type, num_features)
            filepath = os.path.join(out_dir, "heatmap.png")
            fig.savefig(filepath, dpi=150, bbox_inches='tight')
            return filepath
        except ValueError:
            ui.notification_show("Insufficient data to generate heatmap.", type="error")
            return None

    @render.download(filename="volcano_plot.png")
    def download_volcano_plot():
        dataset_a = input.comp_dataset_a.get()
        dataset_b = input.comp_dataset_b.get()
        bait_a = input.comp_bait_a.get()
        bait_b = input.comp_bait_b.get()
        if not all([dataset_a, dataset_b, bait_a, bait_b]):
            ui.notification_show("Select datasets and baits first.", type="error")
            return None
        if dataset_a != dataset_b:
            ui.notification_show("Volcano plot requires baits from the same dataset.", type="error")
            return None
        # Use cached volcano data
        volcano_data = comp_volcano_data_cached()
        if volcano_data.empty:
            ui.notification_show("No data available for volcano plot.", type="error")
            return None
        # Use matplotlib version for export
        fig = create_volcano_plot_matplotlib(volcano_data, bait_a, bait_b)
        filepath = os.path.join(out_dir, "volcano_plot.png")
        fig.savefig(filepath, dpi=150, bbox_inches='tight', facecolor='white')
        return filepath

    @render.download(filename="venn_diagram.png")
    def download_venn_diagram():
        bait_a = input.comp_bait_a.get()
        bait_b = input.comp_bait_b.get()
        if not all([bait_a, bait_b]):
            ui.notification_show("Select baits first.", type="error")
            return None
        # Use cached filtered data to create sets
        data_a = comp_filtered_data_a_cached()
        data_b = comp_filtered_data_b_cached()
        set_a = set(data_a['Prey.ID'].unique()) if len(data_a) > 0 else set()
        set_b = set(data_b['Prey.ID'].unique()) if len(data_b) > 0 else set()
        if len(set_a) == 0 and len(set_b) == 0:
            ui.notification_show("No data available for Venn diagram.", type="error")
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

    _cytoscape_status_msg = reactive.Value("Click a button to test Cytoscape connection...")

    # Feature analysis tab
    @reactive.effect
    @reactive.event(input.feature_analysis)
    def feature_analysis():
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Running protein feature analysis", value=5)

            # Get the selected dataset
            progress.set(message="Running protein feature analysis", detail="Loading dataset...", value=20)
            dataset = pd.read_csv(os.path.join(out_dir, input.feature_dataset.get(), "annotated_scores.csv"))

            progress.set(message="Running protein feature analysis", detail="Processing data...", value=50)
            result = process_refactored(
                dataset,
                columns_for_analysis = ['GO_CC', 'GO_BP', 'GO_MF', 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains'],
                threshold = input.saint_threshold.get()
            )

            progress.set(message="Running protein feature analysis", detail="Generating plots...", value=80)
            # Store the results in the reactive value
            feature_enrichment.set(result)

            progress.set(message="Running protein feature analysis", detail="Saving results...", value=90)
            # Save the results to the dataset directory
            result.to_csv(os.path.join(out_dir, input.feature_dataset.get(), "Feature_enrichment.csv"), index=False)
            
            progress.set(message="Running protein feature analysis", detail="Done!", value=100)

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
        except Exception as e:
            ui.update_select("download_bait_filter", choices=["All"])

    @render.download()
    def download_enrichment():
        """Download filtered enrichment results."""
        dataset = input.feature_dataset.get()
        if not dataset:
            ui.notification_show("No dataset selected.", type="error")
            return None

        feature_file = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(feature_file):
            ui.notification_show("No enrichment results available. Please run feature analysis first.", type="error")
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
            ui.notification_show("No results match the current filters. Try adjusting the thresholds.", type="warning")
            return None

        # Save to temp file and return
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"enrichment_{dataset}_{timestamp}.csv"
        savepath = os.path.join(out_dir, filename)
        filtered_data.to_csv(savepath, index=False)
        return savepath

    @reactive.Calc
    def available_datasets():
        return datasets.get()["Dataset Name"].tolist()

    @reactive.Calc
    def scored_datasets():
        return datasets.get()[datasets.get()['Scored'] == 'Yes']["Dataset Name"].tolist()

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
                ui.update_select("qc_bait", choices=["All"])  # Reset to default if file not found
            except Exception as e:
                ui.update_select("qc_bait", choices=["All"])
        else:
            ui.update_select("qc_bait", choices=["All"])  # Reset to default if no dataset is selected

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
                ui.update_select("comp_bait_a", choices=[])
            except Exception as e:
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
                ui.update_select("comp_bait_b", choices=[])
            except Exception as e:
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

    @render.text
    @reactive.event(input.compare_networks)
    def genes_a_only():
        """Render list of genes only in network A."""
        # Get selections
        dataset_a = input.comp_dataset_a()
        dataset_b = input.comp_dataset_b()
        bait_a = input.comp_bait_a()
        bait_b = input.comp_bait_b()

        if not all([dataset_a, dataset_b, bait_a, bait_b]):
            return "Select datasets and baits to view gene lists"

        # Get cached filtered data
        data_a = comp_filtered_data_a_cached()
        data_b = comp_filtered_data_b_cached()

        # Get sets
        set_a = set(data_a['Prey.ID'].unique()) if len(data_a) > 0 else set()
        set_b = set(data_b['Prey.ID'].unique()) if len(data_b) > 0 else set()
        only_a = set_a - set_b

        if len(only_a) == 0:
            return "No unique genes in Network A"

        # Get gene names
        genes_only_a = data_a[data_a['Prey.ID'].isin(only_a)]['First_Prey_Gene'].unique()
        genes_sorted = sorted(genes_only_a)

        return f"Network A only ({len(genes_sorted)} genes):\n\n" + "\n".join(genes_sorted)

    @render.text
    @reactive.event(input.compare_networks)
    def genes_b_only():
        """Render list of genes only in network B."""
        # Get selections
        dataset_a = input.comp_dataset_a()
        dataset_b = input.comp_dataset_b()
        bait_a = input.comp_bait_a()
        bait_b = input.comp_bait_b()

        if not all([dataset_a, dataset_b, bait_a, bait_b]):
            return "Select datasets and baits to view gene lists"

        # Get cached filtered data
        data_a = comp_filtered_data_a_cached()
        data_b = comp_filtered_data_b_cached()

        # Get sets
        set_a = set(data_a['Prey.ID'].unique()) if len(data_a) > 0 else set()
        set_b = set(data_b['Prey.ID'].unique()) if len(data_b) > 0 else set()
        only_b = set_b - set_a

        if len(only_b) == 0:
            return "No unique genes in Network B"

        # Get gene names
        genes_only_b = data_b[data_b['Prey.ID'].isin(only_b)]['First_Prey_Gene'].unique()
        genes_sorted = sorted(genes_only_b)

        return f"Network B only ({len(genes_sorted)} genes):\n\n" + "\n".join(genes_sorted)

    @render.text
    @reactive.event(input.compare_networks)
    def genes_both():
        """Render list of genes in both networks."""
        # Get selections
        dataset_a = input.comp_dataset_a()
        dataset_b = input.comp_dataset_b()
        bait_a = input.comp_bait_a()
        bait_b = input.comp_bait_b()

        if not all([dataset_a, dataset_b, bait_a, bait_b]):
            return "Select datasets and baits to view gene lists"

        # Get cached filtered data
        data_a = comp_filtered_data_a_cached()
        data_b = comp_filtered_data_b_cached()

        # Get sets
        set_a = set(data_a['Prey.ID'].unique()) if len(data_a) > 0 else set()
        set_b = set(data_b['Prey.ID'].unique()) if len(data_b) > 0 else set()
        both = set_a & set_b

        if len(both) == 0:
            return "No shared genes between networks"

        # Get gene names (use data_a arbitrarily since they're in both)
        genes_both = data_a[data_a['Prey.ID'].isin(both)]['First_Prey_Gene'].unique()
        genes_sorted = sorted(genes_both)

        return f"Both networks ({len(genes_sorted)} genes):\n\n" + "\n".join(genes_sorted)

    @reactive.Effect
    def update_selectize_custom_columns():
        # Update the custom columns selectize input based on the available datasets
        dataset_name = input.download_dataset.get()
        if dataset_name:
            # Read the dataset to get the columns
            dataset_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")
            if os.path.exists(dataset_path):
                df = pd.read_csv(dataset_path)
                available_columns = df.columns.tolist()
                # Sort the columns alphabetically
                available_columns.sort()
                ui.update_selectize("custom_columns", choices=available_columns, selected=["Experiment.ID", "Prey.ID", "SaintScore", "BFDR"])

    custom_dataset = reactive.Value(pd.DataFrame())
    custom_dataset_total = reactive.Value(0)

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

    @render.data_frame
    def custom_table():
        # Use cached data instead of re-reading file on every slider change
        df = cached_download_data()
        if df.empty:
            custom_dataset_total.set(0)
            return pd.DataFrame()

        # Custom cols will come from a user input
        custom_cols = list(input.custom_columns.get())
        if not custom_cols:
            # Default columns if none are selected
            custom_cols = ["Experiment.ID", "Prey.ID", "SaintScore", "BFDR"]

        # Store total row count before filtering
        custom_dataset_total.set(len(df))

        # Apply threshold filtering using centralized function
        thresholds = {
            'SaintScore': input.dl_threshold_saintscore(),
            'BFDR': input.dl_threshold_bfdr(),
            'WD': input.dl_threshold_wd(),
            'WDFDR': input.dl_threshold_wdfdr()
        }
        df = apply_score_thresholds(df, thresholds)

        custom_df = df[custom_cols]
        custom_dataset.set(custom_df)

        return custom_dataset.get()

    @render.text
    def download_row_count():
        """Display the filtered row count."""
        df = custom_dataset.get()
        total = custom_dataset_total.get()

        if df.empty and total == 0:
            return ""

        return f"Showing {len(df)} of {total} interactions"
    
    # Add a download button for the custom dataset
    @render.download()
    def download_custom_dataset():
        if custom_dataset.get().empty:
            ui.notification_show(
                "No custom dataset to download. Please select columns first.",
                type="error",
            )
            return None
        ui.notification_show("Preparing download...", type="message", duration=2)
        # Create a temporary file to save the custom dataset
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"custom_dataset_{timestamp}.csv"
        savepath = os.path.join(out_dir, filename)
        custom_dataset.get().to_csv(savepath, index=False)
        return savepath

    @render.download()
    def download_batch():
        """Download all results for a dataset as a ZIP file."""
        dataset = input.download_dataset.get()
        if not dataset:
            ui.notification_show("No dataset selected.", type="error")
            return None

        ui.notification_show("Preparing ZIP file...", type="message", duration=2)
        dataset_dir = os.path.join(out_dir, dataset)
        if not os.path.exists(dataset_dir):
            ui.notification_show("Dataset directory not found.", type="error")
            return None

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        zip_filename = f"{dataset}_all_{timestamp}.zip"
        zip_path = os.path.join(out_dir, zip_filename)

        # List of files to include in the ZIP
        files_to_include = [
            'merged.csv',
            'annotated_scores.csv',
            'Feature_enrichment.csv',
            'bait.txt',
            'prey.txt',
            'interaction.txt',
            'ED.csv'
        ]

        included_files = []
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
            for filename in files_to_include:
                filepath = os.path.join(dataset_dir, filename)
                if os.path.exists(filepath):
                    zf.write(filepath, filename)
                    included_files.append(filename)

        if not included_files:
            ui.notification_show("No files found to include in ZIP.", type="error")
            os.remove(zip_path)
            return None

        ui.notification_show(f"ZIP created with {len(included_files)} files.", type="message", duration=3)
        return zip_path

    # Cytoscape tab
    @reactive.effect
    @reactive.event(input.test_create_node)
    def create_test_node():
        """Create a test node in Cytoscape."""
        try:
            base_url = get_cytoscape_base_url()

            # Try to connect to Cytoscape
            version = p4c.cytoscape_version_info(base_url=base_url)

            # Create a simple network if none exists
            try:
                current_network = p4c.get_network_name(base_url=base_url)
            except:
                # No network exists, create one
                nodes_df = pd.DataFrame({'id': ['InitialNode']})
                edges_df = pd.DataFrame({'source': [], 'target': []})
                p4c.create_network_from_data_frames(
                    nodes=nodes_df,
                    edges=edges_df,
                    title="ProxiMate Test Network",
                    base_url=base_url
                )

            # Add a test node
            p4c.add_cy_nodes(['TestNode_ProxiMate'], base_url=base_url)

            # Update status
            status_message = f"✓ Successfully created 'TestNode_ProxiMate' in Cytoscape\nCytoscape version: {version['cytoscapeVersion']}"

        except Exception as e:
            status_message = f"✗ Error connecting to Cytoscape:\n{str(e)}\n\nMake sure Cytoscape is running on your host machine."

        # Store status in reactive value for display
        _cytoscape_status_msg.set(status_message)

    @reactive.effect
    @reactive.event(input.test_delete_node)
    def delete_test_node():
        """Delete the test node from Cytoscape."""
        try:
            base_url = get_cytoscape_base_url()

            # Select the test node
            p4c.select_nodes(['TestNode_ProxiMate'], by_col='name', base_url=base_url)

            # Delete selected nodes
            p4c.delete_selected_nodes(base_url=base_url)

            status_message = "✓ Successfully deleted 'TestNode_ProxiMate' from Cytoscape"

        except Exception as e:
            status_message = f"✗ Error deleting node:\n{str(e)}\n\nMake sure the node exists and Cytoscape is running."

        _cytoscape_status_msg.set(status_message)

    @render.text
    def cytoscape_status():
        return _cytoscape_status_msg.get()

app = App(app_ui, server)

if __name__ == "__main__":
    run_app(app, host="0.0.0.0", port=3838)
