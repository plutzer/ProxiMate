from functools import partial
from shiny import App, Inputs, Outputs, Session, reactive, render, ui, run_app
import plotly.graph_objects as go
import plotly.express as px
# from shiny import render_plotly
from shinywidgets import output_widget, render_widget, render_plotly
import pandas as pd
import os
import parse
import subprocess
import zipfile
import tempfile
import datetime
import shutil
from QC_plots import pca_plot, saint_known_retention, roc_plot
from Ann_Enrichment import process_refactored, plot_results


out_dir = "/Outputs"

app_ui = ui.page_navbar(
    ui.nav_spacer(),
    ui.nav_panel("Input Data",
                 "Use either method to parse input data:",
                 ui.layout_columns(
                    ui.card(
                        ui.card_header("Proteomics Data Input"),
                        ui.input_text("mq_dataset_name",
                                    "Dataset Name",
                                    placeholder="No spaces or special characters (/ \\ : * ? \" < > |)"),
                        ui.input_select("input_format", "Input Format",
                                       choices=["MaxQuant", "DIA-NN"],
                                       selected="MaxQuant"),
                        ui.panel_conditional(
                            "input.input_format === 'MaxQuant'",
                            ui.input_file("pg_file", "MaxQuant proteinGroups.txt file"),
                            ui.input_select("mq_quant_type", "Quantification Type",
                                          choices=["Intensity", "LFQ", "Spectral Counts"],
                                          selected="Intensity")
                        ),
                        ui.panel_conditional(
                            "input.input_format === 'DIA-NN'",
                            ui.input_file("diann_matrix_file", "DIA-NN report.pg_matrix.tsv file"),
                            ui.panel_well(
                                "Note: DIA-NN uses intensity-based quantification."
                            )
                        ),
                        ui.input_file("ed_file", "Experimental Design File"),
                        ui.panel_well(
                            "Experimental Design File Format:",
                            "Experiment Name, Type, Bait, Replicate, Bait ID"
                        ),
                        ui.input_action_button("parse_mq", "Parse Data")

                    ),
                    ui.card(
                        ui.card_header('SAINT Input'),
                        ui.input_text("st_dataset_name", "Dataset Name"),
                        ui.layout_columns(
                            ui.panel_well(
                                ui.input_file("bait", "SAINT bait.txt file"),
                                ui.input_file("prey", "SAINT prey.txt file"),
                                ui.input_file("interaction", "SAINT interaction.txt file"),
                                # Add a link to the SAINT documentation for input formats
                                # ui.link("Documentation for SAINT input format", "www.google.com"),# Placeholder
                                ui.input_select("saint_quant_type", "Quantification Type", choices=["Intensity", "LFQ", "Spectral Counts"]),
                                ui.input_action_button("parse_saint", "Parse SAINT Inputs")
                            ),
                            ui.panel_well(
                                "CompPASS Scoring requires a mapping of Experiment IDs to their corresponding uniprot IDs. Edit the values in the Bait ID column for accurate CompPASS scoring.",
                                ui.output_data_frame("bait_table"),
                            ),
                        ),
                    ),
                    col_widths=[4,8]
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
    ui.nav_panel("Score Data",
                    ui.input_select("score_dataset", "Select Dataset", choices=[]), # Need this to be dynamic
                    ui.input_radio_buttons("imputation_method", "Imputation Method", choices={0: "Default", 1: "Prey-specific"}),
                    ui.input_numeric("wdfdr_iterations", "WDFDR Iterations", value=100),
                    ui.input_action_button("score_data", "Score Data"),
                    # ui.output_text_verbatim("log_text"),
    ),
    ui.nav_panel("Quality Controls",
                    ui.layout_columns(
                        ui.input_select("qc_dataset", "Select Dataset", choices=[]), # Need this to be dynamic
                        ui.input_select("qc_bait", "Select Control Bait", choices=["All"]), # Need this to be dynamic
                    ),
                    ui.layout_columns(
                        ui.card(
                            ui.card_header("Raw Data UMAP"),
                            output_widget("raw_pca_plot"),
                        ),
                        ui.card(
                            # Add dropdown to select bait
                            ui.card_header("Known Physical Interactions"),
                            output_widget("known_retention_plot")
                        ),
                        ui.card(
                            ui.card_header("ROC Curve for BioGRID Interactions"),
                            ui.input_radio_buttons("roc_known_type", "True Positive Type", choices={"BioGRID":"All Physical Interactions", "Multivalidated":"Multivalidated Physical Interactions"}),
                            output_widget("roc_curve_plot")
                        ),
                        col_widths=(4,4,4),
                    ),
    ),
    ui.nav_panel("Protein Feature Analysis",
                ui.layout_columns(
                    ui.card(
                        ui.card_header("Parameters for Feature Analysis"),
                        ui.input_select("feature_dataset", "Select Dataset", choices=[]), # Need this to be dynamic
                        ui.input_slider("saint_threshold", "Saint Threshold", min=0.0, max=1.0, value=0.9),
                        ui.input_action_button("feature_analysis", "Run Feature Analysis"),
                    ),
                    ui.card(
                        ui.card_header("Feature Enrichment Analysis"),
                        ui.input_select("feature_type", "Select Feature Type", choices=["GO_CC", "Motifs", "Regions", "Repeats", "Compositions", "Domains"]),
                        ui.input_numeric("num_features", "Number of Features to Display", value=30, min=1, max=100),
                        # Spot for a seaborn heatmap
                        ui.output_plot("feature_enrichment_plot"),
                    ),
                ),
    ),
    ui.nav_panel("Cytoscape",
                    "Page E content"
    ),
    ui.nav_panel("Downloads",
                    ui.input_select("download_dataset", "Select Dataset", choices=[]),
                    ui.layout_columns(
                        ui.card(
                            ui.card_header("Create a Custom Dataset"),
                            "Select Columns for the Custom Dataset",
                            ui.input_selectize("custom_columns", "Select Columns", choices=["Experiment.ID", "Prey.ID", "SaintScore", "BFDR"], multiple=True, selected=["Experiment.ID", "Prey.ID", "SaintScore", "BFDR"]),
                        ),
                        ui.card(
                            ui.card_header("Custom Dataset"),
                            # Render the custom dataset table here
                            ui.output_data_frame("custom_table"), # This dataset should not be editable
                            ui.download_button("download_custom_dataset", "Download Custom Dataset"),
                        ),
                        col_widths=(4,8)
                    ),
                    ui.card(
                        ui.card_header("Download Internal Datasets"),
                        "Premade datasets used internally by ProxiMate.",
                    )
    ),
    sidebar=ui.sidebar(
        "Sidebar content"
    ),
    title="ProxiMate",
)

def server(input: Inputs, output: Outputs, session: Session):
    print("Server started")

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

    @render.data_frame
    def bait_table():
        # Return the bait table for SAINT input
        return render.DataGrid(saint_baits.get(), editable=True)

    @reactive.effect
    @reactive.event(input.bait)
    def update_bait_table():
        # Check to see if the bait file has been uploaded
        if input.bait.get():
            # Read the bait file and set it to the saint_baits reactive value
            print("Bait file uploaded:", input.bait.get()[0]['name'])
            # Read the bait file into a DataFrame

            # Catch bad input files here

            # Check to make sure the columns are correct and don't have any missing values
            saint_baits.set(pd.read_csv(input.bait.get()[0]['datapath'], sep="\t", header=None, index_col=None, names=["Experiment Name", "Bait", "Type"]))
            # Add a new column for Bait ID
            saint_baits.get()['Bait ID'] = 'None'  # Default value for Bait ID
            print("Bait table updated with", len(saint_baits.get()), "rows.")

    @reactive.effect
    @reactive.event(input.parse_mq)
    def parse_mq():
        # Check if the dataset name is valid
        dataset_name = input.mq_dataset_name.get()
        check_result = parse.validate_name(dataset_name, datasets.get()['Dataset Name'].tolist())
        if check_result != 0:
            print(check_result)
            ui.notification_show(
                    f"Parser: {check_result}",
                    type="error",
                )
            return "Error: " + check_result

        # Get the input format (MaxQuant or DIA-NN)
        input_format = input.input_format.get()

        with ui.Progress(min=0, max=1) as progress:
            if input_format == "MaxQuant":
                progress.set(message="Parsing MaxQuant inputs", value=0.25)
                print("Parsing MaxQuant and Experimental Design files")

                # Check if files are uploaded
                if not input.pg_file.get() or not input.ed_file.get():
                    ui.notification_show(
                        "Please upload both proteinGroups.txt and Experimental Design files",
                        type="error"
                    )
                    return "Error: Missing files"

                print(input.pg_file.get()[0]['name'])
                print(input.ed_file.get()[0]['name'])
                curr_dataset = datasets.get()
                progress.set(0.45)
                print("Current directory:", os.getcwd())

                n_exp, n_ctrl = parse.parse_ed_pg(
                    input.pg_file.get()[0]['datapath'],
                    input.ed_file.get()[0]['datapath'],
                    input.mq_quant_type.get(),
                    out_dir + '/' + input.mq_dataset_name.get()
                )
                progress.set(0.85)

                # Update the datasets dataframe
                new_row = pd.DataFrame(
                    [[input.mq_dataset_name.get(), 'MaxQuant', input.mq_quant_type.get(), n_exp, n_ctrl, '', '', '']],
                    columns=datasets.get().columns
                )
                updated_datasets = pd.concat([datasets.get(), new_row], ignore_index=True)
                datasets.set(updated_datasets)
                progress.set(1.0)

            elif input_format == "DIA-NN":
                progress.set(message="Parsing DIA-NN inputs", value=0.25)
                print("Parsing DIA-NN and Experimental Design files")

                # Check if files are uploaded
                if not input.diann_matrix_file.get() or not input.ed_file.get():
                    ui.notification_show(
                        "Please upload both DIA-NN matrix and Experimental Design files",
                        type="error"
                    )
                    return "Error: Missing files"

                print(input.diann_matrix_file.get()[0]['name'])
                print(input.ed_file.get()[0]['name'])
                curr_dataset = datasets.get()
                progress.set(0.45)
                print("Current directory:", os.getcwd())

                n_exp, n_ctrl = parse.parse_diann(
                    input.diann_matrix_file.get()[0]['datapath'],
                    input.ed_file.get()[0]['datapath'],
                    "Intensity",  # DIA-NN always uses intensity
                    out_dir + '/' + input.mq_dataset_name.get()
                )
                progress.set(0.85)

                # Update the datasets dataframe
                new_row = pd.DataFrame(
                    [[input.mq_dataset_name.get(), 'DIA-NN', 'Intensity', n_exp, n_ctrl, '', '', '']],
                    columns=datasets.get().columns
                )
                updated_datasets = pd.concat([datasets.get(), new_row], ignore_index=True)
                datasets.set(updated_datasets)
                progress.set(1.0)

        return "Parsed!"

    @reactive.effect
    @reactive.event(input.parse_saint)
    def parse_saint():
        print("Parsing")
        print(input.bait.get()[0]['datapath'])
        print(input.prey.get()[0]['datapath'])
        print(input.interaction.get()[0]['datapath'])

        # Check if the dataset name is valid
        dataset_name = input.st_dataset_name.get()
        check_result = parse.validate_name(dataset_name, datasets.get()['Dataset Name'].tolist())
        if check_result != 0:
            print(check_result)
            ui.notification_show(
                    f"Parser: {check_result}",
                    type="error",
                )            
            return "Error: " + check_result
        
        # Create the output directory if it doesn't exist
        if not os.path.exists(out_dir + '/' + dataset_name):
            os.makedirs(out_dir + '/' + dataset_name)
        
        n_expts, n_ctrls = parse.parse_from_saint(
            bait_table.data_view(),
            input.prey.get()[0]['datapath'],
            input.interaction.get()[0]['datapath'],
            out_dir + '/' + input.st_dataset_name.get()
        )

        # Copy the bait, prey, and interaction files to the output directory
        shutil.copy(input.bait.get()[0]['datapath'], out_dir + '/' + dataset_name + '/bait.txt')
        shutil.copy(input.prey.get()[0]['datapath'], out_dir + '/' + dataset_name + '/prey.txt')
        shutil.copy(input.interaction.get()[0]['datapath'], out_dir + '/' + dataset_name + '/interaction.txt')

        # Update the datasets dataframe
        new_row = pd.DataFrame([[dataset_name, 'SAINT', input.saint_quant_type.get(), n_expts, n_ctrls, '', '', '']], columns=datasets.get().columns)

        updated_datasets = pd.concat([datasets.get(), new_row], ignore_index=True)
        datasets.set(updated_datasets)
        print("Parsed SAINT inputs. Updated datasets.")

        print("Edited bait table:", bait_table.data_view())
        # Code for parsing goes here



    def do_clear_datasets():
        datasets.set(pd.DataFrame(columns=['Dataset Name', 'Input Type', 'Quant Type', 'Experiments', 'Controls', 'Scored', 'Imputation', 'WDFDR iterations']))
        print("Datasets table cleared.")
        # Clear the output folder
        for root, dirs, files in os.walk(out_dir):
            for file in files:
                abs_file = os.path.join(root, file)
                os.remove(abs_file)
            for dir in dirs:
                abs_dir = os.path.join(root, dir)
                shutil.rmtree(abs_dir)
        print("Output folder cleared.")

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
            print("Session uploaded and extracted to:", session_dest)

            # Load the datasets.csv file into the datasets reactive value
            datasets.set(pd.read_csv(os.path.join(session_dest, "datasets.csv")))


    # Scoring data
    @reactive.effect
    @reactive.event(input.score_data)
    def score_data():
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Scoring data", detail="Gathering inputs...", value=0)
            print("Scoring data")
            print(input.score_dataset.get())
            print(input.imputation_method.get())
            print(input.wdfdr_iterations.get())

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
            print("Scoring result:", result.stdout)
            if result.stderr:
                print("Error:", result.stderr)
            print("Scoring completed.")

            progress.set(message="Scoring data", detail="Adding protein annotation...", value=65)

            # Add script for annotation here
            ann_result = subprocess.run([
                "python3",
                "/Scripts/annotator.py",
                "--scoreFile",
                out_dir + '/' + input.score_dataset.get() + "/merged.csv",
                "--outputDir",
                out_dir + '/' + input.score_dataset.get(),
            ], capture_output=True, text=True)
            print("Annotation result:", ann_result.stdout)
            if ann_result.stderr:
                print("Error:", ann_result.stderr)
            print("Annotation completed.")

            progress.set(message="Scoring data", detail="Updating datasets...", value=90)

            # Update the datasets dataframe
            curr_dataset = datasets.get().copy()
            # Edit the row for the current dataset

            imp_mapping = {0: 'Default', 1: 'Prey-specific'}

            curr_dataset.loc[curr_dataset['Dataset Name'] == input.score_dataset.get(), 'Scored'] = 'Yes'
            curr_dataset.loc[curr_dataset['Dataset Name'] == input.score_dataset.get(), 'Imputation'] = imp_mapping[int(input.imputation_method.get())]
            curr_dataset.loc[curr_dataset['Dataset Name'] == input.score_dataset.get(), 'WDFDR iterations'] = input.wdfdr_iterations.get()
            print(curr_dataset)
            datasets.set(curr_dataset)
            print("Dataset updated.")
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
            print(ctrls)
        
        if not os.path.exists(results_path):
            print(f"File not found for dataset {dataset_name}. Resetting bait selection.")
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
            print(ctrls)
        

        fig = roc_plot(results_path, known_type=input.roc_known_type.get(), ctrl_experiments=ctrls) # Add selected ctrls

        # Return the figure widget
        return fig

    feature_enrichment = reactive.Value(pd.DataFrame())

    # Feature analysis tab
    @reactive.effect
    @reactive.event(input.feature_analysis)
    def feature_analysis():
        with ui.Progress(min=0, max=100) as progress:
            progress.set(message="Running protein feature analysis", value=5)
            print("Running feature analysis")
            print(input.feature_dataset.get())
            print(input.saint_threshold.get())

            # Get the selected dataset
            progress.set(message="Running protein feature analysis", detail="Loading dataset...", value=20)
            dataset = pd.read_csv(os.path.join(out_dir, input.feature_dataset.get(), "annotated_scores.csv"))

            progress.set(message="Running protein feature analysis", detail="Processing data...", value=50)
            result = process_refactored(
                dataset,
                columns_for_analysis = ['GO_CC', 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains'],
                threshold = input.saint_threshold.get()
            )

            progress.set(message="Running protein feature analysis", detail="Generating plots...", value=80)
            # Store the results in the reactive value
            feature_enrichment.set(result)

            progress.set(message="Running protein feature analysis", detail="Saving results...", value=90)
            # Save the results to the dataset directory
            result.to_csv(os.path.join(out_dir, input.feature_dataset.get(), "Feature_enrichment.csv"), index=False)
            
            progress.set(message="Running protein feature analysis", detail="Done!", value=100)
            print("Done running feature analysis. Results saved.")

    @render.plot
    def feature_enrichment_plot():
        
        # Set the feature enrichment to whatever dataset is selected
        dataset = input.feature_dataset.get()
        if not dataset:
            return None
        
        feature_file = os.path.join(out_dir, dataset, "Feature_enrichment.csv")
        if not os.path.exists(feature_file):
            return None

        feature_enrichment.set(pd.read_csv(feature_file))
        
        feature_type = input.feature_type.get()
        if not feature_type:
            feature_type = 'GO_CC'  # Default feature type if none is selected
        # Plot the results for a specific feature type, e.g., 'Domains'

        num_features = input.num_features.get() if input.num_features.get() else 30  # Default to 30 if not set
        if num_features < 1:
            num_features = 30

        # Call the plot_results function to generate the heatmap
        heatmap = plot_results(feature_enrichment.get(), feature_type, num_features)

        return heatmap
    


    @reactive.Calc
    def available_datasets():
        return datasets.get()["Dataset Name"].tolist()

    @reactive.Calc
    def scored_datasets():
        return datasets.get()[datasets.get()['Scored'] == 'Yes']["Dataset Name"].tolist()

    @reactive.Effect
    def update_score_dataset():
        # Update the dropdown choices dynamically
        available_choices = available_datasets()
        scored_choices = scored_datasets()
        ui.update_select("score_dataset", choices=available_choices)
        ui.update_select("qc_dataset", choices=available_choices)
        ui.update_select("download_dataset", choices=available_choices)
        ui.update_select("feature_dataset", choices=scored_choices)

    @reactive.Effect
    def update_qc_bait():
        # Update the bait dropdown choices dynamically based on the selected dataset
        dataset_name = input.qc_dataset.get()
        if dataset_name:
            # Check to see if the dataset has been scored
            try:
                scores = pd.read_csv(os.path.join(out_dir, dataset_name, "annotated_scores.csv"))
                baits = scores['Experiment.ID'].unique().tolist()
                print(baits)
                baits.insert(0, "All")  # Add "All" option
                ui.update_select("qc_bait", choices=baits)
            except FileNotFoundError:
                print(f"File not found for dataset {dataset_name}. Resetting bait selection.")
                ui.update_select("qc_bait", choices=["All"])  # Reset to default if file not found
        else:
            ui.update_select("qc_bait", choices=["All"])  # Reset to default if no dataset is selected

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

    @render.data_frame
    def custom_table():
        dataset = input.download_dataset.get()
        if not dataset:
            return pd.DataFrame()
        
        # Custom cols will come from a user input
        custom_cols = list(input.custom_columns.get())
        print(custom_cols)
        if not custom_cols:
            # Default columns if none are selected
            custom_cols = ["Experiment.ID", "Prey.ID", "SaintScore", "BFDR"]
        
        # Load the dataset from the output directory
        dataset_path = os.path.join(out_dir, dataset, "annotated_scores.csv")
        if not os.path.exists(dataset_path):
            return pd.DataFrame()
        df = pd.read_csv(dataset_path)

        print(df.columns)

        custom_df = df[custom_cols]
        custom_dataset.set(custom_df)

        return custom_dataset.get()
        return custom_dataset.get()
    
    # Add a download button for the custom dataset
    @render.download()
    def download_custom_dataset():
        if custom_dataset.get().empty:
            ui.notification_show(
                "No custom dataset to download. Please select columns first.",
                type="error",
            )
            return None
        # Create a temporary file to save the custom dataset
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"custom_dataset_{timestamp}.csv"
        savepath = os.path.join(out_dir, filename)
        # temp_file = tempfile.NamedTemporaryFile(delete=False, prefix=filename, suffix=".csv")
        custom_dataset.get().to_csv(savepath, index=False)
        return savepath

app = App(app_ui, server)

if __name__ == "__main__":
    run_app(app, host="0.0.0.0", port=3838)
