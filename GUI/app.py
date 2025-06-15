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
from QC_plots import pca_plot, saint_known_retention
from Ann_Enrichment import process_refactored, plot_results


out_dir = "/Outputs"

app_ui = ui.page_navbar(
    ui.nav_spacer(),
    ui.nav_panel("Input Data",
                 "Use either method to parse input data:",
                 ui.layout_columns(
                    ui.card(
                        ui.card_header("MaxQuant & Experimental Design Input"),
                        ui.input_text("mq_dataset_name", "Dataset Name"),
                        ui.input_file("pg_file", "MaxQuant proteinGroups.txt file"),
                        ui.input_file("ed_file", "Experimental Design File"),
                        ui.input_select("mq_quant_type", "Quantification Type", choices=["Intensity", "LFQ", "Spectral Counts"]),
                        "Experimental Design File Format:",
                        "SampleName  Condition BiologicalReplicate",
                        ui.input_action_button("parse_mq", "Parse MaxQuant Inputs")
                        
                    ),
                    ui.card(
                        ui.card_header('SAINT Input'),
                        ui.input_text("st_dataset_name", "Dataset Name"),
                        ui.input_file("bait", "SAINT bait.txt file"),
                        ui.input_file("prey", "SAINT prey.txt file"),
                        ui.input_file("interaction", "SAINT interaction.txt file"),
                        # Add a link to the SAINT documentation for input formats
                        # ui.link("Documentation for SAINT input format", "www.google.com"),# Placeholder
                        ui.input_action_button("parse_saint", "Parse SAINT Inputs")
                    ),
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
                    ui.input_select("qc_dataset", "Select Dataset", choices=[]), # Need this to be dynamic
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
                    ui.input_select("download_dataset", "Select Dataset", choices=[]), # Need this to be dynamic
                    # Add download buttons for all the files
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

    @reactive.effect
    @reactive.event(input.parse_mq)
    def parse_mq():
        print("Parsing")
        print(input.pg_file.get()[0]['name'])
        print(input.ed_file.get()[0]['name'])
        curr_dataset = datasets.get()

        # Code for parsing goes here
        print('parsing...')
        # Print the current directory
        print("Current directory:", os.getcwd())

        n_exp, n_ctrl = parse.parse_ed_pg(input.pg_file.get()[0]['datapath'], input.ed_file.get()[0]['datapath'], input.mq_quant_type.get(), out_dir + '/' + input.mq_dataset_name.get())

        # Update the datasets dataframe
        new_row = pd.DataFrame([[input.mq_dataset_name.get(), f'Protein Groups', input.mq_quant_type.get(), n_exp, n_ctrl, '', '', '']], columns=datasets.get().columns)
        updated_datasets = pd.concat([datasets.get(), new_row], ignore_index=True)
        datasets.set(updated_datasets)

    @reactive.effect
    @reactive.event(input.parse_saint)
    def parse_saint():
        print("Parsing")
        print(input.bait)
        print(input.prey)
        print(input.interaction)

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
        print("Scoring data")
        print(input.score_dataset.get())
        print(input.imputation_method.get())
        print(input.wdfdr_iterations.get())

        # Get the quant type from the datasets dataframe
        quant_type = datasets.get().loc[datasets.get()['Dataset Name'] == input.score_dataset.get(), 'Quant Type'].values[0]

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

    # Quality controls tab
    @render_plotly
    def raw_pca_plot():
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

        if not os.path.exists(results_path):
            return None
        else:
            # Call the saint_known_retention function to generate the plot
            fig = saint_known_retention(results_path)

            # Return the figure widget
            return fig

    feature_enrichment = reactive.Value(pd.DataFrame())

    # Feature analysis tab
    @reactive.effect
    @reactive.event(input.feature_analysis)
    def feature_analysis():
        print("Running feature analysis")
        print(input.feature_dataset.get())
        print(input.saint_threshold.get())

        # Get the selected dataset
        dataset = pd.read_csv(os.path.join(out_dir, input.feature_dataset.get(), "annotated_scores.csv"))

        result = process_refactored(
            dataset,
            columns_for_analysis = ['GO_CC', 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains'],
            threshold = input.saint_threshold.get()
        )

        # Store the results in the reactive value
        feature_enrichment.set(result)

        # Save the results to the dataset directory
        result.to_csv(os.path.join(out_dir, input.feature_dataset.get(), "Feature_enrichment.csv"), index=False)
        
        print("Done running feature analysis. Results saved.")

    @render.plot
    def feature_enrichment_plot():
        # Check if the feature enrichment data is available
        if feature_enrichment.get().empty:
            return None
        
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
        # Extract the "Dataset Name" column
        return datasets.get()["Dataset Name"].tolist()

    @reactive.Effect
    def update_score_dataset():
        # Update the dropdown choices dynamically
        available_choices = available_datasets()
        ui.update_select("score_dataset", choices=available_choices)
        ui.update_select("qc_dataset", choices=available_choices)
        ui.update_select("download_dataset", choices=available_choices)
        ui.update_select("feature_dataset", choices=available_choices)

app = App(app_ui, server)

if __name__ == "__main__":
    run_app(app, host="0.0.0.0", port=3838)
