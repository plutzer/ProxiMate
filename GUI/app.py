from functools import partial
from shiny import App, Inputs, Outputs, Session, reactive, render, ui, run_app
import pandas as pd
import os
import parse
import subprocess

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
                        ui.layout_columns(
                            ui.h1('Parsed Datasets'),
                            ui.input_action_button("clear_datasets", "Clear All Datasets"),
                            ),
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
                    ui.card(
                        ui.card_header("Raw Data UMAP"),
                        # Plotly PCA plot goes here
                    ),
                    ui.card(
                        # Add dropdown to select bait
                        ui.card_header("Known Physical Interactions"),
                        # Plotly line plot goes here
                    )
    ),
    ui.nav_panel("Protein Domains",
                    "Page D content"
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
        columns=['Dataset Name', 'Input Type', 'Experiments', 'Controls', 'Scored', 'Imputation', 'WDFDR iterations']
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
        new_row = pd.DataFrame([[input.mq_dataset_name.get(), 'Protein Groups', n_exp, n_ctrl, '', '', '']], columns=datasets.get().columns)
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

    # Scoring data
    @reactive.effect
    @reactive.event(input.score_data)
    def score_data():
        print("Scoring data")
        print(input.score_dataset.get())
        print(input.imputation_method.get())
        print(input.wdfdr_iterations.get())

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
            input.mq_quant_type.get()
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

app = App(app_ui, server)

if __name__ == "__main__":
    run_app(app, host="0.0.0.0", port=3838)
