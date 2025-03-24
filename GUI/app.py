from functools import partial
from shiny import App, Inputs, Outputs, Session, reactive, render, ui, run_app
import pandas as pd

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
                    ui.input_radio_buttons("imputation_method", "Imputation Method", choices=["Default", "Prey-specific"]),
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
        columns=['Dataset Name', 'Input Type', 'Baits', 'Controls', 'Scored', 'Imputation', 'WDFDR iterations']
    ))

    # Function to render the datasets table
    @render.data_frame
    def render_datasets():
        return render.DataGrid(datasets.get())

    @reactive.effect
    @reactive.event(input.parse_mq)
    def parse_mq():
        print("Parsing")
        print(input.pg_file)
        print(input.ed_file)
        curr_dataset = datasets.get()

        # Code for parsing goes here
        print('parsing...')

        # Update the datasets dataframe
        new_row = pd.DataFrame([[input.mq_dataset_name.get(), 'Protein Groups', '', '', '', '', '']], columns=datasets.get().columns)
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