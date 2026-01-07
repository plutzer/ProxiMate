"""
Network Comparison Module for ProxiMate

This module provides functions for comparing protein-protein interaction networks
between two baits, including data loading, filtering, statistical analysis, and visualization.
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.stats import ttest_ind
import os


def parse_intensity_string(intensity_str):
    """
    Parse intensity string format from ProxiMate data.

    Format: pipe-delimited values with "." representing missing values
    Example: "100.5|200.3|.|150.2" -> [100.5, 200.3, 150.2]

    Parameters:
    -----------
    intensity_str : str or None
        Pipe-delimited intensity string

    Returns:
    --------
    list of float
        Numeric intensity values (missing values excluded)
    """
    if pd.isnull(intensity_str):
        return []

    values = str(intensity_str).split('|')
    numeric_values = []
    for v in values:
        v = v.strip()
        if v != '.' and v != '':
            try:
                numeric_values.append(float(v))
            except ValueError:
                continue

    return numeric_values


def load_and_filter_bait_data(dataset_name, bait_name, thresholds, out_dir="/Outputs"):
    """
    Load and filter data for a specific bait based on thresholds.

    Parameters:
    -----------
    dataset_name : str
        Name of the dataset
    bait_name : str
        Name of the bait (Experiment.ID)
    thresholds : dict
        Dictionary with keys: 'SaintScore', 'BFDR', 'WD', 'WDFDR'
        Example: {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 0.05}
    out_dir : str
        Output directory path (default: "/Outputs")

    Returns:
    --------
    pd.DataFrame
        Filtered dataframe containing only rows passing all thresholds
        Key columns: Prey.ID, First_Prey_Gene, FoldChange, AvgIntensity,
                    SaintScore, BFDR, WD, WDFDR
    """
    results_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

    if not os.path.exists(results_path):
        return pd.DataFrame()

    # Load data
    data = pd.read_csv(results_path)

    # Filter for specific bait
    bait_data = data[data['Experiment.ID'] == bait_name].copy()

    if len(bait_data) == 0:
        return pd.DataFrame()

    # Apply thresholds (AND logic)
    # Handle NaN values in WDFDR (when n_iterations=0)
    filtered = bait_data[
        (bait_data['SaintScore'] >= thresholds['SaintScore']) &
        (bait_data['BFDR'] <= thresholds['BFDR']) &
        (bait_data['WD'] >= thresholds['WD']) &
        (bait_data['WDFDR'].fillna(1.0) <= thresholds['WDFDR'])  # Treat NaN as 1.0 (not passing)
    ]

    return filtered


def calculate_volcano_data(dataset_name, bait_a, bait_b, thresholds_a, thresholds_b, out_dir="/Outputs"):
    """
    Calculate volcano plot data comparing two baits from the same dataset using RAW interaction data.

    For each prey detected in BOTH baits:
    - Calculate mean intensities from raw interaction.txt file
    - Calculate direct log2 fold change: log2(mean_intensity_A / mean_intensity_B)
    - Perform two-sample t-test on test replicate intensities
    - Calculate -log10(p-value)
    - Determine category based on threshold passing from annotated_scores.csv

    Parameters:
    -----------
    dataset_name : str
        Name of the dataset (must be same for both baits)
    bait_a : str
        Bait name for network A
    bait_b : str
        Bait name for network B
    thresholds_a : dict
        Thresholds for bait A {'SaintScore', 'BFDR', 'WD', 'WDFDR'}
    thresholds_b : dict
        Thresholds for bait B
    out_dir : str
        Output directory path (default: "/Outputs")

    Returns:
    --------
    pd.DataFrame
        Columns: Prey.ID, First_Prey_Gene, log2_fc_ratio, neg_log10_pval,
                category ('Both', 'Network A only', 'Network B only', 'Neither'),
                fc_a, fc_b, pval
    """
    interaction_path = os.path.join(out_dir, dataset_name, "interaction.txt")
    ed_path = os.path.join(out_dir, dataset_name, "ED.csv")
    scores_path = os.path.join(out_dir, dataset_name, "annotated_scores.csv")

    if not os.path.exists(interaction_path) or not os.path.exists(ed_path):
        return pd.DataFrame()

    # Load interaction file (tab-delimited, no header)
    # Columns: Experiment Name, Bait, Prey ID, Intensity
    interaction = pd.read_csv(interaction_path, sep="\t", header=None,
                             names=['Experiment', 'Bait', 'Prey.ID', 'Intensity'])

    # Load experimental design
    ed = pd.read_csv(ed_path)

    # Load annotated scores for threshold filtering and gene names
    if os.path.exists(scores_path):
        scores = pd.read_csv(scores_path)
    else:
        scores = None

    # Get test and control experiment names from ED
    test_experiments_a = ed[(ed['Type'] == 'T') & (ed['Bait'] == bait_a)]['Experiment Name'].tolist()
    test_experiments_b = ed[(ed['Type'] == 'T') & (ed['Bait'] == bait_b)]['Experiment Name'].tolist()
    control_experiments = ed[ed['Type'] == 'C']['Experiment Name'].tolist()

    if len(test_experiments_a) == 0 or len(test_experiments_b) == 0:
        return pd.DataFrame()

    # Get interaction data for each bait
    data_a_test = interaction[interaction['Experiment'].isin(test_experiments_a)]
    data_b_test = interaction[interaction['Experiment'].isin(test_experiments_b)]
    data_control = interaction[interaction['Experiment'].isin(control_experiments)]

    # Find common preys between baits A and B
    preys_a = set(data_a_test['Prey.ID'].unique())
    preys_b = set(data_b_test['Prey.ID'].unique())
    common_preys = preys_a & preys_b

    if len(common_preys) == 0:
        return pd.DataFrame()

    results = []
    for prey_id in common_preys:
        # Get test intensities for bait A
        intensities_a_test = data_a_test[data_a_test['Prey.ID'] == prey_id]['Intensity'].values
        # Get test intensities for bait B
        intensities_b_test = data_b_test[data_b_test['Prey.ID'] == prey_id]['Intensity'].values

        # Calculate mean intensities
        mean_a_test = np.mean(intensities_a_test) if len(intensities_a_test) > 0 else 0
        mean_b_test = np.mean(intensities_b_test) if len(intensities_b_test) > 0 else 0

        # Calculate direct fold change (A/B)
        if mean_a_test > 0 and mean_b_test > 0:
            fc_ratio = mean_a_test / mean_b_test
            log2_ratio = np.log2(fc_ratio)
        else:
            fc_ratio = 0
            log2_ratio = 0.0

        # Perform t-test comparing test intensities from A vs B
        if len(intensities_a_test) >= 2 and len(intensities_b_test) >= 2:
            try:
                _, pval = ttest_ind(intensities_a_test, intensities_b_test)
            except:
                pval = 1.0
        else:
            pval = 1.0

        # Calculate -log10(p-value)
        if pval == 0:
            neg_log10_pval = 300
        else:
            neg_log10_pval = -np.log10(pval)

        # Get gene name from scores if available
        gene_name = prey_id  # Default to Prey.ID
        if scores is not None:
            gene_row = scores[scores['Prey.ID'] == prey_id]
            if len(gene_row) > 0 and 'First_Prey_Gene' in gene_row.columns:
                gene_name = gene_row['First_Prey_Gene'].iloc[0]

        # Determine category based on threshold passing from annotated scores
        category = 'Neither'
        if scores is not None:
            # Get scored data for both baits
            score_a = scores[(scores['Experiment.ID'] == bait_a) & (scores['Prey.ID'] == prey_id)]
            score_b = scores[(scores['Experiment.ID'] == bait_b) & (scores['Prey.ID'] == prey_id)]

            passes_a = False
            passes_b = False

            if len(score_a) > 0:
                row_a = score_a.iloc[0]
                passes_a = (
                    row_a['SaintScore'] >= thresholds_a['SaintScore'] and
                    row_a['BFDR'] <= thresholds_a['BFDR'] and
                    row_a['WD'] >= thresholds_a['WD'] and
                    (pd.isna(row_a['WDFDR']) or row_a['WDFDR'] <= thresholds_a['WDFDR'])
                )

            if len(score_b) > 0:
                row_b = score_b.iloc[0]
                passes_b = (
                    row_b['SaintScore'] >= thresholds_b['SaintScore'] and
                    row_b['BFDR'] <= thresholds_b['BFDR'] and
                    row_b['WD'] >= thresholds_b['WD'] and
                    (pd.isna(row_b['WDFDR']) or row_b['WDFDR'] <= thresholds_b['WDFDR'])
                )

            # Determine category
            if passes_a and passes_b:
                category = 'Both'
            elif passes_a:
                category = 'Network A only'
            elif passes_b:
                category = 'Network B only'

        results.append({
            'Prey.ID': prey_id,
            'First_Prey_Gene': gene_name,
            'log2_fc_ratio': log2_ratio,
            'neg_log10_pval': neg_log10_pval,
            'category': category,
            'mean_intensity_a': mean_a_test,
            'mean_intensity_b': mean_b_test,
            'fc_ratio': fc_ratio,
            'pval': pval
        })

    return pd.DataFrame(results)


def create_volcano_plot(volcano_data, bait_a, bait_b):
    """
    Create interactive volcano plot for network comparison.

    Parameters:
    -----------
    volcano_data : pd.DataFrame
        Output from calculate_volcano_data()
    bait_a : str
        Name of bait A
    bait_b : str
        Name of bait B

    Returns:
    --------
    plotly.graph_objects.Figure
        Interactive volcano plot
    """
    if len(volcano_data) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No common prey proteins between selected baits",
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

    fig = go.Figure()

    # Define colors for categories
    category_colors = {
        'Both': '#9467bd',  # Purple
        'Network A only': '#ff7f0e',  # Orange
        'Network B only': '#2ca02c',  # Green
        'Neither': '#d3d3d3'  # Light gray
    }

    # Plot each category separately for legend control
    # Plot in reverse order so important categories are on top
    for category in ['Neither', 'Network B only', 'Network A only', 'Both']:
        cat_data = volcano_data[volcano_data['category'] == category]

        if len(cat_data) == 0:
            continue

        # Create hover text
        hover_text = []
        for _, row in cat_data.iterrows():
            text = (
                f"<b>{row['First_Prey_Gene']}</b><br>"
                f"log2(FC ratio A/B): {row['log2_fc_ratio']:.2f}<br>"
                f"FC ratio (A/B): {row['fc_ratio']:.2f}<br>"
                f"Mean intensity {bait_a}: {row['mean_intensity_a']:.2e}<br>"
                f"Mean intensity {bait_b}: {row['mean_intensity_b']:.2e}<br>"
                f"-log10(p-value): {row['neg_log10_pval']:.2f}<br>"
                f"p-value: {row['pval']:.2e}"
            )
            hover_text.append(text)

        fig.add_trace(go.Scatter(
            x=cat_data['log2_fc_ratio'],
            y=cat_data['neg_log10_pval'],
            mode='markers',
            marker=dict(
                size=8,
                color=category_colors[category],
                line=dict(width=0.5, color='white')
            ),
            text=hover_text,
            hovertemplate='%{text}<extra></extra>',
            name=category
        ))

    # Add reference lines
    # Vertical line at x=0 (no change)
    fig.add_vline(x=0, line_dash="dash", line_color="gray", line_width=1)

    # Horizontal line at p=0.05 (-log10(0.05) = 1.3)
    fig.add_hline(y=-np.log10(0.05), line_dash="dot", line_color="red",
                  line_width=2, annotation_text="p = 0.05",
                  annotation_position="right")

    # Update layout
    fig.update_layout(
        title=f"Network Comparison: {bait_a} vs {bait_b}",
        xaxis_title=f"log2(Mean Intensity {bait_a} / Mean Intensity {bait_b})",
        yaxis_title="-log10(p-value)",
        xaxis=dict(zeroline=True, zerolinewidth=1, zerolinecolor='lightgray'),
        yaxis=dict(zeroline=False),
        hovermode='closest',
        height=500,
        template='plotly_white',
        legend=dict(
            title="Category",
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        )
    )

    return fig


def create_venn_diagram(set_a, set_b, label_a, label_b):
    """
    Create Venn diagram using plotly shapes (no matplotlib_venn dependency).

    Parameters:
    -----------
    set_a : set
        Set of Prey.IDs passing thresholds for bait A
    set_b : set
        Set of Prey.IDs passing thresholds for bait B
    label_a : str
        Label for set A (bait name)
    label_b : str
        Label for set B (bait name)

    Returns:
    --------
    plotly.graph_objects.Figure
        Venn diagram visualization
    """
    # Calculate set sizes
    only_a = len(set_a - set_b)
    only_b = len(set_b - set_a)
    both = len(set_a & set_b)

    fig = go.Figure()

    # Draw two circles using shapes
    # Circle A (left)
    fig.add_shape(
        type="circle",
        xref="x", yref="y",
        x0=0.2, y0=0.2, x1=1.2, y1=1.2,
        line_color="#ff7f0e",  # Orange
        line_width=3,
        fillcolor="rgba(255, 127, 14, 0.2)"
    )

    # Circle B (right)
    fig.add_shape(
        type="circle",
        xref="x", yref="y",
        x0=0.8, y0=0.2, x1=1.8, y1=1.2,
        line_color="#2ca02c",  # Green
        line_width=3,
        fillcolor="rgba(44, 160, 44, 0.2)"
    )

    # Add text annotations for counts
    # Left section (only A)
    fig.add_annotation(
        x=0.5, y=0.7,
        text=f"<b>{only_a}</b>",
        showarrow=False,
        font=dict(size=24, color="black")
    )

    # Intersection
    fig.add_annotation(
        x=1.0, y=0.7,
        text=f"<b>{both}</b>",
        showarrow=False,
        font=dict(size=24, color="black")
    )

    # Right section (only B)
    fig.add_annotation(
        x=1.5, y=0.7,
        text=f"<b>{only_b}</b>",
        showarrow=False,
        font=dict(size=24, color="black")
    )

    # Labels
    fig.add_annotation(
        x=0.5, y=1.5,
        text=f"<b>{label_a}</b>",
        showarrow=False,
        font=dict(size=16, color="#ff7f0e")
    )

    fig.add_annotation(
        x=1.5, y=1.5,
        text=f"<b>{label_b}</b>",
        showarrow=False,
        font=dict(size=16, color="#2ca02c")
    )

    # Update layout
    fig.update_xaxes(range=[0, 2], visible=False)
    fig.update_yaxes(range=[0, 2], visible=False)
    fig.update_layout(
        height=400,
        template='plotly_white',
        showlegend=False,
        margin=dict(l=20, r=20, t=20, b=20)
    )

    return fig
