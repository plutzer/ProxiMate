"""
Matplotlib-based plot export functions for ProxiMate.

These functions create static matplotlib versions of the interactive Plotly plots
for PNG/SVG export, avoiding the kaleido dependency.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


def pca_plot_matplotlib(interaction, experimentalDesign):
    """
    Create a matplotlib PCA plot for export.

    Parameters:
    -----------
    interaction : str
        Path to interaction.txt file
    experimentalDesign : str
        Path to ED.csv file

    Returns:
    --------
    matplotlib.figure.Figure
        PCA scatter plot
    """
    # Load data
    int_df = pd.read_csv(interaction, sep="\t", header=0)
    int_df.columns = ['Experiment', 'BaitName', 'Prey', 'Intensity']
    ed = pd.read_csv(experimentalDesign, sep=",")

    # Make the int table wide
    data = int_df.pivot(index='Prey', columns='Experiment', values='Intensity')

    metadata = int_df[['Experiment', 'BaitName']].drop_duplicates()
    metadata = metadata.merge(ed[['Experiment Name', 'Type']],
                              left_on='Experiment', right_on='Experiment Name', how='left')

    # Clean up the data
    data = data.replace(0, np.nan)
    data = data.dropna(thresh=len(data.columns) * 0.5)
    data = data.apply(lambda row: row.fillna(row.min()), axis=1)
    data = data.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(data.T)

    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df['Experiment'] = data.columns
    pca_df = pca_df.merge(metadata, left_on='Experiment', right_on='Experiment', how='left')

    explained_variance = pca.explained_variance_ratio_

    # Create matplotlib figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Get unique baits and types for coloring/markers
    unique_baits = pca_df['BaitName'].unique()
    unique_types = pca_df['Type'].unique()

    # Color palette
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_baits)))
    color_map = dict(zip(unique_baits, colors))

    # Marker map for types
    markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', 'h']
    marker_map = dict(zip(unique_types, markers[:len(unique_types)]))

    # Plot each combination
    for bait in unique_baits:
        for exp_type in unique_types:
            subset = pca_df[(pca_df['BaitName'] == bait) & (pca_df['Type'] == exp_type)]
            if len(subset) > 0:
                ax.scatter(subset['PC1'], subset['PC2'],
                          c=[color_map[bait]],
                          marker=marker_map[exp_type],
                          s=80, alpha=0.8, edgecolors='white', linewidth=0.5,
                          label=f"{bait} ({exp_type})")

    ax.set_xlabel(f"PC1 ({explained_variance[0]*100:.2f}% variance)", fontsize=12)
    ax.set_ylabel(f"PC2 ({explained_variance[1]*100:.2f}% variance)", fontsize=12)
    ax.set_title("PCA of Interaction Data", fontsize=14, fontweight='bold')

    # Legend below plot
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
              ncol=min(4, len(unique_baits)), fontsize=9)

    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    return fig


def saint_scatter_matplotlib(results_path, bait_name, saintscore_threshold):
    """
    Create a matplotlib scatter plot of SAINT Score vs Fold Change for export.

    Parameters:
    -----------
    results_path : str
        Path to annotated_scores.csv file
    bait_name : str
        Name of the bait to visualize
    saintscore_threshold : float
        Threshold value to draw as horizontal reference line

    Returns:
    --------
    matplotlib.figure.Figure
        Scatter plot
    """
    # Load data
    results = pd.read_csv(results_path, sep=",")

    # Filter for specific bait
    bait_data = results[results['Experiment.ID'] == bait_name].copy()

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))

    if len(bait_data) == 0:
        ax.text(0.5, 0.5, f"No data available for bait: {bait_name}",
                ha='center', va='center', fontsize=14, color='red',
                transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        return fig

    # Check available columns
    has_biogrid = 'In.BioGRID' in bait_data.columns
    has_multivalidated = 'Multivalidated' in bait_data.columns

    # Separate data by BioGRID status
    if has_multivalidated:
        multivalidated = bait_data[bait_data['Multivalidated'] == True].copy()
        in_biogrid = bait_data[(bait_data['In.BioGRID'] == True) & (bait_data['Multivalidated'] != True)].copy()
        not_in_biogrid = bait_data[bait_data['In.BioGRID'] != True].copy()
    elif has_biogrid:
        multivalidated = pd.DataFrame()
        in_biogrid = bait_data[bait_data['In.BioGRID'] == True].copy()
        not_in_biogrid = bait_data[bait_data['In.BioGRID'] != True].copy()
    else:
        multivalidated = pd.DataFrame()
        in_biogrid = pd.DataFrame()
        not_in_biogrid = bait_data.copy()

    # Plot in order: not in BioGRID (background), in BioGRID, multivalidated (foreground)
    if len(not_in_biogrid) > 0:
        ax.scatter(not_in_biogrid['FoldChange'], not_in_biogrid['SaintScore'],
                   c='#1f77b4', s=50, alpha=0.7, edgecolors='white', linewidth=0.5,
                   label='Not in BioGRID', zorder=1)

    if len(in_biogrid) > 0:
        ax.scatter(in_biogrid['FoldChange'], in_biogrid['SaintScore'],
                   c='#ff7f0e', s=50, alpha=0.7, edgecolors='white', linewidth=0.5,
                   label='In BioGRID', zorder=2)

    if len(multivalidated) > 0:
        ax.scatter(multivalidated['FoldChange'], multivalidated['SaintScore'],
                   c='#d62728', s=50, alpha=0.7, edgecolors='white', linewidth=0.5,
                   label='Multivalidated', zorder=3)

    # Add threshold line
    ax.axhline(y=saintscore_threshold, color='red', linestyle='--', linewidth=2,
               label=f'Threshold: {saintscore_threshold}', zorder=4)

    # Add zero line for fold change
    ax.axvline(x=0, color='lightgray', linestyle='-', linewidth=1, zorder=0)

    ax.set_xlabel("Fold Change (log2)", fontsize=12)
    ax.set_ylabel("SAINT Score", fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.set_title(f"SAINT Score vs Fold Change - {bait_name}", fontsize=14, fontweight='bold')

    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig
