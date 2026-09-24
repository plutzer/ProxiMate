"""
Matplotlib-based plot export functions for ProxiMate.

These functions create static matplotlib versions of the interactive Plotly plots
for PNG/SVG export, avoiding the kaleido dependency.
"""

import numpy as np
from matplotlib import colormaps
from matplotlib.figure import Figure

from QC_plots import (KNOWN_STATUS_STYLE, bait_scores, experiment_pca,
                      known_status_split, prey_pca)


def pca_plot_matplotlib(interaction, experimentalDesign, matrix=None):
    """
    Create a matplotlib PCA plot for export.

    Parameters:
    -----------
    interaction : str
        Path to interaction.txt file
    experimentalDesign : str
        Path to ED.csv file
    matrix : pd.DataFrame, optional
        Preprocessed prey x experiment matrix from prepare_pca_matrix; computed
        with default settings when omitted.

    Returns:
    --------
    matplotlib.figure.Figure
        PCA scatter plot
    """
    pca_df, explained_variance = experiment_pca(interaction, experimentalDesign, matrix)

    # Create matplotlib figure
    fig = Figure(figsize=(10, 8))
    ax = fig.subplots()

    # Get unique baits and types for coloring/markers
    unique_baits = pca_df['BaitName'].unique()
    unique_types = pca_df['Type'].unique()

    # Color palette
    colors = colormaps["tab10"](np.linspace(0, 1, len(unique_baits)))
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
    fig.tight_layout()

    return fig


def prey_pca_matplotlib(matrix, color_values=None, color_label=None,
                        color_mode="none", color_threshold=None):
    """
    Create a matplotlib prey-level PCA plot for export.

    Parameters:
    -----------
    matrix : pd.DataFrame
        Preprocessed prey x experiment matrix from prepare_pca_matrix
    color_values : pd.Series, optional
        Per-prey color data indexed by prey ID (numeric for 'continuous',
        labels for 'categorical')
    color_label : str, optional
        Legend / colorbar title
    color_mode : str
        'none', 'continuous', or 'categorical'

    Returns:
    --------
    matplotlib.figure.Figure
        PCA scatter plot
    """
    prey_df, explained_variance = prey_pca(matrix)

    fig = Figure(figsize=(10, 8))
    ax = fig.subplots()

    if color_mode == "continuous":
        values = prey_df['Prey'].map(color_values)
        if color_threshold is not None:
            below = values < color_threshold
            ax.scatter(prey_df.loc[below, 'PC1'], prey_df.loc[below, 'PC2'],
                       c='lightgrey', s=20, alpha=0.7,
                       label=f"{color_label} < {color_threshold:g}")
            sc = ax.scatter(prey_df.loc[~below, 'PC1'], prey_df.loc[~below, 'PC2'],
                            c=values[~below], cmap='viridis', s=20, alpha=0.7)
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), fontsize=9)
        else:
            sc = ax.scatter(prey_df['PC1'], prey_df['PC2'], c=values,
                            cmap='viridis', s=20, alpha=0.7)
        fig.colorbar(sc, ax=ax, label=color_label)
    elif color_mode == "categorical":
        labels = prey_df['Prey'].map(color_values)
        categories = [c for c in labels.dropna().unique()
                      if c not in ("Other", "Unknown")]
        categories = sorted(categories) + ["Other", "Unknown"]
        colors = colormaps["tab20"](np.linspace(0, 1, len(categories)))
        for category, color in zip(categories, colors):
            subset = prey_df[labels == category]
            if len(subset) > 0:
                ax.scatter(subset['PC1'], subset['PC2'], c=[color],
                           s=20, alpha=0.7, label=category)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
                  ncol=4, fontsize=9, title=color_label)
    else:
        ax.scatter(prey_df['PC1'], prey_df['PC2'], s=20, alpha=0.7)

    ax.set_xlabel(f"PC1 ({explained_variance[0]*100:.2f}% variance)", fontsize=12)
    ax.set_ylabel(f"PC2 ({explained_variance[1]*100:.2f}% variance)", fontsize=12)
    ax.set_title("Prey PCA", fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

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
    bait_data = bait_scores(results_path, bait_name)

    fig = Figure(figsize=(10, 7))
    ax = fig.subplots()

    if len(bait_data) == 0:
        ax.text(0.5, 0.5, f"No data available for bait: {bait_name}",
                ha='center', va='center', fontsize=14, color='red',
                transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        return fig

    for zorder, (group, (name, color)) in enumerate(
            zip(known_status_split(bait_data), KNOWN_STATUS_STYLE), start=1):
        if len(group) > 0:
            ax.scatter(group['FoldChange'], group['SaintScore'],
                       c=color, s=50, alpha=0.7, edgecolors='white', linewidth=0.5,
                       label=name, zorder=zorder)

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

    fig.tight_layout()

    return fig
