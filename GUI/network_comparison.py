"""
Network Comparison Module for ProxiMate

This module provides functions for comparing protein-protein interaction networks
between two baits, including data loading, filtering, statistical analysis, and visualization.
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import os
from QC_plots import apply_score_thresholds


# Side-panel BFDR values of 0 are floored here so -log10 stays finite; the cap
# puts them at y = 4, above any BFDR that survives rounding in SAINT output.
BFDR_FLOOR = 1e-4

CATEGORY_COLORS = {
    'Both': '#9467bd',            # Purple
    'Network A only': '#ff7f0e',  # Orange
    'Network B only': '#2ca02c',  # Green
    'Neither': '#d3d3d3'          # Light gray
}

# Draw order: background categories first so hits render on top.
CATEGORY_ORDER = ['Neither', 'Network B only', 'Network A only', 'Both']


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

    # Apply thresholds (AND logic) using centralized function
    filtered = apply_score_thresholds(bait_data, thresholds)

    return filtered


def calculate_volcano_data(dataset_name, bait_a, bait_b, thresholds_a, thresholds_b, out_dir="/Outputs"):
    """
    Calculate volcano plot data comparing two baits from the same dataset using RAW interaction data.

    Every prey detected under either bait's test experiments is returned, tagged by
    a `status` column:

    - 'shared' (nonzero mean intensity under both baits): direct log2 fold change of
      the mean intensities, plus a two-sample t-test on the log2 replicate
      intensities (zero intensities are non-detections and are excluded).  The raw
      p-values are Benjamini-Hochberg adjusted across all tested preys; the plotted
      `neg_log10_pval` is -log10 of the adjusted value.  Preys with fewer than two
      nonzero replicates on either side carry NaN p-values (a placeholder 1.0 would
      distort the BH ranking) and plot at y = 0.
    - 'a_only' / 'b_only' (detected under one bait only): no fold change or p-value
      exists; these rows carry `neg_log10_bfdr`, the -log10 SAINT BFDR under the
      bait where the prey is present, for the flanking presence/absence panels.
      A BFDR of 0 is floored at BFDR_FLOOR; a missing BFDR plots at 0.

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
        Columns: Prey.ID, First_Prey_Gene, status, category
                ('Both', 'Network A only', 'Network B only', 'Neither'),
                mean_intensity_a, mean_intensity_b, fc_ratio, log2_fc_ratio,
                pval, pval_adj, neg_log10_pval, neg_log10_bfdr
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

    # Load annotated scores for threshold filtering, gene names and BFDRs
    if os.path.exists(scores_path):
        scores = pd.read_csv(scores_path)
    else:
        scores = None

    # Get test experiment names from ED
    test_experiments_a = ed[(ed['Type'] == 'T') & (ed['Bait'] == bait_a)]['Experiment Name'].tolist()
    test_experiments_b = ed[(ed['Type'] == 'T') & (ed['Bait'] == bait_b)]['Experiment Name'].tolist()

    if len(test_experiments_a) == 0 or len(test_experiments_b) == 0:
        return pd.DataFrame()

    # Per-prey replicate intensity arrays for each bait
    data_a_test = interaction[interaction['Experiment'].isin(test_experiments_a)]
    data_b_test = interaction[interaction['Experiment'].isin(test_experiments_b)]
    intensities_a = data_a_test.groupby('Prey.ID')['Intensity'].apply(np.asarray).to_dict()
    intensities_b = data_b_test.groupby('Prey.ID')['Intensity'].apply(np.asarray).to_dict()

    all_preys = set(intensities_a) | set(intensities_b)
    if len(all_preys) == 0:
        return pd.DataFrame()

    # Threshold passing comes from the shared helper so the volcano's categories agree
    # with the filtered networks shown elsewhere.  It reads a NaN WDFDR as 1.0, which
    # fails the threshold; scoring with 0 CompPASS iterations produces those NaNs.
    gene_names = {}
    bfdr_a = {}
    bfdr_b = {}
    if scores is not None:
        passing_a = set(apply_score_thresholds(
            scores[scores['Experiment.ID'] == bait_a], thresholds_a)['Prey.ID'])
        passing_b = set(apply_score_thresholds(
            scores[scores['Experiment.ID'] == bait_b], thresholds_b)['Prey.ID'])
        if 'First_Prey_Gene' in scores.columns:
            gene_names = scores.drop_duplicates('Prey.ID').set_index('Prey.ID')['First_Prey_Gene'].to_dict()
        if 'BFDR' in scores.columns:
            bfdr_a = scores[scores['Experiment.ID'] == bait_a].set_index('Prey.ID')['BFDR'].to_dict()
            bfdr_b = scores[scores['Experiment.ID'] == bait_b].set_index('Prey.ID')['BFDR'].to_dict()
    else:
        passing_a = set()
        passing_b = set()

    results = []
    for prey_id in all_preys:
        vals_a = intensities_a.get(prey_id, np.array([]))
        vals_b = intensities_b.get(prey_id, np.array([]))

        mean_a_test = float(np.mean(vals_a)) if len(vals_a) > 0 else 0.0
        mean_b_test = float(np.mean(vals_b)) if len(vals_b) > 0 else 0.0

        if mean_a_test > 0 and mean_b_test > 0:
            status = 'shared'
        elif mean_a_test > 0:
            status = 'a_only'
        elif mean_b_test > 0:
            status = 'b_only'
        else:
            # Only zero intensities under both baits: nothing to plot.
            continue

        fc_ratio = np.nan
        log2_ratio = np.nan
        pval = np.nan
        neg_log10_bfdr = np.nan

        if status == 'shared':
            fc_ratio = mean_a_test / mean_b_test
            log2_ratio = np.log2(fc_ratio)

            # Zero intensities are non-detections, not measurements of zero:
            # exclude them rather than feed log2(0) into the test.
            log_a = np.log2(vals_a[vals_a > 0])
            log_b = np.log2(vals_b[vals_b > 0])
            if len(log_a) >= 2 and len(log_b) >= 2:
                _, pval = ttest_ind(log_a, log_b)
        else:
            bfdr = (bfdr_a if status == 'a_only' else bfdr_b).get(prey_id, np.nan)
            if pd.isna(bfdr):
                neg_log10_bfdr = 0.0
            else:
                neg_log10_bfdr = -np.log10(max(bfdr, BFDR_FLOOR))

        # Determine category based on threshold passing from annotated scores
        passes_a = prey_id in passing_a
        passes_b = prey_id in passing_b
        if passes_a and passes_b:
            category = 'Both'
        elif passes_a:
            category = 'Network A only'
        elif passes_b:
            category = 'Network B only'
        else:
            category = 'Neither'

        results.append({
            'Prey.ID': prey_id,
            'First_Prey_Gene': gene_names.get(prey_id, prey_id),
            'status': status,
            'category': category,
            'mean_intensity_a': mean_a_test,
            'mean_intensity_b': mean_b_test,
            'fc_ratio': fc_ratio,
            'log2_fc_ratio': log2_ratio,
            'pval': pval,
            'neg_log10_bfdr': neg_log10_bfdr
        })

    volcano = pd.DataFrame(results)

    # BH adjustment across the preys that were actually tested; NaN p-values
    # (side-panel rows, too few replicates) stay NaN and do not enter the ranking.
    volcano['pval_adj'] = np.nan
    tested = volcano['pval'].notna()
    if tested.any():
        volcano.loc[tested, 'pval_adj'] = multipletests(
            volcano.loc[tested, 'pval'], method='fdr_bh')[1]

    # Plotted y for the central panel: -log10 adjusted p.  Shared preys without a
    # test sit at 0; an adjusted p of exactly 0 is capped at 300.
    volcano['neg_log10_pval'] = np.nan
    shared = volcano['status'] == 'shared'
    adj = volcano.loc[shared, 'pval_adj']
    volcano.loc[shared, 'neg_log10_pval'] = np.where(
        adj.isna(), 0.0, np.where(adj == 0, 300.0, -np.log10(adj)))

    columns = ['Prey.ID', 'First_Prey_Gene', 'status', 'category',
               'mean_intensity_a', 'mean_intensity_b', 'fc_ratio', 'log2_fc_ratio',
               'pval', 'pval_adj', 'neg_log10_pval', 'neg_log10_bfdr']
    return volcano[columns]


def _shared_hover_text(row, bait_a, bait_b):
    return (
        f"<b>{row['First_Prey_Gene']}</b><br>"
        f"log2(FC ratio A/B): {row['log2_fc_ratio']:.2f}<br>"
        f"FC ratio (A/B): {row['fc_ratio']:.2f}<br>"
        f"Mean intensity {bait_a}: {row['mean_intensity_a']:.2e}<br>"
        f"Mean intensity {bait_b}: {row['mean_intensity_b']:.2e}<br>"
        f"p-value (raw): {row['pval']:.2e}<br>"
        f"adjusted p (BH): {row['pval_adj']:.2e}"
    )


def _side_hover_text(row, bait, mean_column):
    return (
        f"<b>{row['First_Prey_Gene']}</b><br>"
        f"Detected only in {bait}<br>"
        f"Mean intensity {bait}: {row[mean_column]:.2e}<br>"
        f"-log10(BFDR): {row['neg_log10_bfdr']:.2f}<br>"
        f"Category: {row['category']}"
    )


def create_volcano_plot(volcano_data, bait_a, bait_b):
    """
    Create interactive volcano plot for network comparison.

    Three panels: the central volcano holds preys quantified under both baits
    (x = log2 fold change of means, y = -log10 BH-adjusted p); narrow flanking
    jitter strips hold the presence/absence preys (left = only in bait A,
    right = only in bait B) with y = -log10(BFDR) under the bait where present.

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
            text="No prey proteins detected under either bait",
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

    # Strip sides follow the fold-change axis: positive log2 FC means higher in
    # bait A, so A-only preys sit on the right and B-only preys on the left.
    fig = make_subplots(
        rows=1, cols=3, column_widths=[0.15, 0.7, 0.15],
        horizontal_spacing=0.06,
        subplot_titles=(f"Only in {bait_b}", "", f"Only in {bait_a}")
    )

    # Fixed jitter seed: the strip layout is stable across redraws of the same data.
    rng = np.random.default_rng(0)
    categories_in_legend = set()

    # Central volcano: shared preys
    shared = volcano_data[volcano_data['status'] == 'shared']
    for category in CATEGORY_ORDER:
        cat_data = shared[shared['category'] == category]
        if len(cat_data) == 0:
            continue

        hover_text = [_shared_hover_text(row, bait_a, bait_b)
                      for _, row in cat_data.iterrows()]
        categories_in_legend.add(category)
        fig.add_trace(go.Scatter(
            x=cat_data['log2_fc_ratio'],
            y=cat_data['neg_log10_pval'],
            mode='markers',
            marker=dict(
                size=8,
                color=CATEGORY_COLORS[category],
                line=dict(width=0.5, color='white')
            ),
            text=hover_text,
            hovertemplate='%{text}<extra></extra>',
            name=category,
            legendgroup=category
        ), row=1, col=2)

    # Flanking presence/absence strips
    side_panels = [(1, 'b_only', bait_b, 'mean_intensity_b'),
                   (3, 'a_only', bait_a, 'mean_intensity_a')]
    for col, status, bait, mean_column in side_panels:
        panel = volcano_data[volcano_data['status'] == status]
        for category in CATEGORY_ORDER:
            cat_data = panel[panel['category'] == category]
            if len(cat_data) == 0:
                continue

            hover_text = [_side_hover_text(row, bait, mean_column)
                          for _, row in cat_data.iterrows()]
            show_legend = category not in categories_in_legend
            categories_in_legend.add(category)
            fig.add_trace(go.Scatter(
                x=rng.uniform(-0.4, 0.4, len(cat_data)),
                y=cat_data['neg_log10_bfdr'],
                mode='markers',
                marker=dict(
                    size=8,
                    color=CATEGORY_COLORS[category],
                    line=dict(width=0.5, color='white')
                ),
                text=hover_text,
                hovertemplate='%{text}<extra></extra>',
                name=category,
                legendgroup=category,
                showlegend=show_legend
            ), row=1, col=col)

    # Reference lines on the central panel
    fig.add_vline(x=0, line_dash="dash", line_color="gray", line_width=1,
                  row=1, col=2)
    fig.add_hline(y=-np.log10(0.05), line_dash="dot", line_color="red",
                  line_width=2, annotation_text="adj. p = 0.05",
                  annotation_position="top right", row=1, col=2)

    fig.update_xaxes(title_text=f"log2(Mean Intensity {bait_a} / Mean Intensity {bait_b})",
                     zeroline=True, zerolinewidth=1, zerolinecolor='lightgray',
                     row=1, col=2)
    fig.update_yaxes(title_text="-log10(adjusted p-value)", row=1, col=2)
    for col in (1, 3):
        fig.update_xaxes(range=[-1, 1], showticklabels=False, row=1, col=col)
        fig.update_yaxes(title_text="-log10(BFDR)", row=1, col=col)

    fig.update_layout(
        # Centered so the title clears the left strip's subplot title
        title=dict(text=f"Network Comparison: {bait_a} vs {bait_b}",
                   x=0.5, xanchor='center'),
        hovermode='closest',
        height=500,
        template='plotly_white',
        # Horizontal legend below the panels: inside the plot area it would sit on
        # top of the narrow right-hand strip.
        legend=dict(
            title="Category",
            orientation="h",
            yanchor="top",
            y=-0.25,
            xanchor="center",
            x=0.5
        )
    )

    return fig


def create_venn_diagram_matplotlib(set_a, set_b, label_a, label_b):
    """
    Create an area-proportional Venn diagram with matplotlib-venn.

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
    matplotlib.figure.Figure
        Venn diagram visualization
    """
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    fig, ax = plt.subplots(figsize=(8, 6))

    # venn2 cannot draw two empty sets
    if len(set_a) == 0 and len(set_b) == 0:
        ax.text(0.5, 0.5, "No preys pass thresholds in either network",
                ha='center', va='center', fontsize=14, color='#666')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        return fig

    venn = venn2([set_a, set_b], set_labels=(label_a, label_b),
                 set_colors=('#ff7f0e', '#2ca02c'), alpha=0.35, ax=ax)

    # Empty regions have no label
    for region_id in ('10', '01', '11'):
        label = venn.get_label_by_id(region_id)
        if label is not None:
            label.set_fontsize(16)
            label.set_fontweight('bold')

    for label, color in zip(venn.set_labels, ('#ff7f0e', '#2ca02c')):
        if label is not None:
            label.set_fontsize(14)
            label.set_fontweight('bold')
            label.set_color(color)

    plt.tight_layout()
    return fig


def create_volcano_plot_matplotlib(volcano_data, bait_a, bait_b):
    """
    Create a matplotlib volcano plot for export.

    Same three-panel layout as create_volcano_plot: central volcano for shared
    preys, flanking jitter strips for presence/absence preys.

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
    matplotlib.figure.Figure
        Volcano plot
    """
    import matplotlib.pyplot as plt

    if len(volcano_data) == 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, "No prey proteins detected under either bait",
                ha='center', va='center', fontsize=14, color='#666',
                transform=ax.transAxes)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        return fig

    fig = plt.figure(figsize=(13, 8))
    grid = fig.add_gridspec(1, 3, width_ratios=[0.15, 0.7, 0.15], wspace=0.35)
    ax_left = fig.add_subplot(grid[0, 0])
    ax_center = fig.add_subplot(grid[0, 1])
    ax_right = fig.add_subplot(grid[0, 2])

    # Same fixed jitter seed as the interactive plot
    rng = np.random.default_rng(0)

    # Central volcano: shared preys
    shared = volcano_data[volcano_data['status'] == 'shared']
    for category in CATEGORY_ORDER:
        cat_data = shared[shared['category'] == category]
        if len(cat_data) == 0:
            continue

        zorder = CATEGORY_ORDER.index(category) + 1
        ax_center.scatter(cat_data['log2_fc_ratio'], cat_data['neg_log10_pval'],
                          c=CATEGORY_COLORS[category], s=50, alpha=0.7,
                          edgecolors='white', linewidth=0.5,
                          label=category, zorder=zorder)

    ax_center.axvline(x=0, color='gray', linestyle='--', linewidth=1, zorder=0)
    ax_center.axhline(y=-np.log10(0.05), color='red', linestyle=':', linewidth=2, zorder=0)
    ax_center.annotate('adj. p = 0.05', xy=(ax_center.get_xlim()[1], -np.log10(0.05)),
                       xytext=(-5, 5), textcoords='offset points',
                       fontsize=10, color='red', va='bottom', ha='right')

    ax_center.set_xlabel(f"log2(Mean Intensity {bait_a} / Mean Intensity {bait_b})", fontsize=12)
    ax_center.set_ylabel("-log10(adjusted p-value)", fontsize=12)
    ax_center.set_title(f"Network Comparison: {bait_a} vs {bait_b}", fontsize=14, fontweight='bold')
    ax_center.legend(title="Category", loc='upper right', fontsize=9)
    ax_center.grid(True, alpha=0.3)

    # Flanking presence/absence strips, sides matching the fold-change axis:
    # positive log2 FC means higher in bait A, so A-only preys sit on the right.
    side_panels = [(ax_left, 'b_only', bait_b), (ax_right, 'a_only', bait_a)]
    for ax, status, bait in side_panels:
        panel = volcano_data[volcano_data['status'] == status]
        for category in CATEGORY_ORDER:
            cat_data = panel[panel['category'] == category]
            if len(cat_data) == 0:
                continue

            zorder = CATEGORY_ORDER.index(category) + 1
            ax.scatter(rng.uniform(-0.4, 0.4, len(cat_data)),
                       cat_data['neg_log10_bfdr'],
                       c=CATEGORY_COLORS[category], s=50, alpha=0.7,
                       edgecolors='white', linewidth=0.5, zorder=zorder)

        ax.set_xlim(-1, 1)
        ax.set_xticks([])
        ax.set_title(f"Only in {bait}", fontsize=11)
        ax.set_ylabel("-log10(BFDR)", fontsize=12)
        ax.grid(True, alpha=0.3, axis='y')

    return fig
