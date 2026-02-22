import plotly.express as px
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.metrics import roc_curve, auc


def apply_score_thresholds(df, thresholds):
    """
    Filter DataFrame by scoring thresholds.

    Applies combined AND logic: all conditions must pass for a row to be included.
    Handles NaN values in WDFDR (from scoring with 0 iterations) by treating them
    as 1.0 (failing the threshold).

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with columns: SaintScore, BFDR, WD, WDFDR
    thresholds : dict
        Dictionary with keys: 'SaintScore', 'BFDR', 'WD', 'WDFDR'
        Example: {'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 0.05}

    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame containing only rows passing all thresholds
    """
    return df[
        (df['SaintScore'] >= thresholds['SaintScore']) &
        (df['BFDR'] <= thresholds['BFDR']) &
        (df['WD'] >= thresholds['WD']) &
        (df['WDFDR'].fillna(1.0) <= thresholds['WDFDR'])
    ]


# Module-level cache for BioGRID data to avoid repeated file I/O
_biogrid_cache = None
_biogrid_cache_path = None


def _load_biogrid_cached(biogrid_path="/Datasets/biogrid_summary.csv"):
    """
    Load BioGRID data with caching to avoid repeated I/O.

    The BioGRID file is ~2M rows and takes 2-3 seconds to load.
    Cache it at module level so it's only loaded once per session.
    """
    global _biogrid_cache, _biogrid_cache_path

    # Return cached data if path hasn't changed
    if _biogrid_cache is not None and _biogrid_cache_path == biogrid_path:
        return _biogrid_cache

    # Load and cache the data
    try:
        _biogrid_cache = pd.read_csv(biogrid_path)
        _biogrid_cache_path = biogrid_path

        # Convert to string type for consistency (do this once)
        _biogrid_cache['SWISS-PROT Accessions Interactor A'] = _biogrid_cache['SWISS-PROT Accessions Interactor A'].astype(str)
        _biogrid_cache['SWISS-PROT Accessions Interactor B'] = _biogrid_cache['SWISS-PROT Accessions Interactor B'].astype(str)

        return _biogrid_cache
    except FileNotFoundError:
        print(f"Warning: BioGRID file not found at {biogrid_path}")
        return None


def pca_plot(interaction, experimentalDesign):

    int = pd.read_csv(interaction, sep="\t", header=0)
    int.columns = ['Experiment', 'BaitName', 'Prey', 'Intensity']
    ed = pd.read_csv(experimentalDesign, sep=",")

    # ed['shared_id'] = ed['Experiment Name'] + '_' + ed['Replicate'].astype(str)

    # Make the int table wide
    data = int.pivot(index='Prey', columns='Experiment', values='Intensity')

    metadata = int[['Experiment', 'BaitName']].drop_duplicates()

    metadata = metadata.merge(ed[['Experiment Name', 'Type']], left_on='Experiment', right_on='Experiment Name', how='left')

    ### Now make the PCA plot....

    # Clean up the data

    # Convert all 0 values to NaN
    data = data.replace(0, np.nan)

    # Remove rows with more than 75% NaN values
    data = data.dropna(thresh=len(data.columns) * 0.5)

    # Replace NaN values with the minimum of the row
    data = data.apply(lambda row: row.fillna(row.min()), axis=1)    
    
    # Now z-score normalize the data by row
    data = data.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

    # Now make a PCA plot
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(data.T)

    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])

    # Add the original column names to the PCA dataframe
    pca_df['Experiment'] = data.columns

    # Merge with the metadata to get the BaitName and Type
    pca_df = pca_df.merge(metadata, left_on='Experiment', right_on='Experiment', how='left')

    # pca_df['PC1'] = pd.to_numeric(pca_df['PC1'], errors='coerce')
    # pca_df['PC2'] = pd.to_numeric(pca_df['PC2'], errors='coerce')

    # Calculate the explained variance for each component
    explained_variance = pca.explained_variance_ratio_
    

    # Make a plotly scatterplot of the PCA results
    pc1_var = round(explained_variance[0] * 100, 2)
    pc2_var = round(explained_variance[1] * 100, 2)
    fig = px.scatter(
        pca_df,
        x='PC1',
        y='PC2',
        color='BaitName',
        symbol='Type',
        hover_name='Experiment',
        hover_data=['Experiment', 'BaitName'],
        title="PCA of Interaction Data",
        labels={
            'PC1': f'PC1 ({explained_variance[0]*100:.2f}% variance)',
            'PC2': f'PC2 ({explained_variance[1]*100:.2f}% variance)'
        }
    )

    fig.update_layout(
        legend=dict(
            orientation="h",   # horizontal legend
            yanchor="top",
            y=-0.4,            # below the plot area
            xanchor="center",
            x=0.5
        )
    )

    # fig_widget = go.FigureWidget(fig)

    return fig

def saint_known_retention(results_path, ctrl_experiments=None):

    results = pd.read_csv(results_path, sep=",")

    # If ctrl is used, filter the results
    if ctrl_experiments is not None:
        results = results[results['Experiment.ID'].isin(ctrl_experiments)]

    thresholds = np.arange(0, 1.05, 0.05)

    percents = []
    cco_means = []

    for threshold in thresholds:
        subset = results[results['SaintScore'] >= threshold]
        total = len(subset)
        knowns = len(subset[subset['In.BioGRID'] == True])
        percents.append(knowns / total if total > 0 else 0)
        cco_means.append(subset['CCO'].mean())

    fig = px.line(x=thresholds, y=percents, title="Known Retention and Cell Component Similarity by Saint Score Threshold",
                  labels={'x': 'Threshold', 'y': 'Percent Known Retention'},
                  markers=True)
    
    # Give the first trace a name and force it to show in the legend
    fig.data[0].name = "Percent Known Retention"
    fig.data[0].showlegend = True
    
    # Plot the CCO means as well
    fig.add_scatter(x=thresholds, y=cco_means, mode='lines+markers', name='CCO Mean', yaxis='y2')
    fig.update_layout(yaxis2=dict(title='CCO Mean', overlaying='y', side='right', range=[-0.05,1.05]),
                      yaxis=dict(title='Percent Known Retention', range=[-0.05,1.05]))
    fig.update_layout(
        legend=dict(
            orientation="h",   # horizontal legend
            yanchor="top",
            y=-0.4,            # below the plot area
            xanchor="center",
            x=0.5
        )
    )
    return fig

def roc_plot(results_path, known_type, ctrl_experiments=None):

    scores = pd.read_csv(results_path, sep=",")

    # If ctrl is used, filter the results
    if ctrl_experiments is not None:
        scores = scores[scores['Experiment.ID'].isin(ctrl_experiments)]

    # Get the true positives and false positives
    if known_type == 'BioGRID':
        truth = list(scores['In.BioGRID'])
    elif known_type == 'Multivalidated':
        truth = list(scores['Multivalidated'])
    else:
        raise ValueError("Unknown known_type: " + known_type)
    
    # Convert truth to boolean, where True means the interaction is known
    truth = [False if x != x else x for x in truth]

    # Create a plotly figure for the line plot ROC curves
    fig = go.Figure()

    for score in ['SaintScore', 'BFDR', 'WD', 'WDFDR']:
        if 'FDR' in score:
            fpr, tpr, thresholds = roc_curve(truth, -1 * np.array(list(scores[score])))
            thresholds = -1 * thresholds  # Adjust thresholds for FDR scores
        else:
            fpr, tpr, thresholds = roc_curve(truth, list(scores[score]))

        # For arrays (e.g., fpr, tpr, thresholds in ROC)
        mask = ~(np.isnan(fpr) | np.isnan(tpr) | np.isnan(thresholds) |
                np.isinf(fpr) | np.isinf(tpr) | np.isinf(thresholds))
        fpr = fpr[mask]
        tpr = tpr[mask]
        thresholds = thresholds[mask]

        # Prepend (0,0) and append (1,1) if not already present
        if fpr[0] != 0 or tpr[0] != 0:
            fpr = np.insert(fpr, 0, 0.0)
            tpr = np.insert(tpr, 0, 0.0)
            thresholds = np.insert(thresholds, 0, thresholds[0] if len(thresholds) > 0 else 0.0)
        if fpr[-1] != 1 or tpr[-1] != 1:
            fpr = np.append(fpr, 1.0)
            tpr = np.append(tpr, 1.0)
            thresholds = np.append(thresholds, thresholds[-1] if len(thresholds) > 0 else 1.0)

        roc_auc = auc(fpr, tpr)
        # axs.plot(fpr, tpr, label=f'{score} (AUC = {roc_auc:.2f})')
        # Plot the ROC curve for plotly
        customdata = np.stack((thresholds, tpr, fpr), axis=-1)

        fig.add_trace(go.Scatter(
            x=fpr,
            y=tpr,
            mode='lines',
            name=f'{score} (AUC = {roc_auc:.2f})',
            customdata=customdata,
             hovertemplate=(
                'Threshold: %{customdata[0]:.3f}<br>'
                'FPR: %{x:.3f}<br>'
                'TPR: %{y:.3f}<br>'
                '<extra>%{fullData.name}</extra>'
            )
        ))

    # Add the diagonal line that is not interacive
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode='lines', name='Random', line=dict(dash='dash')))

    # Add axes labels and title
    fig.update_layout(title='ROC Curves for Interaction Scores',
                      xaxis_title='False Positive Rate',
                      yaxis_title='True Positive Rate',
                      xaxis=dict(range=[0, 1]),
                      yaxis=dict(range=[0, 1]),
                      legend=dict(
                            title='Scores',
                            orientation="h",   # horizontal legend
                            yanchor="top",
                            y=-0.4,            # just above the plot
                            xanchor="center",
                            x=0.5
                        ))

    return fig


def calculate_network_degrees(passing_interactions, biogrid_path="/Datasets/biogrid_summary.csv"):
    """
    Calculate prey-prey network degree for each prey protein from BioGRID.

    Network degree = number of OTHER prey proteins (in the passing set) that
    this prey interacts with in BioGRID, excluding all bait proteins.

    Parameters:
    -----------
    passing_interactions : pd.DataFrame
        Filtered interactions with columns: First_ID (prey), Bait.ID (bait)
    biogrid_path : str
        Path to biogrid_summary.csv file

    Returns:
    --------
    list of int
        Network degrees for each unique prey (prey-prey interactions only)
    """

    # Handle edge cases
    if len(passing_interactions) == 0:
        return []

    # Load BioGRID data using cache
    biogrid = _load_biogrid_cached(biogrid_path)
    if biogrid is None:
        return [0] * len(passing_interactions)

    # Get unique prey and bait IDs from the passing interactions
    # Use set for O(1) lookup performance in filtering
    prey_ids = set(str(pid) for pid in passing_interactions['First_ID'].unique()
                   if pd.notna(pid) and str(pid) != 'nan')
    bait_ids = set(str(bid) for bid in passing_interactions['Bait.ID'].unique()
                   if pd.notna(bid) and str(bid) != 'nan')

    if len(prey_ids) == 0:
        return []

    # Filter BioGRID for prey-prey interactions only
    # Keep only edges where BOTH proteins are in prey set and NEITHER is in bait set
    col_a = biogrid['SWISS-PROT Accessions Interactor A']
    col_b = biogrid['SWISS-PROT Accessions Interactor B']

    prey_prey_edges = biogrid[
        col_a.isin(prey_ids) & col_b.isin(prey_ids) &
        ~col_a.isin(bait_ids) & ~col_b.isin(bait_ids)
    ]

    # If no prey-prey edges found, return zeros
    if len(prey_prey_edges) == 0:
        return [0] * len(prey_ids)

    # Count degree efficiently using value_counts on concatenated series
    # This is faster than creating separate DataFrames and concatenating
    all_proteins = pd.concat([
        prey_prey_edges['SWISS-PROT Accessions Interactor A'],
        prey_prey_edges['SWISS-PROT Accessions Interactor B']
    ])
    degree_counts = all_proteins.value_counts().to_dict()

    # Calculate degrees for each unique prey
    degrees = [degree_counts.get(prey_id, 0) for prey_id in prey_ids]

    return degrees


def calculate_threshold_metrics(results_path, thresholds, ctrl_experiments=None):
    """
    Calculate metrics for interactions passing thresholds.

    Parameters:
    -----------
    results_path : str
        Path to annotated_scores.csv file
    thresholds : dict
        Dictionary with threshold values: {'SaintScore': float, 'BFDR': float,
                                           'WD': float, 'WDFDR': float}
    ctrl_experiments : list, optional
        List of control experiment IDs to filter by

    Returns:
    --------
    dict
        Dictionary containing:
        - median_network_size: Median number of interactions per bait after filtering
        - enrichment_ratio: Average enrichment of known interactions
        - mean_degree: Mean prey-prey network degree (average number of other
                      passing prey proteins each prey interacts with in BioGRID)
        - total_before: Total interactions before filtering
        - total_after: Total interactions after filtering
    """

    # Load data
    results = pd.read_csv(results_path, sep=",")

    # Filter by control experiments if specified
    if ctrl_experiments is not None:
        results = results[results['Experiment.ID'].isin(ctrl_experiments)]

    # Total interactions before filtering
    total_before = len(results)

    # Known interactions before filtering
    known_before = results['In.BioGRID'].sum() if 'In.BioGRID' in results.columns else 0

    # Combined threshold (all conditions must pass using AND logic)
    passing_all = apply_score_thresholds(results, thresholds)

    # Total interactions after filtering
    total_after = len(passing_all)

    # Known interactions after filtering
    known_after = passing_all['In.BioGRID'].sum() if 'In.BioGRID' in passing_all.columns else 0

    # Calculate median network size (interactions per bait)
    if total_after > 0:
        network_sizes = passing_all.groupby('Experiment.ID').size()
        median_network_size = network_sizes.median()
    else:
        median_network_size = 0

    # Calculate average enrichment of known interactions
    # Enrichment = (known_after / total_after) / (known_before / total_before)
    if total_before > 0 and total_after > 0 and known_before > 0:
        pct_before = known_before / total_before
        pct_after = known_after / total_after
        enrichment_ratio = pct_after / pct_before
    else:
        enrichment_ratio = 0

    # Calculate mean prey-prey network degree from BioGRID
    if total_after > 0:
        degrees = calculate_network_degrees(passing_all)
        mean_degree = np.mean(degrees) if len(degrees) > 0 else 0
    else:
        mean_degree = 0

    return {
        'median_network_size': median_network_size,
        'enrichment_ratio': enrichment_ratio,
        'mean_degree': mean_degree,
        'total_before': total_before,
        'total_after': total_after,
        'known_before': known_before,
        'known_after': known_after
    }


def saint_scatter_plot(results_path, bait_name, saintscore_threshold):
    """
    Create scatter plot of SAINT Score vs Fold Change for a specific bait.

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
    plotly.graph_objects.Figure
        Interactive scatter plot
    """

    # Load data
    results = pd.read_csv(results_path, sep=",")

    # Filter for specific bait
    bait_data = results[results['Experiment.ID'] == bait_name].copy()

    if len(bait_data) == 0:
        # Return empty figure with message
        fig = go.Figure()
        fig.add_annotation(
            text=f"No data available for bait: {bait_name}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16, color="red")
        )
        return fig

    # Helper function to calculate average control intensity
    def calculate_avg_ctrl_intensity(ctrl_intensity_str):
        """
        Parse control intensity string, filter out missing values ('.'),
        and calculate average.
        """
        if pd.isnull(ctrl_intensity_str):
            return np.nan

        # Split by '|' and filter out '.' values
        values = [v.strip() for v in str(ctrl_intensity_str).split('|') if v.strip() != '.']

        if len(values) == 0:
            return np.nan

        try:
            # Convert to float and calculate mean
            numeric_values = [float(v) for v in values]
            return np.mean(numeric_values)
        except (ValueError, TypeError):
            return np.nan

    # Separate data by BioGRID status
    # Handle missing columns gracefully
    has_biogrid = 'In.BioGRID' in bait_data.columns
    has_multivalidated = 'Multivalidated' in bait_data.columns

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

    # Create scatter plot
    fig = go.Figure()

    # Add not in BioGRID (blue) - plot first so it's in the background
    if len(not_in_biogrid) > 0:
        hover_text = []
        for idx, row in not_in_biogrid.iterrows():
            avg_ctrl = calculate_avg_ctrl_intensity(row['ctrlIntensity'])
            if np.isnan(avg_ctrl):
                ctrl_text = "NaN"
            else:
                ctrl_text = f"{avg_ctrl:.2e}"

            text = (
                f"<b>{row['First_Prey_Gene']}</b><br>"
                f"SAINT Score: {row['SaintScore']:.3f}<br>"
                f"BFDR: {row['BFDR']:.3f}<br>"
                f"Fold Change: {row['FoldChange']:.3f}<br>"
                f"Avg Intensity: {row['AvgIntensity']:.2e}<br>"
                f"Avg Ctrl Intensity: {ctrl_text}"
            )
            hover_text.append(text)

        fig.add_trace(go.Scatter(
            x=not_in_biogrid['FoldChange'],
            y=not_in_biogrid['SaintScore'],
            mode='markers',
            marker=dict(
                size=8,
                color='#1f77b4',  # Blue
                line=dict(width=0.5, color='white')
            ),
            text=hover_text,
            hovertemplate='%{text}<extra></extra>',
            name='Not in BioGRID'
        ))

    # Add in BioGRID (orange)
    if len(in_biogrid) > 0:
        hover_text = []
        for idx, row in in_biogrid.iterrows():
            avg_ctrl = calculate_avg_ctrl_intensity(row['ctrlIntensity'])
            if np.isnan(avg_ctrl):
                ctrl_text = "NaN"
            else:
                ctrl_text = f"{avg_ctrl:.2e}"

            text = (
                f"<b>{row['First_Prey_Gene']}</b><br>"
                f"SAINT Score: {row['SaintScore']:.3f}<br>"
                f"BFDR: {row['BFDR']:.3f}<br>"
                f"Fold Change: {row['FoldChange']:.3f}<br>"
                f"Avg Intensity: {row['AvgIntensity']:.2e}<br>"
                f"Avg Ctrl Intensity: {ctrl_text}"
            )
            hover_text.append(text)

        fig.add_trace(go.Scatter(
            x=in_biogrid['FoldChange'],
            y=in_biogrid['SaintScore'],
            mode='markers',
            marker=dict(
                size=8,
                color='#ff7f0e',  # Orange
                line=dict(width=0.5, color='white')
            ),
            text=hover_text,
            hovertemplate='%{text}<extra></extra>',
            name='In BioGRID'
        ))

    # Add multivalidated (red) - plot last so it's on top
    if len(multivalidated) > 0:
        hover_text = []
        for idx, row in multivalidated.iterrows():
            avg_ctrl = calculate_avg_ctrl_intensity(row['ctrlIntensity'])
            if np.isnan(avg_ctrl):
                ctrl_text = "NaN"
            else:
                ctrl_text = f"{avg_ctrl:.2e}"

            text = (
                f"<b>{row['First_Prey_Gene']}</b><br>"
                f"SAINT Score: {row['SaintScore']:.3f}<br>"
                f"BFDR: {row['BFDR']:.3f}<br>"
                f"Fold Change: {row['FoldChange']:.3f}<br>"
                f"Avg Intensity: {row['AvgIntensity']:.2e}<br>"
                f"Avg Ctrl Intensity: {ctrl_text}"
            )
            hover_text.append(text)

        fig.add_trace(go.Scatter(
            x=multivalidated['FoldChange'],
            y=multivalidated['SaintScore'],
            mode='markers',
            marker=dict(
                size=8,
                color='#d62728',  # Red
                line=dict(width=0.5, color='white')
            ),
            text=hover_text,
            hovertemplate='%{text}<extra></extra>',
            name='Multivalidated'
        ))

    # Add horizontal threshold line
    fig.add_hline(
        y=saintscore_threshold,
        line_dash="dot",
        line_color="red",
        line_width=2,
        annotation_text=f"SAINT Score Threshold: {saintscore_threshold}",
        annotation_position="right"
    )

    # Update layout
    fig.update_layout(
        title=f"SAINT Score vs Fold Change - {bait_name}",
        xaxis_title="Fold Change (log2)",
        yaxis_title="SAINT Score",
        xaxis=dict(zeroline=True, zerolinewidth=1, zerolinecolor='lightgray'),
        yaxis=dict(range=[-0.05, 1.05]),
        hovermode='closest',
        height=500,
        template='plotly_white',
        legend=dict(
            yanchor="bottom",
            y=0.01,
            xanchor="right",
            x=0.99
        )
    )

    return fig


if __name__ == "__main__":
    pass