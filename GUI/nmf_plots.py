"""
NMF Clustering Visualization Module

This module provides Non-negative Matrix Factorization (NMF) analysis with t-SNE
visualization for proximity labeling proteomics data. It includes:
- Data loading and pre-filtering by BFDR threshold
- NMF decomposition for dimensionality reduction
- t-SNE visualization of NMF components
- Automatic cluster detection and labeling

Author: ProxiMate GUI Integration
Date: 2025
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from sklearn.decomposition import NMF
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import warnings


def construct_matrices(csv_path, bfdr_threshold):
    """
    Load annotated_scores.csv, pre-filter by BFDR, and create prey × bait matrices.

    Parameters:
    -----------
    csv_path : str
        Path to annotated_scores.csv file
    bfdr_threshold : float
        BFDR threshold for pre-filtering (0.0-1.0)

    Returns:
    --------
    dict containing:
        - V_saint: np.ndarray, prey × bait matrix of SAINT scores
        - V_fc_scaled: np.ndarray, prey × bait matrix of scaled fold changes
        - preys: list of str, prey gene names
        - baits: list of str, bait experiment IDs
        - error: str or None, error message if data insufficient
    """
    try:
        # Load the data
        data = pd.read_csv(csv_path, sep=",")

        # Handle column name variations
        if 'First_Prey_Gene' in data.columns:
            prey_col = 'First_Prey_Gene'
        elif 'PreyGene' in data.columns:
            prey_col = 'PreyGene'
        else:
            return {'error': 'Required columns missing: First_Prey_Gene or PreyGene'}

        # Check for other required columns
        required_cols = ['Experiment.ID', 'SaintScore', 'BFDR', 'FoldChange']
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            return {'error': f'Required columns missing: {", ".join(missing_cols)}'}

        # Filter by BFDR threshold
        data_filtered = data[data['BFDR'] <= bfdr_threshold].copy()

        if len(data_filtered) == 0:
            return {'error': 'Insufficient data after filtering. Try lowering BFDR threshold.'}

        # Get unique preys and baits
        preys_raw = data_filtered[prey_col].unique()
        baits = data_filtered['Experiment.ID'].unique()

        # Simplify prey names (take first identifier before semicolon)
        preys = [str(prey).split(';')[0] for prey in preys_raw]

        # Check for minimum samples
        if len(preys) < 5:
            return {'error': f'Insufficient data after filtering. Only {len(preys)} preys found. Try lowering BFDR threshold.'}
        if len(baits) < 2:
            return {'error': f'Insufficient data after filtering. Only {len(baits)} baits found.'}

        # Add simplified prey column
        data_filtered['Prey_Simple'] = data_filtered[prey_col].str.split(';').str[0]

        # Create SAINT score matrix
        saint_pivot = data_filtered.pivot_table(
            index='Prey_Simple',
            columns='Experiment.ID',
            values='SaintScore',
            aggfunc='first'  # Take first value if duplicates
        )

        # Create fold-change matrix
        fc_pivot = data_filtered.pivot_table(
            index='Prey_Simple',
            columns='Experiment.ID',
            values='FoldChange',
            aggfunc='first'
        )

        # Ensure both matrices have the same prey and bait order
        preys = list(saint_pivot.index)
        baits = list(saint_pivot.columns)

        # Convert to numpy arrays
        V_saint = saint_pivot.values
        V_fc = fc_pivot.values

        # Handle missing values
        # SAINT scores: NaN → 0 (no interaction)
        V_saint = np.nan_to_num(V_saint, nan=0.0)

        # Fold changes: log2 transform and scale
        # First, handle zeros and NaN
        V_fc = np.where(V_fc <= 0, np.nan, V_fc)  # Negative/zero FC → NaN
        V_fc = np.log2(V_fc)  # log2 transform

        # Replace NaN with row median (or 0 if all NaN in row)
        for i in range(len(V_fc)):
            row = V_fc[i]
            row_median = np.nanmedian(row)
            if np.isnan(row_median):
                row_median = 0
            V_fc[i] = np.where(np.isnan(row), row_median, row)

        # Scale each row to [0, 1]
        V_fc_scaled = np.zeros_like(V_fc)
        for i in range(len(V_fc)):
            row = V_fc[i]
            row_min = row.min()
            row_max = row.max()
            if row_max > row_min:
                V_fc_scaled[i] = (row - row_min) / (row_max - row_min)
            else:
                V_fc_scaled[i] = 0.5  # If all same value, set to middle

        return {
            'V_saint': V_saint,
            'V_fc_scaled': V_fc_scaled,
            'preys': preys,
            'baits': baits,
            'error': None
        }

    except FileNotFoundError:
        return {'error': 'Annotated scores file not found'}
    except Exception as e:
        return {'error': f'Error loading data: {str(e)}'}


def perform_nmf(V_matrix, n_components, random_state=42):
    """
    Run NMF decomposition on a matrix.

    Parameters:
    -----------
    V_matrix : np.ndarray
        Input matrix (preys × baits)
    n_components : int
        Number of components for NMF
    random_state : int
        Random seed for reproducibility

    Returns:
    --------
    dict containing:
        - W: np.ndarray, prey loadings (preys × components)
        - H: np.ndarray, bait loadings (components × baits)
        - reconstruction_error: float
        - error: str or None
    """
    try:
        # Suppress convergence warnings
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UserWarning)

            # Initialize NMF model
            model = NMF(
                n_components=n_components,
                init='nndsvda',  # Non-negative double SVD initialization
                solver='cd',  # Coordinate descent solver
                max_iter=500,
                random_state=random_state,
                alpha_W=0.0,
                alpha_H=0.0,
                l1_ratio=0.0
            )

            # Fit and transform
            W = model.fit_transform(V_matrix)
            H = model.components_

            return {
                'W': W,
                'H': H,
                'reconstruction_error': model.reconstruction_err_,
                'error': None
            }

    except Exception as e:
        return {'error': f'NMF failed to converge. Try different parameters.'}


def perform_tsne(matrix, random_state=42):
    """
    Run t-SNE dimensionality reduction.

    Parameters:
    -----------
    matrix : np.ndarray
        NMF loading matrix (samples × components)
    random_state : int
        Random seed for reproducibility

    Returns:
    --------
    np.ndarray: t-SNE coordinates (samples × 2)
    """
    n_samples = len(matrix)

    # Calculate adaptive perplexity
    # Perplexity should be smaller than n_samples
    perplexity = min(30, max(5, n_samples // 4))

    # Initialize and run t-SNE
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        random_state=random_state,
        init='random',
        n_iter=1000
    )

    tsne_coords = tsne.fit_transform(matrix)
    return tsne_coords


def assign_simple_clusters(W_matrix):
    """
    Assign cluster labels using KMeans with silhouette-based selection.

    Parameters:
    -----------
    W_matrix : np.ndarray
        NMF W matrix (samples × components)

    Returns:
    --------
    np.ndarray: cluster label strings (e.g., ["Cluster 1", "Cluster 2", ...])
    """
    n_samples = len(W_matrix)

    # Determine range of clusters to test
    min_clusters = 2
    max_clusters = min(10, n_samples // 2)

    if max_clusters < min_clusters:
        # Too few samples for clustering, assign all to one cluster
        return np.array([f"Cluster 1"] * n_samples)

    # Test different numbers of clusters and find best silhouette score
    best_score = -1
    best_n_clusters = min_clusters
    best_labels = None

    for n_clusters in range(min_clusters, max_clusters + 1):
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        labels = kmeans.fit_predict(W_matrix)

        # Calculate silhouette score
        score = silhouette_score(W_matrix, labels)

        if score > best_score:
            best_score = score
            best_n_clusters = n_clusters
            best_labels = labels

    # Convert numeric labels to string labels
    cluster_labels = np.array([f"Cluster {label + 1}" for label in best_labels])

    return cluster_labels


def create_tsne_plot(tsne_coords, labels, cluster_labels, title):
    """
    Create interactive Plotly scatter plot of t-SNE results.

    Parameters:
    -----------
    tsne_coords : np.ndarray
        t-SNE coordinates (N × 2)
    labels : list of str
        Sample names (prey or bait names)
    cluster_labels : np.ndarray
        Cluster assignments (strings)
    title : str
        Plot title

    Returns:
    --------
    plotly.graph_objects.Figure
    """
    # Create DataFrame for plotting
    plot_df = pd.DataFrame({
        'x': tsne_coords[:, 0],
        'y': tsne_coords[:, 1],
        'label': labels,
        'cluster': cluster_labels
    })

    # Sort by cluster for consistent legend ordering
    plot_df = plot_df.sort_values('cluster')

    # Create scatter plot using plotly express
    fig = px.scatter(
        plot_df,
        x='x',
        y='y',
        color='cluster',
        hover_data={'label': True, 'cluster': True, 'x': False, 'y': False},
        labels={'x': 't-SNE 1', 'y': 't-SNE 2', 'cluster': 'Cluster', 'label': 'Name'},
        title=title
    )

    # Update layout
    fig.update_layout(
        template='plotly_white',
        height=400,
        xaxis_title='t-SNE 1',
        yaxis_title='t-SNE 2',
        legend=dict(
            yanchor="bottom",
            y=0.01,
            xanchor="right",
            x=0.99,
            bgcolor='rgba(255,255,255,0.8)'
        ),
        hovermode='closest'
    )

    # Update traces
    fig.update_traces(
        marker=dict(size=8, line=dict(width=0.5, color='white')),
        hovertemplate='<b>%{customdata[0]}</b><br>Cluster: %{customdata[1]}<extra></extra>'
    )

    return fig


def run_nmf_pipeline(csv_path, bfdr_threshold, n_components):
    """
    Complete NMF pipeline: load data, run NMF, run t-SNE, create plots.

    This is the main entry point called by the GUI.

    Parameters:
    -----------
    csv_path : str
        Path to annotated_scores.csv
    bfdr_threshold : float
        Pre-filtering threshold (0.0-1.0)
    n_components : int
        Number of NMF components

    Returns:
    --------
    dict containing:
        - prey_fig: plotly.graph_objects.Figure or None
        - bait_fig: plotly.graph_objects.Figure or None
        - error: str or None
        - n_preys: int
        - n_baits: int
    """
    # Step 1: Load and construct matrices
    matrix_result = construct_matrices(csv_path, bfdr_threshold)

    if matrix_result.get('error'):
        return {
            'prey_fig': None,
            'bait_fig': None,
            'error': matrix_result['error'],
            'n_preys': 0,
            'n_baits': 0
        }

    V_saint = matrix_result['V_saint']
    V_fc_scaled = matrix_result['V_fc_scaled']
    preys = matrix_result['preys']
    baits = matrix_result['baits']

    n_preys = len(preys)
    n_baits = len(baits)

    # Step 2: Validate n_components
    max_components = min(n_preys, n_baits)
    if n_components > max_components:
        return {
            'prey_fig': None,
            'bait_fig': None,
            'error': f'Too few samples ({n_preys} preys, {n_baits} baits) for {n_components} components. Maximum is {max_components}.',
            'n_preys': n_preys,
            'n_baits': n_baits
        }

    # Step 3: Run NMF on SAINT score matrix
    nmf_result = perform_nmf(V_saint, n_components)

    if nmf_result.get('error'):
        return {
            'prey_fig': None,
            'bait_fig': None,
            'error': nmf_result['error'],
            'n_preys': n_preys,
            'n_baits': n_baits
        }

    W = nmf_result['W']
    H = nmf_result['H']

    # Step 4: Run t-SNE on prey loadings (W matrix)
    prey_tsne_coords = perform_tsne(W)

    # Step 5: Run t-SNE on bait loadings (H transposed)
    bait_tsne_coords = perform_tsne(H.T)

    # Step 6: Assign clusters
    prey_clusters = assign_simple_clusters(W)
    bait_clusters = assign_simple_clusters(H.T)

    # Step 7: Create plots
    prey_fig = create_tsne_plot(
        prey_tsne_coords,
        preys,
        prey_clusters,
        "Prey Clustering"
    )

    bait_fig = create_tsne_plot(
        bait_tsne_coords,
        baits,
        bait_clusters,
        "Bait Clustering"
    )

    return {
        'prey_fig': prey_fig,
        'bait_fig': bait_fig,
        'error': None,
        'n_preys': n_preys,
        'n_baits': n_baits
    }
