import os

import plotly.express as px
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.metrics import roc_curve, auc
from log_config import get_logger
import provenance

logger = get_logger(__name__)


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


def _load_biogrid_cached(biogrid_path):
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
        logger.warning("BioGRID file not found at %s; known-interaction plots "
                       "will have no reference set.", biogrid_path)
        return None


def read_interactions(interaction):
    """The SAINT interaction file as a frame of Experiment, BaitName, Prey, Intensity."""
    df = pd.read_csv(interaction, sep="\t", header=0)
    df.columns = ['Experiment', 'BaitName', 'Prey', 'Intensity']
    return df


def interaction_matrix(interaction):
    """Prey x experiment intensity matrix from a SAINT interaction file; a prey not
    listed under an experiment is NaN."""
    return read_interactions(interaction).pivot(index='Prey', columns='Experiment',
                                                values='Intensity')


def prepare_pca_matrix(interaction, min_detection_frac=0.5,
                       imputation="row_min", normalization="zscore"):
    """Load a SAINT interaction file and return the prey x experiment matrix
    both PCA plots run on.

    Zeros are treated as non-detections. Preys detected in fewer than
    ``min_detection_frac`` of experiments are removed before imputation.

    imputation: 'row_min' (fill with the prey's minimum observed value),
    'zero' (fill with 0), or 'drop' (keep only fully detected preys).

    normalization: 'zscore' (per-prey), 'log2_zscore' (log2(x + 1) then
    per-prey z-score; the pseudocount keeps zero-imputed values finite),
    or 'none'. Preys with zero variance are dropped by the z-score options.

    A prey still missing a value after imputation (one never detected, which 'row_min'
    has nothing to fill from) is dropped whatever the normalization: PCA cannot take NaN.
    """
    data = interaction_matrix(interaction)

    if len(data.columns) < 2:
        raise ValueError(
            f"Only {len(data.columns)} experiments in the interaction data; "
            "PCA needs at least 2.")

    data = data.replace(0, np.nan)
    data = data.dropna(thresh=len(data.columns) * min_detection_frac)

    if imputation == "row_min":
        data = data.apply(lambda row: row.fillna(row.min()), axis=1)
    elif imputation == "zero":
        data = data.fillna(0)
    elif imputation == "drop":
        data = data.dropna()
    else:
        raise ValueError(f"Unknown imputation option: {imputation!r}")
    data = data.dropna(how='any')

    if normalization == "log2_zscore":
        data = np.log2(data + 1)
    if normalization in ("zscore", "log2_zscore"):
        data = data.apply(lambda row: (row - row.mean()) / row.std(), axis=1)
        data = data.dropna(how='any')
    elif normalization != "none":
        raise ValueError(f"Unknown normalization option: {normalization!r}")

    if len(data) < 3:
        raise ValueError(
            f"Only {len(data)} preys remain after filtering "
            f"(min_detection_frac={min_detection_frac}, imputation={imputation!r}); "
            "PCA needs at least 3. Relax the detection or imputation settings.")
    return data


def detection_counts(interaction):
    """Per-prey count of experiments with a nonzero, non-missing intensity."""
    return interaction_matrix(interaction).replace(0, np.nan).notna().sum(axis=1)


def reduce_categorical(series, top_n=12, missing_label="Unknown"):
    """Reduce a possibly multi-valued annotation column to plottable categories.

    Multi-valued entries (semicolon-joined, e.g. HPA 'Main location') are
    reduced to their first value; the top_n most frequent labels are kept and
    the rest collapsed to 'Other'; missing values become ``missing_label``.
    """
    s = series.astype("string").str.split(";").str[0].str.strip()
    top = s.value_counts().head(top_n).index
    s = s.where(s.isin(top) | s.isna(), other="Other")
    return s.fillna(missing_label).astype(str)


def load_pca_metadata(interaction, experimentalDesign):
    """Experiment-level metadata (BaitName from the interaction file, Type from
    the ED file) for labeling PCA points."""
    df = read_interactions(interaction)
    ed = pd.read_csv(experimentalDesign, sep=",")
    metadata = df[['Experiment', 'BaitName']].drop_duplicates()
    return metadata.merge(ed[['Experiment Name', 'Type']],
                          left_on='Experiment', right_on='Experiment Name',
                          how='left')


def experiment_pca(interaction, experimentalDesign, matrix=None):
    """Experiment-level PCA coordinates, one row per experiment.

    Returns ``(frame, explained_variance_ratio)``; the frame carries PC1, PC2,
    Experiment, BaitName and Type.  ``matrix``, if given, is a prepare_pca_matrix
    result; otherwise it is computed with default settings.
    """
    if matrix is None:
        matrix = prepare_pca_matrix(interaction)
    metadata = load_pca_metadata(interaction, experimentalDesign)

    pca = PCA(n_components=2)
    pca_df = pd.DataFrame(data=pca.fit_transform(matrix.T), columns=['PC1', 'PC2'])
    pca_df['Experiment'] = matrix.columns
    pca_df = pca_df.merge(metadata, on='Experiment', how='left')
    return pca_df, pca.explained_variance_ratio_


def prey_pca(matrix):
    """Prey-level PCA coordinates: preys as samples, experiments as features -- an
    independent embedding, not the loadings of the experiment PCA.

    Returns ``(frame, explained_variance_ratio)``; the frame carries PC1, PC2, Prey.
    """
    pca = PCA(n_components=2)
    prey_df = pd.DataFrame(data=pca.fit_transform(matrix), columns=['PC1', 'PC2'])
    prey_df['Prey'] = matrix.index
    return prey_df, pca.explained_variance_ratio_


def pca_plot(interaction, experimentalDesign, matrix=None):
    """Experiment-level PCA (one point per experiment); see ``experiment_pca``."""
    pca_df, explained_variance = experiment_pca(interaction, experimentalDesign, matrix)

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

    return fig


def prey_gene_names(prey_file):
    """Gene name per prey accession from a SAINT prey.txt (no header; accession
    first, gene name last -- the spectral-count layout has sequence length between)."""
    prey = pd.read_csv(prey_file, sep="\t", header=None, dtype=str)
    return pd.Series(prey.iloc[:, -1].values, index=prey.iloc[:, 0].values)


def prey_pca_plot(matrix, color_values=None, color_label=None,
                  color_mode="none", color_threshold=None, gene_names=None):
    """Prey-level PCA figure; see ``prey_pca`` for the embedding.

    color_values: pd.Series indexed by prey (numeric for 'continuous',
    labels for 'categorical'), or None with color_mode 'none'.

    color_threshold: for 'continuous' only -- preys with a value below it are
    drawn grey so the color scale is spent on the informative range.

    gene_names: pd.Series indexed by prey; when given, the hover shows the gene
    name in bold above the accession.
    """
    prey_df, explained_variance = prey_pca(matrix)
    if gene_names is None:
        prey_df['Gene'] = prey_df['Prey']
    else:
        prey_df['Gene'] = prey_df['Prey'].map(gene_names).fillna(prey_df['Prey'])

    labels = {
        'PC1': f'PC1 ({explained_variance[0]*100:.2f}% variance)',
        'PC2': f'PC2 ({explained_variance[1]*100:.2f}% variance)'
    }
    kwargs = dict(x='PC1', y='PC2', hover_name='Gene', hover_data={'Prey': True},
                  title="Prey PCA", labels=labels)

    if color_mode == "continuous":
        prey_df[color_label] = prey_df['Prey'].map(color_values)
        if color_threshold is not None:
            below = prey_df[prey_df[color_label] < color_threshold]
            above = prey_df[~(prey_df[color_label] < color_threshold)]
            fig = px.scatter(above, color=color_label,
                             color_continuous_scale='Viridis', **kwargs)
            fig.add_trace(go.Scatter(
                x=below['PC1'], y=below['PC2'], mode='markers',
                name=f"{color_label} < {color_threshold:g}",
                marker=dict(color='lightgrey'),
                text=below['Gene'], customdata=below[['Prey']],
                hovertemplate="<b>%{text}</b><br>%{customdata[0]}<br>PC1=%{x}<br>PC2=%{y}<extra></extra>",
            ))
            # grey first so scoring preys draw on top of it
            fig.data = fig.data[-1:] + fig.data[:-1]
            fig.update_layout(showlegend=True)
        else:
            fig = px.scatter(prey_df, color=color_label,
                             color_continuous_scale='Viridis', **kwargs)
    elif color_mode == "categorical":
        prey_df[color_label] = prey_df['Prey'].map(color_values)
        categories = [c for c in prey_df[color_label].dropna().unique()
                      if c not in ("Other", "Unknown")]
        categories = sorted(categories) + ["Other", "Unknown"]
        fig = px.scatter(prey_df, color=color_label,
                         category_orders={color_label: categories}, **kwargs)
    else:
        fig = px.scatter(prey_df, **kwargs)

    # Plotly Express lists hover_data after the axes; move the accession up so it
    # sits directly under the bold gene name.
    for trace in fig.data:
        trace.hovertemplate = (trace.hovertemplate
                               .replace("<br>Prey=%{customdata[0]}", "")
                               .replace("<br><br>", "<br>%{customdata[0]}<br>"))
    fig.update_traces(marker=dict(size=5, opacity=0.7))
    fig.update_layout(
        legend=dict(
            orientation="h",
            yanchor="top",
            y=-0.4,
            xanchor="center",
            x=0.5
        )
    )
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


def calculate_network_degrees(passing_interactions, biogrid_path):
    """
    Calculate prey-prey network degree for each prey protein from BioGRID.

    Network degree = number of OTHER prey proteins (in the passing set) that
    this prey interacts with in BioGRID, excluding all bait proteins.

    Parameters:
    -----------
    passing_interactions : pd.DataFrame
        Filtered interactions with columns: First_ID (prey), Bait.ID (bait)
    biogrid_path : str
        Path to the biogrid_summary.csv built for the relevant organism

    Returns:
    --------
    list of int, or None
        Network degrees for each unique prey (prey-prey interactions only).  None when
        the reference set could not be read: a degree of zero is a real result for a
        prey with no published partners, so absent data must not be reported as one.
    """

    # Handle edge cases
    if len(passing_interactions) == 0:
        return []

    # Load BioGRID data using cache
    biogrid = _load_biogrid_cached(biogrid_path)
    if biogrid is None:
        return None

    # Get unique prey and bait IDs from the passing interactions
    # Use set for O(1) lookup performance in filtering
    prey_ids = set(str(pid) for pid in passing_interactions['First_ID'].unique()
                   if pd.notna(pid) and str(pid) != 'nan')
    # Baits are excluded by accession, the form BioGRID uses; results annotated before
    # symbols were resolved carry only the supplied Bait.ID.
    bait_col = 'Bait_Accession' if 'Bait_Accession' in passing_interactions.columns else 'Bait.ID'
    bait_ids = set(str(bid) for bid in passing_interactions[bait_col].unique()
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


def calculate_threshold_metrics(results_path, thresholds, ctrl_experiments=None,
                                biogrid_path=None):
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
    biogrid_path : str, optional
        Reference set for the network degree.  Defaults to the summary built for the
        organism the dataset was annotated against.

    Returns:
    --------
    dict
        Dictionary containing:
        - median_network_size: Median number of interactions per bait after filtering
        - enrichment_ratio: Average enrichment of known interactions
        - mean_degree: Mean prey-prey network degree (average number of other
                      passing prey proteins each prey interacts with in BioGRID),
                      or None when the organism's BioGRID summary is unavailable
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
        if biogrid_path is None:
            biogrid_path = provenance.biogrid_summary_path(
                **provenance.annotation_settings(os.path.dirname(results_path)))
        degrees = calculate_network_degrees(passing_all, biogrid_path)
        if degrees is None:
            mean_degree = None
        else:
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


def bait_scores(results_path, bait_name):
    """The annotated-score rows of one bait."""
    results = pd.read_csv(results_path, sep=",")
    return results[results['Experiment.ID'] == bait_name].copy()


def known_status_split(bait_data):
    """Split a bait's rows into ``(not_in_biogrid, in_biogrid, multivalidated)``.

    Multivalidated rows are left out of the BioGRID group so each row lands in one
    group.  Without a Multivalidated column that group is empty; without an In.BioGRID
    column every row is "not in BioGRID".
    """
    if 'Multivalidated' in bait_data.columns:
        multivalidated = bait_data[bait_data['Multivalidated'] == True].copy()
        in_biogrid = bait_data[(bait_data['In.BioGRID'] == True)
                               & (bait_data['Multivalidated'] != True)].copy()
        not_in_biogrid = bait_data[bait_data['In.BioGRID'] != True].copy()
    elif 'In.BioGRID' in bait_data.columns:
        multivalidated = pd.DataFrame()
        in_biogrid = bait_data[bait_data['In.BioGRID'] == True].copy()
        not_in_biogrid = bait_data[bait_data['In.BioGRID'] != True].copy()
    else:
        multivalidated = pd.DataFrame()
        in_biogrid = pd.DataFrame()
        not_in_biogrid = bait_data.copy()
    return not_in_biogrid, in_biogrid, multivalidated


# Legend name and colour of each known-status group, in draw order: the common group
# goes down first so the rarer groups sit on top.
KNOWN_STATUS_STYLE = (('Not in BioGRID', '#1f77b4'),
                      ('In BioGRID', '#ff7f0e'),
                      ('Multivalidated', '#d62728'))


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

    bait_data = bait_scores(results_path, bait_name)

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

    # SAINTexpress names the quantity columns by input type: its intensity build writes
    # AvgIntensity/ctrlIntensity, its spectral-count build AvgSpec/ctrlCounts.
    if 'AvgIntensity' in bait_data.columns:
        avg_col, ctrl_col, label, fmt = 'AvgIntensity', 'ctrlIntensity', 'Intensity', '.2e'
    elif 'AvgSpec' in bait_data.columns:
        avg_col, ctrl_col, label, fmt = 'AvgSpec', 'ctrlCounts', 'Spec', '.1f'
    else:
        raise KeyError("annotated scores carry neither AvgIntensity nor AvgSpec")

    def calculate_avg_ctrl(ctrl_str):
        """Mean of SAINT's '|'-separated control values, ignoring '.' placeholders."""
        if pd.isnull(ctrl_str):
            return np.nan
        values = [v.strip() for v in str(ctrl_str).split('|') if v.strip() != '.']
        if len(values) == 0:
            return np.nan
        try:
            return np.mean([float(v) for v in values])
        except (ValueError, TypeError):
            return np.nan

    def hover_text(frame):
        texts = []
        for _, row in frame.iterrows():
            avg_ctrl = calculate_avg_ctrl(row[ctrl_col])
            ctrl_text = "NaN" if np.isnan(avg_ctrl) else f"{avg_ctrl:{fmt}}"
            texts.append(
                f"<b>{row['First_Prey_Gene']}</b><br>"
                f"SAINT Score: {row['SaintScore']:.3f}<br>"
                f"BFDR: {row['BFDR']:.3f}<br>"
                f"Fold Change: {row['FoldChange']:.3f}<br>"
                f"Avg {label}: {row[avg_col]:{fmt}}<br>"
                f"Avg Ctrl {label}: {ctrl_text}"
            )
        return texts

    fig = go.Figure()
    for group, (name, color) in zip(known_status_split(bait_data), KNOWN_STATUS_STYLE):
        if len(group) == 0:
            continue
        fig.add_trace(go.Scatter(
            x=group['FoldChange'],
            y=group['SaintScore'],
            mode='markers',
            marker=dict(size=8, color=color, line=dict(width=0.5, color='white')),
            text=hover_text(group),
            hovertemplate='%{text}<extra></extra>',
            name=name
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