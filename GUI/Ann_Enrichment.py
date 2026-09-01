import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import argparse
from collections import Counter
from QC_plots import apply_score_thresholds
import matplotlib.ticker as ticker

# import dash_bio


def enrich_foreground(foreground_ids, all_ids, feature_map):
    """
    foreground_ids: set of proteins deemed 'hits'
    all_ids:      set of all proteins seen in the assay
    feature_map:  dict protein_id -> list of features
    """
    M = len(all_ids)
    n = len(foreground_ids)

    # count K (# proteins in all with each feature) and k (# in foreground with each feature)
    feat_population = Counter(feat
                     for prot, feats in feature_map.items()
                     if prot in all_ids
                     for feat in feats)
    feat_foreground = Counter(feat
                     for prot in foreground_ids
                     for feat in feature_map.get(prot, []))

    rows = []
    for feat, k in feat_foreground.items():
        K = feat_population[feat]
        # skip super‐rare features
        if K < 5 or k < 2:  
            continue
        p = hypergeom.sf(k-1, M, K, n)
        enrich = (k/n) / (K/M)
        rows.append((feat, k, n, K, M, p, enrich))

    df = (pd.DataFrame(rows, 
                       columns=['Feature','k','n','K','M','p_value','enrichment'])
          .sort_values('p_value'))
    # df['adj_p'] = multipletests(df['p_value'], method='fdr_bh')[1]
    return df

def split_and_clean(annotations):
    if isinstance(annotations, str):
        split_anns = annotations.split(';')
        # Remove any leading or trailing whitespace from each annotation
        split_anns = [x.strip() for x in split_anns]
        # Remove any numbers that are at the end of the string
        split_anns = [x.rstrip('0123456789') for x in split_anns]
        # Remove trailing whitespace again
        split_anns = [x.strip() for x in split_anns]
        # Remove any empty strings
        split_anns = [x for x in split_anns if x != '']
        # Remove any annotations that begin with a number
        split_anns = [x for x in split_anns if not x[0].isdigit()]
        return set(split_anns)
    else:
        return set()

RESULT_COLUMNS = ['Bait', 'Feature', 'Feature_type', 'k', 'n', 'K', 'M',
                  'p_value', 'enrichment', 'adj_p']


def process_refactored(data, columns_for_analysis, thresholds):
    """Test each bait's high-confidence preys for enrichment of each feature type.

    ``thresholds`` names the four scores the rest of the application filters on, in the
    form ``apply_score_thresholds`` takes; the foreground is the preys of one bait that
    pass all of them, tested against every prey seen in the run.
    """
    # Get the unique experiments from the data
    experiments = list(data['Experiment.ID'].unique())

    # Get the unique proteins from the data
    all_proteins = set(data['Prey.ID'].unique())

    # Results accumulate in a list and are concatenated once.  Concatenating onto an
    # empty frame instead would leave the count columns as object dtype.
    frames = []

    for column in columns_for_analysis:
        # Create a feature map for this feature type
        feature_df = data[['Prey.ID', column]].copy()
        feature_df.loc[:, 'list'] = feature_df[column].apply(split_and_clean)
        feature_map = dict(zip(feature_df['Prey.ID'], feature_df['list']))

        for experiment in experiments:
            foreground = data[data['Experiment.ID'] == experiment]
            foreground = apply_score_thresholds(foreground, thresholds)
            foreground_ids = set(foreground['Prey.ID'].unique())

            result = enrich_foreground(foreground_ids, all_proteins, feature_map)
            if result.empty:
                continue

            # Add information about the bait and feature type
            result['Bait'] = experiment
            result['Feature_type'] = column

            # Correction is within one bait and feature type, not across them.
            result['adj_p'] = multipletests(result['p_value'], method='fdr_bh')[1]

            frames.append(result)

    if not frames:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    return pd.concat(frames, ignore_index=True)[RESULT_COLUMNS]


MAX_LABEL_CHARS = 42


def _truncate_label(label):
    """Shorten a feature name to something a heatmap row can carry.

    Feature names run to hundreds of characters; drawn in full they take the width the
    map itself needs.
    """
    label = str(label)
    if len(label) <= MAX_LABEL_CHARS:
        return label
    return label[:MAX_LABEL_CHARS - 3] + '...'


def plot_results(results, feature_type, num_features=30):
    # Filter the results for the specific feature type
    filtered_results = results[results['Feature_type'] == feature_type]

    # Get a list of the top features that are passing a threshold
    thresholded_results = filtered_results[filtered_results['adj_p'] <= 0.05]
    thresholded_results = thresholded_results[thresholded_results['enrichment'] >= 2]

    # Get a counts of the number of times each feature is present
    feature_counts = thresholded_results['Feature'].value_counts()

    # Get the top N features
    top_features = feature_counts.head(num_features).index.tolist()

    # Filter the results again for these features
    filtered_results = filtered_results[filtered_results['Feature'].isin(top_features)]

    # Turn this long format into a wide format for plotting
    filtered_results = filtered_results.pivot(index='Feature', columns='Bait', values='enrichment')

    # Fill NaN values with 1 (no enrichment)
    filtered_results = filtered_results.fillna(1)

    # Make anything less than 1 equal to 1 (no enrichment)
    filtered_results[filtered_results < 1] = 1

    # Convert the enrichment to log2 scale
    filtered_results = np.log2(filtered_results)

    grid = sns.clustermap(filtered_results, cmap='viridis', figsize=(12, 8),
                          dendrogram_ratio=(0.18, 0.18), cbar_pos=None)

    ax = grid.ax_heatmap
    ax.set_xlabel("Bait")
    ax.set_ylabel("Feature")

    # data2d holds the data in the order the clustering put it in.
    ordered_columns = grid.data2d.columns
    ordered_rows = grid.data2d.index

    ax.set_xticks(np.arange(len(ordered_columns)) + 0.5)
    ax.set_xticklabels(ordered_columns, rotation=45, ha='right', fontsize=8)

    ax.set_yticks(np.arange(len(ordered_rows)) + 0.5)
    ax.set_yticklabels([_truncate_label(label) for label in ordered_rows], fontsize=7)

    # The colorbar goes in the corner the two dendrograms leave empty, as a cell of the
    # clustermap's own grid rather than as a free-floating inset: only an axes the grid
    # owns is moved when the layout below is recomputed.
    corner = grid.gs[0, 0].subgridspec(2, 1, height_ratios=[2, 1])
    grid.ax_cbar = grid.cax = grid.figure.add_subplot(corner[1])
    grid.figure.colorbar(ax.collections[0], cax=grid.cax, orientation='horizontal')
    grid.cax.xaxis.set_ticks_position('bottom')
    grid.cax.tick_params(labelsize=7)
    # Anchored to the left edge of the bar rather than centered on it: the corner is
    # narrower than the caption, and a centered caption overhangs the canvas.
    grid.cax.set_title('log2 enrichment', loc='left', fontsize=8)

    # seaborn leaves a placeholder engine behind, which freezes the positions it computed
    # for figsize above.  Shiny resizes the figure to the browser card before drawing it
    # and only substitutes an engine of its own when there is none, so without a real one
    # here the margins stay sized for a figure the plot is never drawn at and the labels
    # fall off the canvas.
    grid.figure.set_layout_engine("tight")

    return grid


def main():   
    # Arguments
    parser = argparse.ArgumentParser(description="Enrichment analysis of interaction data")

    parser.add_argument("--input", help="Path to the input file", required=True)
    parser.add_argument("--output", help="Path to the output directory", required=True)
    parser.add_argument("--threshold", help="Minimum SAINT score for the foreground",
                        type=float, default=0.9)
    parser.add_argument("--bfdr", help="Maximum BFDR for the foreground",
                        type=float, default=1.0)
    parser.add_argument("--wd", help="Minimum WD score for the foreground",
                        type=float, default=0.0)
    parser.add_argument("--wdfdr", help="Maximum WDFDR for the foreground",
                        type=float, default=1.0)

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input
    output_dir = args.output
    thresholds = {'SaintScore': args.threshold, 'BFDR': args.bfdr,
                  'WD': args.wd, 'WDFDR': args.wdfdr}

    columns_for_analysis = ['GO_CC', 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains']


    # Load the dataset
    data = pd.read_csv(input_file, sep=",")

    # feature_df, results = process_data(data, columns_for_analysis, threshold)
    results = process_refactored(data, columns_for_analysis, thresholds)

    # test = plot_results(results, 'Domains', num_features=30)
    for feature in columns_for_analysis:
        heatmap = plot_results(results, feature, num_features=30)

        # Save the figure to a file
        plt.savefig(output_dir + f"{feature}_enrichment_analysis.png", dpi=600, bbox_inches='tight')
        plt.close()

    # Show the plot
    plt.show()

    # feature_df.to_csv(output_dir + "Feature_DF_new.csv", index=False)
    results.to_csv(output_dir + "Results_new.csv", index=False)

    # Plot the results

if __name__ == "__main__":
    main()