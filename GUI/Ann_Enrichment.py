import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import argparse
from collections import Counter
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

def process_refactored(data, columns_for_analysis, threshold):
    # Get the unique experiments from the data
    experiments = list(data['Experiment.ID'].unique())

    # Initialize a dataframe for the results
    results = pd.DataFrame(columns=['Bait', 'Feature', 'Feature_type', 'k','n','K','M','p_value','enrichment', 'adj_p'])

    # Get the unique proteins from the data
    all_proteins = set(data['Prey.ID'].unique())

    for column in columns_for_analysis:
        # Create a temporary dataframe to store results for this feature type
        temp_df = pd.DataFrame(columns=['Bait', 'Feature', 'Feature_type', 'k','n','K','M','p_value','enrichment'])

        # Create a feature map for this feature type
        feature_map = {}
        feature_df = data[['Prey.ID', column]].copy()
        feature_df.loc[:, 'list'] = feature_df[column].apply(split_and_clean)
        feature_map = dict(zip(feature_df['Prey.ID'], feature_df['list']))

        print(len(feature_map), " features in feature map")

        for experiment in experiments:
            print('Analyzing experiment:', experiment, 'for feature type:', column)
            foreground = data[data['Experiment.ID'] == experiment]
            foreground = foreground[foreground['SaintScore'] >= threshold]
            foreground_ids = set(foreground['Prey.ID'].unique())

            result = enrich_foreground(foreground_ids, all_proteins, feature_map)

            # Add information about the bait and feature type
            result['Bait'] = experiment
            result['Feature_type'] = column

            # calculate an adjusted p-value for the results
            if len(result) > 0:
                result['adj_p'] = multipletests(result['p_value'], method='fdr_bh')[1]
            else:
                continue

            # Concatenate the results to the temp_df
            if not result.empty:
                temp_df = pd.concat([temp_df, result], ignore_index=True)
        
        # Add the temp_df to the results dataframe
        results = pd.concat([results, temp_df], ignore_index=True)

    return results

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

    # Initialize the plot with custom size
    # plt.figure(figsize=(10, 28))

    heatmap = sns.clustermap(filtered_results, cmap='viridis', cbar_kws={'label': 'Enrichment Score'}, figsize =  (28, 10))
    # plt.title(f"Enrichment Analysis for {feature_type}")
    plt.xlabel("Bait")
    plt.ylabel("Feature")
    heatmap.cax.set_ylabel('log2 Enrichment Score', rotation=270, labelpad=15)

    ax = heatmap.ax_heatmap
    # ax.set_xticklabels(ax.get_xticklabels(), fontsize=5, rotation=45, ha='right')
    # ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)

    # Use the reordered DataFrame from the clustermap (data2d holds the data with proper ordering)
    ordered_columns = heatmap.data2d.columns

    # # Set tick positions corresponding to every column
    ax.set_xticks(np.arange(len(ordered_columns)) + 0.5)
    ax.set_xticklabels(ordered_columns, rotation=45, ha='right', fontsize=5)


    # Force the locator to show every tick
    # ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(len(ordered_columns))))

    # plt.xticks(rotation=45, ha='right', fontsize=8)
    # plt.yticks(fontsize=8)
    # plt.tight_layout()

    return heatmap


def main():   
    # Arguments
    parser = argparse.ArgumentParser(description="Enrichment analysis of interaction data")

    parser.add_argument("--input", help="Path to the input file", required=True)
    parser.add_argument("--output", help="Path to the output directory", required=True)
    parser.add_argument("--threshold", help="Threshold for enrichment analysis", type=float, default=0.9)

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input
    output_dir = args.output
    threshold = args.threshold

    print("Started Script")

    columns_for_analysis = ['GO_CC', 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains']


    # Load the dataset
    data = pd.read_csv(input_file, sep=",")

    # feature_df, results = process_data(data, columns_for_analysis, threshold)
    results = process_refactored(data, columns_for_analysis, threshold)

    # test = plot_results(results, 'Domains', num_features=30)
    for feature in columns_for_analysis:
        print("Plotting results for feature type:", feature)
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