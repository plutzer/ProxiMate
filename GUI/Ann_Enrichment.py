import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import argparse
from collections import Counter
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
        feature_df.loc[:, 'list'] = feature_df[column].apply(lambda x: set(x.split(';')) if isinstance(x, str) else set())
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

            # Concatenate the results to the temp_df
            if not result.empty:
                temp_df = pd.concat([temp_df, result], ignore_index=True)
        
        # Adjust the p-value for this feature type
        temp_df['adj_p'] = multipletests(temp_df['p_value'], method='fdr_bh')[1]
    
        # Add the temp_df to the results dataframe
        results = pd.concat([results, temp_df], ignore_index=True)

    return results
        
def process_data(data, columns_for_analysis, threshold):

    experiments = list(data['Experiment.ID'].unique())

    # Initialize a dataframe with the following columns: Bait, Feature, Feature_type, Foreground, Background, p_value
    results = pd.DataFrame(columns=['Bait', 'Feature', 'Feature_type', 'Foreground', 'Background', 'Total_foreground', 'Total_background', 'p_value'])

    for experiment in experiments:
        print('Analyzing experiment:', experiment)
        # Get the data for that bait
        bait_data = data[data['Experiment.ID'] == experiment]

        # Split the data into foreground and background
        foreground = bait_data[bait_data['SaintScore'] >= threshold]
        background = bait_data[bait_data['SaintScore'] < threshold]

        # Now iterate through the columns of interest
        for column in columns_for_analysis:
            print('Analyzing column:', column)
            # Get the features for the foreground
            foreground_features = []
            for item in list(foreground[column]):
                if item is not np.nan:
                    foreground_features.extend(list(set(item.split(';'))))

            # Now the same for the background
            background_features = []
            for item in list(background[column]):
                if item is not np.nan:
                    background_features.extend(list(set(item.split(';'))))

            # Now statistical tests for each of the features
            foreground_counts = pd.Series(foreground_features).value_counts()
            background_counts = pd.Series(background_features).value_counts()

            features = list(foreground_counts.index)
            counts = list(foreground_counts.values)
            back_counts = [background_counts[feature] if feature in background_counts else 0 for feature in features]

            # Now assemble these lists into a dataframe to perform the statistical test in a vectorized manner
            feature_df = pd.DataFrame({'Feature': features, 'Foreground': counts, 'Background': back_counts})

            # Get the total counts for foreground and background
            total_foreground = len(foreground)
            total_background = len(background)

            # Add these to the dataframe
            feature_df['Total_foreground'] = total_foreground
            feature_df['Total_background'] = total_background

            # Now perform the statistical test for each row in the dataframe
            p_values = []
            for index, row in feature_df.iterrows():
                p_values.append(hypergeom.sf(row['Foreground'], row['Total_foreground'] + row['Total_background'], row['Foreground'] + row['Background'], row['Total_foreground']))

            # Now add the p_values to the dataframe
            feature_df['p_value'] = p_values

            # Add the bait and feature type to the dataframe
            feature_df['Bait'] = experiment
            feature_df['Feature_type'] = column

            # Calculate an enrichment score based on relative abundance in foreground vs background
            epsilon = 1e-10
            feature_df['Enrichment'] = feature_df['Foreground'] / (feature_df['Total_foreground'] + epsilon) / ((feature_df['Background'] + epsilon) / (feature_df['Total_background'] + epsilon))

            # Calculate an adjusted p-value

            # First remove any rows with nan p-values
            feature_df = feature_df.dropna(subset=['p_value'])

            # Sort the dataframe by p-value
            feature_df = feature_df.sort_values('p_value')

            # Add a rank column
            feature_df['Rank'] = range(1, len(feature_df) + 1)

            # Calculate the adjusted p-value
            feature_df['Adjusted_p_value'] = feature_df['p_value'] * len(feature_df) / feature_df['Rank']

            # Now add the results to the main dataframe
            results = pd.concat([results, feature_df])

    return feature_df, results

def plot_results(results, feature_type, num_features=30):
    # Filter the results for the specific feature type
    filtered_results = results[results['Feature_type'] == feature_type]

    # Get a list of the top features that are passing a threshold
    thresholded_results = filtered_results[filtered_results['adj_p'] <= 0.05]

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
    heatmap = sns.clustermap(filtered_results.T, cmap='viridis', cbar_kws={'label': 'Enrichment Score'})
    plt.title(f"Enrichment Analysis for {feature_type}")
    plt.xlabel("Bait")
    plt.ylabel("Feature")

    ax = heatmap.ax_heatmap
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=8, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=5)
    # plt.xticks(rotation=45, ha='right', fontsize=8)
    # plt.yticks(fontsize=8)
    plt.tight_layout()

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

    test = plot_results(results, 'GO_CC')

    # Show the plot
    plt.show()

    # # feature_df.to_csv(output_dir + "Feature_DF_new.csv", index=False)
    # results.to_csv(output_dir + "Results_new.csv", index=False)

    # Plot the results

if __name__ == "__main__":
    main()