import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom

threshold = 0.9

def process_data(data, columns_for_analysis):

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


def main():
    
    print("Started Script")

    out_dir = "C:/Users/isaac/Work/025_MainNetwork/Enrichment/"
    columns_for_analysis = ['GO_CC', 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains']


    # Load the dataset
    data = pd.read_csv("C:/Users/isaac/Work/025_MainNetwork/1_MainNetwork_aggregated_scores.csv")

    feature_df, results = process_data(data, columns_for_analysis)

    feature_df.to_csv(out_dir + "Feature_DF.csv", index=False)
    results.to_csv(out_dir + "Results.csv", index=False)


if __name__ == "__main__":
    main()