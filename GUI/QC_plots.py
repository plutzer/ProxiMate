import plotly.express as px
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import plotly.graph_objects as go



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

    print(len(data), " rows in data")
    # Convert all 0 values to NaN
    data = data.replace(0, np.nan)

    # Remove rows with more than 75% NaN values
    data = data.dropna(thresh=len(data.columns) * 0.5)

    print(len(data), " rows in data after removing rows with more than 75% NaN values")

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

    pca_df['PC1'] = pd.to_numeric(pca_df['PC1'], errors='coerce')
    pca_df['PC2'] = pd.to_numeric(pca_df['PC2'], errors='coerce')

    # Make a plotly scatterplot of the PCA results
    # fig = px.scatter(pca_df, x='PC1', y='PC2', color='BaitName', symbol='Type',
    #                              hover_name='Experiment', hover_data=['Experiment', 'BaitName'],
    #                              title="PCA of Interaction Data")
    print(pca_df)

    # fig_widget = go.FigureWidget(fig)

    return pca_df

def saint_known_retention(results_path, ctrl_experiments=None):

    results = pd.read_csv(results_path, sep=",")

    # If ctrl is used, filter the results
    if ctrl_experiments is not None:
        results = results[results['Experiment'].isin(ctrl_experiments)]

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

    return fig




def main():
    base_dir = 'C:/Users/plutzer/Work/ProxiMate_Testing/'

    # For testing PCA
    # figure = pca_plot(base_dir + "output/interaction.txt", base_dir + "output/ED.csv")

    # For testing known retention
    figure = saint_known_retention(base_dir + "output/annotated_scores.csv")

    figure.show()



if __name__ == "__main__":
    main()