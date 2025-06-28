import plotly.express as px
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import plotly.graph_objects as go
from sklearn.metrics import roc_curve, auc



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
    print(pca_df)
        


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
        scores = scores[scores['Experiment'].isin(ctrl_experiments)]

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




def main():
    # base_dir = 'C:/Users/plutzer/Work/ProxiMate_Testing/'

    # For testing PCA
    # figure = pca_plot(base_dir + "output/interaction.txt", base_dir + "output/ED.csv")

    # For testing known retention
    figure = roc_plot("C:/Users/isaac/Work/025_MainNetwork/1_MainNetwork_aggregated_scores.csv", 'Multivalidated')

    figure.show()
    print('done')


if __name__ == "__main__":
    main()