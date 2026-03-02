import numpy as np
import pandas as pd
from scipy.optimize import minimize
import time
from scipy.stats import lognorm
import os
import math
import matplotlib.pyplot as plt
from patsy import dmatrix
from sklearn.linear_model import LinearRegression
from statsmodels.api import GLM, families


def fit_and_plot_spline(prey_stats, title, df=5):
    """Fit a natural cubic spline with both OLS and logistic regression and plot scatter + fits."""
    x = prey_stats['mean_log2'].values
    y = prey_stats['missingness'].values

    spline_basis = np.asarray(dmatrix("cr(x, df={})".format(df), {"x": x}))
    x_smooth = np.linspace(x.min(), x.max(), 300)
    spline_smooth = np.asarray(dmatrix("cr(x_smooth, df={}, lower_bound={}, upper_bound={})".format(df, x.min(), x.max()),
                                       {"x_smooth": x_smooth}))

    # OLS spline (red)
    ols_model = LinearRegression().fit(spline_basis, y)
    y_ols = ols_model.predict(spline_smooth)

    # Logistic regression spline (blue)
    y_clamped = np.clip(y, 1e-6, 1 - 1e-6)
    glm_model = GLM(y_clamped, spline_basis, family=families.Binomial()).fit()
    y_glm = glm_model.predict(spline_smooth)

    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, alpha=0.5, s=10, color='gray')
    plt.plot(x_smooth, y_ols, color='red', linewidth=2, label='OLS')
    plt.plot(x_smooth, y_glm, color='blue', linewidth=2, label='Logistic')
    plt.xlabel('Mean Log2 Intensity')
    plt.ylabel('Missingness Proportion')
    plt.title(title)
    plt.ylim(-0.05, 1.05)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    # Absolute paths only for testing
    interaction_path = 'Example_datasets\HCM_LFQ\interaction.txt'
    ed_path = 'Example_datasets\HCM_LFQ\ED_clean.csv'


    interaction  = pd.read_csv(interaction_path, sep='\t',header=None)
    # Create column names
    interaction.columns = ['ExperimentID', 'Bait', 'Prey', 'Intensity']

    # Read in ED data
    ed = pd.read_csv(ed_path)

    # Make a dictionary mapping all Baits to their BaitID
    bait_dict = {}
    baits = ed['Bait']
    Bait_ids = ed['Bait ID']
    for i in range(len(baits)):
        bait_dict[baits[i]] = Bait_ids[i]


    # Get control experiments, grouped by Bait from ED file
    controls = ed[ed['Type'] == 'C']
    control_interaction = interaction[interaction['ExperimentID'].isin(controls['Experiment Name'])]

    # Plot for each control bait
    control_baits = controls.groupby('Bait')['Experiment Name'].agg(list)
    for bait, experiments in control_baits.items():
        n_exp = len(experiments)
        if n_exp < 2:
            continue

        bait_interaction = control_interaction[control_interaction['ExperimentID'].isin(experiments)]

        prey_stats = bait_interaction.groupby('Prey')['Intensity'].agg(
            missingness=lambda x: (x == 0).mean(),
            mean_log2=lambda x: np.log2(x[x > 0]).mean()
        ).dropna()

        if len(prey_stats) < 10:
            print('Skipping {} - only {} preys'.format(bait, len(prey_stats)))
            continue
        fit_and_plot_spline(prey_stats, '{} (n={})'.format(bait, n_exp))

    print('Done!')



if __name__ == '__main__':
    main()
