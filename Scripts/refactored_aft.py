import numpy as np
import pandas as pd
from scipy.optimize import minimize
import time
from scipy.stats import lognorm
from scipy.stats import norm
import os
import math
import matplotlib.pyplot as plt
from patsy import dmatrix
from sklearn.linear_model import LinearRegression
from statsmodels.api import GLM, families

# For development, I'll start with hardcoded values for pi
# pi = 0.342 # OLS Avg
# pi = 0.356 # OLS Weighted Avg
# pi = 0.297 # GLM Avg
pi = 0.313 # GLM Weighted Avg

def protein_log_likelihood(prey_intensities, mu, sigma, Tlim, pi):
    y = np.asarray(prey_intensities)
    observed = y > 0
    ll = np.sum(np.log(1/sigma) - (1/(2*sigma**2)) * (y[observed] - mu)**2)
    n_zero = np.sum(~observed)
    if n_zero > 0:
        ll += n_zero * np.log(pi + (1-pi) * norm.cdf((Tlim-mu)/sigma))
    return ll

def neg_likelihood(params, prey_intensities, Tlim, pi):
    mu, sigma = params
    return -protein_log_likelihood(prey_intensities, mu, sigma, Tlim, pi)

def fit_and_plot_spline(prey_stats, title, df=5, eval_type='ols'):
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

    # Evaluate the spline at its maximum intensity
    if eval_type == 'ols':
        return y_ols[-1]
    elif eval_type == 'glm':
        return y_glm[-1]
    else:
        raise ValueError("eval_type must be 'ols' or 'glm'")


def get_initial_params(prey_intensities):
    mu = np.mean(np.log10(prey_intensities.replace(0, np.nan).dropna()))
    sigma = np.std(np.log10(prey_intensities.replace(0, np.nan).dropna()))
    if sigma == 0:
        sigma = 0.5
    return (mu, sigma)


#### MAIN function ####

def filter_impute(prey_path, interaction_path, output_dir, ed_path, impute=False):
    # Read in interaction data
    interaction = pd.read_csv(interaction_path, sep='\t', header=None)
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

    if impute:
        # Read in prey data
        prey_data = pd.read_csv(prey_path, sep='\t', header=None)
        # Change the column names
        prey_data.columns = ['PreyID', 'PreyGene']
        # Add a column for the mu parameter for each prey
        prey_data['mu'] = 0.0

        # Start a timer to track how long the imputation takes
        start = time.time()

        # Precompute nonzero interaction data and log10 intensities once
        interaction_nonzero = interaction[interaction['Intensity'] > 0].copy()
        interaction_nonzero['Intensity_log'] = np.log10(interaction_nonzero['Intensity'])

        # First need an upper-bound estimate for Tlim based on the distribution of intensities
        upper_Tlim = np.log10(np.percentile(interaction_nonzero.groupby('Prey')['Intensity'].min(), 95))
        preys = prey_data['PreyID'].tolist()

        mu_list = []
        sigma_list = []
        tlim_list = []
        imputed = []
        iterations = []
        original_b = []
        pi_list = []

        mu_lower_bound = interaction_nonzero['Intensity_log'].min() - interaction_nonzero['Intensity_log'].std() * 5
        mu_upper_bound = np.log10(np.percentile(interaction_nonzero.groupby('Prey')['Intensity'].mean(), 99))
        sigma_lower_bound = np.percentile(interaction_nonzero.groupby('Prey')['Intensity_log'].std().dropna(), 1)

        # Add BaitID column to full interaction once, before the loop
        interaction['BaitID'] = interaction['Bait'].map(bait_dict)

        # Pre-group interaction by Prey for fast lookup
        grouped = interaction.groupby('Prey')

        # Loop through the preys and impute the intensity values
        for n, prey in enumerate(preys, 1):
            print(f"Processing prey {prey}... ({n}/{len(preys)})")

            # Fast group lookup instead of scanning the full DataFrame
            try:
                prey_interaction = grouped.get_group(prey).copy()
            except KeyError:
                mu_list.append(upper_Tlim)
                sigma_list.append(np.nan)
                tlim_list.append(upper_Tlim)
                iterations.append(0)
                original_b.append(0)
                pi_list.append(pi)
                imputed.append(False)
                prey_data.loc[prey_data['PreyID'] == prey, 'mu'] = 0.0
                continue

            # Remove rows where BaitID is the same as prey (self-interactions)
            prey_interaction = prey_interaction[prey_interaction['BaitID'] != prey]

            prey_intensities = prey_interaction['Intensity'].values
            nonzero_mask = prey_intensities > 0
            nonzero_vals = prey_intensities[nonzero_mask]

            # Compute Tlim
            if len(nonzero_vals) > 0:
                Tlim = min(np.log10(nonzero_vals.min()), upper_Tlim)
            else:
                Tlim = upper_Tlim

            # Compute full_mu from all intensities (before self-interaction filter)
            full_nonzero = prey_interaction['Intensity'].values
            full_nonzero = full_nonzero[full_nonzero > 0]
            full_mu = np.mean(np.log10(full_nonzero)) if len(full_nonzero) > 0 else Tlim

            # If prey intensity is all zeros, skip the optimization
            if len(nonzero_vals) == 0:
                mu_list.append(Tlim)
                sigma_list.append(np.nan)
                tlim_list.append(Tlim)
                iterations.append(0)
                original_b.append(0)
                pi_list.append(pi)
                imputed.append(False)
                prey_data.loc[prey_data['PreyID'] == prey, 'mu'] = 0.0
                continue

            # Initialize parameters to optimize
            log_nonzero = np.log10(nonzero_vals)
            a = np.mean(log_nonzero)
            b = np.std(log_nonzero)
            if b == 0:
                b = 0.5

            obs_sigma = b

            bounds = ((mu_lower_bound, mu_upper_bound), (sigma_lower_bound, sigma_lower_bound * 3 + 3 * obs_sigma))
            # Check to see if the initial parameters are within the bounds
            a = np.clip(a, bounds[0][0], bounds[0][1])
            b = np.clip(b, bounds[1][0], bounds[1][1])

            # Log10-transform non-zero intensities for the likelihood function, keep zeros as 0
            prey_intensities_log = np.where(nonzero_mask, np.log10(np.where(nonzero_mask, prey_intensities, 1)), 0)

            # Optimize the parameters
            res = minimize(neg_likelihood, [a, b], args=(prey_intensities_log, Tlim, pi), method='Nelder-Mead', tol=1e-3, bounds=bounds)

            optimized_mu = res.x[0]
            optimized_sigma = res.x[1]

            if np.isnan(optimized_mu):
                optimized_mu = full_mu

            # Append the optimized parameters to the lists
            mu_list.append(optimized_mu)
            sigma_list.append(optimized_sigma)
            tlim_list.append(Tlim)
            iterations.append(res.nit)
            original_b.append(b)
            pi_list.append(pi)

            # Change the intensity values for the current prey to the imputed values
            if impute:
                prey_data.loc[prey_data['PreyID'] == prey, 'mu'] = 10**optimized_mu
                imputed.append(True)

        # Remove imputation from preys that didn't fit appropriately
        # Calculate the mean and sd of the imputed values
        mean_mu = np.nanmean(mu_list)
        sd_mu = np.nanstd(mu_list)
        # Calculate the mean and sd of the sigma values
        mean_sigma = np.nanmean(sigma_list)
        sd_sigma = np.nanstd(sigma_list)

        # iterate through the prey data
        for i in range(len(mu_list)):
            # Check to see if mu is within 3 standard deviations of the mean
            if (mu_list[i] < mean_mu - 3 * sd_mu) or (mu_list[i] > mean_mu + 3 * sd_mu):
                prey_data.loc[i, 'mu'] = 0.0
                mu_list[i] = 0.0
                imputed[i] = False
            # Check to see if sigma is within 3 standard deviations of the mean
            if (sigma_list[i] < mean_sigma - 3 * sd_sigma) or (sigma_list[i] > mean_sigma + 3 * sd_sigma):
                prey_data.loc[i, 'mu'] = 0.0
                mu_list[i] = 0.0
                imputed[i] = False

        # Write a csv output file using the prey names and the optimized parameters
        output = pd.DataFrame({'Prey': preys[:len(mu_list)], 'mu': mu_list, 'sigma': sigma_list, 'originalSigma': original_b, 'Tlim': tlim_list, 'iterations': iterations, 'imputed': imputed, 'pi': pi_list})
        if impute:
            output.to_csv(output_dir + 'imputed_params.csv', index=False)

    # filter the interaction file to remove zero intensity values
    interaction = interaction[interaction['Intensity'] > 0]

    if impute:
        prey_data.to_csv(output_dir + 'imputed_prey.txt', sep='\t', index=False, header=False)
    interaction.to_csv(output_dir + 'filtered_interaction.txt', sep='\t', index=False, header=False)


def main():
    interaction_path = 'Example_datasets/HCM_LFQ/interaction.txt'
    ed_path = 'Example_datasets/HCM_LFQ/ED_clean.csv'
    prey_path = 'Example_datasets/HCM_LFQ/prey.txt'
    output_dir = 'Example_datasets/HCM_LFQ/'
    filter_impute(prey_path, interaction_path, output_dir, ed_path, impute=True)


def main_spline_test():
    interaction_path = 'Example_datasets/HCM_LFQ/interaction.txt'
    ed_path = 'Example_datasets/HCM_LFQ/ED_clean.csv'

    interaction = pd.read_csv(interaction_path, sep='\t', header=None)
    interaction.columns = ['ExperimentID', 'Bait', 'Prey', 'Intensity']

    ed = pd.read_csv(ed_path)

    # --- Per-bait control spline fits ---
    controls = ed[ed['Type'] == 'C']
    control_interaction = interaction[interaction['ExperimentID'].isin(controls['Experiment Name'])]

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
        fit_and_plot_spline(prey_stats, 'Control: {} (n={})'.format(bait, n_exp))

    # --- Per-bait test spline fits ---
    tests = ed[ed['Type'] == 'T']
    test_interaction = interaction[interaction['ExperimentID'].isin(tests['Experiment Name'])]

    test_baits = tests.groupby('Bait')['Experiment Name'].agg(list)
    for bait, experiments in test_baits.items():
        n_exp = len(experiments)
        if n_exp < 2:
            continue

        bait_interaction = test_interaction[test_interaction['ExperimentID'].isin(experiments)]

        prey_stats = bait_interaction.groupby('Prey')['Intensity'].agg(
            missingness=lambda x: (x == 0).mean(),
            mean_log2=lambda x: np.log2(x[x > 0]).mean()
        ).dropna()

        if len(prey_stats) < 10:
            print('Skipping {} - only {} preys'.format(bait, len(prey_stats)))
            continue
        fit_and_plot_spline(prey_stats, 'Test: {} (n={})'.format(bait, n_exp))

    # --- Full dataset spline fit ---
    prey_stats_all = interaction.groupby('Prey')['Intensity'].agg(
        missingness=lambda x: (x == 0).mean(),
        mean_log2=lambda x: np.log2(x[x > 0]).mean()
    ).dropna()

    fit_and_plot_spline(prey_stats_all, 'Full Dataset (n={})'.format(len(interaction['ExperimentID'].unique())))

    print('Done!')


if __name__ == '__main__':
    main_spline_test()
