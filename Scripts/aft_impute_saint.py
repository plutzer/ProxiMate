import numpy as np
import pandas as pd
from scipy.optimize import minimize
import time
from scipy.stats import lognorm
import os
import cProfile
import math


# Functions for lognormal distribution
def norm_cdf(x, mu, sigma):
    '''
    Compute the cumulative density function of the normal distribution for a given x, mean and standard deviation
    '''
    return (1 + math.erf((x - mu) / np.sqrt(2 * sigma**2))) / 2

# Need a helper function that computes the probability density function of the normal distribution for a given x, mean and standard deviation
def norm_pdf(x, mu, sigma):
    '''
    Compute the probability density function of the normal distribution for a given x, mean and standard deviation
    '''
    return np.exp(-(x - mu)**2 / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi))


def sample_likelihood(intensity, mu, sigma, Tlim):
    '''
    Compute the likelihood of a single sample j
    '''
    if intensity > 0:
        return np.log10(norm_pdf(np.log10(intensity),mu,sigma) + 0.000001)
    else:
        return np.log10(norm_cdf(Tlim,mu,sigma) + 0.000001)

def likelihood(args, prey_intensities, Tlim, n, obs_mu, obs_sigma):
    '''
        Add documentation
    '''
    return -np.sum([sample_likelihood(x, args[0], args[1], Tlim) for x in prey_intensities] + [n*np.log10(norm_pdf(args[0], obs_mu, obs_sigma) + 0.000001)])

def get_initial_params(prey_intensities):
    '''
        Add documentation
    '''
    mu = np.mean(np.log10(prey_intensities.replace(0,np.nan).dropna()))
    sigma = np.std(np.log10(prey_intensities.replace(0,np.nan).dropna()))
    if sigma == 0:
        sigma = 0.5
    return (mu, sigma)
    

#### MAIN  function ####

def filter_impute(prey_path,interaction_path,output_dir,ed_path,impute=False):
    # Read in interaction data
    interaction = pd.read_csv(interaction_path, sep='\t',header=None)
    # Create column names
    interaction.columns = ['ExperimentID', 'Bait', 'Prey', 'Intensity']

    df = interaction[interaction['Intensity'] > 0]
    # y = df.groupby('Prey')['Intensity'].apply(lambda s: np.log10(s).mean())

    # mu0 = float(y.median())
    # mad = float((y - mu0).abs().median())
    # tau  = float(max(1e-3, 1.4826 * mad))   # MAD→SD

    # global_obs_mu    = mu0
    # global_obs_sigma = tau
    # print(f"Global prior (robust): mu0={mu0:.3f}, tau={tau:.3f}")

    ### Code for tunable weight
    prey_counts = interaction.groupby('Prey').agg({'Intensity': 'count'})
    n0 = prey_counts['Intensity'].median()

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
        prey_data = pd.read_csv(prey_path, sep='\t',header=None)
        # Change the column names
        prey_data.columns = ['PreyID', 'PreyGene']
        # Add a column for the mu parameter for each prey
        prey_data['mu'] = 0.0

        # Start a timer to track how long the imputation takes
        start = time.time()

        # First need an upper-bound estimate for Tlim based on the distribution of intensities
        upper_Tlim = np.log10(np.percentile(interaction.replace(0,np.nan).dropna().groupby('Prey').agg({'Intensity': 'min'})['Intensity'],95))
        print("Upper_Tlim: " + str(upper_Tlim))

        preys = prey_data['PreyID'].tolist()

        mu_list = []
        sigma_list = []
        tlim_list = []
        imputed = []
        iterations = []
        original_b = []

        mu_lower_bound = np.log10(np.min(interaction.replace(0,np.nan).dropna()['Intensity'])) - np.std(np.log10(interaction.replace(0,np.nan).dropna()['Intensity']))*3
        mu_upper_bound = np.log10(np.percentile(interaction.replace(0,np.nan).dropna().groupby('Prey').agg({'Intensity': 'mean'})['Intensity'],95))

        interaction_log = interaction.replace(0,np.nan).dropna().copy()
        interaction_log['Intensity_log'] = np.log10(interaction_log['Intensity'])
        sigma_lower_bound = np.percentile(interaction_log.groupby('Prey').agg({'Intensity_log': 'std'})['Intensity_log'].dropna(),5)
        
        #Loop through the preys and impute the intensity values
        for prey in preys:
            # get a vector of the intensities for the current prey
            prey_interaction = interaction.loc[interaction['Prey'] == prey].copy()

            full_mu = np.mean(np.log10(prey_interaction['Intensity'].replace(0,np.nan).dropna()))

            # Add a new column to prey_interaction with the Bait ID from experimentalDesign file
            prey_interaction['BaitID'] = prey_interaction['Bait'].map(bait_dict).copy()

            # Remove rows where BaitID is the same as prey
            prey_interaction = prey_interaction[prey_interaction['BaitID'] != prey]

            prey_intensities = prey_interaction['Intensity']

            Tlim = np.min(pd.concat([np.log10(prey_intensities.replace(0, np.nan).dropna()), pd.Series([upper_Tlim])])) # Change this to a if statement to track how many times upper_Tlim is used

            # If prey intensity is all zeros, skip the optimization
            if len(prey_intensities.replace(0,np.nan).dropna()) == 0:
                print("All zeros: " + prey)
                mu_list.append(Tlim)
                sigma_list.append(np.nan)
                tlim_list.append(Tlim)
                iterations.append(0)
                original_b.append(0)
                imputed.append(False)
                prey_data.loc[prey_data['PreyID'] == prey, 'mu'] = 0.0
                continue
            
            # Initialize parameters to optimize
            (a,b) = get_initial_params(prey_intensities)

            # Change prey_intensities to a list
            prey_intensities = list(prey_intensities)

            n = len(prey_intensities) - prey_intensities.count(0)

            # Tunable weight
            w = n/(n + n0)

            # Calculate the mean and sd of the log10 intensities excluding zeros
            log_prey_intensities = [np.log10(x) for x in prey_intensities if x != 0]
            obs_mu = np.mean(log_prey_intensities)
            obs_sigma = np.std(log_prey_intensities) # TODO: What happens if this is 0?
            if obs_sigma == 0:
                obs_sigma = 0.5

            bounds = ((mu_lower_bound,mu_upper_bound),(sigma_lower_bound,sigma_lower_bound*2 + 2*obs_sigma))
            # Check to see if the initial parameters are within the bounds
            if a < bounds[0][0]:
                print("starting mu too low: " + prey, a, bounds[0][0])
                a = bounds[0][0]
            if a > bounds[0][1]:
                print("starting mu too high: " + prey, a, bounds[0][1])
                a = bounds[0][1]
            if b < bounds[1][0]:
                print("starting sigma too low: " + prey, b, bounds[1][0])
                b = bounds[1][0]
            if b > bounds[1][1]:
                print("starting sigma too high: " + prey, b, bounds[1][1])
                b = bounds[1][1]

            # Optimize the parameters
            res = minimize(likelihood, [a,b], args=(prey_intensities, Tlim, w, obs_mu, obs_sigma), method='Nelder-Mead', tol=1e-3, bounds = bounds)

            optimized_mu = res.x[0]
            optimized_sigma = res.x[1]

            if np.isnan(optimized_mu):
                print("NaN mu: " + prey)
                optimized_mu = full_mu

            # # Make optimized mu the min of Tlim and the optimized mu
            # optimized_mu = min(optimized_mu, Tlim)

            # Append the optimized parameters to the lists
            mu_list.append(optimized_mu)
            sigma_list.append(optimized_sigma)
            tlim_list.append(Tlim)
            iterations.append(res.nit)
            original_b.append(b)

            # Print an update every 100 preys
            if len(mu_list) % 100 == 0:
                print(len(mu_list))

            # Change the intensity values for the current prey to the imputed values for the selected control bait
            if impute:
                prey_data.loc[prey_data['PreyID'] == prey, 'mu'] = 10**optimized_mu
                imputed.append(True)

        # Remove imputation from preys that didn't fit appropriately
        # Calculate the mean and sd of the imputed values
        mean_mu = np.mean(mu_list)
        sd_mu = np.std(mu_list)
        # Calculate the mean and sd of the sigma values
        mean_sigma = np.mean(sigma_list)
        sd_sigma = np.std(sigma_list)
        
        # iterate through the prey data
        for i in range(len(mu_list)):
            # Check to see if mu is within 3 standard deviations of the mean
            if (mu_list[i] < mean_mu - 3*sd_mu) or (mu_list[i] > mean_mu + 3*sd_mu):
                prey_data.loc[i, 'mu'] = 0.0
                mu_list[i] = 0.0
                imputed[i] = False
            # Check to see if sigma is within 3 standard deviations of the mean
            if (sigma_list[i] < mean_sigma - 3*sd_sigma) or (sigma_list[i] > mean_sigma + 3*sd_sigma):
                prey_data.loc[i, 'mu'] = 0.0
                mu_list[i] = 0.0
                imputed[i] = False


        # Stop the timer
        end = time.time()
        print(end - start)

        # Write a csv output file using the prey names and the optimized parameters
        output = pd.DataFrame({'Prey': preys[:len(mu_list)], 'mu': mu_list, 'sigma': sigma_list, 'originalSigma':original_b, 'Tlim': tlim_list, 'iterations': iterations, 'imputed': imputed}) #Edited for debugging
        if impute:
            output.to_csv(output_dir + 'imputed_params.csv', index=False)

    # filter the interaction file to remove zero intensity values
    interaction = interaction[interaction['Intensity'] > 0]

    if impute:
        prey_data.to_csv(output_dir + 'imputed_prey.txt', sep='\t', index=False, header=False)
    interaction.to_csv(output_dir + 'filtered_interaction.txt', sep='\t', index=False, header=False)

def main(): # For testing purposes
    prey_path = 'C:/Users/isaac/Work/20241024_ScoringFixing/Impute_bounds_testing/prey.txt'
    interaction_path = 'C:/Users/isaac/Work/20241024_ScoringFixing/Impute_bounds_testing/interaction.txt'
    output_dir = 'C:/Users/isaac/Work/20241024_ScoringFixing/Impute_bounds_testing/test_output_boundfix/'
    ed_path = 'C:/Users/isaac/Work/20241024_ScoringFixing/Impute_bounds_testing/ED.csv'

    profiler = cProfile.Profile()

    # Check to see if the output directory exists, if not create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    profiler.enable()
    filter_impute(prey_path,interaction_path,output_dir,ed_path,impute=True)
    profiler.disable()

    profiler.dump_stats('C:/Users/isaac/Work/20241024_ScoringFixing/Impute_bounds_testing/profile.prof')


# Add the if name main statement
if __name__ == '__main__':
    main()