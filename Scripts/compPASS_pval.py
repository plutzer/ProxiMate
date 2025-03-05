import numpy as np
import pandas as pd
import argparse
import cProfile
import time
from statsmodels.stats.multitest import multipletests

# Need a function for entropy calculation
def entropy(xs):
    # Convert input to numpy array for easier manipulation
    xs = np.array(xs, dtype=np.float64)
    # Calculate probabilities
    p = (xs + 1/len(xs)) / (np.sum(xs) + 1)
    # Calculate entropy
    ent = np.sum([-x * np.log2(x) for x in p])
    return ent

def get_ave_psm(input):
    # Constructed in the same way as the original R code.
    dups_dropped = input.groupby(["Experiment.ID", "Prey", "Replicate"]).agg({"Bait": "first", "Spectral.Count" : "max"}).reset_index()
    ave_psm = dups_dropped.groupby(["Experiment.ID", "Prey"]).agg(
        Bait=('Bait', lambda x: x.unique()[0]),
        AvePSM = ('Spectral.Count', 'mean'),
        N_Saw = ('Spectral.Count', lambda x: (x > 0).sum()),
        Entropy = ('Spectral.Count', entropy)
    ).reset_index()
    return ave_psm

def permute_prey_matrices(ave_psm_values, prey_means_matrix, prey_sd_matrix, self_interaction, n_saw_values):
    # Returns a single permutation of prey values (counts), not labels
    # Permutes the prey values (spectral counts) within each prey (row), including zero values,
    # ensuring that the observed bait-prey structure and prey-specific variance are preserved.

    permuted_ave_psm_values = np.zeros_like(ave_psm_values)
    permuted_prey_means_matrix = np.zeros_like(prey_means_matrix)
    permuted_prey_sd_matrix = np.zeros_like(prey_sd_matrix)
    permuted_self_interaction = np.zeros_like(self_interaction)
    permuted_n_saw_values = np.zeros_like(n_saw_values)

    # For each prey (row), permute the values across experiments (columns)
    for i in range(ave_psm_values.shape[0]):
        # Get the spectral counts (values) for this prey across all experiments
        values = ave_psm_values[i, :]
        self_interactions = self_interaction[i, :]
        
        # Permute the values for this prey across experiments
        permuted_indices = np.random.permutation(len(values))
        permuted_values = values[permuted_indices]
        permuted_means = prey_means_matrix[i, :][permuted_indices]
        permuted_sd = prey_sd_matrix[i, :][permuted_indices]
        permuted_self_interactions = self_interactions[permuted_indices]
        permuted_n_saw = n_saw_values[i, :][permuted_indices]
        
        # Insert the permuted values back into the permuted matrices
        permuted_ave_psm_values[i, :] = permuted_values
        permuted_prey_means_matrix[i, :] = permuted_means
        permuted_prey_sd_matrix[i, :] = permuted_sd
        permuted_self_interaction[i, :] = permuted_self_interactions
        permuted_n_saw_values[i, :] = permuted_n_saw

    return permuted_ave_psm_values, permuted_prey_means_matrix, permuted_prey_sd_matrix, permuted_self_interaction, permuted_n_saw_values

def normalize_matrix(matrix, quantile):
    if quantile:
        # Collapse the matrix and remove zeros
        data = matrix.flatten()
        data = data[data > 0]
        # Calculate the value at the quantile
        q = np.quantile(data, quantile)
        # Divide the original matrix by the quantile value
        return matrix / q, q
    else:
        return matrix, None

def calculate_wd_matrix(ave_psm_values, prey_means_matrix, prey_sd_matrix, n_exp_with_prey, n_saw_values, n_experiments, norm_factor=None):
    # Calculate the WD score matrix
    inner_terms = (prey_sd_matrix / prey_means_matrix) * np.tile((n_experiments / n_exp_with_prey)[:, np.newaxis], (1, prey_sd_matrix.shape[1]))
    wd_scores = np.sqrt(ave_psm_values * (inner_terms**n_saw_values))

    normalized_wd_scores, q = normalize_matrix(wd_scores, norm_factor)

    return normalized_wd_scores, q

def calculate_wd_pvals_matrix(normalized_wd_scores, ave_psm_values, prey_means_matrix, prey_sd_matrix, self_interaction, n_exp_with_prey, n_saw_values, n_experiments, ave_psm, rows, columns, iterations, q):

    # Do things one permutation at a time to reduce memory usage
    p_values = np.zeros((normalized_wd_scores.shape[0], normalized_wd_scores.shape[1]))

    # Get a matrix of non-zero ave_psm values

    for i in range(iterations):
        print(f"Permutation {i}")
        # if i % 100 == 0:
        #     print(f"Permutation {i}")
        permuted_ave_psm_values, permuted_prey_means_matrix, permuted_prey_sd_matrix, permuted_self_interaction, permuted_n_saw_values = permute_prey_matrices(ave_psm_values, prey_means_matrix, prey_sd_matrix, self_interaction, n_saw_values)
        permuted_normalized_wd_scores, _ = calculate_wd_matrix(permuted_ave_psm_values, permuted_prey_means_matrix, permuted_prey_sd_matrix, n_exp_with_prey, permuted_n_saw_values, n_experiments, norm_factor=0.98)
        #permuted_normalized_wd_scores, _ = calculate_wd_matrix(permuted_ave_psm_values, permuted_self_interaction, n_exp_with_prey, n_saw_values, n_experiments, ave_psm, rows, columns)
        permuted_normalized_wd_scores = permuted_normalized_wd_scores# / q
        
        # Standard incrementation
        p_values += normalized_wd_scores < permuted_normalized_wd_scores

    p_values = p_values / iterations
    return p_values

def score_compPass(input, norm_factor, iterations=None):

    ave_psm = get_ave_psm(input)

    n_experiments = len(ave_psm['Experiment.ID'].unique())

    # Get the unique values for the columns and rows
    columns = ave_psm["Experiment.ID"].unique()
    rows = ave_psm["Prey"].unique()

    # Create matrices for WD calculations
    ave_psm_values = np.zeros((len(rows), len(columns)))
    n_saw_values = np.zeros((len(rows), len(columns)))
    self_interaction = np.zeros((len(rows), len(columns)))
    # prey_means_matrix = np.zeros((len(rows), len(columns)))
    # prey_sd_matrix = np.zeros((len(rows), len(columns)))

    # Fill in the base matrices
    for i, row in ave_psm.iterrows():
        row_index = np.where(rows == row["Prey"])[0][0]
        col_index = np.where(columns == row["Experiment.ID"])[0][0]
        ave_psm_values[row_index, col_index] = row["AvePSM"]
        n_saw_values[row_index, col_index] = row["N_Saw"]
        self_interaction[row_index, col_index] = (row["Bait"] == row["Prey"]) # This matrix is 1 if the bait and prey are the same protein, 0 otherwise

    # Identify edge cases where the prey is only found with itself as bait
    self_interaction_only = list(np.where((np.sum(self_interaction, axis=1) == np.sum(ave_psm_values > 0, axis=1)) == True)[0])

    # Calculate the number of experiments with each prey
    n_exp_with_prey = np.sum(n_saw_values > 0, axis=1)

    # Correct for edge cases where the prey is only found with itself as bait
    n_exp_with_prey[self_interaction_only] = 1

    # Calculate the mean matrix - Yes the prey means needs to be a matrix because of how the original compPASS code handles self-interactions. A prey can have different 'means' for each bait.
    prey_means_matrix = np.tile(((np.sum(ave_psm_values, axis=1) - np.sum(ave_psm_values * self_interaction, axis=1)) / n_experiments)[:, np.newaxis], (1, ave_psm_values.shape[1]))

    # Update the means matrix for the edge cases
    for prey_index in self_interaction_only:
        prey_means_matrix[prey_index, :] = ave_psm_values[prey_index, :]/n_experiments

    # Calculate the base standard deviation matrix
    sum_sq_err = np.tile(np.sum(((ave_psm_values - (ave_psm_values*self_interaction) + (prey_means_matrix*self_interaction)) - prey_means_matrix)**2, axis=1)[:, np.newaxis], (1, ave_psm_values.shape[1]))

    # Update the sum_sq_err matrix for the edge cases
    for prey_index in self_interaction_only:
        sum_sq_err[prey_index, :] = (ave_psm_values[prey_index, :] - prey_means_matrix[prey_index, :])**2 + (prey_means_matrix[prey_index, :]**2)*(n_experiments)

    # Calculate the standard deviation matrix
    prey_sd_matrix = np.sqrt(sum_sq_err/(n_experiments - 1))

    # Clean up the sd and mean matrices by replacing zeros with NaN
    prey_sd_matrix[prey_sd_matrix == 0] = np.nan
    prey_means_matrix[prey_means_matrix == 0] = np.nan

    # Calculate the z-score matrix
    z_scores = (ave_psm_values - prey_means_matrix) / prey_sd_matrix

    # Calculate the WD score matrix
    normalized_wd_scores, norm_val = calculate_wd_matrix(ave_psm_values, prey_means_matrix, prey_sd_matrix, n_exp_with_prey, n_saw_values, n_experiments, norm_factor)

    # Create a new dataframe to combine with the ave_psm dataframe
    experiment_ids = []
    prey_names = []
    bait_names = []
    z_scores_list = []
    normalized_wd_scores_list = []
    sd_list = []
    mean_list = []
    self_interaction_list = []
    self_only_list = []
    n_saw_vals_list = []
    n_exp_with_prey_list = []

    for i, row in ave_psm.iterrows():
        row_index = np.where(rows == row["Prey"])[0][0]
        col_index = np.where(columns == row["Experiment.ID"])[0][0]
        experiment_ids.append(row["Experiment.ID"])
        prey_names.append(row["Prey"])
        bait_names.append(row["Bait"])
        z_scores_list.append(z_scores[row_index, col_index])
        normalized_wd_scores_list.append(normalized_wd_scores[row_index, col_index])
        sd_list.append(prey_sd_matrix[row_index, col_index])
        mean_list.append(prey_means_matrix[row_index, col_index])
        self_interaction_list.append(self_interaction[row_index, col_index])
        self_only_list.append(row_index in self_interaction_only)
        n_saw_vals_list.append(n_saw_values[row_index, col_index])
        n_exp_with_prey_list.append(n_exp_with_prey[row_index])


    # Add the new columns to ave_psm
    ave_psm["Mean"] = mean_list
    ave_psm["SD"] = sd_list
    ave_psm["Z"] = z_scores_list
    ave_psm["WD"] = normalized_wd_scores_list
    ave_psm["Self.Interaction"] = self_interaction_list
    ave_psm["Self.Only"] = self_only_list
    ave_psm["N_Saw"] = n_saw_vals_list
    ave_psm["N_Exp_With_Prey"] = n_exp_with_prey_list

    if iterations:
        print(f"Iterations: {iterations}")
        # Calculate p-values for WD scores
        wd_pvals = calculate_wd_pvals_matrix(normalized_wd_scores, ave_psm_values, prey_means_matrix, prey_sd_matrix, self_interaction, n_exp_with_prey, n_saw_values, n_experiments, ave_psm, rows, columns, iterations, norm_val)

        # Create a vector to add to the dataframe
        wd_pvals_list = []
        for i, row in ave_psm.iterrows():
            row_index = np.where(rows == row["Prey"])[0][0]
            col_index = np.where(columns == row["Experiment.ID"])[0][0]
            wd_pvals_list.append(wd_pvals[row_index, col_index])

        # Perform a benjamini-hochberg correction on the p-values
        wd_pvals_bh = multipletests(wd_pvals_list, method='fdr_bh')[1] # Might want to remove this to avoid confusion

        ave_psm["WD_pval"] = wd_pvals_list
        ave_psm["WDFDR"] = wd_pvals_bh

    return ave_psm


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='FastCompPASS')

    parser.add_argument('--input', type=str, required=True, help='Input file')
    parser.add_argument('--norm_factor', type=str, default=0.98, help='Normalization factor')
    parser.add_argument('--iterations', type=int, default=None, help='Number of iterations for bootstrapping')

    args = parser.parse_args()

    # Get directory name for output files
    output_dir = '/'.join(args.input.split('/')[:-1])
    # print(output_dir)

    # Read input file
    try:
        input = pd.read_csv(args.input, sep='\t', index_col=0).astype({'Prey': str, 'Bait': str})
    except:
        print("Couldn't read input file with tab separator. Trying Commas.")
        input = pd.read_csv(args.input, sep=',', index_col=0).astype({'Prey': str, 'Bait': str})


    result = score_compPass(input, float(args.norm_factor), args.iterations)

    # Write output to file
    result.to_csv(output_dir + '/compPASS.csv', index=False)

if __name__ == '__main__':
    main()