import numpy as np
import pandas as pd
import argparse
import time
from statsmodels.stats.multitest import multipletests
from log_config import get_logger

logger = get_logger(__name__)

# Quantile level at which WD score matrices are normalized.  Observed and permuted
# matrices must share this level for the permutation test to be valid.
DEFAULT_NORM_FACTOR = 0.98

# Seed for the permutation null.  Fixed by default so repeated runs over the same input
# reproduce their WD p-values; any value works, the null only needs an arbitrary stream.
DEFAULT_SEED = 0

# Need a function for entropy calculation
def entropy(xs):
    # Convert input to numpy array for easier manipulation
    xs = np.array(xs, dtype=np.float64)
    # Calculate probabilities
    p = (xs + 1/len(xs)) / (np.sum(xs) + 1) #Pseudocount 1/len to avoid log(0)
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

def permute_prey_matrices(ave_psm_values, prey_means_matrix, prey_sd_matrix, n_saw_values, rng):
    # Returns a single permutation of prey values (counts), not labels
    # Permutes the prey values (spectral counts) within each prey (row), including zero values,
    # ensuring that the observed bait-prey structure and prey-specific variance are preserved.
    # Every matrix is reordered by one shared index vector per row, so the permuted WD
    # matrix holds the same multiset of values as the observed one.
    # rng is a numpy Generator; the caller owns it so a run can be reproduced.

    permuted_ave_psm_values = np.zeros_like(ave_psm_values)
    permuted_prey_means_matrix = np.zeros_like(prey_means_matrix)
    permuted_prey_sd_matrix = np.zeros_like(prey_sd_matrix)
    permuted_n_saw_values = np.zeros_like(n_saw_values)

    # For each prey (row), permute the values across experiments (columns)
    for i in range(ave_psm_values.shape[0]):
        # Get the spectral counts (values) for this prey across all experiments
        values = ave_psm_values[i, :]
        
        # Permute the values for this prey across experiments
        permuted_indices = rng.permutation(len(values))
        permuted_values = values[permuted_indices]
        permuted_means = prey_means_matrix[i, :][permuted_indices]
        permuted_sd = prey_sd_matrix[i, :][permuted_indices]
        permuted_n_saw = n_saw_values[i, :][permuted_indices]
        
        # Insert the permuted values back into the permuted matrices
        permuted_ave_psm_values[i, :] = permuted_values
        permuted_prey_means_matrix[i, :] = permuted_means
        permuted_prey_sd_matrix[i, :] = permuted_sd
        permuted_n_saw_values[i, :] = permuted_n_saw

    return permuted_ave_psm_values, permuted_prey_means_matrix, permuted_prey_sd_matrix, permuted_n_saw_values

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

def calculate_wd_pvals_matrix(normalized_wd_scores, ave_psm_values, prey_means_matrix, prey_sd_matrix, n_exp_with_prey, n_saw_values, n_experiments, iterations, norm_factor, rng):
    # norm_factor must be the level used for the observed scores.  Permutation preserves
    # each row's multiset of WD values, so both matrices share a quantile at any given
    # level and the two normalizations cancel; splitting the levels would rescale only the
    # observed side and bias every p-value.

    # Do things one permutation at a time to reduce memory usage
    p_values = np.zeros((normalized_wd_scores.shape[0], normalized_wd_scores.shape[1]))

    for i in range(iterations):
        permuted_ave_psm_values, permuted_prey_means_matrix, permuted_prey_sd_matrix, permuted_n_saw_values = permute_prey_matrices(ave_psm_values, prey_means_matrix, prey_sd_matrix, n_saw_values, rng)
        permuted_normalized_wd_scores, _ = calculate_wd_matrix(permuted_ave_psm_values, permuted_prey_means_matrix, permuted_prey_sd_matrix, n_exp_with_prey, permuted_n_saw_values, n_experiments, norm_factor=norm_factor)

        # Standard incrementation
        p_values += normalized_wd_scores < permuted_normalized_wd_scores

    p_values = p_values / iterations
    return p_values

def score_compPass(input, norm_factor=DEFAULT_NORM_FACTOR, iterations=None, seed=DEFAULT_SEED):
    """Score a CompPASS input table, optionally with permutation p-values for WD.

    norm_factor is the quantile level used to normalize the WD matrix.  It rescales the
    reported WD scores but leaves WD_pval and WDFDR unchanged, because the permutation
    null is normalized at the same level.

    seed fixes the permutation stream, so two calls on the same input return identical
    WD_pval and WDFDR.  Pass None for a fresh stream on every call.
    """

    rng = np.random.default_rng(seed)

    ave_psm = get_ave_psm(input)

    n_experiments = len(ave_psm['Experiment.ID'].unique())

    # Get the unique values for the columns and rows
    columns = ave_psm["Experiment.ID"].unique()
    rows = ave_psm["Prey"].unique()

    # Prey x experiment matrices; absent pairs are zero.  row_index/col_index give
    # every ave_psm row's position, used to read the matrices back out.
    row_index = pd.Index(rows).get_indexer(ave_psm["Prey"])
    col_index = pd.Index(columns).get_indexer(ave_psm["Experiment.ID"])

    def _matrix(values):
        return (pd.DataFrame({"Prey": ave_psm["Prey"], "Experiment.ID": ave_psm["Experiment.ID"],
                              "value": values})
                .pivot(index="Prey", columns="Experiment.ID", values="value")
                .reindex(index=rows, columns=columns).fillna(0.0).to_numpy(dtype=float))

    ave_psm_values = _matrix(ave_psm["AvePSM"])
    n_saw_values = _matrix(ave_psm["N_Saw"])
    # 1 where the bait and prey are the same protein, 0 otherwise
    self_interaction = _matrix((ave_psm["Bait"] == ave_psm["Prey"]).astype(float))

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
    normalized_wd_scores, _ = calculate_wd_matrix(ave_psm_values, prey_means_matrix, prey_sd_matrix, n_exp_with_prey, n_saw_values, n_experiments, norm_factor)

    # Read each ave_psm row's statistics back out of the matrices
    ave_psm["Mean"] = prey_means_matrix[row_index, col_index]
    ave_psm["SD"] = prey_sd_matrix[row_index, col_index]
    ave_psm["Z"] = z_scores[row_index, col_index]
    ave_psm["WD"] = normalized_wd_scores[row_index, col_index]
    ave_psm["Self.Interaction"] = self_interaction[row_index, col_index]
    ave_psm["Self.Only"] = np.isin(row_index, self_interaction_only)
    ave_psm["N_Saw"] = n_saw_values[row_index, col_index]
    ave_psm["N_Exp_With_Prey"] = n_exp_with_prey[row_index]

    # The p-value columns are always present: downstream filters read WDFDR from
    # every scored table and treat NaN as failing the cutoff.
    ave_psm["WD_pval"] = np.nan
    ave_psm["WDFDR"] = np.nan
    if iterations:
        wd_pvals = calculate_wd_pvals_matrix(normalized_wd_scores, ave_psm_values, prey_means_matrix, prey_sd_matrix, n_exp_with_prey, n_saw_values, n_experiments, iterations, norm_factor, rng)
        wd_pvals_list = wd_pvals[row_index, col_index]

        # Perform a benjamini-hochberg correction on the p-values
        wd_pvals_bh = multipletests(wd_pvals_list, method='fdr_bh')[1] # Might want to remove this to avoid confusion

        ave_psm["WD_pval"] = wd_pvals_list
        ave_psm["WDFDR"] = wd_pvals_bh

    return ave_psm


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description='FastCompPASS')

    parser.add_argument('--input', type=str, required=True, help='Input file')
    parser.add_argument('--norm_factor', type=float, default=DEFAULT_NORM_FACTOR, help='Quantile level for WD score normalization')
    parser.add_argument('--iterations', type=int, default=None, help='Number of iterations for bootstrapping')
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED, help='Seed for the permutation null')

    args = parser.parse_args()

    # Get directory name for output files
    output_dir = '/'.join(args.input.split('/')[:-1])

    # The input is tab- or comma-separated; the header line says which.
    with open(args.input) as handle:
        sep = '\t' if '\t' in handle.readline() else ','
    input = pd.read_csv(args.input, sep=sep, index_col=0).astype({'Prey': str, 'Bait': str})


    result = score_compPass(input, args.norm_factor, args.iterations, seed=args.seed)

    # Write output to file
    result.to_csv(output_dir + '/compPASS.csv', index=False)

if __name__ == '__main__':
    main()