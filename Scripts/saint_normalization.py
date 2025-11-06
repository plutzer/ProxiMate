"""
Functions to extract SAINT normalization parameters and reverse normalization.

SAINT applies z-score normalization to log-transformed intensity values.
These functions help extract the normalization parameters and reverse the transformation.
"""

import pandas as pd
import numpy as np


def extract_saint_normalization_params(interaction_file, bait_file=None):
    """
    Extract normalization parameters used by SAINT from interaction and bait files.

    SAINT collects all non-NaN, valid intensity values from both test and control
    interactions, then applies z-score normalization using the global mean and std.

    Parameters:
    -----------
    interaction_file : str
        Path to SAINT interaction.txt file
        Format: IP_name, Bait_name, Prey_name, Spectral_count
    bait_file : str, optional
        Path to SAINT bait.txt file (to identify control vs test baits)
        Format: IP_name, Bait_name, T/C (test or control)
        If None, assumes all interactions should be included

    Returns:
    --------
    dict with keys:
        - 'mean': global mean of log-transformed intensities
        - 'std': global standard deviation of log-transformed intensities
        - 'n_values': number of values used in calculation
    """
    # Read interaction file
    # Typical format: IP, Bait, Prey, Intensity
    inter_df = pd.read_csv(interaction_file, sep='\t', header=None)

    # Assuming standard SAINT format: columns are IP, Bait, Prey, Intensity
    if inter_df.shape[1] == 4:
        inter_df.columns = ['IP', 'Bait', 'Prey', 'Intensity']
    else:
        raise ValueError(f"Expected 4 columns in interaction file, got {inter_df.shape[1]}")

    # Log-transform intensities (SAINT does this before normalization)
    # Filter out zero or negative values
    valid_intensities = inter_df['Intensity'][inter_df['Intensity'] > 0].copy()
    log_intensities = np.log(valid_intensities)

    # Collect all valid values (filter out -inf from log(0))
    # SAINT uses threshold -1e5 to filter out log(0)
    valid_log_intensities = log_intensities[log_intensities > -1e5]

    # Calculate global mean and std (same as SAINT does)
    global_mean = valid_log_intensities.mean()
    # SAINT uses var1 which is sample variance (divides by N, not N-1)
    # var1(x) = sum((x - mean)^2) / N
    # So std = sqrt(var1) = sqrt(sum((x - mean)^2) / N)
    global_std = valid_log_intensities.std(ddof=0)  # ddof=0 for population std (divides by N)

    return {
        'mean': global_mean,
        'std': global_std,
        'n_values': len(valid_log_intensities)
    }


def reverse_saint_normalization(normalized_intensity, mean, std):
    """
    Convert SAINT normalized intensity back to original intensity.

    SAINT normalization formula:
        normalized = (log(intensity) - mean) / std

    This function reverses both the z-score normalization and log transformation.

    Parameters:
    -----------
    normalized_intensity : float or array-like
        Normalized intensity value(s) from SAINT output
    mean : float
        Global mean used in normalization (from extract_saint_normalization_params)
    std : float
        Global std used in normalization (from extract_saint_normalization_params)

    Returns:
    --------
    float or array-like
        Original intensity value(s) before SAINT normalization
    """
    # Reverse z-score normalization
    log_intensity = (normalized_intensity * std) + mean

    # Reverse log transformation
    original_intensity = np.exp(log_intensity)

    return original_intensity


def normalize_like_saint(intensity, mean, std):
    """
    Apply SAINT-style normalization to intensity values.

    Useful for normalizing new data using previously extracted parameters.

    Parameters:
    -----------
    intensity : float or array-like
        Original intensity value(s)
    mean : float
        Global mean to use for normalization
    std : float
        Global std to use for normalization

    Returns:
    --------
    float or array-like
        Normalized intensity value(s)
    """
    # Log transform
    log_intensity = np.log(intensity)

    # Z-score normalization
    normalized = (log_intensity - mean) / std

    return normalized


# Example usage (modify these variables to test with your data):
if __name__ == "__main__":
    # Test files - modify these paths
    interaction_file = "filtered_interaction.txt"
    bait_file = "bait.dat"  # optional

    # Extract normalization parameters
    print("Extracting normalization parameters from SAINT input files...")
    params = extract_saint_normalization_params(interaction_file)

    print(f"\nNormalization parameters:")
    print(f"  Mean (log scale): {params['mean']:.6f}")
    print(f"  Std (log scale):  {params['std']:.6f}")
    print(f"  N values:         {params['n_values']}")

    # Test reverse normalization
    test_normalized_value = -0.088247
    print(f"\nTest reverse normalization:")
    print(f"  Normalized value: {test_normalized_value}")

    original = reverse_saint_normalization(
        test_normalized_value,
        params['mean'],
        params['std']
    )
    print(f"  Original intensity: {original:.2f}")

    # Verify by re-normalizing
    re_normalized = normalize_like_saint(
        original,
        params['mean'],
        params['std']
    )
    print(f"  Re-normalized value: {re_normalized:.6f}")
    print(f"  Match: {np.isclose(re_normalized, test_normalized_value)}")
