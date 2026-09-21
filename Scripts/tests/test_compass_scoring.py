"""Tests for the CompPASS WD score, its components, and its permutation p-values."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import entropy as scipy_entropy

from compPASS_pval import (
    calculate_wd_matrix,
    entropy,
    get_ave_psm,
    normalize_matrix,
    score_compPass,
)


# --- components ------------------------------------------------------------------

def test_entropy_matches_an_independent_implementation():
    """The (x + 1/n)/(sum + 1) pseudocount yields a proper probability vector."""
    counts = np.array([3.0, 1.0, 0.0, 5.0])
    probabilities = (counts + 1 / len(counts)) / (counts.sum() + 1)

    assert probabilities.sum() == pytest.approx(1.0)
    assert entropy(counts) == pytest.approx(scipy_entropy(probabilities, base=2))


def test_normalize_matrix_takes_the_quantile_over_nonzero_entries_only():
    """Zeros are excluded, so the median of [2, 4, 8] is 4.0 rather than 3.0."""
    matrix = np.array([[0.0, 2.0], [4.0, 8.0]])
    normalized, q = normalize_matrix(matrix, 0.5)

    assert q == pytest.approx(4.0)
    assert np.allclose(normalized, matrix / 4.0)


def test_wd_matrix_matches_the_hand_computed_formula():
    """WD = sqrt(AvePSM * ((sd/mean) * (N_experiments/N_exp_with_prey)) ** N_saw)."""
    ave_psm = np.array([[4.0, 9.0], [16.0, 25.0]])
    means = np.array([[2.0, 2.0], [5.0, 5.0]])
    sds = np.array([[1.0, 1.0], [2.0, 2.0]])
    n_saw = np.array([[1.0, 2.0], [2.0, 1.0]])
    n_exp_with_prey = np.array([2.0, 1.0])

    wd, q = calculate_wd_matrix(ave_psm, means, sds, n_exp_with_prey, n_saw,
                                n_experiments=2, norm_factor=None)

    # Row 0: inner = (1/2) * (2/2) = 0.5.  Row 1: inner = (2/5) * (2/1) = 0.8.
    expected = np.array([
        [np.sqrt(4.0 * 0.5**1), np.sqrt(9.0 * 0.5**2)],
        [np.sqrt(16.0 * 0.8**2), np.sqrt(25.0 * 0.8**1)],
    ])
    assert q is None
    assert np.allclose(wd, expected)


def test_ave_psm_averages_replicates_and_counts_nonzero_ones(comppass_input):
    ave_psm = get_ave_psm(comppass_input).set_index(["Experiment.ID", "Prey"])

    # B1/P1 has replicate counts 10 and 12.
    assert ave_psm.loc[("B1", "P1"), "AvePSM"] == pytest.approx(11.0)
    assert ave_psm.loc[("B1", "P1"), "N_Saw"] == 2
    # B1/P2 has counts 2 and 0, so only one replicate saw it.
    assert ave_psm.loc[("B1", "P2"), "AvePSM"] == pytest.approx(1.0)
    assert ave_psm.loc[("B1", "P2"), "N_Saw"] == 1


def test_ave_psm_takes_the_max_over_duplicate_replicate_rows(comppass_input):
    """Duplicate (Experiment.ID, Prey, Replicate) rows collapse to their maximum."""
    duplicated = comppass_input.copy()
    extra = duplicated[
        (duplicated["Experiment.ID"] == "B1")
        & (duplicated["Prey"] == "P1")
        & (duplicated["Replicate"] == 1)
    ].copy()
    extra["Spectral.Count"] = 4  # lower than the existing 10, so it must be ignored
    duplicated = pd.concat([duplicated, extra], ignore_index=True)

    ave_psm = get_ave_psm(duplicated).set_index(["Experiment.ID", "Prey"])
    assert ave_psm.loc[("B1", "P1"), "AvePSM"] == pytest.approx(11.0)


# --- score_compPass ------------------------------------------------------------

def test_score_compass_statistics_match_hand_computation(comppass_input):
    """Prey P1 has AvePSM 11, 1, 1 across the three baits.

    Mean divides the prey's total by the number of experiments; SD uses the sample
    convention (n_experiments - 1); Z and WD follow from those.  No bait's protein ID
    is a prey here, so nothing is a self-interaction.
    """
    scored = score_compPass(comppass_input, norm_factor=None)
    row = scored.set_index(["Experiment.ID", "Prey"]).loc[("B1", "P1")]

    assert row["Mean"] == pytest.approx(4.333333, abs=1e-6)
    assert row["SD"] == pytest.approx(5.773503, abs=1e-6)
    assert row["Z"] == pytest.approx(1.154701, abs=1e-6)
    assert row["WD"] == pytest.approx(4.418894, abs=1e-6)
    assert not scored["Self.Interaction"].any()
    assert not scored["Self.Only"].any()


def test_prey_seen_only_with_itself_as_bait_is_flagged_self_only(
        comppass_input_with_self_interaction):
    """B4_ID appears solely in bait B4, whose bait protein ID is also B4_ID."""
    scored = score_compPass(comppass_input_with_self_interaction,
                            norm_factor=None).set_index(["Experiment.ID", "Prey"])
    row = scored.loc[("B4", "B4_ID")]

    assert bool(row["Self.Interaction"]) is True
    assert bool(row["Self.Only"]) is True
    # The self-only correction pins the prey's experiment count at 1.
    assert row["N_Exp_With_Prey"] == 1


def test_self_interaction_is_excluded_from_the_prey_mean(
        comppass_input_with_self_interaction):
    """Prey B1_ID has AvePSM 20 (with itself as bait), 2, 4 and 0 across four baits.

    The self-interaction is dropped from the numerator while the denominator stays at
    the full experiment count: (26 - 20) / 4 = 1.5.
    """
    scored = score_compPass(comppass_input_with_self_interaction,
                            norm_factor=None).set_index(["Experiment.ID", "Prey"])

    assert scored.loc[("B1", "B1_ID"), "Mean"] == pytest.approx(1.5)
    assert scored.loc[("B1", "B1_ID"), "SD"] == pytest.approx(1.707825, abs=1e-6)
    # Tiled across experiments, so every row for this prey carries the same mean.
    assert scored.loc[("B2", "B1_ID"), "Mean"] == pytest.approx(1.5)

    assert bool(scored.loc[("B1", "B1_ID"), "Self.Interaction"]) is True
    assert bool(scored.loc[("B2", "B1_ID"), "Self.Interaction"]) is False
    # Seen with other baits too, so the self-only override must not apply.
    assert bool(scored.loc[("B1", "B1_ID"), "Self.Only"]) is False
    assert scored.loc[("B1", "B1_ID"), "N_Exp_With_Prey"] == 3


def test_normalization_rescales_wd_scores_by_one_quantile(comppass_input_large):
    """norm_factor is not inert: a single quantile divides every WD score."""
    unnormalized = score_compPass(comppass_input_large.copy(), norm_factor=None)["WD"]
    at_098 = score_compPass(comppass_input_large.copy(), norm_factor=0.98)["WD"]
    at_050 = score_compPass(comppass_input_large.copy(), norm_factor=0.50)["WD"]

    ratios = (unnormalized.to_numpy() / at_098.to_numpy())
    finite = ratios[np.isfinite(ratios)]
    assert np.allclose(finite, finite[0])
    assert finite[0] > 0
    assert not np.allclose(at_098.to_numpy(), at_050.to_numpy())


# --- the permutation null ------------------------------------------------------

@pytest.mark.parametrize("seed_kwargs", [{}, {"seed": 7}], ids=["default-seed", "explicit"])
def test_wd_pvalues_are_reproducible_across_runs(comppass_input_large, seed_kwargs):
    """Two runs over identical input agree exactly, so WD p-values are reproducible."""
    first = score_compPass(comppass_input_large.copy(), 0.98, iterations=200, **seed_kwargs)
    second = score_compPass(comppass_input_large.copy(), 0.98, iterations=200, **seed_kwargs)

    assert first["WD_pval"].between(0.0, 1.0).all()
    assert np.array_equal(first["WD_pval"].to_numpy(), second["WD_pval"].to_numpy())
    assert np.array_equal(first["WDFDR"].to_numpy(), second["WDFDR"].to_numpy())


def test_different_seeds_draw_different_nulls(comppass_input_large):
    """Guards the reproducibility test: a frozen permutation would satisfy it too."""
    seed_one = score_compPass(comppass_input_large.copy(), 0.98, iterations=200, seed=1)
    seed_two = score_compPass(comppass_input_large.copy(), 0.98, iterations=200, seed=2)

    assert not np.array_equal(seed_one["WD_pval"].to_numpy(),
                              seed_two["WD_pval"].to_numpy())


def test_wd_pvalues_do_not_depend_on_the_normalization_factor(comppass_input_large):
    """Permutation preserves each prey row's multiset of WD values, so observed and null
    matrices share a quantile at any level and the two normalizations cancel exactly."""
    at_098 = score_compPass(comppass_input_large.copy(), 0.98, iterations=100)["WD_pval"]
    at_050 = score_compPass(comppass_input_large.copy(), 0.50, iterations=100)["WD_pval"]
    at_020 = score_compPass(comppass_input_large.copy(), 0.20, iterations=100)["WD_pval"]

    assert np.array_equal(at_098.to_numpy(), at_050.to_numpy())
    assert np.array_equal(at_098.to_numpy(), at_020.to_numpy())


def test_zero_iterations_still_writes_the_pvalue_columns_as_nan(comppass_input):
    """Downstream filters read WDFDR from every scored table; without permutations
    the column exists and is all NaN, which the filters treat as failing."""
    result = score_compPass(comppass_input.copy(), 0.98, iterations=0)
    assert result["WD_pval"].isna().all()
    assert result["WDFDR"].isna().all()
