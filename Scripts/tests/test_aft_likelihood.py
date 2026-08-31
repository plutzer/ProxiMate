"""Tests for the censored-likelihood functions behind AFT imputation.

Covers the two implementations the GUI exposes: --imputation 2 runs refactored_aft
(two-component) and 3 runs one_component_aft.
"""

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize
from scipy.stats import norm

import one_component_aft
import refactored_aft


OBSERVED = np.array([5.1, 4.8, 5.4])
WITH_CENSORED = np.array([5.1, 4.8, 0.0, 5.4, 0.0])
MU, SIGMA, TLIM = 5.0, 0.4, 4.0


# --- agreement between the two Gaussian implementations ------------------------

def test_refactored_with_zero_pi_equals_one_component():
    """refactored_aft's mixture reduces to the one-component form when pi = 0."""
    mixture = refactored_aft.protein_log_likelihood(WITH_CENSORED, MU, SIGMA, TLIM, 0.0)
    one_component = one_component_aft.protein_log_likelihood(
        WITH_CENSORED, MU, SIGMA, TLIM)

    assert mixture == pytest.approx(one_component)


def test_observed_term_is_a_gaussian_log_likelihood_without_the_constant():
    """Both drop the -0.5*log(2*pi) term, so they sit above scipy by n*log(sqrt(2*pi))."""
    expected = (norm.logpdf(OBSERVED, MU, SIGMA).sum()
                + len(OBSERVED) * np.log(np.sqrt(2 * np.pi)))

    assert refactored_aft.protein_log_likelihood(
        OBSERVED, MU, SIGMA, TLIM, 0.0) == pytest.approx(expected)
    assert one_component_aft.protein_log_likelihood(
        OBSERVED, MU, SIGMA, TLIM) == pytest.approx(expected)


# --- the pi mixture ------------------------------------------------------------

def test_pi_of_one_makes_censored_observations_free():
    """With pi = 1 every zero is attributed to the point mass, contributing log(1) = 0."""
    with_zeros = refactored_aft.protein_log_likelihood(
        WITH_CENSORED, MU, SIGMA, TLIM, 1.0)
    without_zeros = refactored_aft.protein_log_likelihood(
        OBSERVED, MU, SIGMA, TLIM, 1.0)

    assert with_zeros == pytest.approx(without_zeros)


def test_censored_observations_lower_the_likelihood_when_pi_is_zero():
    with_zeros = refactored_aft.protein_log_likelihood(
        WITH_CENSORED, MU, SIGMA, TLIM, 0.0)
    without_zeros = refactored_aft.protein_log_likelihood(
        OBSERVED, MU, SIGMA, TLIM, 0.0)

    assert with_zeros < without_zeros


def test_larger_pi_raises_the_likelihood_of_censored_data():
    low = refactored_aft.protein_log_likelihood(WITH_CENSORED, MU, SIGMA, TLIM, 0.1)
    high = refactored_aft.protein_log_likelihood(WITH_CENSORED, MU, SIGMA, TLIM, 0.9)

    assert high > low


# --- optimization behavior -----------------------------------------------------

def test_neg_likelihood_is_minimized_at_the_sample_mean_and_population_std():
    """With no censored values the maximum-likelihood fit is the plain Gaussian MLE."""
    sample = np.array([4.2, 5.0, 5.3, 4.7, 5.8])
    fit = minimize(refactored_aft.neg_likelihood, x0=(4.0, 1.0),
                   args=(sample, 0.0, 0.0), method="Nelder-Mead",
                   options={"xatol": 1e-8, "fatol": 1e-8})

    assert fit.x[0] == pytest.approx(sample.mean(), abs=1e-4)
    assert abs(fit.x[1]) == pytest.approx(sample.std(ddof=0), abs=1e-4)


# --- the one-component floor ---------------------------------------------------

def test_one_component_floors_phi_to_stay_finite_far_above_tlim():
    """mu >> Tlim drives the normal CDF to zero; the 1e-12 clip keeps log() finite."""
    result = one_component_aft.protein_log_likelihood(
        np.array([100.0, 0.0]), mu=100.0, sigma=0.1, Tlim=0.0)

    assert np.isfinite(result)


def test_refactored_without_the_floor_diverges_at_pi_zero():
    """The mixture has no clip, so pi = 0 in the same regime yields -inf."""
    with np.errstate(divide="ignore"):
        result = refactored_aft.protein_log_likelihood(
            np.array([100.0, 0.0]), 100.0, 0.1, 0.0, 0.0)

    assert result == -np.inf


# --- initial parameter estimates ------------------------------------------------

def test_get_initial_params_uses_log10_of_nonzero_values():
    intensities = pd.Series([0.0, 100.0, 1000.0, 10000.0])
    mu, sigma = refactored_aft.get_initial_params(intensities)
    logs = np.log10([100.0, 1000.0, 10000.0])

    assert mu == pytest.approx(logs.mean())
    assert sigma == pytest.approx(logs.std(ddof=0))


def test_get_initial_params_substitutes_a_default_sigma_when_values_are_identical():
    mu, sigma = refactored_aft.get_initial_params(pd.Series([0.0, 500.0, 500.0]))

    assert mu == pytest.approx(np.log10(500.0))
    assert sigma == 0.5
