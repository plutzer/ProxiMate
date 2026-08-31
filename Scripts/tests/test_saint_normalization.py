"""Tests for extracting and reversing SAINT's z-score normalization of log intensities."""

import numpy as np
import pytest

from saint_normalization import (
    extract_saint_normalization_params,
    normalize_like_saint,
    reverse_saint_normalization,
)


def test_params_match_natural_log_of_positive_intensities(saint_interaction_file):
    path, positives = saint_interaction_file
    params = extract_saint_normalization_params(str(path))

    assert params["mean"] == pytest.approx(np.log(positives).mean())
    assert params["n_values"] == len(positives)


def test_std_is_the_population_std(saint_interaction_file):
    """SAINT's var1 divides by N, so the std must use ddof=0, not pandas' default ddof=1."""
    path, positives = saint_interaction_file
    params = extract_saint_normalization_params(str(path))
    log_intensities = np.log(positives)

    assert params["std"] == pytest.approx(log_intensities.std(ddof=0))
    # Pin the distinction: with only 3 values the two conventions differ substantially.
    assert params["std"] != pytest.approx(log_intensities.std(ddof=1))


def test_zero_and_negative_intensities_are_excluded(saint_interaction_file):
    path, positives = saint_interaction_file
    params = extract_saint_normalization_params(str(path))
    # The fixture writes five rows; only three carry positive intensities.
    assert params["n_values"] == 3 == len(positives)


def test_normalization_round_trips():
    intensities = np.array([100.0, 5000.0, 12345.6])
    mean, std = 7.0, 1.3

    normalized = normalize_like_saint(intensities, mean, std)
    assert np.allclose(reverse_saint_normalization(normalized, mean, std), intensities)


def test_normalize_applies_z_score_to_the_log():
    assert normalize_like_saint(np.e**3, 1.0, 2.0) == pytest.approx((3.0 - 1.0) / 2.0)


def test_reverse_is_exp_of_the_unscaled_value():
    assert reverse_saint_normalization(1.0, 1.0, 2.0) == pytest.approx(np.exp(3.0))


def test_wrong_column_count_raises(tmp_path):
    path = tmp_path / "interaction.txt"
    path.write_text("exp\tbait\tprey\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Expected 4 columns"):
        extract_saint_normalization_params(str(path))
