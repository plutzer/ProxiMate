"""Tests for the global BFDR recomputation applied to merged per-group SAINT output."""

import numpy as np
import pandas as pd
import pytest

from bfdr_pool import recompute_bfdr


def _naive_bfdr(scores):
    """Direct O(n^2) transcription of the BFDR definition, used as an oracle."""
    out = []
    for value in scores:
        greater = scores[scores > value]
        if len(greater) == 0:
            out.append(0.0)
        else:
            out.append(min(max(1.0 - greater.sum() / len(greater), 0.0), 1.0))
    return np.array(out)


def test_matches_naive_reference():
    rng = np.random.default_rng(1)
    scores = np.round(rng.random(60), 2)
    result = recompute_bfdr(pd.DataFrame({"AvgP": scores, "BFDR": 0.0}))
    assert np.allclose(result["BFDR"].to_numpy(), _naive_bfdr(scores))


def test_tied_scores_share_a_bfdr_and_do_not_count_each_other():
    df = pd.DataFrame({"AvgP": [0.9, 0.9, 0.5, 0.5, 0.1], "BFDR": 0.0})
    bfdr = recompute_bfdr(df)["BFDR"].to_numpy()

    # Nothing scores above 0.9, so the highest-scoring rows get zero.
    assert bfdr[0] == bfdr[1] == 0.0
    # Above 0.5 sit the two 0.9s: 1 - 1.8/2.  The other 0.5 must not contribute.
    assert bfdr[2] == bfdr[3] == pytest.approx(0.1)
    # Above 0.1 sit 0.9, 0.9, 0.5, 0.5: 1 - 2.8/4.
    assert bfdr[4] == pytest.approx(0.3)


def test_preserves_input_row_order():
    scores = np.array([0.3, 0.95, 0.1, 0.7])
    df = pd.DataFrame({"AvgP": scores, "Prey": list("abcd"), "BFDR": 0.0})
    result = recompute_bfdr(df)

    assert list(result["Prey"]) == list("abcd")
    assert np.allclose(result["AvgP"].to_numpy(), scores)
    assert np.allclose(result["BFDR"].to_numpy(), _naive_bfdr(scores))
