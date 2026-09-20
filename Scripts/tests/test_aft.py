"""Tests for the AFT modules behind --imputation.

`--imputation 0` routes to aft_impute_saint.filter_impute with impute=False, which drops
zero-intensity rows and writes filtered_interaction.txt -- the file SAINTexpress reads.
refactored_aft (--imputation 2) and one_component_aft (--imputation 3) share that tail
through interaction_filter, so the behavior tests run against aft_impute_saint and one
test pins that all three write identical bytes.
"""

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize
from scipy.stats import norm

import aft_impute_saint
import one_component_aft
import refactored_aft
from interaction_filter import write_filtered_interaction


MODULES = [aft_impute_saint, refactored_aft, one_component_aft]


@pytest.fixture
def saint_pipeline_inputs(tmp_path):
    """Write the interaction and ED files filter_impute reads.

    interaction.txt is the 4-column tab-separated SAINT format with no header, carrying
    positive, zero and negative intensities so the filter's boundary is observable.

    The returned output directory ends in a separator: filter_impute builds its output
    paths by string concatenation, and score.py passes ``f"{args.scoreInputs}/"``.
    """
    rows = [
        ("exp_1", "BaitA", "P1", 100.0),
        ("exp_1", "BaitA", "P2", 0.0),
        ("exp_1", "BaitA", "P3", 2500.5),
        ("exp_2", "BaitB", "P1", 0.0),
        ("exp_2", "BaitB", "P2", 75.25),
        ("exp_2", "BaitB", "P3", -3.0),
        ("exp_3", "BaitA", "P1", 640.0),
        ("exp_3", "BaitA", "P2", 0.0),
    ]
    interaction_path = tmp_path / "interaction.txt"
    with open(interaction_path, "w", newline="") as handle:
        for experiment, bait, prey, intensity in rows:
            handle.write(f"{experiment}\t{bait}\t{prey}\t{intensity}\n")

    ed_path = tmp_path / "ED.csv"
    pd.DataFrame([
        {"Experiment Name": "exp_1", "Type": "T", "Bait": "BaitA", "Replicate": 1,
         "Bait ID": "P9"},
        {"Experiment Name": "exp_2", "Type": "C", "Bait": "BaitB", "Replicate": 1,
         "Bait ID": "P8"},
        {"Experiment Name": "exp_3", "Type": "T", "Bait": "BaitA", "Replicate": 2,
         "Bait ID": "P9"},
    ]).to_csv(ed_path, index=False)

    out_dir = tmp_path / "out"
    out_dir.mkdir()

    return {
        "interaction_path": str(interaction_path),
        "ed_path": str(ed_path),
        "out_dir": str(out_dir) + "/",
        "out_path": out_dir,
        "kept": [r for r in rows if r[3] > 0],
    }


def _run(inputs):
    aft_impute_saint.filter_impute("prey.txt", inputs["interaction_path"],
                                   inputs["out_dir"], inputs["ed_path"], impute=False)
    return inputs["out_path"] / "filtered_interaction.txt"


# --- filtering ----------------------------------------------------------------

def test_keeps_only_positive_intensities_in_input_order(saint_pipeline_inputs):
    """The filter is strictly greater than zero, so zero and negative rows both go."""
    filtered = pd.read_csv(_run(saint_pipeline_inputs), sep="\t", header=None)
    expected = saint_pipeline_inputs["kept"]

    assert len(filtered) == len(expected) == 4
    assert (filtered[3] > 0).all()
    assert [tuple(r) for r in filtered.itertuples(index=False)] == expected


def test_writes_no_imputation_outputs(saint_pipeline_inputs):
    _run(saint_pipeline_inputs)
    written = {p.name for p in saint_pipeline_inputs["out_path"].iterdir()}

    assert written == {"filtered_interaction.txt"}


def test_all_three_implementations_agree(saint_pipeline_inputs, tmp_path):
    """Byte-identical filtered_interaction.txt whichever module writes it."""
    contents = []
    for module in MODULES:
        out_dir = tmp_path / f"out_{module.__name__}"
        out_dir.mkdir()
        module.filter_impute("prey.txt", saint_pipeline_inputs["interaction_path"],
                             str(out_dir) + "/", saint_pipeline_inputs["ed_path"],
                             impute=False)
        contents.append((out_dir / "filtered_interaction.txt").read_bytes())

    assert contents[0] == contents[1] == contents[2]


def _interaction_frame():
    return pd.DataFrame([
        {"ExperimentID": "exp_1", "Bait": "BaitA", "Prey": "P1", "Intensity": 100.0},
        {"ExperimentID": "exp_1", "Bait": "BaitA", "Prey": "P2", "Intensity": 0.0},
        {"ExperimentID": "exp_2", "Bait": "BaitB", "Prey": "P3", "Intensity": 250.0},
    ])


def test_helper_drops_working_columns_and_header(tmp_path):
    """The imputation paths attach a BaitID helper to the frame before writing.
    SAINTexpress parses this file positionally, so a fifth column or a header
    would corrupt it."""
    frame = _interaction_frame()
    frame["BaitID"] = ["IDA", "IDA", "IDB"]

    write_filtered_interaction(frame, str(tmp_path) + "/")
    lines = (tmp_path / "filtered_interaction.txt").read_text().strip().splitlines()

    assert all(len(line.split("\t")) == 4 for line in lines)
    assert not any("IDA" in line or "IDB" in line for line in lines)
    assert float(lines[0].split("\t")[3]) == 100.0


def test_helper_writes_columns_in_saint_order(tmp_path):
    """Column order follows the SAINT contract, not the input frame's order."""
    frame = _interaction_frame()[["Intensity", "Prey", "Bait", "ExperimentID"]]

    write_filtered_interaction(frame, str(tmp_path) + "/")
    first = (tmp_path / "filtered_interaction.txt").read_text().splitlines()[0]

    assert first.split("\t") == ["exp_1", "BaitA", "P1", "100.0"]


# --- censored likelihood ------------------------------------------------------

OBSERVED = np.array([5.1, 4.8, 5.4])
WITH_CENSORED = np.array([5.1, 4.8, 0.0, 5.4, 0.0])
MU, SIGMA, TLIM = 5.0, 0.4, 4.0


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


def test_pi_of_one_makes_censored_observations_free():
    """With pi = 1 every zero is attributed to the point mass, contributing log(1) = 0."""
    with_zeros = refactored_aft.protein_log_likelihood(
        WITH_CENSORED, MU, SIGMA, TLIM, 1.0)
    without_zeros = refactored_aft.protein_log_likelihood(
        OBSERVED, MU, SIGMA, TLIM, 1.0)

    assert with_zeros == pytest.approx(without_zeros)


def test_larger_pi_raises_the_likelihood_of_censored_data():
    """Zeros cost likelihood under the Gaussian alone and less as pi grows."""
    without_zeros = refactored_aft.protein_log_likelihood(OBSERVED, MU, SIGMA, TLIM, 0.0)
    none = refactored_aft.protein_log_likelihood(WITH_CENSORED, MU, SIGMA, TLIM, 0.0)
    low = refactored_aft.protein_log_likelihood(WITH_CENSORED, MU, SIGMA, TLIM, 0.1)
    high = refactored_aft.protein_log_likelihood(WITH_CENSORED, MU, SIGMA, TLIM, 0.9)

    assert none < without_zeros
    assert none < low < high


def test_neg_likelihood_is_minimized_at_the_sample_mean_and_population_std():
    """With no censored values the maximum-likelihood fit is the plain Gaussian MLE."""
    sample = np.array([4.2, 5.0, 5.3, 4.7, 5.8])
    fit = minimize(refactored_aft.neg_likelihood, x0=(4.0, 1.0),
                   args=(sample, 0.0, 0.0), method="Nelder-Mead",
                   options={"xatol": 1e-8, "fatol": 1e-8})

    assert fit.x[0] == pytest.approx(sample.mean(), abs=1e-4)
    assert abs(fit.x[1]) == pytest.approx(sample.std(ddof=0), abs=1e-4)


def test_one_component_floors_phi_to_stay_finite_far_above_tlim():
    """mu >> Tlim drives the normal CDF to zero; the 1e-12 clip keeps log() finite."""
    result = one_component_aft.protein_log_likelihood(
        np.array([100.0, 0.0]), mu=100.0, sigma=0.1, Tlim=0.0)

    assert np.isfinite(result)


@pytest.mark.parametrize("values, expected_sigma", [
    ([0.0, 100.0, 1000.0, 10000.0], np.log10([100.0, 1000.0, 10000.0]).std(ddof=0)),
    ([0.0, 500.0, 500.0], 0.5),
], ids=["spread", "identical-values-default-sigma"])
def test_get_initial_params_uses_log10_of_nonzero_values(values, expected_sigma):
    mu, sigma = refactored_aft.get_initial_params(pd.Series(values))
    logs = np.log10([v for v in values if v > 0])

    assert mu == pytest.approx(logs.mean())
    assert sigma == pytest.approx(expected_sigma)


# --- pi estimation --------------------------------------------------------------

def _control_design_and_interaction():
    """Control baits with three, four and one replicates.

    Missingness rises with lower intensity for bait C1 and is uniformly high for C2,
    so the two per-bait spline fits differ and the weighted average is observable.
    """
    rng = np.random.default_rng(7)
    ed_rows, rows = [], []
    for bait, reps, high_missing in [("C1", 3, False), ("C2", 4, True), ("C3", 1, False)]:
        for rep in range(1, reps + 1):
            exp = f"{bait}_r{rep}"
            ed_rows.append({"Experiment Name": exp, "Type": "C", "Bait": bait,
                            "Replicate": rep, "Bait ID": f"{bait}_ID"})
            for p in range(40):
                level = 2.0 + p * 0.1
                missing_prob = 0.8 if high_missing else max(0.0, 0.9 - p * 0.03)
                value = 0.0 if rng.random() < missing_prob else 10 ** (level + rng.normal(0, 0.1))
                rows.append((exp, bait, f"P{p}", value))
    ed_rows.append({"Experiment Name": "T_r1", "Type": "T", "Bait": "BaitT",
                    "Replicate": 1, "Bait ID": "T_ID"})
    rows += [("T_r1", "BaitT", f"P{p}", 1000.0) for p in range(40)]
    interaction = pd.DataFrame(rows, columns=["ExperimentID", "Bait", "Prey", "Intensity"])
    return pd.DataFrame(ed_rows), interaction


def test_estimate_pi_single_bait_and_weighted_average():
    """single_bait fits one control; weighted_average pools controls with enough
    replicates, weighting each fit by its replicate count and skipping the rest."""
    ed, interaction = _control_design_and_interaction()

    pi_c1 = refactored_aft.estimate_pi(interaction, ed, "single_bait", selected_bait="C1")
    pi_c2 = refactored_aft.estimate_pi(interaction, ed, "single_bait", selected_bait="C2")
    pooled = refactored_aft.estimate_pi(interaction, ed, "weighted_average")

    assert 0.0 <= pi_c1 <= 1.0 and 0.0 <= pi_c2 <= 1.0
    assert pi_c1 != pytest.approx(pi_c2)
    assert pooled == pytest.approx((3 * pi_c1 + 4 * pi_c2) / 7)


@pytest.mark.parametrize("kwargs", [
    {"method": "single_bait"},
    {"method": "single_bait", "selected_bait": "missing"},
    {"method": "single_bait", "selected_bait": "C3"},
    {"method": "weighted_average", "min_replicates": 5},
    {"method": "other"},
], ids=["no-bait", "unknown-bait", "too-few-replicates", "no-eligible-controls",
        "unknown-method"])
def test_estimate_pi_rejects_unusable_selections(kwargs):
    ed, interaction = _control_design_and_interaction()

    with pytest.raises(ValueError):
        refactored_aft.estimate_pi(interaction, ed, **kwargs)
