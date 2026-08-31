"""Tests for the non-imputing path of filter_impute.

`--imputation 0` routes to aft_impute_saint.filter_impute with impute=False, which drops
zero-intensity rows and writes filtered_interaction.txt -- the file SAINTexpress reads.
refactored_aft and one_component_aft repeat the same filter-and-write tail, so every test
here runs against all three and one test pins that their outputs agree byte for byte.
"""

import os

import pandas as pd
import pytest

import aft_impute_saint
import one_component_aft
import refactored_aft
from interaction_filter import write_filtered_interaction


MODULES = [aft_impute_saint, refactored_aft, one_component_aft]
MODULE_IDS = ["aft_impute_saint", "refactored_aft", "one_component_aft"]


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
        "rows": rows,
        "kept": [r for r in rows if r[3] > 0],
    }


def _run(module, inputs, prey_path="prey.txt"):
    module.filter_impute(prey_path, inputs["interaction_path"], inputs["out_dir"],
                         inputs["ed_path"], impute=False)
    return inputs["out_path"] / "filtered_interaction.txt"


def _read_filtered(path):
    return pd.read_csv(path, sep="\t", header=None,
                       names=["ExperimentID", "Bait", "Prey", "Intensity"])


@pytest.mark.parametrize("module", MODULES, ids=MODULE_IDS)
def test_writes_filtered_interaction_file(module, saint_pipeline_inputs):
    output = _run(module, saint_pipeline_inputs)
    assert output.exists()


@pytest.mark.parametrize("module", MODULES, ids=MODULE_IDS)
def test_keeps_only_positive_intensities(module, saint_pipeline_inputs):
    """The filter is strictly greater than zero, so zero and negative rows both go."""
    filtered = _read_filtered(_run(module, saint_pipeline_inputs))

    assert len(filtered) == len(saint_pipeline_inputs["kept"]) == 4
    assert (filtered["Intensity"] > 0).all()
    assert 0.0 not in set(filtered["Intensity"])
    assert -3.0 not in set(filtered["Intensity"])


@pytest.mark.parametrize("module", MODULES, ids=MODULE_IDS)
def test_preserves_surviving_rows_in_order(module, saint_pipeline_inputs):
    filtered = _read_filtered(_run(module, saint_pipeline_inputs))
    expected = saint_pipeline_inputs["kept"]

    assert list(filtered["ExperimentID"]) == [r[0] for r in expected]
    assert list(filtered["Bait"]) == [r[1] for r in expected]
    assert list(filtered["Prey"]) == [r[2] for r in expected]
    assert list(filtered["Intensity"]) == [r[3] for r in expected]


@pytest.mark.parametrize("module", MODULES, ids=MODULE_IDS)
def test_output_is_four_columns_without_a_header(module, saint_pipeline_inputs):
    """SAINTexpress reads this file positionally, so an extra column or a header breaks it."""
    output = _run(module, saint_pipeline_inputs)
    lines = output.read_text().strip().splitlines()

    assert all(len(line.split("\t")) == 4 for line in lines)
    # A header would put a non-numeric value in the intensity field of the first row.
    assert float(lines[0].split("\t")[3]) == 100.0


@pytest.mark.parametrize("module", MODULES, ids=MODULE_IDS)
def test_writes_no_imputation_outputs(module, saint_pipeline_inputs):
    _run(module, saint_pipeline_inputs)
    written = {p.name for p in saint_pipeline_inputs["out_path"].iterdir()}

    assert written == {"filtered_interaction.txt"}


@pytest.mark.parametrize("module", MODULES, ids=MODULE_IDS)
def test_prey_file_is_not_read(module, saint_pipeline_inputs):
    """prey.txt is an imputation-only input; the filter path never opens it."""
    missing = os.path.join(saint_pipeline_inputs["out_path"], "does_not_exist.txt")
    assert not os.path.exists(missing)

    output = _run(module, saint_pipeline_inputs, prey_path=missing)
    assert output.exists()


def test_all_three_implementations_agree(saint_pipeline_inputs, tmp_path):
    """Byte-identical output across the three modules, so collapsing the duplicated
    filter-and-write tail into one helper is a safe refactor."""
    contents = []
    for name, module in zip(MODULE_IDS, MODULES):
        out_dir = tmp_path / f"out_{name}"
        out_dir.mkdir()
        module.filter_impute("prey.txt", saint_pipeline_inputs["interaction_path"],
                             str(out_dir) + "/", saint_pipeline_inputs["ed_path"],
                             impute=False)
        contents.append((out_dir / "filtered_interaction.txt").read_bytes())

    assert contents[0] == contents[1] == contents[2]


@pytest.mark.parametrize("module", MODULES, ids=MODULE_IDS)
def test_ed_file_must_carry_bait_id_even_without_imputation(module, saint_pipeline_inputs,
                                                            tmp_path):
    """All three build a bait name -> bait ID dict before checking `impute`, so the
    filter path requires a column it never uses."""
    ed_without_bait_id = tmp_path / "ED_no_bait_id.csv"
    pd.read_csv(saint_pipeline_inputs["ed_path"]).drop(columns=["Bait ID"]).to_csv(
        ed_without_bait_id, index=False)

    with pytest.raises(KeyError, match="Bait ID"):
        module.filter_impute("prey.txt", saint_pipeline_inputs["interaction_path"],
                             saint_pipeline_inputs["out_dir"], str(ed_without_bait_id),
                             impute=False)


# --- write_filtered_interaction, directly ---------------------------------------

def _interaction_frame():
    return pd.DataFrame([
        {"ExperimentID": "exp_1", "Bait": "BaitA", "Prey": "P1", "Intensity": 100.0},
        {"ExperimentID": "exp_1", "Bait": "BaitA", "Prey": "P2", "Intensity": 0.0},
        {"ExperimentID": "exp_2", "Bait": "BaitB", "Prey": "P3", "Intensity": 250.0},
    ])


def test_helper_drops_working_columns(tmp_path):
    """The imputation paths attach a BaitID helper to the frame before writing.
    SAINTexpress parses this file positionally, so a fifth column would corrupt it."""
    frame = _interaction_frame()
    frame["BaitID"] = ["IDA", "IDA", "IDB"]

    write_filtered_interaction(frame, str(tmp_path) + "/")
    lines = (tmp_path / "filtered_interaction.txt").read_text().strip().splitlines()

    assert all(len(line.split("\t")) == 4 for line in lines)
    assert not any("IDA" in line or "IDB" in line for line in lines)


def test_helper_writes_columns_in_saint_order(tmp_path):
    """Column order follows the SAINT contract, not the input frame's order."""
    frame = _interaction_frame()[["Intensity", "Prey", "Bait", "ExperimentID"]]

    write_filtered_interaction(frame, str(tmp_path) + "/")
    first = (tmp_path / "filtered_interaction.txt").read_text().splitlines()[0]

    assert first.split("\t") == ["exp_1", "BaitA", "P1", "100.0"]


def test_helper_reports_kept_and_total_counts(tmp_path):
    kept, total = write_filtered_interaction(_interaction_frame(), str(tmp_path) + "/")

    assert (kept, total) == (2, 3)
