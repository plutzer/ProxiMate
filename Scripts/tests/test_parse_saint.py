"""Tests for parsing SAINT-format inputs.

SAINT's ``bait.txt`` carries no bait protein IDs, which makes this entry point
the odd one out: the experiment counts it reports and the annotations that can
later be derived from it both differ from the other formats.
"""

import pandas as pd
import pytest

import parse
from experimental_design import ExperimentalDesign


@pytest.fixture
def saint_inputs(tmp_path):
    """Write a prey and interaction file for two test baits and one control.

    Six experiments in total, of which two are controls, so a count that
    includes controls is distinguishable from one that does not.
    """
    bait_df = pd.DataFrame({
        "Experiment Name": ["t1_1", "t1_2", "t2_1", "t2_2", "c_1", "c_2"],
        "Bait": ["BaitA", "BaitA", "BaitB", "BaitB", "Ctrl", "Ctrl"],
        "Type": ["T", "T", "T", "T", "C", "C"],
    })
    bait_df["Bait ID"] = "None"

    prey = tmp_path / "prey.txt"
    prey.write_text("P1\tGene1\nP2\tGene2\n")

    interaction = tmp_path / "interaction.txt"
    interaction.write_text("".join(
        f"{e}\t{b}\t{p}\t{v}\n"
        for e, b in zip(bait_df["Experiment Name"], bait_df["Bait"])
        for p, v in (("P1", 10), ("P2", 5))))

    return bait_df, str(prey), str(interaction)


def test_experiment_count_excludes_controls(saint_inputs, tmp_path):
    """`ExperimentalDesign.num_experiments` counts only non-control rows, and the
    other parse entry points return that.  Counting all six here would report a
    different quantity under the same name, in run.json and in the GUI table."""
    bait_df, prey, interaction = saint_inputs
    out = tmp_path / "out"

    n_expts, n_ctrls = parse.parse_from_saint(bait_df, prey, interaction, str(out))

    assert (n_expts, n_ctrls) == (4, 2)


def test_counts_match_the_experimental_design_it_writes(saint_inputs, tmp_path):
    """The reconstructed ED.csv is what scoring reads, so the counts reported at
    parse time must be the counts that file yields."""
    bait_df, prey, interaction = saint_inputs
    out = tmp_path / "out"

    n_expts, n_ctrls = parse.parse_from_saint(bait_df, prey, interaction, str(out))

    ed = ExperimentalDesign(str(out / "ED.csv"))
    assert (n_expts, n_ctrls) == (ed.num_experiments, ed.num_controls)


def test_absent_bait_ids_are_reported(saint_inputs, tmp_path, caplog):
    """bait.txt has no bait protein IDs, so anything keyed on the bait — BioGRID,
    self-interaction — cannot be derived.  Silence there looks like a negative
    result rather than a missing input."""
    bait_df, prey, interaction = saint_inputs
    out = tmp_path / "out"

    with caplog.at_level("WARNING", logger="proximate.parse"):
        parse.parse_from_saint(bait_df, prey, interaction, str(out))

    assert "bait" in caplog.text.lower()
    assert "BioGRID" in caplog.text


def test_real_bait_ids_are_not_reported(saint_inputs, tmp_path, caplog):
    bait_df, prey, interaction = saint_inputs
    bait_df = bait_df.assign(**{"Bait ID": ["Q1"] * 4 + ["Q2"] * 2})
    out = tmp_path / "out"

    with caplog.at_level("WARNING", logger="proximate.parse"):
        parse.parse_from_saint(bait_df, prey, interaction, str(out))

    assert "BioGRID" not in caplog.text
