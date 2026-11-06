"""The GUI's downloadable example dataset parses under every quant type."""

import os

import pandas as pd
import pytest

import parse

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "..", "GUI", "example_data")


@pytest.mark.parametrize("quant_type", ["Intensity", "LFQ", "Spectral Counts"])
def test_example_dataset_parses(quant_type, tmp_path):
    experiments, controls = parse.parse_ed_pg(os.path.join(EXAMPLE, "proteinGroups.txt"),
                                              os.path.join(EXAMPLE, "ED.csv"), quant_type,
                                              str(tmp_path))
    assert (experiments, controls) == (16, 12)


def test_example_baits_have_bait_ids():
    ed = pd.read_csv(os.path.join(EXAMPLE, "ED.csv"))
    tests = ed[ed["Type"] == "T"]
    assert tests["Bait"].nunique() == 8
    assert tests["Bait ID"].notna().all()
