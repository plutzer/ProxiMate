"""Tests for GUI/help_text.py: the tooltip string registry and tip() helper."""

import inspect
import re

import pytest

import help_text
from help_text import TOOLTIPS, tip

# Keys app.py references. Renaming or removing an entry must fail here first.
REQUIRED_KEYS = [
    # Shared score thresholds (used on 4 tabs)
    "saintscore", "bfdr", "wd", "wdfdr",
    # Threshold preset buttons (same values on every tab)
    "preset_stringent", "preset_moderate", "preset_relaxed", "preset_none",
    # Network Scoring
    "dataset_name", "input_format", "quant_type",
    "pg_file", "diann_matrix_file", "fragpipe_file", "msstats_file",
    "ed_file", "ed_file_msstats",
    "saint_bait", "saint_prey", "saint_interaction",
    "organism", "imputation_method", "pi_method", "pi_bait",
    "wdfdr_iterations",
    "clear_datasets", "download_session", "session_file",
    # Data Thresholding
    "pca_imputation", "pca_normalization", "pca_min_detection",
    "experiment_pca", "prey_pca", "prey_pca_color", "qc_bait",
    "metric_network_size", "metric_enrichment", "metric_degree",
    # Protein Feature Analysis
    "feature_analysis", "feature_type", "num_features",
    "download_pvalue_threshold", "download_enrichment_threshold",
    # Network Comparison
    "comp_bait", "volcano_plot", "venn", "gene_lists",
    # Downloads
    "dl_preset", "dl_groups", "custom_columns",
    "dl_genelist_mode", "dl_prohits_abundance",
    "dl_filter_card", "batch_export",
]


def test_required_keys_present():
    missing = [k for k in REQUIRED_KEYS if k not in TOOLTIPS]
    assert not missing, f"TOOLTIPS is missing keys referenced by app.py: {missing}"


def test_values_are_concise_nonempty_strings():
    for key, text in TOOLTIPS.items():
        assert isinstance(text, str), f"{key}: not a string"
        assert text.strip() == text, f"{key}: has leading/trailing whitespace"
        assert text, f"{key}: empty"
        assert len(text) <= 250, f"{key}: {len(text)} chars, tooltips must stay concise"


def test_keys_are_snake_case():
    for key in TOOLTIPS:
        assert re.fullmatch(r"[a-z0-9_]+", key), f"bad key name: {key!r}"


def test_every_entry_carries_the_tooltip_marker():
    # Each string entry in the source ends with a "# tooltip" comment so the
    # text can be found and edited by searching for that marker.
    source = inspect.getsource(help_text)
    assert source.count("# tooltip") >= len(TOOLTIPS)


def test_tip_renders_label_icon_and_text():
    html = str(tip("BFDR (≤)", "bfdr"))
    assert "BFDR" in html
    assert "ⓘ" in html  # the info icon
    assert TOOLTIPS["bfdr"] in html


def test_tip_unknown_key_raises():
    with pytest.raises(KeyError):
        tip("x", "no_such_key")
