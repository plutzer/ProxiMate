"""Tests for the hypergeometric annotation-feature enrichment and its heatmap.

The heatmap is rendered by Shiny's ``render.plot``, which resizes the figure to the
browser card before drawing it, so the layout assertions re-run that resize and measure
the drawn artists rather than reading the source.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from Ann_Enrichment import (MAX_LABEL_CHARS, enrich_foreground, plot_results,
                            process_refactored, split_and_clean)


# --- split_and_clean -----------------------------------------------------------

def test_splits_on_semicolons_and_strips_whitespace():
    assert split_and_clean("Nucleus ; Cytoplasm;  Golgi ") == {
        "Nucleus", "Cytoplasm", "Golgi"}


@pytest.mark.parametrize("value", [np.nan, None, 3.5, ""])
def test_non_string_or_empty_input_yields_an_empty_set(value):
    assert split_and_clean(value) == set()


def test_trailing_digits_are_stripped():
    """UniProt numbers repeated features ("WD 1", "WD 2"); enrichment asks whether a
    protein has WD repeats at all, not which one."""
    assert split_and_clean("WD 1;WD 2;WD 3") == {"WD"}
    assert split_and_clean("Complex1") == {"Complex"}


def test_annotations_beginning_with_a_digit_are_dropped():
    assert split_and_clean("40S Ribosome;Nucleus") == {"Nucleus"}


# --- enrich_foreground ---------------------------------------------------------

@pytest.fixture
def feature_map():
    """Six of ten proteins carry "Nucleus"; two carry the rarer "Vesicle"."""
    mapping = {f"P{i:02d}": set() for i in range(1, 11)}
    for i in range(1, 7):
        mapping[f"P{i:02d}"] = {"Nucleus"}
    mapping["P07"] = {"Vesicle"}
    mapping["P08"] = {"Vesicle"}
    return mapping


def test_pvalue_and_enrichment_match_the_hypergeometric_definition(feature_map):
    all_ids = set(feature_map)                       # M = 10
    foreground = {"P01", "P02", "P03", "P07"}        # n = 4, Nucleus k = 3

    result = enrich_foreground(foreground, all_ids, feature_map).set_index("Feature")
    row = result.loc["Nucleus"]

    assert row["k"] == 3 and row["n"] == 4 and row["K"] == 6 and row["M"] == 10
    # Survival function at k-1 gives P(X >= k), the one-sided over-representation test.
    assert row["p_value"] == pytest.approx(hypergeom.sf(3 - 1, 10, 6, 4))
    assert row["enrichment"] == pytest.approx((3 / 4) / (6 / 10))


def test_rare_and_singly_seen_features_are_dropped(feature_map):
    """Vesicle has K = 2 < 5 so it never reaches the test even with k = 2; Nucleus clears
    K >= 5 but a single foreground hit fails the k >= 2 guard."""
    result = enrich_foreground({"P01", "P02", "P07", "P08"}, set(feature_map), feature_map)
    assert set(result["Feature"]) == {"Nucleus"}

    result = enrich_foreground({"P01"}, set(feature_map), feature_map)
    assert result.empty


# --- process_refactored --------------------------------------------------------

OPEN_THRESHOLDS = {"SaintScore": 0.0, "BFDR": 1.0, "WD": 0.0, "WDFDR": 1.0}


def _thresholds(**overrides):
    """Thresholds that admit everything, narrowed by whichever one a test is about."""
    return {**OPEN_THRESHOLDS, **overrides}


@pytest.fixture
def annotated_scores():
    """One row per (Experiment.ID, Prey.ID), as annotated_scores.csv carries.

    Preys P01-P06 are nuclear, P07-P08 vesicular, P09-P12 cytoplasmic.  In E1 the
    high-scoring preys are P01-P03 plus P07; in E2 they are P01, P02 and P04-P06.
    """
    localization = {}
    for i in range(1, 7):
        localization[f"P{i:02d}"] = "Nucleus"
    for i in (7, 8):
        localization[f"P{i:02d}"] = "Vesicle"
    for i in range(9, 13):
        localization[f"P{i:02d}"] = "Cytoplasm"

    passing = {"E1": {"P01", "P02", "P03", "P07"},
               "E2": {"P01", "P02", "P04", "P05", "P06"}}

    rows = []
    for experiment in ("E1", "E2"):
        for prey, scl in localization.items():
            rows.append({
                "Experiment.ID": experiment,
                "Prey.ID": prey,
                "SaintScore": 0.9 if prey in passing[experiment] else 0.1,
                "BFDR": 0.01,
                "WD": 5.0,
                "WDFDR": 0.01,
                "SCL": scl,
            })
    return pd.DataFrame(rows)


def test_threshold_selects_the_foreground(annotated_scores):
    results = process_refactored(annotated_scores, ["SCL"], _thresholds(SaintScore=0.7))
    e1 = results[results["Bait"] == "E1"].set_index("Feature")

    # E1's foreground is four preys, three of them nuclear, out of twelve overall.
    assert e1.loc["Nucleus", "n"] == 4
    assert e1.loc["Nucleus", "k"] == 3
    assert e1.loc["Nucleus", "K"] == 6
    assert e1.loc["Nucleus", "M"] == 12


def test_adjusted_pvalues_are_corrected_within_each_bait_and_feature_type(annotated_scores):
    """BH runs per (feature type, experiment) group, not across the whole result."""
    results = process_refactored(annotated_scores, ["SCL"], _thresholds(SaintScore=0.7))

    for (_, _), group in results.groupby(["Feature_type", "Bait"]):
        expected = multipletests(group["p_value"].tolist(), method="fdr_bh")[1]
        assert np.allclose(group["adj_p"].to_numpy(dtype=float), expected)


def test_result_columns_are_stable_even_when_nothing_passes(annotated_scores):
    """The GUI and its CSV export read these columns positionally."""
    populated = process_refactored(annotated_scores, ["SCL"], _thresholds(SaintScore=0.7))
    empty = process_refactored(annotated_scores, ["SCL"], _thresholds(SaintScore=1.5))

    assert list(populated.columns) == [
        "Bait", "Feature", "Feature_type", "k", "n", "K", "M",
        "p_value", "enrichment", "adj_p"]
    assert len(empty) == 0
    assert list(empty.columns) == list(populated.columns)


def test_count_columns_are_integers(annotated_scores):
    results = process_refactored(annotated_scores, ["SCL"], _thresholds(SaintScore=0.7))

    for column in ("k", "n", "K", "M"):
        assert results[column].dtype == np.int64, column


def test_every_score_narrows_the_foreground(annotated_scores):
    """The foreground is selected the same way the thresholding and download tabs
    filter, so a bait's enriched features answer to the same four scores."""
    wide = process_refactored(annotated_scores, ["SCL"], _thresholds())
    assert wide[wide["Bait"] == "E1"]["n"].iloc[0] == 12

    for narrowing in ({"BFDR": 0.001}, {"WD": 9.0}, {"WDFDR": 0.001}):
        results = process_refactored(annotated_scores, ["SCL"], _thresholds(**narrowing))
        assert results.empty, f"{narrowing} did not reach the foreground"


def test_conflicting_annotations_for_one_prey_are_an_error():
    """A prey's feature set is a property of the protein, not of the experiment it was
    seen in; two rows disagreeing about it is corrupt input, not a tie to break."""
    rows = []
    for experiment, scl_for_p01 in (("E1", "Nucleus"), ("E2", "Cytoplasm")):
        for i in range(1, 8):
            prey = f"P{i:02d}"
            rows.append({
                "Experiment.ID": experiment,
                "Prey.ID": prey,
                "SaintScore": 0.9,
                "BFDR": 0.01,
                "WD": 5.0,
                "WDFDR": 0.01,
                "SCL": scl_for_p01 if prey == "P01" else "Nucleus",
            })

    with pytest.raises(ValueError):
        process_refactored(pd.DataFrame(rows), ["SCL"], _thresholds(SaintScore=0.7))


# --- plot_results -------------------------------------------------------------

# Shiny renders at 96 ppi.  The card the heatmap sits in is around 8.75 inches wide on a
# maximized window.
PPI = 96.0
CONTAINER_SIZE = (8.75, 4.17)

FEATURES = [
    "COP9 signalosome",
    "Cul2-RING ubiquitin ligase complex",
    "glutamatergic synapse",
    "postsynaptic density",
    "Cul3-RING ubiquitin ligase complex",
    "lysosomal membrane",
    "cytoplasmic ribonucleoprotein granule",
    "nuclear speck",
    "extracellular exosome",
    "focal adhesion",
]
BAITS = ["CUL3-mT", "mT-CUL3", "KLHL12-mT", "mT-KEAP1"]


def _results(features=FEATURES):
    """A long-format enrichment table of the shape ``process_refactored`` returns."""
    rng = np.random.default_rng(0)
    rows = [(bait, feature, "GO_CC", 5, 50, 20, 2000, 0.001, rng.uniform(1, 40), 0.001)
            for bait in BAITS for feature in features]
    return pd.DataFrame(rows, columns=["Bait", "Feature", "Feature_type", "k", "n",
                                       "K", "M", "p_value", "enrichment", "adj_p"])


def _draw(grid, size):
    """Resize and draw the figure the way Shiny does, and return its renderer.

    ``render.plot`` substitutes a tight layout only when the figure carries no layout
    engine of its own, so applying one here mirrors it rather than pre-empting it.
    """
    figure = grid.figure
    figure.set_size_inches(*size)
    figure.set_dpi(PPI)
    if figure.get_layout_engine() is None:
        figure.set_layout_engine(layout="tight")
    figure.canvas.draw()
    return figure.canvas.get_renderer()


@pytest.fixture
def rendered():
    """Draw the heatmap at the size of a maximized window and clean up after."""
    grid = plot_results(_results(), "GO_CC", num_features=30)
    yield grid, _draw(grid, CONTAINER_SIZE)
    plt.close(grid.figure)


def _visible_text(figure):
    for ax in figure.axes:
        for text in [ax.xaxis.label, ax.yaxis.label, ax.title,
                     *ax.get_xticklabels(), *ax.get_yticklabels()]:
            if text.get_visible() and text.get_text():
                yield text


def _clipped(figure, renderer):
    """Text whose bounding box leaves the canvas, and so is drawn only in part."""
    canvas = figure.bbox
    out = []
    for text in _visible_text(figure):
        box = text.get_window_extent(renderer)
        if (box.x0 < canvas.x0 - 0.5 or box.x1 > canvas.x1 + 0.5
                or box.y0 < canvas.y0 - 0.5 or box.y1 > canvas.y1 + 0.5):
            out.append(text.get_text())
    return out


def _overlap(a, b):
    """Area shared by two bounding boxes, in square points."""
    width = min(a.x1, b.x1) - max(a.x0, b.x0)
    height = min(a.y1, b.y1) - max(a.y0, b.y0)
    return max(width, 0) * max(height, 0)


def test_the_layout_is_recomputed_when_the_figure_is_resized(rendered):
    """seaborn leaves behind an engine that does nothing, which would freeze the margins
    at the ones computed for the authored figure size."""
    grid, _ = rendered
    engine = grid.figure.get_layout_engine()

    assert engine is not None
    assert type(engine).__name__ != "PlaceHolderLayoutEngine"
    assert engine.adjust_compatible, "Shiny replaces an engine it cannot adjust and warns"


def test_no_label_is_clipped_by_the_figure_edge(rendered):
    grid, renderer = rendered
    assert _clipped(grid.figure, renderer) == []


def test_the_colorbar_does_not_cover_the_data(rendered):
    grid, _ = rendered
    cbar = grid.cax.get_window_extent()

    assert _overlap(cbar, grid.ax_heatmap.get_window_extent()) == 0
    assert _overlap(cbar, grid.ax_row_dendrogram.get_window_extent()) == 0
    assert _overlap(cbar, grid.ax_col_dendrogram.get_window_extent()) == 0


def test_every_selected_feature_gets_a_labelled_row(rendered):
    grid, _ = rendered

    assert set(grid.data2d.index) == set(FEATURES)
    assert list(grid.ax_heatmap.get_yticks()) == list(np.arange(len(FEATURES)) + 0.5)
    assert [t.get_text() for t in grid.ax_heatmap.get_yticklabels()] == \
        list(grid.data2d.index)


def test_features_are_selected_by_adjusted_p_enrichment_and_count():
    """A feature reaches the map only where adj_p <= 0.05 and enrichment >= 2 in some
    bait; the top-N ranks by how many baits it passes in."""
    results = _results()
    results.loc[results["Feature"] == "COP9 signalosome", "adj_p"] = 0.2
    results.loc[results["Feature"] == "nuclear speck", "enrichment"] = 1.5
    # "focal adhesion" passes in one bait only, so it ranks below every other feature.
    focal = results["Feature"] == "focal adhesion"
    results.loc[focal & (results["Bait"] != BAITS[0]), "adj_p"] = 0.5

    grid = plot_results(results, "GO_CC", num_features=30)
    everything = set(grid.data2d.index)
    plt.close(grid.figure)

    assert "COP9 signalosome" not in everything
    assert "nuclear speck" not in everything
    assert "focal adhesion" in everything

    grid = plot_results(results, "GO_CC", num_features=7)
    top = set(grid.data2d.index)
    plt.close(grid.figure)

    assert len(top) == 7
    assert "focal adhesion" not in top
    assert top < everything


def test_a_long_feature_name_is_truncated_with_an_ellipsis():
    """Drawn in full, one long name takes the width the map itself needs."""
    long_name = "positive regulation of transcription " * 4
    grid = plot_results(_results(FEATURES + [long_name]), "GO_CC", num_features=30)
    labels = [text.get_text() for text in grid.ax_heatmap.get_yticklabels()]
    plt.close(grid.figure)

    assert any(label.endswith("...") for label in labels)
    assert all(len(label) <= MAX_LABEL_CHARS for label in labels)
