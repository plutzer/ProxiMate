"""Layout tests for the feature enrichment heatmap.

The figure is rendered by Shiny's ``render.plot``, which resizes it to whatever the
browser card happens to be before drawing it.  A figure that only looks right at the size
it was authored at is therefore not enough, so the assertions here re-run that resize —
at several widths — before measuring anything.

Everything is measured off the drawn figure rather than read off the source, because the
failure this guards against, a label pushed past the canvas edge or a colorbar laid over
the data, exists only once the artists have been placed.
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from Ann_Enrichment import MAX_LABEL_CHARS, plot_results


# Shiny renders at 96 ppi.  The card the heatmap sits in is around 8.75 inches wide on a
# maximized window; the narrower sizes are the same card on a smaller one.
PPI = 96.0
CONTAINER_SIZES = [(8.75, 4.17), (8.75, 6.8), (5.2, 4.17), (4.2, 3.2)]

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
    yield grid, _draw(grid, CONTAINER_SIZES[0])
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


# --- the figure has to survive Shiny's resize -------------------------------------

def test_the_layout_is_recomputed_when_the_figure_is_resized(rendered):
    """seaborn leaves behind an engine that does nothing, which would freeze the margins
    at the ones computed for the authored figure size."""
    grid, _ = rendered
    engine = grid.figure.get_layout_engine()

    assert engine is not None
    assert type(engine).__name__ != "PlaceHolderLayoutEngine"
    assert engine.adjust_compatible, "Shiny replaces an engine it cannot adjust and warns"


@pytest.mark.parametrize("size", CONTAINER_SIZES)
def test_no_label_is_clipped_by_the_figure_edge(size):
    """The feature names are the point of the plot; a name running off the canvas is the
    failure this module exists for."""
    grid = plot_results(_results(), "GO_CC", num_features=30)
    renderer = _draw(grid, size)

    clipped = _clipped(grid.figure, renderer)
    plt.close(grid.figure)

    assert clipped == []


def test_a_full_map_of_long_names_still_fits():
    """Thirty features of the length GO cellular-component terms actually reach."""
    features = [f"positive regulation of transcription from RNA polymerase II {i}"
                for i in range(30)]
    grid = plot_results(_results(features), "GO_CC", num_features=30)
    renderer = _draw(grid, CONTAINER_SIZES[0])

    clipped = _clipped(grid.figure, renderer)
    plt.close(grid.figure)

    assert clipped == []


# --- the colorbar is placed rather than left where it lands -----------------------

def test_the_colorbar_does_not_cover_the_data(rendered):
    """Left as an inset it keeps a fixed fraction of the figure, which at the card's
    aspect lands across the dendrogram and the first rows of the map."""
    grid, _ = rendered
    cbar = grid.cax.get_window_extent()

    assert _overlap(cbar, grid.ax_heatmap.get_window_extent()) == 0
    assert _overlap(cbar, grid.ax_row_dendrogram.get_window_extent()) == 0
    assert _overlap(cbar, grid.ax_col_dendrogram.get_window_extent()) == 0


def test_the_colorbar_says_what_it_measures(rendered):
    grid, _ = rendered

    assert grid.cax.get_title(loc="left") == "log2 enrichment"


# --- the axis labels belong to the heatmap ----------------------------------------

def test_the_axes_are_labelled_on_the_heatmap(rendered):
    """Named on the heatmap rather than through ``plt``, which labels whichever axes
    happens to be current."""
    grid, _ = rendered

    assert grid.ax_heatmap.get_xlabel() == "Bait"
    assert grid.ax_heatmap.get_ylabel() == "Feature"


def test_the_colorbar_is_not_labelled_as_an_axis_of_the_data(rendered):
    grid, _ = rendered

    assert grid.cax.get_xlabel() == ""
    assert grid.cax.get_ylabel() == ""


# --- rows, columns and their labels -----------------------------------------------

def test_every_selected_feature_gets_a_labelled_row(rendered):
    grid, _ = rendered

    assert set(grid.data2d.index) == set(FEATURES)
    assert list(grid.ax_heatmap.get_yticks()) == list(np.arange(len(FEATURES)) + 0.5)
    assert [t.get_text() for t in grid.ax_heatmap.get_yticklabels()] == \
        list(grid.data2d.index)


def test_a_long_feature_name_is_truncated_with_an_ellipsis():
    """Drawn in full, one long name takes the width the map itself needs."""
    long_name = "positive regulation of transcription " * 4
    grid = plot_results(_results(FEATURES + [long_name]), "GO_CC", num_features=30)
    labels = [text.get_text() for text in grid.ax_heatmap.get_yticklabels()]
    plt.close(grid.figure)

    assert any(label.endswith("...") for label in labels)
    assert all(len(label) <= MAX_LABEL_CHARS for label in labels)
