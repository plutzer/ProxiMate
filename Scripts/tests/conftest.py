"""Shared fixtures for the ProxiMate test suite.

Column convention for CompPASS input
------------------------------------
``ProteinGroups.write_CompPASS`` writes the header
``["Experiment.ID", "Replicate", "Bait", "Prey", "Prey.Name", "Spectral.Count"]``
but fills "Experiment.ID" with the bait *name* and "Bait" with the bait's *protein ID*.
``score_compPass`` relies on this: it flags a self-interaction by comparing the "Bait"
column against the "Prey" column, and both must therefore be protein IDs.  The fixtures
below follow the same convention, so a self-interaction is a row whose bait protein ID
equals its prey ID.
"""

import logging

import matplotlib

# Ann_Enrichment and refactored_aft import pyplot at module level.  conftest is imported
# before any test module, so this is the one place that reliably selects a backend
# needing no display.
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest


def _comppass_frame(counts):
    """Build a CompPASS input frame from {(bait, prey): [count_rep1, count_rep2, ...]}.

    Bait protein IDs are the bait name suffixed with "_ID", so no row is a
    self-interaction unless the prey is itself named "<bait>_ID".
    """
    rows = []
    for (bait, prey), replicate_counts in counts.items():
        for replicate, count in enumerate(replicate_counts, start=1):
            rows.append({
                "Experiment.ID": bait,
                "Replicate": replicate,
                "Bait": f"{bait}_ID",
                "Prey": prey,
                "Prey.Name": prey,
                "Spectral.Count": count,
            })
    return pd.DataFrame(rows)


@pytest.fixture
def comppass_input():
    """Three baits x three preys x two replicates, with no self-interactions.

    Prey P1 has AvePSM 11, 1, 1 across baits B1, B2, B3, which drives the
    hand-computed Mean/SD/Z/WD assertions in test_compass_scoring.py.
    """
    return _comppass_frame({
        ("B1", "P1"): [10, 12], ("B1", "P2"): [2, 0],   ("B1", "P3"): [1, 1],
        ("B2", "P1"): [1, 1],   ("B2", "P2"): [20, 18], ("B2", "P3"): [0, 2],
        ("B3", "P1"): [2, 0],   ("B3", "P2"): [1, 3],   ("B3", "P3"): [30, 30],
    })


@pytest.fixture
def comppass_input_with_self_interaction():
    """Four baits covering both self-interaction branches.

    Prey "B4_ID" is seen only in bait B4, whose protein ID is also "B4_ID" — the
    ``self_interaction_only`` edge case, where the prey mean and SD are overridden.

    Prey "B1_ID" is seen with itself as bait (B1) *and* with baits B2 and B3, so it is
    a self-interaction that is not self-only.  Its AvePSM of 20 in B1 must be excluded
    from the prey mean, which no other fixture exercises.
    """
    return _comppass_frame({
        ("B1", "P1"): [10, 12],   ("B1", "P2"): [2, 0],   ("B1", "B1_ID"): [18, 22],
        ("B2", "P1"): [1, 1],     ("B2", "P2"): [20, 18], ("B2", "B1_ID"): [2, 2],
        ("B3", "P1"): [2, 0],     ("B3", "P2"): [1, 3],   ("B3", "B1_ID"): [4, 4],
        ("B4", "P1"): [1, 1],     ("B4", "B4_ID"): [40, 44],
    })


@pytest.fixture
def comppass_input_large():
    """Six baits x eight preys, sized so two independent permutation runs cannot
    coincidentally agree on every p-value."""
    rng = np.random.default_rng(12345)
    counts = {}
    for b in range(1, 7):
        for p in range(1, 9):
            counts[(f"B{b}", f"P{p}")] = rng.integers(0, 30, size=2).tolist()
    return _comppass_frame(counts)


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def isolate_logging(monkeypatch):
    """Run every test against a pristine ProxiMate logger and restore it afterwards.

    Autouse rather than opt-in: any test that exercises a parse entry point
    attaches a file handler to the process-wide logger pointed at a ``tmp_path``
    that pytest later deletes.  Without restoration those handlers stay attached
    and every later test writes through them.

    The configuration is also torn down up front, and the environment variables
    that would otherwise leak in from whoever launched pytest are cleared, so a
    test sees the logger exactly as a fresh process would.

    Reaches into ``log_config``'s module state deliberately: the one-shot
    ``_initialized`` guard and the handler registry are what have to be undone.
    """
    import log_config

    package = logging.getLogger(log_config.PACKAGE)
    root = logging.getLogger()
    saved = (
        list(package.handlers), package.level, package.propagate,
        list(root.handlers), root.level, log_config._initialized,
        dict(log_config._file_handlers), dict(log_config._file_handler_refs),
    )

    package.handlers.clear()
    log_config._initialized = False
    log_config._file_handlers.clear()
    log_config._file_handler_refs.clear()
    for name in ("PROXIMATE_RUN_ID", "LOG_LEVEL", "PROXIMATE_LOG_DIR"):
        monkeypatch.delenv(name, raising=False)

    yield

    for handler in package.handlers:
        if handler not in saved[0]:
            handler.close()
    (package.handlers[:], package.level, package.propagate,
     root.handlers[:], root.level, log_config._initialized) = saved[:6]
    log_config._file_handlers.clear()
    log_config._file_handlers.update(saved[6])
    log_config._file_handler_refs.clear()
    log_config._file_handler_refs.update(saved[7])
