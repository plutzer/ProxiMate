"""The four score cutoffs every tab reads from the sidebar, and the presets that set them."""

KEYS = ('SaintScore', 'BFDR', 'WD', 'WDFDR')

# The sliders in display order: (id suffix, threshold key, label, maximum, step).  The
# suffix is also the tooltip key.  Every slider runs from 0.
SLIDERS = (
    ('saintscore', 'SaintScore', 'SAINT ≥', 1.0, 0.01),
    ('bfdr', 'BFDR', 'BFDR ≤', 1.0, 0.01),
    ('wd', 'WD', 'WD ≥', 10.0, 0.1),
    ('wdfdr', 'WDFDR', 'WDFDR ≤', 1.0, 0.01),
)

# WD stays at 0 (off) in every preset: the sets are defined on SAINT/BFDR and WDFDR alone.
PRESETS = {
    'SAINT only': {'SaintScore': 0.9, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0},
    'WDFDR only': {'SaintScore': 0.0, 'BFDR': 1.0, 'WD': 0.0, 'WDFDR': 0.05},
    'Both': {'SaintScore': 0.9, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 0.05},
}

# Where the sliders start.  A dataset scored with 0 WDFDR iterations has no WDFDR
# values, so a WDFDR cut would empty it.
DEFAULT = PRESETS['SAINT only']

# Passes every row.
NO_FILTER = {'SaintScore': 0.0, 'BFDR': 1.0, 'WD': 0.0, 'WDFDR': 1.0}


def preset_id(name):
    """The input id of a preset's button, also its tooltip key: 'SAINT only' -> 'preset_saint_only'."""
    return 'preset_' + name.lower().replace(' ', '_')


def active_preset(values):
    """The preset whose four values equal ``values``, or None."""
    for name, preset in PRESETS.items():
        if all(abs(float(values[key]) - preset[key]) < 1e-9 for key in KEYS):
            return name
    return None


def describe(values):
    """One line naming the four cutoffs: 'SAINT ≥ 0.9, BFDR ≤ 0.05, WD ≥ 0, WDFDR ≤ 1'."""
    return ', '.join(f"{label} {float(values[key]):g}" for _, key, label, _, _ in SLIDERS)
