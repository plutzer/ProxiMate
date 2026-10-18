"""The shared threshold presets: their values, matching them to slider values, and
their acceptance by the backend and the help text."""

import backend
import help_text
import thresholds


def test_presets_hold_the_documented_values():
    assert thresholds.PRESETS == {
        'SAINT only': {'SaintScore': 0.9, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0},
        'WDFDR only': {'SaintScore': 0.0, 'BFDR': 1.0, 'WD': 0.0, 'WDFDR': 0.05},
        'Both': {'SaintScore': 0.9, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 0.05},
    }
    assert thresholds.DEFAULT is thresholds.PRESETS['SAINT only']


def test_active_preset_matches_values_or_returns_none():
    for name, values in thresholds.PRESETS.items():
        assert thresholds.active_preset(values) == name
        assert thresholds.active_preset({k: str(v) for k, v in values.items()}) == name
    assert thresholds.active_preset({'SaintScore': 0.7, 'BFDR': 0.05, 'WD': 0.0, 'WDFDR': 1.0}) is None
    assert thresholds.active_preset(thresholds.NO_FILTER) is None


def test_every_preset_and_no_filter_validate():
    for values in list(thresholds.PRESETS.values()) + [thresholds.NO_FILTER, thresholds.DEFAULT]:
        assert backend.validate_thresholds(values) == values


def test_describe_names_the_four_cutoffs():
    assert thresholds.describe(thresholds.DEFAULT) == 'SAINT ≥ 0.9, BFDR ≤ 0.05, WD ≥ 0, WDFDR ≤ 1'


def test_every_preset_button_has_a_tooltip_naming_its_values():
    for name, values in thresholds.PRESETS.items():
        tooltip = help_text.TOOLTIPS[thresholds.preset_id(name)]
        assert thresholds.describe(values) in tooltip
    assert {s[0] for s in thresholds.SLIDERS} <= set(help_text.TOOLTIPS)
