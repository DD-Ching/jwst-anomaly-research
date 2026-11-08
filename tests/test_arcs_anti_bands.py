"""Frozen rules of ``scripts/arcs_anti_bands.py`` (results/arcs/anti-bands/PLAN.md)."""

import importlib.util
from pathlib import Path

import numpy as np

_p = Path(__file__).resolve().parents[1] / "scripts" / "arcs_anti_bands.py"
_spec = importlib.util.spec_from_file_location("arcs_anti_bands", _p)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)


def test_band_gate():
    assert m.band_gate(10, 8)
    assert not m.band_gate(9, 9)  # too few rows
    assert not m.band_gate(20, 15)  # 0.75 < 0.8


def test_band_verdict():
    assert m.band_verdict(True, 15.0) == "agree"
    assert m.band_verdict(True, -14.9) == "agree"
    assert m.band_verdict(True, 15.1) == "differs"
    assert m.band_verdict(False, 2.0) == "unusable"
    assert m.band_verdict(True, np.nan) == "unusable"


def test_row_label_uses_passing_bands_only():
    ok = {"SW": True, "LW": True}
    assert m.row_label({"SW": "agree", "LW": "differs"}, ok) == "blend_flag"
    assert m.row_label({"SW": "agree", "LW": "unusable"}, ok) == "same_axis"
    assert m.row_label({"SW": "unusable", "LW": "unusable"}, ok) == "undetermined"
    assert m.row_label({"SW": "agree", "LW": "differs"}, {"SW": True, "LW": False}) == "same_axis"
    assert m.row_label({"SW": "agree", "LW": "agree"}, {"SW": False, "LW": False}) == (
        "undetermined"
    )


def test_outcome():
    assert m.outcome(["blend_flag", "blend_flag", "same_axis"]) == (
        "open part of the pooled anti excess explained by colour blends or colour structure"
        " (exploratory)"
    )
    assert m.outcome(["same_axis", "same_axis", "undetermined"]).startswith("blend veto passed")
    assert m.outcome(["same_axis", "same_axis", "blend_flag"]) == "mixed / undetermined"
    assert m.outcome(["undetermined"] * 3) == "mixed / undetermined"


def test_row_index_refuses_missing_or_duplicate():
    import pytest
    from astropy.table import Table

    t = Table({"label": [1, 2, 2]})
    assert m.row_index(t, 1) == 0
    for bad in (2, 3):
        with pytest.raises(SystemExit):
            m.row_index(t, bad)
