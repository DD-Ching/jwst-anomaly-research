from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import MaskedColumn, Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location("e1_swift", ROOT / "scripts" / "e1_swift.py")
sw = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sw)


def test_bat_events_selection_and_sigma():
    t = Table(
        {
            "name": ["GRB A", "GRB B", "GRB C", "GRB D"],
            "time": [60000.5, 60001.0, 60002.0, 60003.0],
            "trigger_obs": ["Swift", "Fermi", "Swift", "Swift "],
            "ra": [10.0, 20.0, 30.0, 40.0],
            "dec": [-5.0, 5.0, 15.0, 25.0],
            "error_radius": MaskedColumn(
                [2.146, 1.0, 1.0, 4.292], mask=[False, False, True, False]
            ),
        }
    )
    ev = sw.bat_events(t)
    # Fermi-triggered and missing-error rows dropped; whitespace in trigger_obs ignored
    assert list(ev["name"]) == ["GRB A", "GRB D"]
    assert np.allclose(ev["sigma"], [1 / 60, 2 / 60])
    assert set(ev["cat"]) == {"BAT"}
    assert ev.meta["provenance"] == "derived"


def test_control_rule():
    rng = np.random.default_rng(0)
    null = rng.poisson(2.0, 10000)
    ok, p = sw.control_passes(60, null)
    assert ok and p < 1e-20
    ok, _ = sw.control_passes(40, null)  # too few pairs, however significant
    assert not ok
    ok, p = sw.control_passes(3, null)
    assert not ok and p > 0.1


def test_outcome_rule():
    assert sw.outcome(np.array([0.5, 1e-3]), 35) == "null"
    assert sw.outcome(np.array([0.5, 1.3e-3 / 35]), 35) == "excess"
    assert sw.outcome(np.array([0.5, 1.4e-3 / 35]), 35) == "null"


def test_mask_family_has_35_cells_and_skips_control():
    m = sw.test_mask(sw.P)
    assert int(m.sum()) == 35
    assert not m[sw.CHANNELS.index(sw.CONTROL)].any()
    gw = sw.CHANNELS.index(("BAT", "GW"))
    assert m[gw, :, 2].sum() == 5 and not m[gw, :, :2].any()
