"""Offline tests for scripts/goal_wormhole_timeline.py (carried-wormhole event accounting)."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location(
    "goal_wormhole_timeline", _DIR / "goal_wormhole_timeline.py"
)
gw = importlib.util.module_from_spec(_spec)
sys.modules["goal_wormhole_timeline"] = gw  # dataclasses look the module up
_spec.loader.exec_module(gw)


def test_coast_trip_is_special_relativity():
    tr = gw.coast_trip(20.0, 0.99)
    assert tr.earth_years == pytest.approx(20.0 / 0.99)
    assert tr.proper_years == pytest.approx(20.0 / 0.99 * math.sqrt(1 - 0.99**2))


def test_accel_trip_limits():
    # short trip: Newtonian t = 2 sqrt(D / a), proper time ~ coordinate time
    tr = gw.accel_trip(1e-6, 1.0)
    assert tr.earth_years == pytest.approx(2 * math.sqrt(1e-6), rel=1e-4)
    assert tr.proper_years == pytest.approx(tr.earth_years, rel=1e-4)
    # long trip: Earth time -> D + 2/a (asymptotically light-like); proper time ~ logarithmic
    tr = gw.accel_trip(1e4, 1.0)
    assert tr.earth_years == pytest.approx(1e4 + 2.0, rel=1e-3)
    assert tr.proper_years == pytest.approx(2 * math.log(1e4), rel=0.05)


@pytest.mark.parametrize("d", [1.0, 20.0, 100.0])
def test_one_carried_mouth_never_closes_a_loop_with_light(d):
    for beta in np.linspace(0.01, 0.999999, 400):
        tr = gw.coast_trip(d, float(beta))
        assert gw.loop_margin_years(tr.shift_years, d)["light_plus_throat"] > 0
    assert gw.loop_margin_years(gw.accel_trip(d).shift_years, d)["light_plus_throat"] > 0


def test_a_second_connection_with_another_shift_closes_a_loop():
    tr = gw.coast_trip(20.0, 0.9)
    m = gw.loop_margin_years(tr.shift_years, 20.0, other_shift_years=0.0)
    assert m["two_throats"] < 0
    assert gw.loop_margin_years(5.0, 20.0, other_shift_years=5.0)["two_throats"] == 0.0


def test_throat_maps_invert_and_reproduce_the_mother_s_times():
    tr = gw.coast_trip(20.0, 0.999)
    t_dest = tr.earth_years + 3.0
    assert gw.earth_to_throat(gw.throat_to_earth(t_dest, tr), tr) == pytest.approx(t_dest)
    ev = gw.timeline(tr, stay_years=1.0, ship_return_beta=0.99)
    assert ev["arrive"]["mother_elapsed_via_throat"] == pytest.approx(tr.proper_years)
    assert ev["home_via_throat"]["mother_elapsed"] == pytest.approx(tr.proper_years + 1.0)
    # going home by ship instead costs the mother more than two light-crossing times
    assert ev["home_via_ship"]["mother_elapsed"] > 2 * 20.0
    with pytest.raises(ValueError):
        gw.throat_to_earth(tr.earth_years - 1.0, tr)
