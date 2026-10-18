from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location("e1_icetracks", ROOT / "scripts" / "e1_icetracks.py")
it = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(it)


def _row(run, evt, sub, mjd, le, err, ra, dec):
    return [run, evt, sub, mjd, le, err, ra, dec, 0.0, 0.0]


def test_dr2_events_cuts_and_one_per_readout():
    ev = np.array(
        [
            _row(1, 7, 1, 55000.2, 4.5, 0.5, 10.0, 20.0),  # second subevent of readout (1, 7)
            _row(1, 7, 0, 55000.2, 4.2, 0.4, 11.0, 21.0),  # first subevent: kept
            _row(1, 8, 0, 55000.1, 3.9, 0.4, 12.0, 22.0),  # below the energy cut
            _row(2, 1, 0, 55001.0, 5.0, 0.3, 13.0, -6.0),  # south of the Dec cut
            _row(2, 2, 0, 54999.0, 4.0, 0.2, 14.0, -4.0),  # kept (edge of both cuts)
        ]
    )
    t = it.dr2_events(ev)
    assert list(t["name"]) == ["r2e2", "r1e7"]  # time-sorted
    assert np.allclose(t["sigma"], [0.2, 0.4])
    assert np.allclose(t["ra"], [14.0, 11.0])
    assert set(t["cat"]) == {"DR2"}


def test_control_rule():
    rng = np.random.default_rng(0)
    null = rng.poisson(0.3, 10000)
    ok, p = it.control_passes(40, null)
    assert ok and p <= it.CONTROL_MAX_P
    ok, _ = it.control_passes(it.CONTROL_MIN_PAIRS - 1, null)
    assert not ok
    ok, _ = it.control_passes(2, null)
    assert not ok


def test_outcome_and_mask():
    assert it.outcome(np.array([0.5, 1.3e-3 / 35]), 35) == "excess"
    assert it.outcome(np.array([0.5, 1.4e-3 / 35]), 35) == "null"
    m = it.test_mask(it.P)
    assert int(m.sum()) == 35
    assert not m[it.CHANNELS.index(it.CONTROL)].any()
