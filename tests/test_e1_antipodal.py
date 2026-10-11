"""scripts/e1_antipodal.py: antipodal class, counts and injection (offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

from jwst_anomaly import event_network as en

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("e1_antipodal", _DIR / "e1_antipodal.py")
ap = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ap)


def _s(mjd, ra, dec, sigma=0.5):
    n = len(mjd)
    return en.Sample(
        "X",
        np.asarray(mjd, float),
        np.asarray(ra, float),
        np.asarray(dec, float),
        np.full(n, sigma),
        np.full(n, 2020),
    )


def test_antipodal_threshold_uses_floor_or_errors():
    assert ap.antipodal(np.array([171.0]), 0.5, 0.5)[0]
    assert not ap.antipodal(np.array([169.0]), 0.5, 0.5)[0]
    assert ap.antipodal(np.array([165.0]), 4.0, 4.0)[0]  # 3 * hypot(4, 4) = 17 deg


def test_count_channel_bins_antipodal_pairs_by_lag():
    d = en.DAY
    s = _s([59000, 59000 + 5 / d, 59000 + 600 / d, 59010], [10, 190, 10, 190], [20, -20, 20, -20])
    c = ap.count_channel(s, s, True)
    # (0,1): 5 s antipodal; (1,2): 595 s antipodal; (0,2): same direction; event 3 is > 7 d away
    assert c.tolist() == [1, 0, 1, 0, 0]


def test_injected_antipodal_pairs_are_counted():
    rng = np.random.default_rng(3)
    n = 300
    s = _s(
        59000 + rng.uniform(0, 100, n),
        rng.uniform(0, 360, n),
        np.degrees(np.arcsin(rng.uniform(-1, 1, n))),
    )
    base = ap.count_channel(s, s, True)[3]
    sb = ap.inject_antipodal(s, s, 20, 3600.0, en.DAY, True, rng)
    assert ap.count_channel(sb, sb, True)[3] - base >= 15
