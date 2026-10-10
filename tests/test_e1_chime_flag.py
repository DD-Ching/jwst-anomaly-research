"""scripts/e1_chime_flag.py: the rate-modulated null (offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

from jwst_anomaly import event_network as en

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("e1_chime_flag", _DIR / "e1_chime_flag.py")
cf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cf)


def _sample(n, rng):
    mjd = 59000 + rng.uniform(0, 200, n)
    mjd.sort()
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    return en.Sample("CHIME", mjd, rng.uniform(0, 360, n), dec, np.full(n, 0.2), np.full(n, 2020))


def test_rate_null_keeps_dec_hour_angle_and_optionally_time_of_day():
    rng = np.random.default_rng(1)
    s = _sample(500, rng)
    ha = np.mod(en.gmst_deg(s.mjd) - s.ra, 360)
    for keep_tod in (False, True):
        x = cf.rate_null(s, np.random.default_rng(2), 7, keep_tod)
        np.testing.assert_allclose(x.dec, s.dec)
        dha = np.mod(en.gmst_deg(x.mjd) - x.ra - ha + 180, 360) - 180
        np.testing.assert_allclose(dha, 0, atol=1e-6)
        assert x.mjd.min() >= np.floor(s.mjd.min()) and x.mjd.max() < np.floor(s.mjd.max()) + 1
        tod_same = np.allclose(np.sort(np.mod(x.mjd, 1)), np.sort(np.mod(s.mjd, 1)))
        assert tod_same is keep_tod


def test_count_uses_the_requested_lag_bin():
    p = en.Params()
    mjd = np.array([59000.0, 59000.0 + 600 / en.DAY, 59001.5])
    s = en.Sample(
        "CHIME", mjd, np.array([10.0, 100.0, 200.0]), np.zeros(3), np.full(3, 0.2), np.full(3, 2020)
    )
    assert cf.count(s, p, cf.CELLS["100s-1h"]) == 1  # 600 s, 90 deg apart
    assert cf.count(s, p, cf.CELLS["1h-1d"]) == 0  # 1.5 d and 1.49 d are > 1 d
