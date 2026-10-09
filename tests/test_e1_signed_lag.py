"""scripts/e1_signed_lag.py: signed counts, Skellam tail, one-sided injection (offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

from jwst_anomaly import event_network as en

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("e1_signed_lag", _DIR / "e1_signed_lag.py")
sl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sl)


def _s(cat, mjd, ra, dec, sigma=0.5):
    n = len(mjd)
    return en.Sample(
        cat,
        np.asarray(mjd, float),
        np.asarray(ra, float),
        np.asarray(dec, float),
        np.full(n, sigma),
        np.full(n, 2020),
    )


def test_signed_counts_sign_lag_and_class():
    d = en.DAY
    a = _s("A", [59000.0], [10.0], [0.0])
    b = _s("B", [59000 + 5 / d, 59000 - 600 / d, 59000 + 600 / d], [10.0, 100.0, 100.0], [0.0] * 3)
    c = sl.signed_counts(a, b)
    assert c[0].tolist() == [1, 0, 1]  # 5 s after, same direction
    assert c[2].tolist() == [0, 0, 0]  # 600 s before and after, wide: cancel


def test_skellam_tail_matches_a_single_pair_and_gaussian_for_large_cells():
    rng = np.random.default_rng(0)
    small = rng.poisson(0.03, (4000, 1)) - rng.poisson(0.03, (4000, 1))
    p, _ = sl.two_sided_p(np.array([-1.0]), small)
    assert 0.03 < p[0] < 0.09  # about 2 * 0.03
    big = rng.normal(0, 20, (4000, 1))
    p, z = sl.two_sided_p(np.array([40.0]), big)
    assert abs(z[0] - 2) < 0.1 and abs(p[0] - 0.0455) < 0.01


def test_injection_puts_b_after_a():
    rng = np.random.default_rng(1)
    n = 200
    a = _s("A", 59000 + rng.uniform(0, 50, n), rng.uniform(0, 360, n), rng.uniform(-60, 60, n))
    b = _s("B", 59000 + rng.uniform(0, 50, n), rng.uniform(0, 360, n), rng.uniform(-60, 60, n))
    base = sl.signed_counts(a, b)[3, 1]
    sb = sl.inject_after(a, b, 20, 3600.0, en.DAY, 1, rng)
    assert sl.signed_counts(a, sb)[3, 1] - base >= 10
