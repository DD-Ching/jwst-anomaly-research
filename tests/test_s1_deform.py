"""Offline tests for scripts/s1_deform.py (synthetic light curves only)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
_spec = importlib.util.spec_from_file_location("s1_deform", ROOT / "scripts" / "s1_deform.py")
sd = importlib.util.module_from_spec(_spec)
sys.modules["s1_deform"] = sd
_spec.loader.exec_module(sd)
from jwst_anomaly import burst_twins as bt  # noqa: E402


def _two_pulse(n=200, centres=(60, 110), widths=(5, 12), amps=(1.0, 0.6)):
    t = np.arange(n)
    f = sum(
        a * np.exp(-0.5 * ((t - c) / w) ** 2) for c, w, a in zip(centres, widths, amps, strict=True)
    )
    return np.vstack([f, 0.7 * f])


def test_band_factors_identity_and_direction():
    assert np.allclose(sd.band_factors(1.0, 200.0), 1.0)
    # a copy with energies divided by 2 (s = 2) moves photons from 300-1000 keV into 50-300 keV and
    # from 50-300 into 10-50: the band ratio changes by tens of per cent for E_peak = 100 keV
    f = sd.band_factors(2.0, 100.0)
    assert abs(f[0] / f[1] - 1) > 0.2


def test_stretched_copy_defeats_s1_statistic_but_not_grid():
    x = _two_pulse()
    y, _ = bt.stretch(x, np.ones_like(x), 2.0)
    r0, _ = bt.xcorr_max(x, y)
    # the M1 search stretches y by s from the grid; y was stretched by 2, so s = 0.5 undoes it
    yb, _ = bt.stretch(y, np.ones_like(y), 0.5)
    r1, _ = bt.xcorr_max(x, yb)
    assert r0 < 0.95 < r1


def test_rho_free_ignores_band_ratio_and_bounds_rho():
    x = _two_pulse()
    y = x * np.array([[1.4], [0.6]])  # same shape, different band ratio
    r_tot, _ = bt.xcorr_max(x, y)
    r_free, _ = sd.rho_free(x, y)
    assert r_free > 0.999 and r_free >= r_tot - 1e-12
