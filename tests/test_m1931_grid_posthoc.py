import importlib.util
import sys
from pathlib import Path

import numpy as np

_p = Path(__file__).resolve().parents[1] / "scripts" / "m1931_grid_posthoc.py"
_spec = importlib.util.spec_from_file_location("m1931_grid_posthoc", _p)
mod = importlib.util.module_from_spec(_spec)
sys.modules["m1931_grid_posthoc"] = mod
_spec.loader.exec_module(mod)


def test_peak_fraction_single_pixel_vs_psf():
    hot = np.zeros((5, 5))
    hot[2, 2] = 10.0
    assert mod.peak_fraction(hot, 0.0) == 1.0
    yy, xx = np.mgrid[-2:3, -2:3]
    psf = np.exp(-(xx**2 + yy**2) / (2 * 1.0**2))
    assert mod.peak_fraction(psf + 3.0, 3.0) < 0.25
    assert mod.peak_fraction(np.full((5, 5), np.nan), 0.0) is None
    assert mod.peak_fraction(np.zeros((5, 5)), 0.0) is None


def test_stamp_at_rejects_edges():
    img = np.arange(100.0).reshape(10, 10)
    assert mod.stamp_at(img, 2, 2).shape == (5, 5)
    assert mod.stamp_at(img, 1, 5) is None
    assert mod.stamp_at(img, 5, 8) is None


def test_n_exposures_counts_bits_and_bounds():
    con = np.zeros((2, 4, 4), dtype=np.int32)
    con[0, 1, 2] = 0b1011
    con[1, 1, 2] = 0b1
    con[0, 0, 0] = -1  # all 32 bits set
    assert mod.n_exposures(con, 2, 1) == 4
    assert mod.n_exposures(con, 0, 0) == 32
    assert mod.n_exposures(con, -1, 0) is None
    assert mod.n_exposures(con, 0, 4) is None
    assert mod.n_exposures(con[0], 2, 1) == 3
