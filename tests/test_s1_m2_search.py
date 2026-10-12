"""Offline tests for scripts/s1_m2_search.py (synthetic light curves only)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
_spec = importlib.util.spec_from_file_location("s1_m2_search", ROOT / "scripts" / "s1_m2_search.py")
ms = importlib.util.module_from_spec(_spec)
sys.modules["s1_m2_search"] = ms
_spec.loader.exec_module(ms)
from jwst_anomaly import burst_twins as bt  # noqa: E402


def _pair(s0=1.3, n=160, noise=0.01, seed=0):
    t = np.arange(n)
    f = np.exp(-0.5 * ((t - 40) / 4) ** 2) + 0.6 * np.exp(-0.5 * ((t - 75) / 9) ** 2)
    x = np.vstack([f, 0.7 * f]) * 100
    y, _ = bt.stretch(x, np.ones_like(x), s0)
    rng = np.random.default_rng(seed)
    ex, ey = np.full_like(x, noise * 100), np.full_like(y, noise * 100)
    return x + rng.normal(size=x.shape) * ex, ex, y + rng.normal(size=y.shape) * ey, ey


def test_fine_grid_is_log_uniform_and_symmetric():
    g = np.array(ms.fine_grid())
    assert len(g) == 29 and abs(g[14] - 1.0) < 1e-12
    assert np.allclose(np.diff(np.log(g)), np.log(1.05), atol=2e-4)


def test_refinement_rescues_off_grid_true_pair():
    # a bright copy stretched by 1.3 (1/1.3 = 0.769 lies between grid values 0.7462 and 0.7835)
    x, ex, y, ey = _pair(1.3)
    m, s_grid, p_ref = ms.m2_only(x, ex, 0, y, ey, 0, with_chi2=True)
    ys, eys = bt.stretch(y, ey, s_grid)
    _, lag = ms.sd.rho_free(x, ys)
    p_grid = ms.sd.chi2_free(x, ex, ys, eys, lag)
    assert m > 0.97 and p_ref >= 1e-3 and p_ref > p_grid


def test_m2_skips_resolution_mismatch():
    x, ex, y, ey = _pair(1.0)
    assert np.isnan(ms.m2_only(x, ex, 0, y, ey, 3)[0])
