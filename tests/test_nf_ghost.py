"""Offline tests for the Neutrino Frontier ghost-pair helpers (E-NF1, NF-H04)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly import neutrino_frontier as nf  # noqa: E402

GP = nf.GhostParams()


def _sample(mjd, ra, dec, sigma=1.0):
    mjd = np.asarray(mjd, float)
    return en.Sample(
        "ICECAT",
        mjd,
        np.asarray(ra, float),
        np.asarray(dec, float),
        np.full(len(mjd), sigma),
        en.mjd_year(mjd),
    )


def test_wide_pair_stats_bins_and_weights():
    # pairs: (0,1) 2 h apart, 90 deg -> wide 1h-1d; (0,2) 40 d apart, 1 deg -> same (not wide);
    # (1,2) ~40 d apart, wide -> 30d-180d
    s = _sample([60000.0, 60000.0 + 2 / 24, 60040.0], [0.0, 90.0, 1.0], [0.0, 0.0, 0.0])
    w = np.array([0.5, 0.4, 0.2])
    out = nf.wide_pair_stats(s, w, GP)
    assert out.shape == (7, 2)
    assert out[3, 0] == 1 and out[3, 1] == pytest.approx(0.2)
    assert out[6, 0] == 1 and out[6, 1] == pytest.approx(0.08)
    assert out[:, 0].sum() == 2


def test_scramble_cyclic_stays_in_span_and_keeps_hour_angle():
    rng = np.random.default_rng(0)
    mjd = np.sort(rng.uniform(55000, 60000, 200))
    s = _sample(mjd, rng.uniform(0, 360, 200), rng.uniform(-60, 60, 200))
    sc = nf.scramble_cyclic(s, np.random.default_rng(1), 365.0)
    assert sc.mjd.min() >= mjd.min() and sc.mjd.max() <= mjd.max()
    ha0 = np.mod(s.ra - en.gmst_deg(s.mjd), 360)
    ha1 = np.mod(sc.ra - en.gmst_deg(sc.mjd), 360)
    assert np.allclose(np.mod(ha1 - ha0 + 180, 360) - 180, 0, atol=1e-6)
    assert np.array_equal(sc.dec, s.dec)
    assert not np.allclose(sc.mjd, s.mjd)


def test_injected_ghosts_raise_the_bin_count():
    rng = np.random.default_rng(2)
    n = 300
    s = _sample(
        np.sort(rng.uniform(55000, 60000, n)), rng.uniform(0, 360, n), rng.uniform(-60, 60, n)
    )
    w = np.full(n, 0.5)
    base = nf.wide_pair_stats(s, w, GP)[3, 0]
    si, wi = nf.inject_ghosts(s, w, 1.0, 3600.0, 86400.0, np.random.default_rng(3))
    n_ghost = len(si.mjd) - n
    assert 100 < n_ghost < 200  # about r_g * w * n = 150 parents (edge losses small)
    assert np.all(wi[n:] == 0.5)
    after = nf.wide_pair_stats(si, wi, GP)[3, 0]
    assert after - base >= 0.9 * n_ghost  # each ghost is wide of its parent at 1 h-1 d


def test_pair_ratio_and_empirical_p():
    import nf_ghost as G

    w = np.array([0.5, 0.5, 0.5, 0.5])
    r = G.pair_ratio(w, np.concatenate([w, [0.5, 0.5]]))
    assert r[0] == pytest.approx(30 / 12) and r[1] == pytest.approx(30 / 12)
    null = np.arange(100, dtype=float)[:, None, None] * np.ones((1, 1, 2))
    p = G.empirical_p(np.array([[98.0, 1000.0]]), null)
    assert p[0, 0] == pytest.approx(3 / 101) and p[0, 1] == pytest.approx(1 / 101)
