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


# --- IceTracks-DR2 helpers (E-NF1b) ---


def test_read_icetracks_tab_and_dedupe():
    lines = [
        "#  run  event  subevent  MJD  log10(E/GeV)  AngErr  RA  Dec  Az  Zen",
        '"  110782  744977  255  54562.379  3.31  0.48  203.0  16.7  218.3  106.7"',
        '"  110782  744977  255  54562.379  3.31  0.48  203.0  16.7  218.3  106.7"',
        '"  110782  1620385  255  54562.300  5.17  0.62  75.0  -13.3  349.1  76.6"',
    ]
    a = nf.read_icetracks_tab(lines)
    assert a.shape == (3, 10)
    d = nf.dedupe_events(a)
    assert len(d) == 2 and np.all(np.diff(d[:, 3]) >= 0)


def test_in_uptime_overlapping_intervals():
    start, stop = np.array([0.0, 5.0, 6.0]), np.array([1.0, 10.0, 7.0])
    got = nf.in_uptime(np.array([-1.0, 0.5, 2.0, 6.5, 9.0, 11.0]), start, stop)
    assert got.tolist() == [False, True, False, True, True, False]


def test_jitter_uptime_lands_in_good_runs_and_keeps_hour_angle():
    rng = np.random.default_rng(1)
    t = np.linspace(55000.2, 55009.8, 200)
    s = en.Sample("x", t, np.full(200, 10.0), np.zeros(200), np.ones(200), en.mjd_year(t))
    start = np.arange(55000.0, 55010.0, 1.0)
    stop = start + 0.5
    j = nf.jitter_uptime(s, start, stop, 3.0, rng)
    moved = j.mjd != t
    assert moved.mean() > 0.9
    assert nf.in_uptime(j.mjd[moved], start, stop).all()
    ha0 = en.gmst_deg(t) - s.ra
    ha1 = en.gmst_deg(j.mjd) - j.ra
    assert np.allclose(np.mod(ha1 - ha0 + 180, 360) - 180, 0, atol=1e-6)


def test_wide_pair_counts_drops_same_readout():
    t = 55000 + np.array([0.0, 5.0, 50.0]) / en.DAY
    s = en.Sample(
        "x", t, np.array([0.0, 90.0, 180.0]), np.zeros(3), np.full(3, 0.5), en.mjd_year(t)
    )
    edges = (0.0, 10.0, 100.0)
    assert nf.wide_pair_counts(s, np.array([1, 2, 3]), edges).tolist() == [1.0, 2.0]
    assert nf.wide_pair_counts(s, np.array([1, 1, 3]), edges).tolist() == [0.0, 2.0]


def test_union_days_counts_overlap_once():
    assert nf.union_days(np.array([0.0, 0.5, 3.0]), np.array([1.0, 2.0, 4.0])) == 3.0
