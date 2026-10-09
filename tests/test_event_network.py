"""Offline tests for the E1 event-network statistic (synthetic catalogues only)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from jwst_anomaly import event_network as en


def _sample(cat, mjd, ra, dec, sigma, year=2020):
    n = len(mjd)
    return en.Sample(
        cat,
        np.asarray(mjd, float),
        np.asarray(ra, float),
        np.asarray(dec, float),
        np.full(n, sigma, float) if np.isscalar(sigma) else np.asarray(sigma, float),
        np.full(n, year),
    )


def _random(cat, n, rng, sigma=1.0, span=300.0):
    return _sample(
        cat,
        59000 + rng.uniform(0, span, n),
        rng.uniform(0, 360, n),
        np.degrees(np.arcsin(rng.uniform(-1, 1, n))),
        sigma,
    )


def test_sep_deg_known_values():
    assert en.sep_deg(0, 0, 90, 0) == pytest.approx(90)
    assert en.sep_deg(10, 89.999, 190, 89.999) == pytest.approx(0.002, abs=1e-6)
    assert en.sep_deg(0, 0, 180, 0) == pytest.approx(180)


@pytest.mark.parametrize("same", [True, False])
def test_pairs_within_matches_brute_force(same):
    rng = np.random.default_rng(1)
    ta = 59000 + rng.uniform(0, 30, 200)
    tb = ta if same else 59000 + rng.uniform(0, 30, 150)
    i, j = en.pairs_within(ta, tb, 2 * en.DAY, same)
    got = {(int(a), int(b)) for a, b in zip(i, j, strict=True)}
    want = set()
    for a in range(len(ta)):
        for b in range(len(tb)):
            if abs(ta[a] - tb[b]) <= 2 and (not same or a != b):
                want.add((min(a, b), max(a, b)) if same else (a, b))
    if same:
        got = {(min(a, b), max(a, b)) for a, b in got}
        assert len(i) == len(got)  # each unordered pair once
    assert got == want


def test_sep_class_same_wide_unlocalized():
    p = en.Params()
    cls = en.sep_class(
        np.array([1.0, 10.0, 0.05, 5.0]),
        np.array([1.0, 1.0, 0.001, np.nan]),
        np.array([1.0, 1.0, 0.001, 1.0]),
        p,
    )
    # third pair: beyond 3 sigma but below the 0.1 deg floor -> -1 ("between")
    assert cls.tolist() == [0, 1, -1, 2]


def test_scramble_perm_keeps_times_dec_and_hour_angle():
    rng = np.random.default_rng(2)
    s = _random("CHIME", 300, rng)
    s.year[150:] = 2021
    sc = en.scramble_perm(s, rng)
    for y in (2020, 2021):
        m = s.year == y
        assert np.allclose(np.sort(sc.mjd[m]), np.sort(s.mjd[m]))
    assert np.array_equal(sc.dec, s.dec)
    ha0 = np.mod(en.gmst_deg(s.mjd) - s.ra, 360)
    ha1 = np.mod(en.gmst_deg(sc.mjd) - sc.ra, 360)
    assert np.allclose(np.mod(ha1 - ha0 + 180, 360) - 180, 0, atol=1e-6)


def test_scramble_jit_gbm_keeps_orbit_phase():
    rng = np.random.default_rng(3)
    p = en.Params()
    s = _random("GBM", 400, rng, sigma=5.0)
    perm_rng = np.random.default_rng(9)
    perm = en.scramble_perm(s, np.random.default_rng(9))
    sc = en.scramble_jit(s, perm_rng, p)
    shift = (sc.mjd - perm.mjd) * en.DAY
    resid = np.mod(shift + en.FERMI_ORBIT / 2, en.FERMI_ORBIT) - en.FERMI_ORBIT / 2
    assert np.all(np.abs(resid) <= p.orbit_slop_s + 1e-6)
    assert np.all(np.abs(shift) <= p.jitter_s + p.orbit_slop_s + 1e-6)


def test_count_channel_bins_and_classes():
    p = en.Params()
    a = _sample("ICECAT", [59000.0, 59010.0], [10, 100], [0, 0], 0.5)
    b = _sample("CHIME", [59000.0 + 5 / en.DAY, 59010.0 + 500 / en.DAY], [10.5, 250], [0, 40], 0.2)
    out = en.count_channel(a, b, p, same=False)
    assert out.shape == (len(en.windows_for(p)), 3)
    assert out[0].tolist() == [1, 0, 1]  # 0-10 s: same direction
    assert out[2].tolist() == [0, 1, 1]  # 100 s-1 h: wide
    assert out[:5, 2].sum() == 2


def test_union_window_counts_orbit_multiples():
    p = en.Params()
    a = _sample("GBM", [59000.0, 59000.0 + 3 * en.FERMI_ORBIT / en.DAY], [0, 90], [0, 0], 5.0)
    out = en.count_channel(a, a, p, same=True)
    names = [n for n, _ in en.windows_for(p)]
    assert out[names.index("orbit_k"), 2] == 1
    assert out[names.index("solar_k"), 2] == 0


def test_global_p_trials_correction():
    rng = np.random.default_rng(4)
    null = rng.poisson(5.0, size=(2000, 4, 3, 2))
    mask = np.ones((4, 3, 2), bool)
    typical = null[0]
    _, gp = en.global_p(typical, null, mask)
    assert gp > 0.05
    extreme = typical.copy()
    extreme[1, 1, 0] = 40
    pmin, gp2 = en.global_p(extreme, null, mask)
    assert pmin == pytest.approx(1 / 2001)
    assert gp2 < 0.01


def test_injection_recovered_in_wide_channel():
    rng = np.random.default_rng(5)
    p = en.Params()
    a = _random("GBM", 300, rng, sigma=4.0)
    b = _random("ICECAT", 300, rng, sigma=1.0)
    base = en.count_channel(a, b, p, same=False)
    b2 = en.inject_pairs(a, b, 20, 100.0, 3600.0, p, False, rng)
    assert np.count_nonzero(b2.mjd != b.mjd) == 20
    out = en.count_channel(a, b2, p, same=False)
    # moved events also lose their old pairs, so the net gain is a little below 20
    assert out[2, 1] - base[2, 1] >= 14


def test_loaders():
    ice = pd.DataFrame(
        {
            "NAME": ["IC1", "IC2"],
            "EVENTMJD": [59000.0, 59001.0],
            "RA": [1.0, 2.0],
            "DEC": [0.0, 1.0],
            "RA_ERR_PLUS": [2.146, 1.0],
            "RA_ERR_MINUS": [2.146, 1.0],
            "DEC_ERR_PLUS": [2.146, 1.0],
            "DEC_ERR_MINUS": [2.146, 1.0],
            "CR_VETO": ["FALSE", "TRUE"],
        }
    )
    t = en.icecat_events(ice)
    assert len(t) == 1 and t["sigma"][0] == pytest.approx(1.0)
    assert t.meta["provenance"] == "derived"
    chime = pd.DataFrame(
        {
            "tns_name": ["F1", "F2", "F3", "F3"],
            "repeater_name": [np.nan, "R1", "R1", "R1"],
            "sub_num": [0, 0, 0, 1],
            "mjd_400": [59000.0, 59001.0, 59002.0, 59002.0],
            "ra": [1.0, 2.0, 2.0, 2.0],
            "dec": [0.0, 1.0, 1.0, 1.0],
            "ra_err": [0.1, 0.2, 0.2, 0.2],
            "dec_err": [0.3, 0.1, 0.1, 0.1],
        }
    )
    assert len(en.chime_events(chime)) == 2
    assert len(en.chime_events(chime, one_per_source=False)) == 3
    gw = en.gw_events(pd.DataFrame({"name": ["GW170817"], "gps": [1187008882.4]}))
    assert gw["mjd"][0] == pytest.approx(57982.528524, abs=2e-5)  # 12:41:04.4 UTC
    assert not np.isfinite(gw["sigma"][0])
