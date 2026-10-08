"""Offline tests for the W5 count-map adapter and screen helpers (D-063), on synthetic maps."""

from __future__ import annotations

import math

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import countmap as cm
from jwst_anomaly import exotic_sim as es
from jwst_anomaly import schema, signatures


def test_region_chunks_tile_the_box():
    r = cm.Region("t", 20.0, 25.0, -30.0, -26.0, step=2.0)
    ch = r.chunks()
    assert len(ch) == 3 * 2
    assert sum((b[1] - b[0]) * (b[3] - b[2]) for b in ch) == pytest.approx(5 * 4)
    # spherical area: 5° × (sin(-26°) - sin(-30°)) in deg²
    exp = (180 / math.pi) ** 2 * math.radians(5) * (math.sin(math.radians(-26)) + 0.5)
    assert r.area_deg2() == pytest.approx(exp)


def test_chunk_queries_use_plain_group_by():
    q = cm.chunk_queries((20.0, 22.0, -30.0, -28.0), 23.5)
    assert set(q) == {"gal", "all", "bad"}
    for s in q.values():
        assert s.endswith("GROUP BY nest4096")
        assert "ra>=20.0 AND ra<22.0" in s
        assert "CASE" not in s
    assert "dered_mag_r<23.50" in q["gal"]
    assert "maskbits<>0" in q["bad"]


def test_merge_chunk_and_weight():
    all_ = Table(
        {
            "nest4096": [5, 3, 9],
            "n_all": [10, 20, 4],
            "galdepth_r": [500.0, 500.0, 500.0],
            "nobs_g": [3, 4, 5],
            "nobs_r": [2, 4, 5],
            "nobs_z": [3, 1, 5],
            "ebv": [0.02, 0.03, 0.04],
        }
    )
    gal = Table({"nest4096": [3, 9], "n_gal": [7, 1]})
    bad = Table({"nest4096": [5], "n_bad": [5]})
    t = cm.finalize_map(cm.merge_chunk(gal, all_, bad), source="test")
    assert list(t["pix"]) == [3, 5, 9]
    assert list(t["n_gal"]) == [7, 0, 1]
    assert list(t["n_bad"]) == [0, 5, 0]
    assert list(t["nobs_min"]) == [1, 2, 5]
    assert np.allclose(t["w"], [1.0, 0.5, 1.0])
    # 5 sigma depth for ivar 500: 22.5 - 2.5 log10(5 / sqrt(500))
    assert t["depth_r"][0] == pytest.approx(22.5 - 2.5 * math.log10(5 / math.sqrt(500)))
    assert t.meta["provenance"] == schema.Provenance.OBSERVED.value
    assert t.meta["source"] == "test"


def test_combine_duplicates_sums_split_pixels():
    """A pixel straddling a chunk edge comes back once per chunk; counts must add."""
    t = Table(
        {
            "pix": [7, 3, 7],
            "n_gal": [2, 5, 3],
            "n_all": [4, 9, 6],
            "n_bad": [0, 1, 2],
            "depth_r": [24.5, 24.0, 24.5],
            "nobs_min": [3.0, 2.0, 1.0],
            "ebv": [0.02, 0.05, 0.03],
        }
    )
    c = cm.combine_duplicates(t)
    assert list(c["pix"]) == [3, 7]
    assert list(c["n_gal"]) == [5, 5] and list(c["n_all"]) == [9, 10]
    assert list(c["n_bad"]) == [1, 2]
    assert c["depth_r"][1] == pytest.approx(24.5)
    assert c["nobs_min"][1] == 1.0
    assert c["ebv"][1] == pytest.approx((4 * 0.02 + 6 * 0.03) / 10)


def test_count_map_survey_protocol():
    s = cm.LegacySurveysCountMap(cm.Region("t", 0, 1, 0, 1), cache=None)
    assert isinstance(s, signatures.CountMapSurvey)
    assert not isinstance(s, signatures.CatalogueSurvey)  # pixels are not objects


def test_tabulated_counts_match_power_law():
    m = np.arange(16.0, 23.51, 0.5)
    n = 10 ** (0.4 * (m - 23.5))  # dlogN/dm = 0.4 -> alpha = 1
    f = cm.tabulated_counts(m, n, bright_slope=0.4)
    s = np.array([0.5, 1.0, 10.0, 1e4])
    assert np.allclose(f(s), s ** (-1.0), rtol=1e-9)
    with pytest.raises(ValueError):
        cm.tabulated_counts(m[::-1], n)


def _realistic_counts():
    """Broken power law like r-band galaxy counts: dlogN/dm = 0.6 bright of r = 19, 0.4 fainter."""
    m = np.arange(16.0, 23.51, 0.5)
    lg = np.where(m < 19, 0.6 * (m - 19), 0.4 * (m - 19))
    return cm.tabulated_counts(m, 10**lg, bright_slope=0.6)


def _profile(x):
    return cm.deficit_profile(x, _realistic_counts())


def test_profile_limits():
    """alpha = 1: no change; steep counts: deficit in the core and excess near x = 1."""
    x = np.array([0.1, 0.5, 0.9, 1.1, 3.0])
    assert np.allclose(cm.deficit_profile(x, cm.power_counts(1.0)), 1.0)
    p = cm.deficit_profile(x, cm.power_counts(1.5))
    mu = np.abs(es.image_plane_magnification(x, 1.0, -1))
    assert np.allclose(p, mu**0.5)
    assert p[0] < 0.1 and p[2] > 1.0
    # capped at the critical curve; matches exotic_sim.count_ratio away from it
    assert cm.deficit_profile(np.array([1.0]), cm.power_counts(1.5))[0] == pytest.approx(30**0.5)
    c = _realistic_counts()
    assert np.allclose(_profile(x[:3]), es.count_ratio(x[:3], c, 1.0, n=1.0, sign=-1))
    assert _profile(np.array([0.1]))[0] < 0.1  # realistic counts: an empty core


def test_expected_missing_scales_with_area():
    pf = _profile
    m1, n1 = cm.expected_missing(2.0, 10.0, pf, 1.0)
    m2, n2 = cm.expected_missing(4.0, 10.0, pf, 1.0)
    assert n1 == pytest.approx(math.pi * 4 * 10, rel=1e-3)
    assert m2 == pytest.approx(4 * m1, rel=1e-6)
    assert 0 < m1 < n1


def test_block_sum_and_rasterise():
    img = np.arange(16.0).reshape(4, 4)
    assert cm.block_sum(img, 2).tolist() == [[10.0, 18.0], [42.0, 50.0]]
    idx = np.array([[0, 0, 1], [-1, 1, 1]])
    d = cm.rasterise(np.array([4.0, 9.0]), idx, per_cell=True)
    assert d.sum() == pytest.approx(13.0)
    assert d[0, 0] == pytest.approx(2.0) and d[1, 0] == 0.0


def _synthetic(theta_cells=8.0, n_cell=2.0, size=240, inject=True, seed=3):
    rng = np.random.default_rng(seed)
    pf = _profile
    yy, xx = np.mgrid[:size, :size] + 0.5
    c = size / 2
    lam = np.full((size, size), n_cell)
    if inject:
        lam = lam * pf(np.hypot(xx - c, yy - c) / theta_cells)
    weight = np.ones((size, size))
    weight[:, :10] = 0.0  # a masked strip
    d = rng.poisson(lam) * (weight > 0)
    return d.astype(float), weight, pf, int(c)


def test_matched_filter_recovers_injected_amplitude():
    d, w, pf, c = _synthetic()
    mf = cm.MatchedFilter(w, 8.0, pf, 2.5)
    amp, sig = mf(d)
    assert amp[c, c] == pytest.approx(1.0, abs=0.15)
    assert amp[c, c] / sig[c, c] > 3  # Poisson S/N ≈ 4.4 for 2 galaxies per cell, θ_E = 8 cells
    assert mf.full[c, c] == pytest.approx(1.0) and mf.cover[c, c] == pytest.approx(1.0)
    iy, ix = cm.find_peaks(amp / sig, (mf.full > 0.8) & (mf.cover > 0.8), 8, 3.0)
    assert len(iy) >= 1
    assert np.min(np.hypot(iy - c, ix - c)) <= 3


def test_matched_filter_null_and_masked_edge():
    d, w, pf, c = _synthetic(inject=False)
    mf = cm.MatchedFilter(w, 8.0, pf, 2.5)
    amp, sig = mf(d)
    ok = (mf.full > 0.8) & (mf.cover > 0.8)
    z = (amp / sig)[ok]
    assert abs(np.median(z)) < 0.3
    assert np.std(z) < 2.0
    # the masked strip is not usable (coverage), and is not a deficit
    assert mf.cover[c, 2] < 0.8


def test_inject_deficit_thins_core_and_conserves_far():
    rng = np.random.default_rng(1)
    pf = _profile
    g = np.arange(-40, 41, 0.86)
    x, y = np.meshgrid(g, g)
    x, y = x.ravel(), y.ravel()
    n = np.full(x.size, 50)
    out = cm.inject_deficit(n, (x, y), (0.0, 0.0), 8.0, pf, rng, 2.5)
    r = np.hypot(x, y)
    assert out[r < 2].mean() < 0.35 * 50  # x < 0.25: ratio ≲ 0.3 for these counts
    assert np.all(out[r > 25] == 50)
    assert out.dtype == np.int64 and np.all(out >= 0)
