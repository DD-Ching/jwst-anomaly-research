import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import countmap


def _uniform_field(rng, n_cells=80, mean=4.0):
    gal = rng.poisson(mean, (n_cells, n_cells)).astype(float)
    every = rng.poisson(2.0, gal.shape).astype(float)
    return gal, every.copy(), every


def test_tangent_plane_origin_and_scale():
    x, y = countmap.tangent_plane(np.array([10.0, 10.0]), np.array([0.0, 1.0 / 60]), 10.0, 0.0)
    assert x == pytest.approx([0, 0], abs=1e-9)
    assert y == pytest.approx([0, 1.0], rel=1e-6)


def test_count_grids_split_galaxies_stars_and_coverage():
    tab = Table(
        {
            "ra": [0.001, 0.001, 0.001, 0.001],
            "dec": [0.001, 0.001, 0.001, 0.001],
            "type": ["REX", "PSF", "EXP", "DEV"],
            "mag_r": [20.0, 20.0, 20.0, 23.0],
            "maskbits": [0, 0, 2, 0],
        }
    )
    gal, star, clean, every = countmap.count_grids(tab, 0.0, 0.0, 1.0, 1.0, mag_gal=21.0)
    assert gal.sum() == 1 and star.sum() == 1
    assert clean.sum() == 3 and every.sum() == 4


def test_power_law_counts_extrapolates_with_bright_slope():
    n = countmap.power_law_counts([18, 19, 20], [100, 300, 900], bright_slope=0.6)
    assert n(10 ** (-0.4 * 19)) == pytest.approx(300)
    assert n(10 ** (-0.4 * 17)) == pytest.approx(100 / 10**0.6)


def test_predicted_ratio_is_a_deficit_inside_einstein_radius():
    n = countmap.power_law_counts([16, 18, 20, 21], [10, 120, 1200, 3000])
    r = countmap.predicted_ratio(np.array([0.2, 0.5, 3.0]), n, 21.0)
    assert r[0] < r[1] < 1 and r[2] == pytest.approx(1, abs=0.05)


def test_disk_deficit_null_and_injected_hole():
    rng = np.random.default_rng(1)
    gal, clean, every = _uniform_field(rng)
    _, _, logp = countmap.disk_deficit(gal, clean, every, 6, min_coverage=0.8)
    assert np.nanmin(logp) > -6  # no extreme tail in a Poisson field of 6400 cells

    def ratio(x):
        return np.where(x < 1, 0.05, 1.0)

    holed = countmap.inject_hole(gal, 4.0, (40, 40), 12, ratio, rng)
    obs, exp, logp = countmap.disk_deficit(holed, clean, every, 6, min_coverage=0.8)
    assert np.nanargmin(logp) == pytest.approx(40 * 80 + 40, abs=3 * 80 + 3)
    assert obs[40, 40] < 0.2 * exp[40, 40] and logp[40, 40] < -20


def test_disk_deficit_masks_low_coverage():
    rng = np.random.default_rng(2)
    gal, clean, every = _uniform_field(rng, 40)
    clean[:, :20] = 0  # masked: rows exist, all flagged
    gal[:, :20] = 0
    obs, _, logp = countmap.disk_deficit(gal, clean, every, 4, min_coverage=0.8)
    assert np.isnan(logp[20, 5]) and np.isfinite(logp[20, 35])


def test_clustering_k_recovers_negative_binomial_shape():
    rng = np.random.default_rng(3)
    k_true, mean = 8.0, 20.0
    obs = rng.negative_binomial(k_true, k_true / (k_true + mean), (200, 200)).astype(float)
    k, n_ind = countmap.clustering_k(obs, np.full(obs.shape, mean), 0.5)
    assert n_ind == obs.size and k == pytest.approx(k_true, rel=0.15)
    k_poisson, _ = countmap.clustering_k(
        rng.poisson(mean, (200, 200)).astype(float), np.full((200, 200), mean), 0.5
    )
    assert k_poisson > 100


def test_nb_tail_is_wider_than_poisson():
    lt_nb = countmap.nb_log10_tail([1.0], [17.5], 9.5)
    lt_p = countmap.nb_log10_tail([1.0], [17.5], np.inf)
    assert lt_p[0] < lt_nb[0] < -3  # clustering makes an empty disk much likelier
    assert np.isnan(countmap.nb_log10_tail([np.nan], [1.0], 2.0)[0])
