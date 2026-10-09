"""Offline tests for jwst_anomaly.distance_consistency (D1)."""

import numpy as np
import pytest
from astropy.cosmology import FlatLambdaCDM, FlatwCDM

from jwst_anomaly import distance_consistency as dc

SMALL = dc.Params(n_om=21, n_lns=121, n_pred=20_000, n_obs=20_000)


def test_ratio_model_matches_astropy_distances():
    zd, zs = 0.5, 2.0
    for cosmo, w in (
        (FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0), -1.0),
        (FlatwCDM(H0=70, Om0=0.3, w0=-1.3, Tcmb0=0), -1.3),
    ):
        r_ast = (
            cosmo.angular_diameter_distance(zs) / cosmo.angular_diameter_distance(zd, zs)
        ).value
        assert dc.ratio_model(zd, zs, 0.3, w)[0] == pytest.approx(r_ast, rel=1e-5)
        dd = cosmo.angular_diameter_distance(zd).value
        ddt = (1 + zd) * dd * r_ast
        assert dc.ratio_from_samples(np.array([dd]), np.array([ddt]), zd)[0] == pytest.approx(r_ast)
        hubble = 299792.458 / 70
        assert dc.ddt_model_dimensionless(zd, zs, 0.3, w)[0] * hubble == pytest.approx(
            ddt, rel=1e-5
        )


def test_ratio_is_h0_free_and_weakly_om_dependent():
    r1 = dc.ratio_model(0.6, 1.5, 0.2)[0]
    r2 = dc.ratio_model(0.6, 1.5, 0.4)[0]
    assert abs(np.log(r1 / r2)) < 0.1  # ~7 % over Omega_m 0.2-0.4 (vs ~25 % for D_dt)


def test_pull_recovers_known_offset():
    rng = np.random.default_rng(0)
    obs = rng.normal(0.5, 0.1, 50_000)
    pred = np.zeros(50_000)
    r = dc.pull(obs, pred, rng)
    assert r["z"] == pytest.approx(5.0, abs=0.1)
    assert r["delta_ln"] == pytest.approx(0.5, abs=0.01)
    assert r["z_emp"] > 3.5  # saturates at 1/N


def test_leave_one_out_consistent_and_injected():
    rng = np.random.default_rng(1)
    zd = np.array([0.3, 0.45, 0.6, 0.7, 0.35])
    zs = np.array([0.7, 1.5, 1.4, 1.8, 1.7])
    om, grid = dc.grid_ln_model(zd, zs, "ratio", SMALL)
    truth = np.log(dc.ratio_model(zd, zs, 0.3))
    sig = 0.1
    lns = np.linspace(-0.7, 0.7, SMALL.n_lns)

    def obs(shift):
        return [
            rng.normal(t + (shift if i == 2 else 0.0), sig, 20_000) for i, t in enumerate(truth)
        ]

    res = dc.leave_one_out(obs(0.0), grid, om, lns, rng, SMALL)
    assert max(abs(r["z"]) for r in res) < 2.0
    res = dc.leave_one_out(obs(0.8), grid, om, lns, rng, SMALL)
    assert res[2]["z"] > 4.0  # the shifted object stands out
    assert abs(res[0]["z"]) < res[2]["z"]


def test_prior_predictive_ratio_covers_model():
    rng = np.random.default_rng(2)
    d = dc.prior_predictive_ratio(0.5, 2.0, rng, SMALL)
    assert np.exp(d.min()) <= dc.ratio_model(0.5, 2.0, 0.3)[0] <= np.exp(d.max())


def test_trials_threshold_round_trip():
    thr = dc.local_sigma_threshold(18, 5.0)
    assert thr > 5.0
    assert dc.global_sigma(thr, 18) == pytest.approx(5.0, abs=1e-6)
    assert dc.local_sigma_threshold(1, 5.0) == pytest.approx(5.0)


def test_one_sided_threshold():
    thr1 = dc.local_sigma_threshold(188, 5.0, one_sided=True)
    assert thr1 == pytest.approx(5.81, abs=0.01)
    assert thr1 < dc.local_sigma_threshold(188, 5.0)


def test_macquart_mean_and_delta():
    dm1 = dc.macquart_mean_dm(1.0)[0]
    assert 850 < dm1 < 1050  # ~ 1000 z pc cm^-3 at z ~ 1 (Macquart et al. 2020)
    rng = np.random.default_rng(3)
    d = dc.draw_macquart_delta(0.5, rng, n=50_000)
    assert d.mean() == pytest.approx(1.0, abs=0.03)
    assert d.min() >= 0


def test_frb_grid_matches_monte_carlo():
    p = dc.FRB_DEFAULT
    dm_ism, z = 50.0, 0.4
    pred = dc.FRBPredictive(dm_ism, z, p)
    rng = np.random.default_rng(5)
    n = 200_000
    mc = (
        dm_ism * (1 + p.ism_frac_err * rng.standard_normal(n))
        + rng.uniform(*p.halo_range, n)
        + pred.mean_dm * dc.draw_macquart_delta(z, rng, p, n)
        + np.exp(p.host_mu + p.host_sigma * rng.standard_normal(n)) / (1 + z)
    )
    for d in np.percentile(mc, [1, 10, 50, 90, 99]):
        lo, hi = pred.tails(d)
        assert lo == pytest.approx(np.mean(mc <= d), abs=0.003)
        assert lo + hi == pytest.approx(1.0, abs=1e-6)


def test_frb_tails_reach_beyond_threshold():
    # A Monte Carlo floored at 1/N caps near 4.3-4.8 sigma; the grid must reach the 5.81 sigma flag.
    p = dc.FRB_DEFAULT
    thr = dc.local_sigma_threshold(188, 5.0, one_sided=True)
    low = dc.frb_predictive_tails(20.0, 40.0, 0.8, p)  # far below the Milky-Way floor
    high = dc.frb_predictive_tails(60_000.0, 40.0, 0.8, p)  # beyond any host + cosmic draw
    ok = dc.frb_predictive_tails(900.0, 40.0, 0.8, p)
    assert low["z_low"] > thr and low["p_below_floor"] > 0.5
    assert high["z_high"] > thr
    assert ok["z_low"] < 2.0 and ok["z_high"] < 2.0
    pred = dc.FRBPredictive(40.0, 0.8, p)
    lo, hi = pred.detect_limits(thr)
    assert 0 < lo < ok["dm_pred_median"] < hi
    assert dc.frb_predictive_tails(0.9 * lo, 40.0, 0.8, p)["z_low"] >= thr
    assert dc.frb_predictive_tails(1.1 * hi, 40.0, 0.8, p)["z_high"] >= thr
    assert dc.frb_predictive_tails(1.1 * lo, 40.0, 0.8, p)["z_low"] < thr
