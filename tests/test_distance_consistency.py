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


def test_draw_from_hist_matches_binned_pdf():
    rng = np.random.default_rng(1)
    edges = np.linspace(-0.2, 1.0, 13)
    pdf = np.zeros(12)
    pdf[2], pdf[3] = 1.0, 3.0
    x = dc.draw_from_hist(pdf, edges, 40_000, rng)
    assert x.min() >= edges[2] and x.max() <= edges[4]
    assert np.mean(x >= edges[3]) == pytest.approx(0.75, abs=0.01)


def test_resample_honours_weights():
    rng = np.random.default_rng(2)
    x = dc.resample(np.array([1.0, 2.0]), 20_000, rng, weights=[1.0, 3.0])
    assert np.mean(x == 2.0) == pytest.approx(0.75, abs=0.01)


def test_ddt_with_kext_convention():
    # hierArc/TDCOSMO: D_dt^model = (1 - kappa_ext) D_dt, so a positive kappa_ext raises D_dt
    assert dc.ddt_with_kext(np.array([900.0]), np.array([0.1]))[0] == pytest.approx(1000.0)


def test_load_array_pickle_reads_arrays_and_refuses_other_classes(tmp_path):
    import pickle

    good = tmp_path / "good.pkl"
    good.write_bytes(pickle.dumps([np.arange(3.0), np.ones(2)], protocol=2))
    a, b = dc.load_array_pickle(good)
    assert a.tolist() == [0.0, 1.0, 2.0] and b.tolist() == [1.0, 1.0]

    class Evil:
        def __reduce__(self):
            return (print, ("ran",))

    bad = tmp_path / "bad.pkl"
    bad.write_bytes(pickle.dumps([Evil()], protocol=2))
    with pytest.raises(pickle.UnpicklingError):
        dc.load_array_pickle(bad)


def test_load_data_pickle_reads_containers_and_refuses_other_classes(tmp_path):
    import collections
    import pickle

    good = tmp_path / "lens.pkl"
    good.write_bytes(pickle.dumps([{"z_lens": 0.3, "j_model": np.arange(3.0)}], protocol=2))
    (lens,) = dc.load_data_pickle(good)
    assert lens["z_lens"] == 0.3 and np.array_equal(lens["j_model"], np.arange(3.0))
    bad = tmp_path / "bad.pkl"
    bad.write_bytes(pickle.dumps(collections.OrderedDict(a=1)))
    with pytest.raises(pickle.UnpicklingError):
        dc.load_data_pickle(bad)


def _synthetic_kin_lens(ratio, n_bin=1, gamma_axis=False):
    j = np.full(n_bin, 4e-7)
    sig = np.sqrt(j * ratio) * dc.C_KMS  # exact sigma_v at a_ani-scaling 1
    lens = {
        "j_model": j,
        "sigma_v_measurement": sig,
        "error_cov_measurement": np.diag((0.02 * sig) ** 2),
        "error_cov_j_sqrt": np.zeros((n_bin, n_bin)) if n_bin > 1 else np.array(0.0),
        "kin_scaling_param_list": ["a_ani"],
        "j_kin_scaling_param_axes": [np.linspace(-0.5, 1.0, 7)],
        "j_kin_scaling_grid_list": [np.ones(7)] * n_bin,
        "prior_list": [],
    }
    if gamma_axis:
        lens["kin_scaling_param_list"] = ["a_ani", "gamma_pl"]
        lens["j_kin_scaling_param_axes"] = [np.linspace(-0.5, 1.0, 7), np.linspace(1.5, 2.5, 5)]
        lens["j_kin_scaling_grid_list"] = [np.ones((7, 5))] * n_bin
        lens["prior_list"] = [["gamma_pl", 2.0, 0.05]]
    return lens


@pytest.mark.parametrize("n_bin,gamma_axis", [(1, False), (3, True)])
def test_kin_ln_ratio_recovers_injected_ratio(n_bin, gamma_axis):
    x = np.linspace(np.log(0.1), np.log(20.0), 1500)
    lnl, info = dc.kin_ln_ratio_loglike(_synthetic_kin_lens(1.8, n_bin, gamma_axis), x)
    _, med, hw = dc.grid_moments(x, lnl)
    assert med == pytest.approx(np.log(1.8), abs=0.01)
    # 2 % sigma_v errors per bin -> 4 % on sigma_v^2, sqrt(n_bin) better with independent bins
    assert hw == pytest.approx(0.04 / np.sqrt(n_bin), rel=0.1)
    assert info["n_bins"] == n_bin and info["gamma_pl"] == ("prior" if gamma_axis else None)


def test_grid_pull_and_offset_scatter_fit():
    x = np.linspace(-10, 10, 4001)
    dens = np.exp(-0.5 * (x - 3.0) ** 2) / np.sqrt(2 * np.pi)
    assert dc.grid_pull(x, dens, 0.0, 1e-9) == pytest.approx(3.0, abs=0.01)
    assert dc.grid_pull(x, dens, 6.0, np.sqrt(8.0)) == pytest.approx(-1.0, abs=0.01)
    rng = np.random.default_rng(1)
    s = np.full(400, 0.1)
    d, sd, tau = dc.offset_scatter_fit(0.2 + rng.normal(0, 0.1, 400), s)
    assert d == pytest.approx(0.2, abs=0.03) and tau < 0.05
    d, sd, tau = dc.offset_scatter_fit(-0.1 + rng.normal(0, np.hypot(0.1, 0.3), 400), s)
    assert d == pytest.approx(-0.1, abs=0.06) and tau == pytest.approx(0.3, abs=0.05)
