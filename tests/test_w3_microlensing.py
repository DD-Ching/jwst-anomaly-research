"""Offline tests for scripts/w3_microlensing.py on synthetic light curves (D-057)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from jwst_anomaly import exotic_sim as es

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w3_microlensing", _DIR / "w3_microlensing.py")
w3 = importlib.util.module_from_spec(_spec)
sys.modules["w3_microlensing"] = w3  # dataclasses resolve annotations through sys.modules
_spec.loader.exec_module(w3)

RA, DEC = 270.0, -29.0


def _cadence(n=500, span=1600.0, seed=1):
    rng = np.random.default_rng(seed)
    return np.sort(2456000.0 + rng.uniform(0.0, span, n))


def _synthetic(t, a, fs=0.6, fb=0.4, sigma=0.01, seed=2):
    """Flux light curve fs·A + fb with Gaussian noise of a fixed relative size (simulated)."""
    rng = np.random.default_rng(seed)
    f = fs * a + fb
    sf = sigma * np.sqrt(np.maximum(f, 0.05))
    return w3.LightCurve(t, f + rng.normal(0.0, sf), sf, RA, DEC, "synthetic")


def test_linear_fluxes_recovers_source_and_blend_and_applies_the_published_bounds():
    t = _cadence()
    a = w3.pspl(w3.straight_beta(t, 2456800.0, 20.0, 0.1))
    f = 0.7 * a + 0.3
    fs, fb, chi2 = w3.linear_fluxes(a, f, np.ones_like(f))
    assert fs == pytest.approx(0.7) and fb == pytest.approx(0.3) and chi2 < 1e-15
    fs, fb, _ = w3.linear_fluxes(a, 1.2 * a - 0.5, np.ones_like(f))
    assert fb == pytest.approx(-0.1)  # F_b ≥ −F_min (I = 20.5)
    fs, fb, _ = w3.linear_fluxes(a, 2.0 - 0.1 * a, np.ones_like(f))
    assert fs == 0.0  # no negative source flux: a dip cannot be a negative PSPL
    k = w3.linear_fluxes(np.vstack([a, a]), f, np.ones_like(f))
    assert k[0].shape == (2,)


@pytest.mark.parametrize("sign", [1, -1])
def test_fast_magnification_matches_exotic_sim(sign):
    beta = np.concatenate([np.geomspace(1e-3, 40.0, 400), [1.8901, 2.0, 2.0001, 2.5]])
    for n in (1.0, 2.0):
        ref = es.total_magnification(beta, n, sign)
        got = w3.point_magnification(beta, n, sign)
        ok = np.isfinite(ref) & np.isfinite(got) & (ref > 0)
        assert np.allclose(got[ok], ref[ok], rtol=2e-4)
        fin = np.isfinite(ref)  # exactly on the caustic exotic_sim returns inf
        assert np.all((got[fin] == 0) == (ref[fin] == 0))
        rho = 0.05
        b = np.linspace(0.0, 4.0, 81) if sign == 1 else np.linspace(1.5, 2.6, 81)
        fs_ref = es.finite_source_magnification(b, rho, n, sign)
        fs_got = w3.exotic_magnification(b, n, sign, rho)
        assert np.allclose(fs_got, fs_ref, rtol=5e-3, atol=1e-3)


def test_fitter_recovers_pspl_parameters_and_prefers_ordinary():
    t = _cadence()
    a = w3.pspl(w3.straight_beta(t, 2456800.0, 25.0, 0.15))
    lc = _synthetic(t, a)
    res = w3.fit_event(lc, 2456805.0, 15.0, 0.4, models=("PSPL", "FSPL", "N1neg", "E2pos"))
    ps = res["PSPL"]
    truth = w3.chi2_of("PSPL", lc, {"t0": 2456800.0, "tE": 25.0, "u0": 0.15})[0]
    assert ps["chi2"] <= truth + 1e-3  # the optimiser reaches at least the true minimum
    assert ps["t0"] == pytest.approx(2456800.0, abs=0.3)
    assert ps["tE"] == pytest.approx(25.0, rel=0.1)
    assert ps["u0"] == pytest.approx(0.15, rel=0.1)
    assert ps["fs"] == pytest.approx(0.6, rel=0.1)
    row = w3.summarise("synthetic", res, t.size)
    assert row["dbic_N1neg"] > w3.P.flag_dbic and row["dbic_E2pos"] > w3.P.flag_dbic


def test_injected_negative_mass_event_is_flagged_and_fails_the_published_selection():
    t = _cadence(n=800)
    a = es.light_curve(t, 2456800.0, 30.0, 0.8, 1.0, -1, 0.05)  # umbra between two spikes
    lc = _synthetic(t, a, fs=0.6, fb=0.4)
    res = w3.fit_event(lc, 2456770.0, 10.0, 0.3, models=("PSPL", "FSPL", "N1neg"))
    row = w3.summarise("synthetic", res, t.size)
    assert row["dbic_N1neg"] < -100  # simulated: the exotic model is far better
    assert row["N1neg_u0"] == pytest.approx(0.8, abs=0.1)
    sel = w3.published_selection(lc, res["PSPL"])
    assert not sel["selected"]
    assert sel["n_bump"] == 2 or not sel["passed"]["chi2_fit"]


def test_published_selection_accepts_a_clean_pspl_event():
    t = _cadence(n=900, span=2700.0)
    a = w3.pspl(w3.straight_beta(t, 2457000.0, 20.0, 0.2))
    lc = _synthetic(t, a, fs=0.5, fb=0.1)
    fit = w3.fit_event(lc, 2457000.0, 20.0, 0.2, models=("PSPL",))["PSPL"]
    sel = w3.published_selection(lc, fit)
    assert sel["selected"], sel["passed"]


def test_injection_keeps_cadence_and_vanishes_the_lensed_flux():
    t = _cadence()
    base = w3.LightCurve(t, np.full_like(t, 1.0), np.full_like(t, 0.01), RA, DEC, "b")
    lc = w3.inject_w3(base, 1.0, 0.6, 2456800.0, 30.0, 0.5, 0.01)
    inside = np.abs(t - 2456800.0) < 30.0 * np.sqrt(4 - 0.25) * 0.9
    assert np.allclose(lc.f[inside], 0.4)  # only the unlensed 40 % remains in the umbra
    assert np.all(lc.sf[inside] < 0.01)


def test_rate_mass_conversion_matches_the_stated_geometry():
    # θ_E = sqrt(κ M π_rel), κ = 8.144 mas/M☉, π_rel = 1/4 − 1/8 mas; μ = 5 mas/yr
    te = w3.einstein_time_days(1.0)[0]
    assert te == pytest.approx(np.sqrt(8.144 * 0.125) / 5.0 * 365.25, rel=2e-3)


def test_trajectory_and_magnification_match_mulensmodel():
    mm = pytest.importorskip("MulensModel")
    t = _cadence()
    lc = w3.LightCurve(t, np.ones_like(t), np.ones_like(t), RA, DEC, "x")
    p = {"t0": 2456800.0, "tE": 60.0, "u0": 0.2, "pi_E_N": 0.3, "pi_E_E": -0.2, "t0_par": 2456800}
    x, y = w3.trajectory(lc, p)
    ref = w3.mm_model(lc, p).get_trajectory(t)
    assert np.allclose(x, ref.x, atol=1e-10) and np.allclose(y, ref.y, atol=1e-10)
    m = mm.Model({"t_0": 2456800.0, "u_0": 0.2, "t_E": 60.0})
    q = {"t0": 2456800.0, "tE": 60.0, "u0": 0.2}
    assert np.allclose(w3.magnification("PSPL", lc, q), m.get_magnification(t), rtol=1e-12)
    assert w3.mm_crosscheck("FSPL", lc, {**q, "u0": 0.01, "rho": 0.02}) < 0.01


def test_isolated_outliers_and_season_offsets():
    t = np.concatenate([np.arange(0.0, 100.0), np.arange(365.0, 465.0)])
    f = np.where(t < 200, 1.0, 1.05)  # a season-to-season zero-point jump
    f[50] = 1.2  # one isolated outlier
    f[150:153] = 1.2  # not isolated: three neighbours
    lc = w3.LightCurve(t, f, np.full_like(t, 0.01))
    bad = w3.isolated_outliers(lc, np.where(t < 200, 1.0, 1.05))
    assert list(np.nonzero(bad)[0]) == [50]
    lcs = lc.with_season_offsets()
    assert lcs.n_extra == 1 and lcs.seasons.shape == (2, t.size)
    assert lc.with_season_offsets(trend=True).n_extra == 3
    p = {"t0": 50.0, "tE": 5.0, "u0": 1.0}
    assert w3.chi2_of("PSPL", lcs, p)[0] < w3.chi2_of("PSPL", lc, p)[0]
    sub = lcs.subset(t < 200)
    assert sub.n_extra == 0 and sub.seasons.shape == (1, 100)


def test_spike_pair_starts_bracket_the_umbra():
    t = _cadence(n=900)
    a = es.light_curve(t, 2456800.0, 40.0, 1.0, 1.0, -1, 0.02)
    lc = _synthetic(t, a, sigma=0.002)
    starts = w3.spike_pair_starts(lc, 2.0)
    assert any(abs(s[0] - 2456800.0) < 5 and abs(10 ** s[1] - 40.0) < 15 for s in starts)


def test_fit_checkpoint_drops_a_torn_last_line(tmp_path):
    p = tmp_path / "fits.partial.jsonl"
    assert w3.load_checkpoint(p) == []
    p.write_text('{"event_id": "a", "x": NaN}\n{"event_id": "b"}\n{"event_id": "c", "x"')
    rows = w3.load_checkpoint(p)
    assert [r["event_id"] for r in rows] == ["a", "b"] and np.isnan(rows[0]["x"])


@pytest.mark.parametrize(
    "flag",
    [
        {"event_id": "A", "tests": [], "survives": True},
        {"event_id": "B", "tests": [], "survives": False, "complete": False},
    ],
)
def test_limit_refuses_a_zero_event_limit_unless_vetting_is_a_complete_null(
    tmp_path, monkeypatch, flag
):
    import json

    (tmp_path / "vetting_bulge2019.json").write_text(json.dumps({"flags": [flag]}))
    monkeypatch.setattr(w3, "out_dir", lambda: tmp_path)
    with pytest.raises(SystemExit, match=flag["event_id"]):
        w3.run_limit()
