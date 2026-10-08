"""Offline tests for scripts/w3_moa.py on synthetic difference light curves (D-060)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import exotic_sim as es

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w3_moa", _DIR / "w3_moa.py")
wm = importlib.util.module_from_spec(_spec)
sys.modules["w3_moa"] = wm
_spec.loader.exec_module(wm)


def _cadence(n=900, seed=1):
    """MOA-like cadence: three ~240-d seasons, a few epochs per night (simulated)."""
    rng = np.random.default_rng(seed)
    nights = np.concatenate([np.arange(s, s + 240) for s in (2454000.0, 2454365.0, 2454730.0)])
    t = rng.choice(nights, n) + rng.uniform(0.55, 0.95, n)
    return np.sort(t)


def _noise(t, sigma=100.0, seed=2):
    return np.random.default_rng(seed).normal(0.0, sigma, t.size), np.full(t.size, sigma)


def _w3(t, fs=2000.0, t0=2454480.0, te=10.0, u0=0.8, rho=0.01):
    inj = es.inject_light_curve(t, fs, t0, te, u0, 1.0, -1, rho, 1.0)
    return np.asarray(inj["flux"]) - fs


def test_deficit_scan_finds_an_umbra_and_ignores_white_noise():
    t = _cadence()
    f, sf = _noise(t)
    quiet = wm.deficit_scan(t, f, sf)
    assert quiet["z_min"] > -6.0
    assert not wm.prescreen_pass(quiet["z_min"], quiet["s_min"], quiet["z_min2"])
    sig = _w3(t)
    scan = wm.deficit_scan(t, f + sig, sf)
    assert scan["z_min"] < -wm.P.prescreen_z
    assert 2454480.0 - 20 < 0.5 * (scan["t_lo"] + scan["t_hi"]) < 2454480.0 + 20
    assert wm.prescreen_pass(scan["z_min"], scan["s_min"], scan["z_min2"])


def test_a_pspl_bump_does_not_pass_the_prescreen():
    t = _cadence()
    f, sf = _noise(t)
    inj = es.inject_light_curve(t, 2000.0, 2454480.0, 10.0, 0.2, 1.0, 1, 0.0, 1.0)
    scan = wm.deficit_scan(t, f + np.asarray(inj["flux"]) - 2000.0, sf)
    assert not wm.prescreen_pass(scan["z_min"], scan["s_min"], scan["z_min2"])


def test_repeated_dips_are_rejected_as_a_variable():
    t = _cadence()
    f, sf = _noise(t)
    dips = np.zeros_like(t)
    for c in (2454100.0, 2454480.0, 2454850.0):
        dips[np.abs(t - c) < 6] = -1500.0
    scan = wm.deficit_scan(t, f + dips, sf)
    assert scan["z_min"] < -wm.P.prescreen_z
    assert scan["z_min2"] < -wm.P.prescreen_repeat
    assert not wm.prescreen_pass(scan["z_min"], scan["s_min"], scan["z_min2"])


def test_robust_scale_ignores_a_slow_trend():
    t = _cadence()
    f, sf = _noise(t)
    _b, s = wm.robust_scale(f + 5.0 * (t - t[0]), sf)
    assert s == pytest.approx(1.0, abs=0.15)
    _b, s2 = wm.robust_scale(3.0 * f, sf)
    assert s2 == pytest.approx(3.0, rel=0.15)


def test_cut0_emulation_needs_three_chained_detections():
    t = np.array([0.0, 1.0, 2.0, 20.0, 40.0])
    sf = np.ones(5)
    assert wm.cut0_emulated(t, np.array([3.0, -3.0, 3.0, 0.0, 0.0]), sf)  # either sign counts
    assert not wm.cut0_emulated(t, np.array([3.0, 3.0, 0.0, 3.0, 3.0]), sf)  # gap > 8 d
    assert not wm.cut0_emulated(t, np.array([2.0, 2.0, 2.0, 0.0, 0.0]), sf)  # S/N ≤ 2.7


def test_feature_window_covers_the_umbra_crossing():
    r = {"t0": 100.0, "tE": 10.0, "u0": 1.2, "rho": 0.01}
    lo, hi = wm.feature_window(r, "N1neg")
    half = 10.0 * np.sqrt(4.0 - 1.44)
    assert lo < 100.0 - half and hi > 100.0 + half
    assert wm.notch_in_window(*_dip_lc(), 95.0, 105.0) < -5.0


def _dip_lc():
    t = np.linspace(0.0, 200.0, 801)
    f, sf = _noise(t)
    f[(t > 95) & (t < 105)] -= 500.0
    return t, f, sf


def test_luminosity_function_weights():
    w = wm.lf_weights([15.0, 18.0, 21.0])
    assert np.all(np.diff(w) > 0)
    assert w[2] / w[1] == pytest.approx(10 ** (3 * wm.LF_SLOPE))
    assert 0.9 < wm.lf_fraction_injected() < 1.0


def test_neighbours_and_fit_rows_round_trip():
    ev = Table(
        {
            "event_id": ["a", "b", "c", "d"],
            "chip": [1, 1, 1, 2],
            "subframe": [0, 0, 0, 0],
            "x": [10.0, 15.0, 40.0, 10.0],
            "y": [10.0, 12.0, 10.0, 10.0],
        }
    )
    assert wm.neighbours_of(ev, "a") == ["b"]
    row = {"event_id": "a", "n_points": 100}
    for m, bic in (("PSPL", 50.0), ("N1neg", 20.0), ("E2pos", 60.0), ("E2neg", 30.0)):
        row.update(
            {f"{m}_chi2": bic - 10, f"{m}_dof": 95, f"{m}_bic": bic, f"{m}_t0": 1.0, f"{m}_tE": 5.0}
        )
        row.update({f"{m}_u0": 0.5, f"{m}_fs": 1.0, f"{m}_fb": 0.0})
        if m != "PSPL":
            row[f"{m}_rho"] = 0.01
    tab = Table([row])
    res = wm.res_from_row(tab[0])
    assert set(res) == {"PSPL", "N1neg", "E2pos", "E2neg"}
    assert res["N1neg"]["k"] == 5 and res["N1neg"]["rho"] == 0.01


def test_vetting_removes_a_star_with_repeated_dips():
    """A flagged deficit that repeats elsewhere is a variable star, not one lensing event."""
    t = _cadence(n=500)
    f, sf = _noise(t)
    sig = _w3(t, fs=3000.0)
    for c in (2454100.0, 2454850.0):
        sig[np.abs(t - c) < 15] -= 1500.0
    scan = wm.deficit_scan(t, f + sig, sf)
    scan["t_lo"], scan["t_hi"] = 2454455.0, 2454505.0
    ev = {"event_id": "synthetic", "ra": 279.1, "dec": -23.8}
    out = wm.vet_one((ev, t, f + sig, sf, {}, scan, [], False, None))
    names = [n for n, _ok, _ in out["tests"]]
    assert names[0] == "screen_flag" and out["tests"][0][1]
    assert out["exotic"] in ("N1neg", "E2neg")
    assert ("repeated_deficit", False) in [(n, ok) for n, ok, _ in out["tests"]]
    assert not out["survives"] and not out["complete"]


def test_a_deficit_shared_by_many_objects_is_a_frame_artefact():
    rng = np.random.default_rng(5)
    n = 300
    centres = rng.uniform(2453824.0, 2456970.0, n)
    centres[:12] = 2455000.0 + rng.uniform(-0.5, 0.5, 12)  # twelve objects dim on one night
    pre = Table(
        {
            "event_id": [f"gb22-R-{1 + i % 10}-0-{i}" for i in range(n)],
            "z_min": np.full(n, -12.0),
            "s_min": np.full(n, -20.0),
            "z_min2": np.zeros(n),
            "width": np.full(n, 1.0),
            "t_lo": centres - 0.5,
            "t_hi": centres + 0.5,
            "error": [""] * n,
        }
    )
    pop = wm.deficit_population(pre)
    hot = wm.shared_epoch(pop, 1.0, 2455000.0, "gb22-R-1-0-0")
    assert hot["n"] >= 11 and hot["p"] < wm.COINC_P
    lone = wm.shared_epoch(pop, 1.0, float(centres[100]), "gb22-R-1-0-100")
    assert lone["p"] >= wm.COINC_P
    kept = set(wm.passes(pre)["event_id"])
    assert "gb22-R-1-0-0" not in kept and "gb22-R-1-0-100" in kept


def test_a_flat_bottomed_dip_is_fitted_as_an_eclipse():
    assert wm.trapezoid([0.0, 4.0, 10.0], 0.0, 10.0, 0.2).tolist() == [1.0, 1.0, 0.0]
    t, f, sf = _dip_lc()  # 500-count box dip, 95 < t < 105
    lc = wm.to_lightcurve(t, f, sf)
    ecl = wm.fit_eclipse(lc, 93.0, 107.0)
    assert ecl["tc"] == pytest.approx(100.0, abs=1.0)
    assert ecl["duration"] == pytest.approx(10.0, rel=0.2)
    flat = float(np.sum((f - np.average(f, weights=lc.w)) ** 2 * lc.w))
    assert ecl["chi2"] < flat - 100.0
