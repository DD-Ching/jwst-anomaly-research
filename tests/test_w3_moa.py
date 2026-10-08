"""Offline tests for scripts/w3_moa.py on synthetic difference light curves (D-062)."""

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


def test_merge_chunks_joins_complete_chunks_and_refuses_bad_ones(tmp_path, monkeypatch):
    ids = [f"gb22-R-1-0-{i}" for i in range(7)]
    pre = Table({"event_id": ids, "z_min": [-20.0] * 7, "s_min": [-9.0] * 7, "z_min2": [0.0] * 7})
    pre["width"] = [1.0] * 7
    pre["t_lo"] = 2454000.0 + 300.0 * np.arange(7)  # deficits at unrelated epochs
    pre["t_hi"] = pre["t_lo"] + 1.0
    pre["error"] = ["x"] + [""] * 6  # a pre-screen error is not a pass
    monkeypatch.setattr(wm, "out_dir", lambda: tmp_path)
    monkeypatch.setattr(wm, "results_dir", lambda: tmp_path)
    pre.write(tmp_path / "prescreen_gb22.ecsv")
    passes = ids[1:]

    def chunk(k, n, rows, params=None, pre_params=None):
        tab = Table({"event_id": rows, "dbic_min": [0.0] * len(rows)})
        tab.meta.update(
            fit_params=params or wm.json.dumps(wm.asdict(wm.w3.P)),
            params=pre_params or wm.json.dumps(wm.asdict(wm.P)),
            wall_time_s=1.0,
        )
        wm.w3.write_ecsv_gz(tab, tmp_path / f"{wm.chunk_name((k, n))}.gz")

    chunk(0, 2, passes[0::2])
    with pytest.raises(SystemExit, match="missing chunk 2/2"):
        wm.merge_chunks(2)
    chunk(1, 2, passes[1::2])
    out = Table.read(wm.merge_chunks(2))
    assert list(out["event_id"]) == sorted(passes) and out.meta["wall_time_s"] == 2.0
    chunk(0, 1, passes, params="{}")
    with pytest.raises(SystemExit, match="other Params"):
        wm.merge_chunks(1)
    chunk(0, 1, passes, pre_params="{}")  # other pre-screen Params: other deficit starts
    with pytest.raises(SystemExit, match="other Params"):
        wm.merge_chunks(1)
    chunk(0, 1, passes[:-1])
    with pytest.raises(SystemExit, match="exactly its pre-screen passes"):
        wm.merge_chunks(1)


def test_pool_workers_get_moa_bounds_and_populations(monkeypatch):
    # spawned workers (Windows) do not run main(): the initializer must set both
    monkeypatch.setattr(wm.w3, "P", wm.replace(wm.w3.P, te_bounds=(0.5, 50.0)))
    monkeypatch.setattr(wm, "_POP", None)
    monkeypatch.setattr(wm, "_BASELINE", None)
    monkeypatch.setattr(wm, "FIELD", "gb22")
    wm._init_worker({"pop": ("field", "chips"), "baseline": 3.5, "field": "gb5"})
    assert wm.w3.P.te_bounds == wm.MOA_FIT_BOUNDS["te_bounds"]
    assert wm._POP == ("field", "chips")
    assert wm.baseline_threshold() == 3.5 and wm.FIELD == "gb5"


def test_vet_and_limit_refuse_partial_or_failed_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(wm, "out_dir", lambda: tmp_path)
    ids = [f"gb22-R-1-0-{i}" for i in range(4)]
    pre = Table({"event_id": ids, "z_min": [-20.0] * 4, "s_min": [-9.0] * 4, "z_min2": [0.0] * 4})
    pre["error"] = [""] * 4
    pre["width"] = [1.0] * 4
    pre["t_lo"] = 2454000.0 + 300.0 * np.arange(4)  # deficits at unrelated epochs
    pre["t_hi"] = pre["t_lo"] + 1.0
    pre.write(tmp_path / "prescreen_gb22.ecsv")
    fits = Table({"event_id": ids[:2]})
    fits.meta["chunk"] = "1/2"
    with pytest.raises(SystemExit, match="chunk 1/2 only"):
        wm.check_complete_fits(fits)
    fits.meta["chunk"] = ""
    with pytest.raises(SystemExit, match="exactly the pre-screen passes"):
        wm.check_complete_fits(fits)
    wm.check_complete_fits(Table({"event_id": ids[::-1]}))  # complete: no error
    vet = {"flags": [], "fit_errors": ["gb22-R-1-0-3"]}
    (tmp_path / "vetting_gb22.json").write_text(wm.json.dumps(vet))
    with pytest.raises(SystemExit, match="passes without a fit"):
        wm.run_limit()
    vet["fit_errors"] = []
    (tmp_path / "vetting_gb22.json").write_text(wm.json.dumps(vet))
    Table({"kind": ["W3", "W3"], "error": ["", "boom"]}).write(tmp_path / "injections_gb22.ecsv")
    with pytest.raises(SystemExit, match="1 injections failed"):
        wm.run_limit()


def _scan_table(n=400, seed=7):
    """Pre-screen rows: quiet light curves with red-noise χ²/dof and a few deficits (simulated)."""
    rng = np.random.default_rng(seed)
    tab = Table(
        {
            "event_id": [f"gb21-R-{1 + i % 10}-0-{i}" for i in range(n)],
            "flux_kind": ["difference"] * n,
            "z_min": rng.uniform(-3.5, -0.5, n),
            "s_min": np.full(n, -1.0),
            "z_max": np.full(n, 2.0),
            "z_min2": np.zeros(n),
            "width": np.full(n, 10.0),
            "t_lo": rng.uniform(2453824.0, 2456900.0, n),
            "n_points": np.full(n, 2000),
            "err_scale": np.ones(n),
            "chi2_const": np.exp(rng.normal(0.5, 0.5, n)),
            "offset": 512 * np.arange(n),
            "size": np.full(n, 100),
            "error": np.array([""] * n, dtype="U200"),
        }
    )
    tab["t_hi"] = tab["t_lo"] + 10.0
    tab["z_min"][:20] = -12.0  # deficits: always tracked
    tab["s_min"][:20] = -9.0
    return tab


def test_baseline_calibration_is_a_quantile_of_the_quiet_light_curves():
    tab = _scan_table()
    cal = wm.calibrate_baseline(tab)
    quiet = np.asarray(tab["chi2_const"])[wm.is_quiet(tab)]
    assert cal["n_quiet"] == quiet.size == 380
    assert cal["threshold"] == pytest.approx(np.quantile(quiet, wm.BASELINE_Q))
    tab.meta["quiet_chi2_hist"] = wm.chi2_histogram(quiet)  # the merged-table path
    cal_h = wm.calibrate_baseline(tab)
    assert cal_h["threshold"] == pytest.approx(cal["threshold"], rel=0.02)
    assert cal_h["frac_above_2"] == pytest.approx(np.mean(quiet > 2.0), abs=0.01)
    with pytest.raises(SystemExit, match="too few"):
        few = tab[:60]
        few.meta.pop("quiet_chi2_hist")
        wm.calibrate_baseline(few)


def test_tracked_rows_keep_deficits_errors_and_a_quiet_sample():
    tab = _scan_table()
    tab["error"][25] = "boom"
    keep = wm.tracked_mask(tab)
    assert keep[:20].all() and keep[25]
    frac = keep[wm.is_quiet(tab)].mean()
    assert 0.15 < frac < 0.35  # 1 in QUIET_TRACK_MOD = 4
    assert np.array_equal(keep, wm.tracked_mask(tab))  # deterministic


def test_lf_sampling_gives_equal_weights():
    rng = np.random.default_rng(1)
    mags = wm.sample_magnitudes(rng, 20000, "lf")
    assert mags.min() >= wm.INJ_IS[0] and mags.max() <= wm.INJ_IS[1]
    w = wm.lf_weights(mags) / wm.sampling_density(mags, "lf")
    assert np.ptp(w) / w.mean() < 1e-9  # n_eff = n
    # the sample follows the LF: its density ratio over 3 mag is 10^(3 × slope)
    hist, _ = np.histogram(mags, [17.0, 17.5, 20.0, 20.5])
    assert hist[2] / hist[0] == pytest.approx(10 ** (3 * wm.LF_SLOPE), rel=0.15)
    u = wm.sample_magnitudes(rng, 10, "uniform")
    assert np.allclose(wm.sampling_density(u, "uniform"), 1 / np.ptp(wm.INJ_IS))


def test_merge_prescreen_checks_the_member_count(tmp_path, monkeypatch):
    monkeypatch.setattr(wm, "out_dir", lambda: tmp_path)
    monkeypatch.setattr(wm, "results_dir", lambda: tmp_path)
    monkeypatch.setattr(wm, "FIELD", "gb21")
    monkeypatch.setitem(wm.moa.CUT0_PER_FIELD, 21, 400)
    monkeypatch.setitem(wm.moa.TAR_BYTES, 21, 2 * wm.CHUNK_BYTES - 1)
    tab = _scan_table()
    rows = [dict(zip(tab.colnames, r, strict=True)) for r in tab]
    with pytest.raises(SystemExit, match="missing pre-screen chunk 1/2"):
        wm.merge_prescreen("gb21")
    wm._write_prescreen_chunk("gb21", 0, 2, rows[:150], {"wall_time_s": 1.0})
    assert wm._chunk_done("gb21", 0, 2) and not wm._chunk_done("gb21", 1, 2)
    wm._write_prescreen_chunk("gb21", 1, 2, rows[150:399], {"wall_time_s": 1.0})
    with pytest.raises(SystemExit, match="399 light curves streamed, metadata has 400"):
        wm.merge_prescreen("gb21")
    wm._write_prescreen_chunk("gb21", 1, 2, rows[150:], {"wall_time_s": 1.0})
    pre = Table.read(wm.merge_prescreen("gb21"))
    assert wm.n_light_curves(pre) == 400 and len(pre) < 400
    assert wm.calibrate_baseline(pre)["n_quiet"] == 380  # from the histogram, not the sample
    assert pre.meta["provenance"] == "derived" and pre.meta["source"]


def test_combined_limit_sums_star_years_times_efficiency(tmp_path, monkeypatch):
    monkeypatch.setattr(wm, "results_dir", lambda: tmp_path)
    for field, ns in (("gb21", 1e6), ("gb20", 3e6)):
        rows = [
            {
                "field": field,
                "tE_days": te,
                "rho": rho,
                "n_s": ns,
                "n_s_low": ns / 2,
                "years": 8.0,
                "eff_per_star": 0.1,
                "n_inj": 200,
                "n_recovered": 20,
            }  # fmt: skip
            for te in wm.INJ_TE
            for rho in wm.INJ_RHO
        ]
        Table(rows).write(tmp_path / f"limits_{field}.ecsv")
    out = Table.read(wm.run_combine())
    assert out["rate95_per_star_yr"][0] == pytest.approx(3.0 / (4e6 * 8.0 * 0.1))
    assert out["rate95_conservative"][0] == pytest.approx(2 * out["rate95_per_star_yr"][0])
    assert out["n_fields"][0] == 2 and out["n_inj"][0] == 400


def test_published_star_counts_are_used_where_they_exist():
    pub = wm.field_star_counts("gb21")
    assert pub["n_s"] == wm.moa.NUNOTA_NS[21][1] and "Table 1" in pub["n_s_source"]
    model = wm.field_star_counts("gb22")
    assert model["n_s_low"] < model["n_s"] < model["n_s_high"]
    assert wm.parse_chunk_list("1-3,7") == [0, 1, 2, 6] and wm.parse_chunk_list(None) is None
