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


def test_feature_coverage_finds_the_unobserved_caustic_spike():
    """An exotic fit whose spike falls in an observing gap has no epochs in its feature."""
    t = np.sort(np.concatenate([_cadence(n=300, span=300.0), _cadence(n=300, span=300.0) + 400]))
    lc = _synthetic(t, np.ones_like(t), sigma=0.003)  # flat; the gap is 2456300 … 2456400
    flat = {"t0": 2456350.0, "tE": 1.0, "u0": 50.0, "fs": 0.6, "fb": 0.4}
    in_gap = {"t0": 2456350.0, "tE": 2.0, "u0": 1.0, "rho": 0.001, "fs": 0.6, "fb": 0.4}
    cov = w3.feature_coverage(lc, flat, in_gap, "PSPL", "N1neg")
    assert cov["max_diff_sigma"] > 3 and cov["duration_days"] > 0
    assert cov["n_epochs"] == 0  # the umbra and both spikes sit inside the 100-day gap
    on_data = {**in_gap, "t0": 2456280.0}
    assert w3.feature_coverage(lc, flat, on_data, "PSPL", "N1neg")["n_epochs"] > 3


def test_jackknife_removes_a_preference_built_on_one_epoch():
    t = _cadence(n=400)
    a = w3.pspl(w3.straight_beta(t, 2456800.0, 25.0, 0.2))
    lc = _synthetic(t, a, sigma=0.004)
    i = int(np.argmin(np.abs(lc.t - 2456800.0)))
    lc.f[i] *= 0.5  # one bad measurement at the peak: an exotic dip would love it
    res = w3.fit_event(lc, 2456800.0, 25.0, 0.2, models=("PSPL", "N1neg"))
    jk = w3.jackknife_worst_epochs(lc, res["PSPL"], res["N1neg"], "PSPL", "N1neg")
    assert jk[-1] > res["N1neg"]["bic"] - res["PSPL"]["bic"]  # the preference weakens
    assert jk[-1] > w3.P.flag_dbic  # and no longer flags


def test_jackknife_n_drop_keeps_three_feature_epochs():
    assert [w3.jackknife_n_drop(n) for n in (0, 3, 4, 5, 6, 50)] == [0, 0, 1, 2, 3, 3]


def test_well_sampled_short_w3_event_survives_feature_coverage_and_jackknife():
    """A real t_E = 3 d W3 event sampled by a few epochs must not be vetted away (PR #88 review)."""
    t = _cadence(n=500, span=1600.0)
    truth = {"t0": 2456800.0, "tE": 3.0, "u0": 0.2, "rho": 0.01, "fs": 0.6, "fb": 0.4}
    grid = w3.LightCurve(t, np.ones_like(t), np.ones_like(t), RA, DEC, "x")
    lc = _synthetic(t, (w3.model_flux("N1neg", grid, truth) - 0.4) / 0.6)
    res = w3.fit_event(lc, 2456800.0, 3.0, 0.2, models=("PSPL", "N1neg"))
    assert res["N1neg"]["bic"] - res["PSPL"]["bic"] < w3.P.flag_dbic
    cov = w3.feature_coverage(lc, res["PSPL"], res["N1neg"], "PSPL", "N1neg")
    assert cov["n_epochs"] >= 3 and cov["dchi2_in"] < w3.P.flag_dbic
    n_drop = w3.jackknife_n_drop(cov["n_epochs"])
    jk = w3.jackknife_worst_epochs(lc, res["PSPL"], res["N1neg"], "PSPL", "N1neg", n_drop=n_drop)
    assert jk[-1] < w3.P.flag_dbic
    d2e, _ = w3.two_events_dbic(lc, res["PSPL"], res["N1neg"], "PSPL", "N1neg")
    assert d2e is None or d2e < w3.P.flag_dbic


def test_two_unrelated_bumps_explain_spikes_far_apart():
    """An event plus an unrelated later flare: two PSPL bumps beat any single exotic fit."""
    t = _cadence(n=600, span=1600.0)
    t = np.sort(np.concatenate([t, 2457300.0 + np.array([-0.2, 0.0, 0.3, 0.9])]))
    a1 = w3.pspl(w3.straight_beta(t, 2456300.0, 3.0, 0.5))
    a2 = w3.pspl(w3.straight_beta(t, 2457300.0, 1.0, 0.3))
    lc = _synthetic(t, a1 + 0.3 * a2, sigma=0.005)
    res = w3.fit_event(lc, 2456300.0, 3.0, 0.5, models=("PSPL", "N1neg"))
    d2e, two = w3.two_events_dbic(lc, res["PSPL"], res["N1neg"], "PSPL", "N1neg")
    assert d2e is not None and d2e > w3.P.flag_dbic
    assert abs(two["x"][3] - 2457300.0) < 2.0


@pytest.mark.parametrize(
    "flag",
    [
        {"event_id": "A", "tests": [], "survives": True},
        {"event_id": "B", "tests": [], "survives": False, "complete": False},
        {"event_id": "C", "tests": [], "survives": False},  # older file: completeness unknown
    ],
)
def test_limit_refuses_a_zero_event_limit_unless_vetting_is_a_complete_null(
    tmp_path, monkeypatch, flag
):
    import json

    doc = {"fit_chunk": "", "flags": [flag]}
    (tmp_path / "vetting_bulge2019.json").write_text(json.dumps(doc))
    monkeypatch.setattr(w3, "out_dir", lambda: tmp_path)
    with pytest.raises(SystemExit, match=flag["event_id"]):
        w3.run_limit()


def test_limit_passes_a_complete_null_vetting_to_the_injections(tmp_path, monkeypatch):
    import json

    flag = {"event_id": "A", "tests": [], "survives": False, "complete": True}
    (tmp_path / "vetting_bulge2019.json").write_text(json.dumps({"fit_chunk": "", "flags": [flag]}))
    monkeypatch.setattr(w3, "out_dir", lambda: tmp_path)
    with pytest.raises(FileNotFoundError, match="injections_bulge2019"):
        w3.run_limit()


@pytest.mark.parametrize("doc", [{"fit_chunk": "1/12"}, {}])  # one chunk, or an older file
def test_limit_refuses_vetting_of_a_partial_fit(tmp_path, monkeypatch, doc):
    import json

    (tmp_path / "vetting_bulge2019.json").write_text(json.dumps({**doc, "flags": []}))
    monkeypatch.setattr(w3, "out_dir", lambda: tmp_path)
    with pytest.raises(SystemExit, match="fit chunk"):
        w3.run_limit()


def test_limit_without_vetting_says_so(tmp_path, monkeypatch):
    monkeypatch.setattr(w3, "out_dir", lambda: tmp_path)
    with pytest.raises(SystemExit, match="run `vet` first"):
        w3.run_limit()


def test_parallax_objective_rejects_pi_e_beyond_the_bound():
    t = _cadence()
    lc = _synthetic(t, w3.pspl(w3.straight_beta(t, 2456800.0, 40.0, 0.2)))
    fun = w3._objective("PAR", lc, 2456800.0)
    lt = np.log10(40.0)
    inside = 0.6 * w3.P.pie_max
    assert fun(np.array([2456800.0, lt, 0.2, w3.P.pie_max, 0.1])) == 1e30
    pytest.importorskip("MulensModel")  # the parallax trajectory inside the bound needs it
    assert fun(np.array([2456800.0, lt, 0.2, inside, inside])) < 1e30


def test_fit_drops_checkpoint_rows_fitted_with_other_params(tmp_path, monkeypatch):
    import dataclasses
    import json
    import types

    monkeypatch.setattr(w3, "out_dir", lambda: tmp_path)
    fake = types.SimpleNamespace(
        name="fake", spec=types.SimpleNamespace(reference="x"), events=lambda: None
    )
    monkeypatch.setattr(w3.ogle, "OgleMrozSample", lambda key: fake)
    monkeypatch.setattr(w3, "_jobs", lambda sample, limit, skip, chunk=None: [])
    monkeypatch.setattr(w3, "_join_pub", lambda tab, ev: tab)
    row = {"event_id": "new", "error": "", "best_ordinary": "PSPL", "seconds": 1.0}
    old = {**row, "event_id": "old"}  # written before checkpoints carried a params tag
    (tmp_path / "fits_k.partial.jsonl").write_text(
        json.dumps(old) + "\n" + json.dumps({**row, "params_tag": w3.params_tag()}) + "\n"
    )
    tab = w3.Table.read(w3.run_fit("k", None, 1))
    assert list(tab["event_id"]) == ["new"] and "params_tag" not in tab.colnames
    monkeypatch.setattr(w3, "P", dataclasses.replace(w3.P, pie_max=1e4))
    assert (
        w3.params_tag()
        != json.loads((tmp_path / "fits_k.partial.jsonl").read_text().splitlines()[0])["params_tag"]
    )


def test_parse_chunk_is_one_based_and_chunks_partition_the_sample():
    assert w3.parse_chunk(None) is None
    assert w3.parse_chunk("1/6") == (0, 6)
    with pytest.raises(ValueError):
        w3.parse_chunk("0/6")
    with pytest.raises(ValueError):
        w3.parse_chunk("7/6")
    idx = np.arange(23)
    parts = [idx[k::6] for k, _ in (w3.parse_chunk(f"{i}/6") for i in range(1, 7))]
    assert sorted(np.concatenate(parts).tolist()) == idx.tolist()


def test_chunked_fit_keeps_only_its_chunk_and_writes_a_chunk_table(tmp_path, monkeypatch):
    import json
    import types

    monkeypatch.setattr(w3, "out_dir", lambda: tmp_path)
    (tmp_path / "res").mkdir()
    monkeypatch.setattr(w3, "results_dir", lambda: tmp_path / "res")
    ev = w3.Table({"event_id": ["a", "b", "c", "d"]})
    fake = types.SimpleNamespace(
        name="fake", spec=types.SimpleNamespace(reference="x"), events=lambda: ev
    )
    monkeypatch.setattr(w3.ogle, "OgleMrozSample", lambda key: fake)
    monkeypatch.setattr(w3, "_jobs", lambda sample, limit, skip, chunk=None: [])
    monkeypatch.setattr(w3, "_join_pub", lambda tab, ev: tab)
    row = {"error": "", "best_ordinary": "PSPL", "seconds": 1.0, "params_tag": w3.params_tag()}
    (tmp_path / "fits_k.partial.jsonl").write_text(
        "".join(json.dumps({**row, "event_id": e}) + "\n" for e in "abcd")
    )
    tab = w3.Table.read(w3.run_fit("k", None, 1, chunk=(1, 2)))
    assert sorted(tab["event_id"]) == ["b", "d"]
    assert sorted(w3.Table.read(tmp_path / "fits_k_chunk2of2.ecsv")["event_id"]) == ["b", "d"]
    gz = w3.Table.read(tmp_path / "res" / "fits_k_chunk2of2.ecsv.gz", format="ascii.ecsv")
    assert sorted(gz["event_id"]) == ["b", "d"]  # the tracked copy (D-059)
    assert len((tmp_path / "fits_k.partial.jsonl").read_text().splitlines()) == 4


def _chunk_setup(tmp_path, monkeypatch, n_events=5):
    import types

    monkeypatch.setattr(w3, "out_dir", lambda: tmp_path / "out")
    (tmp_path / "out").mkdir()
    monkeypatch.setattr(w3, "results_dir", lambda: tmp_path)
    ev = w3.Table({"event_id": [f"e{i}" for i in range(n_events)]})
    fake = types.SimpleNamespace(events=lambda: ev)
    monkeypatch.setattr(w3.ogle, "OgleMrozSample", lambda key: fake)
    return ev


def _write_chunk(tmp_path, ids, k, n, params=None):
    import json
    from dataclasses import asdict

    tab = w3.Table({"event_id": ids, "seconds": [1.0] * len(ids), "error": [""] * len(ids)})
    tab.meta.update(params=params or json.dumps(asdict(w3.P)), chunk=f"{k}/{n}")
    w3.write_ecsv_gz(tab, tmp_path / f"fits_k_chunk{k}of{n}.ecsv.gz")


def test_merge_chunks_joins_all_chunks_into_a_whole_sample_table(tmp_path, monkeypatch):
    ev = _chunk_setup(tmp_path, monkeypatch)
    for k in (1, 2):
        _write_chunk(tmp_path, list(ev["event_id"][k - 1 :: 2]), k, 2)
    a = (tmp_path / "fits_k_chunk1of2.ecsv.gz").read_bytes()
    _write_chunk(tmp_path, list(ev["event_id"][0::2]), 1, 2)
    assert (tmp_path / "fits_k_chunk1of2.ecsv.gz").read_bytes() == a  # deterministic bytes
    tab = w3.Table.read(w3.merge_chunks("k", 2))
    assert sorted(tab["event_id"]) == list(ev["event_id"]) and tab.meta["chunk"] == ""
    assert tab.meta["cpu_time_s"] == 5.0


def test_merge_chunks_refuses_missing_foreign_or_misplaced_chunks(tmp_path, monkeypatch):
    ev = _chunk_setup(tmp_path, monkeypatch)
    _write_chunk(tmp_path, list(ev["event_id"][0::2]), 1, 2)
    with pytest.raises(SystemExit, match=r"missing chunks of 2: \[2\]"):
        w3.merge_chunks("k", 2)
    _write_chunk(tmp_path, list(ev["event_id"][1::2]), 2, 2, params='{"other": 1}')
    with pytest.raises(SystemExit, match="other Params"):
        w3.merge_chunks("k", 2)
    _write_chunk(tmp_path, ["e0", "e1", "e3"], 2, 2)  # e0 belongs to chunk 1
    with pytest.raises(SystemExit, match="1 from outside"):
        w3.merge_chunks("k", 2)
    _write_chunk(tmp_path, ["e1", "e3", "e3"], 2, 2)
    with pytest.raises(SystemExit, match="duplicate"):
        w3.merge_chunks("k", 2)
    _write_chunk(tmp_path, ["e1"], 2, 2)  # e3 skipped (no light curve): not the whole sample
    with pytest.raises(SystemExit, match="1 events of chunk 2/2 missing"):
        w3.merge_chunks("k", 2)
