"""Offline tests for jwst_anomaly.gaia_mulens and scripts/w3_gaia.py (D-061)."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from jwst_anomaly import gaia_mulens as gm

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w3_gaia", _DIR / "w3_gaia.py")
wg = importlib.util.module_from_spec(_spec)
sys.modules["w3_gaia"] = wg
_spec.loader.exec_module(wg)
w3 = wg.w3

HEADER = (
    "solution_id,source_id,transit_id,g_transit_time,g_transit_flux,g_transit_flux_error,"
    "g_transit_flux_over_error,g_transit_mag,bp_obs_time,bp_flux,bp_flux_error,bp_flux_over_error,"
    "bp_mag,rp_obs_time,rp_flux,rp_flux_error,rp_flux_over_error,rp_mag,variability_flag_g_reject,"
    "variability_flag_bp_reject,variability_flag_rp_reject,g_other_flags,bp_other_flags,"
    "rp_other_flags,rejected_by_photometry"
)


def _row(t, mag, foe, reject="false", rej_phot="false"):
    return f"1,42,7,{t},1,1,{foe},{mag},,,,,,,,,,,{reject},false,false,1,0,0,{rej_phot}"


def test_parse_epoch_csv_drops_rejected_rows_and_converts_time_and_errors():
    text = "\n".join(
        [
            HEADER,
            _row(1800.5, 17.0, 100.0),
            _row(1700.25, 17.1, 50.0),
            _row(1750.0, 15.0, 80.0, reject="true"),
            _row(1760.0, 15.0, 80.0, rej_phot="true"),
            _row("", "", ""),
        ]
    )
    lc = gm.parse_epoch_csv(text)
    assert len(lc) == 2
    assert lc.meta["n_rejected"] == 3
    np.testing.assert_allclose(lc["time"], [1700.25 + gm.GAIA_T0, 1800.5 + gm.GAIA_T0])
    np.testing.assert_allclose(lc["mag_err"], 2.5 / math.log(10) / np.array([50.0, 100.0]))
    assert set(lc["band"]) == {"G"}


def test_parse_method_table_reads_the_last_column():
    tex = (
        "001 & 6059400613544951552 & 184.43 & -59.02 & 298.6 & 3.5 & 13.58 & A \\\\ \n"
        "002 & 4135677580466393216 & 258.07 & -16.86 & 6.0 & 12.9 & 13.64 & A+B \\\\ \n"
        "003 & 4084984459435880064 & 284.43 & -20.77 & 14.9 & -10.6 & 14.16 & B \\\\ \n"
        "\\hline\n"
    )
    assert gm.parse_method_table(tex) == {
        6059400613544951552: "A",
        4135677580466393216: "A+B",
        4084984459435880064: "B",
    }


def test_rescale_g_errors_follows_eq_9_10():
    e = wg.rescale_g_errors([12.0, 19.0], [0.0, 0.01])
    floor = math.sqrt(30.0) * 10 ** (0.17 * 13.5 - 5.1)
    assert e[0] == floor  # brighter than 13.5 uses 13.5
    assert math.isclose(e[1], math.hypot(0.01, math.sqrt(30.0) * 10 ** (0.17 * 19.0 - 5.1)))


def test_abbe_is_small_for_a_smooth_bump_and_near_one_for_noise():
    rng = np.random.default_rng(3)
    assert abs(wg.abbe(rng.normal(size=4000)) - 1.0) < 0.1
    x = np.exp(-(np.linspace(-3, 3, 200) ** 2))
    assert wg.abbe(x) < 0.01


def _gaia_like(a, fs=1.0, seed=4, n=120):
    rng = np.random.default_rng(seed)
    t = np.sort(gm.T_FIRST + rng.uniform(0.0, gm.T_LAST - gm.T_FIRST, n))
    flux = fs * a(t)
    mag = w3.flux_to_mag(flux)
    err = wg.rescale_g_errors(mag, np.full(n, 0.003))
    mag = mag + rng.normal(0.0, err)
    return w3.LightCurve.from_mag(t, mag, err, event_id="synthetic")


def test_level0_recovers_a_pspl_and_the_emulated_selection_passes_it():
    t0, te, u0 = gm.T_FIRST + 500.0, 60.0, 0.15
    lc = _gaia_like(lambda t: w3.pspl(w3.straight_beta(t, t0, te, u0)), fs=3.0)  # baseline G ≈ 16.8
    l0 = wg.level0_fit(lc, t0 + 5.0, 40.0, 0.3)
    assert abs(l0["t0"] - t0) < 3.0 and abs(l0["tE"] / te - 1) < 0.1 and abs(l0["u0"] - u0) < 0.03
    sel = wg.sample_a_selection(lc, l0, None)
    for cut in ("u0", "te", "g0", "skew", "amp", "t_first", "par_chi2", "abbe"):
        assert sel["passed"][cut], cut


def test_a_pure_dimming_fails_the_skewness_cut():
    t0, te = gm.T_FIRST + 500.0, 60.0

    def dim(t):
        return 1.0 - 0.7 * np.exp(-0.5 * ((t - t0) / te) ** 2)

    lc = _gaia_like(dim)
    l0 = wg.level0_fit(lc, t0, te, 0.3)
    sel = wg.sample_a_selection(lc, l0, None)
    assert sel["skew"] > 0 and not sel["passed"]["skew"] and not sel["selected"]


@pytest.mark.network
def test_live_gaia_events_methods_and_one_light_curve(tmp_path):
    survey = gm.GaiaDR3Microlensing(cache=tmp_path)
    ev = survey.events()
    assert len(ev) == 363
    counts = {m: int(np.sum(ev["method"] == m)) for m in ("A", "B", "A+B")}
    assert counts == {"A": 130, "B": 200, "A+B": 33}  # Wyrzykowski et al. 2023 §4.2
    lc = survey.light_curve(str(ev["event_id"][0]))
    assert len(lc) >= 10 and gm.T_FIRST - 1 <= lc["time"].min() <= lc["time"].max() <= gm.T_LAST + 1


def test_fetch_batch_handles_zip_and_bare_csv(tmp_path, monkeypatch):
    import io
    import zipfile

    csv_text = "\n".join([HEADER, _row(1800.5, 17.0, 100.0)]).encode()
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("EPOCH_PHOTOMETRY-Gaia DR3 111111111.csv", csv_text)
        z.writestr("EPOCH_PHOTOMETRY-Gaia DR3 222222222.csv", csv_text)
    replies = iter([buf.getvalue(), csv_text])

    class _Resp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(gm.urllib.request, "urlopen", lambda *a, **k: _Resp(next(replies)))
    survey = gm.GaiaDR3Microlensing(cache=tmp_path)
    survey._fetch_batch(["111111111", "222222222"])
    survey._fetch_batch(["333333333"])
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "111111111.csv",
        "222222222.csv",
        "333333333.csv",
    ]
    assert len(survey.light_curve("333333333")) == 1


def test_strict_skew_abbe_cut_applies_in_the_first_mission_year(monkeypatch):
    """Table C.1's window (6824.5–7189.5) falls in the DR3 data only as JD − 2450000."""
    t0, te, u0 = gm.T_FIRST + 500.0, 60.0, 0.15
    lc = _gaia_like(lambda t: w3.pspl(w3.straight_beta(t, t0, te, u0)), fs=3.0)
    l0 = wg.level0_fit(lc, t0 + 5.0, 40.0, 0.3)
    # skewness −1, Abbe 0.2: passes log10(abbe) < −0.6 (loose), fails < −0.84 (strict)
    monkeypatch.setattr(wg, "skew", lambda x: -1.0)
    monkeypatch.setattr(wg, "abbe", lambda x: 0.2)
    inside = {**l0, "t0": wg.S.abbe_window_zero + 7000.0}
    outside = {**l0, "t0": wg.S.abbe_window_zero + 7500.0}
    assert gm.T_FIRST <= inside["t0"] <= gm.T_LAST and gm.T_FIRST <= outside["t0"] <= gm.T_LAST
    assert not wg.sample_a_selection(lc, inside, None)["passed"]["abbe"]
    assert wg.sample_a_selection(lc, outside, None)["passed"]["abbe"]
