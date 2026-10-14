from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location(
    "field_followups", ROOT / "scripts" / "field_followups.py"
)
ff = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ff)


def test_dja_photometry_nearest_and_bands():
    t = Table(
        {
            "id": [1, 2],
            "ra": [10.0, 10.001],
            "dec": [-5.0, -5.0],
            "f200w_flux_aper_1": [1.0, 2.0],  # uJy: 1 uJy = 23.9 AB
            "f200w_fluxerr_aper_1": [0.1, 0.1],
            "f444w_flux_aper_1": [-0.2, 1.0],
            "f444w_fluxerr_aper_1": [0.1, 0.1],
            "f090w_flux_aper_1": [np.nan, 1.0],
            "f090w_fluxerr_aper_1": [0.1, 0.1],
        }
    )
    out = ff.dja_photometry(t, 10.0, -5.0 + 0.1 / 3600)
    assert out["id"] == 1 and out["sep_arcsec"] == 0.1
    assert out["bands"]["f200w"] == {"mag": 23.9, "snr": 10.0}
    assert out["bands"]["f444w"] == {"mag": None, "snr": -2.0}
    assert "f090w" not in out["bands"]  # NaN flux is left out


def test_spike_sources_row_and_column():
    t = Table(
        {
            "label": [1, 2, 3, 4, 5],
            "xcentroid": [100.0, 300.0, 103.0, 300.0, 100.0],
            "ycentroid": [100.0, 102.0, 400.0, 300.0, 2000.0],
            "aper50_abmag": [29.0, 16.5, 18.0, 15.0, 15.0],
        }
    )
    out = ff.spike_sources(t, 1, 0.03)
    # 2 on the same row, 3 on the same column; 4 is off both axes; 5 is beyond the reach
    assert [(s["axis"], s["label"]) for s in out] == [("row", 2), ("column", 3)]
    assert out[0]["sep_arcsec"] == round(np.hypot(200, 2) * 0.03, 2)


def test_aperture_mag_point_source():
    n = 101
    pix = 0.03
    jj, ii = np.mgrid[:n, :n].astype(float)
    xx, yy = (ii - 50) * pix, (jj - 50) * pix
    sci = np.full((n, n), 0.5)  # MJy/sr background
    sci[50, 50] += 1000.0
    err = np.full((n, n), 0.01)
    out = ff.aperture_mag(sci, err, xx, yy)
    flux_jy = 1000.0 * (pix / 206264.806) ** 2 * 1e6
    assert out["mag"] == round(-2.5 * np.log10(flux_jy / 3631.0), 2)
    assert out["snr"] > 100


def test_centroid_and_relative_shift():
    n = 61
    pix = 0.03
    jj, ii = np.mgrid[:n, :n].astype(float)
    xx, yy = (ii - 30) * pix, (jj - 30) * pix
    sci = np.exp(-((xx - 0.012) ** 2 + (yy + 0.006) ** 2) / (2 * 0.04**2))
    x0, y0 = ff.centroid(sci, xx, yy)
    assert abs(x0 - 0.012) < 1e-3 and abs(y0 + 0.006) < 1e-3

    rng = np.random.default_rng(1)
    refs = np.array([17.0, -6.0]) + rng.normal(0, 5, (200, 2))
    refs[3] = np.nan
    out = ff.relative_shift(np.array([37.0, -6.0]), refs)
    assert out["n_ref"] == 199
    assert abs(out["shift_mas"][0] - 20.0) < 2.0 and abs(out["shift_mas"][1]) < 2.0
    assert 3.0 < out["shift_sigma"][0] < 5.5


def test_dja_photometry_no_match_beyond_radius():
    t = Table(
        {
            "id": [1],
            "ra": [10.0],
            "dec": [-5.0],
            "f200w_flux_aper_1": [1.0],
            "f200w_fluxerr_aper_1": [0.1],
        }
    )
    out = ff.dja_photometry(t, 10.0, -5.0 + 2.0 / 3600)
    assert out["id"] is None and out["nearest_sep_arcsec"] == 2.0 and out["bands"] == {}


def test_spike_sources_keeps_saturated_nan_magnitude():
    t = Table(
        {
            "label": [1, 2],
            "xcentroid": [100.0, 400.0],
            "ycentroid": [100.0, 101.0],
            "aper50_abmag": [25.0, np.nan],
        }
    )
    out = ff.spike_sources(t, 1, 0.03)
    assert [(s["axis"], s["label"], s["aper50_abmag"]) for s in out] == [("row", 2, None)]


def test_aperture_mag_all_nan_stamp():
    xx, yy = np.meshgrid(np.arange(-5, 6) * 0.03, np.arange(-5, 6) * 0.03)
    nan = np.full(xx.shape, np.nan)
    assert ff.aperture_mag(nan, nan, xx, yy) == {"mag": None, "snr": None}


def test_drop_self_removes_own_match_and_flags_near_source():
    cps = [
        {"sep_from_source_arcsec": 1.3, "f200w_match": 7, "f200w_match_mag": 24.1},
        {"sep_from_source_arcsec": 12.0, "f200w_match": 9, "f200w_match_mag": 25.0},
    ]
    out = ff.drop_self(cps, {"label": 7})
    assert out[0]["near_source"] and out[0]["f200w_match"] is None
    assert not out[1]["near_source"] and out[1]["f200w_match"] == 9
