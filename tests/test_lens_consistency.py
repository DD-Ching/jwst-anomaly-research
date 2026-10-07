"""Tests for scripts/lens_consistency.py (synthetic model, map and catalogs)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from jwst_anomaly.lensmodel import DPIE, LensModel

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("lens_consistency", _DIR / "lens_consistency.py")
lc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lc)

RA0, DEC0 = 110.826989, -73.454723


def _model():
    comp = DPIE("halo", 0.0, 0.0, 0.0, 0.0, 1e-3, 1e4, 900.0, 0.39)  # SIS-like, theta_E(D=1) ~ 23"
    return LensModel([comp], RA0, DEC0, FlatLambdaCDM(H0=70.0, Om0=0.3), source="synthetic")


def test_kappa_map_check_on_own_map(tmp_path):
    model = _model()
    n, scale = 200, 0.1 / 3600
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [RA0, DEC0]
    w.wcs.crpix = [100.5, 100.5]
    w.wcs.cdelt = [-scale, scale]
    jj, ii = np.mgrid[0:n, 0:n]
    ra, dec = w.pixel_to_world_values(ii, jj)
    x, y = model.to_frame(ra, dec)
    path = tmp_path / "k.fits"
    fits.PrimaryHDU(model.kappa_xy(x, y).astype(np.float32), header=w.to_header()).writeto(path)
    res = lc.kappa_map_check(model, path, step=4)
    assert res["n_points"] == 50 * 50
    assert res["median_abs_dkappa"] < 1e-6  # float32 rounding only
    assert abs(res["median_ratio"] - 1) < 1e-6
    assert -10.1 < res["extent_arcsec"][0] < -9.8


def test_chi2pos_from_par(tmp_path):
    par = tmp_path / "best.par"
    par.write_text("#Chi2tot(dof=32): 30.9\n#Chi2pos: 30.913219\nrunmode\n end\n")
    assert lc.chi2pos_from_par(par) == 30.913219
    par.write_text("runmode\n end\n")
    assert lc.chi2pos_from_par(par) is None


def test_orientation_classes():
    model = _model()
    # Three sources 40" North of the centre, where arcs should run East-West (PA 90).
    dec = DEC0 + 40.0 / 3600
    src = Table(
        {
            "label": [1, 2, 3],
            "ra": [RA0, RA0, RA0],
            "dec": [dec, dec, dec],
            "pa_obs": [92.0, 0.0, 45.0],
            "z_phot": [2.0, 3.0, np.nan],
        }
    )
    out = lc.orientation_table(model, src)
    assert list(out["orientation_class"]) == ["aligned", "anti", "mixed"]
    np.testing.assert_allclose(out["offset_z2"], [2.0, 90.0, 45.0], atol=1e-6)
    np.testing.assert_allclose(out["offset_zphot"][:2], [2.0, 90.0], atol=1e-6)
    assert np.isnan(out["offset_zphot"][2])
    assert np.all(out["pa_pred_spread"] < 1e-6)
    assert out.meta["provenance"] == "derived"


def test_load_shapes_and_photoz(tmp_path):
    cat = Table(
        {
            "label": [1, 2],
            "sky_centroid": SkyCoord([RA0, RA0 + 1e-3], [DEC0, DEC0], unit="deg"),
            "sky_orientation": [190.0, 45.0],
            "ellipticity": [0.6, 0.2],
            "semimajor_sigma": [3.0, 1.0],
            "isophotal_area": [50.0, 10.0],
            "isophotal_flux": [100.0, 5.0],
            "isophotal_flux_err": [2.0, 1.0],
            "is_extended": [True, False],
        }
    )
    path = tmp_path / "x_cat.ecsv"
    cat.write(path)
    shapes = lc.load_shapes(path)
    np.testing.assert_allclose(shapes["pa_obs"], [10.0, 45.0])
    np.testing.assert_allclose(shapes["snr"], [50.0, 5.0])
    assert shapes.meta["provenance"] == "observed"
    zout = Table(
        {
            "id": [7],
            "ra": [RA0 + 1e-5],
            "dec": [DEC0],
            "z_phot": [1.5],
            "z160": [1.4],
            "z840": [1.6],
        }
    )
    zpath = tmp_path / "zout.fits"
    zout.write(zpath)
    lc.attach_photoz(shapes, zpath)
    assert shapes["z_phot"][0] == 1.5 and np.isnan(shapes["z_phot"][1])  # second is 13" away
