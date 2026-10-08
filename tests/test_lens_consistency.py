"""Tests for scripts/lens_consistency.py (synthetic model, map and catalogs)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
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
    comp = DPIE(
        "halo", 0.0, 0.0, 0.0, 0.0, 1e-3, 1e4, 900.0, 0.39
    )  # SIS-like, b0 = theta_E(D=1) ~ 35"
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
    # Sources 40" North of the centre, where arcs should run East-West (PA 90) at any z_s.
    dec = DEC0 + 40.0 / 3600
    src = Table(
        {
            "label": [1, 2, 3],
            "ra": [RA0, RA0, RA0],
            "dec": [dec, dec, dec],
            "pa_obs": [92.0, 0.0, 45.0],
            "z_phot": [2.0, 3.0, np.nan],
            "z160": [1.8, 2.7, np.nan],
            "z840": [2.2, 3.3, np.nan],
        }
    )
    out = lc.orientation_table(model, src)
    assert list(out["orientation_class"]) == ["aligned", "anti", "mixed"]
    assert list(out["z_basis"]) == ["photo-z", "photo-z", "grid"]
    np.testing.assert_allclose(out["offset_z2"], [2.0, 90.0, 45.0], atol=1e-6)
    np.testing.assert_allclose(out["offset_zphot"][:2], [2.0, 90.0], atol=1e-6)
    assert np.isnan(out["offset_zphot"][2])
    assert np.all(out["pa_pred_spread"] < 1e-6)
    assert out.meta["provenance"] == "derived"


def test_orientation_uses_the_source_photoz():
    # 7" North of an SIS: at z_s = 1, 2, 4 kappa > 1 (radial stretch, PA 0), but at the source's
    # own z ~ 0.6 kappa < 1 (tangential, PA 90). A tangential source must not be called anti.
    model = _model()
    src = Table(
        {
            "label": [1],
            "ra": [RA0],
            "dec": [DEC0 + 7.0 / 3600],
            "pa_obs": [90.0],
            "z_phot": [0.6],
            "z160": [0.58],
            "z840": [0.62],
        }
    )
    out = lc.orientation_table(model, src)
    assert out["offset_z1"][0] > 85  # the grid would call it anti
    assert out["z_basis"][0] == "photo-z"
    assert out["orientation_class"][0] == "aligned"


def test_class_stats_null_follows_threshold():
    rows = Table({"orientation_class": ["aligned"] * 6 + ["anti"] * 2 + ["mixed"] * 2})
    s30, s45 = lc.class_stats(rows, 30.0), lc.class_stats(rows, 45.0)
    assert (s30["n"], s30["aligned"], s30["anti"]) == (10, 6, 2)
    assert s30["p_random"] == 1 / 3 and s45["p_random"] == 0.5
    assert s30["p_aligned_excess"] < s45["p_aligned_excess"]


def test_read_sigpos(tmp_path):
    par = tmp_path / "input.par"
    par.write_text("image\n    multfile 22 arcs.dat\n    sigposArcsec  0.44\n    end\n")
    assert lc.read_sigpos(par) == 0.44
    par.write_text("image\n    end\n")
    with pytest.raises(SystemExit):
        lc.read_sigpos(par)


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


def test_forced_class_thresholds():
    assert lc.forced_class(np.nan, 50.0) == "no_reference"
    assert lc.forced_class(4.0, 100.0, 50.0) == "undetectable"
    assert lc.forced_class(12.0, 15.0, 1.3) == "recovered"
    assert lc.forced_class(12.0, 600.0, 100.0) == "confused"  # a cluster galaxy dominates
    assert lc.forced_class(12.0, 2.0) == "absent"
    assert lc.forced_class(7.0, 4.0) == "ambiguous"
    assert lc.forced_class(100.0, 6.0, 0.06) == "absent"  # only a far fainter source is there
    assert lc.forced_class(7.0, 6.0, 0.1) == "ambiguous"


def test_aperture_snr_rejects_mostly_invalid_pixels():
    yy, xx = np.mgrid[-50:51, -50:51] * 0.03
    img, err = np.ones(xx.shape), np.full(xx.shape, 0.1)
    f, e = lc.aperture_snr(img, err, xx, yy, 0.0, 0.0)
    assert abs(f) < 1e-9 and e > 0
    img[np.hypot(xx, yy) <= 0.2] = np.nan  # a gap over the aperture
    assert all(np.isnan(lc.aperture_snr(img, err, xx, yy, 0.0, 0.0)))


def _stamp_factory(sources, pix=0.03, half=1.5, noise=0.003):
    """``stamp(ra, dec)`` over a noiseless sky of Gaussian sources given as (ra, dec, flux)."""
    hp = int(np.ceil(half / pix))
    yy, xx = np.mgrid[-hp : hp + 1, -hp : hp + 1] * pix
    cosd = np.cos(np.deg2rad(DEC0))

    def stamp(ra, dec):
        img = np.zeros(xx.shape)
        for sra, sdec, flux in sources:
            dx = -(sra - ra) * cosd * 3600.0  # +x = West, as in the model frame
            dy = (sdec - dec) * 3600.0
            g = np.exp(-((xx - dx) ** 2 + (yy - dy) ** 2) / (2 * 0.06**2))
            img += flux * g / g.sum() if g.sum() > 0 else 0.0
        return img, np.full(xx.shape, noise), xx, yy

    return stamp


def test_forced_check_recovers_offset_image_and_flags_absent_one():
    cosd = np.cos(np.deg2rad(DEC0))
    ref = (RA0, DEC0)
    pred_a = (RA0 + 20 / 3600 / cosd, DEC0)  # an image 0.5" off its prediction
    pred_b = (RA0, DEC0 + 20 / 3600)  # nothing there
    sources = [(ref[0], ref[1], 2.0), (pred_a[0] - 0.5 / 3600 / cosd, pred_a[1], 1.0)]
    stamp = _stamp_factory(sources)
    backtrace = Table(
        {
            "system": ["1", "1"],
            "ra": [ref[0], ref[0] + 40 / 3600 / cosd],
            "dec": [ref[1], ref[1]],
            "magnification": [10.0, 400.0],  # the second is beyond max_ref_mu and has no flux
        }
    )
    table = Table(
        {
            "system": ["1", "1", "1"],
            "ra": [ref[0], pred_a[0], pred_b[0]],
            "dec": [ref[1], pred_a[1], pred_b[1]],
            "magnification": [10.0, 5.0, -5.0],
            "image_class": ["observed", "missing", "missing"],
        }
    )
    table.meta["provenance"] = "derived"
    lc.forced_check(table, backtrace, stamp, search_arcsec=1.0)
    assert list(table["forced_class"]) == ["", "recovered", "absent"]
    assert abs(table["best_dx"][1] - 0.5) < 0.11 and abs(table["best_dy"][1]) < 0.11
    assert 0.8 < table["flux_ratio"][1] < 1.2  # 1.0 observed against 2.0 x 5/10 predicted
    assert table["pred_snr"][2] > 10 and table["best_snr"][2] < 3
    assert table.meta["forced"]["max_ref_mu"] == 50.0

    def gap_stamp(ra, dec):  # valid at the reference only; the predictions fall in a gap
        img, err, xx, yy = stamp(ra, dec)
        if (ra, dec) != ref:
            img = np.full(img.shape, np.nan)
        return img, err, xx, yy

    off = table.copy()
    lc.forced_check(off, backtrace, gap_stamp)
    assert list(off["forced_class"]) == ["", "off_image", "off_image"]


def test_predict_counter_images_reproduces_an_sis_pair():
    from jwst_anomaly import lensmodel

    model = _model()
    grid = lensmodel.DeflectionGrid.compute(model, half_width=60.0, step=0.5)
    z = 2.0
    imgs = lensmodel.find_images(model, grid, 5.0, 3.0, z)
    bright = imgs[np.abs(imgs["magnification"]) > 0.5]
    assert len(bright) == 2
    images = Table(
        {
            "image_id": ["1.1", "1.2"],
            "system": ["1", "1"],
            "ra": bright["ra"],
            "dec": bright["dec"],
            "z": [z, z],
        }
    )
    bt = lensmodel.backtrace_images(model, images, {})
    shapes = Table(
        {
            "label": [1, 2],
            "ra": np.asarray(bright["ra"]),
            "dec": np.asarray(bright["dec"]),
            "mag": [24.0, 25.0],
        }
    )
    table, unpredicted = lc.predict_counter_images(model, grid, bt, shapes, depth=27.0)
    assert unpredicted == []
    classes = list(table["image_class"])
    assert classes.count("observed") == 2
    assert set(classes) <= {"observed", "demagnified"}  # at most a central demagnified image
    assert table.meta["provenance"] == "derived"
