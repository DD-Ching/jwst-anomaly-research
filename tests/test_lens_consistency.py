"""Tests for scripts/lens_consistency.py (synthetic model, map and catalogs)."""

from __future__ import annotations

import importlib.util
import json
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
            "system": ["1", "1", "1", "1"],
            "ra": [ref[0], ref[0] + 40 / 3600 / cosd, pred_a[0], pred_b[0]],
            "dec": [ref[1], ref[1], pred_a[1], pred_b[1]],
            "magnification": [10.0, 400.0, 5.0, -5.0],
            "image_class": ["observed", "observed", "missing", "missing"],
            "sep_image_arcsec": [0.1, 0.2, np.nan, np.nan],
        }
    )
    table.meta["provenance"] = "derived"
    lc.forced_check(table, backtrace, stamp, search_arcsec=1.0)
    assert list(table["forced_class"]) == ["", "", "recovered", "absent"]
    assert abs(table["best_dx"][2] - 0.5) < 0.11 and abs(table["best_dy"][2]) < 0.11
    assert 0.8 < table["flux_ratio"][2] < 1.2  # 1.0 observed against 2.0 x 5/10 predicted
    assert table["pred_snr"][3] > 10 and table["best_snr"][3] < 3
    assert table["search_arcsec"][3] == 1.0  # residuals 0.1-0.2": the base radius
    widened = table.copy()
    lc.forced_check(widened, backtrace, stamp, search_arcsec=1.0, unpredicted=["1.3"])
    assert widened["search_arcsec"][3] == 2.25  # a catalogued image of system 1 is unpredicted
    assert table.meta["forced"]["max_ref_mu"] == 50.0

    def gap_stamp(ra, dec):  # valid at the reference only; the predictions fall in a gap
        img, err, xx, yy = stamp(ra, dec)
        if (ra, dec) != ref:
            img = np.full(img.shape, np.nan)
        return img, err, xx, yy

    off = table.copy()
    lc.forced_check(off, backtrace, gap_stamp)
    assert list(off["forced_class"]) == ["", "", "off_image", "off_image"]

    # a reference whose recentred peak is a brighter neighbour (nothing at the catalogued spot)
    nb = _stamp_factory([(ref[0] + 0.35 / 3600 / cosd, ref[1], 2.0)])
    nbt = table.copy()
    lc.forced_check(nbt, backtrace, nb)
    assert list(nbt["forced_class"])[2:] == ["no_reference", "no_reference"]

    # two usable references whose f/|mu| disagree by more than 3x
    bt2 = Table(
        {
            "system": ["1", "1"],
            "ra": [ref[0], ref[0] - 10 / 3600 / cosd],
            "dec": [ref[1], ref[1]],
            "magnification": [10.0, 5.0],
        }
    )
    st2 = _stamp_factory(
        [(ref[0], ref[1], 2.0), (ref[0] - 10 / 3600 / cosd, ref[1], 0.25)], noise=0.001
    )  # f/|mu| 0.2 against 0.05
    t2 = table.copy()
    t2["magnification"][1] = 5.0
    lc.forced_check(t2, bt2, st2)
    assert list(t2["forced_class"])[2:] == ["inconsistent_reference"] * 2


def test_search_radius_follows_model_residuals():
    t = Table(
        {
            "system": ["1", "1", "1", "2"],
            "image_class": ["observed", "observed", "missing", "observed"],
            "sep_image_arcsec": [0.3, 1.2, np.nan, 0.1],
        }
    )
    assert lc.system_search_radius(t, "1", 1.0, False) == pytest.approx(1.8)
    assert lc.system_search_radius(t, "2", 1.0, False) == 1.0
    assert lc.system_search_radius(t, "2", 1.0, True) == pytest.approx(2.25)  # one unpredicted
    t["sep_image_arcsec"][1] = 5.0
    assert lc.system_search_radius(t, "1", 1.0, False) == lc.MAX_SEARCH_ARCSEC
    assert lc.system_search_radius(t, "2", 4.0, False) == 4.0  # a user radius is never capped


def test_resolved_reference_is_not_used():
    # a reference spread over ~0.5" (resolved): its 0.2" aperture flux does not scale with |mu|
    cosd = np.cos(np.deg2rad(DEC0))
    yy, xx = np.mgrid[-60:61, -60:61] * 0.03

    def stamp(ra, dec):
        dx = -(RA0 - ra) * cosd * 3600.0
        dy = (DEC0 - dec) * 3600.0
        img = np.exp(-((xx - dx) ** 2 + (yy - dy) ** 2) / (2 * 0.35**2))
        return img, np.full(xx.shape, 0.001), xx, yy

    bt = Table({"system": ["1"], "ra": [RA0], "dec": [DEC0], "magnification": [5.0]})
    t = Table(
        {
            "system": ["1", "1"],
            "ra": [RA0, RA0 + 20 / 3600 / cosd],
            "dec": [DEC0, DEC0],
            "magnification": [5.0, 3.0],
            "image_class": ["observed", "missing"],
            "sep_image_arcsec": [0.1, np.nan],
        }
    )
    lc.forced_check(t, bt, stamp)
    assert t["forced_class"][1] == "no_reference"


def test_frame_offset_moves_model_and_images_together():
    model = _model()
    images = Table({"ra": [RA0], "dec": [DEC0]})
    lc.MODELS["_test"] = {"frame_offset_arcsec": (0.3, -0.1)}
    try:
        x0, y0 = model.to_frame(RA0, DEC0)
        lc.apply_frame_offset("_test", model, images)
        x1, y1 = model.to_frame(images["ra"][0], images["dec"][0])
        # the image keeps its model-frame position; on the sky it moved by the offset
        assert abs(x1 - x0) < 1e-6 and abs(y1 - y0) < 1e-6
        dra = (images["ra"][0] - RA0) * np.cos(np.deg2rad(DEC0)) * 3600
        assert abs(dra - 0.3) < 1e-6 and abs((images["dec"][0] - DEC0) * 3600 + 0.1) < 1e-6
    finally:
        del lc.MODELS["_test"]


def test_frame_offset_moves_a_map_model_with_its_maps():
    from jwst_anomaly import lensmodel

    model = _model()
    n, pix = 401, 0.2
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [RA0, DEC0]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-pix / 3600, pix / 3600]
    jj, ii = np.mgrid[0:n, 0:n]
    ax, ay = model.deflection_xy(*model.to_frame(*w.pixel_to_world_values(ii, jj)))
    mm = lensmodel.MapLensModel(ax, ay, w, model.z_lens, model.cosmology, source="synthetic")
    ref = lensmodel.MapLensModel(ax, ay, w.deepcopy(), model.z_lens, model.cosmology, source="ref")
    px, py = np.array([7.3, -11.2]), np.array([4.1, 9.6])
    before = mm.deflection_xy(px, py)
    sky = mm.to_sky(px, py)
    lc.MODELS["_test"] = {"frame_offset_arcsec": (0.6, -0.4)}
    try:
        lc.apply_frame_offset("_test", mm)
    finally:
        del lc.MODELS["_test"]
    # in model-frame coordinates nothing changes: the maps moved with the reference point
    np.testing.assert_allclose(mm.deflection_xy(px, py), before, atol=1e-6)
    # the maps moved +0.6" in RA cos dec and -0.4" in Dec, so a fixed sky position now shows
    # what the unshifted maps had 0.6" lower in RA and 0.4" higher in Dec
    ra_src = sky[0] - 0.6 / 3600.0 / np.cos(np.deg2rad(DEC0))
    dec_src = sky[1] + 0.4 / 3600.0
    np.testing.assert_allclose(
        mm.deflection_xy(*mm.to_frame(*sky)),
        ref.deflection_xy(*ref.to_frame(ra_src, dec_src)),
        atol=2e-3,
    )


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


def test_model_sigpos_rules(tmp_path):
    images = Table({"a": [0.38, 0.57]})
    inp = tmp_path / "input.par"
    inp.write_text("image\n    sigposArcsec 0.44\n    end\n")
    files = {"input.par": inp}
    np.testing.assert_allclose(lc.model_sigpos("smacs0723-iclv2", files, images), [0.44, 0.44])
    np.testing.assert_allclose(
        lc.model_sigpos("abell2744-bergamini23", files, images), [0.38, 0.57]
    )
    np.testing.assert_allclose(lc.model_sigpos("elgordo-caminha23", files, images), [0.621] * 2)


def test_imageplane_check_summary():
    model = _model()
    theta_e = model.components[0].b0 * float(model.dls_ds(2.0))
    grid = lc.lensmodel.DeflectionGrid.compute(model, half_width=1.5 * theta_e, step=0.5)
    ra, dec = model.to_sky(np.array([2.0 + theta_e, 2.0 - theta_e]), np.zeros(2))
    images = Table(
        {"image_id": ["1.1", "1.2"], "system": ["1", "1"], "ra": ra, "dec": dec, "z": [2.0, 2.0]}
    )
    bt = lc.lensmodel.backtrace_images(model, images, {})
    assert lc.grid_half_width(model, images) >= theta_e + 20
    # x = 2 +- theta_E is exact only for a pure SIS; this cut-off profile leaves ~0.05".
    ip, summary = lc.imageplane_check(model, grid, bt, np.array([0.1, 0.1]), 1.0)
    assert summary["n_solved"] == 2 and summary["rms_dtheta_arcsec"] < 0.1
    assert summary["images_over_3sigma"] == [] and summary["chi2_pos_lenstool"] == 1.0


@pytest.mark.network
@pytest.mark.parametrize(
    "name, n_images, n_families, chi2",
    [
        ("elgordo-caminha23", 56, 23, 80.221558),
        ("abell2744-bergamini23", 149, 50, 146.604318),
        ("macs0416-canucs", 303, None, 344.298308),
        ("abell370-canucs", 115, None, 192.611321),
    ],
)
def test_cluster_model_files_load(name, n_images, n_families, chi2):
    files = lc.model_files(name)
    par = lc.lensmodel.parse_lenstool_par(files["best.par"])
    model = lc.lensmodel.LensModel.from_par(par)
    images = lc.lensmodel.load_lenstool_images(files["arcs.dat"])
    assert len(images) == n_images
    if n_families is not None:
        assert len(set(images["system"])) == n_families
    assert lc.chi2pos_from_par(files["best.par"]) == pytest.approx(chi2)
    z = lc.lensmodel.image_redshifts(images, par["z_m_limit"])
    assert np.all(np.isfinite(z))  # every image has a catalogued or fixed redshift
    assert len(model.components) == len(par["potentials"])


def _flux_inputs(offsets_mag=(0.0, 0.0, 0.0), npix=(500, 500, 500), z_lo=(1.0, 1.0, 1.0)):
    model = _model()
    ra = RA0 + np.array([50.0, -20.0, 0.0]) / 3600 / np.cos(np.deg2rad(DEC0))
    dec = DEC0 + np.array([0.0, 0.0, 60.0]) / 3600
    images = Table(
        {
            "image_id": ["1.1", "1.2", "1.3"],
            "system": ["1", "1", "1"],
            "ra": ra,
            "dec": dec,
            "a": [0.5] * 3,
            "z": [2.0] * 3,
        }
    )
    mu = np.asarray(model.evaluate(ra, dec, 2.0)["magnification"], float)
    mag = 25.0 - 2.5 * np.log10(np.abs(mu)) + np.asarray(offsets_mag)
    phot = Table(
        {
            "dja_id": [11, 12, 13],
            "ra": ra,
            "dec": dec + 0.1 / 3600,
            "mag_auto": mag,
            "magerr_auto": [0.02] * 3,
            "npix": list(npix),
            "f150w_mag": mag + 0.5,
            "f150w_mag_err": [0.02] * 3,
            "f444w_mag": mag - 0.5,
            "f444w_mag_err": [0.02] * 3,
            "z_phot": [2.0] * 3,
            "z025": list(z_lo),
            "z975": [3.0] * 3,
        }
    )
    return model, images, phot, mu


def test_flux_ratio_table_consistent_and_injected_outlier():
    model, images, phot, mu = _flux_inputs()
    assert np.all(np.abs(mu) < 20)
    t = lc.flux_ratio_table(model, images, {}, phot)
    np.testing.assert_allclose(t["resid_mag"], 0.0, atol=1e-6)
    assert list(t["flux_class"]) == ["consistent"] * 3
    assert list(t["colour_class"]) == ["consistent"] * 3
    assert t.meta["provenance"] == "derived"
    # A 1.5 mag deficit in one image: leave-one-out keeps it out of its own reference.
    model, images, phot, _ = _flux_inputs(offsets_mag=(1.5, 0.0, 0.0))
    t = lc.flux_ratio_table(model, images, {}, phot)
    assert t["resid_mag"][0] == pytest.approx(1.5)
    assert list(t["flux_class"]) == ["flux_outlier", "consistent", "consistent"]
    assert list(t["colour_class"]) == ["consistent"] * 3  # achromatic: colours still agree
    # 1.6 mag shifts the good images' references by 0.8 mag; only the worst image is flagged.
    model, images, phot, _ = _flux_inputs(offsets_mag=(1.6, 0.0, 0.0))
    t = lc.flux_ratio_table(model, images, {}, phot)
    np.testing.assert_allclose(t["resid_mag"], [1.6, -0.8, -0.8])
    assert list(t["flux_class"]) == ["flux_outlier", "consistent", "consistent"]
    np.testing.assert_allclose(t["resid_err"], np.hypot(0.02, 0.02 / np.sqrt(2)))
    # A chromatic change: only image 1.2 is 0.5 mag redder.
    model, images, phot, _ = _flux_inputs()
    phot["f444w_mag"][1] -= 0.5
    t = lc.flux_ratio_table(model, images, {}, phot)
    assert list(t["colour_class"]) == ["consistent", "colour_outlier", "consistent"]
    assert t["colour_resid"][1] == pytest.approx(0.5)


def test_flux_ratio_table_pair_residual_is_the_full_difference():
    model, images, phot, _ = _flux_inputs(offsets_mag=(1.0, 0.0, 0.0))
    t = lc.flux_ratio_table(model, images[:2], {}, phot)
    np.testing.assert_allclose(t["resid_mag"], [1.0, -1.0])
    assert list(t["flux_class"]) == ["flux_outlier"] * 2
    # A single usable image cannot be tested.
    t = lc.flux_ratio_table(model, images[:1], {}, phot)
    assert t["flux_class"][0] == "untested"


def test_flux_ratio_table_exclusions():
    model, images, phot, _ = _flux_inputs(npix=(50000, 500, 500))
    t = lc.flux_ratio_table(model, images, {}, phot)
    assert list(t["large_segment"]) == [True, False, False]
    assert list(t["flux_class"]) == ["untested", "consistent", "consistent"]
    # z025 = 3.5 excludes z = 2 even with the 0.15 (1 + z) margin; z025 = 2.4 does not.
    model, images, phot, _ = _flux_inputs(z_lo=(3.5, 2.4, 1.0))
    t = lc.flux_ratio_table(model, images, {}, phot)
    assert list(t["photoz_excluded"]) == [True, False, False]
    assert list(t["flux_class"]) == ["untested", "consistent", "consistent"]
    # One DJA source nearest to two images: both are blends.
    model, images, phot, _ = _flux_inputs()
    images["ra"][1], images["dec"][1] = images["ra"][0], images["dec"][0] + 0.05 / 3600
    t = lc.flux_ratio_table(model, images, {}, phot)
    assert list(t["blended"]) == [True, True, False]
    assert list(t["flux_class"]) == ["untested"] * 3
    # Beyond match_arcsec nothing matches; a frame offset brings the matches back.
    model, images, phot, _ = _flux_inputs()
    t = lc.flux_ratio_table(model, images, {}, phot, match_arcsec=0.05)
    assert set(t["flux_class"]) == {"untested"} and set(t["dja_id"]) == {-1}
    t = lc.flux_ratio_table(model, images, {}, phot, match_arcsec=0.05, offset_arcsec=(0.0, 0.1))
    assert list(t["dja_id"]) == [11, 12, 13]


@pytest.mark.network
def test_canucs_macs0416_reproduces_its_lenstool_chi2(tmp_path):
    # D-044: the CANUCS best fit's image-plane chi2pos is 344.30 (sigpos 0.49"); ours is within 10 %
    out = tmp_path / "v"
    lc.main(["--model", "macs0416-canucs", "--out", str(out), "validate"])
    ip = json.loads((out / "macs0416-canucs" / "validate.json").read_text())["image_plane"]
    assert ip["n_solved"] == ip["n_images"] > 250
    assert abs(ip["chi2_pos"] / ip["chi2_pos_lenstool"] - 1) < 0.10


def test_fluxratios_refuses_a_gated_image_list():
    # D-044: abell370-canucs is a source-plane fit whose image list is gated off
    with pytest.raises(SystemExit, match="usable multiple-image list"):
        lc.main(["--model", "abell370-canucs", "fluxratios", "--photometry", "unused.fits"])
