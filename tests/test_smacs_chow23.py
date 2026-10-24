import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import smacs_chow23 as sc  # noqa: E402
import spt0546_inputs as sp  # noqa: E402


def _header(n=101):
    h = fits.Header()
    h["NAXIS"], h["NAXIS1"], h["NAXIS2"] = 2, n, n
    h["CTYPE1"], h["CTYPE2"] = "RA---CAR", "DEC--CAR"
    h["CRPIX1"] = h["CRPIX2"] = (n + 1) / 2
    h["CRVAL1"], h["CRVAL2"] = 110.827212062959, -73.4546612333333
    h["CDELT1"], h["CDELT2"] = -sc.PIXEL_ARCSEC / 3600, sc.PIXEL_ARCSEC / 3600
    return h


def test_readings_cover_unit_and_x_sign():
    assert set(sc.READINGS) == {"A", "A_neg_x", "P", "P_neg_x"}
    assert sc.READINGS["A"] == (1.0, 1.0)
    assert sc.READINGS["P_neg_x"] == (sc.PIXEL_ARCSEC, -1.0)


def test_build_model_applies_scale_and_sign_on_car_grid():
    n = 101
    ax = np.full((n, n), 2.0)
    ay = np.full((n, n), 3.0)
    m = sc.build_model(ax, ay, _header(n), sc.PIXEL_ARCSEC, -1.0, "x", "t")
    assert np.isclose(m.pixel_arcsec, sc.PIXEL_ARCSEC, rtol=1e-4)
    dx, dy = m.deflection_xy(np.array([0.0]), np.array([0.0]))
    assert np.allclose(dx, -2.0 * sc.PIXEL_ARCSEC)
    assert np.allclose(dy, 3.0 * sc.PIXEL_ARCSEC)


def test_unit_check_needs_exactly_one_reading():
    one = {"A": {"1": 0.3, "2": 0.5}, "A_neg_x": {"1": 9.0, "2": 8.0}}
    assert sp.unit_check(one, sc.U_MAX_MEDIAN_ARCSEC)["reading"] == "A"
    two = {"A": {"1": 0.3}, "P": {"1": 0.9}}
    assert not sp.unit_check(two, sc.U_MAX_MEDIAN_ARCSEC)["passes"]


def test_count_differences():
    a = {"1": {"n_predicted": 3}, "2": {"n_predicted": 5}, "9": {"n_predicted": 1}}
    b = {"1": {"n_predicted": 3}, "2": {"n_predicted": 3}}
    assert sc.count_differences(a, b) == ["2"]


def test_fixed_z_images_are_the_16_of_five_systems():
    t = sc.fixed_z_images()
    assert len(t) == 16
    assert sorted(set(t["system"])) == ["1", "2", "3", "5", "7"]
