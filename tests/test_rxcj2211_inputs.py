"""Offline checks of the RXCJ2211-0350 input-check functions."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import rxcj2211_inputs as rj  # noqa: E402


def test_alpha_from_kappa_isothermal_disc():
    # uniform disc kappa0 of radius R: alpha = kappa0 * r inside, kappa0 * R^2 / r outside
    n, pix, k0, R = 201, 0.5, 0.5, 20.0
    y, x = (np.mgrid[:n, :n] - n // 2) * pix
    kappa = np.where(np.hypot(x, y) <= R, k0, 0.0)
    ax, ay = rj.alpha_from_kappa(kappa, pix)
    c = n // 2
    assert ax[c, c + 20] == pytest.approx(k0 * 10.0, rel=0.03)  # inside, r = 10
    assert ay[c + 60, c] == pytest.approx(k0 * R**2 / 30.0, rel=0.03)  # outside, r = 30
    assert ax[c, c - 20] == pytest.approx(-k0 * 10.0, rel=0.03)


def test_check_s_reads_sign_convention():
    ra0, dec0 = 332.94, -3.83
    cosd = np.cos(np.radians(dec0))
    # deflection 10" east, 5" north; x along -RA means deflect_x = -10
    t = Table(
        {
            "ra": [ra0 + 10 / 3600 / cosd],
            "dec": [dec0 + 5 / 3600],
            "ra_source_best": [ra0],
            "dec_source_best": [dec0],
            "deflect_x_best": [-10.0],
            "deflect_y_best": [5.0],
        }
    )
    s = rj.check_s(t)
    assert s["convention"] == "x_along_minus_ra"
    assert s["passes"]
    t["deflect_x_best"] = [10.0]
    assert rj.check_s(t)["convention"] == "x_along_plus_ra"


def test_compare_thresholds():
    z = np.zeros(10)
    assert rj.compare(z + 0.05, z, z, z)["passes"]
    assert not rj.compare(z + 0.2, z, z, z)["passes"]  # median > 0.10
    e = z.copy()
    e[-2:] = 0.5  # 90th percentile > 0.25
    assert not rj.compare(e, z, z, z)["passes"]


def test_affine_fit_recovers_terms():
    rng = np.random.default_rng(1)
    te, tn = rng.normal(0, 30, 50), rng.normal(0, 30, 50)
    re, rn = 0.3 + 0.01 * te - 0.02 * tn, -0.1 + 0.005 * tn
    ce, cn = rj.affine_fit(te, tn, re, rn)
    fe, fn = rj.affine_apply(ce, cn, te, tn)
    assert np.allclose(fe, re) and np.allclose(fn, rn)


def test_dja_check():
    html = "cluster-a-grizli-v7 rxj2211-grizli-v7 abell2744-grizli-v7"
    d = rj.dja_check(html)
    assert d["matches"] == ["rxj2211"] and d["passes"]
    assert not rj.dja_check("abell2744-grizli-v7")["passes"]


def _hdr(n=101, cdelt=0.5 / 3600):
    from astropy.io import fits

    h = fits.Header()
    h.update(
        CTYPE1="RA---CAR", CTYPE2="DEC--CAR", CRVAL1=332.94, CRVAL2=-3.83,
        CRPIX1=(n + 1) / 2, CRPIX2=(n + 1) / 2, CDELT1=-cdelt, CDELT2=cdelt, ZS=3.0,
    )  # fmt: skip
    return h


def test_sample_axes_signs_and_scaling():
    n = 101
    h = _hdr(n)
    ax = np.ones((n, n))  # +1" along array +x (west, CDELT1 < 0)
    ay = np.full((n, n), 2.0)  # +2" north
    e, nn, inside = rj.sample(ax, ay, h, [332.94], [-3.83], [3.0], 3.0)
    assert inside.all()
    assert e[0] == pytest.approx(-1.0) and nn[0] == pytest.approx(2.0)
    e2, _, _ = rj.sample(ax, ay, h, [332.94], [-3.83], [6.0], 3.0)
    assert e2[0] == pytest.approx(-rj.dls_ds(6.0)[0] / rj.dls_ds(3.0)[0])
    # 100" east is outside a 50" map: NaN and not inside
    ra_out = 332.94 + 100 / 3600 / np.cos(np.radians(-3.83))
    e3, _, inside3 = rj.sample(ax, ay, h, [ra_out], [-3.83], [3.0], 3.0)
    assert not inside3[0] and np.isnan(e3[0])


def test_sample_pixel_position():
    # a gradient along array x: value = pixel index; 10" east of centre is 20 px toward lower x
    n = 101
    h = _hdr(n)
    ax = np.tile(np.arange(n, dtype=float), (n, 1))
    ra = 332.94 + 10 / 3600 / np.cos(np.radians(-3.83))
    e, _, _ = rj.sample(ax, np.zeros((n, n)), h, [ra], [-3.83], [3.0], 3.0)
    assert -e[0] == pytest.approx(50 - 20, abs=0.05)


def test_dls_ds():
    assert 0 < rj.dls_ds(3.0)[0] < rj.dls_ds(10.0)[0] < 1
    with pytest.raises(ValueError):
        rj.dls_ds(0.2)


def test_read_csv_keeps_ids_and_nan():
    text = "id,ra,dec,z_lensing,deflect_x_best\n12.10,1.0,2.0,2.0,3.0\n7.1,1.0,2.0,-1.0,\n"
    t = rj.read_csv(text)
    assert list(t["id"]) == ["12.10", "7.1"] and list(t["system"]) == ["12", "7"]
    assert np.isnan(t["deflect_x_best"][1])


def test_readme_header():
    text = (
        "# omega = 0.3000 , hubble = 0.7000\n# zl = 0.3970\n"
        "# chi2 = 174.99 , rms[arcsec] =  0.61\n# MI_REF  1 x\n"
    )
    d = rj.readme_header(text)
    assert d == {"omega": 0.3, "hubble": 0.7, "zl": 0.397, "chi2": 174.99, "rms[arcsec]": 0.61}
