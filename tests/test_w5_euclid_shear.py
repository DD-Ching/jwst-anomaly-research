"""Tests for scripts/w5_euclid_shear.py (synthetic inputs only)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w5_euclid_shear", _DIR / "w5_euclid_shear.py")
wes = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wes)


def test_pa_mapping():
    # position_angle is already PA east of north (MER cutouts); only folded into [0, 180).
    np.testing.assert_allclose(
        wes.pa_east_of_north([0.0, 90.0, -90.0, -45.0]), [0.0, 90.0, 90.0, 135.0]
    )


def test_tangent_plane_axes():
    x, y = wes.tangent_plane([10.01, 10.0], [20.0, 20.01], 10.0, 20.0)
    assert x[0] < 0 and abs(y[0]) < 1e-6 * 3600  # East of the centre is negative x (x West)
    assert abs(x[1]) < 1e-9 and np.isclose(y[1], 36.0, rtol=1e-4)


def _field(n=4000, seed=1):
    rng = np.random.default_rng(seed)
    r = 300.0 * np.sqrt(rng.uniform(0, 1, n))
    a = rng.uniform(0, 2 * np.pi, n)
    e = 0.25 * rng.standard_normal(n) * np.exp(1j * rng.uniform(0, 2 * np.pi, n))
    return r * np.cos(a), r * np.sin(a), np.where(np.abs(e) < 0.9, e, 0.0)


def test_radial_injection_gives_positive_S_and_tangential_negative():
    x, y, e = _field()
    g = wes.radial_shear(x, y, 0.0, 0.0, 60.0, 90.0)
    am = wes.es.ApertureMass(x, y, [0.0], [0.0], 180.0, "pointmass", 90.0, 5)
    s_rad, _ = am.snr(wes.apply_shear(e, g))
    s_tan, _ = am.snr(wes.apply_shear(e, -g))
    s0, _ = am.snr(e)
    assert s_rad[0] > s0[0] + 5 and s_tan[0] < s0[0] - 5


def test_radial_shear_direction_matches_position_angle():
    # A source due North of the centre is stretched North-South (PA 0) by a radial shear.
    g = wes.radial_shear(np.array([0.0]), np.array([100.0]), 0.0, 0.0, 50.0, 10.0)
    assert np.isclose(np.angle(g[0]) / 2, 0.0) and np.isclose(abs(g[0]), 0.25)


def test_star_psf_and_galaxy_shapes():
    stars = Table(
        {"semimajor_axis": [1.0] * 50, "ellipticity": [0.1] * 50, "position_angle": [0.0] * 50}
    )
    psf = wes.star_psf(stars)
    assert psf["n_stars"] == 50 and psf["mean_e1"] > 0  # PA 0 (N-S) is +e1
    gals = Table(
        {"semimajor_axis": [0.5, 4.0], "ellipticity": [0.3, 0.3], "position_angle": [90.0, 90.0]}
    )
    e = wes.galaxy_shapes(gals, psf["sigma_px_median"])
    assert np.isnan(e[0]) and np.isfinite(e[1]) and e[1].real < 0  # PA 90 (E-W) is -e1


def test_moment_pa_with_wcs():
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [10.0, 20.0]
    w.wcs.crpix = [21.0, 21.0]
    w.wcs.cdelt = [-1e-4, 1e-4]  # North up, East left
    yy, xx = np.indices((41, 41)) - 20.0
    # Elongated along +x +y in pixels (toward North-West): PA east of north 135° (i.e. -45°).
    u, v = (xx + yy) / np.sqrt(2), (xx - yy) / np.sqrt(2)
    img = np.exp(-0.5 * (u**2 / 25 + v**2 / 4))
    img[:4] = img[-4:] = 0.0
    img[:4, :] += np.random.default_rng(0).normal(0, 1e-3, (4, 41))
    img[-4:, :] += np.random.default_rng(1).normal(0, 1e-3, (4, 41))
    pa = wes.moment_pa(img, w, 10.0, 20.0, 15.0)
    assert abs(wes.axis_diff(pa, 135.0)) < 2.0


def test_density_limits():
    row = {"flags": [], "injection_efficiency": 0.5, "area_deg2": 0.2}
    fields = {"A": {"theta_e": {f"{t:g}": row for t in wes.THETA_E_ARCSEC}}}
    lim = wes.density_limits(fields)
    assert np.isclose(lim["60"]["n95_deg2"], 30.0)


def test_sheared_catalogue_whole_chain():
    rng = np.random.default_rng(3)
    n = 20000
    gals = Table(
        {
            "semimajor_axis": rng.uniform(3, 6, n),
            "ellipticity": rng.uniform(0, 0.6, n),
            "position_angle": rng.uniform(-90, 90, n),
        }
    )
    g = np.full(n, 0.1 + 0j)  # shear along PA 0 (+e1)
    out = wes.sheared_catalogue(gals, g, rng)
    e = wes.galaxy_shapes(out, 1.0)
    e0 = wes.galaxy_shapes(wes.sheared_catalogue(gals, np.zeros(n, complex), rng), 1.0)
    assert np.nanmean(e.real) > 0.05 and abs(np.nanmean(e0.real)) < 0.01
    # area kept: a * b unchanged
    b_in = gals["semimajor_axis"] * (1 - gals["ellipticity"])
    b_out = out["semimajor_axis"] * (1 - out["ellipticity"])
    np.testing.assert_allclose(out["semimajor_axis"] * b_out, gals["semimajor_axis"] * b_in)


def test_density_limits_gated_by_efficiency():
    blind = {"flags": [], "injection_efficiency": 0.1, "area_deg2": 0.2}
    lim = wes.density_limits({"A": {"theta_e": {f"{t:g}": blind for t in wes.THETA_E_ARCSEC}}})
    assert lim["60"]["n95_deg2"] is None
