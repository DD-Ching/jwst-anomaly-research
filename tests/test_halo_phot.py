"""Synthetic checks of scripts/halo_phot.py (no network)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import halo_phot as hp  # noqa: E402

PIX = 0.031  # arcsec, NIRCam SW i2d


def _grid(half=3.0):
    n = int(np.ceil(half / PIX))
    j, i = np.mgrid[-n : n + 1, -n : n + 1].astype(float)
    return i * PIX, j * PIX


def _halo(xx, yy, x0=1.0, y0=-0.5, amp=5.0, rs=1.5):
    return amp * np.exp(-np.hypot(xx - x0, yy - y0) / rs)


@pytest.mark.parametrize(
    ("route", "amp", "rel"),
    [
        ("medfilt", 0.5, 0.05),  # halo light in the aperture ~ the source's: 3 % bias
        ("annulus", 0.5, 0.10),  # a plane misses the cusp's curvature: 7 %
        ("medfilt", 5.0, 0.40),  # halo ~10x the source: 0.3 mag residual bias, still << raw (12x)
    ],
)
def test_compact_source_on_halo_recovered(route, amp, rel):
    xx, yy = _grid()
    flux = 30.0
    sci = _halo(xx, yy, amp=amp) + hp.gaussian(xx, yy, 0.0, 0.0, flux)
    err = np.full(sci.shape, 0.01)
    f, e, conc = hp.measure(sci, err, xx, yy, route)
    assert f == pytest.approx(flux, rel=rel)
    assert f >= flux  # residual halo light only ever brightens the image
    assert amp > 1 or conc > hp.MIN_CONCENTRATION  # a strong residual halo also lowers conc
    assert e > 0


def test_halo_alone_adds_flux_without_subtraction():
    xx, yy = _grid()
    sci = _halo(xx, yy, amp=0.5) + hp.gaussian(xx, yy, 0.0, 0.0, 30.0)
    raw = float(np.sum(sci[np.hypot(xx, yy) <= hp.R_AP_ARCSEC]))
    assert raw > 1.5 * 30.0  # the halo is a real bias here, so the test above is not trivial


def test_resolved_source_fails_concentration():
    xx, yy = _grid()
    sci = 2.0 * np.exp(-np.hypot(xx, yy) / 0.5)  # exponential disc, scale 0.5"
    _, _, conc = hp.measure(sci, np.full(sci.shape, 0.01), xx, yy, "annulus")
    assert conc < hp.MIN_CONCENTRATION


def test_inject_check_near_unity():
    xx, yy = _grid()
    sci = _halo(xx, yy, amp=0.5) + hp.gaussian(xx, yy, 0.0, 0.0, 30.0)
    err = np.full(sci.shape, 0.01)
    for route in ("medfilt", "annulus"):
        r = hp.inject_check(sci, err, xx, yy, 30.0, route)
        assert r.shape == (8,)
        assert np.all(np.abs(r - 1) < 0.05)


def test_nan_in_aperture_gives_nan():
    xx, yy = _grid()
    sci = hp.gaussian(xx, yy, 0.0, 0.0, 30.0)
    sci[np.hypot(xx, yy) < 0.05] = np.nan
    f, e, c = hp.measure(sci, np.full(sci.shape, 0.01), xx, yy, "medfilt", recentre=0.0)
    assert np.isnan(f) and np.isnan(e) and np.isnan(c)


def test_unknown_route():
    xx, yy = _grid(0.5)
    with pytest.raises(ValueError):
        hp.measure(np.zeros(xx.shape), np.ones(xx.shape), xx, yy, "kron")
