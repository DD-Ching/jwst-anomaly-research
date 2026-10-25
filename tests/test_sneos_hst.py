import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sneos_hst as sh  # noqa: E402

PIX = 0.128


def gauss(shape, x, y, s=0.9):
    yy, xx = np.indices(shape, dtype=float)
    g = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s**2)
    return g / g.sum()


P1, P2 = (23.0, 20.0), (25.0, 26.5)  # ~0.87" apart


def test_aperture_fluxes_recover_ratio_on_background():
    img = 5.0 + 100 * gauss((47, 47), *P1) + 180 * gauss((47, 47), *P2)
    f1, f2 = sh.aperture_fluxes(img, P1, P2, PIX)
    assert f2 / f1 == pytest.approx(1.8, rel=0.03)


def test_random_sigma_matches_white_noise_aperture_scale():
    rng = np.random.default_rng(1)
    img = rng.normal(0, 1.0, (47, 47))
    s = sh.random_sigma(img, P1, P2, PIX)
    n_ap = np.pi * (sh.APER_ARCSEC / PIX) ** 2
    assert 0.5 * np.sqrt(n_ap) < s < 2.0 * np.sqrt(n_ap)


def test_ratio_with_error_and_detection():
    r = sh.ratio_with_error(10.0, 18.0, 1.0, 1.0)
    assert r["ratio"] == pytest.approx(1.8)
    assert r["sratio"] == pytest.approx(1.8 * np.hypot(0.1, 1 / 18))
    assert r["detected"]
    assert not sh.ratio_with_error(4.0, 18.0, 1.0, 1.0)["detected"]


@pytest.mark.parametrize(
    "r,s,want",
    [(1.805, 0.03, "S1"), (1.382, 0.03, "JWST"), (1.58, 0.06, "both"), (1.0, 0.02, "neither")],
)
def test_compare(r, s, want):
    assert sh.compare(np.log(r), s) == want


def test_compare_undetermined():
    assert sh.compare(float("nan"), 0.1) == "undetermined"


def test_injection_trials_recover_on_flat_image():
    rng = np.random.default_rng(2)
    img = rng.normal(0, 0.01, (61, 61))
    psf = gauss((15, 15), 7, 7)
    tr = sh.injection_trials(img, psf, (30.0, 27.0), (32.0, 33.5), 50.0, PIX)
    assert len(tr) == 6
    assert sum(x["recovered"] for x in tr) >= 5


def test_difference_registers_subpixel_offset():
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [292.955, -26.5747]
    w.wcs.crpix = [24, 24]
    w.wcs.cdelt = [-PIX / 3600.0, PIX / 3600.0]
    rng = np.random.default_rng(3)
    yy, xx = np.indices((47, 47), dtype=float)
    blobs = rng.uniform([1, 5, 5], [5, 42, 42], (15, 3))

    def sky(dx=0.0, dy=0.0):
        return sum(
            a * np.exp(-0.5 * ((xx - x - dx) ** 2 + (yy - y - dy) ** 2) / 1.5**2)
            for a, x, y in blobs
        )

    class T:
        pix = PIX

        def section(self, ra, dec, half):
            return sky(0.2, -0.15) + 0.001 * rng.normal(size=xx.shape), None, w

    sci = sky() + 0.001 * rng.normal(size=xx.shape)
    _, info = sh.difference(sci, w, T(), 292.955, -26.5747, np.zeros(sci.shape, bool), PIX)
    assert info["registered"]
    assert info["shift_arcsec"] < 0.01
