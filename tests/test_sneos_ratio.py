import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sneos_ratio as sr  # noqa: E402

PIX = 0.04


def gauss(shape, x, y, s=1.5):
    yy, xx = np.indices(shape, dtype=float)
    g = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s**2)
    return g / g.sum()


def wcs(crpix, scale=PIX):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [292.955, -26.5747]
    w.wcs.crpix = crpix
    w.wcs.cdelt = [-scale / 3600.0, scale / 3600.0]
    return w


def test_aperture_pair_recovers_ratio_with_flat_background():
    shape = (81, 81)
    p1, p2 = (40.0, 34.0), (46.0, 52.0)  # ~0.76" apart at 0.04"/px
    img = 5.0 + 1000.0 * gauss(shape, *p1) + 1162.0 * gauss(shape, *p2)
    out = sr.aperture_pair(img, np.full(shape, 1e-3), p1, p2, PIX)
    assert abs(out["ratio"] - 1.162) < 0.02
    assert out["detected"]


def test_register_offset_recovers_shift():
    shape = (101, 101)
    a = sum(gauss(shape, x, y) for x, y in ((30, 40), (70, 60), (50, 20)))
    b = sum(gauss(shape, x + 0.4, y - 0.7) for x, y in ((30, 40), (70, 60), (50, 20)))
    dx, dy = sr.register_offset(a, b, np.zeros(shape, bool))
    assert abs(dx - 0.4) < 0.15 and abs(dy + 0.7) < 0.15


def test_resample_identity_and_shift():
    src = gauss((60, 60), 30.0, 25.0, 3.0)
    same = sr.resample(src, wcs([30, 30]), wcs([30, 30]), src.shape)
    assert np.allclose(same[1:-1, 1:-1], src[1:-1, 1:-1], atol=1e-12)  # edges may round to NaN
    moved = sr.resample(src, wcs([30, 30]), wcs([32, 30]), src.shape)  # dst grid offset by 2 px
    j, i = np.unravel_index(np.nanargmax(moved), moved.shape)
    assert (i, j) == (32, 25)


def test_inject_pair_and_recover():
    shape = (81, 81)
    psf = gauss((21, 21), 10.0, 10.0)
    img = sr.inject_pair(np.zeros(shape), psf, (35.3, 40.6), (8.0, 18.0), 500.0, 0.70)
    out = sr.aperture_pair(img, np.full(shape, 1e-3), (35.3, 40.6), (43.3, 58.6), PIX)
    assert sr.recovered(out["ratio"], 0.70)
    assert np.isclose(img.sum(), 500.0 * 1.70, rtol=1e-3)


def test_decide_and_combine():
    assert sr.decide(np.log(1.162), 0.05) == "agrees"
    assert sr.decide(np.log(1.162) + 0.5, 0.05) == "differs"
    assert sr.decide(np.log(1.162) + 0.25, 0.01) == "agrees"  # below the 0.30 floor
    assert sr.decide(np.log(1.162) + 0.5, 0.2) == "agrees"  # within 3σ
    assert sr.decide(float("nan"), 0.1) == "undetermined"
    m, s = sr.combine(np.array([0.1, 0.2]), np.array([0.1, 0.1]))
    assert np.isclose(m, 0.15) and s >= 0.1 / np.sqrt(2)
