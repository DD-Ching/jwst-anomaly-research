"""Frozen stage V rules of scripts/m1931_deblend.py on synthetic images (offline)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
pytest.importorskip("photutils")
import m1931_deblend as md  # noqa: E402


def _wcs(rot_deg: float = 0.0, scale: float = 0.0312) -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [64.39, -11.90]
    w.wcs.crpix = [100.0, 100.0]
    r = np.radians(rot_deg)
    s = scale / 3600.0
    w.wcs.cd = [[-s * np.cos(r), s * np.sin(r)], [s * np.sin(r), s * np.cos(r)]]
    return w


def test_sky_pa_north_up_and_rotated():
    w = _wcs()
    # +y is North, +x is West on a north-up east-left image: theta 90 -> PA 0, theta 0 -> PA 90
    assert md.sky_pa([100], [100], [90.0], w)[0] == pytest.approx(0.0, abs=1e-3) or md.sky_pa(
        [100], [100], [90.0], w
    )[0] == pytest.approx(180.0, abs=1e-3)
    assert md.sky_pa([100], [100], [0.0], w)[0] == pytest.approx(90.0, abs=1e-3)
    wr = _wcs(30.0)
    pa0 = md.sky_pa([100], [100], [0.0], w)[0]
    par = md.sky_pa([100], [100], [0.0], wr)[0]
    assert abs(((par - pa0) + 90) % 180 - 90) == pytest.approx(30.0, abs=1e-2)


def test_section_box_margin_and_clip():
    w = _wcs()
    ra, dec = w.pixel_to_world_values([50.5, 150.5], [60.5, 120.5])
    assert md.section_box(ra, dec, w, (1000, 1000), 10) == (50, 132, 40, 162)
    assert md.section_box(ra, dec, w, (100, 100), 10) == (50, 100, 40, 100)


def test_build_deblends_two_blobs_on_a_halo():
    rng = np.random.default_rng(1)
    n = 256
    y, x = np.mgrid[:n, :n]
    halo = 400.0 * np.exp(-np.hypot(x - 128, y - 128) / 300.0)  # halo scale ~ 5 boxes, as on a bCG
    blob = lambda x0, y0, f, s: f * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * s * s))  # noqa: E731
    img = halo + blob(100, 128, 80, 3) + blob(110, 128, 80, 3) + blob(180, 60, 40, 2)
    err = np.full((n, n), 1.0)
    sci = img + rng.normal(0, 1.0, (n, n))
    cat, meta, deb = md.build(sci, err, _wcs())
    assert meta["n_segments_deblended"] >= meta["n_segments_detected"]
    assert len(cat["flux"]) >= 3  # the close pair is split, the third source is found
    assert np.all(cat["flux"] > 0) and np.all(cat["snr"] > 0)
    assert set(cat) == {"label", "ra", "dec", "a", "b", "pa", "flux", "snr", "ext", "ell", "scale"}
    assert np.all((cat["pa"] >= 0) & (cat["pa"] < 180))
