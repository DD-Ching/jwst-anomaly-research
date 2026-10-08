"""Tests for scripts/epoch_difference.py (synthetic two-epoch cutouts)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("epoch_difference", _DIR / "epoch_difference.py")
ed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ed)

RA0, DEC0, SCALE = 24.33, -8.43, 0.031 / 3600
N = 65


def _wcs(dx_pix: float = 0.0) -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [RA0, DEC0]
    w.wcs.crpix = [(N + 1) / 2 + dx_pix, (N + 1) / 2]
    w.wcs.cdelt = [-SCALE, SCALE]
    return w


def _gauss(w: WCS, ra: float, dec: float, amp: float, sigma_pix: float) -> np.ndarray:
    x0, y0 = w.world_to_pixel_values(ra, dec)
    yy, xx = np.mgrid[:N, :N]
    return amp * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sigma_pix**2))


def test_nuclear_brightening_is_centred_and_offset_one_is_not():
    rng = np.random.default_rng(1)
    w1, w2 = _wcs(), _wcs(dx_pix=2.3)  # epoch 2 on a shifted grid
    host = lambda w: _gauss(w, RA0, DEC0, 1.0, 3.0)  # noqa: E731
    e1 = host(w1) + 5.0 + rng.normal(0, 0.002, (N, N))
    nuclear = host(w2) + _gauss(w2, RA0, DEC0, 1.0, 1.2) + 5.0
    res, _ = ed.difference(e1, w1, nuclear, w2, RA0, DEC0)
    assert res["diff_centroid_offset_from_centroid1_arcsec"] < 0.01
    assert res["cumflux_diff"][-1] > 0

    off_ra = RA0 + 0.2 / 3600 / np.cos(np.deg2rad(DEC0))  # 0.2″ east of the nucleus
    offset = host(w2) + _gauss(w2, off_ra, DEC0, 1.0, 1.2) + 5.0
    res2, _ = ed.difference(e1, w1, offset, w2, RA0, DEC0)
    assert 0.15 < res2["diff_centroid_offset_from_centroid1_arcsec"] < 0.25


def test_no_change_gives_small_difference(tmp_path):
    w = _wcs()
    img = _gauss(w, RA0, DEC0, 1.0, 3.0) + 2.0
    paths = []
    for k in (1, 2):
        p = tmp_path / f"e{k}.fits"
        fits.HDUList(
            [fits.PrimaryHDU(), fits.ImageHDU(img, header=w.to_header(), name="SCI")]
        ).writeto(p)
        paths.append(p)
    e1, w1 = ed.read_cutout(paths[0])
    e2, w2 = ed.read_cutout(paths[1])
    res, _ = ed.difference(e1, w1, e2, w2, RA0, DEC0)
    assert abs(res["diff_over_epoch1"][-1]) < 1e-6
    assert (
        ed.main(
            [
                "--ra",
                str(RA0),
                "--dec",
                str(DEC0),
                "--epoch1",
                str(paths[0]),
                "--epoch2",
                str(paths[1]),
            ]
        )
        == 0
    )
