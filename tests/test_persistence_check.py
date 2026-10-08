"""Tests for scripts/persistence_check.py (synthetic dithered _cal exposures)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("persistence_check", _DIR / "persistence_check.py")
pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pc)

RA0, DEC0, SCALE = 24.35, -8.44, 0.031 / 3600
N = 200


def _wcs(dx_pix: float) -> WCS:
    """Pointing moved by dx_pix pixels in RA (a dither)."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [N / 2 + dx_pix, N / 2]
    w.wcs.crval = [RA0, DEC0]
    w.wcs.cdelt = [-SCALE, SCALE]
    return w


def _blob(x: float, y: float, amp: float) -> np.ndarray:
    yy, xx = np.mgrid[0:N, 0:N]
    return amp * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 1.2**2))


def _write(
    path: Path, mjd: float, dx: float, sci: np.ndarray, sat: np.ndarray | None = None
) -> str:
    rng = np.random.default_rng(int(mjd * 1e4) % 2**32)
    sci = sci + rng.normal(0, 0.05, sci.shape)
    dq = (
        np.zeros(sci.shape, dtype=np.uint32)
        if sat is None
        else sat.astype(np.uint32) * pc.SATURATED
    )
    primary = fits.PrimaryHDU()
    primary.header.update(EXPSTART=mjd, DETECTOR="NRCB3", FILTER="F150W", PUPIL="CLEAR")
    hdr = _wcs(dx).to_header()
    hdul = fits.HDUList(
        [
            primary,
            fits.ImageHDU(sci.astype(np.float32), hdr, name="SCI"),
            fits.ImageHDU(np.full(sci.shape, 0.05, np.float32), hdr, name="ERR"),
            fits.ImageHDU(dq, hdr, name="DQ"),
        ]
    )
    hdul.writeto(path)
    return str(path)


def test_classify_and_is_suspect():
    assert pc.classify(np.array([1.0, 2.0]), np.array([False, False])) == "undetected"
    assert pc.classify(np.array([20.0, 1.0]), np.array([True, False])) == "persistence"
    assert pc.classify(np.array([20.0, 9.0]), np.array([False, False])) == "on_sky"
    assert pc.classify(np.array([20.0, 9.0]), np.array([False, True])) == "inconclusive"
    assert pc.is_suspect(1.0, 50.0, 0)
    assert not pc.is_suspect(1.0, 5.0, 0)
    assert pc.is_suspect(1.0, 0.0, 3)  # saturation earlier on the pixel
    assert not pc.is_suspect(-1.0, 50.0, 0)


def test_run_flags_afterimage_and_keeps_real_source(tmp_path):
    # A bright saturated star sits on detector pixel P=(100, 100) in exposure 1 (pointing 0). The
    # telescope then dithers by 40 px; exposure 2 shows a 1 % afterimage at P, which maps to a
    # different sky position. A real faint source is present at another sky position in both.
    p = (100.0, 100.0)
    star = _blob(*p, 500.0)
    sat = np.hypot(*(np.mgrid[0:N, 0:N] - np.array(p)[::-1, None, None])) < 1.5
    real_sky = _wcs(0).all_pix2world(60.0, 140.0, 0)
    x2, y2 = _wcs(40).all_world2pix(*real_sky, 0)
    f1 = _write(tmp_path / "e1_cal.fits", 60000.00, 0, star + _blob(60, 140, 2.0), sat)
    f2 = _write(
        tmp_path / "e2_cal.fits", 60000.01, 40, _blob(*p, 5.0) + _blob(float(x2), float(y2), 2.0)
    )
    ghost_sky = _wcs(40).all_pix2world(*p, 0)
    positions = {
        "ghost": (float(ghost_sky[0]), float(ghost_sky[1])),
        "real": (float(real_sky[0]), float(real_sky[1])),
    }
    per_exp, summ = pc.run(positions, [f1, f2], workers=2)
    verdict = dict(zip(summ["uid"], summ["verdict"], strict=True))
    assert verdict == {"ghost": "persistence", "real": "on_sky"}
    ghost = per_exp[(per_exp["uid"] == "ghost") & per_exp["suspect"]]
    assert len(ghost) == 1 and ghost["prior_sat"][0] > 0
    assert per_exp.meta["provenance"] == "derived"
