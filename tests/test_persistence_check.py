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
    path: Path,
    mjd: float,
    dx: float,
    sci: np.ndarray,
    sat: np.ndarray | None = None,
    substrt: tuple[int, int] = (1, 1),
) -> str:
    rng = np.random.default_rng(int(mjd * 1e4) % 2**32)
    sci = sci + rng.normal(0, 0.05, sci.shape)
    dq = (
        np.zeros(sci.shape, dtype=np.uint32)
        if sat is None
        else sat.astype(np.uint32) * pc.SATURATED
    )
    primary = fits.PrimaryHDU()
    primary.header.update(
        EXPSTART=mjd, EXP_TYPE="NRC_IMAGE", DETECTOR="NRCB3", FILTER="F150W", PUPIL="CLEAR"
    )
    primary.header.update(SUBSTRT1=substrt[0], SUBSTRT2=substrt[1])
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
    one = np.array([1, 1])
    assert pc.classify(np.array([1.0, 2.0]), np.array([False, False]), one) == "undetected"
    assert pc.classify(np.array([20.0, 1.0]), np.array([True, False]), one) == "persistence"
    assert pc.classify(np.array([20.0, 9.0]), np.array([False, False]), one) == "on_sky"
    assert pc.classify(np.array([20.0, 9.0]), np.array([False, True]), one) == "inconclusive"
    # a detection with no earlier exposure checked cannot be cleared
    assert pc.classify(np.array([20.0, 9.0]), np.array([False, False]), np.array([0, 1])) == (
        "inconclusive"
    )
    assert pc.is_suspect(1.0, 50.0, 0)
    assert not pc.is_suspect(1.0, 5.0, 0)
    assert pc.is_suspect(1.0, 0.0, 3)  # saturation earlier on the pixel
    assert not pc.is_suspect(-1.0, 50.0, 0)


def test_run_flags_afterimage_and_keeps_real_source(tmp_path):
    # A bright saturated star sits on detector pixel P=(100, 100) in exposure 1 (pointing 0). The
    # telescope then dithers by 40 px; exposure 2 shows a 1 % afterimage at P, which maps to a
    # different sky position. A real faint source is present at another sky position in all three
    # exposures (the third dithers by another 40 px); the first cannot be cleared (no prior).
    p = (100.0, 100.0)
    star = _blob(*p, 500.0)
    sat = np.hypot(*(np.mgrid[0:N, 0:N] - np.array(p)[::-1, None, None])) < 1.5
    real_sky = _wcs(0).all_pix2world(60.0, 140.0, 0)
    x2, y2 = _wcs(40).all_world2pix(*real_sky, 0)
    f1 = _write(tmp_path / "e1_cal.fits", 60000.00, 0, star + _blob(60, 140, 2.0), sat)
    f2 = _write(
        tmp_path / "e2_cal.fits", 60000.01, 40, _blob(*p, 5.0) + _blob(float(x2), float(y2), 2.0)
    )
    x3, y3 = _wcs(80).all_world2pix(*real_sky, 0)
    f3 = _write(tmp_path / "e3_cal.fits", 60000.02, 80, _blob(float(x3), float(y3), 2.0))
    ghost_sky = _wcs(40).all_pix2world(*p, 0)
    positions = {
        "ghost": (float(ghost_sky[0]), float(ghost_sky[1])),
        "real": (float(real_sky[0]), float(real_sky[1])),
    }
    per_exp, summ = pc.run(positions, [f1, f2, f3], workers=2)
    verdict = dict(zip(summ["uid"], summ["verdict"], strict=True))
    assert verdict == {"ghost": "persistence", "real": "on_sky"}
    ghost = per_exp[(per_exp["uid"] == "ghost") & per_exp["suspect"]]
    assert len(ghost) == 1 and ghost["prior_sat"][0] > 0
    assert per_exp.meta["provenance"] == "derived" and per_exp.meta["failed"] == []
    real = summ[summ["uid"] == "real"][0]
    assert real["n_detected"] == 3 and real["n_clean"] == 2


def test_unreadable_prior_does_not_clear_a_detection(tmp_path, monkeypatch):
    real_sky = _wcs(0).all_pix2world(60.0, 140.0, 0)
    files = []
    for k, dx in enumerate((0, 40, 80)):
        x, y = _wcs(dx).all_world2pix(*real_sky, 0)
        files.append(
            _write(tmp_path / f"e{k}_cal.fits", 60000 + k / 100, dx, _blob(float(x), float(y), 2.0))
        )
    measure = pc.measure_pixel

    def flaky(uri, *args):
        if uri.endswith("e0_cal.fits") and args[0] > 61:  # e0 read as a prior, not as itself
            raise OSError("S3 timeout")
        return measure(uri, *args)

    monkeypatch.setattr(pc, "measure_pixel", flaky)
    per_exp, summ = pc.run({"real": tuple(float(v) for v in real_sky)}, files, workers=1)
    assert summ["verdict"][0] == "inconclusive"  # e1 and e2 each lost a prior: none clean
    assert list(per_exp["n_prior_failed"]) == [0, 1, 1]
    assert len(per_exp.meta["failed"]) == 2


def test_subarray_origin_maps_prior_pixels(tmp_path):
    # Same pointing, but the second exposure is a subarray starting at full-frame pixel (51, 51):
    # its local pixel (x, y) is full-frame (x + 50, y + 50). The bright star at full-frame P in the
    # first exposure must be found under the afterimage at local P - 50 in the second.
    p = (120.0, 120.0)
    f1 = _write(tmp_path / "full_cal.fits", 60000.00, 0, _blob(*p, 500.0))
    sub = np.zeros((N, N))
    sub[:, :] = _blob(p[0] - 50, p[1] - 50, 5.0)
    f2 = _write(tmp_path / "sub_cal.fits", 60000.01, 0, sub, substrt=(51, 51))
    ghost_sky = _wcs(0).all_pix2world(p[0] - 50, p[1] - 50, 0)  # subarray WCS = full WCS here
    per_exp, summ = pc.run({"ghost": tuple(float(v) for v in ghost_sky)}, [f1, f2], workers=1)
    assert summ["verdict"][0] == "persistence"
