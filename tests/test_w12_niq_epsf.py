"""Tests for scripts/w12_niq_epsf.py (synthetic stamps, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.modeling import models

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w12_niq_epsf", _DIR / "w12_niq_epsf.py")
ep = importlib.util.module_from_spec(_spec)
sys.modules["w12_niq_epsf"] = ep
_spec.loader.exec_module(ep)


def _psf(half=12):
    yy, xx = np.indices((2 * half + 1, 2 * half + 1))
    p = models.Moffat2D(1, half, half, gamma=1.8, alpha=2.0)(xx, yy)
    return p / p.sum()


def test_place_conserves_flux_and_centres():
    psf = _psf()
    im = ep.place(psf, (60, 60), 30.3, 25.7)
    assert im.sum() == pytest.approx(1.0, rel=5e-3)  # wings shifted past the stamp edge
    yy, xx = np.indices(im.shape)
    assert (im * xx).sum() / im.sum() == pytest.approx(30.3, abs=0.05)


def test_build_epsf_recovers_the_star_profile():
    psf = _psf()
    img = np.full((200, 200), 0.1)
    xy = [(40.2, 50.7), (120.6, 80.1), (150.0, 160.4)]
    for x, y in xy:
        img += 500 * ep.place(psf, img.shape, x, y)
    built, n = ep.build_epsf(img, np.ones_like(img, bool), xy, 12)
    assert n == 3 and np.abs(built - psf).max() < 0.1 * psf.max()


def test_fit_two_epsf_leaves_no_residual_for_two_point_sources():
    psf = _psf()
    img = (
        0.05 + 300 * ep.place(psf, (90, 90), 30.4, 45.2) + 150 * ep.place(psf, (90, 90), 60.7, 44.6)
    )
    fit, res = ep.fit_two_epsf(img, [(30, 45), (61, 45)], psf)
    assert abs(res).max() < 1e-3 * img.max()
    assert fit.x_0_2.value == pytest.approx(60.7, abs=0.05)
